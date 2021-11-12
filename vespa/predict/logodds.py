#!/usr/bin/env python
"""
 Copyright (C) 2021 Rostlab
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as
 published by the Free Software Foundation, either version 3 of the
 License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.
 
 You should have received a copy of the GNU Affero General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

# -*- coding: utf-8 -*-

import argparse
import time
from pathlib import Path

import numpy as np
import torch
import h5py
from tqdm import tqdm
from transformers import T5ForConditionalGeneration, T5Tokenizer

from vespa.predict.config import (
    CACHE_DIR,
    CSV_SEP,
    MUTANT_ORDER,
    REPLACE_RARE_AA,
    TRANSFORMER_LINK,
    SPIECE_UNDERLINE,
    VERBOSE,
    DEVICE,
)
from vespa.predict import utils


if VERBOSE:
    print("Using device: {}".format(DEVICE))


class T5_condProbas:
    def __init__(self, cache_dir):
        self.cache_dir = cache_dir
        self.tokenizer = self.get_tokenizer()
        self.AAs = MUTANT_ORDER
        self.AA2class = {AA: idx for idx, AA in enumerate(self.AAs)}
        self.class2AA = {idx: AA for idx, AA in enumerate(self.AAs)}
        self.aaIdcs = [
            self.tokenizer.get_vocab()["{}{}".format(SPIECE_UNDERLINE, AA)]
            for AA in self.AAs
        ]
        self.softmax = torch.nn.Softmax(dim=0)

    def get_model(self):
        model = T5ForConditionalGeneration.from_pretrained(
            TRANSFORMER_LINK, cache_dir=self.cache_dir
        )
        model = model.eval()
        model = model.to(DEVICE)
        vocab = T5Tokenizer.from_pretrained(
            TRANSFORMER_LINK, do_lower_case=False, cache_dir=self.cache_dir
        )

        return model, vocab

    def get_tokenizer(self):
        vocab = T5Tokenizer.from_pretrained(
            TRANSFORMER_LINK, do_lower_case=False, cache_dir=self.cache_dir
        )
        return vocab

    def reconstruct_sequence(self, probs):
        return [self.class2AA[yhat] for yhat in probs.argmax(axis=1)]

    def measure_reconstruction_accuracy(self, wt_seq, probs):
        predictions = self.reconstruct_sequence(probs)
        correct = sum(
            [1 for idx, pred in enumerate(predictions) if pred == wt_seq[idx]]
        )
        return correct / len(wt_seq)

    def get_mask_prob(self, masked_seq, mut_idx):
        input_ids = self.tokenizer(masked_seq, return_tensors="pt").input_ids.to(DEVICE)
        try:  # returns a tuple with loss and logits
            with torch.no_grad():  # only logits is useful for us
                output = self.model(input_ids, labels=input_ids).logits
        except RuntimeError as e:
            print("RuntimeError for input with shape: {}".format(input_ids.shape))
            print(e)
            return None
        aa_logits = output[
            0, mut_idx, self.aaIdcs
        ]  # extract for the masked position only logits of 20 std. AAs
        return self.softmax(aa_logits)

    def get_cond_probas(self, seq_dict, mutation_gen: utils.MutationGenerator):
        """
        Compute for all residues in a protein the conditional probabilities for reconstructing single, masked tokens.
        """

        self.model, self.tokenizer = self.get_model()

        result_dict = dict()
        overall_time = time.time()

        for identifier, position_list in tqdm(
            mutation_gen.generate_mutation_mask(),
            desc="Extract Sequence Logodds",
            disable=not VERBOSE,
        ):
            sample_start = time.time()

            wt_seq = seq_dict[identifier]
            seq_len = len(wt_seq)
            result_dict[identifier] = np.ones((seq_len, len(MUTANT_ORDER))) * -1

            for residue_idx, filter_mask in position_list.items():

                # Mask the residue in question
                masked_seq = list(wt_seq)
                masked_seq[
                    residue_idx
                ] = "<extra_id_0>"  # replace AA by special mask token used by T5
                masked_seq = " ".join(masked_seq)

                # Extract probabilities for all mutations
                aa_probas = self.get_mask_prob(masked_seq, residue_idx)
                if aa_probas is None:
                    print(
                        "Skipping {} (L={}) due to RunTimeError.".format(
                            identifier, seq_len
                        )
                    )
                    continue
                result_vec = aa_probas.detach().cpu().numpy().squeeze() * filter_mask
                result_dict[identifier][residue_idx] = result_vec

            accuracy = self.measure_reconstruction_accuracy(
                wt_seq, result_dict[identifier]
            )
            if VERBOSE:
                print(
                    "Acc={:.3f} (time={:.2f}[s] L={})".format(
                        accuracy, time.time() - sample_start, seq_len
                    )
                )

        end = time.time()
        if VERBOSE:
            print("\n############# STATS #############")
            print("Total number of proteins: {}".format(len(result_dict)))
            print(
                "Total time: {:.2f}[s]; time/prot: {:.2f}[s]".format(
                    end - overall_time, (end - overall_time) / len(result_dict)
                )
            )
        return result_dict

    @staticmethod
    def write_single_h5(dmiss_data, h5_path):
        with h5py.File(str(h5_path), "w") as hf:
            for seq_id, log_odds in dmiss_data.items():
                # noinspection PyUnboundLocalVariable
                hf.create_dataset(seq_id, data=log_odds)
        return None

    @staticmethod
    def write_single_csv(
        dmiss_data, csv_path: Path, mutation_gen: utils.MutationGenerator
    ):
        with csv_path.open("w") as out_file:
            out_file.write(f"SeqID_SAV{CSV_SEP}SAV_score\n")
            for id, mutations in T5_condProbas.__gen_mutation_lists(
                dmiss_data, mutation_gen
            ):
                lines = [
                    f"{id}_{sav}{CSV_SEP}{score}\n" for sav, score in mutations.items()
                ]
                out_file.writelines(lines)

    @staticmethod
    def __gen_mutation_lists(dmiss_data, mutation_gen: utils.MutationGenerator):
        AA2idx = mutation_gen._get_aa2idx()
        mutation_dict = dict()
        old_id = None
        for id, pos0idx, fromAA, toAA in mutation_gen.generate_all_mutations():
            if old_id is None:
                old_id = id
            if id != old_id:
                yield old_id, mutation_dict
                old_id = id
                mutation_dict = dict()

            mutation_handle = utils.generate_sav_str(fromAA, pos0idx, toAA)
            mutation_score = dmiss_data[id][pos0idx, AA2idx[toAA]]
            mutation_dict[mutation_handle] = mutation_score
        yield old_id, mutation_dict

    @staticmethod
    def write_csv_dir(dmiss_data, mutation_dir: Path, mutation_gen):
        if mutation_dir.exists() and not mutation_dir.is_dir():
            raise ValueError("Output path must be a directory!")
        else:
            mutation_dir.mkdir(parents=True, exist_ok=True)
        for id, mutations in T5_condProbas.__gen_mutation_lists(
            dmiss_data, mutation_gen
        ):
            csv_path = mutation_dir / f"{id}.csv"
            with csv_path.open("w") as out_file:
                out_file.write(f"SAV{CSV_SEP}SAV_score\n")
                lines = [f"{sav}{CSV_SEP}{score}\n" for sav, score in mutations.items()]
                out_file.writelines(lines)
