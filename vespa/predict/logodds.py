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
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import torch
import h5py
from tqdm import tqdm

from vespa.predict.config import (
    CACHE_DIR,
    CSV_SEP,
    MUTANT_ORDER,
    REPLACE_RARE_AA,
    REVERSE_MUTANT_ORDER,
    TRANSFORMER_LINK,
    SPIECE_UNDERLINE,
    VERBOSE,
    DEVICE,
    LOGODDS
)
from vespa.predict import utils
from vespa.predict.utils_t5 import ProtT5


if VERBOSE:
    print("Using device: {}".format(DEVICE))


@dataclass
class _ProbaVector:
    pos: int
    vector: np.ndarray
    mask: np.ndarray
    wt_vector_pos: int


class T5_condProbas:
    def __init__(self, cache_dir):
        self.prott5 = ProtT5(cache_dir)
        self.tokenizer = self.prott5.get_tokenizer()
        self.AAs = MUTANT_ORDER + "X"
        self.AA2class = {AA: idx for idx, AA in enumerate(self.AAs)}
        self.class2AA = {idx: AA for idx, AA in enumerate(self.AAs)}
        self.aaIdcs = [
            self.tokenizer.get_vocab()["{}{}".format(SPIECE_UNDERLINE, residue)]
            for residue in self.AAs
        ]
        self.softmax = torch.nn.Softmax(dim=0)

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
        ]  # extract for the masked position only AA logits and potentially X if rare_AAs are replaced
        return self.softmax(aa_logits)

    def get_proba_dict(self, seq_dict, mutation_gen: utils.MutationGenerator):
        """
        Compute for all residues in a protein the conditional probabilities for reconstructing single, masked tokens.
        """

        self.model, self.tokenizer = self.prott5.get_model(LOGODDS)

        result_dict = dict()

        for identifier, position_list in tqdm(
            mutation_gen.generate_mutation_mask(include_x_pos=REPLACE_RARE_AA),
            desc="Extract Sequence Logodds",
            disable=not VERBOSE,
        ):

            wt_seq = seq_dict[identifier]
            seq_len = len(wt_seq)
            result_dict[identifier] = seq_len, list()

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
                result_vec = _ProbaVector(
                    pos=residue_idx,
                    vector=aa_probas.detach().cpu().numpy().squeeze(),
                    mask=filter_mask,
                    wt_vector_pos=self.AA2class[wt_seq[residue_idx]],
                )
                result_dict[identifier][1].append(result_vec)

        return result_dict

    @staticmethod
    def get_rec_probabilities(
        probabiliy_dict: Dict[str, Tuple[int, List[_ProbaVector]]]
    ):
        result_dict = dict()
        for id in probabiliy_dict:
            seq_len = probabiliy_dict[id][0]
            array = np.ones((seq_len, len(MUTANT_ORDER))) * -1
            for prob_vec in probabiliy_dict[id][1]:
                # Remove X reconstruction prob as that probability makes little sense
                vec = prob_vec.vector[:-1]
                array[prob_vec.pos] = vec
            result_dict[id] = array
        return result_dict

    @staticmethod
    def get_log_odds(probabiliy_dict: Dict[str, Tuple[int, List[_ProbaVector]]]):
        result_dict = dict()
        for id in probabiliy_dict:
            seq_len = probabiliy_dict[id][0]
            array = np.ones((seq_len, len(MUTANT_ORDER))) * -1
            for prob_vec in probabiliy_dict[id][1]:
                log_vector = np.log(prob_vec.vector)
                log_vector -= log_vector[prob_vec.wt_vector_pos]
                # Remove X reconstruction prob as that probability makes little sense
                vec = log_vector * prob_vec.mask
                vec = vec[:-1]
                array[prob_vec.pos] = vec
            result_dict[id] = array
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
