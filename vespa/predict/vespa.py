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

from typing import Dict, List, Tuple
import pickle as pkl
from collections import defaultdict
from pathlib import Path
import json

from tqdm import tqdm
import h5py
import numpy as np
from Bio.Align import substitution_matrices

from vespa.predict.config import (
    MODEL_PATH_DICT,
    NUM_MUTATION_SCORES,
    OUTPUT_MAP_NAME,
    REVERSE_MUTANT_ORDER,
    VERBOSE,
    VESPAL,
    VESPA,
    CSV_SEP,
    MUTANT_ORDER,
)
from vespa.predict import utils


class VespaPred:
    def __init__(self, vespa=True, vespal=True) -> None:
        self.vespa_pred = vespa
        self.vespal_pred = vespal
        self.subst_matrix = self.__load_subst()
        self.models = self.__load_models()

    def __load_subst(self):
        return substitution_matrices.load("BLOSUM62")

    def write_output_csv(self, output, output_dir: Path):
        if VERBOSE:
            tqdm.write(f"Writing CSV Output to {output_dir.absolute()}")
        if output_dir.exists() and not output_dir.is_dir():
            raise ValueError("Output path must be a directory!")
        else:
            output_dir.mkdir(parents=True, exist_ok=True)
        header_str = f"Mutant{CSV_SEP}{CSV_SEP.join(list(self.models.keys()))}\n"
        map = dict()
        for idx, id in enumerate(output):
            out_path = output_dir / f"{idx}.csv"
            map[idx] = id
            with out_path.open("w") as out_file:
                out_file.write(header_str)
                for mutation in output[id]:
                    info, pred = mutation
                    out_file.write(
                        f"{info}{CSV_SEP}{CSV_SEP.join([utils.pred_to_str(pred[name]) for name in self.models])}\n"
                    )
        with (output_dir / OUTPUT_MAP_NAME).open("w") as map_file:
            json.dump(map, map_file, indent=4)

    def write_output_h5(self, output, seq_2_len: Dict[str, int], output_h5: Path):
        if VERBOSE:
            tqdm.write(f"Writing H5 Output to {output_h5.absolute()}")
        else:
            output_h5.parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(str(output_h5.absolute()), "w") as out_file:
            for seq_idx in output:
                seq_grp = out_file.create_group(seq_idx)
                seq_grp.attrs["residue-order"] = MUTANT_ORDER

                # Get names for models and generate matrices
                results: Dict[str, np.ndarray] = dict.fromkeys(
                    self.models,
                    np.ones((seq_2_len[seq_idx], len(MUTANT_ORDER))) * -1,
                )

                # Write scores to matrix
                for info, pred in output[seq_idx]:
                    fromAA, postion, toAA = utils.parse_sav_str(info)
                    for name in self.models:
                        results[name][postion, REVERSE_MUTANT_ORDER[toAA]] = pred[name]

                # Transfer matrix to file
                for name in self.models:
                    seq_grp.create_dataset(name, data=results[name].T)

    @staticmethod
    def parse_logodds_input(input_file: Path):
        if not input_file.exists() or not input_file.is_file():
            raise ValueError("Logoodds input does not exist.")

        logodds_scores = dict()

        for id, emb in h5py.File(input_file, "r").items():
            id = utils.unify_seq_id(id)
            logodds_scores[id] = np.array(emb)

        return logodds_scores

    @staticmethod
    def parse_cons_input(input_file: Path):
        if not input_file.exists() or not input_file.is_file():
            raise ValueError("Conservation input does not exist.")

        conservation_scores = dict()

        for id, emb in h5py.File(input_file, "r").items():
            id = utils.unify_seq_id(id)
            conservation_scores[id] = np.array(emb)

        return conservation_scores

    @staticmethod
    def load_model(model_name: str):
        if model_name not in MODEL_PATH_DICT:
            raise ValueError(f"Requested model {model_name} not found")

        with open(MODEL_PATH_DICT[model_name], "rb") as file:
            model = pkl.load(file)
        return model

    def __generate_blosum_vec(self, mutation_gen: utils.MutationGenerator):
        vector = np.zeros((mutation_gen.total_mut_num, 1))
        idx = 0
        for (id, pos0idx, fromAA, toAA) in tqdm(
            mutation_gen.generate_all_mutations(),
            total=mutation_gen.total_mut_num,
            desc="Blosum Lookup",
            disable=not VERBOSE,
        ):
            fromAA = (
                fromAA.replace("U", "*")
                .replace("Z", "*")
                .replace("O", "*")
                .replace("X", "*")
            )
            toAA = (
                toAA.replace("U", "*")
                .replace("Z", "*")
                .replace("O", "*")
                .replace("X", "*")
            )
            vector[idx] = self.subst_matrix.get(fromAA).get(toAA)
            idx += 1
        return vector

    def __generate_prob_vec(
        self, logodds_scores, mutation_gen: utils.MutationGenerator
    ):
        vector = np.zeros((mutation_gen.total_mut_num, 1))
        idx = 0
        for (id, pos0idx, fromAA, toAA) in tqdm(
            mutation_gen.generate_all_mutations(),
            total=mutation_gen.total_mut_num,
            desc="Logodds Lookup",
            disable=not VERBOSE,
        ):
            probability = logodds_scores[id][pos0idx, mutation_gen.aa2idx[toAA]]
            vector[idx] = probability
            idx += 1
        return vector

    def __generate_cons_matrix(
        self, conservation_scores, mutation_gen: utils.MutationGenerator
    ):
        vector = np.zeros((mutation_gen.total_mut_num, NUM_MUTATION_SCORES))
        idx = 0
        for (id, pos0idx, fromAA, toAA) in tqdm(
            mutation_gen.generate_all_mutations(),
            total=mutation_gen.total_mut_num,
            desc="Conservation Lookup",
            disable=not VERBOSE,
        ):
            assert id in conservation_scores, f"{id} not in conservation scores"
            assert (
                conservation_scores[id].T.shape[0] > pos0idx
            ), f"Invalid position {pos0idx} for sequence {id}"
            vector[idx] = conservation_scores[id].T[pos0idx]
            idx += 1
        return vector

    def __generate_mutant_info(self, mutation_gen: utils.MutationGenerator):
        info = []
        for (id, pos0idx, fromAA, toAA) in tqdm(
            mutation_gen.generate_all_mutations(),
            total=mutation_gen.total_mut_num,
            desc="Info Generation",
            disable=not VERBOSE,
        ):
            info.append((id, utils.generate_sav_str(fromAA, pos0idx, toAA)))
        return info

    def __generate_input(
        self,
        mutation_gen: utils.MutationGenerator,
        conservation_dict,
        logodds_dict=None,
    ):
        input = dict()
        if self.vespa_pred:
            input["logodds_vec"] = self.__generate_prob_vec(logodds_dict, mutation_gen)
        input["blosum_vec"] = self.__generate_blosum_vec(mutation_gen)
        input["conservation_matrix"] = self.__generate_cons_matrix(
            conservation_dict, mutation_gen
        )
        return self.__generate_input_dict(**input)

    def __generate_input_dict(self, blosum_vec, conservation_matrix, logodds_vec=None):
        input = {}
        if self.vespal_pred:
            input[VESPAL] = np.concatenate((conservation_matrix, blosum_vec), axis=1)
        if self.vespa_pred:
            input[VESPA] = np.concatenate(
                (conservation_matrix, blosum_vec, logodds_vec), axis=1
            )
        return input

    def __load_models(self):
        models = {}
        if self.vespal_pred:
            models[VESPAL] = self.load_model(VESPAL)
        if self.vespa_pred:
            models[VESPA] = self.load_model(VESPA)
        return models

    def generate_predictions(
        self, mutation_gen: utils.MutationGenerator, conservations, logodds=None
    ):
        if self.vespa_pred and logodds is None:
            raise ValueError("Logodds must be present for VESPA predictions")

        input = self.__generate_input(mutation_gen, conservations, logodds)
        # Generate predictions for VESPA or for VESPA and VESPAl
        predictor_out = dict()
        if VERBOSE:
            tqdm.write("Generate Model Predictions")
        for model_name in self.models:
            # Iterate over sub-models (10 per model)
            model_dict = self.models[model_name]
            mod_input = input[model_name]

            # Iterate over data
            model_out = list()
            for fold in model_dict:
                model_out.append(model_dict[fold].predict_proba(mod_input)[:, 1])
            predictor_out[model_name] = np.mean(model_out, axis=0)
        if VERBOSE:
            tqdm.write("Predictions Done; Generate output")
        mutant_info = self.__generate_mutant_info(mutation_gen)
        # Generate output dict
        output = defaultdict(list)
        for idx, (id, mutant_str) in enumerate(mutant_info):
            output[id].append(
                (mutant_str, {name: predictor_out[name][idx] for name in self.models})
            )

        return dict(output)
