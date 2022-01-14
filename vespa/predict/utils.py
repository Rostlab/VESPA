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

from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, Iterator, Tuple
from pyfaidx import Fasta
import numpy as np
import re

NON_WORD_RE = re.compile("[^a-zA-Z0-9]")

from vespa.predict.config import *


def unify_seq_id(id: str) -> str:
    return NON_WORD_RE.sub("_", id)


def parse_fasta_input(input_file: Path, split_char="!", id_field=0):
    def key_function(header):
        if "|" in header:
            header = header.split("|")[1]
        else:
            header = header.split(split_char)[id_field]
        # replace tokens that are mis-interpreted when loading h5
        header = unify_seq_id(header)
        return header

    sequences = dict()

    seqs_file = Fasta(
        str(input_file.absolute()),
        key_function=key_function,
        sequence_always_upper=True,
    )
    for rec in seqs_file:
        seq = str(rec[:])
        sequences[rec.name] = str(
            seq
            if not REPLACE_RARE_AA
            else seq.replace("U", "X").replace("Z", "X").replace("O", "X")
        )
    return sequences


def generate_sav_str(fromAA, pos0idx, toAA):
    return f"{fromAA}{pos0idx}{toAA}"


def generate_mutation_str(id, fromAA, pos0idx, toAA):
    return f"{id}:{generate_sav_str(fromAA, pos0idx, toAA)}"


def parse_mutation_str(mutation: str) -> Tuple[str, str, int, str]:
    parts = NON_WORD_RE.split(mutation.strip())
    sav = parse_sav_str(parts[-1])
    return "_".join(parts[:-1]), *sav


def parse_sav_str(sav: str) -> Tuple[str, int, str]:
    fromAA = sav[0]
    toAA = sav[-1]
    pos0idx = int(sav[1:-1])
    return fromAA, pos0idx, toAA


def pred_to_str(pred: np.float64):
    if RESULT_PREC > 0:
        fmt = f".{RESULT_PREC}f"
        res = f"{pred:fmt}"
    else:
        res = f"{pred}"
    return res


class MutationGenerator:
    def __init__(
        self,
        sequence_dict: Dict[str, str],
        file_path: Path = None,
        order=MUTANT_ORDER,  # ignore X for list of viable mutations
        one_based_file=False,
        rare_aas=False,
    ) -> None:
        self.one_based = one_based_file
        self.sequence_dict = sequence_dict
        self.load_from_file = False
        self.order = order
        if file_path is not None:
            self.load_from_file = True
            self.mutation_file_path = file_path
            self.mutation_dict = self._parse_mutations_from_file(
                self.mutation_file_path, self.one_based
            )
        self.total_mut_num = self.__calc_total_mutation_num()
        self.aa2idx = self._get_aa2idx()

    def _get_aa2idx(self):
        return {aa: idx for idx, aa in enumerate(self.order)}

    def __calc_total_mutation_num(self):
        num = 0
        if self.load_from_file:
            for id in self.mutation_dict:
                num += len(self.mutation_dict[id])
        else:
            for id in self.sequence_dict:
                num += len(self.sequence_dict[id]) * (len(self.order) - 1)
                num += self.sequence_dict[id].count("X")

        return num

    def total_mutation_num(self):
        return self.total_mut_num

    def mutation_num_id(self, id):
        """
        Get the number of requested mutations per sequence id

        Args:
            id (str): Fasta Sequence-Id

        Returns:
            int : number of requested mutations per sequence
        """
        if self.load_from_file:
            seq = self.mutation_dict.get(id, None)
            return len(seq) if seq is not None else 0
        else:
            seq = self.sequence_dict.get(id, None)
            return len(seq) * (len(self.order) - 1) if seq is not None else 0

    @staticmethod
    def _parse_mutations_from_file(mutation_file: Path, one_based_file: bool):
        """
        Read the mutations file into a dictionary

        Args:
            mutation_file (Path): Path to the mutations.txt
            one_based_file (bool): Determines if coordinates in the file start with 1

        Returns:
            Dict[str, List[Tuple[int, str,str ]]]: Dictionary mapping sequence id to a  sequence of mutation positions followed by fromAA and toAA
        """

        mutations_dict = defaultdict(list)
        with mutation_file.open("r") as mutations:
            for line in mutations:
                id, fromAA, pos0idx, toAA = parse_mutation_str(line)
                pos0idx = pos0idx if not one_based_file else (pos0idx - 1)
                mutations_dict[id].append((pos0idx, fromAA, toAA))
        return dict(mutations_dict)

    def __generate_all_mutations_from_file(self) -> Iterator[Tuple[str, ...]]:
        """
        Generator function;
        Yield all mutations for the sequences as listed in the given mutations file.

        Yields:
            Iterator[Tuple[str, ...]]: Generator that yields sequence id, followed by mutation information (see `generate_all_seq_mutations_from_file`).
        """
        for id in self.sequence_dict:
            if id in self.mutation_dict:
                for (
                    pos0idx,
                    fromAA,
                    toAA,
                ) in self.__generate_all_seq_mutations_from_file(id):
                    yield id, pos0idx, fromAA, toAA

    def __generate_all_seq_mutations_from_file(self, id: str):
        for mutation in self.mutation_dict[id]:
            yield mutation

    def __generate_all_mutations(self) -> Iterator[Tuple[str, ...]]:
        """Generate all viable mutations for every sequence in the given fasta file.

        Returns:
            Iterable[Tuple[str, ...]]: sequence-id followed by mutations (see `generate_all_seq_mutations`)

        Yields:
            Iterator[Iterable[Tuple[str, ...]]]: [description]
        """
        for id, sequence in self.sequence_dict.items():
            for pos0idx, fromAA, toAA in self.__generate_all_seq_mutations(id):
                yield id, pos0idx, fromAA, toAA

    def __generate_all_seq_mutations(self, id: str) -> Iterator[Tuple[int, str, str]]:
        """
        Generator function;
        Generate all possible mutations for a given fasta sequence-id.

        Args:
            id (str): the sequence id in the given fasta file.


        Yields:
            Iterator[Tuple[int, str, str]]: Generator for all possible mutations at each sequence postition; format: 0-based index, fromAA, toAA. Does not gernate self mutations.
        """

        sequence = self.sequence_dict[id]
        for pos0idx in range(len(sequence)):
            for toAA in self.order:
                fromAA = (
                    sequence[pos0idx]
                    if not REPLACE_RARE_AA
                    else sequence[pos0idx]
                    .replace("U", "X")
                    .replace("Z", "X")
                    .replace("O", "X")
                )
                if fromAA != toAA:
                    yield pos0idx, fromAA, toAA

    def generate_all_mutations(self):
        """Generate all mutations either from a mutations file if given or all mutations.

        Returns:
            Iterable[Tuple[str, int, str, str]]: Sequence-Id, 0-based sequence postion, wild-type AA, mutant AA
        """
        if self.load_from_file:
            return self.__generate_all_mutations_from_file()
        else:
            return self.__generate_all_mutations()

    def generate_mutations_for_seq(self, seq_id):
        """
        Retrun generator for all mutations for the given sequence id either from file if given or just all possible mutations

        Args:
            seq_id (str): The Fasta id of the required sequence

        Returns:
            Generator[Tuple[int, str, str]] : Generator that generates the relevant sequence 0-based position, fromAA, toAA.
        """
        if self.load_from_file:
            return self.__generate_all_seq_mutations_from_file(seq_id)
        else:
            return self.__generate_all_seq_mutations(seq_id)

    def generate_mutation_mask(self, include_x_pos=False):
        """
        Generator function.
        Generates a mutation mask; essentially a numpy array per amino acid position with 1 on positions where scores should be used. This is mostly relevant if not all mutations per position should be considered.

        Yields:
            Dict[np.array]: dictionary that maps the specific sequence position to the mutations that should be considered.
        """
        old_id = None
        aa2idx = self.aa2idx
        if include_x_pos:
            aa2idx["X"] = len(self.aa2idx)
        idx_filter = defaultdict(lambda: np.zeros(len(aa2idx)))
        # Note: relies on the generators producing mutations in order by id
        # Essentially iterate over mutations until you see new id then return the combined mask.
        for id, pos0idx, fromAA, toAA in self.generate_all_mutations():
            if old_id is None:
                old_id = id
            if id != old_id:
                yield old_id, idx_filter
                old_id = id
                idx_filter = defaultdict(lambda: np.zeros(len(aa2idx)))

            idx_filter[pos0idx][self.aa2idx[toAA]] = 1
        yield old_id, idx_filter

    def generate_mutation_file(self, output_file: Path, percentage=0.2):
        # Sample percentage indices
        indices = np.random.choice(
            self.total_mut_num, (int)(self.total_mut_num * percentage)
        )

        with output_file.open("w") as mutations:
            for idx, (id, pos0idx, fromAA, toAA) in enumerate(
                self.__generate_all_mutations()
            ):
                if idx in indices:
                    mutations.write(
                        f"{generate_mutation_str(id, fromAA, pos0idx, toAA)}\n"
                    )
