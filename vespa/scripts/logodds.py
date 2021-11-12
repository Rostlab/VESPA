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
import argparse

from vespa.predict.config import VERBOSE
from vespa.predict.logodds import *
from vespa.predict import utils


def create_arg_parser():
    """ "Creates and returns the ArgumentParser object."""

    # Instantiate the parser
    parser = argparse.ArgumentParser(
        description=(
            "logodds.py creates T5 logodds for a given text "
            + " file containing sequence(s) in FASTA-format."
        )
    )

    # Required argument
    parser.add_argument(
        "input",
        type=Path,
        metavar="fasta-file",
        help="A path to a fasta-formatted text file containing protein sequence(s).",
    )

    # Optional argument
    parser.add_argument(
        "-o",
        "--output",
        required=False,
        type=Path,
        default=Path("T5_log_odds.h5"),
        help="Output path to store the computed logodds in h5 format. Does not contain the original sequence data (FromAA in SAV format not reconstructable without sequence)! Contains a -1 in the matrix if a specific mutation was not computed.",
    )

    parser.add_argument(
        "--single_csv",
        required=False,
        type=Path,
        help="Output path to store the computed logodds in a single *.csv file. The information is scored as SeqID_SAV;SAVscore.",
    )

    parser.add_argument(
        "--csv_dir",
        required=False,
        type=Path,
        help="Output path to store the computed logodds in multiple csv files (one per sequence). The information is scored as SAV;SAVscore in files named as seqid.csv",
    )

    parser.add_argument(
        "-m",
        "--mutation_file",
        required=False,
        type=Path,
        help="Path to a txt file listing protein ids and their SAVs that should be checked. The file should contain one ProteinID_SAV per line. If not specified all mutations will be considered.",
    )
    parser.add_argument(
        "--one_based_mutations",
        required=False,
        action="store_true",
        help="Flag to determine if mutation file is one based",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        required=False,
        action="store_true",
        help="Control the output of the script.",
    )

    return parser


def create_log_odds(seq_dict, mutation_gen):
    t5_condProbas = T5_condProbas(cache_dir=CACHE_DIR)
    if VERBOSE:
        print("Running DMISS (Deep mutational in-silico scanning.)")
    dmiss_data = t5_condProbas.get_cond_probas(seq_dict, mutation_gen)
    return dmiss_data


def main():
    parser = create_arg_parser()
    args = parser.parse_args()

    VERBOSE = args.verbose

    seq_path = args.input
    out_path = args.output

    seq_dict = utils.parse_fasta_input(seq_path)
    mutation_gen = utils.MutationGenerator(
        sequence_dict=seq_dict,
        file_path=args.mutation_file,
        one_based_file=args.one_based_mutations,
    )

    dmiss_data = create_log_odds(seq_dict, mutation_gen)

    T5_condProbas.write_single_h5(dmiss_data, out_path)

    if args.single_csv:
        T5_condProbas.write_single_csv(dmiss_data, args.single_csv, mutation_gen)
    if args.csv_dir:
        T5_condProbas.write_csv_dir(dmiss_data, args.csv_dir, mutation_gen)