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

from vespa.predict.config import VERBOSE, CACHE_DIR
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
        default=Path("T5_logodds.h5"),
        help="Output path to store the computed logodds in h5 format. Does not contain the original sequence data (FromAA in SAV format not reconstructable without sequence)! Contains a -1 in the matrix if a specific mutation was not computed.",
    )

    parser.add_argument(
        "-r",
        "--reconstruction_probas",
        default=False,
        action="store_true",
        help="Output the reconstrcution probabilites of T5 for a masked language model problem. The output is a 20-dimensional vector for all considered residue positions (all positions with at least one mutation or all positions if no file provided) giving the raw probability for each amino acid according to the MUTANT_ORDER.",
    )

    parser.add_argument(
        "--reconstruction_output",
        type=Path,
        default=Path("T5_reconstruction_probas.h5"),
        help="Output path to store the computed reconstruction probabilities in h5 format. Does not contain the original sequence data (FromAA in SAV format not reconstructable without sequence)!",
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
        help="Flag to determine if mutation file is one based. NOTE: all SAV strings generated bt this tool are one based.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        required=False,
        action="store_true",
        help="Control the output of the script.",
    )

    parser.add_argument(
        "--prott5_weights_cache",
        required=False,
        type=Path,
        default=Path(CACHE_DIR),
        help="Path to directory for storing the weights of the language model ProtT5.",
    )

    return parser


def run_logodds(seq_path, out_path, mutation_file, one_based_mutations, reconstruction_probas,
                reconstruction_output, single_csv, csv_dir, cache_dir):
    if VERBOSE:
        print(f' Start Logodds Computation '.center(80, '#'))
    seq_dict = utils.parse_fasta_input(seq_path)
    mutation_gen = utils.MutationGenerator(
        sequence_dict=seq_dict,
        file_path=mutation_file,
        one_based_file=one_based_mutations,
    )

    dmiss_data, probas = create_log_odds(
        seq_dict, mutation_gen, cache_dir, reconstruction_probas
    )

    # TODO maybe save logodds more frequently?
    T5_condProbas.write_single_h5(dmiss_data, out_path)

    if reconstruction_probas:
        T5_condProbas.write_single_h5(probas, reconstruction_output)

    if single_csv:
        T5_condProbas.write_single_csv(dmiss_data, single_csv, mutation_gen)
    if csv_dir:
        T5_condProbas.write_csv_dir(dmiss_data, csv_dir, mutation_gen)

    if VERBOSE:
        print(f'>> Finished Logodds Computation!')


def create_log_odds(seq_dict, mutation_gen, cache_dir, extract_probas=False):
    t5_condProbas = T5_condProbas(cache_dir=cache_dir)
    if VERBOSE:
        print("Running DMISS (Deep mutational in-silico scanning.)")
    proba_dict = t5_condProbas.get_proba_dict(seq_dict, mutation_gen)
    dmiss_data = t5_condProbas.get_log_odds(proba_dict)
    probas = None
    if extract_probas:
        probas = t5_condProbas.get_rec_probabilities(proba_dict)

    return dmiss_data, probas


def main():
    parser = create_arg_parser()
    args = parser.parse_args()

    VERBOSE = args.verbose

    arguments = {
        'seq_path': args.input,
        'out_path': args.output,
        'mutation_file': args.mutation_file,
        'one_based_mutations': args.one_based_mutations,
        'reconstruction_probas': args.reconstruction_probas,
        'reconstruction_output': args.reconstruction_output,
        'single_csv': args.single_csv,
        'csv_dir': args.csv_dir,
        'cache_dir': args.prott5_weights_cache
    }
    run_logodds(**arguments)


if __name__ == '__main__':
    main()
