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

from vespa.predict.vespa import *
from vespa.predict.config import VERBOSE
from vespa.predict import utils


def setup_argparse():
    parser = argparse.ArgumentParser(
        description="A cli tool for VESPA predictions. Based on the publication 'Embeddings from protein language models predict conservation and variant effects' (https://doi.org/10.21203/rs.3.rs-584804/v1) "
    )
    parser.add_argument(
        "--vespa",
        help="Run the VESPA model",
        default=True,
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "--vespal",
        help="Run the VESPAl model",
        default=True,
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "cons_input",
        type=Path,
        help="Conservation prediction for the relevant sequences.",
    )
    parser.add_argument(
        "--T5_input",
        nargs="?",
        type=Path,
        help="T5 log odds for the relevant sequences. Needs to be present to run VESPA",
    )
    parser.add_argument(
        "fasta_file",
        type=Path,
        help="Fasta sequence file; Needs to be present for mutation generation ",
    )
    parser.add_argument(
        "--output",
        type=Path,
        nargs="?",
        help="Path to output csv file.",
        default=Path("./output/"),
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


def vespa_pred(args):

    seq_path = args.fasta_file
    out_path = args.output

    VERBOSE = args.verbose

    seq_dict = utils.parse_fasta_input(seq_path)
    mutation_gen = utils.MutationGenerator(
        sequence_dict=seq_dict,
        file_path=args.mutation_file,
        one_based_file=args.one_based_mutations,
    )

    vespa = args.vespa
    vespal = args.vespal

    if not vespa and not vespal:
        raise RuntimeError(
            "Requested to run no vespa-model. Please at least allow one model! E.g. `--vespa` or `--vespal`"
        )

    predictor = VespaPred(vespa=vespa, vespal=vespal)

    input = dict()
    input["conservations"] = predictor.parse_cons_input(args.cons_input)
    if vespa:
        input["logodds"] = predictor.parse_logodds_input(args.T5_input)

    output = predictor.generate_predictions(mutation_gen=mutation_gen, **input)
    predictor.write_output(output, out_path)


def main():
    parser = setup_argparse()
    args = parser.parse_args()
    vespa_pred(args)
