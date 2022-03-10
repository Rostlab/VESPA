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

"""
File: conspred
{Description}
"""

# Default Imports
from pathlib import Path
import argparse
import h5py


# Lib Imports
import torch.utils.data
from vespa.predict.config import MODEL_PATH_DICT, VERBOSE

# Module Imports
from vespa.predict.conspred import ProtT5Cons, get_dataloader


def create_arg_parser():
    """ "Creates and returns the ArgumentParser object."""

    # Instantiate the parser
    parser = argparse.ArgumentParser(
        description=(
            "consurf_infer.py predicts conservation in [0,8] for a given text "
            + " file containing sequence(s) in FASTA-format."
        )
    )

    # Required positional argument
    parser.add_argument(
        "input",
        type=Path,
        help="A path to a h5 embedding file, containing per-residue ProtT5 embeddings.",
    )

    # Optional positional argument
    parser.add_argument(
        "-o",
        "--output_prefix",
        required=False,
        default="./conspred",
        type=str,
        help="A prefix for output files. The output will be written at <prefix>_class.fasta and to ,<prefix>_probs.h5 depending on the selected prediction output.",
    )

    # Optional positional argument
    parser.add_argument(
        "-c",
        "--checkpoint",
        required=False,
        type=Path,
        default=None,
        help="A path for the pre-trained checkpoint for the conservation CNN",
    )

    # Optional argument
    parser.add_argument(
        "--output_probs",
        type=bool,
        default=True,
        action=argparse.BooleanOptionalAction,
        help="Output probabilities for all classes, not only class with highest probability. The probabilities are stored in an h5 file with a dataset per-protein of shape Lx20 (L being the protein length). This output is written to <output_prefix>_probs.h5)",
    )

    # Optional argument
    parser.add_argument(
        "--output_classes",
        type=bool,
        default=False,
        action=argparse.BooleanOptionalAction,
        help="Output the conservation class prediction per residue in a fasta-like format with comma-separated per-residue classes. The output is written to <output_prefix>_class.fast)",
    )

    return parser


def run_conspred(seq_path, checkpoint_path, write_probs, write_classes, out_prefix):
    try:
        if VERBOSE:
            print(f" Start Conservation Prediction ".center(80, "#"))
        embeddings = h5py.File(str(seq_path.resolve()), "r")
        data_loader = get_dataloader(embeddings, batch_size=128)
        if VERBOSE:
            print(f" Load model! ")
        conspred = ProtT5Cons(checkpoint_path)
        if VERBOSE:
            print(f" Predict Conservation! ")
        predictions = conspred.conservation_prediction(
            data_loader, prob_return=write_probs, class_return=write_classes
        )

        out_class = Path(str(out_prefix) + "_class.fasta")
        out_probs = Path(str(out_prefix) + "_probs.h5")
        out_class.parent.mkdir(parents=True, exist_ok=True)

        if write_classes:
            conspred.write_cons_class_pred(predictions, out_class)
        if write_probs:
            conspred.write_probabilities(predictions, out_probs)

        if VERBOSE:
            print(f">> Finished Conservation Prediction!")
    finally:
        embeddings.close()


def main():
    parser = create_arg_parser()
    args = parser.parse_args()

    arguments = {
        "seq_path": args.input,
        "checkpoint_path": args.checkpoint
        if args.checkpoint
        else Path(MODEL_PATH_DICT["CONSCNN"]),
        "out_prefix": args.output_prefix,
        "write_probs": args.output_probs,
        "write_classes": args.output_classes,
    }

    run_conspred(**arguments)


if __name__ == "__main__":
    main()
