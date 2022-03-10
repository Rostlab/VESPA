import argparse
from pathlib import Path
import warnings

import sys

vespa_location = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(vespa_location))

from vespa.scripts.embedding import run_embedding
from vespa.scripts.conspred import run_conspred
from vespa.scripts.logodds import run_logodds
from vespa.scripts.vespa_run import run_vespa
from vespa.predict.config import (CACHE_DIR, MODEL_PATH_DICT)


def setup_argparse():
    parser = argparse.ArgumentParser(
        description="A meta script to execute the components of VESPA sequentially to obtian predictions. "
                    "Based on the publication 'Embeddings from protein language models predict conservation and variant effects' (https://doi.org/10.21203/rs.3.rs-584804/v1) "
    )

    parser.add_argument(
        "fasta_input",
        type=Path,
        help="Fasta sequence file; Needs to be present for mutation generation ",
    )

    parser.add_argument(
        "--vespa",
        help="Run the VESPA model",
        default=False,
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "--vespal",
        help="Run the VESPAl model",
        default=True,
        action=argparse.BooleanOptionalAction,
    )

    # Optional argument
    parser.add_argument(
        "--vespa_output_directory",
        required=False,
        type=Path,
        help="Path to directory containing outputs produced by VESPA."
    )

    parser.add_argument(
        "--use_existing_embeddings",
        required=False,
        type=Path,
        help="Path to computed embeddings to use in h5 format. "
             "Does not contain the original sequence data, only the sequence id!",
    )

    # shared arguments
    # between logodds, vespa
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

    # shared arguments
    # between logodds, embedding
    parser.add_argument(
        "--prott5_weights_cache",
        required=False,
        type=Path,
        help="Path to directory for storing the weights of the language model ProtT5.",
    )
    return parser, {'vespa_output_directory': Path("./vespa_run_directory/"),
                    'prott5_weights_cache': CACHE_DIR}


def embedding_argparse(group):
    group.add_argument(
        "--embedding_output",
        required=False,
        type=Path,
        help="Output path to store the computed embeddings in h5 format. "
             "Does not contain the original sequence data, only the sequence id!",
    )

    # return default values for parser
    return {'embedding_output': "T5_embedding.h5"}


def logodds_argparse(group):
    group.add_argument(
        "--logodds_output",
        required=False,
        type=Path,
        help="Output path to store the computed logodds in h5 format. Does not contain the original sequence data (FromAA in SAV format not reconstructable without sequence)! Contains a -1 in the matrix if a specific mutation was not computed.",
    )

    group.add_argument(
        "-r",
        "--reconstruction_probas",
        default=False,
        action="store_true",
        help="Output the reconstrcution probabilites of T5 for a masked language model problem. The output is a 20-dimensional vector for all considered residue positions (all positions with at least one mutation or all positions if no file provided) giving the raw probability for each amino acid according to the MUTANT_ORDER.",
    )

    group.add_argument(
        "--reconstruction_output",
        required=False,
        type=Path,
        help="Output path to store the computed reconstruction probabilities in h5 format. Does not contain the original sequence data (FromAA in SAV format not reconstructable without sequence)!",
    )

    group.add_argument(
        "--single_csv",
        required=False,
        type=Path,
        help="Output path to store the computed logodds in a single *.csv file. The information is scored as SeqID_SAV;SAVscore.",
    )

    group.add_argument(
        "--csv_dir",
        required=False,
        type=Path,
        help="Output path to store the computed logodds in multiple csv files (one per sequence). The information is scored as SAV;SAVscore in files named as seqid.csv",
    )

    return {'logodds_output': "T5_logodds.h5",
            'reconstruction_output': "T5_reconstruction_probas.h5"}


def conspred_argparse(group):
    group.add_argument(
        "--conspred_output_prefix",
        required=False,
        type=Path,
        help="A prefix for output files. The output will be written at <prefix>_class.fasta and to ,<prefix>_probs.h5 depending on the selected prediction output.",
    )

    group.add_argument(
        "-c",
        "--checkpoint",
        required=False,
        type=Path,
        default=None,
        help="A path for the pre-trained checkpoint for the conservation CNN",
    )

    group.add_argument(
        "--output_probs",
        type=bool,
        default=True,
        action=argparse.BooleanOptionalAction,
        help="Output probabilities for all classes, not only class with highest probability. The probabilities are stored in an h5 file with a dataset per-protein of shape Lx20 (L being the protein length). This output is written to <output_prefix>_probs.h5)",
    )

    group.add_argument(
        "--output_classes",
        type=bool,
        default=False,
        action=argparse.BooleanOptionalAction,
        help="Output the conservation class prediction per residue in a fasta-like format with comma-separated per-residue classes. The output is written to <output_prefix>_class.fast)",
    )

    return {'conspred_output_prefix': "conspred"}


def vespa_argparse(group):
    group.add_argument(
        "--vespa_pred_output",
        required=False,
        type=Path,
        nargs="?",
        help="Path to output csv file."
    )

    group.add_argument(
        "--no_csv",
        help="Disable csv output",
        required=False,
        action="store_true",
    )

    group.add_argument(
        "--h5_output",
        type=Path,
        help="Write output to a h5 file that contains Lx20 matrices per sequence per model; -1 means no score computed",
        default=None,
    )

    return {'vespa_pred_output': "output/"}


def main():
    parser, default_args = setup_argparse()
    embedding_group = parser.add_argument_group(title='embedding_args')
    logodds_group = parser.add_argument_group(title='logodss_args')
    conspred_group = parser.add_argument_group(title='conspred_args')
    vespa_group = parser.add_argument_group(title='vespa_args')

    default_emb_args = embedding_argparse(embedding_group)
    default_logodds_args = logodds_argparse(logodds_group)
    default_conspred_args = conspred_argparse(conspred_group)
    default_vespa_args = vespa_argparse(vespa_group)

    args = parser.parse_args()

    vespa_output_directory = args.vespa_output_directory.with_suffix('') if args.vespa_output_directory \
        else default_args['vespa_output_directory']
    vespa_output_directory.mkdir(parents=True, exist_ok=True)
    print(f'Save vespa output to {vespa_output_directory.resolve()}')

    if not args.prott5_weights_cache:
        warnings.warn("Memory Warning! Downloading the weights of ProtT5 due to the absence of a previous download!")
    cache_dir = args.prott5_weights_cache or vespa_output_directory.joinpath(default_args['prott5_weights_cache'])

    embedding_arguments = {
        'seq_path': args.fasta_input,
        'out_path':
            args.embedding_output.with_suffix('.h5') if args.embedding_output
            else vespa_output_directory.joinpath(default_emb_args['embedding_output']),
        'cache_dir': cache_dir
    }

    conspred_arguments = {
        'seq_path':
            args.use_existing_embeddings.with_suffix('.h5') if args.use_existing_embeddings
            else embedding_arguments['out_path'],
        'checkpoint_path':
            args.checkpoint or Path(vespa_location.joinpath(MODEL_PATH_DICT["CONSCNN"])),
        'out_prefix':
            args.conspred_output_prefix.with_suffix('') if args.conspred_output_prefix
            else vespa_output_directory.joinpath(default_conspred_args['conspred_output_prefix']),
        'write_probs': args.output_probs,
        'write_classes': args.output_classes
    }

    logodds_arguments = {
        'seq_path': args.fasta_input,
        'out_path':
            args.logodds_output.with_suffix('.h5') if args.logodds_output
            else vespa_output_directory.joinpath(default_logodds_args['logodds_output']),
        'mutation_file': args.mutation_file,
        'one_based_mutations': args.one_based_mutations,
        'reconstruction_probas': args.reconstruction_probas,
        'reconstruction_output':
            args.reconstruction_output.with_suffix('.h5') if args.reconstruction_output
            else vespa_output_directory.joinpath(default_logodds_args['reconstruction_output']),
        'single_csv': args.single_csv,
        'csv_dir': args.csv_dir,
        'cache_dir': cache_dir
    }

    vespa_arguments = {
        'seq_path': args.fasta_input,
        'out_path':
            args.vespa_pred_output.with_suffix('') if args.vespa_pred_output
            else vespa_output_directory.joinpath(default_vespa_args['vespa_pred_output']),
        'vespa': args.vespa,
        'vespal': args.vespal,
        'cons_input': Path(str(conspred_arguments['out_prefix']) + "_probs.h5"),
        'mutation_file': args.mutation_file,
        'one_based_mutations': args.one_based_mutations,
        'no_csv': args.no_csv,
        'h5_output': args.h5_output,
        'T5_input': logodds_arguments['out_path'],
    }

    if not vespa_arguments['vespa'] and not vespa_arguments['vespal']:
        raise RuntimeError(
            "Requested to run no vespa-model. Please at least allow one model! E.g. `--vespa` or `--vespal`"
        )

    if not args.use_existing_embeddings:
        warnings.warn('Generating embeddings might not be possible due to a GPU without enough memory!')
        run_embedding(**embedding_arguments)
    run_conspred(**conspred_arguments)
    if vespa_arguments['vespa']:
        run_logodds(**logodds_arguments)
    run_vespa(**vespa_arguments)


if __name__ == '__main__':
    main()
