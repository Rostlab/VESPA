import argparse
from vespa.predict.embedding import *
from vespa.predict.config import VERBOSE, CACHE_DIR


def setup_argparse():
    parser = argparse.ArgumentParser(
        description="embedding.py creates T5 embeddings for a given text "
                    + " file containing sequence(s) in FASTA-format."
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
        default=Path("T5_embedding.h5"),
        help="Output path to store the computed embeddings in h5 format. "
             "Does not contain the original sequence data, only the sequence id!",
    )

    parser.add_argument(
        "--prott5_weights_cache",
        required=False,
        type=Path,
        default=Path(CACHE_DIR),
        help="Path to directory for storing the weights of the language model ProtT5.",
    )

    return parser


def run_embedding(seq_path, out_path, cache_dir):
    if VERBOSE:
        print(f' Start Embedding of {seq_path.stem} '.center(80, '#'))
    t5_emb = T5_Embed(cache_dir=cache_dir)
    t5_emb.embed_from_fasta(fasta_path=seq_path, output_path=out_path)
    if VERBOSE:
        print(f'>> Finished Embeddings!')


def main():
    parser = setup_argparse()
    args = parser.parse_args()

    arguments = {
        'seq_path': args.input,
        'out_path': args.output,
        'cache_dir': args.prott5_weights_cache
    }
    run_embedding(**arguments)


if __name__ == '__main__':
    main()
