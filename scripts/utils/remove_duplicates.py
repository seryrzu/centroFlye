from argparse import ArgumentParser
import os
import sys

this_dirname = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(this_dirname, os.path.pardir))

from bio import read_bio_seqs, write_bio_seqs


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    params = parser.parse_args()
    seqs = read_bio_seqs(params.input)
    write_bio_seqs(params.output, seqs)
