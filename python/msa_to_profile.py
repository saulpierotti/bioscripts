#!/usr/bin/python3

"""
Author: Saul Pierotti
Mail: saulpierotti.bioinfo@gmail.com
Last updated: 24/03/2021

This script takes in input a multiple sequence alignment in fasta
format with gaps represented as -. It produces a sequence profile
in output formatted as a numpy array and places it in a joblib
dump. The extension of the input is expected to be .aln
"""

import argparse
import os

import joblib
import numpy as np
from Bio import AlignIO
from Bio.Data import IUPACData
from sklearn.preprocessing import OrdinalEncoder


def get_profile(align_vec, pseudocount=0):
    """
    Takes in input a multiple sequence alignment in a ordinally encoded numpy
    array and produces a sequence profile in output, optionally with a
    pseudocount.
    """
    profile = []

    for column in align_vec.T:
        # add 1 instance of every possible category so that all the counts have the same lenght
        column = np.concatenate([column, range(len(SYMBOLS))])
        _, counts = np.unique(column, return_counts=True)
        # remove the artificial count and add the pseudocount
        counts = counts - 1 + pseudocount
        frequencies = counts / counts.sum()
        profile.append(frequencies)

    return np.array(profile).T


def get_align_vec(align):
    """
    Return a ordinally-encoded vector representing the multiple sequence
    alignment given in input
    """
    categories = np.array(
        [np.array(list(SYMBOLS)) for _ in range(align.get_alignment_length())]
    ).T
    align_vec = OrdinalEncoder().fit(categories).transform(np.array(align))

    return align_vec


def main(args):
    """
    Main function
    """

    assert os.path.isfile(args.i)
    assert args.i.endswith(".aln")
    assert args.o.endswith("profile.joblib.xz")
    assert not os.path.isfile(args.o)

    align = AlignIO.read(args.i, "fasta")
    align_vec = get_align_vec(align)
    profile = get_profile(align_vec)
    joblib.dump(profile, args.o)


def parse_arguments():
    """
    Parse command line arguments.
    """
    description = " ".join(__doc__.splitlines()[4:])

    epilog = ", ".join(__doc__.splitlines()[1:4])
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-i",
        type=str,
        help="""
        the multiple sequence alignment fasta file with extension .aln to
        use in input
        """,
        metavar="<file>",
        required=True,
    )
    parser.add_argument(
        "-o",
        metavar="<file>",
        type=str,
        help="the file where to write the profile numpy array",
        required=True,
    )

    args = parser.parse_args()

    return args


SYMBOLS = IUPACData.extended_protein_letters + "-"

if __name__ == "__main__":
    ARGS = parse_arguments()
    main(ARGS)
