#!/usr/bin/python3

"""
Author: Saul Pierotti
Mail: saulpierotti.bioinfo@gmail.com
Last updated: 22/03/2021

This script takes as input:
- a joblib dump of a dictionary containing a distance matrix for a pdb file and
  other data (output of ./pdb_to_distance_matrix.py)
- a csv file containing sifts mappings from pdb to uniprot

It expects in the current directory:
- a fasta file containing a single uniprot sequence (must be named
  uniprot_id.fasta, where uniprot_id is the id that can be seen in the sifts
  mapping of pdb_id,chain_id)

It produces in output:
- a joblib compressed dump of a numpy matrix containing the pdb distance matrix
  mapped to the uniprot sequence. Missing positions are occupied by numpy.nan

The script can also operate in batch mode using the --id_list flag. In this
case it expects in input:
- a csv file containing pdb_id,chain_id for each of the inputs to process
- a csv file containing sifts mappings from pdb to uniprot

The output is the same (but one output for each input in the csv file).
It expects to find in the current folder the relevant files as for single
inputs.
"""

import argparse
import os

import joblib
import numpy as np
import pandas as pd
from Bio import SeqIO


def map_inner_routine(d_mat_raw, indexes_zip, uniprot_seq, pdb_seq):
    """
    Map the distance matrix d_mat_raw to the sequence uniprot_seq using the
    indexes tuples indexes_zip for marking the represented regions.
    Returns a numpy array containing the mapped distance matrix where
    the indeces correspond to uniprot positions - 1.
    """
    d_mat_mapped = np.empty((len(uniprot_seq), len(uniprot_seq)))
    d_mat_mapped[:, :] = np.nan
    raw_b_row = 0

    for map_b_row, map_e_row in indexes_zip:
        raw_e_row = raw_b_row + (map_e_row - map_b_row)
        raw_b_col = 0

        for map_b_col, map_e_col in indexes_zip:
            raw_e_col = raw_b_col + (map_e_col - map_b_col)
            d_mat_mapped[map_b_row:map_e_row, map_b_col:map_e_col] = d_mat_raw[
                raw_b_row:raw_e_row, raw_b_col:raw_e_col
            ]
            assert (
                uniprot_seq[map_b_col:map_e_col]
                == pdb_seq[raw_b_col:raw_e_col]
            )
            raw_b_col = raw_e_col
        assert uniprot_seq[map_b_row:map_e_row] == pdb_seq[raw_b_row:raw_e_row]
        raw_b_row = raw_e_row

    return d_mat_mapped


def map_distance_matrix(d_mat_file, sifts_df):
    """
    Map a distance matrix contained in d_mat_file using the mapping reported in
    sifts_df. The actual mapping is performed by calling map_inner_routine().
    """
    out = dict()
    d_mat_dict = joblib.load(d_mat_file)
    uniprot_id_set = set(
        sifts_df[sifts_df.PDB == d_mat_dict["pdb_id"]].SP_PRIMARY
    )
    assert len(uniprot_id_set) == 1
    uniprot_id = uniprot_id_set.pop()
    uniprot_seq = next(
        SeqIO.parse(
            uniprot_id + ".fasta",
            "fasta",
        )
    ).seq
    curr_df = sifts_df[
        (sifts_df.PDB == d_mat_dict["pdb_id"])
        & (sifts_df.CHAIN == d_mat_dict["chain_id"])
    ]
    begin_indexes = list(curr_df.SP_BEG - 1)
    end_indexes = list(curr_df.SP_END)
    _ = begin_indexes.sort(), end_indexes.sort()
    out["distance_matrix"] = map_inner_routine(
        d_mat_dict["distance_matrix"],
        # need to call list since zip can be used only for 1 iteration
        list(zip(begin_indexes, end_indexes)),
        uniprot_seq,
        d_mat_dict["sequence"],
    )
    out["uniprot_seq"] = np.array(uniprot_seq)

    return out


def save_out_mat(d_mat_file, sifts_df):
    """
    Call the function that maps the distance matrix and save the output to a
    joblib dump. This function manages the handling of filenames in input and
    output for a single input.
    """
    assert os.path.isfile(d_mat_file)
    assert d_mat_file.endswith(".distance_matrix_dict.joblib.xz")
    d_mat_file_aslist = d_mat_file.split(".")
    outfile = ".".join(
        d_mat_file_aslist[:-3]
        + ["uniprot_distance_matrix"]
        + d_mat_file_aslist[-2:]
    )
    out = map_distance_matrix(d_mat_file, sifts_df)
    joblib.dump(out, outfile)


def main(args):
    """
    Load the sifts mapping dataframe, determine if --id_list has been called
    and call the function that dumps the output to the output file as requested
    on a single file or on all the files on the input list.
    """
    assert os.path.isfile(args.sifts_mapping)
    assert args.sifts_mapping.endswith(".csv")
    sifts_df = pd.read_csv(
        args.sifts_mapping,
        names=[
            "PDB",
            "CHAIN",
            "SP_PRIMARY",
            "RES_BEG",
            "RES_END",
            "PDB_BEG",
            "PDB_END",
            "SP_BEG",
            "SP_END",
        ],
    )

    if args.id_list:
        assert os.path.isfile(args.input_file)
        assert args.input_file.endswith(".csv")

        with open(args.input_file) as handle:
            inputs = [line.rstrip().split(",") for line in handle]

        for pdb_id, chain_id in inputs:
            d_mat_file = (
                pdb_id + "_" + chain_id + ".distance_matrix_dict.joblib.xz"
            )
            save_out_mat(d_mat_file, sifts_df)
    else:
        save_out_mat(args.input_file, sifts_df)


def parse_arguments():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="""
This script takes as input:
- a joblib dump of a dictionary containing a distance matrix for a pdb file and
  other data (output of ./pdb_to_distance_matrix.py)
- a csv file containing sifts mappings from pdb to uniprot

It expects in the current directory:
- a fasta file containing a single uniprot sequence (must be named
  uniprot_id.fasta, where uniprot_id is the id that can be seen in the sifts
  mapping of pdb_id,chain_id)

It produces in output:
- a joblib compressed dump of a numpy matrix containing the pdb distance matrix
  mapped to the uniprot sequence. Missing positions are occupied by numpy.nan

The script can also operate in batch mode using the --id_list flag. In this
case it expects in input:
- a csv file containing pdb_id,chain_id for each of the inputs to process
- a csv file containing sifts mappings from pdb to uniprot

The output is the same (but one output for each input in the csv file).
It expects to find in the current folder the relevant files as for single
inputs.
            """
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="""
        the output of ./pdb_to_distance_matrix.py for pdb_id. If --id_list is
        specified, then this should instead be a csv file containing
        pdb_id,chain_id for each of the inputs to process
        """,
    )
    parser.add_argument(
        "sifts_mapping",
        type=str,
        help="""
        a csv file containing the uniprot to pdb mappings. There should
        be no header lines
        """,
    )
    parser.add_argument(
        "--id_list",
        help="""
            interpret input_file as a csv file containing a list of chains and
            pdb ids instead than as a mmcif file
        """,
        action="store_true",
    )

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    ARGS = parse_arguments()
    main(ARGS)
