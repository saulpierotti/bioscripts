#!/usr/bin/python3

"""
Author: Saul Pierotti
Mail: saulpierotti.bioinfo@gmail.com
Last updated: 22/03/2021

Given the output of the script ./pdb_to_distance_matrix.py, a sifts mapping,
and a uniprot sequence (read implicitly), produces a mapped version of the
distance matrix and saves it in a joblib dump. The dump contains a dictionary
with the distance matrix itself and the uniprot sequence. It expects a file
named uniprot_id.seq for each uniprot_id in the sifts mapping.
"""

import argparse
import os

import joblib
import numpy as np
import pandas as pd
from Bio import SeqIO


def map_d_mat(d_mat_raw, indexes_zip, uniprot_seq, pdb_seq):
    """
    Map the distance matrix d_mat_raw to the sequence uniprot_seq using the
    indexes tuples indexes_zip for marking the represented regions.
    Returns a numpy array containing the mapped distance matrix where
    the indexes correspond to uniprot positions - 1.
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


def get_uniprot_seq(sifts_df, d_mat_dict):
    """
    Reads from sifts_df which uniprot_id corresponds to the pdb chain in
    d_mat_dict and reads it from the current directory. It expects a file with
    the name uniprot_id.fasta in the current directory.
    """
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

    return uniprot_seq


def get_sifts_df(db_file):
    """
    Wrapper for pandas.read_csv with custom parameters
    """
    sifts_df = pd.read_csv(
        db_file,
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

    return sifts_df


def get_out_dict(d_mat, sifts_df, uniprot_seq):
    """
    Map a distance matrix contained in d_mat_file using the mapping reported in
    sifts_df. It expects the uniprot sequences to be in fasta files called .
    """
    out = dict()
    curr_df = sifts_df[
        (sifts_df.PDB == d_mat["pdb_id"])
        & (sifts_df.CHAIN == d_mat["chain_id"])
    ]
    begin_indexes = list(curr_df.SP_BEG - 1)
    end_indexes = list(curr_df.SP_END)
    _ = begin_indexes.sort(), end_indexes.sort()
    out["distance_matrix"] = map_d_mat(
        d_mat["distance_matrix"],
        # need to call list since zip can be used only for 1 iteration
        list(zip(begin_indexes, end_indexes)),
        uniprot_seq,
        d_mat["sequence"],
    )
    out["uniprot_seq"] = np.array(uniprot_seq)

    return out


def main(args):
    """
    Main function
    """
    assert os.path.isfile(args.db)
    assert args.db.endswith(".csv")
    assert os.path.isfile(args.i)
    assert args.i.endswith("pdb_distance_matrix.joblib.xz")
    assert args.o.endswith("uniprot_distance_matrix.joblib.xz")
    assert not os.path.isfile(args.o)
    sifts_df = get_sifts_df(args.db)
    d_mat = joblib.load(args.i)
    uniprot_seq = get_uniprot_seq(sifts_df, d_mat)
    out = get_out_dict(d_mat, sifts_df, uniprot_seq)
    joblib.dump(out, args.o)


def parse_arguments():
    """
    Parse command line arguments.
    """
    description = " ".join(__doc__.splitlines()[4:])

    epilog = ", ".join(__doc__.splitlines()[1:4])
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-i",
        metavar="<file>",
        type=str,
        help="the output of the script ./pdb_to_distance_matrix.py",
        required=True,
    )
    parser.add_argument(
        "-o",
        metavar="<file>",
        type=str,
        help="""the file where to save the joblib dump of the uniprot-mapped
        distance matrix""",
        required=True,
    )
    parser.add_argument(
        "-db",
        metavar="<file>",
        type=str,
        help="""
        a csv file containing the uniprot to pdb mappings. There should
        be no header lines. It can be obtained from
        ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/uniprot_segments_observed.csv.gz
        """,
        required=True,
    )

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    ARGS = parse_arguments()
    main(ARGS)
