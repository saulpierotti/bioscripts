#!/usr/bin/python3

"""
Author: Saul Pierotti
Mail: saulpierotti.bioinfo@gmail.com
Last updated: 19/03/2021

This script takes as input a mmcif file and a chain id.
It computes the pairwise distance matrix of the alpha carbons and
saves it in a xz compressed joblib dump.
It can optionally also take a list of pdb ids and do the same
action for all the files with such names in the current directory.
"""

import argparse
import os

import joblib
import numpy as np
from Bio import PDB, SeqUtils
from scipy import spatial


def get_distance_matrix(pdb_id, chain_id):
    """
    Given a protein structure in mmcif format and a chain id, extract the
    residue type, the coordinate of each residue, and the resseq id. Compute
    all the pairwise euclidean distances among residues. Returns a dictionary
    containing all of these data.
    """
    mmcif_file = pdb_id + ".cif"
    parser = PDB.MMCIFParser()
    structure = parser.get_structure("_", mmcif_file)
    out = {"residue": [], "coordinates": [], "resseq": []}

    matching_chains = 0

    for chain in structure.get_chains():
        if chain.id != chain_id:
            continue
        matching_chains += 1

        for residue in chain.get_residues():
            het_field = residue.id[0]

            if het_field != " ":
                continue
            out["residue"].append(residue.resname)
            out["coordinates"].append(residue["CA"].get_coord())
            out["resseq"].append(residue.id[1])
    assert matching_chains == 1

    out["coordinates"] = np.array(out["coordinates"])
    out["resseq"] = np.array(out["resseq"])
    # the Minkowski 2-norm is the euclidean distance
    out["distance_matrix"] = spatial.distance_matrix(
        out["coordinates"], out["coordinates"], p=2
    )
    out["sequence"] = "".join(
        [
            SeqUtils.IUPACData.protein_letters_3to1[r.capitalize()]
            for r in out["residue"]
        ]
    )
    out["pdb_id"] = pdb_id
    out["chain_id"] = chain_id

    return out


def save_out_dict(pdb_id, chain_id):
    """
    Takes a mmcif file and a chain id. Calculate the distance matrix and other
    information and store the result as a dictionary in a joblib dump.
    """
    mmcif_file = pdb_id + ".cif"
    assert os.path.isfile(mmcif_file)
    assert len(chain_id) == 1
    outfile = pdb_id + "_" + chain_id + ".distance_matrix_dict.joblib.xz"
    out_dict = get_distance_matrix(pdb_id, chain_id)
    joblib.dump(out_dict, outfile)


def main(args):
    """
    Main function
    """

    if args.id_list:
        assert os.path.isfile(args.input_file)
        assert args.input_file.endswith(".csv")

        with open(args.input_file) as handle:
            inputs = [line.rstrip().split(",") for line in handle]

        for pdb_id, chain_id in inputs:
            save_out_dict(pdb_id, chain_id)
    else:
        save_out_dict(args.input_file, args.chain_id)


def parse_arguments():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="""
            This script takes as input a mmcif file and a chain id.
            It computes the pairwise distance matrix of the alpha carbons and
            saves it in a xz compressed joblib dump.
            It can optionally also take a csv list of pdb ids and chain ids and
            do the same action for all the files with such names in the current
            directory.
            """
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="the mmcif file to use or a csv of pdb ids and chain ids",
    )
    parser.add_argument(
        "--chain_id",
        type=str,
        help="the chain id to consider (leave empty if using --id_list)",
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
