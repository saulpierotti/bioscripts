#!/usr/bin/python3

"""
Author: Saul Pierotti
Mail: saulpierotti.bioinfo@gmail.com
Last updated: 19/03/2021
"""

import argparse

import joblib
import numpy as np
from Bio import PDB, SeqUtils
from scipy import spatial


def get_distance_matrix(mmcif_file, chain_id):
    """
    Given a protein structure in mmcif format and a chain id, extract the
    residue type, the coordinate of each residue, and the resseq id. Compute
    all the pairwise euclidean distances among residues. Returns a dictionary
    containing all of these data.
    """
    parser = PDB.MMCIFParser()
    structure = parser.get_structure("_", mmcif_file)
    out = {"residue": [], "coordinates": [], "resseq": []}

    for chain in structure.get_chains():
        if chain.id != chain_id:
            continue

        for residue in chain.get_residues():
            het_field = residue.id[0]

            if het_field != " ":
                continue
            out["residue"].append(residue.resname)
            out["coordinates"].append(residue["CA"].get_coord())
            out["resseq"].append(residue.id[1])

    out["coordinates"] = np.array(out["coordinates"])
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

    return out


def parse_arguments():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="""
            This script takes as input a file containing a sequence profile in
            as a saved numpy array with .profile.npy extension or a text file
            with a list of filenames for sequence profiles. It uses an svm
            model, provided to the script as a joblib dump, for predicting the
            secondary structure of the protein sequence(s) provided. It
            produces in output a fasta-like file with extension .pred.fa
            containing the predicted conformations.
        """
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="""
        a fasta file with a single sequence (default) or a text file containing
        a list of fasta filenames (with --id_list)
        """,
    )
    parser.add_argument(
        "svm_model",
        type=str,
        help="a trained svm model saved as a joblib dump",
    )
    parser.add_argument(
        "--id_list",
        help="""
            interpret input_file as a text file containing a list of fasta
            files (one per row) instead that considering it a fasta file
        """,
        action="store_true",
    )

    arguments = parser.parse_args()

    return arguments


if __name__ == "__main__":
    args = parse_arguments()
    main(args.input_file, args.svm_model, args.id_list)
