#!/usr/bin/python

import sys
import numpy as np


def get_aln(aln_file):
    D = {}
    F = open(aln_file)
    for row in F:
        if row[:3] != "sp|":
            continue
        curr_ID = row.split("|")[1]
        curr_seq = row.split()[1]
        D[curr_ID] = D.get(curr_ID, "") + curr_seq
    F.close()
    return D


def get_profile(D_aln, all_aa):
    n = len(list(D_aln.values())[0])
    M_aln = []
    for i in range(n):
        aln_pos = []
        for key in D_aln.keys():
            aln_pos.append(D_aln[key][i])
        M_aln.append(aln_pos)
    profile = []
    for pos in M_aln:
        v = np.zeros(len(all_aa))
        for aa in pos:
            aa_num = all_aa.find(aa)
            if aa_num > -1:
                v[aa_num] += 1
        sum_aa_pos = v.sum()
        for i in range(len(all_aa)):
            v[i] = v[i] / sum_aa_pos
        profile.append(v)
    return profile


def get_entropy(profile):
    entropy = 0
    for column in profile:
        col_entropy = 0
        for residue in column[:20]:  # I am avoiding to count the gap
            if residue == 0:
                continue
            res_entropy = -residue * np.log(residue)
            col_entropy += res_entropy
        entropy += col_entropy
    return entropy


def get_logo(profile):
    consensus = []
    cons_entropy = []
    cons_frequency = []
    for column in profile:
        residue = column[:20].max()
        residue_index = column[:20].argmax()
        res_entropy = str(
            abs(residue * np.log(residue))
        )  # abs() for avoiding -0.0 in conserved positions
        res_entropy += "0" * (19 - len(res_entropy))  # just for padding in print
        residue = str(residue) + "0" * (
            5 - len(str(residue))
        )  # just for padding in print
        consensus.append(residue_index)
        cons_frequency.append(residue)
        cons_entropy.append(res_entropy)
    return consensus, cons_frequency, cons_entropy


if __name__ == "__main__":
    aas_with_gap = "ARNDCQEGHILKMFPSTWYV-"
    aln_file = sys.argv[1]
    D_aln = get_aln(aln_file)
    profile = get_profile(D_aln, aas_with_gap)
    aln_entropy = get_entropy(profile)
    logo = get_logo(profile)
    zipped_logo = zip(range(1, len(logo[0]) + 1), logo[0], logo[1], logo[2])
    print("Alignment entropy (S):", aln_entropy, "\n")
    for my_tuple in zipped_logo:
        pos = my_tuple[0]
        res = aas_with_gap[:20][my_tuple[1]]
        freq = my_tuple[2]
        entropy = my_tuple[3]
        print(
            "Position:",
            pos,
            "\tResidue:",
            res,
            "\tFrequency:",
            freq,
            "\tEntropy:",
            entropy,
        )
