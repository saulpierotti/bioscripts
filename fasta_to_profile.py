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


def get_profile(D_aln, all_aa="ARNDCQEGHILKMFPSTWYV"):
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


if __name__ == "__main__":
    aln_file = sys.argv[1]
    D_aln = get_aln(aln_file)
    profile = get_profile(D_aln)
    print(profile)
