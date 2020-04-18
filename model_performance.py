#!/usr/bin/python3

"""
Author: Saul Pierotti
Mail: saulpierotti.bioinfo@gmail.com
Last updated: 17/04/2020

Takes as first argument a space-separated file with 3 columns containing:

ID feature class

The input file must contain only unique IDs.
Takes as second argument specifying the desired threshold feature value used
for the classification.
If called from the shell it prints:

confusion_matrix ACC MCC

If called as a library it provides the get_stats() function, which returns the
same information in a tuple.
"""

import sys
import math


def parse_input(filepath):
    with open(filepath) as filein:
        df = []
        for line in filein:
            ID, score, bin_class = line.rstrip().split()
            df.append((ID, float(score), int(bin_class)))
    return df


def get_confusion_mat(df, thr):
    true_pos, true_neg, false_pos, false_neg = 0, 0, 0, 0
    for row in df:
        if row[1] < thr and row[2] == 1:
            true_pos += 1
        elif row[1] < thr and row[2] == 0:
            false_pos += 1
        elif row[1] >= thr and row[2] == 1:
            false_neg += 1
        elif row[1] >= thr and row[2] == 0:
            true_neg += 1
    confusion_mat = ((true_pos, false_pos), (true_neg, false_neg))
    return confusion_mat


def get_ACC(df, thr):
    confusion_mat = get_confusion_mat(df, thr)
    t_pos = confusion_mat[0][0]
    t_neg = confusion_mat[1][0]
    f_pos = confusion_mat[0][1]
    f_neg = confusion_mat[1][1]
    above_frac = t_pos + t_neg
    below_frac = t_pos + t_neg + f_pos + f_neg
    ACC = above_frac / below_frac
    return ACC


def get_MCC(df, thr):
    confusion_mat = get_confusion_mat(df, thr)
    t_pos = confusion_mat[0][0]
    t_neg = confusion_mat[1][0]
    f_pos = confusion_mat[0][1]
    f_neg = confusion_mat[1][1]
    above_frac = t_pos * t_neg - f_pos * f_neg
    below_frac = math.sqrt(
        (t_pos + f_pos) * (t_pos + f_neg) * (t_neg + f_pos) * (t_neg + f_neg)
    )
    MCC = above_frac / max(below_frac, 1)  # otherwise I can divide by 0
    return MCC


def get_stats(filepath, thr, print_res=False):
    df = parse_input(filepath)
    confusion_mat = get_confusion_mat(df, thr)
    ACC = get_ACC(df, thr)
    MCC = get_MCC(df, thr)
    if print_res:
        print(confusion_mat, ACC, MCC)
    else:
        return confusion_mat, ACC, MCC


if __name__ == "__main__":
    print_res = True
    input_file = sys.argv[1]
    threshold_score = float(sys.argv[2])
    get_stats(input_file, threshold_score, print_res)
