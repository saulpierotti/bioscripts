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

If called from the shell it prints a series of statistical scores.
If called as a library it provides the get_stats() function, which returns the
same information in a tuple.
"""

import sys
import math


def parse_input(filepath):
    with open(filepath) as filein:
        d = {}
        for line in filein:
            ID, score, bin_class = line.rstrip().split()
            d[ID] = d.get(ID, [])
            d[ID].append([float(score), int(bin_class)])
        d.values()
        print(df)
    return df


def get_confusion_mat(data, thr):
    # I evaluate only below thresholds since only those are meaningful for blast (above I can have as many as I want!)
    true_pos_arr = p_vals[
        p_vals < thr
    ]  # np.arr<value returns an np.array with [True, False, True, ...]
    false_pos_arr = n_vals[
        n_vals < thr
    ]  # if I use it as an index for the array itself it filters only the values meeting the condition
    # the numerosity of a class is just the lenght of its array
    true_pos = len(true_pos_arr)
    false_pos = len(false_pos_arr)
    # I calculate the remaining class from the differences
    tot_p = len(p_vals)
    tot_n = len(n_vals)
    false_neg = tot_p - true_pos
    true_neg = tot_n - false_pos
    confusion_mat = np.array([[true_pos, false_pos], [true_neg, false_neg]])
    return confusion_mat


def get_ACC(thr, p_vals, n_vals):
    confusion_mat = get_confusion_mat(thr, p_vals, n_vals)
    t_pos = confusion_mat[0][0]
    t_neg = confusion_mat[1][0]
    f_pos = confusion_mat[0][1]
    f_neg = confusion_mat[1][1]
    above_frac = t_pos + t_neg
    below_frac = t_pos + t_neg + f_pos + f_neg
    ACC = above_frac / below_frac
    return ACC


def get_MCC(thr, p_vals, n_vals):
    confusion_mat = get_confusion_mat(thr, p_vals, n_vals)
    t_pos = confusion_mat[0][0]
    t_neg = confusion_mat[1][0]
    f_pos = confusion_mat[0][1]
    f_neg = confusion_mat[1][1]
    above_frac = t_pos * t_neg - f_pos * f_neg
    below_frac = np.sqrt(
        (t_pos + f_pos) * (t_pos + f_neg) * (t_neg + f_pos) * (t_neg + f_neg)
    )
    MCC = above_frac / max(below_frac, 1)  # otherwise I can divide by 0
    return MCC


def get_stats(thr, p_vals, n_vals):
    confusion_mat = get_confusion_mat(thr, p_vals, n_vals)
    ACC = get_ACC(thr, p_vals, n_vals)
    MCC = get_MCC(thr, p_vals, n_vals)
    print("Confusion matrix:\n")
    print(confusion_mat, "\n")
    print("Threshold:", thr)
    print("ACC:", ACC)
    print("MCC:", MCC)
    return confusion_mat, ACC, MCC


if __name__ == "__main__":
    print("model_performance.py called")
    input_file = sys.argv[1]
    threshold_score = sys.argv[2]
    df = parse_input(input_file)
