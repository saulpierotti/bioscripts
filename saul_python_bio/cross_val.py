#!/usr/bin/env python

"""
Author: Saul Pierotti
Mail: saulpierotti.bioinfo@gmail.com
Last updated: 23/04/2020

Takes a dataset and an optional tolerance as arguments
Produces a report file containing the confusion matrix, ACC, MCC and selected
treshold for 10 randomised cross validation runs on a 80/20 split of the
dataset.
It prints a log of the operation to STDOUT.
"""

import numpy as np
import pandas as pd
from sklearn import metrics
from sklearn import utils
import matplotlib.pyplot as plt
import sys
import seaborn as sns

sns.set()


def get_df(path):
    col_names = ["ID", "Feature", "Class"]
    df = pd.read_csv(dataset_path, sep="\s+", names=col_names)
    return df


def get_train_test(df, test_prop=0.2):
    test_size = int(test_prop * len(df))
    df_test, df_train = np.split(df, [test_size])
    return df_train, df_test


def __get_confusion_mat(df, thr):
    feature = df["Feature"].values
    data_class = df["Class"].values
    true_pos, false_pos, true_neg, false_neg = 0, 0, 0, 0
    for i in range(len(feature)):
        if feature[i] < thr and data_class[i] == 1:
            true_pos += 1
        elif feature[i] < thr and data_class[i] == 0:
            false_pos += 1
        elif feature[i] >= thr and data_class[i] == 1:
            false_neg += 1
        elif feature[i] >= thr and data_class[i] == 0:
            true_neg += 1
    confusion_mat = ((true_pos, false_pos), (true_neg, false_neg))
    return confusion_mat


def __get_ACC(df, thr):
    confusion_mat = __get_confusion_mat(df, thr)
    t_pos = confusion_mat[0][0]
    t_neg = confusion_mat[1][0]
    f_pos = confusion_mat[0][1]
    f_neg = confusion_mat[1][1]
    above_frac = t_pos + t_neg
    below_frac = t_pos + t_neg + f_pos + f_neg
    ACC = above_frac / below_frac
    return ACC


def __get_MCC(df, thr):
    confusion_mat = __get_confusion_mat(df, thr)
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


def get_stats(df, thr):
    confusion_mat = __get_confusion_mat(df, thr)
    ACC = __get_ACC(df, thr)
    MCC = __get_MCC(df, thr)
    return confusion_mat, ACC, MCC


def optimize_MCC(
    df, start, stop, tolerance,
):
    print("Optimisation in range", start, stop)
    MCC_list, thr_list = [], []
    for thr in np.logspace(start, stop, 10):
        MCC_list.append(__get_MCC(df, thr))
        thr_list.append(thr)
    best_MCC = max(MCC_list)
    avg_MCC = sum(MCC_list) / len(MCC_list)
    max_i = MCC_list.index(best_MCC)
    diff = abs(avg_MCC - best_MCC)
    print("Tolerance set to", tolerance)
    print("Current difference in MCC", diff)
    if diff > tolerance:
        new_start = np.log10(thr_list[max_i - 1])
        new_stop = np.log10(thr_list[max_i + 1])
        return optimize_MCC(df, new_start, new_stop, tolerance)
    else:
        best_thr = thr_list[max_i]
        return best_thr


def train_and_test(df, tolerance, test_prop=0.2, num_replicates=10):
    train_report = []
    for i in range(num_replicates):
        print("Training replicate", i + 1)
        df = utils.shuffle(df)
        df_train, df_cv = get_train_test(df, test_prop)
        start = np.log10(min(df["Feature"].values)) - 1
        stop = np.log10(max(df["Feature"].values)) + 1
        thr = optimize_MCC(df_train, start, stop, tolerance)
        train_report.append([get_stats(df_cv, thr), thr])
    return train_report


def write_report(train_report):
    sep = "\t"
    with open("cross_val_report.dat", "w") as f:
        f.write(
            "tp"
            + sep
            + "fp"
            + sep
            + "tn"
            + sep
            + "fn"
            + sep
            + "ACC"
            + sep
            + "MCC"
            + sep
            + "thr"
            + "\n"
        )
        for train in train_report:
            stats, thr = train[0], train[1]
            cm, ACC, MCC = stats[0], stats[1], stats[2]
            f.write(
                str(cm[0][0])
                + sep
                + str(cm[0][1])
                + sep
                + str(cm[1][0])
                + sep
                + str(cm[1][1])
                + sep
                + str(ACC)
                + sep
                + str(MCC)
                + sep
                + "{:e}".format(thr)
                + "\n"
            )


def plot_roc_curve(test_set):
    y_true = test_set["Class"].values
    y_score = [-val for val in test_set["Feature"].values]
    fpr, tpr, thresholds = metrics.roc_curve(y_true, y_score, pos_label=1)
    roc_plt = sns.lineplot(x=fpr, y=tpr)
    plt.ylabel("True Positive Rate")
    plt.xlabel("False Positive Rate")
    plt.savefig("roc_plot.png")
    print("AUC for the ROC curve in the complete dataset:", metrics.auc(fpr, tpr))


if __name__ == "__main__":
    dataset_path = sys.argv[1]
    tolerance = float(sys.argv[2])
    df = get_df(dataset_path)
    plot_roc_curve(df)
    train_report = train_and_test(df, tolerance, 0.2, 1)
    write_report(train_report)
