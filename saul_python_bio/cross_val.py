#!/usr/bin/env python

"""
Author: Saul Pierotti
Mail: saulpierotti.bioinfo@gmail.com
Last updated: 24/04/2020

Takes a dataset and an optional tolerance as arguments
Produces a report file containing the confusion matrix, ACC, MCC and selected
treshold for 10 randomised cross validation runs on a 80/20 split of the
dataset.
It produces also a roc curve plot on the entire dataset and a plot of the MCC vs
threshold used for all the training sets.
It prints a log of the operation to STDOUT.
"""

import sys
import numpy as np
import pandas as pd
from sklearn import metrics
from sklearn import utils
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()


def get_df(path):
    print("Inizialising dataset:", path)
    col_names = ["ID", "Feature", "Class"]
    df = pd.read_csv(path, sep="\s+", names=col_names)
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
        if feature[i] > thr and data_class[i] == 1:
            true_pos += 1
        elif feature[i] > thr and data_class[i] == 0:
            false_pos += 1
        elif feature[i] <= thr and data_class[i] == 1:
            false_neg += 1
        elif feature[i] <= thr and data_class[i] == 0:
            true_neg += 1
    confusion_mat = ((true_pos, false_pos), (true_neg, false_neg))
    return confusion_mat


def get_confusion_mat_skl(df, thr):
    y_pred = [1 if (value > thr) else 0 for value in df["Feature"].values]
    y_true = df["Class"].values
    confusion_mat = metrics.confusion_matrix(y_true, y_pred)
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


def get_ACC_skl(df, thr):
    y_pred = [1 if (value > thr) else 0 for value in df["Feature"].values]
    y_true = df["Class"].values
    ACC = metrics.accuracy_score(y_true, y_pred)
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


def get_MCC_skl(df, thr):
    y_pred = [1 if (value > thr) else 0 for value in df["Feature"].values]
    y_true = df["Class"].values
    MCC = metrics.matthews_corrcoef(y_true, y_pred)
    return MCC


def get_stats(df, thr):
    confusion_mat = __get_confusion_mat(df, thr)
    ACC = __get_ACC(df, thr)
    MCC = __get_MCC(df, thr)
    return confusion_mat, ACC, MCC


def thr_explore(df):
    start = 20
    stop = 30
    MCC_list, thr_list = [], []
    for thr in range(start, stop):
        MCC = __get_MCC(df, thr)
        MCC_list.append(MCC)
        thr_list.append(thr)
        print(thr, MCC)
    return MCC_list, thr_list


def get_best_thr(df):
    MCC, thr = thr_explore(df)
    best_MCC = -2
    for i in range(len(MCC)):
        if MCC[i] > best_MCC:
            best_thr = [thr[i]]
            best_MCC = MCC[i]
        elif MCC[i] == best_MCC:
            best_thr.append(thr[i])
    thr_out = sum(best_thr) / len(best_thr)
    print("Selected threshold:", thr_out)
    return thr_out, thr, MCC


def train_and_test(df, test_prop=0.2, num_replicates=10):
    train_report, all_MCCs = [], []
    for i in range(num_replicates):
        print("Training replicate", i + 1)
        df = utils.shuffle(df)
        df_train, df_cv = get_train_test(df, test_prop)
        thr, thr_list, MCC_list = get_best_thr(df_train)
        test_results = get_stats(df_cv, thr)
        print("Performance on the test set", list(test_results))
        train_report.append([test_results, thr])
        all_MCCs += MCC_list
    thr_MCC_report = [thr_list * num_replicates, all_MCCs]
    return train_report, thr_MCC_report


def write_report(train_report_list):
    sep = "\t"
    argv_index = 1
    with open("cross_val_report.dat", "w") as f:
        for train_report in train_report_list:
            f.write("Input file: " + str(sys.argv[argv_index]) + "\n")
            argv_index += 1
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
                + "\n---\n"
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
                    + str(thr)
                    + "\n"
                )
            list_MCC = [train[0][2] for train in train_report]
            avg_MCC = sum(list_MCC) / len(train_report)
            std_MCC = np.std(np.array(list_MCC))
            avg_bitscore = sum([train[1] for train in train_report]) / len(train_report)
            list_cm = [np.array(train[0][0]) for train in train_report]
            avg_cm = [int(val) for row in np.mean(list_cm, axis=0) for val in row]
            f.write(
                "---"
                + "\nAverage MCC: "
                + str(avg_MCC)
                + "\nStandard deviation MCC: "
                + str(std_MCC)
                + "\nAverage bit score threshold: "
                + str(avg_bitscore)
                + "\nAverage confusion matrix: "
                + str(avg_cm)
            )
            f.write("\n### end of dataset\n")


def plot_roc_curve(datasets):
    for dataset in datasets:
        y_true = dataset["Class"].values
        y_score = dataset["Feature"].values
        fpr, tpr, thresholds = metrics.roc_curve(y_true, y_score)
        roc_plt = sns.lineplot(x=fpr, y=tpr, markers=True)
    plt.legend(labels=sys.argv[1:])
    plt.ylabel("True Positive Rate")
    plt.xlabel("False Positive Rate")
    plt.savefig("roc_plot.png")
    plt.clf()
    print("AUC for the ROC curve in the complete dataset:", metrics.auc(fpr, tpr))


def plot_thr_MCC(thr_MCC_report_list):
    for thr_MCC_report in thr_MCC_report_list:
        sns.lineplot(x=thr_MCC_report[0], y=thr_MCC_report[1], markers=True)
    plt.ylabel("MCC")
    plt.xlabel("bit score threshold")
    plt.legend(labels=sys.argv[1:])
    plt.savefig("bitscore_MCC_plot.png")
    plt.clf()


if __name__ == "__main__":
    train_report_list, thr_MCC_report_list = [], []
    dfs = []
    for path in sys.argv[1:]:
        df = get_df(path)
        train_report, thr_MCC_report = train_and_test(df, 0.2, 1)
        train_report_list.append(train_report)
        thr_MCC_report_list.append(thr_MCC_report)
        dfs.append(df)
    plot_roc_curve(dfs)
    plot_thr_MCC(thr_MCC_report_list)
    write_report(train_report_list)
