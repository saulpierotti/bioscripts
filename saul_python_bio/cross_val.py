#!/usr/bin/env python

"""
Author: Saul Pierotti
Mail: saulpierotti.bioinfo@gmail.com
Last updated: 24/04/2020

Takes a dataset and an optional tolerance as arguments
Produces a report file containing the confusion matrix, ACC, MCC and selected
threshold for 10 randomised cross validation runs on a 80/20 split of the
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
from sklearn import model_selection
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()


def get_df(path):
    print("Inizialising dataset:", path)
    col_names = ["ID", "Feature", "Class"]
    df = pd.read_csv(path, sep="\s+", names=col_names)
    return df


def get_train_test(df, set_indeces_iterator, numsplits=5):
    set_indeces = next(set_indeces_iterator)
    df_train = df.iloc[set_indeces[0]]
    df_test = df.iloc[set_indeces[1]]
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
    confusion_mat = ((true_pos, false_pos), (false_neg, true_neg))
    return confusion_mat


def get_confusion_mat_skl(df, thr):
    y_pred = [1 if (value > thr) else 0 for value in df["Feature"].values]
    y_true = df["Class"].values
    confusion_mat = metrics.confusion_matrix(y_true, y_pred)
    return confusion_mat


def __get_ACC(df, thr):
    confusion_mat = __get_confusion_mat(df, thr)
    t_pos = confusion_mat[0][0]
    t_neg = confusion_mat[1][1]
    f_pos = confusion_mat[0][1]
    f_neg = confusion_mat[1][0]
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
    t_neg = confusion_mat[1][1]
    f_pos = confusion_mat[0][1]
    f_neg = confusion_mat[1][0]
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


def thr_explore(df, start=-100, stop=1):
    MCC_list, thr_list = [], []
    for exp in range(start, stop):
        thr = 10 ** exp
        MCC = __get_MCC(df, thr)
        MCC_list.append(MCC)
        thr_list.append(thr)
        print(thr, MCC)
    return MCC_list, thr_list


def get_best_thr(df, start=-100, stop=1):
    MCC, thr = thr_explore(df, start, stop)
    best_MCC = -2
    for i in range(len(MCC)):
        if MCC[i] > best_MCC:
            best_thr = [thr[i]]
            best_MCC = MCC[i]
        elif MCC[i] == best_MCC:
            best_thr.append(thr[i])
    thr_out = 10 ** (sum(np.log10(best_thr)) / len(best_thr))
    print("Selected threshold:", thr_out)
    return thr_out, thr, MCC


def train_and_test(df, numsplits=5, num_replicates=1, start=-100, stop=1):
    train_report, all_MCCs, auc_list, tpr_list, fpr_list = [], [], [], [], []
    for i in range(num_replicates):
        print("\nTraining replicate", i + 1)
        set_indeces_iterator = model_selection.StratifiedKFold(
            n_splits=numsplits, shuffle=True, random_state=1
        ).split(df, y=df["Class"])
        for j in range(numsplits):
            print("\nCross validation fold", j + 1)
            df_train, df_cv = get_train_test(df, set_indeces_iterator, numsplits)
            thr, thr_list, MCC_list = get_best_thr(df_train, start, stop)
            confusion_mat, ACC, MCC = get_stats(df_cv, thr)
            train_report.append((confusion_mat, ACC, MCC, thr))
            print("\nPerformance on the test set:")
            print("ACC", str(ACC))
            print("MCC", str(MCC))
            print("Confusion matrix", confusion_mat)
            all_MCCs += MCC_list
            y_true = df_cv["Class"].values
            y_score = [-value for value in df_cv["Feature"].values]
            fpr, tpr, thresholds = metrics.roc_curve(y_true, y_score)
            curr_auc = metrics.auc(fpr, tpr)
            print("AUC of the test set", curr_auc)
            auc_list.append(curr_auc)
            tpr_list += list(tpr)
            fpr_list += list(fpr)
    thr_MCC_report = (thr_list * num_replicates * numsplits, all_MCCs)
    roc_curve_report = (auc_list, tpr_list, fpr_list)
    return train_report, thr_MCC_report, roc_curve_report


def write_report(train_report_list):
    sep = "\t"
    argv_index = 1
    with open("cross_val_report.dat", "w") as f:
        for train_report in train_report_list:
            f.write("# Input file: " + str(sys.argv[argv_index]) + "\n")
            argv_index += 1
            f.write(
                "tp"
                + sep
                + "fp"
                + sep
                + "fn"
                + sep
                + "tn"
                + sep
                + "ACC"
                + sep
                + "MCC"
                + sep
                + "thr"
                + "\n"
            )
            for train in train_report:
                cm, ACC, MCC, thr = train[0], train[1], train[2], train[3]
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
            list_MCC = [train[2] for train in train_report]
            arr_MCC = np.array(list_MCC)
            avg_MCC = np.mean(arr_MCC)
            std_MCC = np.std(arr_MCC)
            avg_evalue = sum([np.log10(train[3]) for train in train_report]) / len(
                train_report
            )
            list_cm = [np.array(train[0]) for train in train_report]
            tp_arr = np.array([cm[0][0] for cm in list_cm])
            fp_arr = np.array([cm[0][1] for cm in list_cm])
            fn_arr = np.array([cm[1][0] for cm in list_cm])
            tn_arr = np.array([cm[1][1] for cm in list_cm])
            avg_tp = np.mean(tp_arr)
            avg_fp = np.mean(fp_arr)
            avg_fn = np.mean(fn_arr)
            avg_tn = np.mean(tn_arr)
            avg_cm = (avg_tp, avg_fp, avg_fn, avg_tn)
            f.write(
                "---"
                + "\nAverage MCC: "
                + str(avg_MCC)
                + "\nStandard deviation MCC: "
                + str(std_MCC)
                + "\nAverage E value threshold: "
                + str(avg_evalue)
                + "\nAverage confusion matrix: "
                + str(avg_cm)
            )
            if argv_index <= (len(sys.argv) - 1):
                f.write("\n\n")
            else:
                f.write("\n")


def plot_roc_curve(roc_curve_report_list):
    i = 1
    for roc_curve_report in roc_curve_report_list:
        auc_list, tpr_list, fpr_list = roc_curve_report
        auc_arr = np.array(auc_list)
        avg_auc = np.mean(auc_arr)
        auc_sd = np.std(auc_arr)
        sns.lineplot(
            x=fpr_list,
            y=tpr_list,
            markers=True,
            label=sys.argv[i]
            + "\nMean AUC="
            + str(round(avg_auc, 2))
            + ", SD="
            + "{:.2E}".format(auc_sd),
        )
        i += 1
    plt.ylabel("True Positive Rate")
    plt.xlabel("False Positive Rate")
    plt.savefig("roc_plot.png")
    plt.clf()


def plot_thr_MCC(thr_MCC_report_list):
    for thr_MCC_report in thr_MCC_report_list:
        sns.lineplot(x=thr_MCC_report[0], y=thr_MCC_report[1], markers=True)
    plt.semilogx()
    plt.ylabel("MCC")
    plt.xlabel("Threshold (E value)")
    plt.legend(labels=sys.argv[1:])
    plt.savefig("evalue_MCC_plot.png")
    plt.clf()


def main():
    num_replicates = 1
    numsplits = 5
    start_thr = -30
    stop_thr = 1
    dfs, train_report_list, thr_MCC_report_list, roc_curve_report_list = [], [], [], []
    for path in sys.argv[1:]:
        df = get_df(path)
        train_report, thr_MCC_report, roc_curve_report = train_and_test(
            df, numsplits, num_replicates, start_thr, stop_thr
        )
        train_report_list.append(train_report)
        thr_MCC_report_list.append(thr_MCC_report)
        roc_curve_report_list.append(roc_curve_report)
    plot_roc_curve(roc_curve_report_list)
    plot_thr_MCC(thr_MCC_report_list)
    write_report(train_report_list)


if __name__ == "__main__":
    main()
