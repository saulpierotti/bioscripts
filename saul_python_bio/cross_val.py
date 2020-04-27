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


def get_stats(df, thr):
    y_true = df["Class"].values
    y_pred = [1 if (value < thr) else 0 for value in df["Feature"].values]
    confusion_mat = metrics.confusion_matrix(y_true, y_pred, labels=[1, 0])
    ACC = metrics.accuracy_score(y_true, y_pred)
    MCC = metrics.matthews_corrcoef(y_true, y_pred)
    return confusion_mat, ACC, MCC


def get_wrong_predictions(df, thr):
    false_positives = df.loc[(df["Class"] == 0) & (df["Feature"] < thr)]
    false_negatives = df.loc[(df["Class"] == 1) & (df["Feature"] > thr)]
    wrong_pred_report = false_positives, false_negatives
    return wrong_pred_report


def thr_explore(df, start_thr=-100, stop_thr=1, step_thr=1):
    MCC_list, exp_list, ACC_list = [], [], []
    print("Threshold\tACC\tMCC")
    for exp in range(start_thr, stop_thr, step_thr):
        thr = 10 ** exp
        y_true = df["Class"].values
        y_pred = [1 if (value < thr) else 0 for value in df["Feature"].values]
        MCC = metrics.matthews_corrcoef(y_true, y_pred)
        ACC = metrics.accuracy_score(y_true, y_pred)
        MCC_list.append(MCC)
        ACC_list.append(ACC)
        exp_list.append(exp)
        print(thr, ACC, MCC)
    return MCC_list, exp_list, ACC_list


def get_best_thr(df, start_thr=-100, stop_thr=1, step_thr=1):
    MCC_list, exp_list, ACC_list = thr_explore(df, start_thr, stop_thr, step_thr)
    thr_list = [10 ** exp for exp in exp_list]
    best_MCC = -2
    for i in range(len(MCC_list)):
        if MCC_list[i] > best_MCC:
            best_exp_list = [exp_list[i]]
            best_MCC = MCC_list[i]
        elif MCC_list[i] == best_MCC:
            best_exp_list.append(exp_list[i])
    best_exp = np.mean(best_exp_list)
    best_thr = 10 ** best_exp
    print("Selected threshold:", best_thr)
    return best_thr, thr_list, MCC_list, ACC_list


def get_roc_auc(df):
    y_true = df["Class"].values  # ROC curve calculation
    y_score = [-value for value in df["Feature"].values]
    fpr, tpr, thresholds = metrics.roc_curve(y_true, y_score)
    auc = metrics.auc(fpr, tpr)
    return fpr, tpr, auc


def train_and_test(df, numsplits, rand_seed, start_thr=-100, stop_thr=1, step_thr=1):
    (train_report, all_MCCs, all_ACCs, test_auc_list, test_tpr_list, test_fpr_list,) = (
        [],
        [],
        [],
        [],
        [],
        [],
    )
    false_positives, false_negatives = pd.DataFrame(), pd.DataFrame()
    print("\nSplitting the dataset in", numsplits, "Kfolds")
    set_indeces_iterator = model_selection.StratifiedKFold(  # creates an iterator that at each iteration
        n_splits=numsplits,  # returns 2 arrays of indeces for splitting the fold
        shuffle=True,
        random_state=rand_seed,
    ).split(
        df, y=df["Class"]  # the y makes it respect the stratification of the classes
    )
    for j in range(numsplits):
        # training
        print("\nTraining on fold", j + 1)
        df_train, df_test = get_train_test(df, set_indeces_iterator, numsplits)
        thr, thr_list, MCC_list, ACC_list = get_best_thr(
            df_train, start_thr, stop_thr, step_thr
        )
        stats_train = get_stats(df_train, thr)
        stats_test = get_stats(df_test, thr)
        # testing
        train_roc_auc = get_roc_auc(df_train)
        test_roc_auc = get_roc_auc(df_test)
        test_fpr_list += list(test_roc_auc[0])  # for roc plot
        test_tpr_list += list(test_roc_auc[1])  # for roc plot
        test_auc_list.append(test_roc_auc[2])  # for roc plot
        train_report.append(
            (stats_test[0], stats_test[1], stats_test[2], thr, test_roc_auc[2])
        )
        wrong_pred_report = get_wrong_predictions(df_test, thr)
        false_positives = pd.concat([false_positives, wrong_pred_report[0]])
        false_negatives = pd.concat([false_negatives, wrong_pred_report[1]])
        all_MCCs += MCC_list  # for MCC-E value plot
        all_ACCs += ACC_list  # for ACC-E value plot
        print("\nPerformance on fold", j + 1)
        print("    Training set\tTest set")
        print("AUC", train_roc_auc[2], test_roc_auc[2])
        print("ACC", str(stats_train[1]), str(stats_test[1]))
        print("MCC", str(stats_train[2]), str(stats_test[2]))
        print("CM train\n", stats_train[0])
        print("CM test\n", stats_test[0])
        print("\nFalse positives in the test set")
        print(wrong_pred_report[0].to_string(index=False))
        print("\nFalse negatives in the test set")
        print(wrong_pred_report[1].to_string(index=False))
    wrong_pred_report_final = false_positives, false_negatives
    thr_MCC_report = (thr_list * numsplits, all_MCCs, all_ACCs)
    roc_curve_report = (test_auc_list, test_tpr_list, test_fpr_list)
    return train_report, thr_MCC_report, roc_curve_report, wrong_pred_report_final


def get_final_stats(train_report):
    final_stats = {}
    list_AUC = [train[4] for train in train_report]
    arr_AUC = np.array(list_AUC)
    list_ACC = [train[1] for train in train_report]
    arr_ACC = np.array(list_ACC)
    list_MCC = [train[2] for train in train_report]
    arr_MCC = np.array(list_MCC)
    arr_evalue = [train[3] for train in train_report]
    arr_exp = np.log10(arr_evalue)
    avg_exp = np.mean(arr_exp)
    list_cm = [np.array(train[0]) for train in train_report]
    tp_tot = sum([cm[0][0] for cm in list_cm])
    fp_tot = sum([cm[0][1] for cm in list_cm])
    fn_tot = sum([cm[1][0] for cm in list_cm])
    tn_tot = sum([cm[1][1] for cm in list_cm])
    final_stats["avg_AUC"] = np.mean(arr_AUC)
    final_stats["std_AUC"] = np.std(arr_AUC)
    final_stats["avg_ACC"] = np.mean(arr_ACC)
    final_stats["std_ACC"] = np.std(arr_ACC)
    final_stats["avg_MCC"] = np.mean(arr_MCC)
    final_stats["std_MCC"] = np.std(arr_MCC)
    final_stats["avg_evalue"] = 10 ** avg_exp
    final_stats["final_cm"] = (tp_tot, fp_tot, fn_tot, tn_tot)
    return final_stats


def write_report(final_report_list):
    with open("cross_val_report.dat", "w") as f:
        for final_report in final_report_list:
            f.write(final_report)


def get_final_report(train_report, wrong_pred_report, argv_index):
    final_report = ""
    sep = "\t"
    final_report += "# Input file: " + str(sys.argv[argv_index]) + "\n"
    final_report += (
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
        final_report += (
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
    final_stats = get_final_stats(train_report)
    final_report += (
        "---"
        + "\nAverage AUC: "
        + str(final_stats["avg_AUC"])
        + "\nStandard deviation AUC: "
        + str(final_stats["std_AUC"])
        + "\nAverage ACC: "
        + str(final_stats["avg_ACC"])
        + "\nStandard deviation ACC: "
        + str(final_stats["std_ACC"])
        + "\nAverage MCC: "
        + str(final_stats["avg_MCC"])
        + "\nStandard deviation MCC: "
        + str(final_stats["std_MCC"])
        + "\nAverage E value threshold: "
        + str(final_stats["avg_evalue"])
        + "\nFinal confusion matrix (elementwise sum of test confusion matrices): "
        + str(final_stats["final_cm"])
    )
    final_report += (
        "\n\nFalse positives from all the test sets\n"
        + wrong_pred_report[0].drop_duplicates().to_string(index=False)
        + "\n\nFalse negatives from all the test sets\n"
        + wrong_pred_report[1].drop_duplicates().to_string(index=False)
    )
    if argv_index <= (len(sys.argv) - 1):
        final_report += "\n\n"
    else:
        final_report += "\n"
    return final_report


def plot_roc_curve(roc_curve_report_list):
    for i in range(len(sys.argv) - 1):
        auc_list, tpr_list, fpr_list = roc_curve_report_list[i]
        is_beginning = True
        tpr_list = tpr_list[-len(fpr_list) :]
        tpr_list.insert(0, 0.0)
        fpr_list.insert(0, 0.0)
        random_tpr_list = [val for val in fpr_list]
        sns.lineplot(
            x=fpr_list,
            y=tpr_list,
            markers=True,
            color=sns.color_palette()[i],
            estimator=None,
        )
        plt.plot(fpr_list, random_tpr_list, color="r", ls="dashed")
        plt.ylabel("True Positive Rate")
        plt.xlabel("False Positive Rate")
        plt.legend(labels=[sys.argv[i + 1] + " E Value ROC", "Random Classifier"])
        plt.savefig("roc_plot" + sys.argv[i + 1] + ".png")
        plt.clf()


def plot_thr_ACC(thr_MCC_report_list, final_stats_list):
    best_thr, df_list, label_list = [], [], []
    for i in range(len(sys.argv) - 1):
        df = pd.DataFrame(thr_MCC_report_list[i]).T
        df.columns = ["E value", "MCC", "ACC"]
        df["Filename"] = sys.argv[i + 1]
        df_list.append(df)
        best_thr.append(final_stats_list[i]["avg_evalue"])
        label_list.append(sys.argv[i + 1])
    df_all = pd.concat(df_list)
    the_plot = sns.lineplot(
        data=df_all, x="E value", y="ACC", hue="Filename", legend=False
    )
    plt.xlabel("Threshold (E Value)")
    plt.ylabel("ACC")
    plt.legend(labels=label_list)
    plt.semilogx()
    plt.savefig("evalue_ACC_plot.png")
    for i in range(len(sys.argv) - 1):
        plt.axvline(best_thr[i], color=sns.color_palette()[i], ls="dotted")
    plt.savefig("evalue_ACC_plot_thr.png")
    plt.clf()


def plot_thr_MCC(thr_MCC_report_list, final_stats_list):
    best_thr, df_list, label_list = [], [], []
    for i in range(len(sys.argv) - 1):
        df = pd.DataFrame(thr_MCC_report_list[i]).T
        df.columns = ["E value", "MCC", "ACC"]
        df["Filename"] = sys.argv[i + 1]
        df_list.append(df)
        best_thr.append(final_stats_list[i]["avg_evalue"])
        label_list.append(sys.argv[i + 1])
    df_all = pd.concat(df_list)
    the_plot = sns.lineplot(
        x="E value", y="MCC", data=df_all, hue="Filename", legend=False
    )
    plt.semilogx()
    plt.xlabel("Threshold (E Value)")
    plt.ylabel("MCC")
    plt.legend(labels=label_list)
    plt.savefig("evalue_MCC_plot.png")
    label_list = []
    for i in range(len(sys.argv) - 1):
        plt.axvline(best_thr[i], color=sns.color_palette()[i], ls="dotted")
    plt.savefig("evalue_MCC_plot_thr.png")
    plt.clf()


def plot_scatter_all_data(df_list, final_stats_list):
    i = 1
    best_thr = []
    for df in df_list:
        df["Filename"] = sys.argv[i]
        best_thr.append(final_stats_list[i - 1]["avg_evalue"])
        i += 1
    df_all = pd.concat(df_list)
    the_plot = sns.catplot(
        x="Class", y="Feature", data=df_all, dodge=True, hue="Filename"
    )
    the_plot.set_xticklabels(["True negatives", "True positives"])
    the_plot._legend.set_title("")
    plt.ylabel("E Value")
    plt.xlabel("")
    plt.ylim([min(df_all["Feature"].values) / 1e2, max(df_all["Feature"].values) * 1e2])
    plt.semilogy()
    plt.gca().invert_yaxis()
    plt.savefig("scatterplot_all_data.png")
    for i in range(len(sys.argv) - 1):
        plt.axhline(best_thr[i], color=sns.color_palette()[i], ls="dotted")
    plt.savefig("scatterplot_all_data_thr.png")
    plt.clf()
    the_plot_zoomed = sns.catplot(
        x="Class", y="Feature", data=df_all, kind="strip", dodge=True, hue="Filename"
    )
    the_plot_zoomed.set_xticklabels(["True negatives", "True positives"])
    the_plot_zoomed._legend.set_title("")
    plt.ylabel("E Value")
    plt.xlabel("")
    # plt.ylim([1e-15, 1e-5])
    plt.ylim([1e-30, 1e0])
    plt.semilogy()
    plt.gca().invert_yaxis()
    plt.savefig("scatterplot_all_data_zoomed.png")
    for i in range(len(sys.argv) - 1):
        plt.axhline(best_thr[i], color=sns.color_palette()[i], ls="dotted")
    plt.savefig("scatterplot_all_data_zoomed_thr.png")
    plt.clf()


def main():
    numsplits = 5
    rand_seed = 1
    start_thr = -30
    stop_thr = 1
    step_thr = 1
    (
        df_list,
        thr_MCC_report_list,
        roc_curve_report_list,
        final_report_list,
        final_stats_list,
    ) = (
        [],
        [],
        [],
        [],
        [],
    )
    for i in range(1, len(sys.argv)):
        df = get_df(sys.argv[i])
        (
            train_report,
            thr_MCC_report,
            roc_curve_report,
            wrong_pred_report,
        ) = train_and_test(df, numsplits, rand_seed, start_thr, stop_thr, step_thr)
        final_stats = get_final_stats(train_report)
        final_report = get_final_report(train_report, wrong_pred_report, i)
        df_list.append(df)
        thr_MCC_report_list.append(thr_MCC_report)
        roc_curve_report_list.append(roc_curve_report)
        final_report_list.append(final_report)
        final_stats_list.append(final_stats)
    plot_scatter_all_data(df_list, final_stats_list)
    plot_roc_curve(roc_curve_report_list)
    plot_thr_ACC(thr_MCC_report_list, final_stats_list)
    plot_thr_MCC(thr_MCC_report_list, final_stats_list)
    write_report(final_report_list)


if __name__ == "__main__":
    main()
