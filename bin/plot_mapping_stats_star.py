#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 15:51:22 2020

@author: B. Mika-Gospodorz

Input files: tsv file
Output file: tsv and pdf files
Description: Used to plot mapping statistics for STAR
"""

import argparse

import matplotlib
import pandas as pd

# do not use Xwindows backend
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# function to plot mapping statistics
def plot_mapping_stats(df_comb, no_samples, profile, x_lab, xticks_np, percentage, m):
    # define color palette
    my_cmap2 = matplotlib.cm.get_cmap("tab20")
    my_cmap4 = matplotlib.cm.get_cmap("tab20c")
    color = [
        my_cmap4.colors[12],
        my_cmap4.colors[13],
        "#4984b8",
        "#95d0fc",
        my_cmap4.colors[6],
        my_cmap2.colors[15],
        my_cmap2.colors[14],
    ]
    # reverse sample order
    df_comb = df_comb.loc[reversed(df_comb.index)]
    # make plot
    fig = df_comb.plot(
        kind="barh", stacked=True, figsize=(60, no_samples + 3), legend=True, color=color, width=0.8, fontsize=40
    )
    fig.spines["top"].set_visible(False)
    fig.spines["right"].set_visible(False)
    fig.spines["bottom"].set_visible(True)
    fig.spines["left"].set_visible(False)
    # set x label
    plt.xlabel(x_lab, fontsize=40, labelpad=15)
    # set x label ticks
    plt.xticks(xticks_np)
    if percentage == False:  # set scientific format of label ticks if there is no percentage values
        fig.ticklabel_format(style="scientific", axis="x", scilimits=(m, m))
        fig.xaxis.major.formatter._useMathText = True
        tx = fig.xaxis.get_offset_text()
        tx.set_fontsize(40)
        plt.xticks(xticks_np, rotation=90, fontsize=40)
        plt.tick_params(axis="x", which="major", pad=10)
    # set legent position and format
    fig.legend(
        loc="upper center",
        ncol=4,
        bbox_to_anchor=(0.5, 1.4),
        bbox_transform=plt.gcf().transFigure,
        frameon=False,
        prop={"size": 50},
    )
    # set sample names as y labels
    fig.set_yticklabels(df_comb.index)
    # save used table
    df_comb.to_csv("mapping_stats_" + profile + ".tsv", sep="\t")
    # save plot
    plt.savefig(
        "mapping_stats_star_" + profile + ".pdf",
        dpi=300,
        orientation="landscape",
        transparent=False,
        bbox_inches="tight",
    )


parser = argparse.ArgumentParser(description="""plot mapping statistics""")
parser.add_argument(
    "-i",
    "--input_files",
    metavar="<input_files>",
    default="*.csv",
    help="Path to the outputfiles from mapping statistics ",
)
args = parser.parse_args()


# read mapping statistic table as data frame
df_stats = pd.read_csv(args.input_files, sep="\t", index_col=0)

# extract number of samples
no_samples = df_stats.shape[0]


# calculate percentage of mapping staistics
if "trimmed_reads" in df_stats.columns:
    pathogen_uniquely_mapped_percent = (df_stats["pathogen_uniquely_mapped_reads"] / df_stats["total_raw_reads"]) * 100
    host_uniquely_mapped_percent = (df_stats["host_uniquely_mapped_reads"] / df_stats["total_raw_reads"]) * 100
    pathogen_multi_mapped_percent = (df_stats["pathogen_multi_mapped_reads"] / df_stats["total_raw_reads"]) * 100
    host_multi_mapped_percent = (df_stats["host_multi_mapped_reads"] / df_stats["total_raw_reads"]) * 100
    cross_mapped_percent = (df_stats["cross_mapped_reads"] / df_stats["total_raw_reads"]) * 100
    unmapped_percent = (df_stats["unmapped_reads"] / df_stats["total_raw_reads"]) * 100
    trimmed_percent = (df_stats["trimmed_reads"] / df_stats["total_raw_reads"]) * 100

    df_comb = pd.concat(
        [
            host_uniquely_mapped_percent,
            host_multi_mapped_percent,
            pathogen_uniquely_mapped_percent,
            pathogen_multi_mapped_percent,
            cross_mapped_percent,
            unmapped_percent,
            trimmed_percent,
        ],
        axis=1,
    )
    df_comb.columns = [
        "host uniquely mapped reads",
        "host multi-mapped reads",
        "pathogen uniquely mapped reads",
        "pathogen multi-mapped reads",
        "cross-mapped reads",
        "unmapped reads",
        "trimmed reads",
    ]
else:
    pathogen_uniquely_mapped_percent = (df_stats["pathogen_uniquely_mapped_reads"] / df_stats["processed_reads"]) * 100
    host_uniquely_mapped_percent = (df_stats["host_uniquely_mapped_reads"] / df_stats["processed_reads"]) * 100
    pathogen_multi_mapped_percent = (df_stats["pathogen_multi_mapped_reads"] / df_stats["processed_reads"]) * 100
    host_multi_mapped_percent = (df_stats["host_multi_mapped_reads"] / df_stats["processed_reads"]) * 100
    cross_mapped_percent = (df_stats["cross_mapped_reads"] / df_stats["processed_reads"]) * 100
    unmapped_percent = (df_stats["unmapped_reads"] / df_stats["processed_reads"]) * 100

    df_comb = pd.concat(
        [
            host_uniquely_mapped_percent,
            host_multi_mapped_percent,
            pathogen_uniquely_mapped_percent,
            pathogen_multi_mapped_percent,
            cross_mapped_percent,
            unmapped_percent,
        ],
        axis=1,
    )
    df_comb.columns = [
        "host uniquely mapped reads",
        "host multi-mapped reads",
        "pathogen uniquely mapped reads",
        "pathogen multi-mapped reads",
        "cross-mapped reads",
        "unmapped reads",
    ]

# plot mapping statistics expressed as percentages
plot_mapping_stats(df_comb, no_samples, "samples_percentage", "[ % ]", np.arange(0, 105, step=5), percentage=True, m=0)


# plot mapping statistics
if "trimmed_reads" in df_stats.columns:  # if no. of trimmed reads is defined
    df_comb = pd.concat(
        [
            df_stats["host_uniquely_mapped_reads"],
            df_stats["host_multi_mapped_reads"],
            df_stats["pathogen_uniquely_mapped_reads"],
            df_stats["pathogen_multi_mapped_reads"],
            df_stats["cross_mapped_reads"],
            df_stats["unmapped_reads"],
            df_stats["trimmed_reads"],
        ],
        axis=1,
    )
    df_comb.columns = [
        "host uniquely mapped reads",
        "host multi-mapped reads",
        "pathogen uniquely mapped reads",
        "pathogen multi-mapped reads",
        "cross-mapped reads",
        "unmapped reads",
        "trimmed reads",
    ]
    # calculate max total no. of reads and define x label values, round to 7 digits
    max_limit = round(df_stats["total_raw_reads"].max(), -6)
    if max_limit > 0:  # if total number of reads is higher than 10^6 specify label ticks adjusted to this magnitude
        # divide max_limit by 20 (number of label ticks) to get step
        step2 = int(max_limit / 20)
        # define label ticks
        array_labels = np.arange(0, max_limit + step2, step=step2)
        array_labels2 = np.around(array_labels)
        # set m value to specify format of scientific format of x label ticks
        m = 6
    else:  # define label ticks if number of reads is smaller than 10^6
        step = int(df_stats["total_raw_reads"].max() / 20)
        array_labels = np.arange(0, df_stats["total_raw_reads"].max() + step, step=step)
        array_labels2 = np.around(array_labels)
        # set m value to specify format of scientific format of x label ticks
        m = 0
else:
    df_comb = pd.concat(
        [
            df_stats["host_uniquely_mapped_reads"],
            df_stats["host_multi_mapped_reads"],
            df_stats["pathogen_uniquely_mapped_reads"],
            df_stats["pathogen_multi_mapped_reads"],
            df_stats["cross_mapped_reads"],
            df_stats["unmapped_reads"],
        ],
        axis=1,
    )
    df_comb.columns = [
        "host uniquely mapped reads",
        "host multi-mapped reads",
        "pathogen uniquely mapped reads",
        "pathogen multi-mapped reads",
        "cross-mapped reads",
        "unmapped reads",
    ]
    # calculate max total no. of processed reads and define x label values, round to 7 digits
    max_limit = round(df_stats["processed_reads"].max(), -6)
    if max_limit > 0:  # if total number of reads is higher than 10^6 specify label ticks adjusted to this magnitude
        # divide max_limit by 20 (number of label ticks) to get step
        step2 = int(max_limit / 20)
        # define label ticks
        array_labels = np.arange(0, max_limit + step2, step=step2)
        array_labels2 = np.around(array_labels)
        # set m value to specify format of scientific format of x label ticks
        m = 6
    else:  # define label ticks if number of reads is smaller than 10^6
        step = int(df_stats["processed_reads"].max() / 20)
        array_labels = np.arange(0, df_stats["processed_reads"].max() + step, step=step)
        array_labels2 = np.around(array_labels)
        # set m value to specify format of scientific format of x label ticks
        m = 0

# plot mapping statistics
plot_mapping_stats(df_comb, no_samples, "samples_total_reads", "Number of reads", array_labels2, percentage=False, m=m)
