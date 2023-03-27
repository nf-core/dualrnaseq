#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 15:00:57 2019

@author: B. Mika-Gospodorz

Input file: list of quantification tables
Output files: quantification_stats_*.tsv, quantification_results_*.tsv (HTSeq option) or
Description: Used to merge quantification results from all samples
"""

import argparse

import pandas as pd


# function to merge HTSeq quantification results
def collect_quantification_data_HTseq(input_files: list, gene_attribute):
    # initiate merged quantification data frame
    quant_merged_table = pd.DataFrame()
    # iterate over sample results
    for input_file in input_files:
        # read quantification data
        quant_table = pd.read_csv(input_file, sep="\t", header=0)
        if input_file == input_files[0]:  # initiate first column of quant_merged_table data frame
            quant_merged_table = quant_table
        else:  # extend quant_merged_table data frame with results of new sample
            quant_merged_table = pd.concat([quant_merged_table, quant_table], axis=1)
    quant_merged_table.index.names = [gene_attribute]
    # sort quant_merged_table by column labels
    quant_merged_table = quant_merged_table.sort_index(axis=1)

    # extract last 5 rows of HTSeq quantification table that contain statistics and save results
    alignment_stats = quant_merged_table["__no_feature":"__alignment_not_unique"]
    alignment_stats.to_csv("quantification_stats_htseq.tsv", sep="\t")
    # remove statistics from quantification results and save quant_merged_table
    quant_merged_table = quant_merged_table.drop(
        ["__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique"]
    )
    quant_merged_table.to_csv("quantification_results_htseq.tsv", sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Combine counts from all samples""")
    parser.add_argument(
        "-i", "--input_files", metavar="<input_files>", nargs="+", help="Path to quantification results "
    )
    parser.add_argument("-a", "--gene_attribute", metavar="<gene_attribute>", help="gene attribute")
    args = parser.parse_args()

    # collect either Salmon or HTSeq quantification results
    collect_quantification_data_HTseq(args.input_files, args.gene_attribute)
