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


def get_quant_table(input_file: list, organism: str):
    ORGANISM_TO_FILE = {
        "both": input_file + "/quant.sf",
        "pathogen": input_file + "/pathogen_quant.sf",
        "host": input_file + "/host_quant.sf",
    }
    quant_table = pd.read_csv(
        ORGANISM_TO_FILE[organism],
        sep="\t",
        header=0,
        index_col=0,
        dtype={"Name": str, "Length": int, "EffectiveLength": float, "TPM": float, "NumReads": float},
    )
    return quant_table


# function to merge Salmon quantification results
def collect_quantification_data_salmon(input_files, gene_attribute, organism):
    # iterate over sample results
    for input_file in input_files:
        # read tximport results
        if organism == "host_gene_level":
            quant_table = pd.read_csv(input_file, sep="\t", header=0, index_col=0)
        else:
            quant_table = get_quant_table(input_file, organism)
            # create new column names that keep sample name
            name_file = input_file.split("/")
            new_col = [name_file[0] + "_" + column for column in quant_table.columns]
            quant_table.columns = new_col
        if input_file == input_files[0]:  # initiate first column
            quant_merged_table = quant_table
        else:  # extend quant_merged_table data frame with results of new sample
            quant_merged_table = pd.concat([quant_merged_table, quant_table], axis=1)
    quant_merged_table.index.names = [gene_attribute]
    # sort quant_merged_table by column labels
    quant_merged_table = quant_merged_table.sort_index(axis=1)
    # save results
    quant_merged_table.to_csv(f"combined_{organism}.tsv", sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Combine counts from all samples""")
    parser.add_argument(
        "-i", "--input_files", metavar="<input_files>", nargs="+", help="Path to quantification results "
    )
    parser.add_argument("-a", "--gene_attribute", metavar="<gene_attribute>", help="gene attribute")
    parser.add_argument(
        "-org",
        "--organism",
        metavar="<organism>",
        help="host, pathogen, both, host_gene_level - option avaiable for Salmon",
    )

    args = parser.parse_args()
    collect_quantification_data_salmon(args.input_files, args.gene_attribute, args.organism)
