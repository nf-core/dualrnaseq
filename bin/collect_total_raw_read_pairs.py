#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 14:45:52 2020

@author: B. Mika-Gospodorz

Input file: tsv file containing number of reads in fastq files created with count_total_reads.sh
Output file: total_raw_read_pairs_fastq.tsv file that contains number of raw read pairs in each sample
Description: Used to collect number of total raw read paires in each sample
"""

import argparse
import pandas as pd


parser = argparse.ArgumentParser(description="""collect total raw read paires""")
parser.add_argument("-i", "--input_files", metavar='<input_files>', default='*.tsv', help="Path to the total_raw_reads_fastq.tsv")
args = parser.parse_args()


# read total_raw_reads_fastq.tsv file
df_stats = pd.read_csv(args.input_files,sep="\t",index_col=0,header=None)

# extract fastq file names
index_names = df_stats.index

# remove '_1' and '_2' suffixes from fastq file names
sample_names = [name.rsplit('_', 1)[0] for name in index_names]

# define new sample names and remove duplicated results (each sample contains results from both *_1 and *_2 fastq files
df_stats.index = sample_names
df_stats = df_stats.reset_index().drop_duplicates(subset='index', keep='last').set_index('index')

 # save results
df_stats.to_csv("total_raw_read_pairs_fastq.tsv", sep='\t',header=False)
