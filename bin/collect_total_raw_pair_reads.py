#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 14:45:52 2020

@author: bozena
"""
import argparse
import pandas as pd




parser = argparse.ArgumentParser(description="""collect total raw read paires""")
parser.add_argument("-i", "--input_files", metavar='<input_files>', default='*.csv', help="Path to the total_raq_reads_fastq.csv")
args = parser.parse_args()


df_stats = pd.read_csv(args.input_files,sep="\t",index_col=0,header=None)

index_names = df_stats.index

sample_names = [name.rsplit('_', 1)[0] for name in index_names]



df_stats.index = sample_names

df_stats = df_stats.reset_index().drop_duplicates(subset='index', keep='last').set_index('index')


df_stats.to_csv("total_raw_reads_fastq2.csv", sep='\t',header=False)