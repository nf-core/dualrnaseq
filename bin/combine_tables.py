#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 09:26:36 2019

@author: B.Mika-Gospodorz


Input files: list of .txt files with mapping statistics for different samples generated with count_multi_mapped_reads.sh, count_uniquely_mapped_read_pairs, count_uniquely_mapped_read_pairs.sh or count_multi_mapped_read_pairs.sh
Output file: .tsv file with combined statistics from all samples
Description: Used to collect mapping statistics from all samples
"""


import argparse
import pandas as pd

# function to combine quantification results
def combine_stat_tables(input_files, suffix):
    combined_tables = pd.DataFrame() # initialize data frame
    for input_file in input_files: # iterate over input files
        # read file 
        table_read = pd.read_csv(input_file,sep=" ", header=None, index_col=1)
        # extract sample name
        sample_name = str(table_read[0][0])
        table_read = table_read.drop(0,1)
        table_read.columns = [sample_name]
        # add new table to combined_tables data frame
        combined_tables = pd.concat([combined_tables, table_read], axis=1) 
    # transpose table to have 'host' and 'pathogen' labels as column names
    combined_tables = combined_tables.transpose()
    combined_tables.columns = ['host_'+ suffix,'pathogen_' + suffix]
    return combined_tables
        

parser = argparse.ArgumentParser(description="""Combine mapping statistis from all samples""")
parser.add_argument("-i", "--input_files", metavar='<input_files>', nargs='+', default='*.txt', help="Path to mapping statistic results ")
parser.add_argument("-o", "--output_dir", metavar='<output>', help="output dir",default='.')
parser.add_argument("-s", "--suffix", metavar='<suffix>', help="uniquely_mapped/multi_mapped",default='.')
args = parser.parse_args()

# combine results
combineTables = combine_stat_tables(args.input_files, args.suffix)

#save results
combineTables.to_csv(args.output_dir,sep='\t')  

