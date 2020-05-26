#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 09:26:36 2019

@author: bozena
"""


import argparse
import pandas as pd


def combine_stat_tables(input_files, suffix):
    combined_tables = pd.DataFrame()
    for input_file in input_files:
        table_read = pd.read_csv(input_file,sep=" ", header=None, index_col=1)
        sample_name = str(table_read[0][0])
        table_read = table_read.drop(0,1)
        table_read.columns = [sample_name]
        combined_tables = pd.concat([combined_tables, table_read], axis=1) 
    combined_tables = combined_tables.transpose()
    combined_tables.columns = ['host_'+ suffix,'pathogen_' + suffix]
    return combined_tables
        

parser = argparse.ArgumentParser(description="""Merges mapping stats for all the samples in a project""")
parser.add_argument("-i", "--input_files", metavar='<input_files>', nargs='+', default='*.txt', help="Path to the outputfiles from mapping statistics ")
parser.add_argument("-o", "--output_dir", metavar='<output>', help="output dir",default='.')
parser.add_argument("-s", "--suffix", metavar='<suffix>', help="uniquely_mapped/multi_mapped",default='.')
args = parser.parse_args()

combineTables = combine_stat_tables(args.input_files, args.suffix)
combineTables.to_csv(args.output_dir,sep='\t')  

