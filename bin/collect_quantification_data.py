#!/usr/bin/env python3
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
def collect_quantification_data_HTseq(input_files, profile, gene_attribute):
    # initiate merged quantification data frame
    quant_merged_table = pd.DataFrame()
    for input_file in input_files: #iterate over sample results
        # read quantification data
        quant_table = pd.read_csv(input_file,sep='\t',header=0)
        if input_file == input_files[0]:        # initiate first column of quant_merged_table data frame
            quant_merged_table = quant_table
        else:                                   # extend quant_merged_table data frame with results of new sample
            quant_merged_table = pd.concat([quant_merged_table, quant_table], axis=1)
    quant_merged_table.index.names = [gene_attribute]
    # sort quant_merged_table by column labels 
    quant_merged_table = quant_merged_table.sort_index(axis=1)

    # extract last 5 rows of HTSeq quantification table that contain statistics and save results
    alignment_stats = quant_merged_table['__no_feature':'__alignment_not_unique']
    alignment_stats.to_csv("quantification_stats_" + profile + ".csv", sep='\t')
    # remove statistics from quantification results and save quant_merged_table
    quant_merged_table = quant_merged_table.drop(['__no_feature','__ambiguous','__too_low_aQual','__not_aligned','__alignment_not_unique'])
    quant_merged_table.to_csv("quantification_results_" + profile + ".csv",sep='\t')


# function to merge Salmon quantification results
def collect_quantification_data_salmon(input_files, gene_attribute, organism):
    for input_file in input_files: #iterate over sample results
            if organism == 'host_gene_level': #read tximport results
                quant_table = pd.read_csv(input_file,sep='\t',header=0,index_col=0)
            else: 
                if organism == 'both': # read quantification table containing both host and pathogen transcripts
                    quant_table = pd.read_csv(input_file + '/quant.sf',sep='\t',header=0,index_col=0,dtype={'Name':str,'Length': int,'EffectiveLength':float,'TPM':float,'NumReads': float})
                elif organism == 'pathogen': # read quantification table containing pathogen transcripts
                    quant_table = pd.read_csv(input_file + '/pathogen_quant.sf',sep='\t',header=0,index_col=0,dtype={'Name':str,'Length': int,'EffectiveLength':float,'TPM':float,'NumReads': float})
                elif organism == 'host': # read quantification table containing host transcripts
                      quant_table = pd.read_csv(input_file + '/host_quant.sf',sep='\t',header=0,index_col=0,dtype={'Name':str,'Length': int,'EffectiveLength':float,'TPM':float,'NumReads': float})               
                # create new column names that keep sample name
                name_file = input_file.split('/')
                new_col = [name_file[0] + '_' + column for column in quant_table.columns]
                quant_table.columns = new_col
            if input_file == input_files[0]: #initiate first column
                quant_merged_table = quant_table
            else:                            # extend quant_merged_table data frame with results of new sample
                quant_merged_table = pd.concat([quant_merged_table, quant_table], axis=1)
    quant_merged_table.index.names = [gene_attribute]
    # sort quant_merged_table by column labels 
    quant_merged_table = quant_merged_table.sort_index(axis=1)
    # save results
    if organism == 'both':
        quant_merged_table.to_csv("combined_quant.tsv",sep='\t')        
    elif organism == 'pathogen':
        quant_merged_table.to_csv("pathogen_combined_quant.tsv",sep='\t')
    elif organism == 'host':
        quant_merged_table.to_csv("host_combined_quant.tsv",sep='\t')
    elif organism == 'host_gene_level':
        quant_merged_table.to_csv("host_combined_gene_level.tsv",sep='\t')



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Merges counts for all samples""")
    parser.add_argument("-i", "--input_files", metavar='<input_files>', nargs='+',help="Path to quantification results ")
    parser.add_argument("-q", "--quantifier", metavar='<quantifier>',help="name of quantifier")
    parser.add_argument("-a", "--gene_attribute", help="gene attribute")
    parser.add_argument("-org", "--organism", help="host, pathogen, both, host_gene_level - option avaiable for Salmon")

    args = parser.parse_args()

    # collect either Salmon or HTSeq quantification results
    if args.quantifier == 'htseq':
        collect_quantification_data_HTseq(args.input_files, args.quantifier, args.gene_attribute)
    if args.quantifier == 'salmon':
        collect_quantification_data_salmon(args.input_files, args.gene_attribute, args.organism)

