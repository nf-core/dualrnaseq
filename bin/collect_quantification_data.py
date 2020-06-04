#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 15:00:57 2019

@author: bozena
"""

import argparse
import pandas as pd



def collect_quantification_data_HTseq(input_files, profile, gene_attribute):
    quant_table = pd.DataFrame()
    quant_merged_table = pd.DataFrame()
    alignment_stats = pd.DataFrame()
    for input_file in input_files:
        quant_table = pd.read_csv(input_file,sep='\t',header=0)
		#initiate a  first  column
        if input_file == input_files[0]:
            quant_merged_table = quant_table
        else:
            quant_merged_table = pd.concat([quant_merged_table, quant_table], axis=1)
    quant_merged_table.index.names = [gene_attribute]
    quant_merged_table = quant_merged_table.sort_index(axis=1)
    alignment_stats = quant_merged_table['__no_feature':'__alignment_not_unique']
    alignment_stats.to_csv("quantification_stats_" + profile + ".csv", sep='\t')
    quant_merged_table = quant_merged_table.drop(['__no_feature','__ambiguous','__too_low_aQual','__not_aligned','__alignment_not_unique'])
    quant_merged_table.to_csv("quantification_results_" + profile + ".csv",sep='\t')


def collect_quantification_data_salmon(input_files, gene_attribute, organism):
    for input_file in input_files:
            if organism == 'both':
                quant_table = pd.read_csv(input_file + '/quant.sf',sep='\t',header=0,index_col=0,dtype={'Name':str,'Length': int,'EffectiveLength':float,'TPM':float,'NumReads': float})
            elif organism == 'pathogen':
                quant_table = pd.read_csv(input_file + '/pathogen_quant.sf',sep='\t',header=0,index_col=0,dtype={'Name':str,'Length': int,'EffectiveLength':float,'TPM':float,'NumReads': float})
            elif organism == 'host':
                  quant_table = pd.read_csv(input_file + '/host_quant.sf',sep='\t',header=0,index_col=0,dtype={'Name':str,'Length': int,'EffectiveLength':float,'TPM':float,'NumReads': float})               
            elif organism == 'host_gene_level':
                  quant_table = pd.read_csv(input_file + '/host_quant_gene_level.sf',sep='\t',header=0,index_col=0,dtype={'Name':str,'Length': int,'TPM':float,'NumReads': float})               

            name_file = input_file.split('/')
            #create new column names that keep sample name
            new_col = [name_file[0] + '_' + column for column in quant_table.columns]
            quant_table.columns = new_col
            if input_file == input_files[0]:
                #initiate first columns
                  quant_merged_table = quant_table
            else:
                  quant_merged_table = pd.concat([quant_merged_table, quant_table], axis=1)
    quant_merged_table.index.names = [gene_attribute]
    quant_merged_table = quant_merged_table.sort_index(axis=1)
    if organism == 'both':
        quant_merged_table.to_csv("combined_quant.csv",sep='\t')        
    elif organism == 'pathogen':
        quant_merged_table.to_csv("pathogen_combined_quant.csv",sep='\t')
    elif organism == 'host':
        quant_merged_table.to_csv("host_combined_quant.csv",sep='\t')
    elif organism == 'host_gene_level':
        quant_merged_table.to_csv("host_combined_gene_level.csv",sep='\t')
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Merges the counts for all the samples in a project""")
    parser.add_argument("-i", "--input_files", metavar='<input_files>', nargs='+', default='*.count_u_m',
                                   help="Path to the outputfiles from HTseq ")
    parser.add_argument("-q", "--quantifier", metavar='<quantifier>',help="name of the quantifier")
    parser.add_argument("-a", "--gene_attribute", help="gene attribute")
    parser.add_argument("-org", "--organism", help="host, pathogen or both")
    parser.add_argument("-p", "--profile", help="uniquely_mapped/multi_mapped",default='')    
        
    args = parser.parse_args()
    
    if args.quantifier == 'htseq':
        collect_quantification_data_HTseq(args.input_files, args.profile, args.gene_attribute)
    if args.quantifier == 'salmon':
        collect_quantification_data_salmon(args.input_files, args.gene_attribute,args.organism)

