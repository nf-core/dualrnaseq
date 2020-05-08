#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 10:25:47 2019

@author: bozena
"""

import argparse
import pandas as pd

    
def mapping_stats_salmon(quantification_table_path,gene_attribute,organism):
    col_names = pd.read_csv(quantification_table_path, sep = '\t', nrows=0).columns
    types_dict = {gene_attribute: str}
    types_dict.update({col: float for col in col_names if col not in types_dict})
    quantification_table = pd.read_csv(quantification_table_path, sep = '\t',index_col=0,dtype=types_dict)
    quantification_table = quantification_table.fillna(0)
    quant_sum = {}
    for sample_name in quantification_table:
        if 'NumReads' in sample_name:
            quant_sum.update({sample_name:sum(quantification_table[sample_name])})
    total_counts_pd = pd.DataFrame.from_dict(quant_sum,orient='index')
    total_counts_pd.columns = [organism]
    return total_counts_pd

def mapping_stats_htseq(quantification_table_path,gene_attribute,organism):
    quantification_table = pd.read_csv(quantification_table_path, sep = '\t',index_col=0,dtype={'gene_ID':str})
    quant_sum = {}
    for sample_name in quantification_table:
        quant_sum.update({sample_name:sum(quantification_table[sample_name])})
    total_counts_pd = pd.DataFrame.from_dict(quant_sum,orient='index')
    total_counts_pd.columns = [organism]
    return total_counts_pd

   



parser = argparse.ArgumentParser(description="""salmon - mapping statistics""")
parser.add_argument("-q_p", "--quantification_table_pathogen", metavar='<quantification_table_pathogen>', help="Path to the quantification table")
parser.add_argument("-q_h", "--quantification_table_host", metavar='<quantification_table_host>', help="Path to the quantification table")
parser.add_argument("-total_raw", "--total_no_raw_reads", help="Path to the table with total number of reads for each sample")
parser.add_argument("-total_processed", "--total_no_processed_reads", help="Path to the table with total number of reads for each sample")
parser.add_argument("-a", "--gene_attribute", help="gene attribute")
parser.add_argument("-t", "--tool", metavar='<tool>')
parser.add_argument("-o", "--output_dir", metavar='<output>', help="output dir",default='.')
args = parser.parse_args()

if args.tool == 'salmon':
    #mapped reads
    pathogen_total_counts = mapping_stats_salmon(args.quantification_table_pathogen,args.gene_attribute,'pathogen')
    host_total_counts = mapping_stats_salmon(args.quantification_table_host,args.gene_attribute,'host')
    combined_total_mapped_reads = pd.concat([pathogen_total_counts, host_total_counts], axis=1)
    new_index1 = [sample.split('_NumReads')[0] for sample in combined_total_mapped_reads.index]
    combined_total_mapped_reads.index = new_index1
    combined_total_mapped_reads['total_mapped_reads'] = combined_total_mapped_reads.sum(axis=1)
    
    #read total number of reads
    processed_reads_salmon = pd.read_csv(args.total_no_processed_reads,sep="\t",index_col=0, names=['processed_reads'])
    
    total_reads = pd.read_csv(args.total_no_raw_reads,sep="\t",index_col=0, names=['total_raw_reads'])
    new_index2 = [sample.split('.fastq')[0] for sample in total_reads.index]
    total_reads.index = new_index2
   # total_reads.columns = ['total_no_reads']

    
    results_df = pd.concat([combined_total_mapped_reads, processed_reads_salmon, total_reads], axis=1)
    #unmapped reads
    results_df['unmapped_reads'] = results_df['processed_reads'] - results_df['total_mapped_reads']
    results_df['trimmed_reads'] = results_df['total_raw_reads'] - results_df['processed_reads']

    results_df.to_csv(args.output_dir, sep='\t')
    
elif args.tool == 'htseq':
    pathogen_total_counts = mapping_stats_htseq(args.quantification_table_pathogen,args.gene_attribute,'pathogen')
    host_total_counts = mapping_stats_htseq(args.quantification_table_host,args.gene_attribute,'host')
    combined_total_mapped_reads = pd.concat([pathogen_total_counts, host_total_counts], axis=1)
    combined_total_mapped_reads['total'] = combined_total_mapped_reads.sum(axis=1)
    combined_total_mapped_reads.to_csv(args.output_dir, sep='\t')

