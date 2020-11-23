#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 10:25:47 2019

@author: B. Mika-Gospodorz

Input files: host and pathogen quantification tables (e.g. pathogen_quant_salmon.tsv, host_quantification_uniquely_mapped_htseq.tsv), raw read statistics, star statistics 
Output file: tsv file
Description: Used to collect mapping and quantifications statistics for STAR, HTSeq, Salmon or HTSeq results. 
"""

import argparse
import pandas as pd


# function to sum up number of mapped reads from quantification table
def mapping_stats(quantification_table_path,gene_attribute,organism):
    # read quantification table
    col_names = pd.read_csv(quantification_table_path, sep = '\t', nrows=0).columns
    types_dict = {gene_attribute: str}
    types_dict.update({col: float for col in col_names if col not in types_dict})
    quantification_table = pd.read_csv(quantification_table_path, sep = '\t',index_col=0,dtype=types_dict)
    quantification_table = quantification_table.fillna(0)
    # initialize dict. 
    quant_sum = {}
    for sample_name in quantification_table: # iterate over columns of quantification_table  
        if 'NumReads' in sample_name:        # for each column (sample) with 'NumReads' sum up the number of reads and add into quant_sum dict. 
            quant_sum.update({sample_name:sum(quantification_table[sample_name])})
    # create data frame from dict.
    total_counts_pd = pd.DataFrame.from_dict(quant_sum,orient='index')
    total_counts_pd.columns = [organism]
    return total_counts_pd


parser = argparse.ArgumentParser(description="""collects and generates mapping statistics""")
parser.add_argument("-q_p", "--quantification_table_pathogen", metavar='<quantification_table_pathogen>', default='.', help="Path to pathogen quantification table; Salmon and HTSeq")
parser.add_argument("-q_h", "--quantification_table_host", metavar='<quantification_table_host>', default='.', help="Path to host quantification table; Salmon and HTSeq")
parser.add_argument("-total_raw", "--total_no_raw_reads",metavar='<total_no_raw_reads>',default='.',help="Path to table with total number of raw reads for each sample; Salmon and STAR")
parser.add_argument("-total_processed", "--total_no_processed_reads", metavar='<total_no_processed_reads>', help="Path to table with total number of processed reads by STAR or Salmon")
parser.add_argument("-m_u", "--mapped_uniquely", metavar='<stats mapped uniquely >', default='.', help="Path to STAR mapping stats of uniquely mapped reads; STAR")
parser.add_argument("-m_m", "--mapped_multi", metavar='<stats multi mapped >', default='.', help="Path to STAR mapping stats of multi mapped reads; STAR")
parser.add_argument("-c_m", "--cross_mapped", metavar='<stats cross mapped >', default='.', help="Path to STAR mapping stats of cross_mapped reads; STAR")
parser.add_argument("-star", "--star_stats", metavar='<stats star >', default='.', help="Path to mapping statistics of STAR; HTSeq") 
parser.add_argument("-star_pr", "--star_processed", metavar='<stats star_processed >', default='.',help="Path to STAR stats of processed reads; Salmon in alignment-based mode") 
parser.add_argument("-a", "--gene_attribute", default='.', help="gene attribute used in quantification; Salmon and HTSeq")
parser.add_argument("-t", "--tool", metavar='<tool>', help="salmon, salmon_alignment, htseq or star")
parser.add_argument("-o", "--output_dir", metavar='<output>', help="output dir",default='.')
args = parser.parse_args()

# collect statistics for Salmon Selective Alignment mode
if args.tool == 'salmon': 
    # collect assigned pathogen reads
    pathogen_total_counts = mapping_stats(args.quantification_table_pathogen,args.gene_attribute,'pathogen')
    # collect assigned host reads
    host_total_counts = mapping_stats(args.quantification_table_host,args.gene_attribute,'host')
    # combine host and pathogen mapped reads
    combined_total_mapped_reads = pd.concat([pathogen_total_counts, host_total_counts], axis=1)
    # rename colnames - remove '_NumReads' suffix
    new_index1 = [sample.split('_NumReads')[0] for sample in combined_total_mapped_reads.index]
    combined_total_mapped_reads.index = new_index1
    # calculate total mapped reads 
    combined_total_mapped_reads['total_mapped_reads'] = combined_total_mapped_reads.sum(axis=1)
    if args.total_no_raw_reads.endswith('.tsv'): # if tsv table is defined in total_no_raw_reads argument
       # read total number of raw reads
       total_reads = pd.read_csv(args.total_no_raw_reads,sep="\t",index_col=0, names=['total_raw_reads'])
       # read total number of reads processed by Salmon
       processed_reads_salmon = pd.read_csv(args.total_no_processed_reads,sep="\t",index_col=0, names=['processed_reads'])
       # combine statistics
       results_df = pd.concat([combined_total_mapped_reads, processed_reads_salmon, total_reads], axis=1) 
       # calculate unmapped reads
       results_df['unmapped_reads'] = results_df['processed_reads'] - results_df['total_mapped_reads']
       # calculate trimmed reads
       results_df['trimmed_reads'] = results_df['total_raw_reads'] - results_df['processed_reads'] 
    else: # if tsv table is not defined in total_no_raw_reads argument 
       # read total number of reads processed by Salmon
       processed_reads_salmon = pd.read_csv(args.total_no_processed_reads,sep="\t",index_col=0, names=['processed_reads'])
       results_df = pd.concat([combined_total_mapped_reads, processed_reads_salmon], axis=1)
       # calculate unmapped reads
       results_df['unmapped_reads'] = results_df['processed_reads'] - results_df['total_mapped_reads']
    # save results
    results_df2 = results_df.sort_index()
    results_df2.to_csv(args.output_dir, sep='\t')
    
# collect statistics for Salmon alignment-based mode
elif  args.tool == 'salmon_alignment':
    # collect assigned pathogen reads
    pathogen_total_counts = mapping_stats(args.quantification_table_pathogen,args.gene_attribute,'pathogen')
    # collect assigned host reads
    host_total_counts = mapping_stats(args.quantification_table_host,args.gene_attribute,'host')
    # combine host and pathogen mapped reads
    combined_total_assigned_reads = pd.concat([pathogen_total_counts, host_total_counts], axis=1)
    # rename colnames - remove '_NumReads' suffix
    new_index1 = [sample.split('_NumReads')[0] for sample in combined_total_assigned_reads.index]
    combined_total_assigned_reads.index = new_index1
    # calculate total assigned reads 
    combined_total_assigned_reads['total_assigned_reads'] = combined_total_assigned_reads.sum(axis=1)
    # extracted mapped reads from salmon log file
    combined_total_mapped_reads = pd.read_csv(args.total_no_processed_reads,sep="\t",index_col=0, names=['total_mapped_reads'])
    #read total number of reads processed by STAR
    processed_reads_star = pd.read_csv(args.star_processed,sep="\t",index_col=0, names=['processed_reads'])
    if args.total_no_raw_reads.endswith('.tsv'):
          #read total number of raw reads
          total_reads = pd.read_csv(args.total_no_raw_reads,sep="\t",index_col=0, names=['total_raw_reads'])
          results_df = pd.concat([processed_reads_star, total_reads, combined_total_assigned_reads, combined_total_mapped_reads], axis=1) 
         # calculate unmapped reads
          results_df['unmapped_reads'] = results_df['processed_reads'] - results_df['total_mapped_reads']
         # calculate unassigned reads
          results_df['unassigned_reads'] = results_df['total_mapped_reads'] - results_df['total_assigned_reads']
         # calculate trimmed reads
          results_df['trimmed_reads'] = results_df['total_raw_reads'] - results_df['processed_reads'] 
    else:
          results_df = pd.concat([processed_reads_star, combined_total_assigned_reads,  combined_total_mapped_reads], axis=1)
         # calculate unmapped reads
          results_df['unmapped_reads'] = results_df['processed_reads'] - results_df['total_mapped_reads']
         # calculate unassigned reads
          results_df['unassigned_reads'] = results_df['total_mapped_reads'] - results_df['total_assigned_reads']
    # save results
    results_df2 = results_df.sort_index()
    results_df2.to_csv(args.output_dir, sep='\t')
    
# collect statistics for HTSeq
elif args.tool == 'htseq':
    # collect assigned pathogen reads
    pathogen_total_counts = mapping_stats(args.quantification_table_pathogen,args.gene_attribute,'pathogen_assigned_reads')
    # collect assigned host reads
    host_total_counts = mapping_stats(args.quantification_table_host,args.gene_attribute,'host_assigned_reads')
    combined_total_mapped_reads = pd.concat([pathogen_total_counts, host_total_counts], axis=1)
    # rename colnames - remove '_NumReads' suffix
    new_index1 = [sample.split('_NumReads')[0] for sample in combined_total_mapped_reads.index]
    combined_total_mapped_reads.index = new_index1
    # calculate total assigned reads 
    combined_total_mapped_reads['total_assigned_reads'] = combined_total_mapped_reads.sum(axis=1)
    # read alignment statistics
    star_stats = pd.read_csv(args.star_stats,sep="\t",index_col=0 )
    # combine statistics
    results_df = pd.concat([star_stats,combined_total_mapped_reads], axis=1)
    # calculate unassigned host reads
    results_df['unassigned_host_reads'] = results_df['host_uniquely_mapped_reads'] - results_df['host_assigned_reads']
    # calculate unassigned pathogen reads
    results_df['unassigned_pathogen_reads'] = results_df['pathogen_uniquely_mapped_reads'] - results_df['pathogen_assigned_reads']
    # save results
    results_df2 = results_df.sort_index()
    results_df2.to_csv(args.output_dir, sep='\t')
   
# collect statistics for STAR
elif args.tool == 'star':
    # read total number of uniquely mapped reads
     mapped_uniquely = pd.read_csv(args.mapped_uniquely,sep="\t",index_col=0)
    # read total number of multi-mapped reads
     mapped_multi = pd.read_csv(args.mapped_multi,sep="\t",index_col=0)
    # read total number of cross-mapped reads
     cross_mapped = pd.read_csv(args.cross_mapped,sep="\t", header=None, index_col=0, names=['cross_mapped_reads'])
     cross_mapped.index = [sample.replace("_cross_mapped_reads.txt", "") for sample in cross_mapped.index]
    # combine statistics
     combined_total_mapped_reads = pd.concat([mapped_uniquely, mapped_multi, cross_mapped ], axis=1)
    # calculate total mapped reads 
     combined_total_mapped_reads['total_mapped_reads'] = combined_total_mapped_reads.sum(axis=1)
     if args.total_no_raw_reads.endswith('.tsv'):
          # read total number of raw reads
          total_reads = pd.read_csv(args.total_no_raw_reads,sep="\t",index_col=0, names=['total_raw_reads'])
          # read total number of processed reads 
          processed_reads_star = pd.read_csv(args.total_no_processed_reads,sep="\t",index_col=0, names=['processed_reads'])
          # combine statistics
          results_df = pd.concat([processed_reads_star, total_reads, combined_total_mapped_reads], axis=1) 
          # calculate unmapped reads
          results_df['unmapped_reads'] = results_df['processed_reads'] - results_df['total_mapped_reads']
          # calculate trimmed reads
          results_df['trimmed_reads'] = results_df['total_raw_reads'] - results_df['processed_reads'] 
     else:
          #read total number of processed reads 
          processed_reads_star = pd.read_csv(args.total_no_processed_reads,sep="\t",index_col=0, names=['processed_reads'])
          # combine statistics
          results_df = pd.concat([processed_reads_star, combined_total_mapped_reads], axis=1)
          # calculate unmapped reads
          results_df['unmapped_reads'] = results_df['processed_reads'] - results_df['total_mapped_reads']
     #save results
     results_df2 = results_df.sort_index()
     results_df2.to_csv(args.output_dir, sep='\t')
