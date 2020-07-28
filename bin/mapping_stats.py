#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 10:25:47 2019

@author: bozena
"""

import argparse
import pandas as pd

    
def mapping_stats(quantification_table_path,gene_attribute,organism):
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
    col_names = pd.read_csv(quantification_table_path, sep = '\t', nrows=0).columns
    types_dict = {gene_attribute: str}
    types_dict.update({col: float for col in col_names if col not in types_dict})
    quantification_table = pd.read_csv(quantification_table_path, sep = '\t',index_col=0,dtype=types_dict)
    quant_sum = {}
    for sample_name in quantification_table:
        quant_sum.update({sample_name:sum(quantification_table[sample_name])})
    total_counts_pd = pd.DataFrame.from_dict(quant_sum,orient='index')
    total_counts_pd.columns = [organism]
    return total_counts_pd

   

parser = argparse.ArgumentParser(description="""salmon - mapping statistics""")
parser.add_argument("-q_p", "--quantification_table_pathogen", metavar='<quantification_table_pathogen>', default='.', help="Path to the quantification table")
parser.add_argument("-q_h", "--quantification_table_host", metavar='<quantification_table_host>', default='.', help="Path to the quantification table")
parser.add_argument("-total_raw", "--total_no_raw_reads", help="Path to the table with total number of reads for each sample")
parser.add_argument("-total_processed", "--total_no_processed_reads", help="Path to the table with total number of reads for each sample")
parser.add_argument("-m_u", "--mapped_uniquely", metavar='<stats mapped uniquely >', default='.', help="Path to the STAr mapping stats uniquely mapped")
parser.add_argument("-m_m", "--mapped_multi", metavar='<stats multi mapped >', default='.', help="Path to the STAr mapping stats multi mapped")
parser.add_argument("-c_m", "--cross_mapped", metavar='<stats cross mapped >', default='.', help="Path to the STAr mapping stats cross_mapped")
parser.add_argument("-star", "--star_stats", metavar='<stats star >', default='.', help="Path to the STAR mapping /  for htseq") 
parser.add_argument("-star_pr", "--star_processed", metavar='<stats star >', default='.') 
parser.add_argument("-a", "--gene_attribute", default='.', help="gene attribute")
parser.add_argument("-t", "--tool", metavar='<tool>')
parser.add_argument("-o", "--output_dir", metavar='<output>', help="output dir",default='.')
args = parser.parse_args()


if args.tool == 'salmon':
    #mapped reads
    pathogen_total_counts = mapping_stats(args.quantification_table_pathogen,args.gene_attribute,'pathogen')
    host_total_counts = mapping_stats(args.quantification_table_host,args.gene_attribute,'host')
    combined_total_mapped_reads = pd.concat([pathogen_total_counts, host_total_counts], axis=1)
    new_index1 = [sample.split('_NumReads')[0] for sample in combined_total_mapped_reads.index]
    combined_total_mapped_reads.index = new_index1
    combined_total_mapped_reads['total_assigned_reads'] = combined_total_mapped_reads.sum(axis=1)
    

    if args.total_no_raw_reads.endswith('.csv'):
       #read total number of raw reads
       total_reads = pd.read_csv(args.total_no_raw_reads,sep="\t",index_col=0, names=['total_raw_reads'])
       #read total number of processed reads 
       processed_reads_salmon = pd.read_csv(args.total_no_processed_reads,sep="\t",index_col=0, names=['processed_reads'])
       results_df = pd.concat([combined_total_mapped_reads, processed_reads_salmon, total_reads], axis=1) 
       #unmapped reads
       results_df['unassigned_reads'] = results_df['processed_reads'] - results_df['total_assigned_reads']
       results_df['trimmed_reads'] = results_df['total_raw_reads'] - results_df['processed_reads'] 
    else:
       #read total number of processed reads 
       processed_reads_salmon = pd.read_csv(args.total_no_processed_reads,sep="\t",index_col=0, names=['processed_reads'])
       results_df = pd.concat([combined_total_mapped_reads, processed_reads_salmon], axis=1)
       results_df['unassigned_reads'] = results_df['processed_reads'] - results_df['total_assigned_reads']


    results_df2 = results_df.sort_index()
    results_df.to_csv(args.output_dir, sep='\t')
    
    
elif  args.tool == 'salmon_alignment':
    #mapped reads
    pathogen_total_counts = mapping_stats(args.quantification_table_pathogen,args.gene_attribute,'pathogen')
    host_total_counts = mapping_stats(args.quantification_table_host,args.gene_attribute,'host')
    combined_total_assigned_reads = pd.concat([pathogen_total_counts, host_total_counts], axis=1)
    new_index1 = [sample.split('_NumReads')[0] for sample in combined_total_assigned_reads.index]
    combined_total_assigned_reads.index = new_index1
    combined_total_assigned_reads['total_assigned_reads'] = combined_total_assigned_reads.sum(axis=1)
    #mapped reads - extracted from salmon log file
    combined_total_mapped_reads = pd.read_csv(args.total_no_processed_reads,sep="\t",index_col=0, names=['total_mapped_reads'])
    print(combined_total_mapped_reads)
    #read total number of processed reads 
    processed_reads_star = pd.read_csv(args.star_processed,sep="\t",index_col=0, names=['processed_reads'])
    if args.total_no_raw_reads.endswith('.csv'):
          #read total number of raw reads
          total_reads = pd.read_csv(args.total_no_raw_reads,sep="\t",index_col=0, names=['total_raw_reads'])


          results_df = pd.concat([processed_reads_star, total_reads, combined_total_assigned_reads, combined_total_mapped_reads], axis=1) 
       #   unmapped reads
          results_df['unmapped_reads'] = results_df['processed_reads'] - results_df['total_mapped_reads']
          results_df['unassigned_reads'] = results_df['total_mapped_reads'] - results_df['total_assigned_reads']
          results_df['trimmed_reads'] = results_df['total_raw_reads'] - results_df['processed_reads'] 
 
    else:
          results_df = pd.concat([processed_reads_star, combined_total_assigned_reads,  combined_total_mapped_reads], axis=1)
          results_df['unmapped_reads'] = results_df['processed_reads'] - results_df['total_mapped_reads']
          results_df['unassigned_reads'] = results_df['total_mapped_reads'] - results_df['total_assigned_reads']

    results_df2 = results_df.sort_index()
    results_df.to_csv(args.output_dir, sep='\t')
    
    
elif args.tool == 'htseq':
    pathogen_total_counts = mapping_stats(args.quantification_table_pathogen,args.gene_attribute,'pathogen_assigned_reads')
    host_total_counts = mapping_stats(args.quantification_table_host,args.gene_attribute,'host_assigned_reads')
    combined_total_mapped_reads = pd.concat([pathogen_total_counts, host_total_counts], axis=1)
    new_index1 = [sample.split('_NumReads')[0] for sample in combined_total_mapped_reads.index]
    combined_total_mapped_reads.index = new_index1
    combined_total_mapped_reads['total_assigned_reads'] = combined_total_mapped_reads.sum(axis=1)

    #alignment stats
    star_stats = pd.read_csv(args.star_stats,sep="\t",index_col=0 )
    
    results_df = pd.concat([star_stats,combined_total_mapped_reads], axis=1)

    results_df['unassigned_host_reads'] = results_df['host_uniquely_mapped_reads'] - results_df['host_assigned_reads']
    results_df['unassigned_pathogen_reads'] = results_df['pathogen_uniquely_mapped_reads'] - results_df['pathogen_assigned_reads']

    results_df2 = results_df.sort_index()
    results_df.to_csv(args.output_dir, sep='\t')
   
elif args.tool == 'star':
    #mapped reads
     mapped_uniquely = pd.read_csv(args.mapped_uniquely,sep="\t",index_col=0)
     mapped_multi = pd.read_csv(args.mapped_multi,sep="\t",index_col=0)
     cross_mapped = pd.read_csv(args.cross_mapped,sep="\t", header=None, index_col=0, names=['cross_mapped_reads'])
     cross_mapped.index = [sample.replace("_cross_mapped_reads.txt", "") for sample in cross_mapped.index]
     combined_total_mapped_reads = pd.concat([mapped_uniquely, mapped_multi, cross_mapped ], axis=1)
     combined_total_mapped_reads['total_mapped_reads'] = combined_total_mapped_reads.sum(axis=1)

     if args.total_no_raw_reads.endswith('.csv'):
          #read total number of raw reads
          total_reads = pd.read_csv(args.total_no_raw_reads,sep="\t",index_col=0, names=['total_raw_reads'])
          #read total number of processed reads 
          processed_reads_star = pd.read_csv(args.total_no_processed_reads,sep="\t",index_col=0, names=['processed_reads'])
          results_df = pd.concat([processed_reads_star, total_reads, combined_total_mapped_reads], axis=1) 
       #   unmapped reads
          results_df['unmapped_reads'] = results_df['processed_reads'] - results_df['total_mapped_reads']
          results_df['trimmed_reads'] = results_df['total_raw_reads'] - results_df['processed_reads'] 
     else:
          #read total number of processed reads 
          processed_reads_star = pd.read_csv(args.total_no_processed_reads,sep="\t",index_col=0, names=['processed_reads'])
          results_df = pd.concat([processed_reads_star, combined_total_mapped_reads], axis=1)
          results_df['unmapped_reads'] = results_df['processed_reads'] - results_df['total_mapped_reads']
     
     results_df2 = results_df.sort_index()
     results_df2.to_csv(args.output_dir, sep='\t')