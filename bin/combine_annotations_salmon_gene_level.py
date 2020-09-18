#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 16:22:50 2020

@author: B.Mika-Gospodrz

Input file: host_combined_gene_level.tsv file that contains gene-level salmon estimates
            .tsv file with transcript annotations extracted from gff using extract_annotations_from_gff.py
Output files: *combined_quant_gene_level_annotations.tsv with gene annotations and quantification results
Description: Used to combine annotations with gene-level salmon quantification results
"""

import argparse
import pandas as pd


def combine_annotations_quant_salmon_gene(quantification_table, annotations_table, gene_attribute, organism):
    # read quantification results
    col_names = pd.read_csv(quantification_table, sep = '\t', nrows=0).columns
    types_dict = {gene_attribute: str}
    types_dict.update({col: float for col in col_names if col not in types_dict})
    quantification = pd.read_csv(quantification_table,sep="\t",index_col=0, dtype = types_dict)

    # read annotations 
    annotations = pd.read_csv(annotations_table,sep="\t",index_col=0, dtype='str')
    # keep gene annotations for gene level estimates 
    genes = annotations[['gene_id','gene_name','gene_type']]
    unique_genes = genes.drop_duplicates(subset='gene_id', keep="first")
    unique_genes = unique_genes.set_index('gene_id')

    # combine gene annotations and gene -level quantification results
    quant_merged_table = pd.concat([unique_genes, quantification], axis=1, join = 'inner').sort_index()
    quant_merged_table.index.names = [gene_attribute] 
    # save the results
    if organism == 'pathogen':
         quant_merged_table.to_csv("pathogen_combined_quant_gene_level_annotations.tsv",sep='\t')
    elif organism == 'host':
         quant_merged_table.to_csv("host_combined_quant_gene_level_annotations.tsv",sep='\t')
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Merges the counts for all the samples""")
    parser.add_argument("-q", "--quantification_table", metavar='<input_files>', help="Path to the quantification results ")
    parser.add_argument("-annotations", "--annotations", metavar='<input_files>', help="Path to the annotations extracted from gff file")
    parser.add_argument("-a", "--gene_attribute", help="gene attribute")
    parser.add_argument("-org", "--organism", help="host, pathogen or both")
    
        
    args = parser.parse_args()
    
    # combine reasults of all samples
    combine_annotations_quant_salmon_gene(args.quantification_table,args.annotations,args.gene_attribute,args.organism)
