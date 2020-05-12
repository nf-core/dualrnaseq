#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 12:07:56 2020

@author: bozena
"""
import argparse
import pandas as pd


def combine_annotations_quant(quantification_table, annotations_table, gene_attribute, organism):
    col_names = pd.read_csv(quantification_table, sep = '\t', nrows=0).columns
    types_dict = {gene_attribute: str}
    types_dict.update({col: float for col in col_names if col not in types_dict})
    quantification = pd.read_csv(quantification_table,sep="\t",index_col=0, dtype = types_dict)
    annotations = pd.read_csv(annotations_table,sep="\t",index_col=0, dtype='str')
    quant_merged_table = pd.concat([annotations, quantification], axis=1, join = 'inner').sort_index()
    quant_merged_table.index.names = [gene_attribute] 
    if organism == 'pathogen':
        quant_merged_table.to_csv("pathogen_combined_quant_annotations.csv",sep='\t')
    elif organism == 'host':
        quant_merged_table.to_csv("host_combined_quant_annotations.csv",sep='\t')
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Merges the counts for all the samples in a project""")
    parser.add_argument("-q", "--quantification_table", metavar='<input_files>', help="Path to the quantification results ")
    parser.add_argument("-annotations", "--annotations", metavar='<input_files>', help="Path to the annotations extracted from gff file")
    parser.add_argument("-a", "--gene_attribute", help="gene attribute")
    parser.add_argument("-org", "--organism", help="host, pathogen or both")
    
        
    args = parser.parse_args()
    
    
    combine_annotations_quant(args.quantification_table,args.annotations,args.gene_attribute,args.organism)
