#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 10:20:19 2019

@author: B. Mika-Gospodorz

Input files: tsv quantification table, tsv annotation table created with extract_annotations_from_gff.py , optional: RNA_classes_to_replace.tsv to group subclasses into main RNA classes for calculating the statistics (available for host)
Output files: RNA_classes_percentage_*.tsv, RNA_classes_sum_counts*.tsv, gene_types_groups_*.tsv
Description: Used to calculate mapping statistics of RNA classes. 

"""

import argparse
import pandas as pd
from collections import defaultdict



# function to create distionary
def create_dictionary_of_RNA_classes_from_set(keys,values):
    zipped = list(zip(keys,values))
    return dict(zipped)
        

# function to percentage of RNA classes
def calculate_percentage_RNA_classes(RNA_classes_sum_counts):
    return((RNA_classes_sum_counts/RNA_classes_sum_counts.sum())*100)


# function to add counts to 'gene_types' dict. with gene annotations
def add_counts_to_genes(quantification_df,gene_types,feature):  
    # add counts from quantification for each gene/transcripts present in annotation dict., if gene_id/transcript_id is missing in quantification, set counts to 0.
    for id_gene in gene_types:
        if id_gene[feature] in quantification_df.index:
            id_gene['counts'] = quantification_df.at[id_gene[feature]]
        else:
            id_gene['counts'] = 0 
    return(gene_types)

# function to create table with transcript/gene annotations, RNA class assigned to each entry, and counts
def get_gene_type_host(gene_RNAtype_counts, set_RNA_classes,gene_attribute):
    # extract main RNA class and subclasses of RNA from dict. 
    for RNA_class_group, RNA_classes in list(set_RNA_classes.items()):
        # identify main RNA class for each gene/transcript
        for gene in gene_RNAtype_counts:
            if gene['gene_type'] in RNA_classes:
                gene['gene_type'] = RNA_class_group
    # create data frame from dict. and save it
    gene_RNAtype_counts_df = pd.DataFrame.from_dict(gene_RNAtype_counts) 
    gene_RNAtype_counts_df.to_csv("host_gene_types_groups_" + gene_attribute.split('_')[0] + ".tsv",sep = '\t', index = False)
    
# function to sum up number of reads for each RNA class
def sum_counts_for_each_RNA_class(gene_RNAtype_counts, set_RNA_classes):
    RNA_classes_sum_counts = []
    # iterate over RNA classes {'main RNA class': ['RNA subclasses']}
    for RNA_class_group, RNA_classes in list(set_RNA_classes.items()):
        # extract genes from particular RNA class
        all_genes_in_RNA_class = ([gene for gene in gene_RNAtype_counts if gene['gene_type'] in RNA_class_group])
        # sum counts for all genes from the same RNA class
        class_sum_counts = sum([gene_counts['counts'] for gene_counts in all_genes_in_RNA_class])
        # add dictionary of RNA class and summed number of reads
        RNA_classes_sum_counts.append({'gene_type': RNA_class_group, 'sum_counts':float(class_sum_counts)})
        # create data frame from dict. of RNA classes and counts
        RNA_classes_sum_counts_df = pd.DataFrame.from_dict(RNA_classes_sum_counts)
        RNA_classes_sum_counts_df.index = RNA_classes_sum_counts_df['gene_type']
        # remove 'gene_type' header 
        RNA_classes_sum_counts_df = RNA_classes_sum_counts_df.drop('gene_type', axis=1)
    return(RNA_classes_sum_counts_df)


# function to calculate RNA class statistics for each sample
def RNA_classes_for_each_sample(quantification_table_path,gene_types,set_RNA_classes,organism, quantifier, gene_attribute):
    # extract column names
    col_names = pd.read_csv(quantification_table_path, sep = '\t', nrows=0).columns
    # define data types in each column of quantification table
    types_dict = {gene_attribute: str}
    types_dict.update({col: float for col in col_names if col not in types_dict})
    # read quantification table
    quantification_table = pd.read_csv(quantification_table_path, sep = '\t',index_col=0, dtype=types_dict)
    quantification_table = quantification_table.fillna(0)
    # initialize RNA class statistic data frames 
    df_percentage_RNA_all_samples  = pd.DataFrame(columns = set_RNA_classes,index=['NA'])
    df_RNA_classes_sum_counts_all_samples  = pd.DataFrame(columns = set_RNA_classes,index=['NA'])
    # calculate statistics for each sample 
    for sample_name in quantification_table:
        # extract columns with Number of Reads
        if 'NumReads' in sample_name:
            # add counts for each gene in dict. with annotations 
            gene_RNA_type_counts = add_counts_to_genes(quantification_table[sample_name],gene_types,gene_attribute)
            # create table with annotations and RNA classes assigned to genes/transcript 
            if organism == 'host':
                get_gene_type_host(gene_RNA_type_counts,set_RNA_classes,gene_attribute)
            # create data frame with summed number of reads for each RNA class
            sum_counts_RNA_class_df = sum_counts_for_each_RNA_class(gene_RNA_type_counts,set_RNA_classes)
            df_RNA_classes_sum_counts = sum_counts_RNA_class_df.sort_values('sum_counts', ascending=0)
            # define sample name as column
            df_RNA_classes_sum_counts = pd.DataFrame(df_RNA_classes_sum_counts).T
            df_RNA_classes_sum_counts.index = [sample_name]
            df_RNA_classes_sum_counts = df_RNA_classes_sum_counts.T
            # calculate percentage of RNA class
            percentage_RNA_class = calculate_percentage_RNA_classes(df_RNA_classes_sum_counts)
            # add results of RNA class to data frames with mapping statistics
            df_RNA_classes_sum_counts_all_samples = pd.concat([df_RNA_classes_sum_counts_all_samples, df_RNA_classes_sum_counts.T],sort=True)  
            df_percentage_RNA_all_samples = pd.concat([df_percentage_RNA_all_samples, percentage_RNA_class.T],sort=True)
    # remove missing values from data frames
    df_percentage_RNA_all_samples = df_percentage_RNA_all_samples.dropna()
    df_RNA_classes_sum_counts_all_samples = df_RNA_classes_sum_counts_all_samples.dropna()
    # check it data frame contains any non zero value. If so, teturn 'true' and plot RNA class statistics in another process
    non_zero = df_RNA_classes_sum_counts_all_samples.any()
    if non_zero.any():
        print('true')
    else: 
        print('false')
    # save results
    df_RNA_classes_sum_counts_all_samples.to_csv(organism + "_RNA_classes_sum_counts_" + quantifier + ".tsv", sep = '\t')  
    df_percentage_RNA_all_samples.to_csv(organism + "_RNA_classes_percentage_" + quantifier + ".tsv", sep = '\t')  
       
        
    
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

parser = argparse.ArgumentParser(description="""RNA class statistics""")
parser.add_argument("-q", "--quantification_table", metavar='<quantification_table>', help="Path to the quantification table")
parser.add_argument("-a", "--gene_attribute", help="gene attribute")
parser.add_argument("-annotations", "--gene_annotations", help="gene annotations")
parser.add_argument("-q_tool", "--quantifier", metavar='<quantifier>', help="name of the quantifier")
parser.add_argument("-org", "--organism", help = 'host or pathogen')
parser.add_argument("-rna", "--rna_classes", help = 'group of RNA classes to be replaced', default='')
parser.add_argument("-o", "--output_dir", metavar='<output>', help="output dir",default='.')


args = parser.parse_args()


if args.organism == 'pathogen': 
    #read annotation file 
    gene_types_df = pd.read_csv(args.gene_annotations, sep = '\t', dtype='str',index_col=None)
    #replace pathogen attribute with host attribute which is present in quantification results
    pathogen_attribute = gene_types_df.columns[0]
    gene_types_df.rename(columns = {pathogen_attribute : args.gene_attribute}, inplace = True)
    # create dictionary with information from annotations for each transcript/gene 
    gene_types = gene_types_df.to_dict('records')
    # list of unique gene types present in annotations 
    set_gene_types_pathogen = list(set([gene["gene_type"] for gene in gene_types]))
    # create dictionary {'name of main RNA class': 'corresponding RNA class present in annotations'}, e.g. {'sRNA':'sRNA'}
    dict_set_RNA_classes = create_dictionary_of_RNA_classes_from_set(set_gene_types_pathogen,set_gene_types_pathogen)

elif args.organism == 'host':
    #read annotation file 
    gene_types_df = pd.read_csv(args.gene_annotations, sep = '\t', dtype='str',index_col=None)
    # create dictionary from annotation file
    gene_types = gene_types_df.to_dict('records') 
    # create list of avaiable gene types
    set_gene_types_host = list(set([gene['gene_type'] for gene in gene_types]))   
    #read groups of rna classes
    RNA_classes_to_replce_df = pd.read_csv(args.rna_classes, sep = '\t', dtype='str',index_col=None)
    #create dictionary of RNA classes for each entry (line of the RNA_classes_to_replace.tsv file)
    RNA_classes_to_replce = RNA_classes_to_replce_df.to_dict('records')

    # create dictionary with grouped subclasses of RNA and their main RNA class, e.g. {'CDS': ['IG_LV_gene', 'IG_C_gene'],'pseudogene': ['pseudogene', 'transcribed_unprocessed_pseudogene']}
    grouped_RNA_classes = defaultdict(list) 
    for d in RNA_classes_to_replce:
          for key, value in d.items():
              grouped_RNA_classes[key].append(value)

    #remove nan from values of grouped_RNA_classes dict.
    rna_to_replace = {}
    for k,v in grouped_RNA_classes.items():
          v = [x for x in v if str(x) != 'nan']
          rna_to_replace.update({k:v})

    #remove RNA subclasses defined in RNA_classes_to_replace.tsv from set_gene_types_host before creating final dictionary
    for RNA_class, RNA_class_to_replace in list(rna_to_replace .items()):
          set_gene_types_host = [elem for elem in set_gene_types_host if elem not in RNA_class_to_replace]

    # create dictionary {'main RNA class': 'corresponding RNA class present in annotations'}, e.g. {'sRNA':'sRNA'}
    dict_set_RNA_classes = create_dictionary_of_RNA_classes_from_set(set_gene_types_host,set_gene_types_host)

    # add grouped_RNA_classes dict. into dict_set_RNA_classes dict. 
    dict_set_RNA_classes.update(rna_to_replace)


# calculate RNA class statistics for either host or pathogen
RNA_classes_for_each_sample(args.quantification_table,gene_types,dict_set_RNA_classes,args.organism, args.quantifier, args.gene_attribute)


