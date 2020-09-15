#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 10:20:19 2019

@author: bozena
"""
import argparse
import pandas as pd
from collections import defaultdict


# if not os.path.exists("host"):
#     os.makedirs("host")
# if not os.path.exists("pathogen"):
#     os.makedirs("pathogen")
    
       
def create_dictionary_of_RNA_classes_from_set(keys,values):
    zipped = list(zip(keys,values))
    return dict(zipped)
        

####~~~~~~~~~~~~~~~~~~~~~~~ get statistics  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    
def calculate_percentage_RNA_classes(RNA_classes_sum_counts):
    return((RNA_classes_sum_counts/RNA_classes_sum_counts.sum())*100)


def add_counts_to_genes(quantification_df,gene_types,feature):
    for id_gene in gene_types:
        if id_gene[feature] in quantification_df.index:
            id_gene['counts'] = quantification_df.at[id_gene[feature]]
        else:
            id_gene['counts'] = 0 
    return(gene_types)


def get_gene_type_host(gene_RNAtype_counts, set_RNA_classes,gene_attribute):
    for RNA_class_group, RNA_classes in list(set_RNA_classes.items()):
        for gene in gene_RNAtype_counts:
            if gene['gene_type'] in RNA_classes:
                gene['gene_type'] = RNA_class_group
    gene_RNAtype_counts_df = pd.DataFrame.from_dict(gene_RNAtype_counts) 
    gene_RNAtype_counts_df.to_csv("host_gene_types_groups_" + gene_attribute.split('_')[0] + ".csv",index = False)
    
    
def sum_counts_for_each_RNA_class(gene_RNAtype_counts, set_RNA_classes):
    RNA_classes_sum_counts = []
    for RNA_class_group, RNA_classes in list(set_RNA_classes.items()):  
        all_genes_in_RNA_class = ([gene for gene in gene_RNAtype_counts if gene['gene_type'] in RNA_class_group]) 
        class_sum_counts = sum([gene_counts['counts'] for gene_counts in all_genes_in_RNA_class])
        RNA_classes_sum_counts.append({'gene_type': RNA_class_group, 'sum_counts':float(class_sum_counts)})
        RNA_classes_sum_counts_df = pd.DataFrame.from_dict(RNA_classes_sum_counts)
        RNA_classes_sum_counts_df.index = RNA_classes_sum_counts_df['gene_type']
        RNA_classes_sum_counts_df = RNA_classes_sum_counts_df.drop('gene_type', axis=1)
    return(RNA_classes_sum_counts_df)    

        
                                                                          
# def RNA_classes_for_each_sample_htseq(quantification_table_path,gene_types,set_RNA_classes,organism, profile, gene_attribute):
#     col_names = pd.read_csv(quantification_table_path, sep = '\t', nrows=0).columns
#     types_dict = {gene_attribute: str}
#     types_dict.update({col: float for col in col_names if col not in types_dict})  
#     quantification_table = pd.read_csv(quantification_table_path, sep = '\t',index_col=0,dtype=types_dict)
#     df_percentage_RNA_all_samples  = pd.DataFrame(columns = set_RNA_classes,index=['NA'])
#     df_RNA_classes_sum_counts_all_samples  = pd.DataFrame(columns = set_RNA_classes,index=['NA'])
#     for sample_name in quantification_table:
#         gene_RNA_type_counts = add_counts_to_genes(quantification_table[sample_name],gene_types,gene_attribute)
#         # get gene types for host
#         if organism == 'host':
#             get_gene_type_host(gene_RNA_type_counts,set_RNA_classes,gene_attribute)
#         sum_counts_RNA_class_df = sum_counts_for_each_RNA_class(gene_RNA_type_counts,set_RNA_classes)
#         #sample name as an index
#         df_RNA_classes_sum_counts = sum_counts_RNA_class_df.sort_values('sum_counts', ascending=0)
#         df_RNA_classes_sum_counts = pd.DataFrame(df_RNA_classes_sum_counts).T        
#         df_RNA_classes_sum_counts.index = [sample_name]
#         df_RNA_classes_sum_counts = df_RNA_classes_sum_counts.T
#         percentage_RNA_class = calculate_percentage_RNA_classes(df_RNA_classes_sum_counts)
#         df_RNA_classes_sum_counts_all_samples = pd.concat([df_RNA_classes_sum_counts_all_samples, df_RNA_classes_sum_counts.T],sort=True)  
#         df_percentage_RNA_all_samples = pd.concat([df_percentage_RNA_all_samples, percentage_RNA_class.T],sort=True)
#     df_percentage_RNA_all_samples = df_percentage_RNA_all_samples.dropna()
#     df_RNA_classes_sum_counts_all_samples = df_RNA_classes_sum_counts_all_samples.dropna()
#     df_RNA_classes_sum_counts_all_samples.to_csv(organism + "_RNA_classes_sum_counts_" + profile + ".csv", sep = '\t')  
#     df_percentage_RNA_all_samples.to_csv(organism + "_RNA_classes_percentage_" + profile + ".csv", sep = '\t')  
        
    
      
def RNA_classes_for_each_sample(quantification_table_path,gene_types,set_RNA_classes,organism, profile, gene_attribute):
    col_names = pd.read_csv(quantification_table_path, sep = '\t', nrows=0).columns
    types_dict = {gene_attribute: str}
    types_dict.update({col: float for col in col_names if col not in types_dict})
    quantification_table = pd.read_csv(quantification_table_path, sep = '\t',index_col=0, dtype=types_dict)
    quantification_table = quantification_table.fillna(0)
    df_percentage_RNA_all_samples  = pd.DataFrame(columns = set_RNA_classes,index=['NA'])
    df_RNA_classes_sum_counts_all_samples  = pd.DataFrame(columns = set_RNA_classes,index=['NA'])
    for sample_name in quantification_table:
        if 'NumReads' in sample_name:
            gene_RNA_type_counts = add_counts_to_genes(quantification_table[sample_name],gene_types,gene_attribute)
            # get gene types for host
            if organism == 'host':
                get_gene_type_host(gene_RNA_type_counts,set_RNA_classes,gene_attribute)  
            sum_counts_RNA_class_df = sum_counts_for_each_RNA_class(gene_RNA_type_counts,set_RNA_classes)
            #sample name as an index
            df_RNA_classes_sum_counts = sum_counts_RNA_class_df.sort_values('sum_counts', ascending=0)
            df_RNA_classes_sum_counts = pd.DataFrame(df_RNA_classes_sum_counts).T        
            df_RNA_classes_sum_counts.index = [sample_name]
            df_RNA_classes_sum_counts = df_RNA_classes_sum_counts.T
            percentage_RNA_class = calculate_percentage_RNA_classes(df_RNA_classes_sum_counts)
            df_RNA_classes_sum_counts_all_samples = pd.concat([df_RNA_classes_sum_counts_all_samples, df_RNA_classes_sum_counts.T],sort=True)  
            df_percentage_RNA_all_samples = pd.concat([df_percentage_RNA_all_samples, percentage_RNA_class.T],sort=True)
    df_percentage_RNA_all_samples = df_percentage_RNA_all_samples.dropna()
    df_RNA_classes_sum_counts_all_samples = df_RNA_classes_sum_counts_all_samples.dropna()
    non_zero = df_RNA_classes_sum_counts_all_samples.any()
    if non_zero.any():
        print(True)
    else: 
        print(False) 
    df_RNA_classes_sum_counts_all_samples.to_csv(organism + "_RNA_classes_sum_counts_" + profile + ".csv", sep = '\t')  
    df_percentage_RNA_all_samples.to_csv(organism + "_RNA_classes_percentage_" + profile + ".csv", sep = '\t')  
       
        
    
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

parser = argparse.ArgumentParser(description="""RNA class statistics""")
parser.add_argument("-q", "--quantification_table", metavar='<quantification_table_host>', help="Path to the quantification table")
parser.add_argument("-a", "--gene_attribute", help="gene attribute")
parser.add_argument("-annotations", "--gene_annotations", help="gene annotations")
parser.add_argument("-q_tool", "--quantifier", metavar='<quantifier>', help="name of the quantifier")
parser.add_argument("-org", "--organism", help = 'host or pathogen')
parser.add_argument("-rna", "--rna_classes", help = 'group of RNA classes to be replaced', default='')
parser.add_argument("-o", "--output_dir", metavar='<output>', help="output dir",default='.')
parser.add_argument("-p", "--profile", default='', metavar='<profile_of_mapping>', help="multi/uniquely -mapped")


args = parser.parse_args()


if args.organism == 'pathogen': 
    
    gene_types_df = pd.read_csv(args.gene_annotations, sep = '\t', dtype='str',index_col=None)
    #replace pathogen attribute with host attribute which is present in the quantification results
    pathogen_attribute = gene_types_df.columns[0]
    gene_types_df.rename(columns = {pathogen_attribute : args.gene_attribute}, inplace = True)
    gene_types = gene_types_df.to_dict('records')
    set_gene_types_pathogen = list(set([gene["gene_type"] for gene in gene_types]))
    dict_set_RNA_classes = create_dictionary_of_RNA_classes_from_set(set_gene_types_pathogen,set_gene_types_pathogen)


elif args.organism == 'host':
    #dictionary of annotations
    
    gene_types_df = pd.read_csv(args.gene_annotations, sep = '\t', dtype='str',index_col=None)
    gene_types = gene_types_df.to_dict('records')

    #create a list of avaiable gene types
    set_gene_types_host = list(set([gene['gene_type'] for gene in gene_types]))   

    #read groups of rna class
    RNA_classes_to_replce_df = pd.read_csv(args.rna_classes, sep = '\t', dtype='str',index_col=None)
    RNA_classes_to_replce = RNA_classes_to_replce_df.to_dict('records')

    # create a dictionary of rna classes to be replaced 
    dd = defaultdict(list)
    for d in RNA_classes_to_replce:
          for key, value in d.items():
              dd[key].append(value)

    #remove nan
    rna_to_replace = {}
    for k,v in dd.items():
          v = [x for x in v if str(x) != 'nan']
          rna_to_replace.update({k:v})

    #remove RNA classes to replace from set_gene_types_host before creating a dictionary
    for RNA_class, RNA_class_to_replace in list(rna_to_replace .items()):
          set_gene_types_host = [elem for elem in set_gene_types_host if elem not in RNA_class_to_replace]
    dict_set_RNA_classes = create_dictionary_of_RNA_classes_from_set(set_gene_types_host,set_gene_types_host)
    dict_set_RNA_classes.update(rna_to_replace)


#if args.quantifier == 'htseq':
RNA_classes_for_each_sample(args.quantification_table,gene_types,dict_set_RNA_classes,args.organism, args.profile, args.gene_attribute)
#elif args.quantifier == 'salmon':
#RNA_classes_for_each_sample(args.quantification_table,gene_types,dict_set_RNA_classes,args.organism, args.profile, args.gene_attribute)

