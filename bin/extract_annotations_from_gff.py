#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 16:11:46 2020

@author: bozena
"""

import argparse
import pandas as pd
import os

if not os.path.exists("host"):
    os.makedirs("host")
if not os.path.exists("pathogen"):
    os.makedirs("pathogen")
      


####~~~~~~~~~~~~~~~~~~~~~~~ Extract annotations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
def extract_gene_types_host(gff,gene_feature,feature):
    gene_type_host = []
    r_pos_host = []
    for line in open(gff):
        d = line.rstrip()  #remove '\n'
        if ((d[0] != '#') and (d != '')):
            d = d.split('\t')
            if d[2] in gene_feature:
                      d1 = d[8].split(';')
                      #in case of tRNAs transcript_id was replaced with 'Parent' in gff
                      ID_pos_host = [pos.split('=') for pos in d1 if pos.startswith(feature)][0][1]
                      transcript_name_host = [pos.split('=') for pos in d1 if pos.startswith('transcript_name')][0][1]
                      r_pos_host = [pos.split('=') for pos in d1 if pos.startswith('gene_type')][0][1]
                      if d[0] == 'chrM':
                          g_type = 'mitoRNA'
                      elif r_pos_host.startswith('Pseudo_'):
                          g_type = 'pseudogene'
                      elif 'tRNA' in r_pos_host:
                          g_type = 'tRNA'
                      else:
                          g_type = r_pos_host     
                      # if feature.startswith('gene'):
                      #    gene_type_host.append({feature + '_ID':ID_pos_host, 'gene_type':g_type, feature + '_name':feature_name_host})   
                      # else:
                      ID_gene_pos_host = [pos.split(' ') for pos in d1 if pos.startswith('gene_id')][0][0].split("=")[1]
                      gene_name_host = [pos.split(' ') for pos in d1 if pos.startswith('gene_name')][0][0].split("=")[1]
                      gene_type_host.append({'transcript_id' :ID_pos_host, 'transcript_name':transcript_name_host, 'gene_ID':ID_gene_pos_host,'gene_name':gene_name_host, 'gene_type':g_type, })                    
    return gene_type_host





def extract_gene_types_pathogen(gff,gene_feature,feature):
    gene_type_pathogen = []
    for line in open(gff):
        d = line.rstrip()  #remove '\n'
        if (not d.startswith('#') and (d != '')):
              d = d.split('\t')
              if d[2] in gene_feature:
                d1 = d[8].split(';')
                ID_pos_pathogen = [pos.split('=') for pos in d1 if pos.startswith(feature)][0][1]
                if ('name' or 'Name' in d1):
                    feature_name = [pos.split('=') for pos in d1 if pos.startswith(tuple(['name','Name']))]
                    if feature_name:
                        feature_name_pathogen = feature_name[0][1]
                    else:
                        feature_name_pathogen = ''
                else:
                    feature_name_pathogen = ''
                     
                if d[2] == 'sRNA':
                    g_type = 'sRNA'
                elif d[2] == 'ncRNA':
                    g_type = 'ncRNA'
                elif d[2] == 'tRNA':
                    g_type = 'tRNA'
                elif d[2] == 'rRNA':
                    g_type = 'rRNA'
                elif ('gene_biotype' in d1):
                    r_pos_host = [pos.split('=') for pos in d1 if pos.startswith('gene_biotype')]
                    g_type = r_pos_host[0][1]
                else:
                    g_type = d[2]
                gene_type_pathogen.append({feature: ID_pos_pathogen, 'name':feature_name_pathogen,'gene_type':g_type})   
    return gene_type_pathogen
    



parser = argparse.ArgumentParser(description="""RNA class statistics""")
#parser.add_argument("-q", "--quantification_table", metavar='<quantification_table_host>', help="Path to the quantification table")
parser.add_argument("-gff", "--gff", metavar='<gtf_host_genome>', help="Path to the host reference gtf, multiple gtfs (read as a list)")
parser.add_argument("-f", "--gene_feature", nargs='+', help="gene feature defined int the 3rd column of gff file")
parser.add_argument("-a", "--gene_attribute", help="gene attribute")
parser.add_argument("-q_tool", "--quantifier", metavar='<quantifier>', help="name of the quantifier")
parser.add_argument("-org", "--organism", help = 'host or pathogen')
parser.add_argument("-o", "--output", metavar='<output>', help="output name")
#parser.add_argument("-p", "--profile", default='', metavar='<profile_of_mapping>', help="multi/uniquely -mapped")

args = parser.parse_args()


gene_features = [feature.replace('[' , '').replace(']','').replace(',','')  for feature in args.gene_feature]
if args.organism == 'pathogen': 
      gene_types = extract_gene_types_pathogen(args.gff,gene_features, args.gene_attribute)
      gene_annotations_pathogen_df = pd.DataFrame.from_dict(gene_types)
      gene_annotations_pathogen_df.to_csv(args.output + '_' + args.gene_attribute + '_' + args.quantifier + ".csv",index =False, sep = '\t')     
elif args.organism == 'host':
    #dictionary of annotations
    gene_types = extract_gene_types_host(args.gff,gene_features, args.gene_attribute)
    #there are duplicates for exons
    gene_types_unique = pd.DataFrame(gene_types).drop_duplicates().to_dict('records')
    gene_annotations_pathogen_df = pd.DataFrame.from_dict(gene_types_unique)
    gene_annotations_pathogen_df.to_csv(args.output + '_' + args.gene_attribute + '_' + args.quantifier + ".csv",index=False, sep = '\t')
    

    
