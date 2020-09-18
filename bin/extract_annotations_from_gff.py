#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 16:11:46 2020

@author: B. Mika-Gospodorz and R. Hayward

Input file: gff file
Output files: tsv file
Description: Used to extract annotations from gff file for each gene/transcript (transcript_id, transcript_name, gene_id, gene_name, gene_type) 
"""

import argparse
import pandas as pd
import os


# function to extract annotations from host gff file
def extract_gene_types_host(gff,gene_feature,gene_attribute):
    gene_type_host = []  # initialize list of dictionaries with genes/transcripts and their annotations
    for line in open(gff): #looping through GFF file
        d = line.rstrip()  #remove '\n'
        if ((d[0] != '#') and (d != '')): # skip both commented and empty lines
            d = d.split('\t') # separating tabbed words into a list
            if d[2] in gene_feature: # check if gene feature (3rd column of gff file) is present in list of features used in quantification
                      d1 = d[8].split(';') # split column 8 by ;
                      # extract ID of gene attribute
                      ID_pos_host = [pos.split('=') for pos in d1 if pos.startswith(gene_attribute)][0][1]
                      # extract gene type 
                      gt_host = [pos.split('=') for pos in d1 if pos.startswith('gene_type')][0][1]
                      if d[0] == 'chrM':    # if gene is located on Mitochondrial chromosome, set gene type as mitoRNA
                          g_type = 'mitoRNA'
                      elif gt_host.startswith('Pseudo_'): # if gene type contains Pseudo_ prefix, set gene type as pseudogene 
                          g_type = 'pseudogene'
                      elif 'tRNA' in gt_host:  # if gene type contains tRNA, set gene type as tRNA 
                          g_type = 'tRNA' 
                      else:
                          g_type = gt_host # set original gene type
                      # extract gene_id
                      ID_gene_pos_host = [pos.split(' ') for pos in d1 if pos.startswith('gene_id')][0][0].split("=")[1]
                      # extract gene_name
                      gene_name_host = [pos.split(' ') for pos in d1 if pos.startswith('gene_name')][0][0].split("=")[1]
                      # if quantification is performed at gene-level, store results for genes (in case of HTSeq, attribute is set to 'gene_id')
                      if gene_attribute.startswith('gene'): 
                          gene_type_host.append({'gene_id':ID_gene_pos_host,'gene_name':gene_name_host, 'gene_type':g_type})   
                      # store rusults for transcripts - in case of Salmon, attribute is set to 'parent'
                      else:  
                          # extract transcript_name
                          transcript_name_host = [pos.split('=') for pos in d1 if pos.startswith('transcript_name')][0][1]
                          gene_type_host.append({'transcript_id' :ID_pos_host, 'transcript_name':transcript_name_host, 'gene_id':ID_gene_pos_host,'gene_name':gene_name_host, 'gene_type':g_type})                    
    return gene_type_host




# function to extract annotations from pathogen gff file
def extract_gene_types_pathogen(gff,gene_feature,gene_attribute):
    gene_type_pathogen = [] #create new list
    for line in open(gff): #looping through GFF file
        d = line.rstrip()  #remove '\n'
        if (not d.startswith('#') and (d != '')): #ignoring comments and blank lines
              d = d.split('\t') #separating tabbed words into a list
              if d[2] in gene_feature: #if values from 3rd col of GFF are in gene_feature
                d1 = d[8].split(';') #split column 8 by ;
                
                #Error handler to ignore rows that don't contain the same format as other rows
				#This is a common issue in bacterial GTFs/GFFs which are often composed with non-uniform rows
                try:
                    #further split contents from col 8 by '=', taking first occurance and 2nd value. 
					#ie, ID=1234 becomes 'ID', '1234', and '1234' is stored
                    ID_pos_pathogen = [pos.split('=') for pos in d1 if pos.startswith(gene_attribute)][0][1]
                except Exception:
                    continue
                    
                #Search for field 'Name/name' and store contents
                if ('name' or 'Name' in d1):
                    feature_name = [pos.split('=') for pos in d1 if pos.startswith(tuple(['name','Name']))]
                    if feature_name:
                        feature_name_pathogen = feature_name[0][1]
                    else:
                        feature_name_pathogen = ''
                else:
                    feature_name_pathogen = ''
                
                #Capture biotypes
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
    


parser = argparse.ArgumentParser(description="""extract annotations from gff""")
parser.add_argument("-gff", "--gff", metavar='<gff_annotations>', help="Path to gff file")
parser.add_argument("-f", "--gene_feature", nargs='+', help="gene features from 3rd column of gff file used for quantification")
parser.add_argument("-a", "--gene_attribute", help="gene attribute from 9th column of gff file")
parser.add_argument("-q_tool", "--quantifier", metavar='<quantifier>', help="name of quantifier")
parser.add_argument("-org", "--organism", help = 'host or pathogen')
parser.add_argument("-o", "--output", metavar='<output>', help="output file name")

args = parser.parse_args()

# create list of gene features used in quantification
gene_features = [feature.replace('[' , '').replace(']','').replace(',','')  for feature in args.gene_feature]

if args.organism == 'pathogen': 
      # create dictionary of annotations for each gff entry with desired gene feature
      gene_types = extract_gene_types_pathogen(args.gff,gene_features, args.gene_attribute)
      # create data frame of quantified genes/transcripts with annotations (e.g. gene type and gene name)
      gene_annotations_pathogen_df = pd.DataFrame.from_dict(gene_types)
      # check if data frame returns any values  
      if gene_annotations_pathogen_df.empty:
        print('No features matched the input criteria of: ', gene_features, ' and ', args.gene_attribute)
      else:
      # save results
        gene_annotations_pathogen_df.to_csv(args.output + '_' + args.gene_attribute + '_' + args.quantifier + ".tsv",index =False, sep = '\t')     
        
elif args.organism == 'host':
    # dictionary of annotations for gff entries with desired gene feature
    gene_types = extract_gene_types_host(args.gff,gene_features, args.gene_attribute)
    # remove duplicates that occur in case of exons. Exons of the same isoforms have the same 'transcript_id'.
    gene_types_unique = pd.DataFrame(gene_types).drop_duplicates().to_dict('records')
    # create data frame of quantified genes/transcripts with annotations and save it
    gene_annotations_host_df = pd.DataFrame.from_dict(gene_types_unique)
    gene_annotations_host_df.to_csv(args.output + '_annotations_' + args.quantifier + ".tsv",index=False, sep = '\t')

