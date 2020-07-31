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
                      ID_gene_pos_host = [pos.split(' ') for pos in d1 if pos.startswith('gene_id')][0][0].split("=")[1]
                      gene_name_host = [pos.split(' ') for pos in d1 if pos.startswith('gene_name')][0][0].split("=")[1]
                      if feature.startswith('gene'):
                          gene_type_host.append({'gene_id':ID_gene_pos_host,'gene_name':gene_name_host, 'gene_type':g_type})   
                      else:
                          gene_type_host.append({'transcript_id' :ID_pos_host, 'transcript_name':transcript_name_host, 'gene_id':ID_gene_pos_host,'gene_name':gene_name_host, 'gene_type':g_type})                    
    return gene_type_host





def extract_gene_types_pathogen(gff,gene_feature,feature):
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
                    ID_pos_pathogen = [pos.split('=') for pos in d1 if pos.startswith(feature)][0][1]
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
      #store identified gene types
      gene_annotations_pathogen_df = pd.DataFrame.from_dict(gene_types)
      #Check if df returned any values  
      if gene_annotations_pathogen_df.empty:
        print('No features matched the input criteria of: ', gene_features, ' and ', args.gene_attribute)
      else:
        gene_annotations_pathogen_df.to_csv(args.output + '_' + args.gene_attribute + '_' + args.quantifier + ".csv",index =False, sep = '\t')     
        
elif args.organism == 'host':
    #dictionary of annotations
    gene_types = extract_gene_types_host(args.gff,gene_features, args.gene_attribute)
    #there are duplicates for exons
    gene_types_unique = pd.DataFrame(gene_types).drop_duplicates().to_dict('records')
    gene_annotations_pathogen_df = pd.DataFrame.from_dict(gene_types_unique)
    gene_annotations_pathogen_df.to_csv(args.output + '_annotations_' + args.quantifier + ".csv",index=False, sep = '\t')
    

    
