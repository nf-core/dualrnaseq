#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 11:26:46 2019

@author: bozena
"""

import argparse
from Bio import SeqIO

def create_transcriptome(fasta_records_dict, gff_files,feature, gene_attribute, output_file_name):
    with open(output_file_name, 'a') as out_name:
        for gff_file in gff_files:  
            for line in open(gff_file):
                 d = line.rstrip()  #remove '\n'
                 if ((d[0] != '#') and (d != '')):
                    d_list = d.split('\t')
                    if d_list[2] in feature:
                        reference_name = d_list[0]
                        split_8 = d_list[8].split(';')
                        #find index of id of interest
                        index_feature = [split_8.index(el) for el in split_8 if gene_attribute in el]
                        description = [s.split("=")[1] for s in split_8]
                        #re-order ids from column 8
                        if not index_feature:
                            print('lack of ' + gene_attribute + ' attribute for gene ID:' + split_8[0].split('=')[1] )
                        else: 
                            #seq_col_8 = list(range(0,len(split_8)))
                            #indexes_without_feature = [x for x in seq_col_8 if x != index_feature[0]]
                            #new_description = ([description[index_feature[0]]] + [description[i] for i in indexes_without_feature])
                            #out_name.write('>' + '|'.join(new_description) + '\n')
                            out_name.write('>' + description[index_feature[0]] + '\n')
                            if d_list[6] == '+' :
                                out_name.write(str(fasta_records_dict[reference_name].seq[int(d_list[3])-1:int(d_list[4])]) + '\n') 
                            elif d_list[6] == '-' :
                                out_name.write(str(fasta_records_dict[reference_name].seq[int(d_list[3])-1:int(d_list[4])].reverse_complement()) + '\n') 
                        

parser = argparse.ArgumentParser()
parser.add_argument("-fasta",nargs='+',help="genome fasta file")
parser.add_argument("-gff", nargs='+', help="gff file")
parser.add_argument("-f", "--gene_feature", nargs='+', help="gene feature defined int the 3rd column of gff file")
parser.add_argument("-a", "--gene_attribute", help="gene attribute")
parser.add_argument("-o", help="output file name")

args = parser.parse_args()

gene_features = [feature.replace('[' , '').replace(']','').replace(',','') for feature in args.gene_feature ]
fasta_records_dict =  dict()
for fasta_file in  args.fasta:
        fasta_records_dict.update(SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta')))
create_transcriptome(fasta_records_dict, args.gff, gene_features, args.gene_attribute, args.o)
    

