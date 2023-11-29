#!/bin/bash

#-------------------------
#
# Description: 
# Script to split HTSeq quantification tables into host and pathogen results
#
# Created by B. Mika-Gospodorz
# Date: 3rd June 2020
#
# Input files:  $1 HTSeq quantification table 
# 		$2 tsv file that contains host annotations extracted with extract_annotations_from_gff.py
# 		$3 tsv file that contains pathogen annotations extracted with extract_annotations_from_gff.py
# 		$4 output file name
# Output: host_* and pathogen_* tsv files defined in the 4th argument
#
#-------------------------

quant_table=$1
host_annotations=$2
pathogen_annotations=$3
out_name=$4

# extract host quantification results based on gene annotation from tsv file (extract first column from host annotation table with gene names and extract quantification results for them)
awk -F"\t" '!(NR<=1) {print $1}' $host_annotations | awk 'NR==FNR{a[$0]=$0}NR>FNR{if($1==a[$1])print $0}' - $quant_table > host_quant
# copy header from quantification table to host quantification table
awk 'NR==1 {print; exit}' $quant_table | cat - host_quant > host_$out_name


# extract pathogen quantification results based on gene annotation from tsv file
awk -F"\t" '!(NR<=1) {print $1}' $pathogen_annotations | awk 'NR==FNR{a[$0]=$0}NR>FNR{if($1==a[$1])print $0}' - $quant_table > pathogen_quant
# copy header from quantification table to pathogen quantification table
awk 'NR==1 {print; exit}' $quant_table | cat - pathogen_quant > pathogen_$out_name