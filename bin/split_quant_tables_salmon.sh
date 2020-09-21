#!/bin/bash

#-------------------------
#
# Description: 
# Script to split Salmon quantification tables into host and pathogen results
#
# Created by B. Mika-Gospodorz
# Date: 4th June 2020
#
# Input files:  $1 pathogen transcriptome fasta file
# 		$2 host transcriptome fasta file
# 		$3 Salmon quantification table 
# 		$4 output file name
# Output: host_* and pathogen_* tsv files defined in the 4th argument
#
#-------------------------

transcriptome_pathogen=$1
transcriptome_host=$2
quant_table=$3
suffix=$4

# extract pathogen quantification results based on transcript names defined in transcriptome file (extract transcript names from transcriptome data file and extract quantification results for them)
grep ">" $transcriptome_pathogen | awk -F ">" '{ print $2 }' | awk 'NR==FNR{a[$0]=$0}NR>FNR{if($1==a[$1])print $0}' - $quant_table > pathogen_quant
# copy header from quantification table to pathone quantification table
awk 'NR==1 {print; exit}' $quant_table | cat - pathogen_quant > pathogen_$suffix

# extract host quantification results based on transcript names defined in transcriptome file
grep ">" $transcriptome_host | awk -F ">" '{ print $2 }' | awk -F ' ' '{print $1}' | awk 'NR==FNR{a[$0]=$0}NR>FNR{if($1==a[$1])print $0}' - $quant_table > host_quant
# copy header from quantification table to host quantification table
awk 'NR==1 {print; exit}' $quant_table | cat - host_quant > host_$suffix


