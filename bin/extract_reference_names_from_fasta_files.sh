#!/bin/bash

#-------------------------
#
# Description: 
# Script to extract chromosome and plasmid names from genome fasta files
#
# Created by B. Mika-Gospodorz
# Date: 5th September 2019
#
# Input files:  $1 name of output file, eg. reference_host_names.txt
# 		$@ list of fasta files
# Output: file with a given name defined in the first argument that contains names of chromosomes or plasmids defined in the fasta files
#
#-------------------------

# assign output file name to variable 
output_name=$1
shift

for file in "$@" #iterate over fasta files
do   
   grep ">" ${file} | awk -F" " '{print $1}' | sed 's/>//' >> ${output_name}  # extract name followed '>' character and append to output file 
done
