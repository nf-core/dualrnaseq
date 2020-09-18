#!/bin/bash

#-------------------------
#
# Description: 
# Script to count cross mapped reads in *_cross_mapped_reads.txt file created by remove_crossmapped_reads_BAM.sh
#
# Created by B. Mika-Gospodorz
# Date: 27th May 2020
#
# Input files:  list of *_cross_mapped_reads.txt files
# Output: cross_mapped_reads_sum.txt
#
#-------------------------



for i in "$@" # iterate over .txt files
do
lines=$(wc -l $i | awk '{print $1}')  # count number of lines
echo -n -e "$i\t$lines\n" >> cross_mapped_reads_sum.txt  # appends to cross_mapped_reads_sum.txt
done


