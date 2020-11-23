#!/bin/bash

#-------------------------
#
# Description: 
# Script to count number of reads in fastq files
#
# Created by B. Mika-Gospodorz
# Date: 2nd April 2020
#
# Input files:  list of fastq files
# Output: name of sample followed by number of reads
#
#-------------------------


for i in "$@" # iterate over fastq files
do
lines=$(zcat -f $i| wc -l) # count lines
count=$(($lines / 4))  # divide number of lines by 4 (fastq file uses four lines per read)
# modify sample name
name=$(echo "$i" | tr -s \: _) 
name2=$(echo "$name" | tr -s - _)
name2=$(echo "${name/.fastq.gz/}")
name2=$(echo "${name2/.fq.gz/}")
name2=$(echo "${name2/.fastq/}")
name2=$(echo "${name2/.fq/}")
name2=$(echo "${name2/_trimmed/}")
echo -n -e "$name2\t$count\n" # return sample name followed by number of reads
done

