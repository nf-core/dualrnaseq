#!/bin/bash

#-------------------------
#
# Description: 
# Script to count multi mapped reads (paired-end reads) from .bam file
#
# Created by B. Mika-Gospodorz
# Date: 17th July 2020
#
# Input files:  $1 bam file 
# 		$2 reference_host_names.txt file that contains host references extracted with extract_reference_names_from_fasta_files.sh
# 		$3 reference_pathogen_names.txt file that contains pathogen references extracted with extract_reference_names_from_fasta_files.sh
# 		$4 sample name
# 		$5 output file name *_multi_mapped.txt
# Output: *_uniquely_mapped.txt file with list of uniquely mapped reads
#-------------------------

bam_file=$1
host_reference_names=$2
pathogen_reference_names=$3
sample_name=$4
out_name=$5

# count reads that mapped uniquely onto host genome (-f 66 to extract read mapped in proper pair and first in pair) 
host=$(samtools view -f 66 -h $bam_file | grep -f $host_reference_names | fgrep -w NH:i:1 | echo "$sample_name host `wc -l`")
# count reads that mapped uniquely onto pathogen genome 
pathogen=$(samtools view -f 66 -h $bam_file | grep -f $pathogen_reference_names | fgrep -w NH:i:1 | echo "$sample_name pathogen `wc -l`")
# combine information about host and pathogen uniquely mapped reads
counts="${host}\n${pathogen}" 
# save into file
echo -e $counts > $out_name
