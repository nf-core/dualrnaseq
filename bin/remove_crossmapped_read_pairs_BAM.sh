#!/bin/bash

#-------------------------
#
# Description: 
# Script to remove cross-mapped reads from bam file 
#
# Created by B. Mika-Gospodorz
# Date: 17th July 2020
#
# Input files:  $1 bam file 
#		$2 directory where the extract_crossmapped_reads.py script is stored
# 		$3 reference_host_names.txt file that contains host references extracted with extract_reference_names_from_fasta_files.sh
# 		$4 reference_pathogen_names.txt file that contains pathogen references extracted with extract_reference_names_from_fasta_files.sh
# 		$5 output file name for list of cross-mapped reads, e.g. *_cross_mapped_reads.txt
# 		$6 output file name for bam file with removed cross-mapped reads, e.g. *_no_crossmapped.bam
# Output: *_cross_mapped_reads.txt, *_no_crossmapped.bam
#
#-------------------------


alignment=$1
extract_crossmapped_reads_script_path=$2
host_reference=$3
pathogen_reference=$4
out_name=$5

# extract multi-mapped reads from bam file (-f 0x2 - read mapped in proper pair) and run extract_crossmapped_reads.py to extract reads that mapped onto both host and pathogen genomes simultaneously
samtools view -f 0x2 -h $alignment | fgrep -vw NH:i:1 | python $extract_crossmapped_reads_script_path/extract_crossmapped_reads.py -h_ref $host_reference -p_ref $pathogen_reference -o $out_name

cross_mapped_fragments=$out_name
out_bam_name=$6

# remove cross-mapped reads from bam file and create a new one
samtools view -h $alignment | fgrep -wvf $cross_mapped_fragments | samtools view -bS -o $out_bam_name -
