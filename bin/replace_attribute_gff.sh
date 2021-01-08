#!/bin/bash

#-------------------------
#
# Description: 
# Script to replace gene attributes in 9th column of gff file, e.g. locus_tag with transcript_id, transcript_id with parent
#
# Created by B. Mika-Gospodorz
# Date: 20th April 2020
#
# Input files:  $1 gff file
# 		$2 output file name
# 		$3 new attribute
# 		$4 current attribute 
# Output: gff file defined in the 2nd argument
#
#-------------------------


input_gff=$1
out_file=$2

# replace current attribute in 9th column of gff file with a new one and create new gff file
awk -v new_atr="$3" -v current_atr="$4" -F '\t' '{gsub(current_atr,new_atr,$9);print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' "${input_gff}" > ${out_file}





