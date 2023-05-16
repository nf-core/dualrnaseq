#!/bin/bash

#-------------------------
#
# Description:
# Script to replace feature in 3rd column of gff file with "quant"
#
# Created by B. Mika-Gospodorz
# Date: 19th March 2020
#
# Input files:  $1 gff file
# 		$2 output file name
# 		$3 list of features to replace with "quant", e.g. ["exon", "tRNA"]
# Output: gff file defined in the 2nd argument
#
#-------------------------

input_gff=$1
out_file=$2
shift
shift

# read features
feature2=`echo "$@" | tr -d [`
feature2=`echo "$feature2" | tr -d ]`
feature2=`echo "$feature2" | tr -d ,`

# create regular expression for list of features, e.g. exon|tRNA
pat=`echo "$feature2" | tr " " "|"`

# replace given features in the 3rd column of gff file with "quant" and create new gff file
awk -v orig="${pat}" -F '\t' '{gsub(orig,"quant",$3);print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' "${input_gff}" > ${out_file}



