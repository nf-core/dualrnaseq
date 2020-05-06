#!/bin/bash
input_gff=$1
out_file=$2

awk -v pathogen="$4" -v host="$3" -F '\t' '{gsub(pathogen,host,$9);print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' "${input_gff}" > ${out_file}


# OFS='\t'
