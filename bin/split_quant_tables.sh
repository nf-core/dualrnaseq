#!/bin/bash
quant_table=$1
host_annotations=$2
pathogen_annotations=$3
out_name=$4


host_quant=$(awk -F"\t" '!(NR<=1) {print $1}' $host_annotations| awk 'NR==FNR{a[$0]=$0}NR>FNR{if($1==a[$1])print $0}' - $quant_table) 
awk 'NR==1 {print; exit}' $quant_table | cat - $host_quant > host_$out_name

pathogen_quant=$(awk -F"\t" '!(NR<=1) {print $1}' $pathogen_annotations| awk 'NR==FNR{a[$0]=$0}NR>FNR{if($1==a[$1])print $0}' - $quant_table) 
awk 'NR==1 {print; exit}' $quant_table | cat - $pathogen_quant > pathogen_$out_name
