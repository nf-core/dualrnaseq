#!/bin/bash
transcriptome_pathogen=$1
transcriptome_host=$2
quant_table=$3
suffix=$4

grep ">" $transcriptome_pathogen | awk -F ">" '{ print $2 }' | awk 'NR==FNR{a[$0]=$0}NR>FNR{if($1==a[$1])print $0}' - $quant_table > pathogen_quant
awk 'NR==1 {print; exit}' $quant_table | cat - pathogen_quant > pathogen_quant_$suffix

grep ">" $transcriptome_host | awk -F ">" '{ print $2 }' | awk -F ' ' '{print $1}' | awk 'NR==FNR{a[$0]=$0}NR>FNR{if($1==a[$1])print $0}' - $quant_table > host_quant
awk 'NR==1 {print; exit}' $quant_table | cat - host_quant > host_quant_$suffix


