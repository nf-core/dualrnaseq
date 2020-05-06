#!/bin/bash
input_gff=$1
out_file=$2
shift
shift

feature2=`echo "$@" | tr --delete [`
feature2=`echo "$feature2" | tr --delete ]`
feature2=`echo "$feature2" | tr --delete ,`
#echo $feature2

pat=`echo "$feature2" | tr " " "|"`

echo $pat

awk -v orig="${pat}" -F '\t' '{gsub(orig,"quant",$3);print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' "${input_gff}" > ${out_file}



