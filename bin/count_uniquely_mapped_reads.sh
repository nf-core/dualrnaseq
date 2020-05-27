#!/bin/bash

bam_file=$1
host_reference_names=$2
pathogen_reference_names=$3
sample_name=$4
out_name=$5


host=$(samtools view -F 4 -h $bam_file | grep -f $host_reference_names | fgrep -w NH:i:1 | echo "$sample_name host `wc -l`")
pathogen=$(samtools view -F 4 -h $bam_file | grep -f $pathogen_reference_names | fgrep -w NH:i:1 | echo "$sample_name pathogen `wc -l`")
counts="${host}\n${pathogen}" 
echo -e $counts > $out_name


