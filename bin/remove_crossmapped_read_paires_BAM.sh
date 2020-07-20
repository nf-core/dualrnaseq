#!/bin/bash
alignment=$1
extract_crossmapped_reads_script_path=$2
host_reference=$3
pathogen_reference=$4
out_name=$5



samtools view -f 0x2 -h $alignment | fgrep -vw NH:i:1 | python $extract_crossmapped_reads_script_path/extract_crossmapped_reads.py -h_ref $host_reference -p_ref $pathogen_reference -o $out_name

cross_mapped_fragments=$out_name
out_bam_name=$6

samtools view -h $alignment | fgrep -wvf $cross_mapped_fragments | samtools view -bS -o $out_bam_name -
