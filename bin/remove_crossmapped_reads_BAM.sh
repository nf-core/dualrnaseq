#!/bin/bash

alignment=$1
echo $alignment

cross_mapped_reads=$2
echo $cross_mapped_reads

out_name=$3

echo $out_name
samtools view -h $alignment | fgrep -wvf $cross_mapped_reads | samtools view -bS -o $out_name -
