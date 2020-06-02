#!/bin/bash

alignment=$1
cross_mapped_reads=$2
out_name=$3

samtools view -h $alignment | fgrep -wvf $cross_mapped_reads | samtools view -bS -o $out_name -
