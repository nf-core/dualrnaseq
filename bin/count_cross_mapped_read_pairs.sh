#!/bin/bash

for i in "$@"
do
lines=$(wc -l $i | awk '{print $1}')
count_reads=$(($lines / 2))

echo -n -e "$i\t$count_reads\n" >> cross_mapped_reads_sum.txt
done


