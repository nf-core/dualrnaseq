#!/bin/bash

for i in "$@"
do
lines=$(wc -l $i | awk '{print $1}')

echo -n -e "$i\t$lines\n" >> cross_mapped_reads_sum.txt
done


