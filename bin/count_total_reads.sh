#!/bin/bash

for i in "$@"
do
lines=$(zcat -f $i| wc -l)
count=$(($lines / 4))
name=$(echo "$i" | tr \: _)
name2=$(echo "${name/.fastq.gz/}")
name2=$(echo "${name2/.fq.gz/}")
name2=$(echo "${name2/.fastq/}")
name2=$(echo "${name2/.fq/}")
echo -n -e "$name2\t$count\n"
done

