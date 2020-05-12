#!/bin/bash

for i in "$@"
do
lines=$(zcat -f $i| wc -l)
count=$(($lines / 4))

name=$(echo "$i" | tr \: _)
echo -n -e "$name\t$count\n"
done

