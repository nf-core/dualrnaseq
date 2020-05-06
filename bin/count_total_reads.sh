#!/bin/bash

for i in "$@"
do
lines=$(zcat -f $i| wc -l)
count=$(($lines / 4))
echo -n -e "$i\t$count\n"
done

