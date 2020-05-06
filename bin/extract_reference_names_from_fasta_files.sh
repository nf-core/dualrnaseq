#!/bin/bash
output_name=$1
shift

for file in "$@"
do
   grep ">" ${file} | awk -F" " '{print $1}' | sed 's/>//' >> ${output_name}

done
