#!/bin/bash

if [ $3 = "salmon" ]; then
	processed=$(grep "num_processed" $1 | sed 's/num_processed//g'| sed 's/[^a-zA-Z0-9]//g') 
	echo -e "$2\t${processed}" > $2.txt
elif [ $3 = "salmon_alignment" ]; then
	processed=$(grep "num_mapped" meta_info.json | sed 's/num_mapped//g'| sed 's/[^a-zA-Z0-9]//g') 
	echo -e "$2\t${processed}" > $2.txt
elif [ $3 == "star" ]; then
	processed=$(grep "Number of input reads" $1 | sed 's/Number of input reads//g'| sed 's/[^a-zA-Z0-9]//g')
	echo -e "$2\t${processed}" > $2.txt
fi




