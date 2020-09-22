#!/bin/bash

#-------------------------
#
# Description: 
# Script to extract information on processed reads by Salmon and STAR
#
# Created by B. Mika-Gospodorz
# Date: 8th April 2020
#
# Input files:  $1 log file which contains information about processed reads by a tool: meta_info.json for Salmon or *Log.final.out for STAR
# 		$2 sample name
# 		$3 name of tool: salmon, salmon_alignment, or star
# Output: {sample_name}.txt file with number of processed reads by tool
#-------------------------

if [ $3 == "salmon" ]; then # for Salmon extract 'num_processed' from meta_info.json file
	processed=$(grep "num_processed" $1 | sed 's/num_processed//g'| sed 's/[^a-zA-Z0-9]//g') 
	echo -e "$2\t${processed}" > $2.txt
elif [ $3 == "salmon_alignment" ]; then # for Salmon alignment-based mode extract 'num_mapped' from meta_info.json file
	processed=$(grep "num_mapped" $1 | sed 's/num_mapped//g'| sed 's/[^a-zA-Z0-9]//g') 
	echo -e "$2\t${processed}" > $2.txt
elif [ $3 == "star" ]; then  # for STAR extract "Number of input reads" from *Log.final.out file
	processed=$(grep "Number of input reads" $1 | sed 's/Number of input reads//g'| sed 's/[^a-zA-Z0-9]//g')
	echo -e "$2\t${processed}" > $2.txt
fi




