#!/bin/bash

processed=$(grep "num_processed" $1 | sed 's/num_processed//g'| sed 's/[^a-zA-Z0-9]//g') 

echo -e "$2\t${processed}" > $2.txt

