#!/usr/bin/bash

echo "Give the name of your gene of interest as first argument"
echo "Give the name of the fasta file containing sequences as second argument"

"""Extract fasta sequences with given IDs/names
^-- get rid of line separator
"""

grep -A1 -i "$1" $2 | grep -v "^--" > $1.fasta 



## --no-group-separator should be used instead of -v but  --no-group-separator is not working in mac grep
## grep -w -A1 --no-group-separator -i "$1" $2  > $1.fasta 
