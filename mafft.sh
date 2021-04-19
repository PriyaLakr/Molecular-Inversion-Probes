#!/usr/bin/bash

echo "Give name of your fasta file as first argument, path of ouput folder as second argument, and path of input files as third argument. Run -> bash mafft.sh filename outpath inpath" 

name=$(basename $1 .fasta) 
in_data=$2
out="$2/${name}_out"

mkdir $in_data 
cd $in_data 
cp $3/$1 $in_data
mkdir $out



for i in $1; do mafft  --auto --leavegappyregion $i > $out/${name}_aligned.fasta; done

