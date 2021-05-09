#!/usr/bin/bash

echo "This is script for successfuly installing MIPgen and generating BWA index"

# path_to_store_mipgen"
inpath=$1

cd $inpath

git clone https://github.com/shendurelab/MIPGEN
cd MIPGEN
make

echo "MIPGEN is successfully download in $input folder"


tar -xf mipgen_example.tgz

cd mipgen_practice

## bwa index write here!! 


