# Molecular-Inversion-Probes


# Designing MIPs for human genes

MIPgen pipeline by Shendurlab 

Dependencies:
Bwa 0.7.17
samtools 1.11


1. Obtain MIPgen

  git clone https://github.com/shendurelab/MIPGEN
  cd MIPGEN
  make

  tar -xf mipgen_example.tgz

  cd mipgen_practice

  mkdir index

2. Generate index from Human reference file 

  bwa index -p GRCh38.primary_assembly.genome.fa  /data1/priya/new_data/reference/bwa/GRCh38.primary_assembly.genome.fa  

Keep the name of index and fasta file same

3. Running MIPgen

  /data1/priya/mipgen/MIPGEN/mipgen -regions_to_scan /data1/priya/mipgen/MIPGEN/mipgen_practice/practice_genes.bed -project_name /data1/priya/mipgen/MIPGEN/tools/prac_design -min_capture_size 162 -max_capture_size 162 -bwa_genome_index /data1/priya/new_data/reference/GRCh38.primary_assembly.genome.fa



# Designing MIPs for non-human genes

This workflow will generate molecular inversion probes for non-human genes, for instance, bacterial, viral and fungal genes

Requirements:

Python
NCBI-BLAST+
Primer3 

Scripts are added in the scripts folder

