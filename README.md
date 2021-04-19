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


2. Generate index from Human reference genome file (GRCh38.primary_assembly.genome.fa)
        
        mkdir index
        bwa index -p GRCh38.primary_assembly.genome.fa  /data1/priya/new_data/reference/bwa/GRCh38.primary_assembly.genome.fa  

Keep the name of index and fasta file same

3. Running MIPgen
        
        input="location where MIPGEN is stored"
        cd $input
        ./mipgen -regions_to_scan /data1/priya/mipgen/MIPGEN/mipgen_practice/practice_genes.bed -project_name /data1/priya/mipgen/MIPGEN/tools/prac_design -min_capture_size 162 -max_capture_size 162 -bwa_genome_index /data1/priya/new_data/reference/GRCh38.primary_assembly.genome.fa



# Designing MIPs for non-human genes (Under development)
 

This workflow will generate molecular inversion probes for non-human genes, for instance, bacterial, viral and fungal genes

Requirements:

Python
MAFFT

Scripts are located in pl_path_mips.py

Sequences for AMR genes were obtained from https://card.mcmaster.ca/ontology/ database
Sequences were aligned using MAFFT.sh script
Multiple sequence alignment was used to identify regions containing extensive gaps, and such regions were avoided in downstream analysis.


Design of molecular inversion probes (MIPs).
The general format for the molecular inversion probes (MIPs) used here is a common NNN bp linker flanked by an extension arm of 18 to 20 bp and a ligation arm of 20 to 24 bp.
The unique arms of each MIP target a specific 100 bp genomic region. 
Linker can be the NGS sequencing adaptors and sample-specific barcodes. 
MIPs for a given target region were chosen iteratively from the 5' to 3' end of the target region. Each successive MIP is chosen to satisfy the following criteria in order of priority: 1) resides on the opposite strand as the previous MIP, 2) doesn’t overlap with previous MIP, 3) avoids gapped regions as obtained from MSA. 
Good probe features: lack polyNs (6) GC : 40-70 tm >= 60
If no suitable MIPs can be picked, 2nd condition could be avoided.

run —help flag for inputs in pl_path_mips.py

