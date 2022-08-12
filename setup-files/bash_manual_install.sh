#!/bin/bash

cd /home/chutter/

eval "$(conda shell.bash hook)"

conda create --name PhyloCap

conda activate PhyloCap

conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority flexible

conda install -c conda-forge r-base=4.0.2 r-devtools r-ape r-stringr r-data.table r-seqinr r-foreach r-doparallel r-rdrop2 r-biomartr

conda install -c bioconda bioconductor-rsamtools bioconductor-genomicranges bioconductor-biostrings fastp spades mafft bwa samtools=1.1.0 gatk4 trimal iqtree
