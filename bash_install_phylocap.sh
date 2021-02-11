#!/bin/bash

mkdir conda

eval "$(conda shell.bash hook)"

conda create phylocap

conda activate phylocap

conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority flexible

conda install -c conda-forge r-base=3.5 r-devtools r-ape r-stringr r-data.table r-seqinr r-foreach r-doparallel r-rdrop2

conda install -c bioconda bioconductor-rsamtools bioconductor-genomicranges bioconductor-biostrings fastp spades bbmap mafft trnascan-se bwa samtools=1.10 gatk4 trimal cap3 iqtree
