# PHYLOCAP

R package For processing high-throughput sequencing data from targeted sequence capture.


# Installation of R package

PHYLOCAP is an R package that has been tested on R version 3.5. To install PHYLOCAP, you can use the R package devtools. Here are step-by-step instructions for installation:

1) Install devtools by typing in your R console: install.packages("devtools")

2) Install PHYLOCAP by typing in your R console: devtools::install_github("chutter/PHYLOCAP", update = "never")

3) Devtools will ask you to install the package dependecies, select "Yes". If devtools asks you to update packages, you may choose to do so. I would recommend not to install packages from source if devtools asks you and to include the update parameter set to "never" in install_github. 

4) Devtools should finish and say the package loaded properly. Load the package with library(PHYLOCAP). 

And installation should be done! 


# Installation of prerequisites 

PhyloCap uses several R packages and other outside programs for certain functions. 

1. R packages

From CRAN
- devtools
- ape
- stringr
- data.table
- seqinr
- foreach
- doparallel
- rdrop2

From BioConductor
- rsamtools
- genomicranges
- biostrings

2. Outside programs

- fastp: adaptor trimming and paired-end read merging
- bwa: read mapping
- spades: assembly
- BLAST: matching assembled contigs to targets, other utilities
- mafft: creating alignments
- trimal: trimming alignments
- IqTree: gene tree and concatenation trees
- GATK4: variant calling functions
- SamTools: variant calling and read mapping tools

The R packages and outside programs can be installed manually. The easiest and quickest way is to use the Anaconda package manager. Anaconda can be downloaded and installed for different operating systems from https://anaconda.org. Miniconda is recommended. Once a version of Anaconda is installed or loaded, a new clean environment should be created for PhyloCap. Anaconda can be set up and configured through the following steps: 

1. Create new environment:

```
conda create phylocap

conda activate phylocap
```


2. Next, you should configure anaconda by adding and ordering channels as follows:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority flexible
```


3. Finally, the required programs can be installed in two commands. The programs should be installed together so that the correct package dependency versions can be sorted out. 

```
conda install -c conda-forge r-base=3.5 r-devtools r-ape r-stringr r-data.table r-seqinr r-foreach r-doparallel r-rdrop2

conda install -c bioconda bioconductor-rsamtools bioconductor-genomicranges bioconductor-biostrings fastp spades mafft bwa samtools=1.1.0 gatk4 trimal iqtree
```

All the functions for PhyloCap should be ready to go! 

< coming soon a function to test if they can found >


# PhyloCap pipeline tutorials 

[Installation: detailed installation instructions and trouble-shooting ](https://github.com/chutter/PhyloCap/wiki/Installation)

[Tutorial 1: Data preprocessing (cleaning raw reads, decontamination, assembly)](https://github.com/chutter/PhyloCap/wiki/Tutorial-1)

[Tutorial 2: Data alignment (target marker matching, alignment, alignment filtering)](https://github.com/chutter/PhyloCap/wiki/Tutorial-2)

Tutorial 3: Tree estimation (gene trees, concatenation, species trees) COMING SOON

[Troubleshooting: random issues that arise](https://github.com/chutter/PhyloCap/wiki/Troubleshooting)



