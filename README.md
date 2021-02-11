# PHYLOCAP

R package For processing high-throughput sequencing data from targeted sequence capture.


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

The program dependencies can be installed using the provided bash script ("bash_install_phylocap.sh") or PBS cluster script ("qsub_install_phylocap.sh"). More detailed manual installation directions are provided in the <b>Installation</b> tutorial below if those scripts do not work. Note that for the PBS script, usernames and file paths need to be modified to match yours. 

All the functions for PhyloCap should be ready to go! 

< coming soon a function to test if they can found >


# Installation of R package

The main functions of PhyloCap are contained in an R package that has been tested on R version 3.5 and use the listed programs above along with custom scripts. To install PhyloCap from GitHub, you can use the R package devtools. Here are step-by-step instructions for installation:

1) Install devtools by typing in your R console: install.packages("devtools", dependencies = TRUE)

2) Install PHYLOCAP by typing in your R console: devtools::install_github("chutter/PHYLOCAP", update = "never", dependencies = FALSE)

3) If Devtools asks you to install the package dependencies, select "No", because the dependencies are installed above. If devtools asks you to update packages, you may choose to do so. I would recommend not to install packages from source if devtools asks you and to include the update parameter set to "never" in install_github. 

4) Devtools should finish and say the package loaded properly. Load the package with library(PHYLOCAP) in your R script. 

And installation should be done! 


# PhyloCap pipeline tutorials 

[Installation: detailed installation instructions and trouble-shooting ](https://github.com/chutter/PhyloCap/wiki/Installation)

[Tutorial 1: Data preprocessing (cleaning raw reads, decontamination, assembly)](https://github.com/chutter/PhyloCap/wiki/Tutorial-1)

[Tutorial 2: Data alignment (target marker matching, alignment, alignment filtering)](https://github.com/chutter/PhyloCap/wiki/Tutorial-2)

Tutorial 3: Tree estimation (gene trees, concatenation, species trees) COMING SOON

[Troubleshooting: random issues that arise](https://github.com/chutter/PhyloCap/wiki/Troubleshooting)



