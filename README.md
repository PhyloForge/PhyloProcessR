# PHYLOCAP

R package For processing high-throughput sequencing data from targeted sequence capture

Main functionality coming soon! Some implemented functions and their purpose are as follows: 

Coming Soon!

# Installation of R package

PHYLOCAP is an R package that has been tested on R version 3.5 or greater. To install PHYLOCAP, you can use the R package devtools. Here are step-by-step instructions for installation:

1) Install devtools by typing in your R console: install.packages("devtools")

2) Install PHYLOCAP by typing in your R console: devtools::install_github("chutter/PHYLOCAP", update = "never")

3) Devtools will ask you to install the package dependecies, select "Yes". If devtools asks you to update packages, you may choose to do so. I would recommend not to install packages from source if devtools asks you and to include the update parameter set to "never" in install_github. Ape is problemic from source and I could not get it to install on my machine. If devtools is still giving you trouble, you can install the dependencies with "install.packages(c("ape", "stringr"))". Then rerun Step 2 and skip package updates. 

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

The R packages and outside programs can be installed manually. The easiest and quickest way is to use the Anaconda package manager. Anaconda can be downloaded and installed for different operating systems from https://anaconda.org. Miniconda is recommended. Once a version of Anaconda is installed or loaded, a new clean environment should be created for PhyloCap. Anaconda can be configured through the following steps: 

1. Create new environment:

```
conda create phylocap

conda activate phylocap
```


2. Next, you should configure anaconda by adding and ordering channels as follows:

```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority flexible
```


3. Finally, the required programs can be installed in two commands. The programs should be installed together so that the correct package dependency versions can be sorted out. 

```
conda install -c conda-forge r-base=3.5 r-devtools r-ape r-stringr r-data.table r-seqinr r-foreach r-doparallel r-rdrop2

conda install -c bioconda bioconductor-rsamtools bioconductor-genomicranges bioconductor-biostrings fastp spades mafft bwa samtools gatk4 trimal iqtree
```

All the functions for PhyloCap should be ready to go! 

< coming soon a function to test if they can found >


# PhyloCap pipeline tutorials 

1. Preprocessing raw reads and assembly 
2. Assembly contig matching 
3. Alignment and trimming 





1) first install and load the R package. Its a good idea to install new every time as this package is being updated frequently. The warning that the package has no updates can be ignored. 

```r
devtools::install_github("chutter/PHYLOCAP")
library(PHYLOCAP)

```
\
