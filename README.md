# PhyloProcessR

R package for processing high-throughput sequencing data from raw reads to alignments for many samples from targeted sequence capture for use in phylogenomic/phylogenetic analyses. 

The R package and pipeline does the following:

1) Organize raw read data
2) Remove adaptor contamination and merge paired-end reads
3) Decontaminate reads from other organisms
4) Assemble cleaned reads into contigs 
5) Use a sample-based iterative mapping approach to call SNPs and export for popular programs
6) Match contigs to design targets for sequence capture
7) Align and trim contigs from samples
8) Concatenate all targets or only targets from the same gene


# PhyloProcessR prerequisites 

PhyloProcessR uses several R packages and other outside programs for certain functions, which will be installed all at once below using an anaconda environment file.

1. R packages (R version 4.2.2 tested)
- From CRAN: devtools, ape, stringr, data.table, seqinr, foreach, doparallel, rdrop2, biomartr
- From BioConductor: rsamtools, genomicranges, biostrings

2. Stand alone programs
- fastp: adaptor trimming and paired-end read merging
- ORNA: read normalization
- bwa: read mapping
- hisat2: alternative mapper
- spades: assembly
- BLAST: matching assembled contigs to targets, other utilities
- mafft: creating alignments
- trimal: trimming alignments
- IqTree: gene tree and concatenation trees
- GATK4: variant calling functions
- SamTools: variant calling and read mapping tools


# Quick installation instructions

First, clone this repository to your computer to obtain the setup files. Or alternatively go to the green "Code" button in top right of this repository and select "download ZIP".

```bash
git clone https://github.com/chutter/PhyloProcessR.git
```

Second, change your working directory in the terminal to the downloaded repository. The key file here is the "environment.yml" anaconda environment file, which must be present in the working directory being used. 

```bash
cd PhyloProcessR/setup-files/
```

The R packages and outside programs can be installed manually or more easily through the anaconda environment file provided (version numbers are provided in environment file for reporting and exact replication). To install with the environment file, the easiest and quickest way is to first install the Anaconda package manager. Anaconda can be downloaded and installed for different operating systems from https://anaconda.org. Miniconda is recommended as it has a smaller footprint (smaller size and fewer files). Once installed, you can create a new environment for PhyloProcessR by: 

```bash
conda env create -f environment.yml -n PhyloProcessR
```

**** WARNING: It is possible that the environment file may fail, however, it has been tested on Linux and MacOS on April 3 2023 and installed fine. For MacOS, you must use the X84 (not M1) version of anaconda as most packages are not available for M1 but can be emulated through X84. Occasionally things break and there are manual installation methods in the Wiki (the first tutorial). 

And finally, the cloned GitHub directory may be deleted after installing the prerequisites through the conda env file that manually installs the anaconda environment. There are some useful example files (also in the tutorial here), which could be saved.   

To use the environment, it must first be activated in your current terminal session or placed in your cluster job script. 

```bash
conda activate PhyloProcessR
```

# Installation of R package

The main functions of PhyloProcessR are contained in an R package that has been tested on R version 4.0.2 and use the listed programs above along with custom scripts. To install PhyloProcessR from GitHub, you can use the R package devtools included in the environment above. When running in a cluster environment, the code for installation here should be included at the top of your R script with your selected PhyloProcessR functions. Here are step-by-step instructions for installation:

1) Install PhyloProcessR by typing in your R console: 

```R
devtools::install_github("chutter/PhyloProcessR", upgrade = "never", dependencies = FALSE)
```

The update = "never" flag ensures that packages already installed via the anaconda environment are not changed, which will often break things. Additionally, dependencies = FALSE is set for the same reason. 


2) Devtools should finish and say the package loaded properly with no errors. Load the package in your R script with:

```R
library(PhyloProcessR)
```

And installation should be done! All the functions for PhyloProcessR should be ready to go! It is recommended to keep the install line above in your R script as the package is frequently updated for bugs and other features. 


3) You can run the following function to see if PhyloProcessR can find the dependencies: 

< coming soon a function to test if they can found >


# PhyloProcessR workflows

PhyloProcessR is organised into a series of workflows, each covering a distinct stage of the pipeline. Configuration files and R scripts for each workflow are provided in the `workflows/` directory.

| Workflow | Script | Description |
|---|---|---|
| **Workflow 1** | `workflow-1_preprocess.R` | Organise raw reads, remove adaptors, decontaminate, normalise, and merge paired-end reads |
| **Workflow 2** | `workflow-2_assembly.R` | De novo assembly with SPAdes; match contigs to target markers |
| **Workflow 3** | `workflow-3_variantCalling.R` | SNP calling and IUPAC/haplotype consensus generation |
| **Workflow 4** | `workflow-4_alignment.R` | Align target markers across samples |
| **Workflow 5** | `workflow-5_trimming.R` | Trim alignments, concatenate genes, build unlinked dataset; optionally include novel markers from Workflow X4 |
| **Workflow X3** | `workflow-X3_legacy-integration.R` | Integrate Sanger/GenBank legacy alignments into the capture dataset; supports NEXUS conversion and mitochondrial loci |
| **Workflow X4** | `workflow-X4_novel-loci.R` | Discover novel shared genomic regions from unmapped reads, assemble and align them as new loci |

Each workflow has a matching configuration file (e.g. `workflow-1_configuration-file.R`) where all parameters are set. See the tutorials below for detailed guidance.


# PhyloProcessR pipeline tutorials

[Installation: detailed installation instructions and trouble-shooting](https://github.com/chutter/PhyloProcessR/wiki/Installation:-detailed-installation-instructions-and-trouble-shooting)

[Tutorial 1: PhyloProcessR configuration](https://github.com/chutter/PhyloProcessR/wiki/Tutorial-1:-PhyloProcessR-configuration)
— Setting up working directories, renaming files, and configuring the decontamination database.

[Tutorial 2: PhyloProcessR pipeline workflows](https://github.com/chutter/PhyloProcessR/wiki/Tutorial-2:-PhyloProcessR-pipeline-workflows)
— Step-by-step guide to running Workflows 1–5, X3, and X4, including expected outputs and directory structures.

[Tutorial 3: Assess sequence capture results](https://github.com/chutter/PhyloProcessR/wiki/Tutorial-3:-Assess-results)
— Summarise capture success across samples and loci.

[Tutorial 4: Legacy data integration (Workflow X3)](https://github.com/chutter/PhyloProcessR/wiki/Tutorial-4:-Legacy-Integration)
— Full guide to integrating Sanger or GenBank alignments into a sequence-capture dataset, including NEXUS conversion, name-matching strategies, mitochondrial loci, and gene concatenation.

