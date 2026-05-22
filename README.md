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


# PhyloProcessR pipeline tutorials 

[Installation: detailed installation instructions and trouble-shooting ](https://github.com/chutter/PhyloProcessR/wiki/Installation:-detailed-installation-instructions-and-trouble-shooting)

[Tutorial 1: PhyloProcessR configuration](https://github.com/chutter/PhyloProcessR/wiki/Tutorial-1:-PhyloProcessR-configuration)

[Tutorial 2: PhyloProcessR pipeline workflows](https://github.com/chutter/PhyloProcessR/wiki/Tutorial-2:-PhyloProcessR-pipeline-workflows)

[Tutorial 3: Assess sequence capture results](https://github.com/chutter/PhyloProcessR/wiki/Tutorial-3:-Assess-results)

[Tutorial 4: Combine legacy genbank data with sequence capture](https://github.com/chutter/PhyloProcessR/wiki/Tutorial-4:-Legacy-Integration)


# Legacy data integration (`integrateLegacy`)

Workflow X3 integrates Sanger or GenBank legacy sequence alignments into an existing sequence-capture dataset. For each legacy alignment, a representative sequence is BLASTed against the target marker file to identify the matching capture locus. The legacy sequences are then added to that locus alignment via MAFFT, and optionally merged with any capture sequences from the same sample.

### Mitochondrial loci

Mitochondrial loci (e.g. 12S, 16S, ND1) are not present in a typical nuclear probe set and will return no BLAST hit against the nuclear target file. Set `include.mitochondrial = TRUE` and point `mito.alignment.directory` at a set of existing mitochondrial capture alignments (e.g. produced by [MitoTrawlR](https://github.com/chutter/MitoTrawlR)) to enable a second BLAST step that matches these loci automatically.

### Controlling how sample names are matched: `name.match`

Because legacy (Sanger/GenBank) datasets often use different voucher IDs or formatting conventions from sequence-capture datasets, three matching strategies are available:

| Mode | Merge identical names? | Merge same species? | Handles `-` vs `_` differences? | Pre-selects one legacy seq per species? | Unmatched legacy sequences |
|---|:---:|:---:|:---:|:---:|---|
| `"exact"` | ✅ | ❌ | ❌ | ❌ | Added as separate rows |
| `"fuzzy"` | ✅ | ❌ | ✅ | ❌ | Added as separate rows |
| `"species"` | ✅ | ✅ | ❌ | ✅ | Excluded (one per species only) |

**`"exact"`** — sequence names must be identical for merging. All legacy sequences are added to the alignment; those with a name that exactly matches a capture sequence are merged column-by-column (preferring non-gap characters at each site).

**`"fuzzy"`** — before comparing names, hyphens, dots, and spaces are converted to underscores and names are lowercased. This handles common formatting differences such as `MZUTI-2436` vs `MZUTI_2436`. Sequences whose normalised names match are merged and the merged sequence retains the original capture alignment name. All unmatched legacy sequences are added as separate rows with their original names preserved.

**`"species"`** — strips the trailing specimen/voucher ID (the last underscore-delimited field) before matching, so that e.g. `Centrolene_bacatum_MZUTI-2436` (capture) and `Centrolene_bacatum_KU12345` (legacy) are treated as the same taxon and merged into a single sequence named `Centrolene_bacatum`. When the legacy alignment contains multiple sequences for the same species, one representative is pre-selected: the sequence whose full name matches a capture specimen is preferred; otherwise the sequence with the most informative (non-gap, non-missing) bases is used.

> **Which mode to use?** If your capture and legacy datasets use the same voucher IDs, use `"exact"`. If they share species but with different voucher formatting, use `"fuzzy"`. If you want to collapse everything to species-level regardless of voucher ID, use `"species"`.

