# PhyloProcessR

[![R CMD check](https://github.com/chutter/PhyloProcessR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/chutter/PhyloProcessR/actions/workflows/R-CMD-check.yaml)
[![License: GPL v3+](https://img.shields.io/badge/License-GPL_v3%2B-blue.svg)](LICENSE)

PhyloProcessR is an R-based workflow framework for converting raw targeted
sequence-capture reads into curated, analysis-ready phylogenomic datasets. It
combines domain-specific R functions with established bioinformatics programs
for stages that require external executables.

The package and supplied workflows support:

1. Organizing raw read data.
2. Removing adapters, quality filtering, normalizing, and merging paired-end reads.
3. Removing contaminant reads.
4. Assembling cleaned reads and recovering targeted loci.
5. Mapping reads, calling variants, and generating IUPAC or haplotype consensuses.
6. Aligning, filtering, and trimming recovered loci.
7. Assessing capture, missing data, depth, paralogy, and alignment quality.
8. Integrating legacy sequences and discovering shared novel loci.
9. Concatenating loci by target or gene for downstream analyses.

The R package can be installed and tested without every external program.
Individual workflow stages require only the command-line tools they invoke; the
provided container and Conda environment install the complete environment.

## Prerequisites

PhyloProcessR requires R 4.0 or later. Its R package dependencies are declared
in `DESCRIPTION` and installed by standard R package installers. Depending on
the workflow, external tools can include fastp, BWA, HISAT2, SPAdes, BLAST,
MAFFT, trimAl, IQ-TREE, GATK, and Samtools. See the workflow configuration files
for the programs used by each stage.

## Complete workflow environment

Clone the repository to obtain the container recipes, Conda environment, and
workflow configuration files:

```bash
git clone https://github.com/chutter/PhyloProcessR.git
```

Then change to the setup directory:

```bash
cd PhyloProcessR/setup-files/
```

The recommended way to run the full pipeline is with Docker for local Linux or
macOS use, or Apptainer/Singularity on an HPC system. This keeps the R and
bioinformatics dependencies in a reproducible environment.

### Option 1: Run a pre-built container (recommended)

Pull the published image from Docker Hub.

For Docker:
```bash
docker pull chutter/phyloprocessr:latest
```

Mount the current working directory and run a workflow script:
```bash
docker run -v $(pwd):/app -w /app chutter/phyloprocessr:latest Rscript workflow-5_trimming.R
```

For Apptainer or Singularity, convert the same image to a local `.sif` file:
```bash
apptainer pull phyloprocessr.sif docker://chutter/phyloprocessr:latest
```

Then run a workflow:
```bash
apptainer exec phyloprocessr.sif Rscript workflow-5_trimming.R
```

Inside either container, tools are on `PATH`; workflow configurations can use
`conda.env = ""` or `/opt/conda/bin/` where a tool-directory argument is
required.

### Option 2: Build a container from the recipes

Use this option when you need to modify the supplied environment.

Build Docker from the `setup-files` directory. On Apple Silicon, the
`linux/amd64` platform maintains compatibility with tools distributed only for
x86_64 Linux.
```bash
docker build --platform linux/amd64 -t chutter/phyloprocessr:1.0.0 -t chutter/phyloprocessr:latest .
```

Build an Apptainer image from the supplied definition:
```bash
apptainer build phyloprocessr.sif Apptainer.def
```

### Option 3: Create the Conda environment

After installing a Conda-compatible package manager, create the environment
from `setup-files/environment.yml`:

```bash
conda env create -f environment.yml -n PhyloProcessR
```

Activate it before running a workflow:

```bash
conda activate PhyloProcessR
```

## Install the R package

Install the development version from GitHub:

```R
install.packages("remotes")
remotes::install_github("chutter/PhyloProcessR")
```

When working inside the supplied Conda environment or container, its pinned R
dependencies are already present. Install a local checkout without upgrading
them:

```R
remotes::install_local(".", upgrade = "never", dependencies = FALSE)
```

Load the package with:

```R
library(PhyloProcessR)
```

For reproducible analyses, record the installed package version or Git commit
rather than reinstalling the moving development branch in every script.

To check external programs in the active Conda environment:

```R
setupCheck(anaconda.environment = Sys.getenv("CONDA_PREFIX"))
```


## Workflows

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


## Tutorials

[Installation: detailed installation instructions and trouble-shooting](https://github.com/chutter/PhyloProcessR/wiki/Installation:-detailed-installation-instructions-and-trouble-shooting)

[Tutorial 1: PhyloProcessR configuration](https://github.com/chutter/PhyloProcessR/wiki/Tutorial-1:-PhyloProcessR-configuration)
— Setting up working directories, renaming files, and configuring the decontamination database.

[Tutorial 2: PhyloProcessR pipeline workflows](https://github.com/chutter/PhyloProcessR/wiki/Tutorial-2:-PhyloProcessR-pipeline-workflows)
— Step-by-step guide to running Workflows 1–5, X3, and X4, including expected outputs and directory structures.

[Tutorial 3: Assess sequence capture results](https://github.com/chutter/PhyloProcessR/wiki/Tutorial-3:-Assess-results)
— Summarise capture success across samples and loci.

[Tutorial 4: Legacy data integration (Workflow X3)](https://github.com/chutter/PhyloProcessR/wiki/Tutorial-4:-Legacy-Integration)
— Full guide to integrating Sanger or GenBank alignments into a sequence-capture dataset, including NEXUS conversion, name-matching strategies, mitochondrial loci, and gene concatenation.

## Citation, contributing, and license

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). Contributions
are welcome; see [`CONTRIBUTING.md`](CONTRIBUTING.md) and
[`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md). PhyloProcessR is distributed under
the [GNU General Public License, version 3 or later](LICENSE).
