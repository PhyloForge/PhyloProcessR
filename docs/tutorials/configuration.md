# Configure a PhyloProcessR project

PhyloProcessR supplies complete reference workflows. Each operation is also an
independently callable R function. You can omit stages or combine functions in
a project-specific script.

## 1. Copy the workflow files

Copy the required script and its configuration file from `workflows/` into
your project directory. Keep each pair together.

For example, Workflow 1 uses:

```text
workflow-1_preprocess.R
workflow-1_configuration-file.R
```

The script reads the configuration file from the current directory.

## 2. Set the main project values

Edit the first part of each configuration file. Use full paths for input files
and external program directories.

Typical Workflow 1 values are:

```r
install.latest.github = FALSE
working.directory = "/full/path/to/project"
dataset.name = "project-name"
read.directory = "/full/path/to/raw-reads"
sample.file = "file_rename.csv"
target.fasta = "/full/path/to/target-markers.fa"
threads = 8
memory = 32
overwrite = FALSE
quiet = TRUE
conda.env = "/full/path/to/envs/PhyloProcessR/bin"
```

Keep `overwrite = FALSE` for the first run. The workflow can then keep completed
results. Set `overwrite = TRUE` only when you intend to replace output.

Keep `install.latest.github = FALSE` for normal work. Set it to `TRUE` only to
test the current development version.

## 3. Create the sample-name file

Use a comma-separated file with the columns `File` and `Sample`.

```csv
File,Sample
CRH111,Spinomantis_elegans_CRH111
CRH1644,Aglyptodactylus_securifer_CRH1644
CRH0481,Boophis_burgeri_CRH0481
```

The `File` value is a unique part of the input file name. Do not include the
read-pair or lane suffix.

For these files:

```text
CRH111_AX1212_L001_R1.fastq.gz
CRH111_AX1212_L001_R2.fastq.gz
```

you can use `CRH111_AX1212` as the `File` value.

The `Sample` value is the name used in contigs, alignments, and summary files.
Use a unique value for each biological sample. Use underscores instead of
spaces. Do not use one sample name as the prefix of another sample name.

Use one row for each lane or library. Give each input a different `File` value.
Use the same `Sample` value for inputs from one biological sample:

```csv
File,Sample
CRH111_L001,Spinomantis_elegans_CRH111
CRH111_L002,Spinomantis_elegans_CRH111
CRH111_LIB03,Spinomantis_elegans_CRH111
```

PhyloProcessR processes the lanes separately and combines them for assembly.

## 4. Configure contamination removal

Set `decontamination = TRUE` to remove read pairs that map to contaminant
references. Use the example file at
`setup-files/decontamination_database.csv`.

Set these values when PhyloProcessR must download the references:

```r
decontamination = TRUE
contaminant.genome.list = "decontamination_database.csv"
download.contaminant.genomes = TRUE
decontamination.path = NULL
```

Set `decontamination.path` to an existing reference directory when you cannot
download files during the analysis.

## 5. Set external program paths

The `*.path` parameters must identify the directory that contains each
executable. A Conda environment usually uses one directory for all programs:

```r
conda.env = "/full/path/to/envs/PhyloProcessR/bin"
fastp.path = conda.env
samtools.path = conda.env
bwa.path = conda.env
spades.path = conda.env
blast.path = conda.env
```

Each workflow checks only the programs that it uses.

## 6. Use functions in a custom workflow

You do not have to run the complete scripts. Load the package and call only the
required functions:

```r
library(PhyloProcessR)
tool_bin <- file.path(Sys.getenv("CONDA_PREFIX"), "bin")

fastqStats(
  read.directory = "raw-reads",
  output.name = "raw-read-summary",
  threads = 2,
  mem = 4
)

fastpComplete(
  input.reads = "raw-reads",
  output.directory = "processed-reads/cleaned-reads",
  fastp.path = tool_bin,
  threads = 2,
  mem = 4
)
```

Keep the custom R script with the analysis. Record all parameter values in the
script. This makes the project-specific workflow reproducible.

## Next step

Continue with [Run the workflows](workflows.md).
