# Assess PhyloProcessR results

Use summary functions to find failed samples, low capture efficiency, high
missing data, or possible paralogs. Run each function on a separate output
directory. Keep `overwrite = FALSE` until you intend to replace a summary.

When you use the supplied Conda environment, define its executable directory
once:

```r
tool_bin <- file.path(Sys.getenv("CONDA_PREFIX"), "bin")
```

## Summarize paired reads

`readStats()` counts lanes, read pairs, reads, nucleotides, and gigabases for
each sample. The input directory must contain one subdirectory for each sample.

```r
readStats(
  read.directory = "processed-reads/cleaned-reads",
  output.directory = "data-analysis/read-stats",
  overwrite = FALSE
)
```

The function writes `sample-read-summary.txt` and per-lane sample tables.
Compare read counts across samples. Investigate samples with unexpectedly low
counts.

## Summarize assemblies

`assemblyStats()` reports total bases, total contigs, and contig-length
statistics. The input directory must contain one FASTA file for each sample.

```r
assemblyStats(
  assembly.directory = "data-analysis/contigs/1_draft-contigs",
  output.directory = "data-analysis/assembly-stats",
  overwrite = FALSE
)
```

The function writes `sample-assembly-summary.txt`. Compare total bases and
contig counts across samples. A large difference can indicate low read depth,
contamination, or an assembly problem.

## Assess putative paralogs

`paralogStats()` uses BLAST to find target loci with more than one matching
contig in a sample.

```r
paralogStats(
  assembly.directory = "data-analysis/contigs/1_draft-contigs",
  target.file = "target-markers.fa",
  output.directory = "data-analysis/paralog-stats",
  min.match.percent = 65,
  min.match.length = 60,
  min.match.coverage = 35,
  threads = 4,
  memory = 8,
  overwrite = FALSE,
  quiet = TRUE,
  blast.path = tool_bin
)
```

The function writes `sample-paralog-summary.txt`. A second hit is a putative
paralog, not proof of gene duplication. Review sequence length, similarity,
coverage, and the locus alignment before you remove data.

## Calculate capture specificity

`sampleSpecificity()` calculates the proportion of cleaned read pairs that map
to the target markers. It uses BWA, GATK, and Samtools.

```r
sampleSpecificity(
  read.directory = "processed-reads/cleaned-reads",
  target.file = "target-markers.fa",
  output.directory = "data-analysis/sample-specificity",
  threads = 4,
  memory = 8,
  overwrite = FALSE,
  quiet = TRUE,
  bwa.path = tool_bin,
  gatk4.path = tool_bin,
  samtools.path = tool_bin
)
```

The function writes `sample-specificity_summary.txt` and a raw table for each
sample. Low specificity means that many cleaned reads do not map to the target
set. It does not identify the cause by itself.

## Calculate capture sensitivity

`sampleSensitivity()` calculates the recovered proportion of each target for
each sample. The input must contain PHYLIP alignments whose file names match
the target names.

```r
sampleSensitivity(
  alignment.directory = "data-analysis/alignments/untrimmed_all-markers",
  target.file = "target-markers.fa",
  output.directory = "data-analysis/sample-sensitivity",
  threads = 4,
  memory = 8,
  mafft.path = tool_bin,
  overwrite = FALSE,
  quiet = TRUE
)
```

The function writes raw per-marker values and
`sample-sensitivity_summary.txt`. Low sensitivity means that a sample recovered
only a small part of the target regions.

## Assess alignments

Workflow 4 writes per-locus and per-sample alignment summaries in `logs/`.
Workflow 5 can also run `alignmentAssess()` through `superTrimmer()`.

Review these values:

- number of taxa in each alignment;
- alignment length;
- missing bases and gaps;
- samples removed by the similarity filter;
- samples removed by the coverage filter;
- loci that fail the minimum-taxon or minimum-length threshold.

Do not use one threshold for all projects without review. Select thresholds
that are applicable to the taxonomic scale, target design, and downstream
method.

## Interpret specificity and sensitivity together

Specificity and sensitivity measure different properties:

- High specificity means that many reads map to the targets.
- High sensitivity means that much of each target region is recovered.

A sample can have high specificity and low sensitivity when reads concentrate
on a small number of loci. A sample can have lower specificity and good
sensitivity when sufficient reads cover most targets.

Use both measures with the read, assembly, and alignment summaries.
