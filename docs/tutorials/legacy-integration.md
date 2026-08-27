# Integrate legacy sequence data

Workflow X3 adds Sanger, GenBank, or other legacy sequences to an existing
sequence-capture dataset. Use the workflow to increase taxon or locus coverage
without new sequencing.

The workflow identifies a matching capture locus with BLAST. It then uses
MAFFT to add the legacy sequences to that alignment.

## Required input

Prepare these files:

1. A directory of untrimmed sequence-capture alignments.
2. A directory of legacy alignments, with one alignment for each locus.
3. The target-marker FASTA file used for the capture project.
4. The Workflow X3 script and configuration file.

Use untrimmed capture alignments when possible. Trimmed alignments can exclude
flanking sequence that overlaps a legacy marker.

Legacy alignments can use PHYLIP or FASTA format. All files in one run must use
the format set by `legacy.format`.

Optional input includes:

- a partitioned NEXUS matrix;
- mitochondrial capture alignments;
- a `Marker` and `Gene` metadata table for exon concatenation.

## Set the principal paths

Edit `workflow-X3_configuration-file.R`:

```r
working.directory = "/full/path/to/project"
target.file = "/full/path/to/target-markers.fa"
alignment.directory = "data-analysis/alignments/untrimmed_all-markers"
alignment.format = "phylip"
legacy.directory = "/full/path/to/legacy-alignments"
legacy.format = "phylip"
threads = 8
memory = 40
overwrite = FALSE
quiet = TRUE
```

Run the workflow:

```bash
Rscript workflow-X3_legacy-integration.R
```

## Convert a partitioned NEXUS matrix

Enable NEXUS conversion when the legacy data is one concatenated matrix with
a `BEGIN SETS` block and `charset` definitions:

```r
convert.nexus = TRUE
nexus.file = "/full/path/to/legacy-matrix.nex"
nexus.output.directory = "data-analysis/alignments/legacy-alignments"
max.missing.percent = 100
```

`convertNexusPartitions()` writes one PHYLIP file for each charset. It joins
multiple ranges for the same charset in their listed order. It converts `?` to
`N` for downstream compatibility.

Reduce `max.missing.percent` to remove a sample from a locus when its missing
percentage is too high.

## Select a name-matching method

Set `name.match` to `"exact"`, `"fuzzy"`, or `"species"`.

### Exact matching

Use exact matching when both datasets use identical sample names:

```r
name.match = "exact"
```

For example, `MZUTI-2436` matches only `MZUTI-2436`.

### Fuzzy matching

Use fuzzy matching when both datasets use the same voucher but different
separators:

```r
name.match = "fuzzy"
```

The comparison removes hyphens, underscores, periods, and spaces. It also uses
lowercase characters. Thus, `MZUTI-2436`, `MZUTI_2436`, and `MZUTI2436` match.

### Species matching

Use species matching when the datasets contain the same species but different
voucher identifiers:

```r
name.match = "species"
```

The comparison removes the last underscore-delimited name field. For example,
`Centrolene_bacatum_MZUTI-2436` and
`Centrolene_bacatum_KU12345` become `Centrolene_bacatum`.

Review names before you use this method. The method assumes that the last field
is a voucher identifier.

## Merge sequences from the same sample

Keep the default value to merge matched capture and legacy sequences:

```r
combine.same.sample = TRUE
```

The merge examines each alignment column. It uses an available base when the
other sequence has a gap or `N`. When both sequences have a base, it keeps the
capture base.

Set `combine.same.sample = FALSE` to keep both sequences as separate rows.

## Include unmatched legacy loci

The default behavior skips a legacy locus when BLAST finds no capture target:

```r
include.uncaptured.legacy = FALSE
```

Set the value to `TRUE` to save an unmatched legacy locus as a separate
alignment.

Set this value to write a complete dataset that contains unchanged capture
alignments and integrated alignments:

```r
include.all.together = TRUE
```

## Include mitochondrial loci

Nuclear target files usually do not contain mitochondrial markers. Supply
mitochondrial capture alignments when the legacy dataset contains loci such as
12S, 16S, or ND1:

```r
include.mitochondrial = TRUE
mito.alignment.directory = "/full/path/to/mitochondrial-alignments"
mito.alignment.format = "phylip"
```

Workflow X3 first searches the nuclear targets. If that search has no hit, it
searches consensus sequences from the mitochondrial alignments.

[MitoTrawlR](https://github.com/chutter/MitoTrawlR) can produce mitochondrial
capture alignments from sequence-capture reads.

## Concatenate exons and gather unlinked loci

Enable gene concatenation when the target set contains multiple exons from the
same gene:

```r
feature.gene.names = "data-analysis/gene_metadata.txt"
concatenate.genes = TRUE
minimum.exons = 2
concatenate.legacy.genes = TRUE
gather.unlinked = TRUE
```

The metadata file must contain the columns `Marker` and `Gene`.

Set `concatenate.legacy.genes = TRUE` to include integrated legacy loci in the
gene alignments. Set it to `FALSE` to concatenate only capture exons.

`gatherUnlinked()` uses one concatenated alignment for each multi-exon gene. It
also keeps single-exon genes and other unlinked loci.

## Trim the integrated dataset

Enable trimming after you inspect the untrimmed integrated output:

```r
trim.alignments = TRUE
min.taxa.alignment = 4
min.alignment.length = 100
run.TrimAl = TRUE
trim.column = TRUE
min.column.gap.percent = 50
trim.external = TRUE
min.external.percent = 50
trim.coverage = TRUE
min.coverage.percent = 35
min.coverage.bp = 60
```

These values are starting points. Confirm that they are applicable to the
project before the final analysis.

## Output directories

Workflow X3 writes output under `data-analysis/legacy-integration/`.

| Directory | Content |
|---|---|
| `untrimmed_legacy-only` | Loci that received legacy data |
| `untrimmed_legacy-all` | Complete capture and legacy dataset |
| `untrimmed_legacy-genes` | Optional concatenated genes |
| `untrimmed_legacy-unlinked` | Optional unlinked dataset |
| `trimmed_legacy-unlinked` | Optional final trimmed dataset |

The workflow also writes
`untrimmed_legacy-integration_summary.txt`. The table reports locus counts and
informative base pairs for each taxon.

Review these columns:

| Column | Meaning |
|---|---|
| `Taxon` | Sample or species name |
| `Total_Loci` | Number of alignments that contain the taxon |
| `Capture_Loci` | Capture loci without added legacy data |
| `Legacy_Nuclear_Loci` | Integrated nuclear loci |
| `Legacy_Mito_Loci` | Integrated mitochondrial loci |
| `Legacy_Total` | Total integrated nuclear and mitochondrial loci |
| `Total_BP` | Informative bases, without gaps or `N` |
| `Pct_Loci` | Percentage of available loci that contain the taxon |

Use this summary to find taxa with low legacy coverage or no matching capture
sample.

## Resume a run

Keep `overwrite = FALSE` to preserve existing output. Workflow X3 skips output
that is already complete where supported.

Set `overwrite = TRUE` only after you confirm the input and output paths. The
workflow can remove and recreate output directories.
