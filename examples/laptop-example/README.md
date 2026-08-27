# PhyloProcessR seven-sample laptop example

This compact example demonstrates how independently callable PhyloProcessR
functions can be composed into a reproducible target-capture workflow. It
includes seven biological samples: five single-lane samples (`CRH01`-`CRH05`),
one two-lane sample (`Test_dataset_1`), and one four-lane sample
(`Test_dataset_2`). The lane structure is retained so input discovery and the
combination of technical lanes are tested without treating lanes as taxa.

The fixture is a functional example, not a biological benchmark or a
performance comparison with another pipeline. The full source read files and
the original two-sample fixture remain unchanged.

## Data license

The reduced FASTQ fixture, 40-locus target panel, and expected biological
outputs are released under the Creative Commons CC0 1.0 Universal public-domain
dedication; see `DATA_LICENSE.md`. R and shell scripts remain covered by the
repository's GPL-3.0-or-later software license.

## Compact marker and read selection

The five CRH libraries were screened against the 13,515 sequences in
`Ranoidea-V1_Markers.fa`. There were 5,124 gene-unique markers with mapped
reads in all five samples. Candidate markers were ranked by support in the
weakest sample, with one marker retained per gene. The two loci previously
recovered from both original samples were forced into a 40-marker candidate
panel so the example contains both broad five-taxon alignments and two
seven-taxon alignments.

Both mates were retained whenever either mate mapped to the selected panel.
Each CRH sample was capped at 10,000 paired reads; the original two samples
retain 2,138 pairs across their six lanes. The complete fixture therefore
contains 52,138 read pairs in 22 compressed FASTQ files (about 9.2 MB).
Selection tables and preparation details are stored in `provenance/`.

## Requirements

Install PhyloProcessR and create its Conda environment as described in the
main repository documentation. The preprocessing and assembly example uses
`fastp`, SPAdes, CD-HIT-EST, and BLAST. The alignment example additionally
uses MAFFT. Defaults are two CPU threads and a 4 GB SPAdes memory limit. No
cluster scheduler, network access, or full FrogCap marker file is required.

## Run preprocessing and assembly

```bash
export PHYLOPROCESSR_BIN="$CONDA_PREFIX/bin"
Rscript run_laptop_example.R
```

During package development, load a checkout directly:

```bash
export PHYLOPROCESSR_BIN="$CONDA_PREFIX/bin"
export PHYLOPROCESSR_SOURCE=/path/to/PhyloProcessR
Rscript run_laptop_example.R
```

The script composes `fastqStats()`, `fastpComplete()`, `assembleSpades()`,
`reduceRedundancy()`, and `removeOffTargetContigs()`. Outputs are written under
`results/`, including per-sample contigs, count summaries, target recovery,
step timings, and the software environment.

On the validation machine, the complete preprocessing and assembly run took
about 16.9 minutes with two threads and a 4 GB memory cap. The five CRH
samples recovered 39-44 target contigs, including paralog-suffixed matches;
both original samples recovered the same two expected loci. Thirty-nine of
the 40 base loci occurred in at least four samples.

## Run annotation and alignment

After preprocessing and assembly completes, run:

```bash
Rscript run_alignment_example.R
```

This composes `annotateTargets()` and `alignTargets()`, uses a four-taxon
minimum, and ranks completed alignments by taxon occupancy, mean missingness,
alignment length, and locus name. All 38 completed candidate alignments are
retained in `results/alignments/forty-candidate/`; the selected review set is
written to `results/alignments/best-20/` with statistics in
`results/best_20_alignment_summary.tsv`.

Validation produced 38 candidate alignments. The selected 20 contain five to
seven taxa and span 983-12,146 bp. The two seven-taxon loci are
`ranoidea-01644_kif2a-ex2` and `ranoidea-05399_ckap5-ex3`.

## Run trimming and QC

After creating the best-20 alignment set, run:

```bash
Rscript run_trimming_qc_example.R
```

This composes `trimAlignmentTargets()` and `superTrimmer()`. Alignments are
restricted to the reference target region, processed with TrimAl, screened for
samples at least 0.40 distant from the consensus, trimmed at poorly covered
edges and columns, filtered for sample coverage, and assessed against the
four-taxon, 100-bp, and 50%-gap thresholds.

All 20 alignments passed the validation run, and no samples were removed. The
final alignments contain five to seven taxa, span 205-11,589 bp, and contain
0-13.13% gaps. Target-region trimming and QC took approximately 28 seconds
with two threads. Results are written to
`results/alignments/best-20_trimmed/`, with per-locus and per-sample summaries
in `results/final_20_qc_summary.tsv` and
`results/final_20_sample_occupancy.tsv`.

## Interpretation

This example is designed to exercise modular workflow construction, mixed
single- and multilane inputs, missing data, target recovery, annotation,
alignment, trimming, and QC on a laptop. Runtime is informational rather than a pass/fail
criterion. Exact contig sequences can vary with dependency versions; sample
counts, recovered locus identities, and completed alignment counts are the
primary checks documented in `expected/`.
