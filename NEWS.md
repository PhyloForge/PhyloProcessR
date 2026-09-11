# PhyloProcessR (development version)

## Workflow X3 legacy integration

- Workflow X3 now saves all results below the configurable
  `output.directory`, which defaults to `data-analysis/legacy-integration`.
- Workflow 4 can copy the complete integrated X3 alignment set into
  `untrimmed_all-markers` when `include.legacy = TRUE`. These alignments are
  then available to workflow 5 trimming and dataset construction.
- Workflow 4 accepts an optional CSV, TSV, TXT, XLS, or XLSX legacy-name map.
  It renames legacy samples and merges mapped legacy and sequence-capture rows.
- Workflow 4 can add per-genome target sequences produced by
  `extractGenomeTarget()` before it aligns each capture locus.

## Workflow X5 paralog analysis

- Workflow X5 collects saved workflow 4 candidates, expands untrimmed marker
  alignments, applies copy-aware trimming, infers checked IQ-TREE gene trees,
  and exports retained markers or both groups from one supported two-copy split.
- The workflow writes complete target, copy, split, membership, and output
  reports without changing workflow 4 or workflow 5 outputs. Split markers need
  explicit biological metadata before they can enter an unlinked dataset.
- X5 locus files and tree directories use locus names. Tree and diagnostic
  alignment labels keep sample names and add `_1`, `_2`, and later suffixes
  only when a sample has multiple copies at that locus.
- X5 does not treat a deep branch as a paralog concern when every sample has one
  sequence at the locus.
- When X5 cannot support a two-copy split, it keeps the best available copy for
  each sample. Candidate rank selects the best copy. Aligned coverage and stable
  sequence metadata resolve ties. X5 re-trims the single-copy alignment and
  retains the locus only when the standard final thresholds pass.

## Workflow 1 preprocessing

### Breaking changes

- `fastpComplete()` is renamed `fastpClean()`. The new function builds the fastp
  command from logical arguments, so one pass over the reads runs every step
  that is TRUE. It takes `remove.adaptors`, `remove.duplicate.reads`,
  `error.correction`, `quality.trim.reads`, `quality.filter`,
  `low.complexity.filter`, `trim.poly.x` and `min.read.length`. The summary CSV
  is now `logs/fastpClean_summary.csv`.
- The workflow 1 configuration file no longer has `fastp.complete`. Use
  `clean.reads` to run or skip the fastp step. The separate steps are no longer
  a different code path. They are settings on the one fastp command.

### Behaviour changes

- `fastpClean()` uses `--compression 6`, the value that every other fastp step
  uses. `fastpComplete()` used `--compression 8`.
- Workflow 1 no longer calls `removeAdaptors()`, `removeDuplicateReads()`,
  `readErrorCorrection()` or `qualityTrimReads()`. The functions are still
  exported and still work on their own.

## Workflow 2 assembly

### Deprecations

- `expandMissingAssembly()` is deprecated. It warns and still runs. Use
  `assembleBinnedTargets()` instead. The new function recovers the same missing
  targets, extends the targets a sample already has, assembles one target at a
  time so a low-coverage locus keeps its own coverage distribution, and can run
  more than one round.

### Behaviour changes

- Every function that matches contigs to targets now uses the same defaults:
  `min.match.percent = 60`, `min.match.length = 50`, `min.match.coverage = 30`.
  This covers `removeOffTargetContigs()`, `curateTargetContigs()`,
  `assembleBinnedTargets()` and `annotateTargets()`. The workflow 2 and workflow
  4 configuration files were changed to match.
- `curateTargetContigs()` and `annotateTargets()` no longer count N padding as
  sequence in the length test. Two fragments of one target are joined with Ns,
  and the padding is not recovered sequence.
- `assembleBinnedTargets()` gained `rescue.failed.divergent`, now default
  `TRUE`. After round 1 it uses LAST to recruit reads for the targets that no
  bin produced, because bwa needs about 90 percent identity and a divergent
  target recruits nothing. It costs one more pass over the reads, about 8
  minutes for a 5 million read sample, and recovered 916 and 656 targets on two
  test runs, about 9 percent of the output.
- `assembleBinnedTargets()` changed `rescue.missing` to default `FALSE`. Setting
  it takes the targets with no sequence away from `rescue.failed.divergent` and
  sends them down a stricter path, where a seed must bait a bin, gate at
  `min.pairs` and clear `min.contig.length`. On a test sample `TRUE` gave 9,884
  targets and `FALSE` gave 11,111. The extra targets are short, a median of
  135 bp, and the alignment steps can drop them.
- `assembleBinnedTargets()` gained `parallel.samples`, default `1`. It divides
  `threads` and `memory` between the samples that run at the same time, the way
  `assembleSpades` does, and `binned.parallel.samples` sets it from the workflow
  2 configuration file. One sample already uses every thread it is given, since
  the bin assembly is thousands of single-threaded megahit jobs, so raise this
  to fill a node across a batch rather than to make one sample faster. About 87
  percent of a run scales with cores.
- `assembleBinnedTargets()` writes the summary row of each sample to
  `log.directory/sample_summaries/<sample>.csv` and joins them into
  `assembleBinnedTargets_summary.csv`. The joined file is written through a
  temporary file and renamed, so concurrent samples cannot lose a row or leave a
  half-written table.
- `assembleBinnedTargets()` adds `medianPreviousLength` and
  `medianBinnedLength` to the summary, because `percentExtended` is measured
  against whatever sat in `assembly.directory` and is comparable only between
  runs that began from the same contigs. It also adds `targetsFromContig`,
  `targetsFromDraft`, `targetsFromReference` and `targetsFromRescue`, which are
  the targets each bait source returned rather than the baits it was given, and
  `targetsUnderMinLength`, the contigs that reached the assembly below
  `min.contig.length`.
- `assembleBinnedTargets()` documents `mapping.reads`. Unmerged reads recover
  more targets than merged reads, 1 to 2 percent across three samples, because a
  capture insert straddles the target edge and merging makes one half-off-target
  query out of two mates. The default was already `"decontaminated-reads"`.
- `assembleBinnedTargets()` assembles each bin with megahit instead of SPAdes.
  SPAdes returns nothing for a bin below about 30 read pairs, which is a third
  to a half of all bins. The function takes `megahit.path` for this.
- The rescue step of `assembleBinnedTargets()` assembles one target at a time
  with cap3, in place of one pooled SPAdes run. It takes `cap3.path`, and
  `spades.path` is gone. A full run rescues tens of thousands of targets, which
  no single assembly can hold.
- The rescue seeds are filtered on `min.match.length` rather than
  `min.contig.length`, because cap3 seeds are short. `min.contig.length` now
  applies only to the binned contigs.

## Workflow 1 read preprocessing

### Behaviour changes

- `removeContamination()` now uses `map.match` to decide which read pairs are
  removed. A read pair is a contaminant when either mate aligns at or above the
  identity threshold. The pair is removed and counted. A pair that aligns below
  the threshold is kept. Before this change, `map.match` only set the
  contamination report, and every mapped pair was removed at any identity. Read
  sets made with an earlier version are more strongly filtered than the
  configured threshold states.
- `assessCaptureEfficiency()` counts primary alignments only, so
  `pctReadsOnTarget` can no longer go above 100.
- `removeContamination()` rebuilds the BWA index on its own when the contaminant
  reference files change. An old index can no longer be used with a new
  contaminant list.
- Workflow 1 now uses the `decontamination.path` and
  `download.contaminant.genomes` settings. Both were ignored before.
- The contaminant database is controlled by the new
  `overwrite.contaminant.database` setting, so a new read run no longer
  downloads every contaminant genome again.

### Bug fixes

- `assembleBinnedTargets()` counts a merged read as a whole insert. The bin gate
  counts BAM records and divides by two, and a merged READ3 read is one record
  that spans the whole insert, so it scored as half a pair. A library with 27
  percent merged reads lost 13.5 percent of its inserts at the gate. The read cap
  is corrected by the same change.
- `quiet = TRUE` now silences every stage of a piped command. R appends its
  redirection to the end of the string, where a shell binds it to the last stage
  only, so `samtools` was quiet and `bwa` was not.
- `mergePairedEndReads()` no longer fails with a missing argument error when the
  first search for a lane returns no files.
- `removeDuplicateReads()` passes `--dup_calc_accuracy` to fastp correctly. The
  flag was malformed, so the setting had no effect.
- Workflow 1 passes the Dropbox token to `dropboxDownload()` and no longer calls
  the `rdrop2` package, which is not a dependency.
- `dropboxDownload()` treats a lane as complete only when both read files are
  present, so a download that stopped between the two mates is finished.
- `assessCaptureEfficiency()` no longer fails when every sample is skipped.
- Sample matching uses fixed strings and a name separator, so `Sample1` no
  longer matches `Sample10`, and a sample name that holds a regular expression
  character is safe.
- All external commands are quoted, so a path that holds a space works.

### Improvements

- `assembleBinnedTargets()` skips a sample that has no contigs of its own in
  either `assembly.directory` or `draft.assembly.directory`, with a warning that
  names the fix. Such a sample puts every target in the rescue pool. One outgroup
  of 34,058 targets and 12.6 GB of reads ran for 14 hours and was about 26
  percent through the cap3 step when it was stopped. The check is a file test and
  runs before the lanes are joined, so the skip costs nothing.
- `assembleBinnedTargets()` writes
  `logs/assembleBinnedTargets_summary.csv`, one row per sample. The row is added
  as each sample finishes, so a batch that stops early keeps the rows it earned,
  and a rerun of one sample replaces its row. It records the targets recovered,
  the targets extended and the base pairs added, the bait sources, both rescue
  steps, the bins and targets of round 1, and the run time. The new
  `log.directory` argument sets where it goes. Default: `"logs"`.
- Every external command is checked before the first sample is processed, and a
  failure now stops the run with the command and its exit status.
- Read counts come from the fastp JSON report and from samtools instead of a
  second pass over each fastq file with gzip.
- `removeContamination()` maps, filters, counts, and writes the clean reads in
  one streaming pass. It no longer sorts by coordinate or keeps a BAM file.
- Resuming skips only the lanes that are complete. An interrupted lane is
  processed again instead of being skipped for good.
- fastp writes its HTML and JSON reports straight into `logs/sample_logs`. The
  JSON report is kept for tools such as MultiQC.
- `createContaminantDB()` downloads each genome once, checks that every download
  is FASTA, and continues an interrupted database.
- `sraDownload()` reads the ENA file report for the true file paths and their
  MD5 sums, so runs with an unusual file layout work and every download is
  verified.
- `fastqStats()` counts the read files in parallel and accepts uncompressed
  fastq files.
- `organizeReads()` can link the reads instead of copying them.
- `mergePairedEndReads()` records the number of merged reads.

# PhyloProcessR 1.0.0

## Initial public release

- Provides independently callable R functions and complete reference workflows
  for constructing reproducible, project-specific target-capture phylogenomic
  pipelines.
- Supports raw-read processing, assembly and target recovery, variant-aware
  consensus generation, alignment, quality control and filtering, paralog
  assessment, dataset construction, legacy-data integration, and recovery of
  novel or shared loci.
- Includes configurable workflow scripts with optional installation of the
  current GitHub development version; reproducible use of the installed package
  remains the default.
- Adds a redistributable seven-sample laptop example with single- and multilane
  libraries, a reduced 40-marker target panel, expected outputs, and provenance
  records.
- Adds versioned installation, configuration, workflow, assessment, and
  legacy-data tutorials in `docs/tutorials/`.
- Adds automated `testthat` regression coverage and GitHub Actions package
  checks.
