# Run the PhyloProcessR workflows

The supplied workflows are reference analyses. You can run a complete script,
run selected sections, or call the same functions from a custom R script.

## Before you run a workflow

1. Activate the software environment.
2. Put the workflow script and its configuration file in the same directory.
3. Set `working.directory` and all input paths.
4. Set the external program paths.
5. Keep `overwrite = FALSE` until you intend to replace output.

Run a workflow from a terminal:

```bash
conda activate PhyloProcessR
Rscript workflow-1_preprocess.R
```

You can also open the script in an R development environment and run selected
function calls.

## Standard workflows

| Workflow | Purpose | Principal functions |
|---|---|---|
| 1 | Process raw reads | `organizeReads`, `fastqStats`, `fastpClean`, `removeContamination`, `mergePairedEndReads` |
| 2 | Assemble reads and recover targets | `assembleSpades`, `reduceRedundancy`, `removeOffTargetContigs`, `expandMissingAssembly` |
| 3 | Call variants and make consensus contigs | `prepareBAM`, `mapReferenceSample`, `haplotypeCaller`, `genotypeSamples`, `VCFtoContigs` |
| 4 | Annotate contigs and align targets | `filterHeterozygosity`, `annotateTargets`, `alignTargets` |
| 5 | Trim and construct datasets | `trimAlignmentTargets`, `alignMACSE`, `concatenateGenes`, `gatherUnlinked`, `superTrimmer` |

## Workflow 1: Process raw reads

Use `workflow-1_configuration-file.R` with `workflow-1_preprocess.R`.

Workflow 1 can use local reads, Dropbox reads, or NCBI SRA reads. Do not enable
Dropbox and SRA input at the same time.

The script can perform these operations:

1. Organize and rename read files.
2. Calculate FASTQ statistics.
3. estimate capture efficiency.
4. Remove adapters and duplicates with fastp.
5. Correct read errors and trim low-quality bases.
6. Remove reads that map to contaminant references.
7. Merge paired reads.

The selected output directories occur under `processed-reads/`. Keep the last
paired-read output for mapping. Keep the selected assembly input for Workflow
2.

Use this command:

```bash
Rscript workflow-1_preprocess.R
```

## Workflow 2: Assemble and recover targets

Use `workflow-2_configuration-file.R` with `workflow-2_assembly.R`.

Set `assembly.reads` to the applicable Workflow 1 output. The default
configuration uses `pe-merged-reads` for assembly.

Workflow 2 performs these operations:

1. Assemble each sample with SPAdes.
2. Reduce redundant contigs with CD-HIT-EST.
3. Keep contigs that match the target markers.
4. Optionally recover missing loci with `expandMissingAssembly`.

The main contig outputs occur under `data-analysis/contigs/`.

If `expand.missing = TRUE`, use paired, unmerged reads for `mapping.reads`.
Do not use `pe-merged-reads` for this parameter.

```bash
Rscript workflow-2_assembly.R
```

## Workflow 3: Call variants

Use `workflow-3_configuration-file.R` with
`workflow-3_variant-calling.R`.

Workflow 3 maps reads to each sample's assembly, calls variants, and makes
alternate-reference contigs. It can produce ordinary sequences and IUPAC codes
at supported heterozygous SNPs; these are not phased haplotypes or random allele
samples. Reference bases are retained outside passing variants. The workflow
requires BWA, Samtools, and GATK.

Variant hard filters and depth filters answer different questions. QD, QUAL,
SOR, FS, MQ, and rank-sum settings label individual variant records; QD 2 is not
a 2x coverage cutoff, QUAL is not per-base quality or genotype GQ, and a
high-depth variant does not rescue a failed hard filter. A rejected correction
leaves the assembly base unless the independent depth rule masks that position.

Depth processing then operates across the original full contig span, including
invariant and zero-coverage positions. In `site` mode, depths below
`min.site.depth` become N. In `mean` mode, sample-contigs below the full-span mean
are removed. `both` applies both rules, while `none` preserves the earlier
behavior. `max.n.proportion` optionally removes a contig after masking. Site
masking currently supports SNP output only because indels shift coordinates.

For example, with a 10x site cutoff, both a passing SNP and an invariant base at
2x become N, while at 12x the passing SNP is applied and the invariant assembly
base remains. A SNP failing QUAL at 12x leaves the assembly base. A 100-base
contig with 50 bases at 20x and 50 at zero has mean depth 10x: it passes a 10x
mean rule but has half its sequence masked by a 10x site rule.

The effective order is optional BQSR, genotyping and record hard filtering,
application of passing variants, depth masking, mean/N contig filtering, and
then workflow 4's heterozygosity and target filters. N is missing sequence;
IUPAC codes are heterozygosity, and the two are reported separately. Workflow
4's current heterozygosity denominator is total length (including Ns): a
100-base contig with 80 Ns and 10 IUPAC bases is 10% by that definition, but 50%
among its 20 non-N bases. Its target match coverage is alignment coverage (%),
not sequencing depth (x), and a depth-passing contig can still fail assignment.

```bash
Rscript workflow-3_variant-calling.R
```

Skip this workflow when draft contigs are sufficient for the project.

## Workflow 4: Annotate and align targets

Use `workflow-4_configuration-file.R` with `workflow-4_alignment.R`.

Set `contig.directory` to the contig set that you want to analyze. Set
`target.file` to the target-marker FASTA file.

Workflow 4 can filter contigs with high IUPAC ambiguity. It then matches
contigs to target loci without collapsing similar copies. A candidate must pass
the configured identity, supported-length, and target-coverage filters. The
default target-coverage floor remains 30 percent.

`paralog.action = "exclude"` removes a sample-target assignment when another
qualifying copy has at least 80 percent of the best score, at least 80 percent
of its target coverage, and is within five identity percentage points. The
`"best"` setting always keeps the top candidate. Both settings save all
qualifying copies from multi-copy targets in
`data-analysis/contigs/9_paralog-contigs/`. Candidate decisions and raw search
hits are written under `logs/sample_logs/` for later review.

The principal alignment output is
`data-analysis/alignments/untrimmed_all-markers/`. A sample sequence with no
comparable reference positions is removed before the final alignment is saved.

To include legacy samples, first run workflow 4 for the capture-only alignments,
and then run workflow X3 with `include.all.together = TRUE`. Run workflow 4
again with `align.targets = FALSE` and `include.legacy = TRUE`. Set
`legacy.alignment.directory` to the X3 `-all` output. Workflow 4 copies those
alignments into `untrimmed_all-markers`, where they replace the capture-only
version of each matching locus. Workflow 5 can then trim the capture and legacy
samples together.

Use `legacy.rename.file` when a legacy sample and its sequence-capture sample
have different names. CSV, TSV, TXT, XLS, and XLSX files are accepted. Put the
legacy name in the first column and the sequence-capture name in the second
column. The preferred headings are `Legacy_Name` and `SeqCap_Name`; workflow 4
uses the first two columns when the headings differ. When both names are present
in an alignment, workflow 4 merges the rows and keeps the
sequence-capture base at conflicting sites. When only the legacy name is
present, workflow 4 renames that row.

```r
legacy.rename.file = "data-analysis/legacy-name-map.tsv"
```

To add target sequences extracted from genome assemblies, set
`include.genomes = TRUE` and set `genome.target.directory` to the top-level
output directory from `extractGenomeTarget()`. Workflow 4 finds each genome's
`*_target-matches.fa` file recursively. It adds those sequences to the capture
sequences before MAFFT aligns each locus.

```r
include.genomes = TRUE
genome.target.directory = "data-analysis/genome-targets"
```

```bash
Rscript workflow-4_alignment.R
```

## Workflow 5: Trim and construct datasets

Use `workflow-5_configuration-file.R` with `workflow-5_trimming.R`.

Workflow 5 can perform these operations:

1. Trim alignments to target regions.
2. Extract flanking regions.
3. Optionally refine no-flank alignments with MACSE when all targets are coding.
4. Concatenate exons from the same gene.
5. Gather one unlinked alignment for each gene or marker.
6. Remove poor samples, columns, edges, or alignments.
7. Make a named alignment subset.

Common output directories include:

| Directory | Content |
|---|---|
| `untrimmed_all-markers` | Initial per-locus alignments |
| `untrimmed_genes` | Exons concatenated by gene |
| `untrimmed_all-unlinked` | Concatenated genes and single markers |
| `trimmed_all-unlinked` | Trimmed unlinked dataset |
| `untrimmed_no-flanks` | Target regions without flanks |
| `untrimmed_only-flanks` | Flanking regions without targets |

The `marker` and `gene` columns in the gene metadata file connect exon files to
genes. The capitalized aliases `Marker` and `Gene` are also accepted. Workflow 5
stops before replacing output when the metadata columns are missing, one marker
maps to conflicting genes, or no alignment names match the table.

When `include.novel.markers = TRUE`, the workflow requires the Workflow X4 novel
alignment directory. It rejects marker-name collisions and builds genes from the
combined ordinary and novel input. Novel markers that are absent from gene
metadata remain separate in the unlinked dataset. These markers are not thereby
shown to be biologically independent. Set `overwrite = TRUE` after changing
dataset composition so that derived outputs are rebuilt.

`trim.alignments` controls filtering of the full-marker dataset. The target-only
and flank-only construction steps use `trim.to.targets` and `trim.to.flanks`.
The subset step reads `subset.alignment.directory`. Its default directory exists
only after a per-marker trimming run, so a gene-based dataset must be selected
explicitly and matched with gene IDs.

The configured minimum taxa and final minimum length rules are exclusive: an
alignment at the configured value is rejected. Column trimming removes columns
at or above the configured gap percentage. For example, a value of 30 removes a
column with 30 percent gaps. Sample percentage coverage is measured against the
longest sample. Ambiguity conversion uses a deterministic A/T-priority mapping.
Some low-level helpers leave alignments with three or fewer taxa unchanged, but
the final alignment assessment still applies its configured thresholds.

MACSE runs only within `trim.to.targets` and receives every no-flank alignment.
Enable it only for targets known to be coding and in frame. Its `trimmed_exons`
and `trimmed_genes` directories are separate products and do not replace the
standard no-flank unlinked dataset.

```bash
Rscript workflow-5_trimming.R
```

Set `include.novel.markers = TRUE` only after Workflow X4 produces
`untrimmed_novel-markers`.

## Additional workflows

### Workflow X1: Joint genotyping

Workflow X1 maps all samples to a common consensus reference and makes a joint
VCF file. Use it when the analysis requires variants in a common coordinate
system.

```bash
Rscript workflow-X1_joint-genotype_VCF.R
```

### Workflow X2: Capture assessment

Workflow X2 processes one sample at a time. It calculates raw-read, cleaned-read,
and target-mapping statistics. This design limits temporary disk use for a
large sample set.

```bash
Rscript workflow-X2_capture-assessment.R
```

The final table is `logs/X2_capture-assessment_FINAL.csv`.

### Workflow X3: Legacy-data integration

Workflow X3 adds Sanger or GenBank sequences to capture alignments. It can also
convert a partitioned NEXUS matrix and integrate mitochondrial alignments.

See [Integrate legacy data](legacy-integration.md).

### Workflow X4: Recover novel shared loci

Workflow X4 finds regions that are captured in multiple samples but are absent
from the target set. It performs these main operations:

1. Find shared covered regions against a reference genome.
2. Assemble reads for each candidate region.
3. Filter and collect novel contigs.
4. Align and trim the novel loci.

The final per-locus alignments occur in
`data-analysis/alignments/untrimmed_novel-markers/`.

```bash
Rscript workflow-X4_novel-loci.R
```

After completion, set `include.novel.markers = TRUE` in the Workflow 5
configuration.

### Workflow X5: Assess and separate candidate copies

Workflow X5 starts from workflow 4 alignments and candidate-copy records. It
adds saved copies to the applicable untrimmed alignment. It also creates an
alignment for a saved-copy target that has no workflow 4 alignment. Existing
alignment rows stay fixed when MAFFT adds copies, and the workflow checks this
condition before it accepts the expanded alignment.

The workflow performs these operations:

1. Match each alignment row and saved sequence to an exact workflow 4 candidate.
2. Add recognizable copy labels and expand each untrimmed marker alignment.
3. Trim with unique biological samples as the occupancy units.
4. Infer one checked IQ-TREE gene tree for each informative marker.
5. Retain a marker, separate one supported two-group split, or exclude it with a reason.
6. Export reports and optional downstream gene and unlinked datasets.

Use `workflow-X5_configuration-file.R` with
`workflow-X5_paralog-analysis.R`. Test the proposed thresholds with a small
target subset before you run the full target union.

Locus directories, alignment files, and tree files use the locus name. Sequence
labels keep the workflow 4 sample name. If a sample has multiple copies for one
locus, the labels use `_1`, `_2`, and later consecutive suffixes. A sample with
one copy keeps its original name.

```bash
Rscript workflow-X5_paralog-analysis.R
```

Accepted alignments occur in
`data-analysis/paralog-analysis/trimmed_all-markers/`. A supported split writes
both `Target_copyA` and `Target_copyB`. These names are local labels for putative
loci. They do not confirm orthology or physical independence. Split markers do
not enter an unlinked dataset unless curated metadata permits that use.

The first version separates two groups only. It reports families with more
complex copy patterns as unresolved. When a supported split is not available,
the workflow keeps the best available copy for each sample. Candidate rank
selects the best copy. Aligned coverage and stable sequence metadata resolve
ties. It then trims the single-copy alignment again and applies the standard
alignment thresholds. The locus is retained only when this final alignment
passes. A tree split can also reflect species history, alleles, or assembly
error. Repeated samples on both sides supply the automatic copy evidence. A
split between disjoint sample sets is not automatic evidence of duplication.
Long branches in a locus with one sequence per sample do not receive a paralog
review flag. Reciprocal copy loss can be difficult to distinguish with these
data.

The main reports occur in `data-analysis/paralog-analysis/tables/`. They include
the target and copy maps, quality decisions, split evidence, group membership,
accepted markers, and one final outcome for each source target. Original
workflow 4 and workflow 5 directories are not changed.

## Resume or replace output

Use `overwrite = FALSE` to keep completed output and resume supported stages.
Use `overwrite = TRUE` only after you confirm the output paths. Some functions
remove and recreate their output directory when this value is `TRUE`.

## Record the analysis

Keep these items with the results:

- all configuration files;
- all custom R scripts;
- the PhyloProcessR version or Git commit;
- the Conda environment file or container tag;
- log files and summary tables.

Continue with [Assess the results](assess-results.md).
