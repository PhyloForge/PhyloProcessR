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
| 1 | Process raw reads | `organizeReads`, `fastqStats`, `fastpComplete`, `removeContamination`, `mergePairedEndReads` |
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

Workflow 3 maps reads to sample assemblies, calls variants, and makes consensus
contigs. It can produce IUPAC consensus contigs and separate haplotype contigs.
The workflow requires BWA, Samtools, and GATK.

```bash
Rscript workflow-3_variant-calling.R
```

Skip this workflow when draft contigs are sufficient for the project.

## Workflow 4: Annotate and align targets

Use `workflow-4_configuration-file.R` with `workflow-4_alignment.R`.

Set `contig.directory` to the contig set that you want to analyze. Set
`target.file` to the target-marker FASTA file.

Workflow 4 can filter contigs with high IUPAC ambiguity. It then matches
contigs to target loci and creates one alignment for each locus. The principal
output is `data-analysis/alignments/untrimmed_all-markers/`.

```bash
Rscript workflow-4_alignment.R
```

## Workflow 5: Trim and construct datasets

Use `workflow-5_configuration-file.R` with `workflow-5_trimming.R`.

Workflow 5 can perform these operations:

1. Trim alignments to target regions.
2. Extract flanking regions.
3. Refine coding alignments with MACSE.
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
