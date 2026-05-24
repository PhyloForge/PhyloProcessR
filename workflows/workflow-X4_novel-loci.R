# Installs updated package version
devtools::install_github("chutter/PhyloProcessR", upgrade = "never", dependencies = FALSE)
library(PhyloProcessR)

#source("/Users/chutter/Dropbox/Research/0_Github/R_Projects/PhyloProcessR/workflows/workflow-X4_configuration-file.R")

source("workflow-X4_configuration-file.R")
setwd(working.directory)

##################################################################################################
##################################################################################################
## Step 1: Discover novel genomic regions
##
## Builds a consensus reference from existing capture alignments, maps reads to it,
## extracts the unmapped reads, maps those to the reference genome, and identifies
## regions covered by reads in >= min.samples samples.
##
## Outputs to data-analysis/novel-loci-discovery/:
##   sample-bams/      — per-sample genome-mapped BAM files (input to step 2)
##   novel_regions.bed — coordinates of shared novel regions
##   novel_targets.fa  — genome sequences at those regions (target file for steps 4-5)
##################################################################################################

discoverSharedRegions(
  alignment.directory = alignment.directory,
  alignment.format    = alignment.format,
  read.directory      = read.directory,
  genome.file         = genome.file,
  output.directory    = "data-analysis/novel-loci-discovery",
  min.samples         = min.samples,
  min.coverage        = min.coverage,
  min.region.length   = min.region.length,
  max.merge.distance  = max.merge.distance,
  min.mapping.quality = min.mapping.quality,
  threads             = threads,
  memory              = memory,
  overwrite           = overwrite,
  quiet               = quiet,
  hisat2.path         = hisat2.path,
  samtools.path       = samtools.path,
  bedtools.path       = bedtools.path
)

##################################################################################################
##################################################################################################
## Step 2: Assemble per-sample contigs for each novel region
##
## For each sample and each discovered region, extracts the reads that mapped to that
## region and runs SPAdes to assemble contigs. One FASTA per sample is written to
## data-analysis/contigs/9_genome-contigs/.
##################################################################################################

if (file.exists("data-analysis/novel-loci-discovery/novel_regions.bed")) {

  assembleSharedRegions(
    discover.directory  = "data-analysis/novel-loci-discovery",
    output.directory    = "data-analysis/contigs/9_genome-contigs",
    min.reads.assemble  = min.reads.assemble,
    kmer.values         = spades.kmer.values,
    threads             = threads,
    memory              = memory,
    overwrite           = overwrite,
    quiet               = quiet,
    spades.path         = spades.path,
    samtools.path       = samtools.path,
    bedtools.path       = bedtools.path,
    blast.path          = blast.path
  )

}# end if

##################################################################################################
##################################################################################################
## Step 3: Filter contigs for heterozygosity
##
## Removes contigs with an excessive proportion of IUPAC ambiguity bases, which can
## indicate chimeric assembly or mis-assembled paralogs.
##################################################################################################

if (file.exists("data-analysis/contigs") == FALSE) { dir.create("data-analysis/contigs") }

if (heterozygote.filter == TRUE) {
  filterHeterozygosity(
    iupac.directory    = "data-analysis/contigs/9_genome-contigs",
    output.directory   = "data-analysis/contigs/10_genome-filtered",
    removed.directory  = "data-analysis/contigs/9b_genome-removed",
    threshold          = heterozygote.filter.threshold,
    min.length         = heterozygote.min.length,
    threads            = threads,
    memory             = memory,
    overwrite          = overwrite
  )
  input.contigs = "data-analysis/contigs/10_genome-filtered"
} else {
  input.contigs = "data-analysis/contigs/9_genome-contigs"
}

##################################################################################################
##################################################################################################
## Step 4: Collect novel contigs
##
## Because assembleSharedRegions names every contig by its source genomic region
## (e.g. chr3_450000_450800_contig_1), the locus assignment is already encoded in
## the name. collectNovelContigs groups by region, picks the longest contig per
## (region, sample) pair, and writes the to-align FASTA — no BLAST or probe set
## involved.
##################################################################################################

if (annotate.targets == TRUE) {
  collectNovelContigs(
    contig.directory  = input.contigs,
    output.name       = paste0("data-analysis/", dataset.name),
    min.contig.length = min.match.length,
    min.taxa          = min.taxa.alignment,
    threads           = threads,
    overwrite         = overwrite
  )
}# end annotate.targets

##################################################################################################
##################################################################################################
## Step 5: Align novel loci across samples
##
## Runs MAFFT to produce multi-sample alignments for each novel locus, using the
## novel_targets.fa as the reference. Output alignments can be merged with the existing
## untrimmed_all-markers alignments or kept separate for downstream trimming (workflow 5).
##################################################################################################

dir.create("data-analysis/alignments", showWarnings = FALSE)

if (align.targets == TRUE) {
  alignTargets(
    targets.to.align  = paste0("data-analysis/", dataset.name, "_to-align.fa"),
    target.file       = "data-analysis/novel-loci-discovery/novel_targets.fa",
    output.directory  = "data-analysis/alignments/untrimmed_novel-markers",
    min.taxa          = min.taxa.alignment,
    removal.threshold = removal.threshold,
    algorithm         = alignment.algorithm,
    threads           = threads,
    memory            = memory,
    overwrite         = overwrite,
    quiet             = quiet,
    mafft.path        = mafft.path
  )
}# end align.targets

### End workflow
