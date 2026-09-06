source("workflow-2_configuration-file.R")
if (isTRUE(get0("install.latest.github", ifnotfound = FALSE))) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    stop("Install the remotes package to use install.latest.github = TRUE.")
  }
  remotes::install_github("PhyloForge/PhyloProcessR", upgrade = "never",
                          dependencies = FALSE)
}
library(PhyloProcessR)
setwd(working.directory)

##################################################################################################
##################################################################################################
#################################################
## Step 1: Assemble reads
##################

# Begins by creating processed read directory
dir.create("data-analysis", showWarnings = FALSE)
dir.create("data-analysis/contigs", showWarnings = FALSE)

# Assembles merged paired end reads with spades
assembleSpades(
  input.reads = paste0(processed.reads, "/", assembly.reads),
  output.directory = "data-analysis/spades-assembly-raw",
  assembly.directory = "data-analysis/contigs/1_draft-contigs",
  mismatch.corrector = spades.mismatch.corrector,
  isolate = spades.isolate,
  kmer.values = spades.kmer.values,
  threads = threads,
  memory = memory,
  overwrite = overwrite,
  save.corrected.reads = save.corrected.reads,
  clean.up.spades = clean.up.spades,
  quiet = quiet,
  spades.path = spades.path
)

#Reduces contig redundancy by removing the shortest contig in a set of similar contigs using cd-hit-est
reduceRedundancy(
  assembly.directory = "data-analysis/contigs/1_draft-contigs",
  output.directory = "data-analysis/contigs/2_reduced-redundancy",
  similarity = similarity,
  cdhit.path = cdhit.path,
  memory = memory,
  threads = threads,
  overwrite = overwrite,
  quiet = quiet
)

removeOffTargetContigs(
  assembly.directory = "data-analysis/contigs/2_reduced-redundancy",
  target.markers = target.markers,
  output.directory = "data-analysis/contigs/3_target-contigs",
  min.match.length = target.match.length,
  min.match.percent = target.match.percent,
  min.match.coverage = target.match.coverage,
  blast.path = blast.path,
  last.path = last.path,
  search.method = target.search.method,
  memory = memory,
  threads = threads,
  overwrite = overwrite,
  quiet = quiet
)

##################################################################################################
##################################################################################################
#################################################
## Step 1b: Curate the target contigs
## Joins the fragments of one target that sit on separate contigs, and cuts apart
## a contig that spans more than one target. This runs here as well as at the end
## because the steps below use these contigs as baits. A chimeric bait recruits
## the reads of two loci into one bin, and a fragment of a target only recruits
## across that fragment.
##################

if (isTRUE(get0("curate.contigs", ifnotfound = TRUE))) {
  curateTargetContigs(
    assembly.directory = "data-analysis/contigs/3_target-contigs",
    target.file        = target.markers,
    output.directory   = "data-analysis/contigs/3a_curated-contigs",
    min.match.percent  = curate.match.percent,
    min.match.length   = curate.match.length,
    min.match.coverage = curate.match.coverage,
    search.method      = curate.search.method,
    threads            = threads,
    memory             = memory,
    blast.path         = blast.path,
    last.path          = last.path,
    cdhit.path         = cdhit.path,
    overwrite          = overwrite,
    quiet              = quiet
  )
}#end curate.contigs

# The steps below start from the curated contigs when that step ran
contig.start = "data-analysis/contigs/3_target-contigs"
if (isTRUE(get0("curate.contigs", ifnotfound = TRUE))) {
  contig.start = "data-analysis/contigs/3a_curated-contigs"
}

##################################################################################################
##################################################################################################
#################################################
## Step 2 (optional): Recover missing loci
## Attempts to assemble loci absent from each sample by mapping reads to
## cross-sample contigs and/or the original reference, then running SPAdes.
## Enable with expand.missing = TRUE in the configuration file.
## The recovered contigs are saved to 3b_expanded-contigs. To use them, set
## expanded.contig.directory in the workflow 3 configuration file.
##################

if (isTRUE(get0("expand.missing", ifnotfound = FALSE))) {
  expandMissingAssembly(
    assembly.directory  = contig.start,
    read.directory      = processed.reads,
    mapping.reads       = mapping.reads,
    reference           = target.markers,
    output.directory    = "data-analysis/expand-missing-assembly",
    expanded.directory  = "data-analysis/contigs/3b_expanded-contigs",
    phase2.reference    = phase2.reference,
    recover.all.missing = recover.all.missing,
    min.match.percent   = expand.match.percent,
    min.match.length    = expand.match.length,
    min.match.coverage  = expand.match.coverage,
    mismatch.corrector  = spades.mismatch.corrector,
    kmer.values         = spades.kmer.values,
    memory              = memory,
    threads             = threads,
    spades.path         = spades.path,
    hisat2.path         = hisat2.path,
    samtools.path       = samtools.path,
    fastp.path          = fastp.path,
    blast.path          = blast.path,
    overwrite           = overwrite,
    quiet               = quiet
  )
} # end expand.missing

##################################################################################################
##################################################################################################
#################################################
## Step 3 (optional): Binned per-locus assembly
## Assembles every target on its own. Reads are binned by their best matching
## target, and each bin is assembled separately. This recovers targets the
## whole-library assembly lost, and it extends the targets it found, because the
## mate of an anchored read reaches into the flanking sequence.
## Enable with binned.assembly = TRUE in the configuration file.
## The merged contigs are saved to 3c_binned-contigs. To use them, set
## expanded.contig.directory in the workflow 3 configuration file.
##################

if (isTRUE(get0("binned.assembly", ifnotfound = FALSE))) {

  # Builds on the expanded contigs when step 2 ran, on the target contigs if not
  binned.input = contig.start
  if (isTRUE(get0("expand.missing", ifnotfound = FALSE))) {
    binned.input = "data-analysis/contigs/3b_expanded-contigs"
  }

  # An empty setting turns off the LAST search of the draft assembly
  binned.draft = get0("binned.draft.directory", ifnotfound = "")
  if (is.null(binned.draft) || nchar(binned.draft) == 0) binned.draft = NULL

  assembleBinnedTargets(
    read.directory     = processed.reads,
    mapping.reads      = binning.reads,
    target.markers     = target.markers,
    assembly.directory = binned.input,
    draft.assembly.directory = binned.draft,
    output.directory   = "data-analysis/binned-target-assembly",
    binned.directory   = "data-analysis/contigs/3c_binned-contigs",
    locus.set          = binned.locus.set,
    bait.source        = binned.bait.source,
    min.bait.coverage  = binned.min.bait.coverage,
    rescue.missing     = binned.rescue.missing,
    iterations         = binned.iterations,
    min.pairs          = binned.min.pairs,
    max.pairs          = binned.max.pairs,
    min.match.length   = binned.match.length,
    min.match.percent  = binned.match.percent,
    min.match.coverage = binned.match.coverage,
    max.extension      = binned.max.extension,
    max.target.hits    = binned.max.target.hits,
    multi.copy         = binned.multi.copy,
    kmer.values        = binned.kmer.values,
    memory             = memory,
    threads            = threads,
    bwa.path           = bwa.path,
    samtools.path      = samtools.path,
    megahit.path       = megahit.path,
    cap3.path          = cap3.path,
    last.path          = last.path,
    overwrite          = overwrite,
    quiet              = quiet
  )
} # end binned.assembly

##################################################################################################
##################################################################################################
#################################################
## Step 4 (optional): Curate the contigs
## Joins the fragments of one target that sit on separate contigs, and cuts
## apart a contig that spans more than one target. This runs before variant
## calling, because the variant caller maps the reads back to these contigs. A
## contig that spans two targets collects the reads of both loci and gives wrong
## genotypes, and no later step repairs that.
## Enable with curate.contigs = TRUE in the configuration file.
## The curated contigs are saved to 3d_curated-contigs. Use them in workflow 3.
##################

if (isTRUE(get0("curate.contigs", ifnotfound = FALSE))) {

  # Takes the last contig set that the steps above produced
  curate.input = contig.start
  if (isTRUE(get0("expand.missing", ifnotfound = FALSE))) {
    curate.input = "data-analysis/contigs/3b_expanded-contigs"
  }
  if (isTRUE(get0("binned.assembly", ifnotfound = FALSE))) {
    curate.input = "data-analysis/contigs/3c_binned-contigs"
  }

  curateTargetContigs(
    assembly.directory = curate.input,
    target.file        = target.markers,
    output.directory   = "data-analysis/contigs/3d_curated-contigs",
    min.match.percent  = curate.match.percent,
    min.match.length   = curate.match.length,
    min.match.coverage = curate.match.coverage,
    search.method      = curate.search.method,
    threads            = threads,
    memory             = memory,
    blast.path         = blast.path,
    last.path          = last.path,
    cdhit.path         = cdhit.path,
    overwrite          = overwrite,
    quiet              = quiet
  )
} # end curate.contigs

#End script
