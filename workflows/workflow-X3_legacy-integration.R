# Installs updated package version
devtools::install_github("chutter/PhyloProcessR", upgrade = "never", dependencies = FALSE)
library(PhyloProcessR)

#source("/Users/chutter/Dropbox/Research/0_Github/R_Projects/PhyloProcessR/workflows/workflow-X3_configuration-file.R")

source("workflow-X3_configuration-file.R")
setwd(working.directory)

##################################################################################################
##################################################################################################
## Step 0 (optional): Convert a concatenated NEXUS file into per-locus phylip files
##
## If your legacy data is a single NEXUS matrix with a BEGIN SETS / charset block
## (e.g. exported from PAUP*, MrBayes, or FigTree), set convert.nexus = TRUE and
## provide nexus.file. The file is split by charset into separate phylip files written
## to nexus.output.directory, which then becomes the legacy.directory for step 1.
##################################################################################################

if (convert.nexus == TRUE) {
  if (!dir.exists(nexus.output.directory) || overwrite == TRUE) {
    convertNexusPartitions(
      nexus.file          = nexus.file,
      output.directory    = nexus.output.directory,
      output.format       = "phylip",
      min.taxa.alignment  = min.taxa.alignment,
      max.missing.percent = max.missing.percent,
      overwrite           = overwrite,
      quiet               = quiet
    )
  } else {
    print(paste0("Nexus output directory already exists, skipping conversion: ",
                 nexus.output.directory))
  }
  # Use the converted output as the legacy alignment source
  legacy.directory = nexus.output.directory
  legacy.format    = "phylip"
}# end convert.nexus

##################################################################################################
##################################################################################################
## Step 1: Integrate legacy alignments into sequence-capture alignments
##################################################################################################

if (length(list.files("data-analysis/legacy-integration/untrimmed_legacy-only")) == 0 || overwrite == TRUE) {
  integrateLegacy(
    alignment.directory = alignment.directory,
    alignment.format = alignment.format,
    output.directory = "data-analysis/legacy-integration/untrimmed_legacy",
    legacy.directory = legacy.directory,
    legacy.format = legacy.format,
    target.markers = target.file,
    combine.same.sample = combine.same.sample,
    name.match = name.match,
    include.uncaptured.legacy = include.uncaptured.legacy,
    include.all.together = include.all.together,
    include.mitochondrial = include.mitochondrial,
    mito.alignment.directory = mito.alignment.directory,
    mito.alignment.format = mito.alignment.format,
    threads = threads,
    memory = memory,
    overwrite = overwrite,
    quiet = quiet,
    mafft.path = mafft.path,
    blast.path = blast.path
  )
} else {
  print("Legacy integration output already exists and is non-empty, skipping: data-analysis/legacy-integration/untrimmed_legacy-only")
}

# Select working directory for downstream steps:
# -all contains the full dataset (capture + legacy); -only contains only the integrated files
if (include.all.together == TRUE) {
  integrated.dir = "data-analysis/legacy-integration/untrimmed_legacy-all"
} else {
  integrated.dir = "data-analysis/legacy-integration/untrimmed_legacy-only"
}

##################################################################################################
##################################################################################################
## Step 2: Trim the legacy-only alignments
##
## Trims the raw legacy-integrated exon alignments (untrimmed_legacy-only) directly,
## producing a trimmed_legacy-only set. Useful for inspecting which legacy loci passed
## filters before downstream concatenation or gene-tree inference.
##################################################################################################

if (trim.alignments == TRUE) {
  superTrimmer(
    alignment.dir          = "data-analysis/legacy-integration/untrimmed_legacy-only",
    alignment.format       = "phylip",
    output.dir             = "data-analysis/legacy-integration/trimmed_legacy-only",
    overwrite              = overwrite,
    TrimAl                 = run.TrimAl,
    TrimAl.path            = trimAl.path,
    trim.similarity        = trim.similarity,
    similarity.threshold   = similarity.threshold,
    mafft.path             = mafft.path,
    trim.column            = trim.column,
    convert.ambiguous.sites = convert.ambiguous.sites,
    alignment.assess       = FALSE,
    trim.external          = trim.external,
    trim.coverage          = trim.coverage,
    min.coverage.percent   = min.coverage.percent,
    min.external.percent   = min.external.percent,
    min.column.gap.percent = min.column.gap.percent,
    min.alignment.length   = min.alignment.length,
    min.taxa.alignment     = min.taxa.alignment,
    min.coverage.bp        = min.coverage.bp,
    threads                = threads,
    memory                 = memory
  )
}# end trim.alignments

##################################################################################################
##################################################################################################
## Step 3: Concatenate exons into genes and gather unlinked dataset
##################################################################################################

if (concatenate.genes == TRUE) {

  if (concatenate.legacy.genes == TRUE) {
    # Concatenate genes from the full integrated dataset (capture + legacy)
    concatenateGenes(
      alignment.folder = integrated.dir,
      output.folder = "data-analysis/legacy-integration/untrimmed_legacy-genes",
      feature.gene.names = feature.gene.names,
      input.format = "phylip",
      output.format = "phylip",
      minimum.exons = minimum.exons,
      remove.reverse = FALSE,
      overwrite = overwrite,
      threads = threads,
      memory = memory
    )
  } else {
    # Concatenate genes from the original capture alignments only;
    # legacy loci remain as separate alignments and are picked up by gatherUnlinked
    concatenateGenes(
      alignment.folder = alignment.directory,
      output.folder = "data-analysis/legacy-integration/untrimmed_legacy-genes",
      feature.gene.names = feature.gene.names,
      input.format = "phylip",
      output.format = "phylip",
      minimum.exons = minimum.exons,
      remove.reverse = FALSE,
      overwrite = overwrite,
      threads = threads,
      memory = memory
    )
  }# end concatenate.legacy.genes

  if (gather.unlinked == TRUE) {
    # Exon directory is always the full integrated dataset so legacy stand-alone
    # loci are included regardless of whether they were concatenated
    gatherUnlinked(
      gene.alignment.directory = "data-analysis/legacy-integration/untrimmed_legacy-genes",
      exon.alignment.directory = integrated.dir,
      output.directory = "data-analysis/legacy-integration/untrimmed_legacy-unlinked",
      feature.gene.names = feature.gene.names,
      overwrite = overwrite
    )
  }# end gather.unlinked

  if (trim.alignments == TRUE) {
    superTrimmer(
      alignment.dir = "data-analysis/legacy-integration/untrimmed_legacy-unlinked",
      alignment.format = "phylip",
      output.dir = "data-analysis/legacy-integration/trimmed_legacy-unlinked",
      overwrite = overwrite,
      TrimAl = run.TrimAl,
      TrimAl.path = trimAl.path,
      trim.similarity = trim.similarity,
      similarity.threshold = similarity.threshold,
      mafft.path = mafft.path,
      trim.column = trim.column,
      convert.ambiguous.sites = convert.ambiguous.sites,
      alignment.assess = FALSE,
      trim.external = trim.external,
      trim.coverage = trim.coverage,
      min.coverage.percent = min.coverage.percent,
      min.external.percent = min.external.percent,
      min.column.gap.percent = min.column.gap.percent,
      min.alignment.length = min.alignment.length,
      min.taxa.alignment = min.taxa.alignment,
      min.coverage.bp = min.coverage.bp,
      threads = threads,
      memory = memory
    )
  }# end trim.alignments

}# end concatenate.genes

##################################################################################################
## If no gene concatenation is needed
##################################################################################################

if (concatenate.genes == FALSE) {

  if (trim.alignments == TRUE) {
    superTrimmer(
      alignment.dir = integrated.dir,
      alignment.format = "phylip",
      output.dir = "data-analysis/legacy-integration/trimmed_legacy",
      overwrite = overwrite,
      TrimAl = run.TrimAl,
      TrimAl.path = trimAl.path,
      trim.similarity = trim.similarity,
      similarity.threshold = similarity.threshold,
      mafft.path = mafft.path,
      trim.column = trim.column,
      convert.ambiguous.sites = convert.ambiguous.sites,
      alignment.assess = FALSE,
      trim.external = trim.external,
      trim.coverage = trim.coverage,
      min.coverage.percent = min.coverage.percent,
      min.external.percent = min.external.percent,
      min.column.gap.percent = min.column.gap.percent,
      min.alignment.length = min.alignment.length,
      min.taxa.alignment = min.taxa.alignment,
      min.coverage.bp = min.coverage.bp,
      threads = threads,
      memory = memory
    )
  }# end trim.alignments

}# end concatenate.genes == FALSE

### End workflow
