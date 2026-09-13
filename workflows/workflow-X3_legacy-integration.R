source("workflow-X3_configuration-file.R")
if (isTRUE(get0("install.latest.github", ifnotfound = FALSE))) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    stop("Install the remotes package to use install.latest.github = TRUE.")
  }
  remotes::install_github("PhyloForge/PhyloProcessR", upgrade = "never",
                          dependencies = FALSE)
}
library(PhyloProcessR)
setwd(working.directory)

legacy.output.base = file.path(output.directory, "untrimmed_legacy")
legacy.only.directory = paste0(legacy.output.base, "-only")
legacy.all.directory = paste0(legacy.output.base, "-all")
legacy.trimmed.only.directory = file.path(output.directory, "trimmed_legacy-only")
legacy.gene.directory = file.path(output.directory, "untrimmed_legacy-genes")
legacy.trimmed.gene.directory = file.path(output.directory, "trimmed_legacy-genes")
legacy.unlinked.directory = file.path(output.directory, "untrimmed_legacy-unlinked")
legacy.trimmed.unlinked.directory = file.path(output.directory, "trimmed_legacy-unlinked")
legacy.trimmed.directory = file.path(output.directory, "trimmed_legacy")
# The integration summary is written only after a complete run, so its presence
# marks a finished integration for resume.
integration.summary = paste0(legacy.output.base, "-integration_summary.txt")

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
      nexus.file = nexus.file,
      output.directory = nexus.output.directory,
      output.format = "phylip",
      min.taxa.alignment = min.taxa.alignment,
      max.missing.percent = max.missing.percent,
      overwrite = overwrite,
      quiet = quiet
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

if (!file.exists(integration.summary) || overwrite == TRUE) {
  addLegacyAlignments(
    alignment.directory = alignment.directory,
    alignment.format = alignment.format,
    output.directory = legacy.output.base,
    legacy.directory = legacy.directory,
    legacy.format = legacy.format,
    target.markers = target.file,
    merge = merge,
    rename.file = rename.file,
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
  print(paste0("Legacy integration already completed (summary present), skipping: ",
               integration.summary))
}

# Select working directory for downstream steps:
# -all contains the full dataset (capture + legacy); -only contains only the integrated files
if (include.all.together == TRUE) {
  integrated.dir = legacy.all.directory
} else {
  integrated.dir = legacy.only.directory
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
  if (length(list.files(legacy.only.directory)) == 0) {
    print("Integration retained no legacy-only alignments; skipping the legacy-only trim step.")
  } else {
    superTrimmer(
      alignment.dir = legacy.only.directory,
      alignment.format = "phylip",
      output.dir = legacy.trimmed.only.directory,
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
  }
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
      output.folder = legacy.gene.directory,
      feature.gene.names = feature.gene.names,
      input.format = "phylip",
      output.format = "phylip",
      minimum.exons = minimum.exons,
      remove.reverse = FALSE,
      remove.duplicates = remove.duplicates,
      overwrite = overwrite,
      threads = threads,
      memory = memory
    )
  } else {
    # Concatenate genes from the original capture alignments only; legacy loci
    # remain as separate alignments and are picked up by gatherUnlinked.
    # Legacy sequences added to a locus whose gene is a multi-exon gene are
    # represented only by the capture-only gene here, so that legacy data does
    # not enter the gene dataset. Set concatenate.legacy.genes = TRUE to include
    # legacy sequences in the gene concatenation.
    print(paste0("concatenate.legacy.genes = FALSE: legacy sequences in multi-exon ",
                 "genes are not in the gene dataset; set TRUE to include them."))
    concatenateGenes(
      alignment.folder = alignment.directory,
      output.folder = legacy.gene.directory,
      feature.gene.names = feature.gene.names,
      input.format = "phylip",
      output.format = "phylip",
      minimum.exons = minimum.exons,
      remove.reverse = FALSE,
      remove.duplicates = remove.duplicates,
      overwrite = overwrite,
      threads = threads,
      memory = memory
    )
  }# end concatenate.legacy.genes

  if (gather.unlinked == TRUE) {
    # Exon directory is always the full integrated dataset so legacy stand-alone
    # loci are included regardless of whether they were concatenated
    gatherUnlinked(
      gene.alignment.directory = legacy.gene.directory,
      exon.alignment.directory = integrated.dir,
      output.directory = legacy.unlinked.directory,
      feature.gene.names = feature.gene.names,
      overwrite = overwrite
    )
  }# end gather.unlinked

  if (trim.alignments == TRUE) {
    # Trim the dataset the enabled stages actually produced: the unlinked set when
    # gathering is on, otherwise the concatenated gene set.
    if (gather.unlinked == TRUE) {
      trim.input  = legacy.unlinked.directory
      trim.output = legacy.trimmed.unlinked.directory
    } else {
      trim.input  = legacy.gene.directory
      trim.output = legacy.trimmed.gene.directory
    }
    if (length(list.files(trim.input)) == 0) {
      stop("No alignments to trim in ", trim.input,
           ". Check the output of the preceding concatenation/gathering stage.")
    }
    superTrimmer(
      alignment.dir = trim.input,
      alignment.format = "phylip",
      output.dir = trim.output,
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
    if (length(list.files(integrated.dir)) == 0) {
      stop("No alignments to trim in ", integrated.dir,
           ". The integration stage produced no usable output.")
    }
    superTrimmer(
      alignment.dir = integrated.dir,
      alignment.format = "phylip",
      output.dir = legacy.trimmed.directory,
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
