source("workflow-5_configuration-file.R")
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
## Trimming to targets
##################################################################################################

# Trims alignments to target sequence leaving out flanks
if (trim.to.targets == TRUE) {
  trimAlignmentTargets(
    alignment.directory = "data-analysis/alignments/untrimmed_all-markers",
    alignment.format = "phylip",
    target.file = target.file,
    target.direction = TRUE,
    output.directory = "data-analysis/alignments/untrimmed_no-flanks",
    min.alignment.length = min.alignment.length,
    min.taxa.alignment = min.taxa.alignment,
    threads = threads,
    memory = memory,
    overwrite = overwrite,
    mafft.path = mafft.path
  )

  if (run.macse == TRUE) {
    alignMACSE(
      alignment.folder = "data-analysis/alignments/untrimmed_no-flanks",
      output.folder = "data-analysis/alignments/trimmed_exons",
      alignment.format = "phylip",
      output.format = "phylip",
      macse.path = macse.path,
      genetic.code = macse.genetic.code,
      threads = threads,
      memory = memory,
      overwrite = overwrite,
      quiet = quiet
    )

    # Concatenates genes from the MACSE refined exons
    concatenateGenes(
      alignment.folder = "data-analysis/alignments/trimmed_exons",
      output.folder = "data-analysis/alignments/trimmed_genes",
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
  }

  # Concatenates genes from the untrimmed-markers original alignments
  concatenateGenes(
    alignment.folder = "data-analysis/alignments/untrimmed_no-flanks",
    output.folder = "data-analysis/alignments/untrimmed_genes_no-flanks",
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

  # Gathers the unlinked markers (genes and single exons / UCEs)
  gatherUnlinked(
    gene.alignment.directory = "data-analysis/alignments/untrimmed_genes_no-flanks",
    exon.alignment.directory = "data-analysis/alignments/untrimmed_no-flanks",
    output.directory = "data-analysis/alignments/untrimmed_no-flanks-unlinked",
    feature.gene.names = feature.gene.names,
    overwrite = overwrite
  )

  if (trim.alignments == TRUE) {
    superTrimmer(
      alignment.dir = "data-analysis/alignments/untrimmed_no-flanks-unlinked",
      alignment.format = "phylip",
      output.dir = "data-analysis/alignments/trimmed_no-flanks-unlinked",
      overwrite = overwrite,
      TrimAl = run.TrimAl,
      TrimAl.path = trimAl.path,
      trim.similarity = trim.similarity,
      similarity.threshold = similarity.threshold,
      mafft.path = mafft.path,
      trim.column = trim.column,
      convert.ambiguous.sites = convert.ambiguous.sites,
      alignment.assess = alignment.assess,
      trim.external = trim.external,
      trim.coverage = trim.coverage,
      min.coverage.percent = min.coverage.percent,
      min.external.percent = min.external.percent,
      min.column.gap.percent = min.column.gap.percent,
      min.alignment.length = min.alignment.length,
      min.taxa.alignment = min.taxa.alignment,
      max.alignment.gap.percent = max.alignment.gap.percent,
      min.coverage.bp = min.coverage.bp,
      threads = threads,
      memory = memory
    )
  } # end if

}# end trim to targets

##################################################################################################
##################################################################################################
## Trimming to flanks
##################################################################################################

# Trims alignments to only the flanks, leaving out the target marker
if (trim.to.flanks == TRUE) {
  # Trim out the target region leaving only the flanks.
  makeFlankAlignments(
    alignment.directory = "data-analysis/alignments/untrimmed_all-markers",
    alignment.format = "phylip",
    output.directory = "data-analysis/alignments/untrimmed_only-flanks",
    reference.type = "target",
    reference.path = target.file,
    target.direction = TRUE,
    concatenate.intron.flanks = TRUE,
    threads = threads,
    memory = memory,
    overwrite = overwrite,
    mafft.path = mafft.path
  )

  # Concatenates genes from the untrimmed-markers original alignments
  concatenateGenes(
    alignment.folder = "data-analysis/alignments/untrimmed_only-flanks",
    output.folder = "data-analysis/alignments/untrimmed_genes_only-flanks",
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

  # Gathers the unlinked markers (genes and single exons / UCEs)
  gatherUnlinked(
    gene.alignment.directory = "data-analysis/alignments/untrimmed_genes_only-flanks",
    exon.alignment.directory = "data-analysis/alignments/untrimmed_only-flanks",
    output.directory = "data-analysis/alignments/untrimmed_only-flanks-unlinked",
    feature.gene.names = feature.gene.names,
    overwrite = overwrite
  )

  if (trim.alignments == TRUE) {
    superTrimmer(
      alignment.dir = "data-analysis/alignments/untrimmed_only-flanks-unlinked",
      alignment.format = "phylip",
      output.dir = "data-analysis/alignments/trimmed_only-flanks-unlinked",
      overwrite = overwrite,
      TrimAl = run.TrimAl,
      TrimAl.path = trimAl.path,
      trim.similarity = trim.similarity,
      similarity.threshold = similarity.threshold,
      mafft.path = mafft.path,
      trim.column = trim.column,
      convert.ambiguous.sites = convert.ambiguous.sites,
      alignment.assess = alignment.assess,
      trim.external = trim.external,
      trim.coverage = trim.coverage,
      min.coverage.percent = min.coverage.percent,
      min.external.percent = min.external.percent,
      min.column.gap.percent = min.column.gap.percent,
      min.alignment.length = min.alignment.length,
      min.taxa.alignment = min.taxa.alignment,
      max.alignment.gap.percent = max.alignment.gap.percent,
      min.coverage.bp = min.coverage.bp,
      threads = threads,
      memory = memory
    )
  } # end if
}# end trim to flanks

##################################################################################################
##################################################################################################
## Select the exon alignments to concatenate and trim
##################################################################################################

# When novel markers from workflow X4 are included, merge them with the
# target-marker alignments into one directory first, then work from that.
if (include.novel.markers == TRUE) {
  combineMarkerAlignments(
    alignment.directories = c("data-analysis/alignments/untrimmed_all-markers",
                              "data-analysis/alignments/untrimmed_novel-markers"),
    output.directory = "data-analysis/alignments/untrimmed_all-plus-novel",
    overwrite = overwrite
  )
  marker.input = "data-analysis/alignments/untrimmed_all-plus-novel"
} else {
  marker.input = "data-analysis/alignments/untrimmed_all-markers"
}

##################################################################################################
##################################################################################################
## Create concatenated genes and unlinked datasets
##################################################################################################

if (concatenate.genes == TRUE) {
  # Concatenates genes from the untrimmed-markers original alignments
  concatenateGenes(
    alignment.folder = marker.input,
    output.folder = "data-analysis/alignments/untrimmed_genes",
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

  if (gather.unlinked == TRUE){
    # Gathers the unlinked markers (genes and single exons / UCEs)
    gatherUnlinked(
      gene.alignment.directory = "data-analysis/alignments/untrimmed_genes",
      exon.alignment.directory = marker.input,
      output.directory = "data-analysis/alignments/untrimmed_all-unlinked",
      feature.gene.names = feature.gene.names,
      overwrite = overwrite
    )

    # Novel markers without a gene assignment are added back individually
    if (include.novel.markers == TRUE) {
      copyUnmatchedMarkers(
        alignment.directory = "data-analysis/alignments/untrimmed_novel-markers",
        output.directory = "data-analysis/alignments/untrimmed_all-unlinked",
        feature.gene.names = feature.gene.names,
        overwrite = overwrite
      )
    }
  }#end if

  # Trims the unlinked
  if (trim.alignments == TRUE && gather.unlinked == TRUE) {
    superTrimmer(
      alignment.dir = "data-analysis/alignments/untrimmed_all-unlinked",
      alignment.format = "phylip",
      output.dir = "data-analysis/alignments/trimmed_all-unlinked",
      overwrite = overwrite,
      TrimAl = run.TrimAl,
      TrimAl.path = trimAl.path,
      trim.similarity = trim.similarity,
      similarity.threshold = similarity.threshold,
      mafft.path = mafft.path,
      trim.column = trim.column,
      convert.ambiguous.sites = convert.ambiguous.sites,
      alignment.assess = alignment.assess,
      trim.external = trim.external,
      trim.coverage = trim.coverage,
      min.coverage.percent = min.coverage.percent,
      min.external.percent = min.external.percent,
      min.column.gap.percent = min.column.gap.percent,
      min.alignment.length = min.alignment.length,
      min.taxa.alignment = min.taxa.alignment,
      max.alignment.gap.percent = max.alignment.gap.percent,
      min.coverage.bp = min.coverage.bp,
      threads = threads,
      memory = memory
    )
  }
} # end concatenate genes


##################################################################################################
## If no concatenated genes are needed
##################################################################################################

if (concatenate.genes == FALSE) {
  # Trims the selected marker alignments directly
  if (trim.alignments == TRUE) {
    superTrimmer(
      alignment.dir = marker.input,
      alignment.format = "phylip",
      output.dir = "data-analysis/alignments/trimmed_all-markers",
      overwrite = overwrite,
      TrimAl = run.TrimAl,
      TrimAl.path = trimAl.path,
      trim.similarity = trim.similarity,
      similarity.threshold = similarity.threshold,
      mafft.path = mafft.path,
      trim.column = trim.column,
      convert.ambiguous.sites = convert.ambiguous.sites,
      alignment.assess = alignment.assess,
      trim.external = trim.external,
      trim.coverage = trim.coverage,
      min.coverage.percent = min.coverage.percent,
      min.external.percent = min.external.percent,
      min.column.gap.percent = min.column.gap.percent,
      min.alignment.length = min.alignment.length,
      min.taxa.alignment = min.taxa.alignment,
      max.alignment.gap.percent = max.alignment.gap.percent,
      min.coverage.bp = min.coverage.bp,
      threads = threads,
      memory = memory
    )
  }
} # end concatenate genes


##################################################################################################
## Alignment subset
##################################################################################################

if (run.subset == TRUE) {
  if (!dir.exists(subset.alignment.directory)) {
    stop("The configured subset alignment directory does not exist: ",
         subset.alignment.directory,
         ". Check concatenate.genes and trim.alignments.")
  }
  makeAlignmentSubset(
    alignment.directory = subset.alignment.directory,
    alignment.format = "phylip",
    output.directory = paste0("data-analysis/alignments/", subset.name),
    subset.reference = subset.reference,
    subset.fasta.file = subset.fasta,
    subset.grep.string = subset.grep.string,
    subset.blast.targets = target.file,
    blast.path = blast.path,
    threads = threads,
    memory = memory,
    overwrite = overwrite
  )
} # end subset


### End workflow
