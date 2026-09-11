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

prepareCombinedAlignments = function(include.novel.markers = FALSE,
                                     overwrite = FALSE) {
  ordinary.directory = "data-analysis/alignments/untrimmed_all-markers"
  if (!include.novel.markers) return(ordinary.directory)

  novel.directory = "data-analysis/alignments/untrimmed_novel-markers"
  if (!dir.exists(novel.directory)) {
    stop("include.novel.markers is TRUE, but the novel-marker directory is missing: ",
         novel.directory)
  }
  pattern = "\\.(phy|phylip|fa|fas|fasta|nex|nexus)$"
  ordinary.files = list.files(ordinary.directory, pattern = pattern,
                              full.names = TRUE, ignore.case = TRUE)
  novel.files = list.files(novel.directory, pattern = pattern,
                           full.names = TRUE, ignore.case = TRUE)
  if (length(novel.files) == 0) stop("No novel-marker alignments were found.")
  markerId = function(path = NULL) {
    sub(pattern, "", basename(path), ignore.case = TRUE)
  }
  duplicate.ids = intersect(markerId(ordinary.files), markerId(novel.files))
  if (length(duplicate.ids) > 0) {
    stop("Duplicate marker IDs occur in ordinary and novel inputs: ",
         paste(utils::head(duplicate.ids, 5), collapse = ", "), ".")
  }

  combined.directory = "data-analysis/alignments/untrimmed_all-plus-novel"
  if (dir.exists(combined.directory)) unlink(combined.directory, recursive = TRUE)
  dir.create(combined.directory, recursive = TRUE)
  for (source in c(ordinary.files, novel.files)) {
    copied = file.copy(source, file.path(combined.directory, basename(source)),
                       overwrite = overwrite)
    if (!isTRUE(copied)) stop("Could not stage alignment: ", source)
  }
  combined.directory
}

copyUnmappedNovelMarkers = function(output.directory = NULL,
                                    metadata.file = NULL,
                                    overwrite = FALSE) {
  if (!include.novel.markers) return(invisible(NULL))
  metadata = data.table::fread(metadata.file, header = TRUE)
  if (!"marker" %in% names(metadata) && "Marker" %in% names(metadata)) {
    data.table::setnames(metadata, "Marker", "marker")
  }
  novel.directory = "data-analysis/alignments/untrimmed_novel-markers"
  pattern = "\\.(phy|phylip|fa|fas|fasta|nex|nexus)$"
  novel.files = list.files(novel.directory, pattern = pattern,
                           full.names = TRUE, ignore.case = TRUE)
  markerId = function(path = NULL) {
    sub(pattern, "", basename(path), ignore.case = TRUE)
  }
  unmapped = novel.files[!markerId(novel.files) %in% metadata$marker]
  for (source in unmapped) {
    destination = file.path(output.directory, basename(source))
    copied = file.copy(source, destination, overwrite = overwrite)
    if (!isTRUE(copied) && !(file.exists(destination) && !overwrite)) {
      stop("Could not copy novel alignment: ", source)
    }
  }
  message(length(unmapped), " unmapped novel marker(s) added individually.")
}

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
    # Fix the installs for this
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
## Create concatenated genes and unlinked datasets
##################################################################################################

if (concatenate.genes == TRUE) {
  exon.dir = prepareCombinedAlignments(include.novel.markers, overwrite)
  # Concatenates genes from the untrimmed-markers original alignments
  concatenateGenes(
    alignment.folder = exon.dir,
    output.folder = "data-analysis/alignments/untrimmed_genes",
    feature.gene.names = feature.gene.names,
    input.format = "phylip",
    output.format = "phylip",
    minimum.exons = minimum.exons,
    remove.reverse = FALSE,
    overwrite = overwrite,
    threads = threads,
    memory = memory
  )

  if (gather.unlinked == TRUE){
    # Gathers the unlinked markers (genes and single exons / UCEs)
    gatherUnlinked(
      gene.alignment.directory = "data-analysis/alignments/untrimmed_genes",
      exon.alignment.directory = exon.dir,
      output.directory = "data-analysis/alignments/untrimmed_all-unlinked",
      feature.gene.names = feature.gene.names,
      overwrite = overwrite
    )
    copyUnmappedNovelMarkers("data-analysis/alignments/untrimmed_all-unlinked",
                             feature.gene.names, overwrite)
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
  trim.input.dir = prepareCombinedAlignments(include.novel.markers, overwrite)

  # Trims the unlinked
  if (trim.alignments == TRUE) {
    superTrimmer(
      alignment.dir = trim.input.dir,
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
