#Installs updated package version
devtools::install_github("chutter/PhyloProcessR", upgrade = "never", dependencies = FALSE)
library(PhyloProcessR)

#source("/Users/chutter/Dropbox/Research/0_Github/R_Projects/PhyloProcessR/workflows/workflow-5_configuration-file.R")
source("workflow-5_configuration-file.R")
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
    output.format = "phylip",
    min.alignment.length = min.alignment.length,
    min.taxa.alignment = 4,
    threads = threads,
    memory = memory,
    overwrite = overwrite,
    mafft.path = mafft.path
  )

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
    # Fix the installs for this
    superTrimmer(
      alignment.dir = "data-analysis/alignments/untrimmed_no-flanks-unlinked",
      alignment.format = "phylip",
      output.dir = "data-analysis/alignments/trimmed_no-flanks-unlinked",
      output.format = "phylip",
      overwrite = overwrite,
      TrimAl = FALSE,
      TrimAl.path = trimAl.path,
      trim.column = FALSE,
      convert.ambiguous.sites = convert.ambiguous.sites,
      alignment.assess = FALSE,
      trim.external = FALSE,
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
    output.format = "phylip",
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
      output.format = "phylip",
      overwrite = overwrite,
      TrimAl = run.TrimAl,
      TrimAl.path = trimAl.path,
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
  } # end if
}# end trim to flanks

##################################################################################################
##################################################################################################
## Create concatenated genes and unlinked datasets
##################################################################################################

if (concatenate.genes == TRUE) {
  # Concatenates genes from the untrimmed-markers original alignments
  concatenateGenes(
    alignment.folder = "data-analysis/alignments/untrimmed_all-markers",
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
      exon.alignment.directory = "data-analysis/alignments/untrimmed_all-markers",
      output.directory = "data-analysis/alignments/untrimmed_all-unlinked",
      feature.gene.names = feature.gene.names,
      overwrite = overwrite
    )
  }#end if

  # Trims the unlinked
  superTrimmer(
    alignment.dir = "data-analysis/alignments/untrimmed_all-unlinked",
    alignment.format = "phylip",
    output.dir = "data-analysis/alignments/trimmed_all-unlinked",
    output.format = "phylip",
    overwrite = overwrite,
    TrimAl = run.TrimAl,
    TrimAl.path = trimAl.path,
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
} # end concatenate genes


##################################################################################################
## If no concatenated genes are needed
##################################################################################################

if (concatenate.genes == FALSE) {
  # Trims the unlinked
  superTrimmer(
    alignment.dir = "data-analysis/alignments/untrimmed_all-markers",
    alignment.format = "phylip",
    output.dir = "data-analysis/alignments/trimmed_all-markers",
    output.format = "phylip",
    overwrite = overwrite,
    TrimAl = run.TrimAl,
    TrimAl.path = trimAl.path,
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
} # end concatenate genes


# makeAlignmentSubset(alignment.directory = "data-analysis/alignments/untrimmed_all-markers",
#                     alignment.format = "phylip",
#                     output.directory = "data-analysis/alignments/untrimmed_UCE",
#                     output.format = "phylip",
#                     subset.reference = "blast",
#                     subset.fasta.file = "subset_uce-consensus.fa",
#                     subset.grep.string = NULL,
#                     subset.blast.targets = target.file,
#                     blast.path = blast.path,
#                     threads = threads,
#                     memory = memory,
#                     overwrite = overwrite)


### End workflow
