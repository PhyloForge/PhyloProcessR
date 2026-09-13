source("workflow-4_configuration-file.R")
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
## Annotation and paralog filtering
##################################################################################################

if (file.exists("data-analysis/contigs") == FALSE){ dir.create("data-analysis/contigs") }

if (heterozygote.filter == TRUE){
  #Remove contigs with too much heterozygosity
  filterHeterozygosity(
    iupac.directory = contig.directory,
    output.directory = "data-analysis/contigs/7_filtered-contigs",
    removed.directory = "data-analysis/contigs/6_removed-contigs",
    threshold = heterozygote.filter.threshold,
    min.length = heterozygote.min.length,
    threads = threads,
    memory = memory,
    overwrite = overwrite
  )

  input.contigs = "data-analysis/contigs/7_filtered-contigs"

} else { input.contigs = contig.directory }

if (annotate.targets == TRUE) {
  # annotates targets
  annotateTargets(
    assembly.directory = input.contigs,
    target.file = target.file,
    alignment.contig.name = paste0("data-analysis/", dataset.name),
    output.directory = "data-analysis/contigs/8_annotated-contigs",
    min.match.percent = min.match.percent,
    min.match.length = min.match.length,
    min.match.coverage = min.match.coverage,
    paralog.action = paralog.action,
    paralog.score.ratio = paralog.score.ratio,
    paralog.coverage.ratio = paralog.coverage.ratio,
    paralog.identity.delta = paralog.identity.delta,
    paralog.directory = "data-analysis/contigs/9_paralog-contigs",
    threads = threads,
    memory = memory,
    overwrite = overwrite,
    quiet = quiet,
    blast.path = blast.path,
    last.path = last.path,
    search.method = annotate.search.method
  )
}#end if

# Create alignments folder
dir.create("data-analysis/alignments", recursive = TRUE, showWarnings = FALSE)

if (align.targets == TRUE) {
  # Aligns target markers from annotation files
  alignTargets(
    targets.to.align = paste0("data-analysis/", dataset.name, "_to-align.fa"),
    target.file = target.file,
    additional.sequence.directory = if (include.genomes == TRUE) {
      genome.target.directory
    } else {
      NULL
    },
    output.directory = "data-analysis/alignments/untrimmed_all-markers",
    min.taxa = min.taxa.alignment,
    removal.threshold = removal.threshold,
    algorithm = alignment.algorithm,
    subset.start = subset.start,
    subset.end = subset.end,
    threads = threads,
    memory = memory,
    overwrite = overwrite,
    quiet = quiet,
    mafft.path = mafft.path
  )
}#end if

# Legacy integration runs in workflow X3 (addLegacyAlignments). It is not part of
# this workflow.
