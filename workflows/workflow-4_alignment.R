#Installs updated package version
devtools::install_github("chutter/PhyloProcessR", upgrade = "never", dependencies = FALSE)
library(PhyloProcessR)

#source("/Users/chutter/Dropbox/Research/0_Github/R_Projects/PhyloProcessR/workflows/workflow-4_configuration-file.R")

source("workflow-4_configuration-file.R")
setwd(working.directory)

##################################################################################################
##################################################################################################
## Annotation and paralog filtering
##################################################################################################

#Remove contigs with too much heterozygosity
filterHeterozygosity(
  iupac.directory = contig.directory,
  output.directory = "data-analysis/contigs/7_filtered-contigs",
  removed.directory = "data-analysis/contigs/6_removed-contigs",
  threshold = heterozyote.filter.threshold,
  threads = threads,
  memory = memory,
  overwrite = overwrite
)

if (annotate.targets == TRUE) {
  # annotates targets
  annotateTargets(
    assembly.directory = "data-analysis/contigs/7_filtered-contigs",
    target.file = target.file,
    alignment.contig.name = paste0("data-analysis/", dataset.name),
    output.directory = "data-analysis/contigs/8_annotated-contigs",
    min.match.percent = min.match.percent,
    min.match.length = min.match.length,
    min.match.coverage = min.match.coverage,
    threads = threads,
    memory = memory,
    trim.target = trim.to.targets,
    overwrite = overwrite,
    quiet = quiet,
    blast.path = blast.path,
    cdhit.path = cdhit.path
  )
}#end if

# Create alignments folder
dir.create("data-analysis/alignments")

if (align.targets == TRUE) {
  # Aligns target markers from annotation files
  alignTargets(
    targets.to.align = paste0("data-analysis/", dataset.name, "_to-align.fa"),
    output.directory = "data-analysis/alignments/untrimmed_all-markers",
    min.taxa = min.taxa.alignment,
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