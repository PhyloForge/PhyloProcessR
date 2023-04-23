#Installs updated package version
devtools::install_github("chutter/PhyloProcessR", upgrade = "never", dependencies = FALSE)
library(PhyloProcessR)

#setwd("/Users/chutter/Dropbox/Research/0_Github/Test-dataset")
#source("/Users/chutter/Dropbox/Research/0_Github/PhyloCap/work-flows/workflow-1_configuration-file.R")
source("workflow-2_configuration-file.R")
setwd(working.directory)

##################################################################################################
##################################################################################################
#################################################
## Step 1: Assemble reads
##################

# Begins by creating processed read directory
dir.create("data-analysis")
dir.create("data-analysis/contigs")

# Assembles merged paired end reads with spades
assembleSpades(
  input.reads = paste0(processed.reads, "/", assembly.reads),
  output.directory = "data-analysis/spades-assembly-raw",
  assembly.directory = "data-analysis/contigs/1_draft-contigs",
  mismatch.corrector = spades.mismatch.corrector,
  kmer.values = spades.kmer.values,
  threads = threads,
  memory = memory,
  overwrite = overwrite,
  save.corrected.reads = save.corrected.reads,
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
  blast.path = blast.path,
  memory = memory,
  threads = threads,
  overwrite = overwrite,
  quiet = quiet
)
#End script
