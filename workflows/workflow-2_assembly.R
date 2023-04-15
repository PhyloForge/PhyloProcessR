#Installs updated package version
devtools::install_github("chutter/PhyloProcessR", upgrade = "never", dependencies = FALSE)
library(PhyloProcessR)
library(foreach)

#setwd("/Users/chutter/Dropbox/Research/0_Github/Test-dataset")
#source("/Users/chutter/Dropbox/Research/0_Github/PhyloCap/work-flows/workflow-1_configuration-file.R")
source("workflow-2_configuration-file.R")
setwd(working.directory)

##################################################################################################
##################################################################################################
#################################################
## Step 1: Assemble reads
##################

#Begins by creating processed read directory
dir.create("data-analysis")

if (denovo.assembly == TRUE){

input.reads = paste0(processed.reads, "/", assembly.reads)

  #Assembles merged paired end reads with spades
  assembleSpades(input.reads = input.reads,
                 output.directory = paste0("data-analysis/spades-assembly"),
                 assembly.directory = "data-analysis/draft-assemblies",
                 mismatch.corrector = spades.mismatch.corrector,
                 kmer.values = spades.kmer.values,
                 threads = threads,
                 memory = memory,
                 overwrite = overwrite,
                 save.corrected.reads = save.corrected.reads,
                 quiet = quiet,
                 spades.path = spades.path)
}#end if

#End script
