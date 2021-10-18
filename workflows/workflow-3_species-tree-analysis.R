#Installs updated package version
devtools::install_github("chutter/PhyloCap", upgrade = "never", force = TRUE)
library(PhyloCap)
library(foreach)

source("configuration-file.R")

##################################################################################################
##################################################################################################
#################################################
## Step 0: Pre-checks before running
##################

#Checks if everything is installed
pass.fail = setupCheck(anaconda.environment =  NULL,
                       fastp.path = fastp.path,
                       samtools.path = samtools.path,
                       bwa.path = bwa.path,
                       spades.path = spades.path,
                       bbmap.path = bbmap.path,
                       blast.path = blast.path,
                       mafft.path = mafft.path,
                       iqtree.path = iqtree.path,
                       trimAl.path = trimAl.path,
                       julia.path = julia.path,
                       taper.path = taper.path)

if (pass.fail == FALSE){ print("Some required programs are missing") } else {
  print("all required programs are found, PhyloCap pipeline continuing...")
}

##################################################################################################
## Step 1: Preprocess reads
##################################################################################################

setwd(work.dir)
dir.create("processed-reads")

## Step 1: estimate gene trees
##################

if (estimate.gene.trees == TRUE) {
  estimateGeneTrees(alignment.directory = "data-analysis/alignments-trimmed",
                    output.directory = "data-analysis/gene-trees",
                    min.taxa = min.taxa.tree,
                    subset.start = 0,
                    subset.end = 1,
                    threads = threads,
                    memory = memory,
                    overwrite = overwrite,
                    quiet = quiet,
                    cleanup.files = cleanup.genetrees,
                    iqtree.path = iqtree.path)

}#end

