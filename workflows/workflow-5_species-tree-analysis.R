#Installs updated package version
devtools::install_github("chutter/PhyloProcessR", upgrade = "never", dependencies = FALSE)
library(PhyloProcessR)
library(foreach)

source("workflow-3_configuration-file.R")

##################################################################################################
##################################################################################################
#################################################
## Step 0: Pre-checks before running
##################

#Checks if everything is installed
pass.fail = setupCheck(anaconda.environment =  NULL,
                       iqtree.path = iqtree.path)

if (pass.fail == FALSE){ print("Some required programs are missing") } else {
  print("all required programs are found, PhyloProcessR pipeline continuing...")
}

##################################################################################################
## Step 1: Preprocess reads
##################################################################################################

setwd(working.directory)
dir.create("data-analysis/gene-trees")

## Step 1: estimate gene trees
##################

estimateGeneTrees(alignment.directory = alignment.directory,
                  output.directory = paste0("data-analysis/gene-trees/", dataset.name),
                  min.taxa = min.taxa.tree,
                  subset.start = 0,
                  subset.end = 1,
                  threads = threads,
                  memory = memory,
                  overwrite = overwrite,
                  quiet = quiet,
                  cleanup.files = cleanup.genetrees,
                  iqtree.path = iqtree.path)

