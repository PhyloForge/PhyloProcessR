#Installs updated package version
devtools::install_github("chutter/PhyloCap", upgrade = "never", dependencies = FALSE)
library(PhyloCap)
#Load in packages
devtools::install_github("chutter/AstralPlane", upgrade = "never", force = TRUE)
library(AstralPlane)
#Load in packages
devtools::install_github("chutter/FilterZone", upgrade = "never", force = TRUE)
library(FilterZone)
library(foreach)

#working.directory = "/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/Guibemantis-liber"
#setwd(working.directory)

#source("/Users/chutter/Dropbox/Research/0_Github/PhyloCap/work-flows/workflow-2_configuration-file.R")
source("workflow-4_configuration-file.R")

##################################################################################################
##################################################################################################
#################################################
## Step 0: Pre-checks before running
##################

#Checks if everything is installed
pass.fail = setupCheck(anaconda.environment =  NULL,
                       iqtree.path = iqtree.path)

if (pass.fail == FALSE){ print("Some required programs are missing") } else {
  print("all required programs are found, PhyloCap pipeline continuing...")
}

##################################################################################################
## Step 1: Preprocess reads
##################################################################################################

setwd(working.directory)
dir.create("data-analysis/concatenation")

## Step 1: estimate gene trees
##################
#Alignment summary function [20 minutes]
align.summary = FilterZone::summarizeAlignments(alignment.path = alignment.directory,
                                                file.export = "alignment_stats",
                                                alignment.format = "phylip",
                                                dataset.name = dataset.name,
                                                overwrite = overwrite)

#Apply filters and create summary table of filters [< 1 min]
filt.summary = FilterZone::filterSummary(alignment.data = align.summary,
                                         alignment.folder = alignment.directory,
                                         dataset.name  = dataset.name,
                                         file.out = "filter_summary",
                                         length.filters = filter.length,
                                         sample.filters = filter.sample,
                                         prop.pis.filters = filter.prop.pis,
                                         count.pis.filters = filter.count.pis,
                                         overwrite = overwrite)

#Make filtered alignments datasets [10 minutes]
FilterZone::filterAlignments(filter.summary = filt.summary,
                             alignment.data = align.summary,
                             alignment.folder = alignment.directory,
                             format = "concatenated",
                             min.alignments = min.alignments,
                             overwrite = overwrite)

#Gathers alignment files and runs a tree for each
alignment.files = list.files("filtered-alignments-concatenated", full.names = T)
alignment.files = alignment.files[grep(".phy$", alignment.files)]
out.path = paste0(working.directory, "/data-analysis/concatenation-trees")
dir.create(out.path)

### Somehow make a loop that submits jobs all at once
for (i in 1:length(alignment.files)){

  align.name = gsub(".phy$", "", alignment.files[i])
  align.name = gsub(".*\\/", "", align.name)

  analysis.concatenationTree(alignment.file = alignment.files[i],
                             output.directory = out.path,
                             output.name = align.name,
                             partition.file = NULL,
                             partition.scheme = partition.scheme,
                             codon.partition = FALSE,
                             program = "IQTREE",
                             msub.type = "nuclear",
                             uf.bootstrap = uf.bootstrap,
                             rcluster = rcluster,
                             threads = threads,
                             memory = memory,
                             iqtree.path = iqtree.path,
                             overwrite = overwrite)

}#end i loop


