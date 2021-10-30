#Installs updated package version
devtools::install_github("chutter/PhyloCap", upgrade = "never", dependencies = FALSE)
library(PhyloCap)
library(foreach)

source("workflow-2_configuration-file.R")
setwd(work.directory)

##################################################################################################
##################################################################################################
#################################################
## Step 1: Preprocess reads
##################

#Checks if everything is installed
pass.fail = setupCheck(anaconda.environment =  NULL,
                       blast.path = blast.path,
                       mafft.path = mafft.path,
                       trimAl.path = trimAl.path,
                       julia.path = julia.path,
                       taper.path = taper.path)

if (pass.fail == FALSE){ print("Some required programs are missing") } else {
  print("all required programs are found, PhyloCap pipeline continuing...")
}


if (match.targets == TRUE){
  #match targets
  matchTargets(assembly.directory = "data-analysis/draft-assemblies",
               target.file = target.file,
               alignment.contig.name = paste0("data-analysis/", dataset.name),
               output.directory = "data-analysis/match-targets",
               min.match.percent = min.match.percent,
               min.match.length = min.match.length,
               min.match.coverage = min.match.coverage,
               threads = threads,
               memory = memory,
               trim.target = trim.to.targets,
               overwrite = overwrite,
               quiet = quiet,
               blast.path = blast.path,
               bbmap.path = bbmap.path)
}#end if

#################################################
## Step 3: Align targets and trim
##################

dir.create("data-analysis/alignments")

if (align.matching.targets == TRUE){
  #align targets
  alignTargets(targets.to.align = paste0("data-analysis/", dataset.name, "_to-align.fa"),
               output.directory = "data-analysis/alignments/untrimmed_all-markers",
               min.taxa = min.taxa.alignment,
               subset.start = 0,
               subset.end = 1,
               threads = threads,
               memory = memory,
               overwrite = overwrite,
               quiet = quiet,
               mafft.path = mafft.path)
}#end if

#### Functions for separating into data types
trimAlignmentTargets(alignment.directory = "data-analysis/alignments/untrimmed_all-markers",
                     alignment.format = "phylip",
                     target.file = target.file,
                     target.direction = TRUE,
                     output.directory = "data-analysis/alignments/untrimmed_no-flank",
                     output.format = "phylip",
                     min.alignment.length = 100,
                     min.taxa.alignment = min.taxa.alignment,
                     threads = threads,
                     memory = memory,
                     overwrite = overwrite,
                     mafft.path = mafft.path)

makeIntronAlignments(alignment.directory = "data-analysis/alignments/untrimmed_all-markers",
                     alignment.format = "phylip",
                     output.directory = "data-analysis/alignments/untrimmed_introns",
                     output.format = "phylip",
                     reference.type = "target",
                     reference.path = target.file,
                     target.direction = TRUE,
                     concatenate.intron.flanks = TRUE,
                     threads = threads,
                     memory = memory,
                     overwrite = overwrite,
                     mafft.path = mafft.path)


makeAlignmentSubset(alignment.directory = "data-analysis/alignments/untrimmed_all-markers",
                    alignment.format = "phylip",
                    output.directory = "data-analysis/alignments/untrimmed_UCE",
                    output.format = "phylip",
                    subset.reference = "blast",
                    subset.fasta.file = "subset_uce-consensus.fa",
                    subset.grep.string = NULL,
                    subset.blast.targets = target.file,
                    blast.path = blast.path,
                    threads = threads,
                    memory = memory,
                    overwrite = overwrite)


if (trim.alignments == TRUE){
  #Fix the installs for this
  batchTrimAlignments(alignment.dir = "data-analysis/alignments/untrimmed_all-markers",
                      alignment.format = "phylip",
                      output.dir = "data-analysis/alignments/trimmed_all-markers",
                      output.format = "phylip",
                      overwrite = overwrite,
                      TAPER = FALSE,
                      TAPER.path = taper.path,
                      julia.path = julia.path,
                      TrimAl = run.TrimAl,
                      TrimAl.path = trimAl.path,
                      trim.column = trim.column,
                      convert.ambiguous.sites = convert.ambiguous.sites,
                      alignment.assess = F,
                      trim.external = trim.external,
                      trim.coverage = trim.coverage,
                      min.coverage.percent = min.coverage.percent,
                      min.external.percent = min.external.percent,
                      min.column.gap.percent = min.column.gap.percent,
                      min.alignment.length = min.alignment.length,
                      min.taxa.alignment = min.taxa.alignment,
                      min.coverage.bp = min.coverage.bp,
                      threads = threads,
                      memory = memory)
}#end if


### End workflow
