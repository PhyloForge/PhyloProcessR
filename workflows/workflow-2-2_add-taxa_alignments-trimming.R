#Installs updated package version
devtools::install_github("chutter/PhyloCap", upgrade = "never", dependencies = FALSE)
library(PhyloCap)
library(foreach)

##################################################################################################
##################################################################################################
#################################################
## Step 1: Preprocess reads
##################

#################################################
## Step 3: Align targets and trim
##################

setwd("/Volumes/LaCie/Anolis/data-analysis")
dir.create("new_alignments")

matchTargets(assembly.directory = "/Volumes/LaCie/Anolis/data-analysis/input-samples",
             target.file = "/Volumes/LaCie/Anolis/data-analysis/Hutter_uce5k_loci.fa",
             alignment.contig.name = "anolis_out",
             output.directory = "match-targets",
             min.match.percent = 50,
             min.match.length = 40,
             min.match.coverage = 30,
             threads = 1,
             memory = 1,
             trim.target = FALSE,
             overwrite = FALSE,
             resume = TRUE,
             quiet = TRUE,
             blast.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin",
             bbmap.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin")

#Funciotn for adding taxa to alignments
addTaxaAlignments(alignment.directory = "/Volumes/LaCie/Anolis/data-analysis/original_alignments",
                  alignment.format = "phylip",
                  output.directory = "/Volumes/LaCie/Anolis/data-analysis/new_alignments",
                  output.format = "phylip",
                  sample.markers = "/Volumes/LaCie/Anolis/data-analysis/anolis_out_to-align.fa",
                  threads = 10,
                  memory = 40,
                  overwrite = FALSE,
                  mafft.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin")


#Fix the installs for this
batchTrimAlignments(alignment.dir = "new_alignments",
                    alignment.format = "phylip",
                    output.dir = "trimmed_alignments",
                    output.format = "phylip",
                    overwrite = FALSE,
                    TAPER = FALSE,
                    TAPER.path = NULL,
                    julia.path = NULL,
                    TrimAl = TRUE,
                    TrimAl.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin",
                    trim.column = TRUE,
                    convert.ambiguous.sites = TRUE,
                    alignment.assess = F,
                    trim.external = TRUE,
                    trim.coverage = TRUE,
                    min.coverage.percent = 60,
                    min.external.percent = 60,
                    min.column.gap.percent = 60,
                    min.alignment.length = 200,
                    min.taxa.alignment = 4,
                    min.coverage.bp = 60,
                    threads = 8,
                    memory = 24)




### End workflow
