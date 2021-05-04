#Installs updated package version
devtools::install_github("chutter/PhyloCap", upgrade = "never", force = TRUE)
library(PhyloCap)
library(foreach)

#################################################
## Configuration file
##################

#Directories
work.dir = "/Users/chutter/Dropbox/Research/0_Github/Test-dataset"
read.dir = "raw-reads"
file.rename = "/Users/chutter/Dropbox/Research/0_Github/Test-dataset/file_rename.csv"
target.file = "/Users/chutter/Dropbox/Research/0_Github/PhyloCap/setup-configuration_files/Ranoidea_All-Markers_Apr21-2019.fa"
dataset.name = "Test"

#Basic settings
threads = 4
memory = 8
overwrite = FALSE
resume = TRUE
quiet = TRUE

#Pre-processing settings
decontamination = TRUE
download.contaminant.genomes = TRUE
contaminant.genome.list = "/Users/chutter/Dropbox/Research/0_Github/PhyloCap/setup-configuration_files/decontamination_database.csv"
decontamination.path = NULL
include.human = TRUE  #Include human genome, unless human is study organism or UCEs in mammals are used
include.mouse = TRUE  #Include human genome, unless human is study organism or UCEs in mammals are used
decontamination.match = 0.99 #
spades.kmer.values = c(21,33,55,77,99,127)
spades.mismatch.corrector = TRUE
save.corrected.reads = FALSE

#Target matching and alignment settings
trim.to.targets = FALSE #whether to trim to the targets or keep the entire contig
min.percent.id = 0.5
min.match.length = 40
min.match.coverage = 0.5
min.taxa.alignment = 4 #minimum taxa for an alignment
min.taxa.tree = 4 #minimum taxa for an alignment

#Alignement trimming parameters
run.TAPER = TRUE
run.TrimAl = TRUE
trim.column = TRUE
convert.ambiguous.sites = TRUE
alignment.assess = TRUE
trim.external = TRUE
trim.coverage = TRUE
min.coverage.percent = 35
min.external.percent = 50
min.column.gap.percent = 50
min.align.length = 100
min.taxa.count = 5
min.gap.percent = 50
min.sample.bp = 60
cleanup.genetrees = TRUE #Only saves the ML tree; FALSE saves all IQTree files for each gene tree

#Program paths
fastp.path = "/Users/chutter/conda/PhyloCap/bin"
samtools.path = "/Users/chutter/conda/PhyloCap/bin"
bwa.path = "/Users/chutter/conda/PhyloCap/bin"
spades.path = "/Users/chutter/conda/PhyloCap/bin"
bbmap.path = "/Users/chutter/conda/PhyloCap/bin"
blast.path = "/Users/chutter/conda/PhyloCap/bin"
mafft.path = "/Users/chutter/conda/PhyloCap/bin"
iqtree.path = "/Users/chutter/conda/PhyloCap/bin"
trimAl.path = "/Users/chutter/conda/PhyloCap/bin"
taper.path = "/Users/chutter/conda/PhyloCap/bin"
julia.path = "/Users/chutter/conda/PhyloCap/bin"


##################################################################################################
##################################################################################################
#################################################
## Step 1: Preprocess reads
##################

setwd(work.dir)
dir.create("processed-reads")

organizeReads(read.directory = read.dir,
              output.dir = "processed-reads/organized-reads",
              rename.file = file.rename,
              overwrite = overwrite)

removeAdaptors(input.reads = "processed-reads/organized-reads",
               output.directory = "processed-reads/adaptor-removed-reads",
               fastp.path = fastp.path,
               threads = threads,
               mem = memory,
               resume = resume,
               overwrite = overwrite,
               quiet = quiet)

#Creates the database by downloading
createContaminantDB(decontamination.list = contaminant.genome.list,
                    output.directory = "contaminant-references",
                    include.human = TRUE,
                    include.univec = TRUE,
                    overwrite = overwrite)

## remove external contamination
removeContamination(input.reads = "processed-reads/adaptor-removed-reads",
                    output.directory = "processed-reads/decontaminated-reads",
                    decontamination.path = "contaminant-references",
                    map.match = decontamination.match,
                    samtools.path = samtools.path,
                    bwa.path = bwa.path,
                    threads = threads,
                    mem = memory,
                    resume = resume,
                    overwrite = overwrite,
                    overwrite.reference = overwrite,
                    quiet = quiet)

#merge paired end reads
mergePairedEndReads(input.reads = "processed-reads/decontaminated-reads",
                    output.directory = "processed-reads/pe-merged-reads",
                    fastp.path = fastp.path,
                    threads = threads,
                    mem = memory,
                    resume = resume,
                    overwrite = overwrite,
                    quiet = quiet)

dir.create("data-analysis")

#Assembles merged paired end reads with spades
assembleSpades(input.reads = "processed-reads/pe-merged-reads",
               output.directory = "processed-reads/spades-assembly",
               assembly.directory = "data-analysis/draft-assemblies",
               mismatch.corrector = spades.mismatch.corrector,
               kmer.values = spades.kmer.values,
               threads = threads,
               memory = memory,
               overwrite = overwrite,
               resume = resume,
               save.corrected.reads = save.corrected.reads,
               quiet = quiet,
               spades.path = spades.path)

#################################################
## Step 2: Match targets and annotate contigs
##################

#match targets
matchTargets(assembly.directory = "data-analysis/draft-assemblies",
             target.file = target.file,
             alignment.contig.name = paste0("data-analysis/", dataset.name),
             output.directory = "data-analysis/match-targets",
             min.percent.id = min.percent.id,
             min.match.length = min.match.length,
             min.match.coverage = min.match.coverage,
             threads = threads,
             memory = memory,
             trim.target = trim.to.targets,
             overwrite = overwrite,
             resume = resume,
             quiet = quiet,
             blast.path = blast.path,
             bbmap.path = bbmap.path)

#################################################
## Step 3: Align targets and trim
##################

#align targets
alignTargets(targets.to.align = paste0("data-analysis/", dataset.name, "_to-align.fa"),
             output.directory = "data-analysis/alignments",
             min.taxa = min.taxa.alignment,
             subset.start = 0,
             subset.end = 1,
             threads = threads,
             memory = memory,
             overwrite = overwrite,
             resume = resume,
             quiet = quiet,
             mafft.path = mafft.path)

#Fix the installs for this
batchTrimAlignments(alignment.dir = "data-analysis/alignments",
                    alignment.format = "phylip",
                    output.dir = "data-analysis/alignments-trimmed",
                    output.format = "phylip",
                    overwrite = overwrite,
                    resume = resume,
                    TAPER = run.TAPER,
                    TAPER.path = taper.path,
                    julia.path = julia.path,
                    TrimAl = run.TrimAl,
                    TrimAl.path = trimAl.path,
                    trim.column = trim.column,
                    convert.ambiguous.sites = convert.ambiguous.sites,
                    alignment.assess = alignment.assess,
                    trim.external = trim.external,
                    trim.coverage = trim.coverage,
                    min.coverage.percent = min.coverage.percent,
                    min.external.percent = min.external.percent,
                    min.column.gap.percent = min.column.gap.percent,
                    min.align.length = min.align.length,
                    min.taxa.count = min.taxa.count,
                    min.gap.percent = min.gap.percent,
                    min.sample.bp = min.sample.bp,
                    threads = threads,
                    memory = memory)

#################################################
## Step 4: Tree stuff
##################

estimateGeneTrees(alignment.directory = "data-analysis/alignments-trimmed",
                  output.directory = "data-analysis/gene-trees",
                  min.taxa = min.taxa.tree,
                  subset.start = 0,
                  subset.end = 1,
                  threads = threads,
                  memory = memory,
                  overwrite = overwrite,
                  resume = resume,
                  quiet = quiet,
                  cleanup.files = cleanup.genetrees,
                  iqtree.path = iqtree.path)


