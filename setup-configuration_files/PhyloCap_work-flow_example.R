#Installs updated package version
devtools::install_github("chutter/PhyloCap", upgrade = "never", force = TRUE)
library(PhyloCap)
library(foreach)

source("/Users/chutter/Dropbox/Research/0_Github/PhyloCap/setup-configuration_files/configuration-file.R")

##################################################################################################
##################################################################################################
#################################################
## Step 1: Preprocess reads
##################

#Checks if everything is installed
pass.fail = setupCheck(anaconda.environment = "/Users/chutter/conda/PhyloCap",
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

if (pass.fail == FALSE){ stop("Some required programs are missing") } else {
  print("all required programs are found, PhyloCap pipeline continuing...")
}

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
                    include.human = include.human,
                    include.univec = include.univec,
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
             min.match.percent = min.match.percent,
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
                    min.alignment.length = min.alignment.length,
                    min.taxa.alignment = min.taxa.alignment.trim,
                    min.alignment.gap.percent = min.alignment.gap.percent,
                    min.coverage.bp = min.coverage.bp,
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


