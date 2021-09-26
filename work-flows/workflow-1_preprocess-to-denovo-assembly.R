#Installs updated package version
devtools::install_github("chutter/PhyloCap", upgrade = "never", force = TRUE)
library(PhyloCap)
library(foreach)

source("workflow-1_configuration-file.R")
setwd(working.directory)

##################################################################################################
##################################################################################################
#################################################
## Step 1: Preprocess reads
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

if (pass.fail == FALSE){ stop("Some required programs are missing") } else {
  print("all required programs are found, PhyloCap pipeline continuing...")
}

#Begins by creating processed read directory
dir.create("processed-reads")

if (dropbox.download == TRUE){
  #Authorizes token
  rdrop2::drop_auth(rdstoken = dropbox.token)
  #Run download function
  dropboxDownload(sample.spreadsheet = sample.file,
                  dropbox.directory = dropbox.directory,
                  out.directory = paste0(processed.reads, "/raw-reads"),
                  overwrite = TRUE,
                  resume = FALSE)

  read.dir = paste0(processed.reads, "/raw-reads")
  organize.reads = FALSE
}#end if

#Organizes reads if scattered elsewhere i.e. creates a sub-dataset
if (organize.reads == TRUE) {
  organizeReads(read.directory = read.dir,
                output.dir = paste0(processed.reads, "/organized-reads"),
                rename.file = sample.file,
                overwrite = overwrite)
  input.reads = paste0(processed.reads, "/organized-reads")
} else {input.reads = read.dir }

if (summary.fastq == TRUE){
  #This function creates a summary of the fastq files per sample for number of reads
  summary.fastqStats(read.directory = input.reads,
                     output.name = "fastq-stats",
                     read.length = 150,
                     threads = threads,
                     mem = memory,
                     overwrite = overwrite)
}#end summary.fastq if

if (remove.adaptors == TRUE) {
  removeAdaptors(input.reads = input.reads,
                 output.directory = paste0(processed.reads, "/adaptor-removed-reads"),
                 fastp.path = fastp.path,
                 threads = threads,
                 mem = memory,
                 resume = resume,
                 overwrite = overwrite,
                 quiet = quiet)
  input.reads = paste0(processed.reads, "/adaptor-removed-reads")
}

#Runs decontamination of reads
if (decontamination == TRUE){
  #Creates the database by downloading
  createContaminantDB(decontamination.list = contaminant.genome.list,
                      output.directory = "contaminant-references",
                      include.human = include.human,
                      include.univec = include.univec,
                      overwrite = overwrite)

  ## remove external contamination
  removeContamination(input.reads = input.reads,
                      output.directory = paste0(processed.reads, "/decontaminated-reads"),
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
  input.reads = paste0(processed.reads, "/decontaminated-reads")
}


if (merge.pe.reads == TRUE){
  #merge paired end reads
  mergePairedEndReads(input.reads = input.reads,
                      output.directory =  paste0(processed.reads, "/pe-merged-reads"),
                      fastp.path = fastp.path,
                      threads = threads,
                      mem = memory,
                      resume = resume,
                      overwrite = overwrite,
                      quiet = quiet)
  input.reads = paste0(processed.reads, "/pe-merged-reads")
} #end decontamination

dir.create("data-analysis")

if (denovo.assembly == TRUE){
  #Assembles merged paired end reads with spades
  assembleSpades(input.reads = input.reads,
                 output.directory = paste0(processed.reads, "/spades-assembly"),
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
}#end if

#End script
