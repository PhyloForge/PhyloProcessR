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

setwd(work.dir)
dir.create("processed-reads")

if (organize.reads == TRUE) {
  organizeReads(read.directory = read.dir,
                output.dir = paste0(processed.reads, "/organized-reads"),
                rename.file = file.rename,
                overwrite = overwrite)
  input.reads = paste0(processed.reads, "/organized-reads")
} else {input.reads = read.dir }

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

#################################################
#################################################
#################################################
#################################################
## Step 2: Match targets and annotate contigs
##################

work.dir<-"/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/Microhylidae_SeqCap/New_Work_2021" #Your main project directory
setwd(work.dir)
dir.create("data-analysis")
assembly.directory<-"/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/Microhylidae_SeqCap/New_Work_2021/Assembled_Contigs"
target.file<-"/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/Microhylidae_SeqCap/New_Work_2021/Master_Ranoidea_All-Markers_Apr21-2019.fa"
subset.fasta.file = "/Users/chutter/Dropbox/Research/0_Github/FrogCap_Pipeline/Source_Files/Hutter_uce5k_loci.fa"

if (match.targets == TRUE){
  #match targets
  matchTargets(assembly.directory = assembly.directory,
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
}#end if

#################################################
## Step 3: Align targets and trim
##################

if (align.matching.targets == TRUE){
  #align targets
  alignTargets(targets.to.align = paste0("data-analysis/", dataset.name, "_to-align.fa"),
               output.directory = "alignments/untrimmed_all-markers",
               min.taxa = min.taxa.alignment,
               subset.start = 0,
               subset.end = 1,
               threads = threads,
               memory = memory,
               overwrite = overwrite,
               resume = resume,
               quiet = quiet,
               mafft.path = mafft.path)
}#end if

#### Functions for separating into data types
trimAlignmentTargets(alignment.directory = "alignments/untrimmed_all-markers",
                     alignment.format = "phylip",
                     target.file = target.file,
                     target.direction = TRUE,
                     output.directory = "alignments/untrimmed_no-flank",
                     output.format = "phylip",
                     min.alignment.length = 100,
                     min.taxa.alignment = min.taxa.alignment,
                     threads = threads,
                     memory = memory,
                     overwrite = overwrite,
                     resume = resume,
                     mafft.path = mafft.path)

makeIntronAlignments(alignment.directory = "alignments/untrimmed_all-markers",
                     alignment.format = "phylip",
                     output.directory = "alignments/untrimmed_introns",
                     output.format = "phylip",
                     reference.type = "target",
                     reference.path = target.file,
                     target.direction = TRUE,
                     concatenate.intron.flanks = TRUE,
                     threads = threads,
                     memory = memory,
                     overwrite = overwrite,
                     resume = resume,
                     mafft.path = mafft.path)


makeAlignmentSubset(alignment.directory = "alignments/untrimmed_all-markers",
                    alignment.format = "phylip",
                    output.directory = "alignments/untrimmed_UCE",
                    output.format = "phylip",
                    subset.reference = "blast",
                    subset.fasta.file = "/Users/chutter/Dropbox/Research/0_Github/FrogCap_Pipeline/Source_Files/Hutter_uce5k_loci.fa",
                    subset.grep.string = NULL,
                    subset.blast.targets = "/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/Microhylidae_SeqCap/New_Work_2021/Master_Ranoidea_All-Markers_Apr21-2019.fa",
                    blast.path = blast.path,
                    threads = threads,
                    memory = memory,
                    overwrite = overwrite)


if (trim.alignments == TRUE){
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
}#end if

#################################################
## Step 4: Tree stuff
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
                    resume = resume,
                    quiet = quiet,
                    cleanup.files = cleanup.genetrees,
                    iqtree.path = iqtree.path)

}#end
