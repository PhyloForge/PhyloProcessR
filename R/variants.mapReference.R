#' @title variants.mapReference
#'
#' @description Function for running the program spades to assemble short read sequencing data
#'
#' @param read.directory directory of processed reads
#'
#' @param output.directory save name for the output directory
#'
#' @param full.path.spades contigs are added into existing alignment if algorithm is "add"
#'
#' @param mismatch.corrector algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param kmer.values if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param threads number of computation processing threads
#'
#' @param mem amount of system memory to use
#'
#' @param resume TRUE to skip samples already completed
#'
#' @param overwrite TRUE to overwrite a folder of samples with output.dir
#'
#' @param quiet TRUE to supress screen output

#' @return an alignment of provided sequences in DNAStringSet format. Also can save alignment as a file with save.name
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
#'
#'
#' @export

variants.mapReference = function(bam.directory = NULL,
                                 output.directory = "variant-discovery",
                                 reference.file = "all",
                                 subreference.name = NULL,
                                 samtools.path = NULL,
                                 bwa.path = NULL,
                                 picard.path = NULL,
                                 threads = 1,
                                 memory = 1,
                                 resume = TRUE,
                                 overwrite.reference = TRUE,
                                 quiet = TRUE) {

    #Debugging
  #Home directoroies
  # library(doParallel)
  # work.dir = "/Volumes/Armored/Test/variant-calling" #Your main project directory
  # dir.create(work.dir)
  # setwd(work.dir)
  #
  # bam.directory = "/Volumes/Armored/Test/variant-calling/variant-discovery"
  # reference.file = "/Users/chutter/Dropbox/Research/0_Github/FrogCap-Sequence-Capture/Probe_Sets/Ranoidea-V2_Markers.fa"
  # subreference.name = "rag1"
  # output.directory = "variant-discovery"
  # auto.readgroup = T #Keep it T unless it crashes.
  # threads = 4 #number of threads, 8-10 is a recommended amount
  # memory = 8
  # samtools.path = "/Users/chutter/miniconda3/bin"
  # bwa.path = "/usr/local/bin"
  # picard.path = "/Users/chutter/miniconda3/bin"
  # resume = TRUE
  # quiet = FALSE
  # overwrite.reference = TRUE

  #Same adds to bbmap path
  if (is.null(samtools.path) == FALSE){
    b.string = unlist(strsplit(samtools.path, ""))
    if (b.string[length(b.string)] != "/") {
      samtools.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { samtools.path = "" }

  #Same adds to bbmap path
  if (is.null(bwa.path) == FALSE){
    b.string = unlist(strsplit(bwa.path, ""))
    if (b.string[length(b.string)] != "/") {
      bwa.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { bwa.path = "" }

  #Same adds to bbmap path
  if (is.null(picard.path) == FALSE){
    b.string = unlist(strsplit(picard.path, ""))
    if (b.string[length(b.string)] != "/") {
      picard.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { picard.path = "" }

  #Quick checks
  if (is.null(bam.directory) == TRUE){ stop("Please provide the bam directory.") }
  if (file.exists(bam.directory) == F){ stop("BAM folder not found.") }

  #Creates output directory
  if (dir.exists("logs") == F){ dir.create("logs") }

  #Read in sample data
  bam.files = list.files(bam.directory, recursive = T, full.names = T)
  bam.files = bam.files[grep("all_reads.bam$", bam.files)]
  sample.names = list.dirs(bam.directory, recursive = F, full.names = F)

  #Resumes file download
  if (resume == TRUE && overwrite == FALSE){
    done.files = list.files(output.directory, full.names = T, recursive = T)
    done.files = done.files[grep("final-mapped-all.bam", done.files)]
    done.names = gsub("/Lane_.*", "", done.files)
    done.names = gsub(".*\\/", "", done.names)
    sample.names = sample.names[!sample.names %in% done.names]
  }

  if (length(sample.names) == 0){ return("no samples remain to analyze.") }
  if (subreference.name == "all") { subreference.name = NULL }

  ############################################################################################
  ########### Step 0 #########################################################################
  ##### Set up reference
  ############################################################################################

  #Create combined and indexed reference
  reference.location = paste0("ref-index/reference.fa")
  if (dir.exists("ref-index") == TRUE && overwrite.reference == TRUE){ system(paste0("rm -rf ref-index")) }
  if (dir.exists("ref-index") == FALSE){
    #Makes new directory
    dir.create("ref-index")

    if (is.null(subreference.name) != TRUE){

      reference.seq = Biostrings::readDNAStringSet(filepath = reference.file)
      sub.ref = reference.seq[grep(subreference.name, names(reference.seq))]
      #Saves final set
      final.loci = as.list(as.character(sub.ref))
      writeFasta(sequences = final.loci, names = names(final.loci),
                 paste0("ref-index/reference.fa"),
                 nbchar = 1000000, as.string = T)
    } else {
      #If no subreference copy actual reference
      system(paste0("cp ", reference.file, " ref-index/reference.fa"))
    }

    #Indexes the reference
    system(paste0(bwa.path, "bwa index -p ref-index/reference ", reference.location),
           ignore.stderr = quiet, ignore.stdout = quiet)

    system(paste0(samtools.path, "samtools faidx ", reference.location))
    system(paste0(picard.path, "picard -Xmx", memory, "G",
                  " CreateSequenceDictionary REFERENCE=", reference.location,
                  " OUTPUT=ref-index/reference.dict",
                  " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"))
  }#end if

  ############################################################################################
  ########### Step 1 #########################################################################
  ##### Start up loop for each sample
  ############################################################################################

  #Runs through each sample
  for (i in 1:length(sample.names)) {
    #################################################
    ### Part A: prepare for loading and checks
    #################################################
    sample.dir = paste0(output.directory, "/", sample.names[i])
    if (file.exists(sample.dir) == FALSE) { dir.create(sample.dir) }

    #Gets the reads for the sample
    sample.bams = bam.files[grep(pattern = paste0(sample.names[i], "/"), x = bam.files)]
    #Checks the Sample column in case already renamed
    if (length(sample.bams) == 0){ sample.bams = sample.bams[grep(pattern = sample.names[i], x = sample.bams)] }

    sample.bams = unique(gsub("_1.f.*|_2.f.*|_3.f.*|-1.f.*|-2.f.*|-3.f.*|_R1_.*|_R2_.*|_R3_.*|_READ1_.*|_READ2_.*|_READ3_.*|_R1.f.*|_R2.f.*|_R3.f.*|-R1.f.*|-R2.f.*|-R3.f.*|_READ1.f.*|_READ2.f.*|_READ3.f.*|-READ1.f.*|-READ2.f.*|-READ3.f.*|_singleton.*|-singleton.*|READ-singleton.*|READ_singleton.*|_READ-singleton.*|-READ_singleton.*|-READ-singleton.*|_READ_singleton.*", "", sample.bams))

    #Returns an error if reads are not found
    if (length(sample.bams) == 0 ){
      stop(sample.names[i], " does not have any reads present for files ")
    } #end if statement

    #CReates new directory
    report.path = paste0("logs/", sample.names[i])
    if (file.exists(report.path) == FALSE) { dir.create(report.path) }

    for (j in 1:length(sample.bams)){

      lane.bams = sample.bams[grep(pattern = paste0(sample.bams[j], "_"), x = sample.bams)]

      #Checks the Sample column in case already renamed
      if (length(lane.bams) == 0){ lane.bams = sample.bams[grep(pattern = sample.bams[j], x = sample.bams)] }
      #Returns an error if reads are not found
      if (length(lane.bams) == 0 ){
        stop(sample.bams[j], " does not have any reads present for files ")
      } #end if statement

      #Gets lane names
      lane.name = paste0("Lane_", j)
      lane.dir = paste0(sample.dir, "/", lane.name)
      dir.create(lane.dir)

      #Run BWA Mem
      system("set -o pipefail")

      #Piped verison
      tmp.dir = paste0(lane.dir, "/tmp")
      system(paste0(picard.path, "picard -Xmx", memory, "G",
                    " SamToFastq I=", lane.dir, "/all_reads.bam FASTQ=/dev/stdout TMP_DIR=", tmp.dir,
                    " CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true",
                    " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true | ",
                    bwa.path, "bwa mem -M -p -t ", threads, " ref-index/reference /dev/stdin | ",
                    picard.path, "picard -Xmx", memory, "G",
                    " MergeBamAlignment ALIGNED_BAM=/dev/stdin UNMAPPED_BAM=", lane.dir, "/all_reads.bam",
                    " OUTPUT=", lane.dir, "/cleaned_final.bam",
                    " R=", reference.location, " CREATE_INDEX=true ADD_MATE_CIGAR=true",
                    " CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true",
                    " MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS",
                    " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true",
                    " TMP_DIR=", tmp.dir))

      system(paste0("rm -r ", tmp.dir))
      #system(paste0("samtools view -H cleaned_final.bam | grep '@RG'"))

      ############################################################################################
      ########### Step 4 #########################################################################
      ##### Sort Sam, Mark duplicates, Sort Sam again, finalize mapped set ofreads
      ############################################################################################

      if (file.exists(paste0(lane.dir, "/cleaned_final.bam")) != TRUE){ stop("Stopped. Something failed after mapping before sorting.")}


      #Sort by coordinate for input into MarkDuplicates
      system(paste0(picard.path, "picard -Xmx", memory, "G",
                    " SortSam INPUT=", lane.dir, "/cleaned_final.bam",
                    " OUTPUT=", lane.dir, "/cleaned_final_sort.bam",
                    " CREATE_INDEX=true SORT_ORDER=coordinate",
                    " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"))

      #Marks duplicate reads
      system(paste0(picard.path, "picard -Xmx", memory, "G",
                    " MarkDuplicates INPUT=", lane.dir, "/cleaned_final_sort.bam",
                    " OUTPUT=", lane.dir, "/cleaned_final_md.bam",
                    " CREATE_INDEX=true METRICS_FILE=", report.path, "/duplicate_metrics.txt",
                    " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"))

      #Sorts and stuff
      system(paste0(picard.path, "picard -Xmx", memory, "G",
                    " SortSam INPUT=", lane.dir, "/cleaned_final_md.bam",
                    " OUTPUT=/dev/stdout SORT_ORDER=coordinate | ",
                    picard.path, "picard -Xmx", memory, "G",
                    " SetNmAndUqTags INPUT=/dev/stdin OUTPUT=", lane.dir, "/final-mapped-all.bam",
                    " CREATE_INDEX=true R=", reference.location,
                    " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"))

      if (file.exists(paste0(lane.dir, "/final-mapped-all.bam")) != TRUE){ stop("Stopped. Something failed after mapping during sorting.")}

      #Intermediate files are deleted
      system(paste0("rm ", lane.dir, "/cleaned_final_sort*"))
      system(paste0("rm ", lane.dir, "/cleaned_final*"))

      print(paste0(sample.names[i], " ", lane.name, " completed read mapping to reference!"))

    }#end j loop

    print(paste0(sample.names[i], " completed read mapping to reference!"))

  }#end i loop

}#end function

