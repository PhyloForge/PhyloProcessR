#' @title haplotypeCaller
#'
#' @description Runs GATK4 HaplotypeCaller in GVCF mode (-ERC GVCF) on
#'   per-sample BAM files. When a sample has multiple lanes, the lane BAMs are
#'   first merged, sorted, and deduplicated with GATK MergeSamFiles,
#'   MarkDuplicates, and SetNmAndUqTags before calling haplotypes. The reference
#'   can be a per-sample assembly or a shared consensus reference. Samples are
#'   processed in parallel.
#'
#' @param mapping.directory path to the directory of per-sample BAM files and
#'   reference indices (output of mapReferenceSample() or
#'   mapReferenceConsensus()).
#'
#' @param output.directory path to the directory where per-sample GVCF files
#'   will be saved.
#'
#' @param reference.type character; "sample" to use each sample's own reference
#'   FASTA (located at mapping.directory/sample/index/reference.fa), or
#'   "consensus" to use a shared reference at index/reference.fa in the working
#'   directory.
#'
#' @param gatk4.path system path to the directory containing the gatk
#'   executable; NULL searches the system PATH.
#'
#' @param temp.directory path to a temporary directory for GATK JVM temp files;
#'   NULL uses the current working directory.
#'
#' @param ploidy integer ploidy to pass to HaplotypeCaller (-ploidy).
#'
#' @param threads number of parallel samples to process simultaneously.
#'
#' @param memory total RAM in GB to allocate as the JVM heap (-Xmx).
#'
#' @param overwrite logical; if TRUE the output directory is deleted and
#'   recreated before processing.
#'
#' @param quiet logical; currently unused.
#'
#' @return invisibly; writes per-sample GVCF (.g.vcf.gz) files and realigned
#'   BAM files to output.directory.
#'
#' @export

haplotypeCaller = function(mapping.directory = NULL,
                          output.directory = "haplotype-caller",
                          reference.type = c("sample", "consensus"),
                          gatk4.path = NULL,
                          temp.directory = NULL,
                          ploidy = 2,
                          threads = 1,
                          memory = 1,
                          overwrite = FALSE,
                          quiet = TRUE) {

  #Debugging
  #Home directoroies
  # Debugging
  # Home directoroies
  # library(PhyloCap)
  # library(doParallel)
  # setwd("/Volumes/LaCie/Mantellidae")
  # output.directory <- "variant-discovery/haplotype-caller"
  # mapping.directory <- "/Volumes/LaCie/Mantellidae/variant-discovery/sample-mapping"

  # gatk4.path <- "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"
  # samtools.path <- "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"
  # bwa.path <- "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"

  # threads <- 4
  # memory <- 8
  # quiet <- FALSE
  # overwrite <- TRUE

  # Same adds to bbmap path
  require(foreach)

  # Same adds to bbmap path
  if (is.null(gatk4.path) == FALSE) {
    b.string <- unlist(strsplit(gatk4.path, ""))
    if (b.string[length(b.string)] != "/") {
      gatk4.path <- paste0(append(b.string, "/"), collapse = "")
    } # end if
  } else {
    gatk4.path <- ""
  }

  #Quick checks
  if (is.null(mapping.directory) == TRUE){ stop("Please provide the bam directory.") }
  if (file.exists(mapping.directory) == F){ stop("BAM folder not found.") }

  # Creates output directory
  if (dir.exists("logs") == F) {
    dir.create("logs")
  }

  # Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.directory) == F) {
    dir.create(output.directory)
  } else {
    if (overwrite == TRUE) {
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } # end else

  if (is.null(temp.directory) == TRUE){
    temp.directory = getwd()
  }


  #Read in sample data
  bam.files = list.files(mapping.directory, recursive = T, full.names = T)
  bam.files = bam.files[grep("final-mapped-all.bam$", bam.files)]
  sample.names <- list.dirs(mapping.directory, recursive = F, full.names = F)

  # Resumes file download
  if (overwrite == FALSE) {
    done.files <- list.files(output.directory, full.names = T, recursive = T)
    done.files <- done.files[grep("gatk4-haplotype-caller.g.vcf.gz$", done.files)]
    done.names <- gsub("/gatk4-haplotype-caller.g.vcf.gz$", "", done.files)
    done.names <- gsub(".*\\/", "", done.names)
    sample.names <- sample.names[!sample.names %in% done.names]
  }

  if (length(sample.names) == 0){ return("no samples remain to analyze.") }

  ############################################################################################
  ########### Step 1 #########################################################################
  ##### Start up loop for each sample
  ############################################################################################
  #Sets up multiprocessing
  cl <- parallel::makeCluster(threads, outfile = "")
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  mem.cl <- floor(memory / threads)

  #Loops through each locus and does operations on them
  foreach(i=1:length(sample.names), .packages = c("foreach")) %dopar% {
    #Runs through each sample
    #for (i in 1:length(sample.names)) {
    #################################################
    ### Part A: prepare for loading and checks
    #################################################
    sample.dir = paste0(output.directory, "/", sample.names[i])
    if (file.exists(sample.dir) == FALSE) { dir.create(sample.dir) }

    #Gets the reads for the sample
    sample.bams = bam.files[grep(pattern = paste0(sample.names[i], "/"), x = bam.files)]
    if (length(sample.bams) == 0){ sample.bams = bam.files[grep(pattern = sample.names[i], x = bam.files)] }

    #Returns an error if reads are not found
    if (length(sample.bams) == 0 ){
      stop(sample.names[i], " does not have any reads present for files ")
    } #end if statement

    #CReates new directory
    report.path = paste0("logs/", sample.names[i])
    if (file.exists(report.path) == FALSE) { dir.create(report.path) }

    # Sets up merging of bams from different lanes
    input.string = paste0("-I ", sample.bams, collapse = " ")
    if (reference.type == "sample") {
      reference.path = paste0(mapping.directory, "/", sample.names[i], "/index/reference.fa")
    }

    if (reference.type == "consensus") {
      reference.path = paste0("index/reference.fa")
    }

    if (length(sample.bams) != 1) {

      merge.dir = paste0(mapping.directory, "/", sample.names[i], "/Lane_Merge")
      dir.create(merge.dir)

      # Next combine .bam files together!
      system(paste0(
        gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", mem.cl, "G\"",
        " MergeSamFiles",
        " ", input.string, " -O ", merge.dir, "/final-mapped-merge.bam",
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
      ))

      # Sort by coordinate for input into MarkDuplicates
      system(paste0(
        gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", mem.cl, "G\"",
        " SortSam",
        " -INPUT ", merge.dir, "/final-mapped-merge.bam",
        " -OUTPUT ", merge.dir, "/final-mapped-sort.bam",
        " -CREATE_INDEX true -SORT_ORDER coordinate",
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
      ))

      # Marks duplicate reads
      system(paste0(
        gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", mem.cl, "G\"",
        " MarkDuplicates",
        " -INPUT ", merge.dir, "/final-mapped-sort.bam",
        " -OUTPUT ", merge.dir, "/final-mapped-dup.bam",
        " -CREATE_INDEX true -METRICS_FILE logs/", sample.names[i], "/duplicate_metrics.txt",
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
      ))

      # Sorts and stuff
      system(paste0(
        gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", mem.cl, "G\"",
        " SortSam",
        " -INPUT ", merge.dir, "/final-mapped-dup.bam",
        " -OUTPUT /dev/stdout -SORT_ORDER coordinate",
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true | ",
        gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", mem.cl, "G\"",
        " SetNmAndUqTags",
        " -INPUT /dev/stdin -OUTPUT ", merge.dir, "/final-mapped-all.bam",
        " -CREATE_INDEX true -R ", reference.path,
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
      ))

      input.bam = paste0(merge.dir, "/final-mapped-all.bam")

      # Delete old files to make more space
      system(paste0("rm ", merge.dir, "/final-mapped-dup.bam"))
      system(paste0("rm ", merge.dir, "/final-mapped-sort.bam"))
      system(paste0("rm ", merge.dir, "/final-mapped-merge.bam"))
    } else {
      # lane 1 if thats all there is
      input.bam = paste0(mapping.directory, "/", sample.names[i], "/Lane_1/final-mapped-all.bam")
    } # end else

    #Starts to finally look for Haplotypes! *here
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", mem.cl, "G\"",
      " HaplotypeCaller",
      " -R ", reference.path, " -O ", sample.dir, "/gatk4-haplotype-caller.g.vcf.gz",
      " -I ", input.bam,
      " -ERC GVCF",
      " -ploidy ", ploidy,
      " -bamout ", sample.dir, "/gatk4-haplotype-caller.bam"
    ))

    print(paste0(sample.names[i], " completed GATK4 haplotype caller!"))

  }# end i loop

  parallel::stopCluster(cl)


}#end function

# END SCRIPT
