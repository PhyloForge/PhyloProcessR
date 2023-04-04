#' @title variants.haplotypeCallerGATK4
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

variants.haplotypeCallerGATK4 = function(bam.directory = NULL,
                                         output.directory = "haplotype-caller",
                                         samtools.path = NULL,
                                         gatk4.path = NULL,
                                         threads = 1,
                                         memory = 1,
                                         overwrite = TRUE,
                                         quiet = TRUE) {

  #Debugging
  #Home directoroies
  # Debugging
  # Home directoroies
  library(PhyloCap)
  library(doParallel)
  setwd("/Volumes/LaCie/Mantellidae")
  output.directory <- "variant-discovery/haplotype-caller"
  bam.directory <- "/Volumes/LaCie/Mantellidae/variant-discovery/sample-mapping"

  gatk4.path <- "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"
  samtools.path <- "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"
  bwa.path <- "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"

  threads <- 4
  memory <- 8
  quiet <- FALSE
  overwrite <- TRUE


  # Same adds to bbmap path
  require(doParallel)
  
  # Same adds to bbmap path
  if (is.null(samtools.path) == FALSE) {
    b.string <- unlist(strsplit(samtools.path, ""))
    if (b.string[length(b.string)] != "/") {
      samtools.path <- paste0(append(b.string, "/"), collapse = "")
    } # end if
  } else {
    samtools.path <- ""
  }

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
  if (is.null(bam.directory) == TRUE){ stop("Please provide the bam directory.") }
  if (file.exists(bam.directory) == F){ stop("BAM folder not found.") }

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

  #Read in sample data
  bam.files = list.files(bam.directory, recursive = T, full.names = T)
  bam.files = bam.files[grep("final-mapped-all.bam$", bam.files)]
  sample.names <- list.dirs(bam.directory, recursive = F, full.names = F)

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
  cl = makeCluster(threads)
  registerDoParallel(cl)
  mem.cl = floor(memory/threads)

  #Loops through each locus and does operations on them
  foreach(i=1:length(sample.names), .packages = c("foreach", "ShortRead")) %dopar% {
    #Runs through each sample
    #for (i in 1:length(sample.names)) {
    #################################################
    ### Part A: prepare for loading and checks
    #################################################
    sample.dir = paste0(output.directory, "/", sample.names[i])
    if (file.exists(sample.dir) == FALSE) { dir.create(sample.dir) }

    #Gets the reads for the sample
    sample.bams = bam.files[grep(pattern = paste0(sample.names[i], "/"), x = bam.files)]
    #Checks the Sample column in case already renamed
    if (length(sample.bams) == 0){ sample.bams = sample.bams[grep(pattern = sample.names[i], x = sample.bams)] }

    #Returns an error if reads are not found
    if (length(sample.bams) == 0 ){
      stop(sample.names[i], " does not have any reads present for files ")
    } #end if statement

    #CReates new directory
    report.path = paste0("logs/", sample.names[i])
    if (file.exists(report.path) == FALSE) { dir.create(report.path) }

    # Sets up merging of bams from different lanes
    input.string = paste0("-I ", sample.bams, collapse = " ")
    reference.path = paste0(bam.directory, "/", sample.names[i], "/index/reference.fa")

    if (length(sample.bams) != 1) {

      merge.dir = paste0(bam.directory, "/", sample.names[i], "/Lane_Merge")
      dir.create(merge.dir)

      # Next combine .bam files together!
      system(paste0(
        gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
        " MergeSamFiles",
        " ", input.string, " -O ", merge.dir, "/final-mapped-merge.bam",
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
      ))

      # Sort by coordinate for input into MarkDuplicates
      system(paste0(
        gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
        " SortSam",
        " -INPUT ", merge.dir, "/final-mapped-merge.bam",
        " -OUTPUT ", merge.dir, "/final-mapped-sort.bam",
        " -CREATE_INDEX true -SORT_ORDER coordinate",
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
      ))

      # Marks duplicate reads
      system(paste0(
        gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
        " MarkDuplicates",
        " -INPUT ", merge.dir, "/final-mapped-sort.bam",
        " -OUTPUT ", merge.dir, "/final-mapped-dup.bam",
        " -CREATE_INDEX true -METRICS_FILE logs/", sample.names[i], "/duplicate_metrics.txt",
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
      ))

      # Sorts and stuff
      system(paste0(
        gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
        " SortSam",
        " -INPUT ", merge.dir, "/final-mapped-dup.bam",
        " -OUTPUT /dev/stdout -SORT_ORDER coordinate",
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true | ",
        gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
        " SetNmAndUqTags",
        " -INPUT /dev/stdin -OUTPUT ", merge.dir, "/final-mapped-all.bam",
        " -CREATE_INDEX true -R ", reference.path, "/reference.fa",
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
      ))

      input.bam = paste0(merge.dir, "/final-mapped-all.bam")

      # Delete old files to make more space
      system(paste0("rm ", merge.dir, "/final-mapped-dup.bam"))
      system(paste0("rm ", merge.dir, "/final-mapped-sort.bam"))
      system(paste0("rm ", merge.dir, "/final-mapped-merge.bam"))
    } else {
      # lane 1 if thats all there is
      input.bam = paste0(bam.directory, "/", sample.names[i], "/Lane_1/final-mapped-all.bam")
    } # end else

    #Starts to finally look for Haplotypes! *here
    system(paste0(
      gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
      " HaplotypeCaller",
      " -R ", reference.path, " -O ", sample.dir, "/gatk4-haplotype-caller.g.vcf.gz",
      " -I ", input.bam,
      " -ERC GVCF --max-alternate-alleles 3",
      " -bamout ", sample.dir, "/gatk4-haplotype-caller.bam"
    ))

    print(paste0(sample.names[i], " completed GATK4 haplotype caller!"))

  }# end i loop

  stopCluster(cl)


}#end function

# END SCRIPT
