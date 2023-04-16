#' @title variants.prepareBAM
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
#' your.tree <- ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data <- astralPlane(
#'   astral.tree = your.tree,
#'   outgroups = c("species_one", "species_two"),
#'   tip.length = 1
#' )
#'
#' @export

variants.prepareBAM <- function(read.directory = NULL,
                                output.directory = "variant-discovery/sample-mapping",
                                samtools.path = NULL,
                                bwa.path = NULL,
                                gatk4.path = NULL,
                                threads = 1,
                                memory = 1,
                                overwrite = FALSE,
                                auto.readgroup = TRUE,
                                quiet = TRUE) {
 
  # library(PhyloCap)
  # library(doParallel)
  # setwd("/Volumes/LaCie/Mantellidae")
  # assembly.directory <- "/Volumes/LaCie/Mantellidae/expanded-assemblies"
  # output.directory <- "variant-discovery/sample-mapping"
  # read.directory <- "/Volumes/LaCie/Mantellidae/reads"

  # gatk4.path <- "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"
  # samtools.path <- "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"
  # bwa.path <- "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"

  # auto.readgroup <- T
  # threads <- 4
  # memory <- 8
  # quiet <- FALSE
  # overwrite <- TRUE

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
  if (is.null(bwa.path) == FALSE) {
    b.string <- unlist(strsplit(bwa.path, ""))
    if (b.string[length(b.string)] != "/") {
      bwa.path <- paste0(append(b.string, "/"), collapse = "")
    } # end if
  } else {
    bwa.path <- ""
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

  # Quick checks
  if (is.null(read.directory) == TRUE) {
    stop("Please provide input reads.")
  }

  if (file.exists(read.directory) == F) {
    stop("Input reads not found.")
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

  # Creates output directory
  if (dir.exists("logs") == F) {
    dir.create("logs")
  }

  # Read in sample data **** sample is run twice?!
  reads <- list.files(read.directory, recursive = T, full.names = T)
  sample.names <- list.files(read.directory, recursive = F, full.names = F)

  # Resumes file download
  if (overwrite == FALSE) {
    done.files <- list.files(output.directory)
    sample.names <- sample.names[!sample.names %in% done.files]
  }

  if (length(sample.names) == 0) {
    return("no samples remain to analyze.")
  }

  ############################################################################################
  ########### Step 1 #########################################################################
  ##### Start up loop for each sample
  ############################################################################################

  # Sets up multiprocessing
  cl <- makeCluster(threads)
  registerDoParallel(cl)
  mem.cl <- floor(memory / threads)

  # Loops through each locus and does operations on them
  foreach(i = 1:length(sample.names), .packages = c("foreach", "ShortRead")) %dopar% {
    # Runs through each sample
    # for (i in 1:length(sample.names)) {
    #################################################
    ### Part A: prepare for loading and checks
    #################################################
    sample.dir <- paste0(output.directory, "/", sample.names[i])
    if (file.exists(sample.dir) == FALSE) {
      dir.create(sample.dir)
    }

    # Gets the reads for the sample
    sample.reads <- reads[grep(pattern = paste0(sample.names[i], "_"), x = reads)]
    # Checks the Sample column in case already renamed
    if (length(sample.reads) == 0) {
      sample.reads <- reads[grep(pattern = sample.names[i], x = reads)]
    }

    sample.reads <- unique(gsub("_1.f.*|_2.f.*|_3.f.*|-1.f.*|-2.f.*|-3.f.*|_R1_.*|_R2_.*|_R3_.*|_READ1_.*|_READ2_.*|_READ3_.*|_R1.f.*|_R2.f.*|_R3.f.*|-R1.f.*|-R2.f.*|-R3.f.*|_READ1.f.*|_READ2.f.*|_READ3.f.*|-READ1.f.*|-READ2.f.*|-READ3.f.*|_singleton.*|-singleton.*|READ-singleton.*|READ_singleton.*|_READ-singleton.*|-READ_singleton.*|-READ-singleton.*|_READ_singleton.*", "", sample.reads))

    # Returns an error if reads are not found
    if (length(sample.reads) == 0) {
      stop(sample.names[i], " does not have any reads present for files ")
    } # end if statement

    # CReates new directory
    report.path <- paste0("logs/", sample.names[i])
    if (file.exists(report.path) == FALSE) {
      dir.create(report.path)
    }

    for (j in 1:length(sample.reads)) {

      #Gets the reads from an individual lane
      lane.reads <- reads[grep(pattern = paste0(sample.reads[j], "_"), x = reads)]

      # Checks the Sample column in case already renamed
      if (length(lane.reads) == 0) {
        lane.reads <- reads[grep(pattern = sample.reads[j], x = reads)]
      }
      # Returns an error if reads are not found
      if (length(lane.reads) == 0) {
        stop(sample.reads[j], " does not have any reads present for files ")
      } # end if statement

      # Gets lane names
      lane.name <- paste0("Lane_", j)
      lane.dir <- paste0(sample.dir, "/", lane.name)
      dir.create(lane.dir)

      ############################################################################################
      ########### Step 2 #########################################################################
      ##### Create unmapped reads set
      ############################################################################################

      ############################
      # convert fastqs to a sam file
      system(paste0(
        gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
        " FastqToSam -FASTQ ", lane.reads[1], " -FASTQ2 ", lane.reads[2],
        " -OUTPUT ", lane.dir, "/fastqsam.bam",
        " -SAMPLE_NAME ", sample.names[i],
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
      ))

      # Revert the sam to a bam file. Cleans and compresses
      system(paste0(
        gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
        " RevertSam -I ", lane.dir, "/fastqsam.bam -O ", lane.dir, "/revertsam.bam",
        " -SANITIZE true -MAX_DISCARD_FRACTION 0.005",
        " -ATTRIBUTE_TO_CLEAR XT -ATTRIBUTE_TO_CLEAR XN -ATTRIBUTE_TO_CLEAR AS",
        " -ATTRIBUTE_TO_CLEAR OP -SORT_ORDER queryname",
        " -RESTORE_ORIGINAL_QUALITIES true -REMOVE_DUPLICATE_INFORMATION true",
        " -REMOVE_ALIGNMENT_INFORMATION true",
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
      ))

      # Tries to automatically find the read groups from the fasta headers
      if (auto.readgroup == T) {
        # Illumina machines generally
        # @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<sample number | barcode1'+barcode2'>
        all.data <- ShortRead::id(ShortRead::readFastq(lane.reads[1]))[1]
        tmp.data <- as.character(all.data)
        header.data <- unlist(strsplit(tmp.data, ":"))
        header.data <- header.data[1:4]
        names(header.data) <- c("instrument", "run", "flowcell", "lane")

        RGID <- paste0(header.data["flowcell"], ".", header.data["lane"])
        RGPU <- paste0(header.data["flowcell"], ".", header.data["lane"], ".", sample.names[i])

        # Read groups are assigned
        system(paste0(
          gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
          " AddOrReplaceReadGroups -I ", lane.dir, "/revertsam.bam -O ", lane.dir, "/all_reads.bam",
          " -RGSM ", sample.names[i], " -RGPU ", RGPU, " -RGID ", RGID,
          " -RGLB LIB-", sample.names[i], " -RGPL ILLUMINA",
          " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
        ))
      } else {
        # Assign read groups all the same
        system(paste0(
          gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
          " AddOrReplaceReadGroups -I ", lane.dir, "/revertsam.bam -O ", lane.dir, "/all_reads.bam",
          " -RGSM ", sample.names[i], " -RGPU FLOWCELL1.LANE", j, ".", sample.names[i], " -RGID FLOWCELL1.LANE", j,
          " -RGLB LIB-", sample.names[i], " -RGPL ILLUMINA",
          " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
        ))
      } # end if

      # Intermediate files are deleted
      system(paste0("rm ", lane.dir, "/revertsam.bam"))
      system(paste0("rm ", lane.dir, "/fastqsam.bam"))

      print(paste0(sample.names[i], " ", lane.name, " completed BAM creation!"))
    } # end j loop

    print(paste0(sample.names[i], " completed BAM creation!"))
  } # end i loop

  stopCluster(cl)
} # end function


#### END SCRIPT
