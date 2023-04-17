#' @title mapReference
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

mapReference = function(bam.directory = NULL,
                        output.directory = "variant-calling/sample-mapping",
                        assembly.directory = NULL,
                        check.assemblies = TRUE,
                        reference.file = NULL,
                        samtools.path = NULL,
                        bwa.path = NULL,
                        gatk4.path = NULL,
                        threads = 1,
                        memory = 1,
                        overwrite = FALSE,
                        quiet = TRUE) {
  # Debugging
  # library(PhyloCap)
  # library(foreach)
  # setwd("/Volumes/LaCie/Mantellidae/data-analysis")
  # assembly.directory <- "/Volumes/LaCie/Mantellidae/expanded-assemblies"
  # output.directory <- "variant-calling/sample-mapping"
  # reference.file <- "/Volumes/LaCie/Ultimate_FrogCap/Final_Files/FINAL_marker-seqs_Mar14-2023.fa"
  # bam.directory <- "/Volumes/LaCie/Mantellidae/data-analysis/variant-calling/sample-mapping"

  # gatk4.path <- "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  # samtools.path <- "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  # bwa.path <- "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"

  # check.assemblies = FALSE
  # auto.readgroup <- T
  # threads <- 4
  # memory <- 8
  # quiet <- FALSE
  # overwrite <- FALSE

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
  if (is.null(bam.directory) == TRUE) {
    stop("Please provide the bam directory.")
  }
  if (file.exists(bam.directory) == FALSE) {
    stop("BAM folder not found.")
  }

  # Creates output directory
  if (dir.exists("logs") == FALSE) {
    dir.create("logs")
  }

  # Read in sample data
  bam.files <- list.files(bam.directory, recursive = TRUE, full.names = TRUE)
  bam.files <- bam.files[grep("all_reads.bam$", bam.files)]
  sample.names <- list.dirs(bam.directory, recursive = FALSE, full.names = FALSE)

  #Checks to see if assemblies match to reads
  sample.files <- list.files(assembly.directory)

  #Stops if TRUE
  if (check.assemblies == TRUE) {
    if (length(sample.files) != length(sample.names)) {
      stop("Assembly count does not match raw read count. Ensure that all samples have been assembled or set check.assemblies == FALSE")
    }
  }

  #Removes reads for missing assemblies if FALSE
  if (check.assemblies == FALSE) {
    sample.names = sample.names[sample.names %in% gsub(".fa$|.fasta$", "", sample.files)]
  }

  # Resumes file download
  if (overwrite == FALSE) {
    done.files <- list.files(output.directory, full.names = TRUE, recursive = TRUE)
    done.files <- done.files[grep("final-mapped-all.bam", done.files)]
    done.names <- gsub("/Lane_.*", "", done.files)
    done.names <- unique(gsub(".*\\/", "", done.names))
    sample.names <- sample.names[!sample.names %in% done.names]
  }

  if (length(sample.names) == 0) {
    return("no samples remain to analyze.")
  }

  ############################################################################################
  ########### Step 1 #########################################################################
  ##### Index reference
  ############################################################################################

  dir.create(output.directory)

  for (i in seq_along(sample.files)) {
    # Creates output directory for the sample
    sample.name <- gsub(".fa|.fasta", "", sample.files[i])
    dir.create(paste0(output.directory, "/", sample.name))
    dir.create(paste0(output.directory, "/", sample.name, "/index"))

    system(paste0(
      "cp ", assembly.directory, "/", sample.files[i], " ",
      output.directory, "/", sample.name, "/index/reference.fa"
    ))

    reference.location <- paste0(output.directory, "/", sample.name, "/index/reference.fa")

    # Indexes the reference
    system(paste0(bwa.path, "bwa index -a bwtsw ", reference.location),
      ignore.stderr = quiet, ignore.stdout = quiet
    )

    # Also creates a samtools index
    system(paste0(samtools.path, "samtools faidx ", reference.location))

    system(paste0(
      gatk4.path, "gatk CreateSequenceDictionary --REFERENCE ", reference.location,
      " --OUTPUT ", output.directory, "/", sample.name, "/index/reference.dict",
      " --USE_JDK_DEFLATER true --USE_JDK_INFLATER true"
    ))
  } # end i loop for indexing reference

  ############################################################################################
  ########### Step 1 #########################################################################
  ##### Start up loop for each sample
  ############################################################################################

  # Runs through each sample
  for (i in seq_along(sample.names)) {
    #################################################
    ### Part A: prepare for loading and checks
    #################################################
    sample.dir <- paste0(output.directory, "/", sample.names[i])
    
    # Gets the reads for the sample
    sample.bams <- bam.files[grep(pattern = paste0(sample.names[i], "/"), x = bam.files)]
    
    # Checks the Sample column in case already renamed
    if (length(sample.bams) == 0) {
      sample.bams <- sample.bams[grep(pattern = sample.names[i], x = sample.bams)]
    }

    # Returns an error if reads are not found
    if (length(sample.bams) == 0) {
      stop(sample.names[i], " does not have any reads present for files ")
    } # end if statement

    # CReates new directory
    report.path <- paste0("logs/", sample.names[i])
    if (file.exists(report.path) == FALSE) {
      dir.create(report.path)
    }

    for (j in seq_along(sample.bams)) {
      
      # Gets lane names
      lane.name <- paste0("Lane_", j)
      lane.dir <- paste0(sample.dir, "/", lane.name)

      # Run BWA Mem
      system("set -o pipefail")

      # Piped verison
      tmp.dir <- paste0(lane.dir, "/tmp")
      reference.location = paste0(output.directory, "/", sample.names[i], "/index/reference.fa")
      system(paste0(
        gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
        " SamToFastq -I ", lane.dir, "/all_reads.bam -FASTQ /dev/stdout -TMP_DIR ", tmp.dir,
        " -CLIPPING_ATTRIBUTE XT -CLIPPING_ACTION 2 -INTERLEAVE true -NON_PF true",
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true | ",
        bwa.path, "bwa mem -M -p -t ", threads, " ",
        reference.location, " /dev/stdin | ",
        gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
        " MergeBamAlignment -ALIGNED_BAM /dev/stdin -UNMAPPED_BAM ", lane.dir, "/all_reads.bam",
        " -OUTPUT ", lane.dir, "/cleaned_final.bam",
        " -R ", reference.location, " -CREATE_INDEX true -ADD_MATE_CIGAR true",
        " -CLIP_ADAPTERS false -CLIP_OVERLAPPING_READS true -INCLUDE_SECONDARY_ALIGNMENTS true",
        " -MAX_INSERTIONS_OR_DELETIONS -1 -PRIMARY_ALIGNMENT_STRATEGY MostDistant -ATTRIBUTES_TO_RETAIN XS",
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true",
        " -TMP_DIR ", tmp.dir
      ))

      system(paste0("rm -r ", tmp.dir))
      #system(paste0(samtools.path, "samtools view -H ", lane.dir, "/cleaned_final.bam | grep '@RG'"))

      ############################################################################################
      ########### Step 4 #########################################################################
      ##### Sort Sam, Mark duplicates, Sort Sam again, finalize mapped set ofreads
      ############################################################################################

      if (file.exists(paste0(lane.dir, "/cleaned_final.bam")) != TRUE) {
        stop("Stopped. Something failed after mapping before sorting.")
      }

      # Sort by coordinate for input into MarkDuplicates
      system(paste0(
        gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
        " SortSam -INPUT ", lane.dir, "/cleaned_final.bam",
        " -OUTPUT ", lane.dir, "/cleaned_final_sort.bam",
        " -CREATE_INDEX true -SORT_ORDER coordinate",
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
      ))

      # Marks duplicate reads
      system(paste0(
        gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
        " MarkDuplicates -INPUT ", lane.dir, "/cleaned_final_sort.bam",
        " -OUTPUT ", lane.dir, "/cleaned_final_md.bam",
        " -CREATE_INDEX true -METRICS_FILE ", report.path, "/duplicate_metrics.txt",
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
      ))

      # Sorts and stuff
      system(paste0(
        gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
        " SortSam -INPUT ", lane.dir, "/cleaned_final_md.bam",
        " -OUTPUT /dev/stdout -SORT_ORDER coordinate | ",
        gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
        " SetNmAndUqTags -INPUT /dev/stdin -OUTPUT ", lane.dir, "/final-mapped-all.bam",
        " -CREATE_INDEX true -R ", reference.location,
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
      ))

      if (file.exists(paste0(lane.dir, "/final-mapped-all.bam")) != TRUE) {
        stop("Stopped. Something failed after mapping during sorting.")
      }

      # Intermediate files are deleted
      system(paste0("rm ", lane.dir, "/cleaned_final_sort*"))
      system(paste0("rm ", lane.dir, "/cleaned_final*"))

      print(paste0(sample.names[i], " ", lane.name, " completed read mapping to reference!"))
    } # end j loop

    print(paste0(sample.names[i], " completed read mapping to reference!"))
  } # end i loop
  
} # end function
