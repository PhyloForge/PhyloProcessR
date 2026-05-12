#' @title mapReferenceSample
#'
#' @description Indexes a per-sample reference assembly with BWA and GATK, then
#'   maps per-sample pre-processed BAM files (from prepareBAM()) against that
#'   sample's own reference using the GATK best-practices pipeline (SamToFastq
#'   | bwa mem | MergeBamAlignment | SortSam | MarkDuplicates |
#'   SetNmAndUqTags). Used for per-sample variant calling workflows where each
#'   sample is mapped to its own draft assembly.
#'
#' @param mapping.directory path to the directory containing per-sample
#'   sub-directories with pre-processed BAM files (all_reads.bam) created by
#'   prepareBAM().
#'
#' @param assembly.directory path to the directory of per-sample draft assembly
#'   FASTA files (.fa or .fasta), one file per sample.
#'
#' @param check.assemblies logical; if TRUE an error is raised when the number
#'   of assemblies does not match the number of sample BAM directories. If FALSE
#'   samples without a matching assembly are silently skipped.
#'
#' @param samtools.path system path to the directory containing samtools; NULL
#'   searches the system PATH.
#'
#' @param bwa.path system path to the directory containing bwa; NULL searches
#'   the system PATH.
#'
#' @param gatk4.path system path to the directory containing the gatk
#'   executable; NULL searches the system PATH.
#'
#' @param threads number of CPU threads for BWA and GATK operations.
#'
#' @param memory total RAM in GB to allocate as the JVM heap (-Xmx).
#'
#' @param temp.directory path to a GATK JVM temp directory; NULL uses the
#'   current working directory.
#'
#' @param overwrite logical; if FALSE samples that already have a
#'   final-mapped-all.bam are skipped.
#'
#' @param quiet logical; if TRUE BWA and samtools stdout/stderr are suppressed.
#'
#' @return invisibly; writes final-mapped-all.bam files to per-sample lane
#'   sub-directories in mapping.directory, and per-sample BWA indices to
#'   mapping.directory/sample/index/.
#'
#' @export

mapReferenceSample = function(mapping.directory = NULL,
                              assembly.directory = NULL,
                              check.assemblies = TRUE,
                              samtools.path = NULL,
                              bwa.path = NULL,
                              gatk4.path = NULL,
                              threads = 1,
                              memory = 1,
                              temp.directory = NULL,
                              overwrite = FALSE,
                              quiet = TRUE) {
  # Debugging
  # library(PhyloCap)
  # library(foreach)
  # setwd("/Volumes/LaCie/Mantellidae/data-analysis")
  # assembly.directory <- "/Volumes/LaCie/Mantellidae/expanded-assemblies"
  # mapping.directory <- "variant-calling/sample-mapping"
  # mapping.directory <- "/Volumes/LaCie/Mantellidae/data-analysis/variant-calling/sample-mapping"

  # gatk4.path <- "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  # samtools.path <- "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  # bwa.path <- "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"

  # check.assemblies = FALSE
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
  if (is.null(mapping.directory) == TRUE) {
    stop("Please provide the bam directory.")
  }
  if (file.exists(mapping.directory) == FALSE) {
    stop("BAM folder not found.")
  }

  if (is.null(temp.directory) == TRUE){
    temp.directory = getwd()
  }


  # Creates output directory
  if (dir.exists("logs/sample_logs") == F){ dir.create("logs/sample_logs", recursive = TRUE) }

  # Read in sample data
  bam.files <- list.files(mapping.directory, recursive = TRUE, full.names = TRUE)
  bam.files <- bam.files[grep("all_reads.bam$", bam.files)]
  sample.names <- list.dirs(mapping.directory, recursive = FALSE, full.names = FALSE)

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
    done.files <- list.files(mapping.directory, full.names = TRUE, recursive = TRUE)
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

  dir.create(mapping.directory)

  for (i in seq_along(sample.files)) {
    sample.name <- gsub(".fa$|.fasta$", "", sample.files[i])
    # Skip samples already completed in resume mode
    if (!sample.name %in% sample.names) next
    # Skip if index already exists
    if (dir.exists(paste0(mapping.directory, "/", sample.name, "/index"))) next
    dir.create(paste0(mapping.directory, "/", sample.name), showWarnings = FALSE)
    dir.create(paste0(mapping.directory, "/", sample.name, "/index"))

    system(paste0(
      "cp ", assembly.directory, "/", sample.files[i], " ",
      mapping.directory, "/", sample.name, "/index/reference.fa"
    ))

    reference.location <- paste0(mapping.directory, "/", sample.name, "/index/reference.fa")

    # Indexes the reference
    system(paste0(bwa.path, "bwa index -a bwtsw ", reference.location),
      ignore.stderr = quiet, ignore.stdout = quiet
    )

    # Also creates a samtools index
    system(paste0(samtools.path, "samtools faidx ", reference.location))

    system(paste0(
      gatk4.path, "gatk CreateSequenceDictionary --REFERENCE ", reference.location,
      " --OUTPUT ", mapping.directory, "/", sample.name, "/index/reference.dict",
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
    sample.dir <- paste0(mapping.directory, "/", sample.names[i])

    # Gets the reads for the sample
    sample.bams <- bam.files[grep(pattern = paste0(sample.names[i], "/"), x = bam.files)]

    if (length(sample.bams) == 0) {
      sample.bams <- bam.files[grep(pattern = sample.names[i], x = bam.files)]
    }

    if (length(sample.bams) == 0) {
      stop(sample.names[i], " does not have any reads present for files ")
    } # end if statement

    # CReates new directory
    report.path <- paste0("logs/sample_logs/", sample.names[i])
    if (file.exists(report.path) == FALSE) {
      dir.create(report.path)
    }

    for (j in seq_along(sample.bams)) {

      # Gets lane names
      lane.name <- paste0("Lane_", j)
      lane.dir <- paste0(sample.dir, "/", lane.name)

      tmp.dir <- paste0(lane.dir, "/tmp")
      dir.create(tmp.dir, showWarnings = FALSE)
      reference.location = paste0(mapping.directory, "/", sample.names[i], "/index/reference.fa")
      system(paste0(
        gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
        " SamToFastq -I ", lane.dir, "/all_reads.bam -FASTQ /dev/stdout -TMP_DIR ", tmp.dir,
        " -CLIPPING_ATTRIBUTE XT -CLIPPING_ACTION 2 -INTERLEAVE true -NON_PF true",
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true | ",
        bwa.path, "bwa mem -M -p -t ", threads, " ",
        reference.location, " /dev/stdin | ",
        gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
        " MergeBamAlignment -ALIGNED_BAM /dev/stdin -UNMAPPED_BAM ", lane.dir, "/all_reads.bam",
        " -OUTPUT ", lane.dir, "/cleaned_final.bam",
        " -R ", reference.location, " -CREATE_INDEX true -ADD_MATE_CIGAR true",
        " -CLIP_ADAPTERS false -CLIP_OVERLAPPING_READS true -INCLUDE_SECONDARY_ALIGNMENTS true",
        " -MAX_INSERTIONS_OR_DELETIONS -1 -PRIMARY_ALIGNMENT_STRATEGY MostDistant -ATTRIBUTES_TO_RETAIN XS",
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true",
        " -TMP_DIR ", tmp.dir
      ))

      system(paste0("rm -rf ", tmp.dir))
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
        gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
        " SortSam -INPUT ", lane.dir, "/cleaned_final.bam",
        " -OUTPUT ", lane.dir, "/cleaned_final_sort.bam",
        " -CREATE_INDEX true -SORT_ORDER coordinate",
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
      ))

      # Marks duplicate reads
      system(paste0(
        gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
        " MarkDuplicates -INPUT ", lane.dir, "/cleaned_final_sort.bam",
        " -OUTPUT ", lane.dir, "/cleaned_final_md.bam",
        " -CREATE_INDEX true -METRICS_FILE ", report.path, "/duplicate_metrics.txt",
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
      ))

      # Sorts and stuff
      system(paste0(
        gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
        " SortSam -INPUT ", lane.dir, "/cleaned_final_md.bam",
        " -OUTPUT /dev/stdout -SORT_ORDER coordinate | ",
        gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
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
