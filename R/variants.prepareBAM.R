#' @title prepareBAM
#'
#' @description Converts processed paired-end fastq.gz read files into
#'   unmapped BAM files with correctly assigned read groups, following the GATK
#'   best-practices pre-processing pipeline. For each sample lane, GATK
#'   FastqToSam creates an unmapped BAM, RevertSam cleans it, and
#'   AddOrReplaceReadGroups assigns read group metadata. When auto.readgroup is
#'   TRUE, read group information (flowcell ID and lane) is parsed automatically
#'   from the Illumina fastq headers. Samples are processed in parallel.
#'   The resulting all_reads.bam files are inputs to mapReferenceSample() or
#'   mapReferenceConsensus().
#'
#' @param read.directory path to a directory of processed paired-end fastq.gz
#'   read files, organised in per-sample sub-directories.
#'
#' @param output.directory path to the directory where per-sample per-lane
#'   unmapped BAM sub-directories will be created.
#'
#' @param auto.readgroup logical; if TRUE read group fields (RGID, RGPU) are
#'   automatically extracted from the Illumina fastq header of the first read.
#'   If FALSE generic read group labels are assigned using the lane index.
#'
#' @param samtools.path system path to the directory containing samtools; NULL
#'   searches the system PATH.
#'
#' @param bwa.path system path to the directory containing bwa; NULL searches
#'   the system PATH (currently unused in this step but reserved for pipeline
#'   consistency).
#'
#' @param gatk4.path system path to the directory containing the gatk
#'   executable; NULL searches the system PATH.
#'
#' @param threads number of parallel samples to process simultaneously.
#'
#' @param memory total RAM in GB to allocate as the JVM heap (-Xmx).
#'
#' @param temp.directory path to a GATK JVM temp directory; NULL uses the
#'   current working directory.
#'
#' @param overwrite logical; if TRUE the output directory is deleted and
#'   recreated. If FALSE samples that already have a final-mapped-all.bam are
#'   skipped.
#'
#' @param quiet logical; currently unused.
#'
#' @return invisibly; writes per-sample per-lane all_reads.bam files to
#'   output.directory.
#'
#' @export

prepareBAM = function(read.directory = NULL,
                      output.directory = "sample-mapping",
                      auto.readgroup = TRUE,
                      samtools.path = NULL,
                      bwa.path = NULL,
                      gatk4.path = NULL,
                      threads = 1,
                      memory = 1,
                      temp.directory = NULL,
                      overwrite = FALSE,
                      quiet = TRUE) {

  # library(PhyloCap)
  # library(doParallel)
  # setwd("/Volumes/LaCie/Anax")
  # output.directory <- "data-analysis/joint-genotyping/sample-mapping"
  # read.directory <- "/Volumes/LaCie/Anax/reads"

  # gatk4.path <- "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  # samtools.path <- "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  # bwa.path <- "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"

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

  if (file.exists(read.directory) == FALSE) {
    stop("Input reads not found.")
  }

  if (is.null(temp.directory) == TRUE){
    temp.directory = getwd()
  }

  # Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.directory) == FALSE) {
    dir.create(output.directory)
  } else {
    if (overwrite == TRUE) {
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } # end else

  # Creates output directory
  if (dir.exists("logs/sample_logs") == F){ dir.create("logs/sample_logs", recursive = TRUE) }

  # Read in sample data **** sample is run twice?!
  reads <- list.files(read.directory, recursive = TRUE, full.names = TRUE)
  sample.names <- list.files(read.directory, recursive = FALSE, full.names = FALSE)

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
  ##### Start up loop for each sample
  ############################################################################################

  mem.cl <- floor(memory / threads)

  # Helper: run a GATK command, capture stderr to a log file, and stop with a
  # descriptive message if the exit code is non-zero.  system() does not throw
  # R errors on GATK failures — exit code checking is the only reliable way to
  # detect a crash.
  run.gatk = function(cmd, stderr.log) {
    exit.code = system(paste0(cmd, " 2>> ", stderr.log))
    if (exit.code != 0) {
      gatk.msg = if (file.exists(stderr.log)) {
        tail.lines = tail(readLines(stderr.log), 20)
        paste(tail.lines, collapse = "\n")
      } else { "(no stderr captured)" }
      stop("GATK exited with code ", exit.code, ".\nLast lines of stderr:\n", gatk.msg)
    }
    invisible(exit.code)
  }

  parallel::mclapply(seq_along(sample.names), function(i) {

    sample.id  = sample.names[i]
    sample.dir = paste0(output.directory, "/", sample.id)
    log.file   = paste0("logs/sample_logs/FAILURE_", sample.id, "_prepareBAM.txt")

    tryCatch({

      if (file.exists(sample.dir) == FALSE) { dir.create(sample.dir) }

      sample.reads <- reads[grep(pattern = paste0(sample.id, "_"), x = reads)]
      if (length(sample.reads) == 0) {
        sample.reads <- reads[grep(pattern = sample.id, x = reads)]
      }
      sample.reads <- unique(gsub("_1.f.*|_2.f.*|_3.f.*|-1.f.*|-2.f.*|-3.f.*|_R1_.*|_R2_.*|_R3_.*|_READ1_.*|_READ2_.*|_READ3_.*|_R1.f.*|_R2.f.*|_R3.f.*|-R1.f.*|-R2.f.*|-R3.f.*|_READ1.f.*|_READ2.f.*|_READ3.f.*|-READ1.f.*|-READ2.f.*|-READ3.f.*|_singleton.*|-singleton.*|READ-singleton.*|READ_singleton.*|_READ-singleton.*|-READ_singleton.*|-READ-singleton.*|_READ_singleton.*", "", sample.reads))

      if (length(sample.reads) == 0) {
        stop(sample.id, " has no reads in ", read.directory)
      }

      report.path <- paste0("logs/sample_logs/", sample.id)
      if (file.exists(report.path) == FALSE) { dir.create(report.path) }

      for (j in 1:length(sample.reads)) {

        lane.reads <- reads[grep(pattern = paste0(sample.reads[j], "_"), x = reads)]
        if (length(lane.reads) == 0) {
          lane.reads <- reads[grep(pattern = sample.reads[j], x = reads)]
        }
        if (length(lane.reads) == 0) {
          stop("Lane ", j, " of ", sample.id, " has no reads")
        }

        lane.name  = paste0("Lane_", j)
        lane.dir   = paste0(sample.dir, "/", lane.name)
        stderr.log = paste0(report.path, "/", lane.name, "_gatk-stderr.txt")

        dir.create(lane.dir, showWarnings = FALSE)

        gatk = paste0(gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=",
                      temp.directory, " -Xmx", mem.cl, "G\"")

        # Step 1: FastqToSam
        run.gatk(paste0(gatk, " FastqToSam",
                        " -FASTQ ", lane.reads[1], " -FASTQ2 ", lane.reads[2],
                        " -OUTPUT ", lane.dir, "/fastqsam.bam",
                        " -SAMPLE_NAME ", sample.id,
                        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"),
                 stderr.log)

        # Step 2: RevertSam
        run.gatk(paste0(gatk, " RevertSam",
                        " -I ", lane.dir, "/fastqsam.bam",
                        " -O ", lane.dir, "/revertsam.bam",
                        " -SANITIZE true -MAX_DISCARD_FRACTION 0.005",
                        " -ATTRIBUTE_TO_CLEAR XT -ATTRIBUTE_TO_CLEAR XN -ATTRIBUTE_TO_CLEAR AS",
                        " -ATTRIBUTE_TO_CLEAR OP -SORT_ORDER queryname",
                        " -RESTORE_ORIGINAL_QUALITIES true -REMOVE_DUPLICATE_INFORMATION true",
                        " -REMOVE_ALIGNMENT_INFORMATION true",
                        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"),
                 stderr.log)

        file.remove(paste0(lane.dir, "/fastqsam.bam"))

        # Step 3: AddOrReplaceReadGroups
        if (auto.readgroup == TRUE) {
          all.data    <- ShortRead::id(ShortRead::readFastq(lane.reads[1]))[1]
          header.data <- unlist(strsplit(as.character(all.data), ":"))
          header.data <- header.data[1:4]
          names(header.data) <- c("instrument", "run", "flowcell", "lane")
          RGID <- paste0(header.data["flowcell"], ".", header.data["lane"])
          RGPU <- paste0(RGID, ".", sample.id)
        } else {
          RGID <- paste0("FLOWCELL1.LANE", j)
          RGPU <- paste0(RGID, ".", sample.id)
        }

        run.gatk(paste0(gatk, " AddOrReplaceReadGroups",
                        " -I ", lane.dir, "/revertsam.bam",
                        " -O ", lane.dir, "/all_reads.bam",
                        " -RGSM ", sample.id,
                        " -RGPU ", RGPU, " -RGID ", RGID,
                        " -RGLB LIB-", sample.id, " -RGPL ILLUMINA",
                        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"),
                 stderr.log)

        file.remove(paste0(lane.dir, "/revertsam.bam"))

        # Keep the stderr log only on failure; clean it up on success
        if (file.exists(stderr.log)) { file.remove(stderr.log) }

        print(paste0(sample.id, " ", lane.name, " completed BAM creation!"))

      } # end j loop

      print(paste0(sample.id, " completed BAM creation!"))

    }, error = function(e) {
      # Clean up any partial intermediate files so resume starts clean
      for (f in c("fastqsam.bam", "revertsam.bam")) {
        partial = list.files(sample.dir, pattern = f, recursive = TRUE, full.names = TRUE)
        if (length(partial) > 0) { file.remove(partial) }
      }
      msg = conditionMessage(e)
      writeLines(msg, log.file)
      warning(sample.id, ": prepareBAM failed — see ", log.file)
    })

  }, mc.cores = threads) # end i loop

} # end function


#### END SCRIPT
