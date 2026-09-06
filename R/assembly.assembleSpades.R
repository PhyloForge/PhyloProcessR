#' @title assembleSpades
#'
#' @description Runs SPAdes (\code{spades.py}) on a directory of processed reads
#'   to produce de novo genome assemblies. Each sample subdirectory under
#'   \code{input.reads} is assembled independently with support for multi-lane
#'   and mixed paired/single-end/merged read configurations. The
#'   \code{scaffolds.fasta} output from each sample is copied to
#'   \code{assembly.directory} as \code{<sample>.fa}. Samples for which a
#'   \code{.fa} file already exists in \code{assembly.directory} are skipped
#'   unless \code{overwrite = TRUE}. The \code{--careful} and \code{--isolate}
#'   SPAdes modes cannot be used simultaneously.
#'
#' @param input.reads path to a directory of processed reads. Each sample must
#'   occupy its own subdirectory containing FASTQ files whose names encode read
#'   number and pair identity. The subdirectory name is the sample name and is
#'   matched exactly.
#'
#' @param output.directory path to the directory where per-sample SPAdes working
#'   directories will be written. Default:
#'   \code{"processed-reads/spades-assembly"}.
#'
#' @param assembly.directory path to the directory where the final scaffold
#'   FASTA files (.fa) are copied after assembly. Default:
#'   \code{"draft-assemblies"}.
#'
#' @param spades.path path to the directory containing \code{spades.py}. If
#'   \code{NULL} expected on the system PATH. Default: \code{NULL}.
#'
#' @param mismatch.corrector logical; if \code{TRUE} passes \code{--careful} to
#'   SPAdes to reduce mismatches and indels in the assembly. Cannot be \code{TRUE}
#'   when \code{isolate = TRUE}. Default: \code{TRUE}.
#'
#' @param isolate logical; if \code{TRUE} passes \code{--isolate} to SPAdes,
#'   recommended for highly covered isolate genomes. Cannot be \code{TRUE} when
#'   \code{mismatch.corrector = TRUE}. Default: \code{FALSE}.
#'
#' @param kmer.values integer vector of k-mer sizes passed to SPAdes with
#'   \code{-k}. Default: \code{c(33, 55, 77, 99, 127)}.
#'
#' @param threads number of CPU threads passed to SPAdes with \code{-t}.
#'   Default: \code{1}.
#'
#' @param memory RAM in GB passed to SPAdes with \code{-m}. Default: \code{4}.
#'
#' @param overwrite logical; if \code{TRUE} existing output and assembly
#'   directories are deleted and recreated and all samples are rerun. Default:
#'   \code{FALSE}.
#'
#' @param save.corrected.reads logical; if \code{FALSE} (default) the
#'   \code{corrected/} subdirectory produced by SPAdes is deleted after assembly
#'   to save disk space. Ignored when \code{clean.up.spades = TRUE}. Default:
#'   \code{FALSE}.
#'
#' @param clean.up.spades logical; if \code{TRUE} the entire SPAdes working
#'   directory for each sample is deleted after scaffolds are copied, keeping
#'   only the final \code{.fa} file in \code{assembly.directory}. Overrides
#'   \code{save.corrected.reads}. Default: \code{FALSE}.
#'
#' @param temp.directory path to a dedicated directory for temporary files; NULL defaults to R tempdir().
#' @param quiet logical; if \code{TRUE} SPAdes screen output is suppressed.
#'   Default: \code{TRUE}.
#'
#' @return Invisibly returns nothing. Assembled scaffolds for each sample are
#'   saved as \code{<assembly.directory>/<sample>.fa}.
#'
#' @export

assembleSpades = function(input.reads = NULL,
                          output.directory = "processed-reads/spades-assembly",
                          assembly.directory = "draft-assemblies",
                          spades.path = NULL,
                          mismatch.corrector = TRUE,
                          isolate = FALSE,
                          kmer.values = c(33,55,77,99,127),
                          threads = 1,
                          memory = 4,
                          overwrite = FALSE,
                          save.corrected.reads = FALSE,
                          temp.directory = NULL,
                          clean.up.spades = FALSE,
                          quiet = TRUE) {

  # #debug
  # setwd("/Users/chutter/Dropbox/Research/0_Github/Test-dataset")
  # input.reads = "/Users/chutter/Dropbox/Research/0_Github/Test-dataset/processed-reads/pe-merged-reads"
  # spades.path = "/Users/chutter/miniconda3/bin/spades.py"
  # output.directory = "processed-reads/spades-assembly"
  # assembly.directory = "draft-assemblies"
  # mismatch.corrector = F
  # kmer.values = c(21,33,55,77,99,127)
  # threads = 1
  # memory = 4
  # overwrite = FALSE
  # quiet = TRUE
  # resume = TRUE

  # Same adds to bbmap path
  if (is.null(spades.path) == FALSE) {
    b.string = unlist(strsplit(spades.path, ""))
    if (b.string[length(b.string)] != "/") {
      spades.path = paste0(append(b.string, "/"), collapse = "")
    } # end if
  } else {
    spades.path = ""
  }

  if (is.null(temp.directory) == TRUE) {
    temp.directory = tempdir()
  }

  # Quick checks
  if (is.null(input.reads) == TRUE) {
    stop("Please provide input reads.")
  }
  if (dir.exists(input.reads) == F) {
    stop("Input reads not found.")
  }
  if (is.null(output.directory) == TRUE) {
    stop("Please provide an output directory.")
  }
  if (is.null(assembly.directory) == TRUE) {
    stop("Please provide an contig save directory.")
  }

  # Sets directory and reads
  if (dir.exists(output.directory) == F) {
    dir.create(output.directory, recursive = TRUE)
  } else {
    if (overwrite == TRUE) {
      unlink(output.directory, recursive = TRUE)
      dir.create(output.directory, recursive = TRUE)
    }
  } # end else

  # Sets directory and reads
  if (dir.exists(assembly.directory) == F) {
    dir.create(assembly.directory, recursive = TRUE)
  } else {
    if (overwrite == TRUE) {
      unlink(assembly.directory, recursive = TRUE)
      dir.create(assembly.directory, recursive = TRUE)
    }
  } # end else

  # Creates the log directory. Failed samples keep their spades.log here.
  if (dir.exists("logs/sample_logs") == F){ dir.create("logs/sample_logs", recursive = TRUE) }

  if (isolate == TRUE && mismatch.corrector == TRUE) {
    stop("Both --isolate or --careful (mismatch corrector) can not be used together. Please choose only one.")
  }

  if (mismatch.corrector == FALSE && isolate == FALSE) {
    mismatch.string = ""
  }
  
  if (isolate == TRUE) {
    mismatch.string = "--isolate "
  }
  
  if (mismatch.corrector == TRUE) {
    mismatch.string = "--careful "
  }

  # Sets up the reads. The extension is matched on the file name only, so a
  # directory called "fastq" does not select every file below it.
  files <- list.files(path = input.reads, full.names = T, recursive = T)
  reads <- files[grep(pattern = "fastq|fq|clustS", x = basename(files))]

  # The sample is the first directory below input.reads. The name is kept for
  # each read so samples are matched exactly and not by a regular expression.
  read.samples <- gsub(paste0(input.reads, "/"), "", reads, fixed = TRUE)
  read.samples <- gsub("/.*", "", read.samples)
  samples <- unique(read.samples)

  # Skips samples already finished
  if (overwrite == FALSE) {
    done.names <- list.files(assembly.directory, pattern = "\\.fa$")
    samples <- samples[!samples %in% gsub("\\.fa$", "", done.names)]
  }

  if (length(samples) == 0) {
    stop("No samples to run or incorrect directory.")
  }
  #Header data for features and whatnot
  for (i in seq_along(samples)){

    sample.reads = reads[read.samples == samples[i]]

    #Returns an error if reads are not found
    if (length(sample.reads) == 0 ){
      warning(samples[i], " does not have any reads present. Skipping.")
      next
    } #end if statement

    #Skip samples with empty or near-empty read files (e.g. all reads removed by decontamination)
    file.sizes = file.info(sample.reads)$size
    if (any(is.na(file.sizes)) || max(file.sizes, na.rm = TRUE) < 1000) {
      warning(samples[i], " read files are empty or near-empty (max file size: ",
              max(file.sizes, na.rm = TRUE), " bytes). Skipping.")
      next
    }

    #Run SPADES on sample
    k.val = paste(kmer.values, collapse = ",")

    #Creates assembly reads folder if not present
    save.assem = paste0(output.directory, "/", samples[i])
    dir.create(save.assem, showWarnings = FALSE, recursive = TRUE)

    #Sorts reads
    sample.lanes = unique(gsub("_1.f.*|_2.f.*|_3.f.*|-1.f.*|-2.f.*|-3.f.*|_R1_.*|_R2_.*|_R3_.*|_READ1_.*|_READ2_.*|_READ3_.*|_R1.f.*|_R2.f.*|_R3.f.*|-R1.f.*|-R2.f.*|-R3.f.*|_READ1.f.*|_READ2.f.*|_READ3.f.*|-READ1.f.*|-READ2.f.*|-READ3.f.*|_singleton.*|-singleton.*|READ-singleton.*|READ_singleton.*|_READ-singleton.*|-READ_singleton.*|-READ-singleton.*|_READ_singleton.*", "", sample.reads))

    #Creates a spades character string to run different reads and library configurations.
    #Spades numbers each library type separately, so paired and single-end
    #libraries get their own counters.
    final.read.string = c()
    pe.count = 0
    se.count = 0
    for (j in seq_along(sample.lanes)) {
      #Gets the sample reads
      lib.reads = sample.reads[grep(sample.lanes[j], sample.reads, fixed = TRUE)]

      #Concatenate together
      lib.read1 = lib.reads[grep("_1.f.*|-1.f.*|_R1_.*|-R1_.*|_R1-.*|-R1-.*|READ1.*|_R1.fast.*|-R1.fast.*", lib.reads)]
      lib.read2 = lib.reads[grep("_2.f.*|-2.f.*|_R2_.*|-R2_.*|_R2-.*|-R2-.*|READ2.*|_R2.fast.*|-R2.fast.*", lib.reads)]
      lib.read3 = lib.reads[grep("_3.f.*|-3.f.*|_R3_.*|-R3_.*|_R3-.*|-R3-.*|READ3.*|_R3.fast.*|-R3.fast.*|_READ3.fast.*|-READ3.fast.*|_singleton.*|-singleton.*|READ-singleton.*|READ_singleton.*|_READ-singleton.*|-READ_singleton.*|-READ-singleton.*|_READ_singleton.*", lib.reads)]

      #Checks for different read lengths. Spades needs a paired library before it
      #accepts merged reads, so a lone merged file is given as single-end.
      read.string = ""
      if (length(lib.read1) == 1 && length(lib.read2) == 1) {
        pe.count = pe.count + 1
        read.string = paste0("--pe", pe.count, "-1 ", shQuote(lib.read1),
                             " --pe", pe.count, "-2 ", shQuote(lib.read2), " ")
        if (length(lib.read3) == 1) {
          read.string = paste0(read.string, "--pe", pe.count, "-m ", shQuote(lib.read3), " ")
        }
      } else if (length(lib.read1) == 1 && length(lib.read2) == 0) {
        se.count = se.count + 1
        read.string = paste0("--s", se.count, " ", shQuote(lib.read1), " ")
      } else if (length(lib.read1) == 0 && length(lib.read2) == 0 && length(lib.read3) == 1) {
        se.count = se.count + 1
        read.string = paste0("--s", se.count, " ", shQuote(lib.read3), " ")
      }

      #Warns when a library does not match a supported layout. Spades would drop
      #these reads without a message.
      if (read.string == "" || length(lib.read3) > 1) {
        warning(samples[i], ", library ", basename(sample.lanes[j]),
                ": unsupported read layout (", length(lib.read1), " read1, ",
                length(lib.read2), " read2, ", length(lib.read3),
                " merged or singleton files). These reads are not assembled.")
      }

      final.read.string = paste0(final.read.string, read.string)
    }#end j loop

#     pe.read1 = sample.reads[grep("_READ1", sample.reads)]
#     if (length(pe.read1) != 0){ pe.read1.string = paste0("--pe", rep(1:length(pe.read1)), "-1 ", pe.read1, collapse = " ") }
#     pe.read2 = sample.reads[grep("_READ2", sample.reads)]
#     if (length(pe.read2) != 0){ pe.read2.string = paste0("--pe", rep(1:length(pe.read2)), "-2 ", pe.read2, collapse = " ") }
#     mg.read3 = sample.reads[grep("_READ3", sample.reads)]
#     if (length(mg.read3) != 0){ mg.read3.string = paste0("--pe", rep(1:length(mg.read3)), "-m ", mg.read3, collapse = " ") }

    #Skips the sample when no library produced a usable spades string
    if (length(final.read.string) == 0 || final.read.string == "") {
      warning(samples[i], " has no reads in a layout spades accepts. Skipping.")
      next
    }

    tmp.dir <- paste0(temp.directory, "/spades_", samples[i])
    dir.create(tmp.dir, showWarnings = FALSE, recursive = TRUE)

    #Runs spades command
    system(paste0(spades.path, "spades.py ", final.read.string,
                  "--tmp-dir ", shQuote(tmp.dir), " -o ", shQuote(save.assem),
                  " -k ", k.val, " ", mismatch.string,
                  "-t ", threads, " -m ", memory),
           ignore.stdout = quiet, ignore.stderr = quiet)

    #Warns if spades failed, also copies new assemblies to assembly.directory
    if (file.exists(paste0(save.assem, "/scaffolds.fasta")) == TRUE ){
      file.copy(paste0(save.assem, "/scaffolds.fasta"),
                paste0(assembly.directory, "/", samples[i], ".fa"), overwrite = TRUE)
    } else {
      #Keeps the spades log so the failure can be inspected after clean up.
      if (file.exists(paste0(save.assem, "/spades.log")) == TRUE) {
        file.copy(paste0(save.assem, "/spades.log"),
                  paste0("logs/sample_logs/FAILURE_", samples[i], "_spades.log"),
                  overwrite = TRUE)
      }
      unlink(tmp.dir, recursive = TRUE)
      warning(paste0("spades error for ", samples[i],
                     ", check logs/sample_logs/FAILURE_", samples[i], "_spades.log."))
      next
    }

    if (clean.up.spades == TRUE) {
      unlink(save.assem, recursive = TRUE)
    } else {
      if (save.corrected.reads == FALSE) {
        unlink(paste0(save.assem, "/corrected"), recursive = TRUE)
      }
    }
    unlink(tmp.dir, recursive = TRUE)
    print(paste0(samples[i], " Completed Spades asssembly!"))

  }#end sample loop

}#end function

