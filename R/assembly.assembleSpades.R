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
#'   SPAdes. SPAdes then runs MismatchCorrector after the assembly. This step
#'   adds about a third to the run time and finds almost no extra targets, so
#'   the default is \code{FALSE}. Cannot be \code{TRUE} when
#'   \code{isolate = TRUE}. Numbers in HANDOFF-assembly-speed.md. Default:
#'   \code{FALSE}.

#' @param error.correction logical; if \code{TRUE} SPAdes runs BayesHammer to
#'   correct read errors before it assembles. If \code{FALSE} the function
#'   passes \code{--only-assembler} and SPAdes skips that step. BayesHammer
#'   also takes about a third of the run time. The draft contigs are baits for
#'   \code{assembleBinnedTargets}, which tolerates high divergence, and later
#'   steps map the reads back to call sites. The small loss of base accuracy
#'   therefore does not reach the final sequences. Default: \code{FALSE}.
#'
#' @param isolate logical; if \code{TRUE} passes \code{--isolate} to SPAdes,
#'   recommended for highly covered isolate genomes. Cannot be \code{TRUE} when
#'   \code{mismatch.corrector = TRUE}. Default: \code{FALSE}.
#'
#' @param kmer.values integer vector of k-mer sizes passed to SPAdes with
#'   \code{-k}. Default: \code{c(33, 55, 77, 99, 127)}.
#'
#' @param threads total number of CPU threads. The function divides them
#'   between the samples that run at the same time, so each SPAdes receives
#'   \code{floor(threads / parallel.samples)}. Default: \code{1}.

#' @param retry.failed logical; if \code{TRUE} the function assembles again any
#'   sample that failed, one at a time, with all of \code{threads} and
#'   \code{memory}. A large sample can exhaust its share and still assemble
#'   when it has everything. The retry runs only when
#'   \code{parallel.samples > 1}, because a serial first pass already used
#'   every resource. A skipped sample is not retried: empty or unusable reads
#'   fail whatever the memory. Default: \code{TRUE}.

#' @param parallel.samples number of samples to assemble at the same time.
#'   SPAdes scales poorly above about 8 threads, so several small runs finish a
#'   set sooner than one large run. \code{threads} and \code{memory} are
#'   divided between them. Raise this only when the memory of one run allows it:
#'   every concurrent run holds its own peak. Default: \code{1}, which
#'   assembles the samples one after another. Numbers in
#'   HANDOFF-assembly-speed.md.
#'
#' @param memory total RAM in GB. The function divides it between the samples
#'   that run at the same time, so each SPAdes receives
#'   \code{floor(memory / parallel.samples)}. Default: \code{4}.
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
#'   saved as \code{<assembly.directory>/<sample>.fa}. The SPAdes log of each
#'   sample is copied to \code{logs/sample_logs/<sample>/spades.log}. The log
#'   holds the run time of each stage and the peak memory, and it stays after
#'   \code{clean.up.spades} deletes the working directory.
#'
#' @export

assembleSpades = function(input.reads = NULL,
                          output.directory = "processed-reads/spades-assembly",
                          assembly.directory = "draft-assemblies",
                          spades.path = NULL,
                          mismatch.corrector = FALSE,
                          error.correction = FALSE,
                          isolate = FALSE,
                          kmer.values = c(33,55,77,99,127),
                          threads = 1,
                          parallel.samples = 1,
                          retry.failed = TRUE,
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

  # Creates the log directory. Every sample keeps its spades.log below this.
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

  # --only-assembler turns off BayesHammer. SPAdes then writes no corrected
  # reads, so save.corrected.reads has nothing to keep.
  if (error.correction == TRUE) {
    correction.string = ""
  } else {
    correction.string = "--only-assembler "
    save.corrected.reads = FALSE
  }

  # Sets up the reads. The extension is matched on the file name only, so a
  # directory called "fastq" does not select every file below it.
  reads <- .listFastqFiles(input.reads, recursive = TRUE)

  # The sample is the first directory below input.reads. The name is kept for
  # each read so samples are matched exactly and not by a regular expression.
  read.samples <- gsub(paste0(input.reads, "/"), "", reads, fixed = TRUE)
  read.samples <- gsub("/.*", "", read.samples)
  samples <- unique(read.samples)

  # Skips samples already finished
  if (overwrite == FALSE) {
    done.names <- list.files(assembly.directory, pattern = "\\.fa$")
    done.paths <- file.path(assembly.directory, done.names)
    done.names <- done.names[file.exists(done.paths) & file.info(done.paths)$size > 0]
    samples <- samples[!samples %in% gsub("\\.fa$", "", done.names)]
  }

  if (length(samples) == 0) {
    return(invisible(NULL))
  }

  spades.command = .toolCommand("spades.py", spades.path)

  # Divides the resources between the samples that run at the same time. SPAdes
  # scales poorly above about 8 threads, so several small runs finish a set
  # sooner than one large run. Each concurrent run holds its own peak memory.
  if (!is.numeric(threads) || length(threads) != 1 || !is.finite(threads) || threads < 1 ||
      !is.numeric(memory) || length(memory) != 1 || !is.finite(memory) || memory < 1 ||
      !is.numeric(parallel.samples) || length(parallel.samples) != 1 ||
      !is.finite(parallel.samples) || parallel.samples < 1) {
    stop("threads, memory, and parallel.samples must be positive finite values.")
  }
  threads = floor(threads)
  parallel.samples = min(floor(parallel.samples), length(samples), threads,
                         max(1, floor(memory)))
  thread.cl = floor(threads / parallel.samples)
  mem.cl = floor(memory / parallel.samples)

  if (parallel.samples > 1) {
    print(paste0("Assembling ", parallel.samples, " samples at a time with ",
                 thread.cl, " threads and ", mem.cl, "GB each."))
  }

  #Header data for features and whatnot
  # Assembles one sample. The thread and memory counts are arguments so the
  # retry below can give a failed sample everything.
  assemble.one = function(i, use.threads, use.memory) {
  tryCatch({

    sample.reads = reads[read.samples == samples[i]]

    #Returns an error if reads are not found
    if (length(sample.reads) == 0 ){
      print(paste0(samples[i], " does not have any reads present. Skipping."))
      return("skipped")
    } #end if statement

    #Skip samples with empty or near-empty read files (e.g. all reads removed by decontamination)
    file.sizes = file.info(sample.reads)$size
    if (any(is.na(file.sizes)) || max(file.sizes, na.rm = TRUE) < 1000) {
      print(paste0(samples[i], " read files are empty or near-empty (max file size: ",
                   max(file.sizes, na.rm = TRUE), " bytes). Skipping."))
      return("skipped")
    }

    #Run SPADES on sample
    k.val = paste(kmer.values, collapse = ",")

    #Creates assembly reads folder if not present
    save.assem = paste0(output.directory, "/", samples[i])
    dir.create(save.assem, showWarnings = FALSE, recursive = TRUE)

    #Sorts reads
    read.base = basename(sample.reads)
    read.number = rep(NA_integer_, length(sample.reads))
    read.number[grepl("(_R1|-R1|_READ1|-READ1|READ1|_1|-1)([_.-]|$)", read.base,
                      ignore.case = TRUE)] = 1L
    read.number[grepl("(_R2|-R2|_READ2|-READ2|READ2|_2|-2)([_.-]|$)", read.base,
                      ignore.case = TRUE)] = 2L
    read.number[grepl("(_R3|-R3|_READ3|-READ3|READ3|_3|-3|singleton)", read.base,
                      ignore.case = TRUE)] = 3L
    key.pattern = paste0("(_R[123]|-R[123]|_READ[123]|-READ[123]|READ[123]|",
                         "_[123]|-[123]|_singleton|-singleton|READ.singleton).*$")
    read.key = sub(key.pattern, "", read.base, ignore.case = TRUE)
    read.key = sub("\\.(fastq|fq)(\\.gz)?$", "", read.key, ignore.case = TRUE)
    read.key = sub("[_.-]+$", "", read.key)
    if (any(is.na(read.number)) || any(nchar(read.key) == 0)) {
      stop(samples[i], " has FASTQ files without a supported read-number token: ",
           paste(read.base[is.na(read.number) | nchar(read.key) == 0], collapse = ", "))
    }
    sample.lanes = unique(read.key)

    #Creates a spades character string to run different reads and library configurations.
    #Spades numbers each library type separately, so paired and single-end
    #libraries get their own counters.
    final.read.string = c()
    pe.count = 0
    se.count = 0
    for (j in seq_along(sample.lanes)) {
      #Gets the sample reads
      lib.index = read.key == sample.lanes[j]
      lib.reads = sample.reads[lib.index]

      #Concatenate together
      lib.read1 = lib.reads[read.number[lib.index] == 1]
      lib.read2 = lib.reads[read.number[lib.index] == 2]
      lib.read3 = lib.reads[read.number[lib.index] == 3]

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
        print(paste0(samples[i], ", library ", basename(sample.lanes[j]),
                     ": unsupported read layout (", length(lib.read1), " read1, ",
                     length(lib.read2), " read2, ", length(lib.read3),
                     " merged or singleton files). These reads are not assembled."))
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
      print(paste0(samples[i], " has no reads in a layout spades accepts. Skipping."))
      return("skipped")
    }

    tmp.dir <- paste0(temp.directory, "/spades_", samples[i])
    dir.create(tmp.dir, showWarnings = FALSE, recursive = TRUE)

    #Runs spades command
    spades.status = system(paste0(spades.command, " ", final.read.string,
                  "--tmp-dir ", shQuote(tmp.dir), " -o ", shQuote(save.assem),
                  " -k ", k.val, " ", mismatch.string, correction.string,
                  "-t ", use.threads, " -m ", use.memory),
           ignore.stdout = quiet, ignore.stderr = quiet)

    # Saves the spades log for every sample, not only for a failure. The log
    # gives the run time of each stage and the peak memory, and the working
    # directory may be deleted below.
    sample.log.directory = paste0("logs/sample_logs/", samples[i])
    if (dir.exists(sample.log.directory) == FALSE) {
      dir.create(sample.log.directory, recursive = TRUE, showWarnings = FALSE)
    }
    if (file.exists(paste0(save.assem, "/spades.log")) == TRUE) {
      file.copy(paste0(save.assem, "/spades.log"),
                paste0(sample.log.directory, "/spades.log"), overwrite = TRUE)
    }

    #Warns if spades failed, also copies new assemblies to assembly.directory
    scaffold.file = paste0(save.assem, "/scaffolds.fasta")
    if (spades.status == 0 && file.exists(scaffold.file) == TRUE &&
        file.access(scaffold.file, 4) == 0 && file.size(scaffold.file) > 0) {
      final.file = paste0(assembly.directory, "/", samples[i], ".fa")
      temp.file = paste0(final.file, ".tmp-", Sys.getpid())
      if (file.copy(scaffold.file, temp.file, overwrite = TRUE) == FALSE ||
          file.size(temp.file) == 0 || file.rename(temp.file, final.file) == FALSE) {
        unlink(temp.file)
        stop("Could not publish the SPAdes scaffolds for ", samples[i], ".")
      }
    } else {
      unlink(tmp.dir, recursive = TRUE)
      print(paste0("spades error for ", samples[i],
                   ", check ", sample.log.directory, "/spades.log."))
      return("failed")
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
    return("done")

  }, error = function(e) {
    print(paste0(samples[i], " failed: ", conditionMessage(e)))
    return("failed")
  })
  }#end assemble.one

  results = parallel::mclapply(seq_along(samples),
                               function(i) assemble.one(i, thread.cl, mem.cl),
                               mc.cores = parallel.samples) # end i loop

  # A large sample can exhaust its share of the memory and still assemble when
  # it has all of it. The retry runs the failed samples one at a time with
  # every thread and gigabyte. It runs only when the first pass divided the
  # resources, because a serial first pass already used everything. A skipped
  # sample is not retried: empty or unusable reads fail whatever the memory.
  retry.index = which(vapply(results, function(x) identical(x, "failed"), logical(1)))
  if (retry.failed == TRUE && parallel.samples > 1 && length(retry.index) != 0) {
    print(paste0("Retrying ", length(retry.index), " failed sample(s) one at a time with ",
                 threads, " threads and ", memory, "GB."))
    for (i in retry.index) {
      # Clears the partial working directory so spades starts clean.
      unlink(paste0(output.directory, "/", samples[i]), recursive = TRUE)
      results[[i]] = assemble.one(i, threads, memory)
    }
  }

  # A warning raised in a forked child never reaches the parent, so a failed
  # sample must be counted here. The count is taken after the retry.
  fail.index = which(vapply(results, function(x) identical(x, "failed"), logical(1)))
  if (length(fail.index) != 0) {
    print(paste0(length(fail.index), " of ", length(samples),
                 " samples failed: ", paste(samples[fail.index], collapse = ", "),
                 ". Each has a log in logs/sample_logs/<sample>/spades.log."))
  }

}#end function
