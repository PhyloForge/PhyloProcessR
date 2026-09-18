# Internal functions used by screenReads.


.screenSampleNames = function(read.directory = NULL,
                              sample.file = NULL,
                              use.dropbox = FALSE) {

  if (use.dropbox == TRUE) {
    sample.data = read.csv(sample.file, stringsAsFactors = FALSE)

    required.columns = c("File", "Sample")
    if (all(required.columns %in% colnames(sample.data)) == FALSE) {
      stop("sample.file must contain File and Sample columns: ", sample.file)
    }

    sample.names = unique(sample.data$Sample)
  } else {
    sample.data = NULL
    sample.names = .listSampleNames(read.directory)
  }

  sample.names = sample.names[is.na(sample.names) == FALSE]
  sample.names = sample.names[nchar(sample.names) > 0]

  if (length(sample.names) == 0) {
    stop("No sample names were found in the configured read source.")
  }

  return(list(names = sample.names, data = sample.data))
}


.screenLocalReads = function(read.directory = NULL,
                             sample.name = NULL,
                             raw.directory = NULL) {

  sample.directory = file.path(read.directory, sample.name)

  if (dir.exists(sample.directory) == TRUE) {
    input.files = .listFastqFiles(sample.directory)
    return(list(directory = sample.directory,
                files = input.files,
                staged = FALSE))
  }

  local.files = .listFastqFiles(read.directory, recursive = FALSE)
  local.names = basename(local.files)
  input.files = .matchPrefix(local.files, local.names, sample.name)

  input.directory = file.path(raw.directory, sample.name)
  if (dir.exists(input.directory) == TRUE) {
    unlink(input.directory, recursive = TRUE)
  }
  dir.create(input.directory, recursive = TRUE, showWarnings = FALSE)

  output.files = file.path(input.directory, basename(input.files))
  linked = file.symlink(input.files, output.files)

  if (length(linked) > 0 && all(linked) == FALSE) {
    stop("Could not stage the local reads for ", sample.name, ".")
  }

  input.files = .listFastqFiles(input.directory)
  return(list(directory = input.directory,
              files = input.files,
              staged = TRUE))
}


.screenDropboxReads = function(sample.data = NULL,
                               sample.name = NULL,
                               dropbox.directory = NULL,
                               dropbox.token = NULL,
                               raw.directory = NULL) {

  input.directory = file.path(raw.directory, sample.name)
  dir.create(input.directory, recursive = TRUE, showWarnings = FALSE)

  sample.rows = sample.data[sample.data$Sample == sample.name, ]
  temporary.file = tempfile(fileext = ".csv")
  write.csv(sample.rows, temporary.file, row.names = FALSE)

  on.exit(unlink(temporary.file), add = TRUE)

  dropboxDownload(sample.spreadsheet = temporary.file,
                  dropbox.directory = dropbox.directory,
                  dropbox.token = dropbox.token,
                  output.directory = input.directory,
                  overwrite = FALSE,
                  skip.not.found = TRUE)

  input.files = .listFastqFiles(input.directory)
  return(list(directory = input.directory,
              files = input.files,
              staged = FALSE))
}


.screenGzipFiles = function(input.files = NULL) {

  gzip.files = input.files[grepl("\\.gz$", input.files, ignore.case = TRUE)]
  bad.files = character()

  for (gzip.file in gzip.files) {
    command = paste0("gzip -t ", shQuote(gzip.file))
    status = system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)

    if (status != 0) {
      bad.files = c(bad.files, gzip.file)
    }
  }

  return(bad.files)
}


.screenCompletePairs = function(input.files = NULL) {

  lane.prefixes = .stripReadSuffix(input.files)

  if (length(lane.prefixes) == 0) {
    return(FALSE)
  }

  for (lane.prefix in lane.prefixes) {
    lane.files = .matchPrefix(input.files, input.files, lane.prefix)
    ordered.files = tryCatch(.orderReadPair(lane.files),
                             error = function(e) NULL)

    if (is.null(ordered.files) == TRUE) {
      return(FALSE)
    }
  }

  return(TRUE)
}


.screenCleanLocalStage = function(input.directory = NULL,
                                   staged = FALSE) {

  if (staged == TRUE && dir.exists(input.directory) == TRUE) {
    unlink(input.directory, recursive = TRUE)
  }

  return(invisible(NULL))
}


.screenRemoveSummarySample = function(summary.file = NULL,
                                       sample.name = NULL) {

  if (file.exists(summary.file) == FALSE) {
    return(invisible(NULL))
  }

  summary.data = tryCatch(read.csv(summary.file, stringsAsFactors = FALSE),
                          error = function(e) NULL)

  if (is.null(summary.data) == TRUE ||
      "Sample" %in% colnames(summary.data) == FALSE) {
    return(invisible(NULL))
  }

  summary.data = summary.data[summary.data$Sample != sample.name, , drop = FALSE]
  write.csv(summary.data, summary.file, row.names = FALSE)

  return(invisible(NULL))
}


.screenResetIncompleteSample = function(sample.name = NULL,
                                         cleaned.directory = NULL,
                                         log.directory = NULL,
                                         completion.file = NULL) {

  cleaned.sample.directory = file.path(cleaned.directory, sample.name)
  if (dir.exists(cleaned.sample.directory) == TRUE) {
    unlink(cleaned.sample.directory, recursive = TRUE)
  }

  fastp.files = list.files(log.directory,
                           pattern = "_fastp-clean",
                           full.names = TRUE)
  if (length(fastp.files) > 0) {
    unlink(fastp.files)
  }

  capture.metadata = list.files(log.directory,
                                pattern = "_capture-metadata\\.csv$",
                                full.names = TRUE)
  if (length(capture.metadata) > 0) {
    unlink(capture.metadata)
  }

  if (file.exists(completion.file) == TRUE) {
    unlink(completion.file)
  }

  summary.files = c("logs/fastp_summary.csv",
                    "logs/sample-capture-assessment_summary.csv",
                    "logs/barcodeSampleScan_summary.csv",
                    "logs/X0_fastq-stats_rolling.csv")

  for (summary.file in summary.files) {
    .screenRemoveSummarySample(summary.file, sample.name)
  }

  return(invisible(NULL))
}


.screenAggregateFastq = function(fastq.data = NULL) {

  if (nrow(fastq.data) == 0) {
    return(data.frame())
  }

  summary.data = aggregate(cbind(Read1_Count,
                                 Read2_Count,
                                 Read3_Count,
                                 Total_Reads,
                                 Read_Pairs,
                                 MegaBasePairs) ~ Sample,
                           data = fastq.data,
                           FUN = sum)

  summary.data$Read_Length = fastq.data$Read_Length[1]
  total.reads = summary.data$Read1_Count
  total.reads = total.reads + summary.data$Read2_Count
  total.reads = total.reads + summary.data$Read3_Count
  summary.data$Reads_Per_Million = total.reads / 1000000

  return(summary.data)
}


.screenAggregateFastp = function(summary.file = NULL) {

  if (file.exists(summary.file) == FALSE) {
    return(data.frame())
  }

  fastp.data = read.csv(summary.file, stringsAsFactors = FALSE)
  summary.data = aggregate(cbind(startPairs,
                                 removePairs,
                                 endPairs) ~ Sample,
                           data = fastp.data,
                           FUN = sum)

  removed.fraction = summary.data$removePairs / summary.data$startPairs
  summary.data$pctRemovedByFastp = round(removed.fraction * 100, 2)

  columns = c("Sample",
              "startPairs",
              "removePairs",
              "endPairs",
              "pctRemovedByFastp")
  summary.data = summary.data[, columns]

  return(summary.data)
}


.screenMergeTables = function(table.list = NULL) {

  keep.tables = list()

  for (table.data in table.list) {
    if (is.data.frame(table.data) == TRUE && nrow(table.data) > 0) {
      keep.tables[[length(keep.tables) + 1]] = table.data
    }
  }

  merged.data = keep.tables[[1]]

  if (length(keep.tables) > 1) {
    for (i in 2:length(keep.tables)) {
      merged.data = merge(merged.data,
                          keep.tables[[i]],
                          by = "Sample",
                          all = TRUE)
    }
  }

  return(merged.data)
}


#' @title screenReads
#'
#' @description Screens raw sequence-capture reads one sample at a time. The
#'   function counts raw reads, cleans them with fastp, measures capture
#'   efficiency, and writes one combined sample summary. It can use local reads
#'   or download each sample from Dropbox. Processing one sample at a time
#'   limits temporary disk use.
#'
#' @param read.directory directory that contains local FASTQ files. Files can
#'   be flat or stored in one directory per sample.
#' @param processed.reads directory used for temporary raw and cleaned reads.
#' @param use.dropbox logical; download reads from Dropbox when TRUE.
#' @param sample.file CSV file with File and Sample columns for Dropbox input.
#' @param dropbox.directory Dropbox directory that contains the reads.
#' @param dropbox.token saved Dropbox token file.
#' @param delete.raw.reads logical; delete downloaded raw reads after a sample
#'   completes. Local source reads are never deleted.
#' @param delete.cleaned.reads logical; delete cleaned reads after a sample
#'   completes.
#' @param target.fasta target probe or marker FASTA file.
#' @param read.length expected read length used by code{fastqStats}.
#' @param remove.adaptors logical passed to code{fastpClean}.
#' @param remove.duplicate.reads logical passed to code{fastpClean}.
#' @param error.correction logical passed to code{fastpClean}.
#' @param quality.trim.reads logical passed to code{fastpClean}.
#' @param quality.filter logical passed to code{fastpClean}.
#' @param low.complexity.filter logical passed to code{fastpClean}.
#' @param trim.poly.x logical passed to code{fastpClean}.
#' @param min.read.length minimum cleaned read length passed to
#'   code{fastpClean}.
#' @param run.barcode.scan logical; run MItoTrawlR barcode identification.
#' @param barcode.fasta barcode reference FASTA file.
#' @param barcode.database.fasta local barcode database, or NULL for a remote
#'   search.
#' @param barcode.hits.per.sample number of barcode hits to keep.
#' @param barcode.min.iterations minimum barcode assembly iterations.
#' @param barcode.max.iterations maximum barcode assembly iterations.
#' @param barcode.min.ref.id starting barcode recruitment identity.
#' @param barcode.per.max.length barcode maximum-length allowance.
#' @param fastp.path path to fastp or its directory.
#' @param bwa.path path to BWA or its directory.
#' @param samtools.path path to samtools or its directory.
#' @param bbmap.path path to BBMap or its directory.
#' @param spades.path path to SPAdes or its directory.
#' @param cap3.path path to CAP3 or its directory.
#' @param blast.path path to BLAST or its directory.
#' @param threads number of CPU threads.
#' @param memory memory allocation in GB.
#' @param quiet logical; hide external program output when TRUE.
#'
#' @return Invisibly returns the final sample summary. The function writes the
#'   summary to code{logs/X0_read-screening_FINAL.csv}.
#'
#' @export

screenReads = function(read.directory = NULL,
                       processed.reads = "processed-reads",
                       use.dropbox = FALSE,
                       sample.file = NULL,
                       dropbox.directory = NULL,
                       dropbox.token = NULL,
                       delete.raw.reads = TRUE,
                       delete.cleaned.reads = TRUE,
                       target.fasta = NULL,
                       read.length = 150,
                       remove.adaptors = TRUE,
                       remove.duplicate.reads = TRUE,
                       error.correction = TRUE,
                       quality.trim.reads = FALSE,
                       quality.filter = TRUE,
                       low.complexity.filter = TRUE,
                       trim.poly.x = TRUE,
                       min.read.length = 60,
                       run.barcode.scan = FALSE,
                       barcode.fasta = NULL,
                       barcode.database.fasta = NULL,
                       barcode.hits.per.sample = 5,
                       barcode.min.iterations = 3,
                       barcode.max.iterations = 10,
                       barcode.min.ref.id = 0.70,
                       barcode.per.max.length = 0.50,
                       fastp.path = NULL,
                       bwa.path = NULL,
                       samtools.path = NULL,
                       bbmap.path = NULL,
                       spades.path = NULL,
                       cap3.path = NULL,
                       blast.path = NULL,
                       threads = 1,
                       memory = 8,
                       quiet = TRUE) {

  if (is.null(target.fasta) == TRUE || file.exists(target.fasta) == FALSE) {
    stop("Target FASTA file not found: ", target.fasta)
  }

  if (use.dropbox == TRUE &&
      (is.null(sample.file) == TRUE || file.exists(sample.file) == FALSE)) {
    stop("Sample file not found: ", sample.file)
  }

  if (use.dropbox == FALSE &&
      (is.null(read.directory) == TRUE || dir.exists(read.directory) == FALSE)) {
    stop("Read directory not found: ", read.directory)
  }

  dir.create(processed.reads, showWarnings = FALSE)
  dir.create("logs/sample_logs", recursive = TRUE, showWarnings = FALSE)
  dir.create("sample-capture-assessment", showWarnings = FALSE)

  if (run.barcode.scan == TRUE) {
    dir.create("barcode-assessment", showWarnings = FALSE)
  }

  raw.directory = file.path(processed.reads, "raw-reads")
  cleaned.directory = file.path(processed.reads, "cleaned-reads")
  dir.create(raw.directory, recursive = TRUE, showWarnings = FALSE)
  dir.create(cleaned.directory, recursive = TRUE, showWarnings = FALSE)

  sample.source = .screenSampleNames(read.directory = read.directory,
                                     sample.file = sample.file,
                                     use.dropbox = use.dropbox)
  sample.names = sample.source$names
  sample.data = sample.source$data

  cat("Found", length(sample.names), "unique samples.\n")

  sample.status = data.frame(Sample = sample.names,
                             Status = rep("pending", length(sample.names)),
                             stringsAsFactors = FALSE)
  all.fastq.stats = data.frame()
  target.md5 = unname(tools::md5sum(target.fasta))

  for (i in seq_along(sample.names)) {
    sample.name = sample.names[i]
    cat("\nSample", i, "of", length(sample.names), ":", sample.name, "\n")

    assessment.directory = file.path("sample-capture-assessment", sample.name)
    log.directory = file.path("logs", "sample_logs", sample.name)
    completion.file = file.path(log.directory,
                                paste0(sample.name, "_X0-complete.csv"))

    completion.matches = FALSE

    if (file.exists(completion.file) == TRUE) {
      completion.data = tryCatch(read.csv(completion.file,
                                          stringsAsFactors = FALSE),
                                 error = function(e) NULL)

      if (is.null(completion.data) == FALSE) {
        completion.matches = identical(as.character(completion.data$TargetMD5),
                                       target.md5)
      }
    }

    target.files = list.files(assessment.directory,
                              pattern = "_per-target-counts\\.csv$",
                              full.names = TRUE)

    if (completion.matches == TRUE && length(target.files) > 0) {
      cat("Already complete. Skipping sample.\n")
      sample.status$Status[sample.status$Sample == sample.name] = "complete"

      rolling.file = "logs/X0_fastq-stats_rolling.csv"
      if (file.exists(rolling.file) == TRUE) {
        rolling.data = read.csv(rolling.file, stringsAsFactors = FALSE)
        sample.rows = rolling.data[rolling.data$Sample == sample.name, ]
        all.fastq.stats = rbind(all.fastq.stats, sample.rows)
      }

      next
    }

    if (dir.exists(assessment.directory) == TRUE &&
        completion.matches == FALSE) {
      unlink(assessment.directory, recursive = TRUE)

      metadata.files = list.files(log.directory,
                                  pattern = "_capture-metadata\\.csv$",
                                  full.names = TRUE)
      if (length(metadata.files) > 0) {
        unlink(metadata.files)
      }
    }

    # A sample without a valid X0 completion file is incomplete. Remove only
    # that sample's temporary cleaning results before the retry. This prevents
    # stale fastp metadata from stopping a resumed workflow.
    .screenResetIncompleteSample(sample.name = sample.name,
                                 cleaned.directory = cleaned.directory,
                                 log.directory = log.directory,
                                 completion.file = completion.file)

    if (use.dropbox == TRUE) {
      input.data = .screenDropboxReads(sample.data = sample.data,
                                       sample.name = sample.name,
                                       dropbox.directory = dropbox.directory,
                                       dropbox.token = dropbox.token,
                                       raw.directory = raw.directory)
    } else {
      input.data = .screenLocalReads(read.directory = read.directory,
                                     sample.name = sample.name,
                                     raw.directory = raw.directory)
    }

    input.directory = input.data$directory
    input.files = input.data$files
    local.stage = input.data$staged

    if (length(input.files) == 0) {
      warning(sample.name, " does not have read files. Skipping.")
      sample.status$Status[sample.status$Sample == sample.name] =
        "failed: no reads"
      .screenCleanLocalStage(input.directory, local.stage)
      next
    }

    bad.files = .screenGzipFiles(input.files)
    if (length(bad.files) > 0) {
      warning(sample.name, " has a corrupt or incomplete gzip file. Skipping.")
      dir.create(log.directory, recursive = TRUE, showWarnings = FALSE)
      failure.file = file.path(log.directory,
                               paste0("FAILURE_", sample.name,
                                      "_corrupted-download.txt"))
      failure.message = paste("Corrupt or incomplete gzip file(s):",
                              paste(basename(bad.files), collapse = ", "))
      writeLines(failure.message, failure.file)

      if (use.dropbox == TRUE) {
        file.remove(bad.files)
      }

      sample.status$Status[sample.status$Sample == sample.name] =
        "failed: corrupt reads"
      .screenCleanLocalStage(input.directory, local.stage)
      next
    }

    complete.pairs = .screenCompletePairs(input.files)
    if (complete.pairs == FALSE) {
      warning(sample.name, " has an incomplete read pair. Skipping.")
      sample.status$Status[sample.status$Sample == sample.name] =
        "failed: incomplete read pair"
      .screenCleanLocalStage(input.directory, local.stage)
      next
    }

    fastq.data = fastqStats(read.directory = input.directory,
                            output.name = "fastq-stats-temp",
                            read.length = read.length,
                            threads = threads,
                            mem = memory,
                            overwrite = TRUE)

    fastq.data = fastq.data[fastq.data$Sample == sample.name, ]
    all.fastq.stats = rbind(all.fastq.stats, fastq.data)

    if (file.exists("fastq-stats-temp.csv") == TRUE) {
      unlink("fastq-stats-temp.csv")
    }

    fastpClean(input.reads = input.directory,
               output.directory = cleaned.directory,
               remove.adaptors = remove.adaptors,
               remove.duplicate.reads = remove.duplicate.reads,
               error.correction = error.correction,
               quality.trim.reads = quality.trim.reads,
               quality.filter = quality.filter,
               low.complexity.filter = low.complexity.filter,
               trim.poly.x = trim.poly.x,
               min.read.length = min.read.length,
               fastp.path = fastp.path,
               threads = threads,
               mem = memory,
               overwrite = FALSE,
               quiet = quiet)

    cleaned.sample.directory = file.path(cleaned.directory, sample.name)
    cleaned.files = .listFastqFiles(cleaned.sample.directory)

    if (length(cleaned.files) < 2) {
      warning(sample.name, " did not produce a cleaned read pair. Skipping.")
      sample.status$Status[sample.status$Sample == sample.name] =
        "failed: incomplete cleaned read pair"
      .screenCleanLocalStage(input.directory, local.stage)
      next
    }

    if (run.barcode.scan == TRUE) {
      if (requireNamespace("MItoTrawlR", quietly = TRUE) == FALSE) {
        stop("run.barcode.scan = TRUE requires the MItoTrawlR package.")
      }

      MItoTrawlR::barcodeSampleScan(
        input.reads = cleaned.sample.directory,
        output.directory = "barcode-assessment",
        barcode.fasta = barcode.fasta,
        database.fasta = barcode.database.fasta,
        hits.per.sample = barcode.hits.per.sample,
        per.max.length = barcode.per.max.length,
        min.iterations = barcode.min.iterations,
        max.iterations = barcode.max.iterations,
        min.ref.id = barcode.min.ref.id,
        bbmap.path = bbmap.path,
        spades.path = spades.path,
        cap3.path = cap3.path,
        blast.path = blast.path,
        memory = memory,
        threads = threads,
        overwrite = FALSE,
        quiet = quiet
      )
    }

    assessCaptureEfficiency(
      input.reads = cleaned.sample.directory,
      output.directory = "sample-capture-assessment",
      target.fasta = target.fasta,
      bwa.path = bwa.path,
      samtools.path = samtools.path,
      threads = threads,
      mem = memory,
      overwrite = FALSE,
      quiet = quiet
    )

    if (use.dropbox == TRUE && delete.raw.reads == TRUE) {
      unlink(input.directory, recursive = TRUE)
    }

    dir.create(log.directory, recursive = TRUE, showWarnings = FALSE)
    completion.data = data.frame(Sample = sample.name,
                                 TargetMD5 = target.md5,
                                 stringsAsFactors = FALSE)
    write.csv(completion.data, completion.file, row.names = FALSE)
    sample.status$Status[sample.status$Sample == sample.name] = "complete"

    if (delete.cleaned.reads == TRUE) {
      unlink(cleaned.sample.directory, recursive = TRUE)
    }

    .screenCleanLocalStage(input.directory, local.stage)

    write.csv(all.fastq.stats,
              "logs/X0_fastq-stats_rolling.csv",
              row.names = FALSE)
    cat("Sample", sample.name, "complete.\n")
  }

  fastq.summary = .screenAggregateFastq(all.fastq.stats)
  fastp.summary = .screenAggregateFastp("logs/fastp_summary.csv")

  capture.summary = data.frame()
  capture.file = "logs/sample-capture-assessment_summary.csv"
  if (file.exists(capture.file) == TRUE) {
    capture.summary = read.csv(capture.file, stringsAsFactors = FALSE)
  }

  barcode.summary = data.frame()
  barcode.file = "logs/barcodeSampleScan_summary.csv"
  if (run.barcode.scan == TRUE && file.exists(barcode.file) == TRUE) {
    barcode.summary = read.csv(barcode.file, stringsAsFactors = FALSE)
  }

  summary.tables = list(sample.status,
                        fastq.summary,
                        fastp.summary,
                        capture.summary,
                        barcode.summary)
  final.summary = .screenMergeTables(summary.tables)

  wanted.columns = c("Sample",
                     "Status",
                     "Read_Pairs",
                     "MegaBasePairs",
                     "Reads_Per_Million",
                     "Read_Length",
                     "startPairs",
                     "removePairs",
                     "endPairs",
                     "pctRemovedByFastp",
                     "readPairs",
                     "mappedReads",
                     "targetsHit",
                     "totalTargets",
                     "pctTargetsHit",
                     "pctReadsOnTarget",
                     "MappedReads",
                     "ContigLength",
                     "BestMatch",
                     "Pident",
                     "AlignLength",
                     "Evalue",
                     "Bitscore")
  wanted.columns = wanted.columns[wanted.columns %in% colnames(final.summary)]
  final.summary = final.summary[, wanted.columns, drop = FALSE]

  output.file = "logs/X0_read-screening_FINAL.csv"
  write.csv(final.summary, output.file, row.names = FALSE)
  cat("Final summary written to", output.file, "\n")

  return(invisible(final.summary))
}
