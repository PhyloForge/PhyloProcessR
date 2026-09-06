# Internal helpers for the workflow 1 preprocess functions. These functions are
# not exported. They hold the logic that was previously repeated in each
# preprocess function: program lookup, safe shell calls, sample matching, read
# counting, and summary log writing.


# Builds the shell command for an external program. program.path can be the
# directory that holds the executable, the full path to the executable, or NULL
# to use the system PATH. The function stops when the program is not found.
# Length of a sequence without its Ns. Two fragments of one target are joined with
# N padding, so the padding would otherwise count as recovered sequence.
.baseWidth = function(seqs = NULL) {
  if (length(seqs) == 0) return(integer(0))
  return(as.integer(Biostrings::width(seqs) -
                    Biostrings::letterFrequency(seqs, "N")))
}#end .baseWidth


.toolCommand = function(program = NULL,
                        program.path = NULL) {

  if (is.null(program.path) == TRUE || nchar(program.path) == 0) {
    found.path = Sys.which(program)
    if (nchar(found.path) == 0) {
      stop(program, " was not found on the system PATH. Set the program path in the configuration file.")
    }
    return(shQuote(as.character(found.path)))
  }

  # Accepts either the executable itself or the directory that holds it
  program.path = sub("/+$", "", program.path)
  if (basename(program.path) == program) {
    command.path = program.path
  } else {
    command.path = file.path(program.path, program)
  }

  if (file.exists(command.path) == FALSE) {
    stop(program, " was not found at ", command.path, ". Check the program path in the configuration file.")
  }

  return(shQuote(command.path))
}#end .toolCommand


# Runs a shell command and stops when it fails. A failed external tool used to
# pass unnoticed and cause a confusing error in a later step.
# Set keep.stdout to TRUE when the command sends its own output to a file. R adds
# "> /dev/null" for ignore.stdout, and that redirection wins over the one in the
# command, which leaves an empty output file.
.runCommand = function(command = NULL,
                       quiet = TRUE,
                       task = "external command",
                       keep.stdout = FALSE) {

  status = system(command,
                  ignore.stdout = quiet && keep.stdout == FALSE,
                  ignore.stderr = quiet)
  if (status != 0) {
    stop("The ", task, " step failed with exit status ", status, ".\nCommand: ", command)
  }
  return(invisible(status))
}#end .runCommand


# Runs a shell command and returns the standard output. Stops when the command
# fails.
.runCommandOutput = function(command = NULL,
                             task = "external command") {

  result = suppressWarnings(system(command, intern = TRUE))
  status = attr(result, "status")
  if (is.null(status) == FALSE && status != 0) {
    stop("The ", task, " step failed with exit status ", status, ".\nCommand: ", command)
  }
  return(result)
}#end .runCommandOutput


# Returns the file paths relative to a base directory.
.relativePaths = function(file.paths = NULL,
                          base.dir = NULL) {

  base.dir = sub("/+$", "", base.dir)
  return(sub(paste0("^", base.dir, "/+"), "", file.paths, fixed = FALSE))
}#end .relativePaths


# Selects the files whose match string starts with a prefix and a separator.
# The comparison uses fixed strings, so a regular expression character in a
# sample name is safe. The separator stops a short name from matching a longer
# one, for example Sample1 and Sample10.
.matchPrefix = function(file.paths = NULL,
                        match.strings = NULL,
                        prefix = NULL) {

  keep = startsWith(match.strings, paste0(prefix, "_")) |
         startsWith(match.strings, paste0(prefix, "-")) |
         startsWith(match.strings, paste0(prefix, ".")) |
         startsWith(match.strings, paste0(prefix, "/"))

  return(file.paths[keep])
}#end .matchPrefix


# Removes the read and mate suffix from a read file path. The result is the lane
# prefix, for example Sample_L001. Multiple lanes give multiple prefixes.
.stripReadSuffix = function(file.paths = NULL) {

  read.suffix = paste0("_1.f.*|_2.f.*|_3.f.*|-1.f.*|-2.f.*|-3.f.*|_R1_.*|_R2_.*|_R3_.*|",
                       "_READ1_.*|_READ2_.*|_READ3_.*|_R1.f.*|_R2.f.*|_R3.f.*|",
                       "-R1.f.*|-R2.f.*|-R3.f.*|_READ1.f.*|_READ2.f.*|_READ3.f.*|",
                       "-READ1.f.*|-READ2.f.*|-READ3.f.*|_singleton.*|-singleton.*|",
                       "READ-singleton.*|READ_singleton.*|_READ-singleton.*|",
                       "-READ_singleton.*|-READ-singleton.*|_READ_singleton.*")

  return(unique(gsub(read.suffix, "", file.paths)))
}#end .stripReadSuffix


# Lists the sample names in a read directory. Samples are sub-directories when
# they exist, otherwise the file names give the sample names.
.listSampleNames = function(read.directory = NULL) {

  sample.names = list.dirs(read.directory, recursive = FALSE, full.names = FALSE)

  if (length(sample.names) == 0) {
    sample.names = list.files(read.directory, recursive = FALSE, full.names = FALSE)
    sample.names = unique(gsub("_L00.*", "", sample.names))
    sample.names = sample.names[nchar(sample.names) > 0]
  }

  return(sample.names)
}#end .listSampleNames


# Counts the reads in a fastq file. The file can be compressed or uncompressed.
.countFastqReads = function(read.file = NULL) {

  if (is.na(read.file) == TRUE) { return(NA_real_) }
  if (file.exists(read.file) == FALSE) { return(NA_real_) }

  if (grepl("\\.gz$", read.file) == TRUE) {
    command = paste0("gzip -cd ", shQuote(read.file), " | wc -l")
  } else {
    command = paste0("wc -l < ", shQuote(read.file))
  }

  line.count = as.numeric(.runCommandOutput(command, task = "read counting"))
  return(line.count / 4)
}#end .countFastqReads


# Reads the read counts from a fastp JSON report. This replaces a second pass
# over each fastq file with gzip, which doubled the read and write load.
# fastp counts both mates, so the totals are divided by 2 to give read pairs.
# Merged reads are single sequences and are not divided.
.fastpReadCounts = function(json.file = NULL) {

  if (file.exists(json.file) == FALSE) { return(NULL) }

  report = jsonlite::fromJSON(json.file)
  if (is.null(report$summary$before_filtering$total_reads) == TRUE) { return(NULL) }

  counts = list(startPairs = report$summary$before_filtering$total_reads / 2,
                endPairs = report$summary$after_filtering$total_reads / 2,
                mergedReads = NA_real_)

  # fastp reports merged reads in a separate block when --merge is used
  if (is.null(report$merged_and_filtered$total_reads) == FALSE) {
    counts$mergedReads = report$merged_and_filtered$total_reads
    counts$endPairs = counts$startPairs - counts$mergedReads
  }

  return(counts)
}#end .fastpReadCounts


# Writes a summary log and keeps the rows from earlier runs. Rows for the
# samples in this run replace the older rows. This lets a resumed run build one
# complete log instead of a log that holds only the last batch of samples.
.appendSummary = function(summary.data = NULL,
                          out.csv = NULL) {

  if (is.null(summary.data) == TRUE || nrow(summary.data) == 0) { return(invisible(NULL)) }

  if (file.exists(out.csv) == TRUE) {
    existing = read.csv(out.csv, stringsAsFactors = FALSE)
    if (identical(sort(names(existing)), sort(names(summary.data))) == TRUE) {
      existing = existing[!existing$Sample %in% summary.data$Sample, , drop = FALSE]
      summary.data = rbind(existing[, names(summary.data), drop = FALSE], summary.data)
    }
  }

  write.csv(summary.data, file = out.csv, row.names = FALSE)
  return(invisible(summary.data))
}#end .appendSummary


# Replaces the characters that are unsafe in a file name or a CSV field. Sample
# names reach file paths, so a space or a comma breaks the rename table.
.sanitizeName = function(sample.name = NULL) {

  clean.name = gsub("[^A-Za-z0-9._-]+", "_", trimws(as.character(sample.name)))
  clean.name = gsub("_+", "_", clean.name)
  clean.name = gsub("^_|_$", "", clean.name)
  return(clean.name)
}#end .sanitizeName


# Deletes a directory and creates it again. This replaces system("rm -r"), which
# failed on a path that holds a space.
.resetDirectory = function(output.directory = NULL) {

  if (dir.exists(output.directory) == TRUE) { unlink(output.directory, recursive = TRUE) }
  dir.create(output.directory, recursive = TRUE, showWarnings = FALSE)
  return(invisible(output.directory))
}#end .resetDirectory


# Reports whether every output file of a lane is present. An interrupted lane is
# processed again instead of being skipped for good.
# require.size adds the test that each file holds data. Set it to FALSE where an
# empty output file is a valid result, for example the merged read file of a
# sample whose read pairs do not overlap.
.laneComplete = function(output.files = NULL,
                         require.size = TRUE) {

  if (length(output.files) == 0) { return(FALSE) }
  if (all(file.exists(output.files)) == FALSE) { return(FALSE) }
  if (require.size == FALSE) { return(TRUE) }
  return(all(file.info(output.files)$size > 0))
}#end .laneComplete


# Runs one fastp cleaning step over every sample lane in a read directory. The
# fastp preprocess functions differ only in the fastp arguments they use, so
# they all call this helper. The HTML and JSON reports are written straight into
# the sample log directory. An earlier version wrote them to the working
# directory under a fixed name, which two runs in the same directory could
# overwrite.
.runFastpStep = function(input.reads = NULL,
                         output.directory = NULL,
                         fastp.path = NULL,
                         fastp.args = "",
                         task = NULL,
                         report.tag = NULL,
                         summary.csv = NULL,
                         merge.reads = FALSE,
                         threads = 1,
                         mem = 8,
                         overwrite = FALSE,
                         quiet = TRUE) {

  #Quick checks
  options(stringsAsFactors = FALSE)
  if (is.null(input.reads) == TRUE){ stop("Please provide raw reads.") }
  if (file.exists(input.reads) == F){ stop("Input reads not found.") }
  if (is.null(output.directory) == TRUE){ stop("Please provide an output directory.") }

  #Checks that fastp is installed before any sample is processed
  fastp.command = .toolCommand("fastp", fastp.path)

  #Sets up the output directory
  if (dir.exists(output.directory) == F){
    dir.create(output.directory, recursive = TRUE)
  } else {
    if (overwrite == TRUE){ .resetDirectory(output.directory) }
  }#end else

  #Creates output directory
  if (dir.exists("logs/sample_logs") == F){ dir.create("logs/sample_logs", recursive = TRUE) }

  #Read in sample data
  input.reads = sub("/+$", "", input.reads)
  reads = list.files(input.reads, recursive = T, full.names = T)
  read.names = .relativePaths(reads, input.reads)
  sample.names = .listSampleNames(input.reads)

  if (length(sample.names) == 0){ return("no samples remain to analyze.") }

  #Creates the summary log
  summary.data = data.frame(Sample = as.character(),
                            Lane = as.character(),
                            Task = as.character(),
                            Program = as.character(),
                            startPairs = as.numeric(),
                            removePairs = as.numeric(),
                            endPairs = as.numeric())

  if (merge.reads == TRUE){ summary.data$mergedReads = as.numeric() }

  for (i in seq_along(sample.names)) {
    #################################################
    ### Part A: prepare for loading and checks
    #################################################
    sample.reads = .matchPrefix(reads, read.names, sample.names[i])

    #Returns a warning if reads are not found
    if (length(sample.reads) == 0 ){
      warning(sample.names[i], " does not have any reads present. Skipping.")
      next
    } #end if statement

    #Check for empty or near-empty input files (sequencing failures)
    file.sizes = file.info(sample.reads)$size
    file.sizes = file.sizes[is.na(file.sizes) == FALSE]
    if (length(file.sizes) == 0 || max(file.sizes) < 1000) {
      largest.size = if (length(file.sizes) == 0) 0 else max(file.sizes)
      failure.msg = paste0("Sample failed: input read files are empty or near-empty",
                           " (max file size: ", largest.size, " bytes).",
                           " This indicates a sequencing or library preparation failure.")
      writeLines(failure.msg, paste0("logs/sample_logs/FAILURE_", sample.names[i], ".txt"))
      warning(sample.names[i], " has empty input read files. Skipping.")
      next
    }

    #Creates new directory
    out.path = paste0(output.directory, "/", sample.names[i])
    report.path = paste0("logs/sample_logs/", sample.names[i])
    if (file.exists(out.path) == FALSE) { dir.create(out.path, recursive = TRUE) }
    if (file.exists(report.path) == FALSE) { dir.create(report.path, recursive = TRUE) }

    lane.prefixes = .stripReadSuffix(sample.reads)

    for (j in seq_along(lane.prefixes)){

      lane.reads = sort(.matchPrefix(reads, reads, lane.prefixes[j]))

      #Returns a warning if reads are not found
      if (length(lane.reads) < 2 ){
        warning(lane.prefixes[j], " does not have a read pair present. Skipping.")
        next
      } #end if statement

      lane.name = basename(lane.prefixes[j])

      #################################################
      ### Part B: Runs fastp
      #################################################
      #sets up output reads
      outreads = c(paste0(out.path, "/", lane.name, "_READ1.fastq.gz"),
                   paste0(out.path, "/", lane.name, "_READ2.fastq.gz"))
      if (merge.reads == TRUE){ outread.m = paste0(out.path, "/", lane.name, "_READ3.fastq.gz") }

      expected.files = outreads
      if (merge.reads == TRUE){ expected.files = c(outreads, outread.m) }

      html.report = paste0(report.path, "/", lane.name, "_", report.tag, ".html")
      json.report = paste0(report.path, "/", lane.name, "_", report.tag, ".json")

      # Skips a lane only when every output file is present and fastp wrote its
      # report, which it does at the end of a successful run. An earlier version
      # skipped a sample as soon as its directory existed, which left an
      # interrupted sample unfinished for good.
      if (overwrite == FALSE && .laneComplete(expected.files, require.size = FALSE) == TRUE &&
          is.null(.fastpReadCounts(json.report)) == FALSE) {
        print(paste0(lane.name, " is already complete. Skipping."))
        next
      }

      merge.arg = ""
      if (merge.reads == TRUE){ merge.arg = paste0(" --merged_out ", shQuote(outread.m)) }

      .runCommand(paste0(fastp.command,
                         " --in1 ", shQuote(lane.reads[1]), " --in2 ", shQuote(lane.reads[2]),
                         " --out1 ", shQuote(outreads[1]), " --out2 ", shQuote(outreads[2]),
                         merge.arg, " ", fastp.args,
                         " --html ", shQuote(html.report), " --json ", shQuote(json.report),
                         " --report_title ", shQuote(lane.name),
                         " --thread ", threads),
                  quiet = quiet, task = paste0("fastp ", task))

      #################################################
      ### Part C: Gathers stats from the fastp report
      #################################################
      # fastp already counts the reads before and after filtering, so the counts
      # come from its JSON report. Counting the fastq files again with gzip
      # doubled the read and write load of every step.
      lane.counts = .fastpReadCounts(json.report)
      if (is.null(lane.counts) == TRUE){
        lane.counts = list(startPairs = .countFastqReads(lane.reads[1]),
                           endPairs = .countFastqReads(outreads[1]),
                           mergedReads = NA_real_)
      }

      temp.remove = data.frame(Sample = sample.names[i],
                               Lane = gsub(".*_", "", lane.name),
                               Task = task,
                               Program = "fastp",
                               startPairs = lane.counts$startPairs,
                               removePairs = lane.counts$startPairs - lane.counts$endPairs,
                               endPairs = lane.counts$endPairs)

      if (merge.reads == TRUE){ temp.remove$mergedReads = lane.counts$mergedReads }

      summary.data = rbind(summary.data, temp.remove)

      print(paste0(lane.name, " completed ", task, "!"))
    }#end lane j loop

    print(paste0(sample.names[i], " Completed ", task, "!"))
  }#end sample i loop

  .appendSummary(summary.data, summary.csv)

  return(invisible(summary.data))
}#end .runFastpStep


# Runs a shell pipeline and stops when any stage fails. A plain pipeline returns
# the status of its last command only, so a failure in an earlier stage, for
# example bwa, would otherwise pass unnoticed.
.runPipeline = function(command = NULL,
                        quiet = TRUE,
                        task = "external command",
                        keep.stdout = FALSE) {

  wrapped = paste0("bash -c ", shQuote(paste0("set -o pipefail; ", command)))
  return(.runCommand(wrapped, quiet = quiet, task = task, keep.stdout = keep.stdout))
}#end .runPipeline


# Writes the awk program that splits a mapped read stream into contaminant pairs
# and clean pairs. The program reads collated SAM records, so both mates of a
# pair arrive together. A pair is a contaminant when either mate aligns at or
# above the identity threshold. Clean pairs go to standard output, the per
# reference counts go to the COUNTS file, and the pair totals go to the STATS
# file. The same test decides both the removal and the count, so the read set
# and the contamination report always agree.
.writeContaminantAwk = function(script.file = NULL) {

  writeLines(c(
    'BEGIN { FS = "\\t"; OFS = "\\t"; name = ""; n = 0; hit = 0; removed = 0; kept = 0 }',
    '/^@/ { print; next }',
    '{',
    '  if ($1 != name) {',
    '    if (n > 0) { emit() }',
    '    name = $1; n = 0; hit = 0',
    '  }',
    '  n = n + 1',
    '  buf[n] = $0',
    '  if (int($2 / 4) % 2 == 0) {',
    '    nm = -1',
    '    for (i = 12; i <= NF; i++) {',
    '      if (substr($i, 1, 5) == "NM:i:") { nm = substr($i, 6) + 0; break }',
    '    }',
    '    qlen = length($10)',
    '    if (nm >= 0 && qlen > 0 && nm / qlen <= MAXMM) { hit = 1; count[$3] = count[$3] + 1 }',
    '  }',
    '}',
    'END {',
    '  if (n > 0) { emit() }',
    '  for (r in count) { line = r "\\t" count[r]; print line > COUNTS }',
    '  line = removed "\\t" kept; print line > STATS',
    '}',
    'function emit(   i) {',
    '  if (hit == 1) { removed = removed + 1 }',
    '  else { kept = kept + 1; for (i = 1; i <= n; i++) { print buf[i] } }',
    '}'
  ), script.file)

  return(invisible(script.file))
}#end .writeContaminantAwk
