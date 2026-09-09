# Shared internal helpers for workflow 1 preprocessing: program lookup, safe
# shell calls, sample matching, read counting, and summary log writing.


# Length of a sequence without its Ns. Two fragments of one target are joined with
# N padding, so the padding would otherwise count as recovered sequence.
.baseWidth = function(seqs = NULL) {
  if (length(seqs) == 0) return(integer(0))
  return(as.integer(Biostrings::width(seqs) -
                    Biostrings::letterFrequency(seqs, "N")))
}#end .baseWidth


# Builds the shell command for an external program. program.path can be the
# directory that holds the executable, the full path to the executable, or NULL
# to use the system PATH. The function stops when the program is not found.
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
                       keep.stdout = FALSE,
                       stderr.log = NULL) {

  drop.stdout = quiet && keep.stdout == FALSE
  drop.stderr = quiet && is.null(stderr.log)

  if (!is.null(stderr.log)) {
    dir.create(dirname(stderr.log), recursive = TRUE, showWarnings = FALSE)
    command = paste0(command, " 2>> ", shQuote(stderr.log))
  }

  # The default /bin/sh does not report a failed early pipeline stage. Run
  # pipelines with pipefail so a mapper or converter failure reaches the caller.
  if (grepl("|", command, fixed = TRUE)) {
    command = paste0("bash -o pipefail -c ", shQuote(command))
  }

  status = system(command,
                  ignore.stdout = drop.stdout,
                  ignore.stderr = drop.stderr)
  if (status != 0) {
    stop("The ", task, " step failed with exit status ", status, ".\nCommand: ", command)
  }
  return(invisible(status))
}#end .runCommand


# Explicit pipeline entry point. .runCommand also detects pipes for backwards
# compatibility, but callers use this name when pipefail is part of the contract.
.runPipeline = function(command = NULL, quiet = TRUE,
                        task = "external pipeline", stderr.log = NULL) {
  .runCommand(command, quiet = quiet, task = task, keep.stdout = TRUE,
              stderr.log = stderr.log)
}


.validateResources = function(threads = 1, memory = 1, samples = 1,
                              simultaneous.jvms = 1) {
  if (length(threads) != 1 || !is.finite(threads) || threads < 1 ||
      threads != as.integer(threads)) stop("threads must be a positive integer.")
  if (length(memory) != 1 || !is.finite(memory) || memory <= 0)
    stop("memory must be a positive number of GB.")
  workers = min(as.integer(threads), max(1L, as.integer(samples)))
  heap.mb = floor(memory * 1024 / workers / simultaneous.jvms)
  if (heap.mb < 1) stop("The requested memory budget provides no usable JVM heap.")
  list(workers = workers, heap.mb = heap.mb)
}


.ensureDirectory = function(path, label = "directory") {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop("Could not create ", label, ": ", path)
  invisible(path)
}


.stageMarker = function(directory, stage) file.path(directory, paste0(".", stage, ".complete"))

.stageComplete = function(directory, stage, outputs) {
  file.exists(.stageMarker(directory, stage)) && length(outputs) > 0 &&
    all(file.exists(outputs)) && all(file.info(outputs)$size > 0)
}

.markStageComplete = function(directory, stage, details = character()) {
  writeLines(c("complete=true", details), .stageMarker(directory, stage))
}

.invalidateStage = function(directory, stage) {
  marker = .stageMarker(directory, stage)
  if (file.exists(marker)) file.remove(marker)
  invisible(NULL)
}

.collectWorkers = function(results, sample.names, stage) {
  failed = vapply(seq_along(sample.names), function(i) {
    result = if (i <= length(results)) results[[i]] else NULL
    is.null(result) || inherits(result, "try-error") || !isTRUE(result$success)
  }, logical(1))
  if (any(failed)) {
    messages = vapply(which(failed), function(i) {
      result = if (i <= length(results)) results[[i]] else NULL
      if (is.list(result) && !is.null(result$message)) result$message else "worker did not return success"
    }, character(1))
    stop(stage, " failed for sample(s): ",
         paste(paste0(sample.names[failed], " (", messages, ")"), collapse = ", "))
  }
  invisible(results)
}

.gatkCommand = function(gatk, temp.directory, heap.mb) {
  java.options = paste0("-Djava.io.tmpdir=", temp.directory, " -Xmx", heap.mb, "m")
  paste(gatk, "--java-options", shQuote(java.options))
}

.laneDirectories = function(sample.directory) {
  dirs = list.dirs(sample.directory, recursive = FALSE, full.names = TRUE)
  dirs[grepl("^Lane_[0-9]+$", basename(dirs))]
}

.selectedSampleBam = function(mapping.directory, sample, use.base.recalibration = FALSE) {
  if (use.base.recalibration) {
    bqsr = file.path(mapping.directory, sample, "Lane_Merge", "bqsr-mapped-all.bam")
    if (file.exists(bqsr)) return(bqsr)
    bqsr = list.files(file.path(mapping.directory, sample), "^bqsr-mapped-all\\.bam$",
                      recursive = TRUE, full.names = TRUE)
    if (length(bqsr) == 1) return(bqsr)
  }
  merged = file.path(mapping.directory, sample, "Lane_Merge", "final-mapped-all.bam")
  if (file.exists(merged)) return(merged)
  lanes = .laneDirectories(file.path(mapping.directory, sample))
  bams = file.path(lanes, "final-mapped-all.bam")
  bams = bams[file.exists(bams)]
  if (length(bams) == 1) return(bams)
  stop("Could not select exactly one calling BAM for sample ", sample, ".")
}


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

  base.dir = paste0(sub("/+$", "", base.dir), "/")
  relative.paths = file.paths
  inside.base = startsWith(file.paths, base.dir)
  relative.paths[inside.base] = substring(file.paths[inside.base], nchar(base.dir) + 1)
  return(relative.paths)
}#end .relativePaths


# Lists only supported FASTQ files. Reports, sentinels, and other files in a
# read directory must not take part in sample or lane discovery.
.listFastqFiles = function(read.directory = NULL,
                           recursive = TRUE) {

  reads = list.files(read.directory, recursive = recursive, full.names = TRUE)
  reads = reads[grepl("\\.(fastq|fq)(\\.gz)?$", reads, ignore.case = TRUE)]
  return(reads)
}#end .listFastqFiles


# Selects the files whose match string starts with a prefix and a recognized
# separator or bare READ label. Fixed-string comparisons keep regular
# expression characters safe and stop Sample1 from matching Sample10.
.matchPrefix = function(file.paths = NULL,
                        match.strings = NULL,
                        prefix = NULL) {

  keep = startsWith(match.strings, paste0(prefix, "_")) |
         startsWith(match.strings, paste0(prefix, "-")) |
         startsWith(match.strings, paste0(prefix, ".")) |
         startsWith(match.strings, paste0(prefix, "/")) |
         startsWith(toupper(match.strings), paste0(toupper(prefix), "READ"))

  return(file.paths[keep])
}#end .matchPrefix


# Removes the read and mate suffix from a read file path. The result is the lane
# prefix, for example Sample_L001. Multiple lanes give multiple prefixes.
.stripReadSuffix = function(file.paths = NULL) {

  named.suffix = paste0("(_R[123]|-R[123]|_READ[123]|-READ[123]|READ[123])",
                        ".*\\.(fastq|fq)(\\.gz)?$")
  lane.prefixes = sub(named.suffix, "", file.paths, ignore.case = TRUE)
  numeric.suffix = paste0("(_[123]|-[123])([_.-].*)?",
                          "\\.(fastq|fq)(\\.gz)?$")
  lane.prefixes = sub(numeric.suffix, "", lane.prefixes, ignore.case = TRUE)
  singleton.suffix = paste0("(_singleton|-singleton|READ-singleton|READ_singleton|",
                            "_READ-singleton|-READ_singleton|-READ-singleton|",
                            "_READ_singleton).*\\.(fastq|fq)(\\.gz)?$")
  lane.prefixes = sub(singleton.suffix, "", lane.prefixes, ignore.case = TRUE)
  return(unique(lane.prefixes))
}#end .stripReadSuffix


# Puts read files in first-mate, second-mate, optional third-read order. Every
# naming form accepted by the older preprocessing functions remains supported,
# including extra instrument text around the mate label.
.orderReadFiles = function(read.files = NULL,
                           allow.third = TRUE) {

  read.files = read.files[grepl("\\.(fastq|fq)(\\.gz)?$", read.files, ignore.case = TRUE)]
  base.names = basename(read.files)
  mate.pattern = function(number) {
    paste0("((_R", number, "|-R", number, "|_READ", number, "|-READ", number,
           "|READ", number, ").*|(_", number, "|-", number, ")([_.-].*)?)",
           "\\.(fastq|fq)(\\.gz)?$")
  }
  mate.matches = sapply(1:3, function(number) {
    grepl(mate.pattern(number), base.names, ignore.case = TRUE)
  })
  if (length(read.files) == 1) { mate.matches = matrix(mate.matches, nrow = 1) }
  singleton.pattern = paste0("(_singleton|-singleton|READ-singleton|READ_singleton|",
                             "_READ-singleton|-READ_singleton|-READ-singleton|",
                             "_READ_singleton)")
  singleton = grepl(singleton.pattern, base.names, ignore.case = TRUE)
  mate.matches[singleton, ] = FALSE
  mate.matches[singleton, 3] = TRUE
  mate.number = max.col(mate.matches, ties.method = "first")
  mate.number[rowSums(mate.matches) != 1] = NA_integer_

  valid.length = length(read.files) == 2 || (allow.third == TRUE && length(read.files) == 3)
  if (sum(mate.number == 1, na.rm = TRUE) != 1 ||
      sum(mate.number == 2, na.rm = TRUE) != 1 ||
      sum(mate.number == 3, na.rm = TRUE) > as.integer(allow.third) ||
      any(is.na(mate.number)) || valid.length == FALSE) {
    stop("Expected exactly one READ1 and one READ2 FASTQ file, but found: ",
         paste(base.names, collapse = ", "))
  }

  return(read.files[order(mate.number)])
}#end .orderReadFiles


.orderReadPair = function(read.files = NULL) {

  return(.orderReadFiles(read.files, allow.third = FALSE))
}#end .orderReadPair


# Lists the sample names in a read directory. Samples are sub-directories when
# they exist, otherwise the file names give the sample names.
.listSampleNames = function(read.directory = NULL) {

  sample.names = list.dirs(read.directory, recursive = FALSE, full.names = FALSE)

  if (length(sample.names) == 0) {
    read.files = .listFastqFiles(read.directory, recursive = FALSE)
    sample.names = basename(.stripReadSuffix(read.files))
    sample.names = unique(sub("_L[0-9]+$", "", sample.names))
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

  report = tryCatch(jsonlite::fromJSON(json.file), error = function(e) NULL)
  if (is.null(report) == TRUE ||
      is.null(report$summary$before_filtering$total_reads) == TRUE ||
      is.null(report$summary$after_filtering$total_reads) == TRUE) { return(NULL) }

  before.reads = suppressWarnings(as.numeric(report$summary$before_filtering$total_reads))
  after.reads = suppressWarnings(as.numeric(report$summary$after_filtering$total_reads))
  if (length(before.reads) != 1 || length(after.reads) != 1 ||
      is.finite(before.reads) == FALSE || is.finite(after.reads) == FALSE ||
      before.reads < 0 || after.reads < 0) { return(NULL) }

  counts = list(startPairs = before.reads / 2,
                endPairs = after.reads / 2,
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
                          out.csv = NULL,
                          replace.samples = NULL) {

  if (is.null(summary.data) == TRUE) { return(invisible(NULL)) }
  if (is.null(replace.samples) == TRUE && nrow(summary.data) > 0) {
    replace.samples = unique(summary.data$Sample)
  }

  if (file.exists(out.csv) == TRUE) {
    existing = read.csv(out.csv, stringsAsFactors = FALSE)
    if (identical(sort(names(existing)), sort(names(summary.data))) == TRUE) {
      existing = existing[!existing$Sample %in% replace.samples, , drop = FALSE]
      summary.data = rbind(existing[, names(summary.data), drop = FALSE], summary.data)
    }
  }

  dir.create(dirname(out.csv), recursive = TRUE, showWarnings = FALSE)
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

  output.path = normalizePath(output.directory, mustWork = FALSE)
  protected.paths = c(normalizePath("/", mustWork = TRUE),
                      normalizePath(path.expand("~"), mustWork = TRUE),
                      normalizePath(getwd(), mustWork = TRUE))
  if (nchar(output.directory) == 0 || output.path %in% protected.paths) {
    stop("Refusing to reset a protected or empty output directory.")
  }
  if (dir.exists(output.directory) == TRUE) { unlink(output.directory, recursive = TRUE) }
  dir.create(output.directory, recursive = TRUE, showWarnings = FALSE)
  return(invisible(output.directory))
}#end .resetDirectory


# Rejects an output directory that is the input directory or contains it. It
# also rejects an output below the input because recursive read discovery would
# otherwise find prior outputs as new inputs.
.normalizedPath = function(path = NULL) {

  missing.parts = character()
  existing.path = path
  while (file.exists(existing.path) == FALSE) {
    missing.parts = c(basename(existing.path), missing.parts)
    parent.path = dirname(existing.path)
    if (identical(parent.path, existing.path)) { break }
    existing.path = parent.path
  }
  normalized.path = normalizePath(existing.path, mustWork = TRUE)
  if (length(missing.parts) > 0) {
    normalized.path = do.call(file.path, as.list(c(normalized.path, missing.parts)))
  }
  return(normalized.path)
}#end .normalizedPath


.checkDirectoryOverlap = function(input.directory = NULL,
                                  output.directory = NULL) {

  input.path = .normalizedPath(input.directory)
  output.path = .normalizedPath(output.directory)

  nested = identical(input.path, output.path) ||
           startsWith(paste0(input.path, "/"), paste0(output.path, "/")) ||
           startsWith(paste0(output.path, "/"), paste0(input.path, "/"))
  if (nested == TRUE) {
    stop("The input and output directories must be separate and cannot contain one another.")
  }
  return(invisible(TRUE))
}#end .checkDirectoryOverlap


.checkFileOutsideOutput = function(input.file = NULL,
                                   output.directory = NULL) {

  input.path = .normalizedPath(input.file)
  output.path = .normalizedPath(output.directory)
  if (identical(input.path, output.path) ||
      startsWith(input.path, paste0(output.path, "/"))) {
    stop("The output directory cannot contain an input file that it may overwrite.")
  }
  return(invisible(TRUE))
}#end .checkFileOutsideOutput


# Checks the File/Sample tables used by organization and download functions.
# Sample names become directory names and therefore cannot contain path or
# control components. Distinct names that clean to the same value are rejected.
.validateRenameTable = function(sample.data = NULL,
                                sanitize.samples = FALSE) {

  required.columns = c("File", "Sample")
  if (all(required.columns %in% names(sample.data)) == FALSE) {
    stop("The rename table must contain File and Sample columns.")
  }

  file.names = trimws(as.character(sample.data$File))
  sample.names = trimws(as.character(sample.data$Sample))
  if (any(is.na(file.names)) || any(nchar(file.names) == 0)) {
    stop("The File column cannot contain missing or blank values.")
  }
  if (any(is.na(sample.names)) || any(nchar(sample.names) == 0)) {
    stop("The Sample column cannot contain missing or blank values.")
  }
  if (any(grepl("[/\\\\]", sample.names)) || any(sample.names %in% c(".", ".."))) {
    stop("Sample names cannot contain directory separators or path components.")
  }

  clean.names = .sanitizeName(sample.names)
  name.map = unique(data.frame(original = sample.names, clean = clean.names,
                               stringsAsFactors = FALSE))
  if (sanitize.samples == TRUE && any(duplicated(name.map$clean))) {
    stop("Distinct sample names become identical after filename cleaning.")
  }

  sample.data$File = file.names
  sample.data$Sample = if (sanitize.samples == TRUE) clean.names else sample.names
  return(sample.data)
}#end .validateRenameTable


# Records the lightweight identity of large read files. Paths, sizes, and
# modification times catch ordinary replacement or reordering without hashing
# every FASTQ file on every resume.
.laneMetadata = function(input.files = NULL,
                         parameters = list()) {

  input.files = normalizePath(input.files, mustWork = TRUE)
  file.data = file.info(input.files)
  metadata = c()
  for (i in seq_along(input.files)) {
    metadata[paste0("input.", i, ".path")] = input.files[i]
    metadata[paste0("input.", i, ".size")] = format(file.data$size[i], scientific = FALSE)
    metadata[paste0("input.", i, ".modified")] = format(as.numeric(file.data$mtime[i]),
                                                          scientific = FALSE)
  }
  if (length(parameters) > 0) {
    parameter.values = vapply(parameters, function(value) paste(value, collapse = ";"), character(1))
    names(parameter.values) = paste0("parameter.", names(parameters))
    metadata = c(metadata, parameter.values)
  }
  return(metadata)
}#end .laneMetadata


.metadataMatches = function(metadata.file = NULL,
                            metadata = NULL) {

  saved = .readLaneMetadata(metadata.file)
  if (is.null(saved) == TRUE) { return(FALSE) }
  current = data.frame(Field = names(metadata), Value = unname(as.character(metadata)),
                       stringsAsFactors = FALSE)
  return(identical(saved, current))
}#end .metadataMatches


.readLaneMetadata = function(metadata.file = NULL) {

  if (file.exists(metadata.file) == FALSE) { return(NULL) }
  saved = tryCatch(read.csv(metadata.file, stringsAsFactors = FALSE,
                            colClasses = "character"),
                   error = function(e) NULL)
  if (is.null(saved) == TRUE || !identical(names(saved), c("Field", "Value"))) {
    return(NULL)
  }
  return(saved)
}#end .readLaneMetadata


# A readable metadata file with different values means the completed output
# belongs to other inputs or settings. Missing or malformed metadata is treated
# as an interrupted write and the lane is rebuilt.
.metadataConflicts = function(metadata.file = NULL,
                              metadata = NULL) {

  if (is.null(.readLaneMetadata(metadata.file)) == TRUE) { return(FALSE) }
  return(.metadataMatches(metadata.file, metadata) == FALSE)
}#end .metadataConflicts


.writeLaneMetadata = function(metadata = NULL,
                              metadata.file = NULL) {

  dir.create(dirname(metadata.file), recursive = TRUE, showWarnings = FALSE)
  temp.file = tempfile(pattern = paste0(basename(metadata.file), "-"),
                       tmpdir = dirname(metadata.file))
  write.csv(data.frame(Field = names(metadata), Value = unname(as.character(metadata)),
                       stringsAsFactors = FALSE), temp.file, row.names = FALSE)
  if (file.rename(temp.file, metadata.file) == FALSE) {
    unlink(temp.file)
    stop("Could not publish completion metadata at ", metadata.file, ".")
  }
  return(invisible(metadata.file))
}#end .writeLaneMetadata


# Moves completed temporary files into place. Temporary files are created in
# the destination directory so rename is atomic on the local file system.
.publishFiles = function(temp.files = NULL,
                         output.files = NULL) {

  if (length(temp.files) != length(output.files)) {
    stop("Temporary and output file lists have different lengths.")
  }
  for (i in seq_along(output.files)) {
    if (file.exists(temp.files[i]) == FALSE) {
      stop("Expected temporary output was not created: ", temp.files[i])
    }
    if (file.exists(output.files[i]) == TRUE) { unlink(output.files[i]) }
    if (file.rename(temp.files[i], output.files[i]) == FALSE) {
      stop("Could not publish output file: ", output.files[i])
    }
  }
  return(invisible(output.files))
}#end .publishFiles


# Reports whether every output file of a lane is present. An interrupted lane is
# processed again instead of being skipped for good.
# require.size adds the test that each file holds data. Set it to FALSE where an
# empty output file is a valid result, for example the merged read file of a
# sample whose read pairs do not overlap.
.laneComplete = function(output.files = NULL,
                         require.size = TRUE,
                         metadata.file = NULL,
                         metadata = NULL) {

  if (length(output.files) == 0) { return(FALSE) }
  if (all(file.exists(output.files)) == FALSE) { return(FALSE) }
  if (require.size == TRUE && all(file.info(output.files)$size > 0) == FALSE) { return(FALSE) }
  if (is.null(metadata.file) == FALSE &&
      .metadataMatches(metadata.file, metadata) == FALSE) { return(FALSE) }
  return(TRUE)
}#end .laneComplete


# Runs one fastp cleaning step over every sample lane in a read directory. The
# fastp preprocess functions differ only in the fastp arguments they use, so
# they all call this helper. The HTML and JSON reports are written straight into
# the sample log directory so simultaneous runs cannot overwrite one another.
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
  if (length(threads) != 1 || is.numeric(threads) == FALSE ||
      is.finite(threads) == FALSE || threads < 1) {
    stop("threads must be one positive number.")
  }
  if (length(overwrite) != 1 || is.logical(overwrite) == FALSE || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE.")
  }
  .checkDirectoryOverlap(input.reads, output.directory)

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
  reads = .listFastqFiles(input.reads)
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

    # A zero-byte file is incomplete. Small compressed files can still contain
    # valid reads, including a valid empty FASTQ from an earlier filter.
    failure.file = paste0("logs/sample_logs/FAILURE_", sample.names[i], ".txt")
    file.sizes = file.info(sample.reads)$size
    if (any(is.na(file.sizes)) || any(file.sizes == 0)) {
      writeLines("Sample skipped because at least one input read file is missing or zero bytes.",
                 failure.file)
      warning(sample.names[i], " has a missing or zero-byte input read file. Skipping.")
      next
    }
    if (file.exists(failure.file) == TRUE) { unlink(failure.file) }

    #Creates new directory
    out.path = paste0(output.directory, "/", sample.names[i])
    report.path = paste0("logs/sample_logs/", sample.names[i])
    if (file.exists(out.path) == FALSE) { dir.create(out.path, recursive = TRUE) }
    if (file.exists(report.path) == FALSE) { dir.create(report.path, recursive = TRUE) }

    lane.prefixes = .stripReadSuffix(sample.reads)

    for (j in seq_along(lane.prefixes)){

      lane.reads = .matchPrefix(reads, reads, lane.prefixes[j])

      #Returns a warning if reads are not found
      if (length(lane.reads) != 2 ){
        warning(lane.prefixes[j], " does not have a read pair present. Skipping.")
        next
      } #end if statement
      lane.reads = tryCatch(.orderReadPair(lane.reads), error = function(e) {
        warning(conditionMessage(e))
        return(NULL)
      })
      if (is.null(lane.reads) == TRUE) { next }

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
      metadata.file = paste0(report.path, "/", lane.name, "_", report.tag, "-metadata.csv")
      metadata = .laneMetadata(lane.reads,
                               list(task = task, fastp.arguments = fastp.args,
                                    merge.reads = merge.reads))

      # The metadata file is written last, after outputs and reports. It is the
      # completion marker and records the input identity and relevant settings.
      if (overwrite == FALSE && .laneComplete(expected.files, require.size = FALSE) == TRUE &&
          file.exists(html.report) == TRUE &&
          is.null(.fastpReadCounts(json.report)) == FALSE &&
          .metadataMatches(metadata.file, metadata) == TRUE) {
        lane.counts = .fastpReadCounts(json.report)
        temp.remove = data.frame(Sample = sample.names[i],
                                 Lane = gsub(".*_", "", lane.name),
                                 Task = task,
                                 Program = "fastp",
                                 startPairs = lane.counts$startPairs,
                                 removePairs = lane.counts$startPairs - lane.counts$endPairs,
                                 endPairs = lane.counts$endPairs)
        if (merge.reads == TRUE){ temp.remove$mergedReads = lane.counts$mergedReads }
        summary.data = rbind(summary.data, temp.remove)
        print(paste0(lane.name, " is already complete. Skipping."))
        next
      }
      if (overwrite == FALSE && .metadataConflicts(metadata.file, metadata) == TRUE) {
        stop(lane.name, " was completed with different inputs or settings. ",
             "Use overwrite = TRUE to replace it.")
      }

      unlink(c(expected.files, html.report, json.report, metadata.file))
      temp.reads = vapply(expected.files, function(output.file) {
        tempfile(pattern = paste0(basename(output.file), "-"),
                 tmpdir = dirname(output.file), fileext = ".fastq.gz")
      }, character(1))
      temp.html = tempfile(pattern = paste0(basename(html.report), "-"),
                           tmpdir = report.path, fileext = ".html")
      temp.json = tempfile(pattern = paste0(basename(json.report), "-"),
                           tmpdir = report.path, fileext = ".json")
      on.exit(unlink(c(temp.reads, temp.html, temp.json)), add = TRUE)

      merge.arg = ""
      if (merge.reads == TRUE){ merge.arg = paste0(" --merged_out ", shQuote(temp.reads[3])) }

      .runCommand(paste0(fastp.command,
                         " --in1 ", shQuote(lane.reads[1]), " --in2 ", shQuote(lane.reads[2]),
                         " --out1 ", shQuote(temp.reads[1]), " --out2 ", shQuote(temp.reads[2]),
                         merge.arg, " ", fastp.args,
                         " --html ", shQuote(temp.html), " --json ", shQuote(temp.json),
                         " --report_title ", shQuote(lane.name),
                         " --thread ", threads),
                  quiet = quiet, task = paste0("fastp ", task))

      #################################################
      ### Part C: Gathers stats from the fastp report
      #################################################
      # fastp already counts the reads before and after filtering, so the counts
      # come from its JSON report. Counting the fastq files again with gzip
      # doubled the read and write load of every step.
      lane.counts = .fastpReadCounts(temp.json)
      if (is.null(lane.counts) == TRUE){
        stop("fastp did not write a complete JSON report for ", lane.name, ".")
      }

      .publishFiles(c(temp.reads, temp.html, temp.json),
                    c(expected.files, html.report, json.report))
      .writeLaneMetadata(metadata, metadata.file)

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
    '    alen = alignedLength($6)',
    '    if (nm >= 0 && alen > 0 && nm / alen <= MAXMM) { hit = 1; count[$3] = count[$3] + 1 }',
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
    '}',
    'function alignedLength(cigar,   rest, token, amount, operation, total) {',
    '  rest = cigar; total = 0',
    '  while (match(rest, /^[0-9]+[MIDNSHP=X]/)) {',
    '    token = substr(rest, RSTART, RLENGTH)',
    '    amount = substr(token, 1, length(token) - 1) + 0',
    '    operation = substr(token, length(token), 1)',
    '    if (operation == "M" || operation == "I" || operation == "D" ||',
    '        operation == "=" || operation == "X") { total = total + amount }',
    '    rest = substr(rest, RLENGTH + 1)',
    '  }',
    '  return(total)',
    '}'
  ), script.file)

  return(invisible(script.file))
}#end .writeContaminantAwk
