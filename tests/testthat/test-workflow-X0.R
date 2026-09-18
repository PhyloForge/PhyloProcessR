# Regression tests for workflow-X0_read-screening.R. The end-to-end tests
# run the real workflow with local reads and are skipped when bwa, fastp, or
# samtools are not on the PATH. Barcode identification is left disabled, so
# MItoTrawlR is not required.

with_x0_working_directory <- function(path, code) {
  old.directory <- setwd(path)
  on.exit(setwd(old.directory), add = TRUE)
  force(code)
}

find_x0_workflow <- function(start = getwd()) {
  directory <- normalizePath(start, mustWork = TRUE)
  repeat {
    candidate <- file.path(directory, "workflows", "workflow-X0_read-screening.R")
    if (file.exists(candidate)) return(candidate)
    parent <- dirname(directory)
    if (identical(parent, directory)) return(character())
    directory <- parent
  }
}

x0_random_seq <- function(n) paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")

x0_revcomp <- function(s) {
  chartr("ACGT", "TGCA", paste(rev(strsplit(s, "")[[1]]), collapse = ""))
}

# Writes one gzipped paired FASTQ lane whose reads derive from one target.
x0_write_lane <- function(prefix, source.seq, npairs = 40) {
  r1 <- gzfile(paste0(prefix, "_READ1.fastq.gz"), "w")
  r2 <- gzfile(paste0(prefix, "_READ2.fastq.gz"), "w")
  left  <- substr(source.seq, 1, 100)
  right <- x0_revcomp(substr(source.seq, 101, 200))
  for (i in seq_len(npairs)) {
    id <- paste0("@read", i)
    writeLines(c(id, left,  "+", strrep("I", nchar(left))),  r1)
    writeLines(c(id, right, "+", strrep("I", nchar(right))), r2)
  }
  close(r1); close(r2)
}

# Builds a self-contained X0 project: three targets, a flat local read set with
# Sample1 (two lanes from t1 and t2) and Sample10 (one lane from t1). The names
# check that Sample1 selection never captures Sample10.
x0_build_project <- function(root, tools, n.targets = 3, delete.reads = FALSE,
                             write.reads = TRUE) {
  read.directory <- file.path(root, "raw-input")
  dir.create(read.directory, recursive = TRUE, showWarnings = FALSE)

  set.seed(42)
  targets <- vapply(seq_len(3), function(i) x0_random_seq(200), character(1))
  names(targets) <- paste0("t", seq_len(3))
  if (n.targets > 3) {
    set.seed(7)
    extra <- vapply(seq_len(n.targets - 3), function(i) x0_random_seq(200), character(1))
    names(extra) <- paste0("t", seq(4, n.targets))
    targets <- c(targets, extra)
  }
  target.fasta <- file.path(root, "targets.fa")
  writeLines(as.vector(rbind(paste0(">", names(targets)), unname(targets))), target.fasta)

  # Read files are written once. Rewriting them would change their modification
  # times, which the cleaning step treats as changed inputs.
  if (write.reads == TRUE) {
    x0_write_lane(file.path(read.directory, "Sample1_L001"),  targets["t1"])
    x0_write_lane(file.path(read.directory, "Sample1_L002"),  targets["t2"])
    x0_write_lane(file.path(read.directory, "Sample10_L001"), targets["t1"])
  }

  config <- c(
    "install.latest.github = FALSE",
    paste0("working.directory = ", deparse(root)),
    "threads = 1",
    "memory = 1",
    "quiet = TRUE",
    "use.dropbox = FALSE",
    "use.sra = FALSE",
    paste0("read.directory = ", deparse(read.directory)),
    "sample.file = 'unused.csv'",
    "read.length = 100",
    "processed.reads = 'processed-reads'",
    "dropbox.token = 'unused'",
    "dropbox.directory = 'unused'",
    "sra.info.file = 'unused.csv'",
    "sra.sample.name.column = NULL",
    "sra.filter.strategy = NULL",
    "sra.max.retries = 3",
    "sra.retry.delay = 10",
    paste0("delete.raw.reads = ", delete.reads),
    paste0("delete.cleaned.reads = ", delete.reads),
    paste0("target.fasta = ", deparse(target.fasta)),
    "run.barcode.scan = FALSE",
    "barcode.fasta = 'unused'",
    "barcode.database.fasta = NULL",
    "barcode.min.iterations = 3",
    "barcode.max.iterations = 10",
    "barcode.min.ref.id = 0.70",
    "barcode.per.max.length = 0.50",
    "barcode.hits.per.sample = 5",
    "fastp.remove.adaptors = TRUE",
    "fastp.remove.duplicate.reads = TRUE",
    "fastp.error.correction = TRUE",
    "fastp.quality.trim.reads = FALSE",
    "fastp.quality.filter = TRUE",
    "fastp.low.complexity.filter = TRUE",
    "fastp.trim.poly.x = TRUE",
    "fastp.min.read.length = 60",
    paste0("fastp.path = ", deparse(unname(tools["fastp"]))),
    paste0("samtools.path = ", deparse(unname(tools["samtools"]))),
    paste0("bwa.path = ", deparse(unname(tools["bwa"]))),
    "bbmap.path = 'unused'",
    "spades.path = 'unused'",
    "cap3.path = 'unused'",
    "blast.path = 'unused'"
  )
  writeLines(config, file.path(root, "workflow-X0_configuration-file.R"))

  workflow.file <- find_x0_workflow()
  file.copy(workflow.file, file.path(root, basename(workflow.file)), overwrite = TRUE)
  list(read.directory = read.directory, target.fasta = target.fasta,
       workflow = file.path(root, basename(workflow.file)))
}

x0_run <- function(root) {
  project.file <- basename(find_x0_workflow())
  with_x0_working_directory(root, {
    source(project.file, local = new.env())
  })
}

x0_tools <- function() {
  tools <- c(bwa = unname(Sys.which("bwa")),
             fastp = unname(Sys.which("fastp")),
             samtools = unname(Sys.which("samtools")))
  tools
}


test_that("X0 keeps processing code inside screenReads", {
  workflow.file <- find_x0_workflow()
  workflow.lines <- readLines(workflow.file, warn = FALSE)

  expect_true(any(grepl("screenReads\\(", workflow.lines)))
  expect_false(any(grepl("for \\(", workflow.lines)))
  expect_false(any(grepl("write.csv\\(", workflow.lines)))
  expect_false(any(grepl("fastqStats\\(", workflow.lines)))
  expect_false(any(grepl("fastpClean\\(", workflow.lines)))
  expect_false(any(grepl("assessCaptureEfficiency\\(", workflow.lines)))
})


test_that("X0 processes each sample separately and reports the target union", {
  tools <- x0_tools()
  skip_if(any(tools == ""), "bwa, fastp, or samtools is not on the PATH")

  root <- tempfile("x2-endtoend-")
  dir.create(root)
  x0_build_project(root, tools)
  x0_run(root)

  final <- read.csv(file.path(root, "logs", "X0_read-screening_FINAL.csv"),
                    stringsAsFactors = FALSE)
  final$Sample <- as.character(final$Sample)

  # Item 1: Sample1 and Sample10 stay separate. Item 7: one row per sample.
  expect_setequal(final$Sample, c("Sample1", "Sample10"))
  expect_true(all(final$Status == "complete"))

  # Item 5/2: read pairs are per sample, summed across the sample's own lanes.
  expect_equal(final$Read_Pairs[final$Sample == "Sample1"], 80)
  expect_equal(final$Read_Pairs[final$Sample == "Sample10"], 40)

  # Item 6: Sample1 hits t1 and t2 across two lanes (union of two of three).
  expect_equal(final$targetsHit[final$Sample == "Sample1"], 2)
  expect_equal(final$targetsHit[final$Sample == "Sample10"], 1)
  expect_equal(final$totalTargets[final$Sample == "Sample1"], 3)
})


test_that("X0 keeps local reads and resumes completed samples", {
  tools <- x0_tools()
  skip_if(any(tools == ""), "bwa, fastp, or samtools is not on the PATH")

  root <- tempfile("x2-resume-")
  dir.create(root)
  project <- x0_build_project(root, tools, delete.reads = FALSE)
  x0_run(root)

  # Item 2: local reads and cleaned reads are retained when the delete flags are FALSE.
  expect_equal(length(list.files(project$read.directory)), 6)
  expect_setequal(list.files(file.path(root, "processed-reads", "cleaned-reads")),
                  c("Sample1", "Sample10"))

  # Completion records were written after the results.
  expect_true(file.exists(file.path(root, "logs", "sample_logs", "Sample1",
                                    "Sample1_X0-complete.csv")))

  # Item 4: a second run recognizes both samples as complete and does not remap.
  messages <- capture.output(x0_run(root))
  expect_equal(sum(grepl("Already complete", messages)), 2)
  expect_false(any(grepl("capture assessment complete", messages)))
})


test_that("X0 reprocesses a sample when the target reference changes", {
  tools <- x0_tools()
  skip_if(any(tools == ""), "bwa, fastp, or samtools is not on the PATH")

  root <- tempfile("x2-target-change-")
  dir.create(root)
  project <- x0_build_project(root, tools, n.targets = 3)
  x0_run(root)

  # Rebuild only the target FASTA with an extra locus, changing its content
  # identity while leaving the read files unchanged.
  x0_build_project(root, tools, n.targets = 4, write.reads = FALSE)
  messages <- capture.output(x0_run(root))

  # Item 4: the changed target prevents a stale skip and remaps both samples.
  expect_false(any(grepl("Already complete", messages)))
  final <- read.csv(file.path(root, "logs", "X0_read-screening_FINAL.csv"),
                    stringsAsFactors = FALSE)
  # Item 6: the recovery percentage uses the new total, not a cached value.
  expect_true(all(final$totalTargets == 4))
})


test_that("X0 fails a sample with an incomplete lane without deleting reads", {
  tools <- x0_tools()
  skip_if(any(tools == ""), "bwa, fastp, or samtools is not on the PATH")

  root <- tempfile("x2-incomplete-")
  dir.create(root)
  project <- x0_build_project(root, tools)
  # Remove one mate of Sample1 lane 2 so the sample has an incomplete pair.
  file.remove(file.path(project$read.directory, "Sample1_L002_READ2.fastq.gz"))

  x0_run(root)

  final <- read.csv(file.path(root, "logs", "X0_read-screening_FINAL.csv"),
                    stringsAsFactors = FALSE)
  final$Sample <- as.character(final$Sample)

  # Item 3: Sample1 is reported as failed, Sample10 still completes.
  expect_true(grepl("failed", final$Status[final$Sample == "Sample1"]))
  expect_equal(final$Status[final$Sample == "Sample10"], "complete")
  # No completion record for the failed sample.
  expect_false(file.exists(file.path(root, "logs", "sample_logs", "Sample1",
                                     "Sample1_X0-complete.csv")))
})


# Dropbox source validation is tested separately from the external service.
