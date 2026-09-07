with_preprocess_test_directory = function(path, code) {
  old.directory = setwd(path)
  on.exit(setwd(old.directory), add = TRUE)
  force(code)
}


write_test_fastq = function(path) {
  writeLines(c("@read", "ACGT", "+", "IIII"), path)
}


test_that("read matching preserves supported names and extensions", {
  first.reads = c(
    "Sample_L001_1.fastq.gz", "Sample_L001-1.fq",
    "Sample_L001_R1_001.fastq.gz", "Sample_L001-R1-extra.machine.fastq",
    "Sample_L001_READ1_random.fq.gz", "Sample_L001-READ1.random.fastq.gz",
    "Sample_L001READ1instrument.fq"
  )
  second.reads = c(
    "Sample_L001_2.fastq.gz", "Sample_L001-2.fq",
    "Sample_L001_R2_001.fastq.gz", "Sample_L001-R2-extra.machine.fastq",
    "Sample_L001_READ2_random.fq.gz", "Sample_L001-READ2.random.fastq.gz",
    "Sample_L001READ2instrument.fq"
  )

  for (i in seq_along(first.reads)) {
    expect_identical(
      basename(PhyloProcessR:::.stripReadSuffix(c(first.reads[i], second.reads[i]))),
      "Sample_L001"
    )
    expect_identical(
      PhyloProcessR:::.orderReadPair(c(second.reads[i], first.reads[i])),
      c(first.reads[i], second.reads[i])
    )
    lane.prefix = PhyloProcessR:::.stripReadSuffix(c(first.reads[i], second.reads[i]))[1]
    expect_identical(
      PhyloProcessR:::.matchPrefix(c(first.reads[i], second.reads[i]),
                                   c(first.reads[i], second.reads[i]), lane.prefix),
      c(first.reads[i], second.reads[i])
    )
  }

  third.read = "Sample_L001_READ3_extra.fastq.gz"
  ordered = PhyloProcessR:::.orderReadFiles(
    c(third.read, "Sample_L001_READ2_x.fastq.gz", "Sample_L001_READ1_x.fastq.gz")
  )
  expect_true(all(mapply(grepl, c("READ1", "READ2", "READ3"), ordered)))
})


test_that("relative read paths do not treat directories as regular expressions", {
  base.directory = "/tmp/reads[old]+set"
  paths = file.path(base.directory, c("Sample1_R1.fastq.gz", "Sample1_R2.fastq.gz"))
  expect_identical(
    PhyloProcessR:::.relativePaths(paths, base.directory),
    c("Sample1_R1.fastq.gz", "Sample1_R2.fastq.gz")
  )
})


test_that("lane metadata detects changed inputs and malformed reports", {
  root = tempfile("preprocess-metadata-")
  dir.create(root)
  input.file = file.path(root, "read.fastq")
  metadata.file = file.path(root, "logs", "lane.csv")
  write_test_fastq(input.file)

  metadata = PhyloProcessR:::.laneMetadata(input.file, list(setting = "one"))
  PhyloProcessR:::.writeLaneMetadata(metadata, metadata.file)
  expect_true(PhyloProcessR:::.metadataMatches(metadata.file, metadata))
  expect_false(PhyloProcessR:::.metadataMatches(
    metadata.file,
    PhyloProcessR:::.laneMetadata(input.file, list(setting = "two"))
  ))

  bad.json = file.path(root, "bad.json")
  writeLines("{not complete", bad.json)
  expect_null(PhyloProcessR:::.fastpReadCounts(bad.json))
})


test_that("fastp resume rebuilds complete summaries and rejects changed settings", {
  root = tempfile("preprocess-fastp-")
  dir.create(root)
  input.directory = file.path(root, "input", "Sample")
  output.directory = file.path(root, "output")
  dir.create(input.directory, recursive = TRUE)
  for (lane in c("L001", "L002")) {
    write_test_fastq(file.path(input.directory, paste0("Sample_", lane, "_READ1.fastq")))
    write_test_fastq(file.path(input.directory, paste0("Sample_", lane, "_READ2.fastq")))
  }

  fake.fastp = file.path(root, "fastp")
  writeLines(c(
    "#!/bin/sh",
    "while [ $# -gt 0 ]; do",
    "  case \"$1\" in",
    "    --out1) shift; out1=\"$1\" ;;",
    "    --out2) shift; out2=\"$1\" ;;",
    "    --html) shift; html=\"$1\" ;;",
    "    --json) shift; json=\"$1\" ;;",
    "  esac",
    "  shift",
    "done",
    "printf '@read\\nACGT\\n+\\nIIII\\n' > \"$out1\"",
    "printf '@read\\nACGT\\n+\\nIIII\\n' > \"$out2\"",
    "printf '<html></html>\\n' > \"$html\"",
    "printf '{\"summary\":{\"before_filtering\":{\"total_reads\":2},\"after_filtering\":{\"total_reads\":2}}}\\n' > \"$json\""
  ), fake.fastp)
  Sys.chmod(fake.fastp, mode = "0755")
  summary.file = file.path(root, "logs", "fastp-summary.csv")

  run.fastp = function(arguments = "--setting-one") {
    PhyloProcessR:::.runFastpStep(
      input.reads = dirname(input.directory), output.directory = output.directory,
      fastp.path = fake.fastp, fastp.args = arguments, task = "cleaning",
      report.tag = "cleaning", summary.csv = summary.file,
      threads = 1, overwrite = FALSE, quiet = TRUE
    )
  }

  with_preprocess_test_directory(root, run.fastp())
  expect_equal(nrow(read.csv(summary.file)), 2)
  unlink(summary.file)
  with_preprocess_test_directory(root, run.fastp())
  expect_equal(nrow(read.csv(summary.file)), 2)
  expect_error(with_preprocess_test_directory(root, run.fastp("--setting-two")),
               "different inputs or settings")
})


test_that("deduplication uses the fastp default accuracy", {
  fastp.calls = list()
  local_mocked_bindings(
    .runFastpStep = function(...) {
      fastp.calls[[length(fastp.calls) + 1]] <<- list(...)
    },
    .package = "PhyloProcessR"
  )

  fastpClean(input.reads = "input", output.directory = "output")
  removeDuplicateReads(input.reads = "input", output.directory = "output")

  expect_length(fastp.calls, 2)
  expect_identical(fastp.calls[[1]]$summary.csv, "logs/fastp_summary.csv")
  for (fastp.call in fastp.calls) {
    expect_match(fastp.call$fastp.args, "(^| )--dedup($| )")
    expect_false(grepl("--dup_calc_accuracy", fastp.call$fastp.args,
                       fixed = TRUE))
  }
})


test_that("the environment pins the tested fastp version", {
  environment.file = test_path("..", "..", "setup-files", "environment.yml")
  expect_true(file.exists(environment.file))
  environment.lines = readLines(environment.file)
  expect_identical(grep("^  - fastp=", environment.lines, value = TRUE),
                   "  - fastp=1.3.6")
})


test_that("empty summaries remove obsolete sample rows", {
  root = tempfile("preprocess-summary-")
  dir.create(root)
  summary.file = file.path(root, "contaminants.csv")
  old.summary = data.frame(Sample = c("A", "B"), Lane = c("L001", "L001"),
                           Contaminant = c("human", "mouse"), Reads = c(2, 3))
  write.csv(old.summary, summary.file, row.names = FALSE)
  empty.summary = old.summary[0, ]

  PhyloProcessR:::.appendSummary(empty.summary, summary.file, replace.samples = "A")
  updated = read.csv(summary.file, stringsAsFactors = FALSE)
  expect_identical(updated$Sample, "B")
})


test_that("organizeReads compresses plain FASTQ files and keeps their sources", {
  root = tempfile("preprocess-organize-")
  dir.create(root)
  input.directory = file.path(root, "input")
  output.directory = file.path(root, "output")
  dir.create(input.directory)
  read1 = file.path(input.directory, "source.random_R1_001.fastq")
  read2 = file.path(input.directory, "source.random_R2_001.fq")
  write_test_fastq(read1)
  write_test_fastq(read2)
  rename.file = file.path(root, "rename.csv")
  write.csv(data.frame(File = "source", Sample = "Sample A"),
            rename.file, row.names = FALSE)

  with_preprocess_test_directory(root, {
    organizeReads(input.directory, output.directory, rename.file,
                  link.reads = TRUE, overwrite = FALSE)
  })

  outputs = list.files(output.directory, recursive = TRUE, full.names = TRUE)
  expect_length(outputs, 2)
  expect_true(all(grepl("[.]fastq[.]gz$", outputs)))
  expect_equal(PhyloProcessR:::.countFastqReads(outputs[1]), 1)
  expect_true(file.exists(read1))
  expect_true(file.exists(read2))
})


test_that("rename tables reject unsafe names and output overlap", {
  expect_error(
    PhyloProcessR:::.validateRenameTable(data.frame(File = "a", Sample = "../sample")),
    "path components"
  )
  expect_error(
    PhyloProcessR:::.validateRenameTable(
      data.frame(File = c("a", "b"), Sample = c("A B", "A_B")),
      sanitize.samples = TRUE
    ),
    "identical"
  )

  root = tempfile("preprocess-overlap-")
  input.directory = file.path(root, "input")
  dir.create(input.directory, recursive = TRUE)
  expect_error(
    PhyloProcessR:::.checkDirectoryOverlap(input.directory,
                                           file.path(input.directory, "output")),
    "cannot contain"
  )
})


test_that("SRA filters fail clearly when metadata columns are absent", {
  root = tempfile("preprocess-sra-")
  dir.create(root)
  sra.file = file.path(root, "runs.csv")
  write.csv(data.frame(Run = "SRR000001"), sra.file, row.names = FALSE)

  expect_error(
    sraDownload(sra.info.file = sra.file,
                output.directory = file.path(root, "reads"),
                filter.library.layout = "PAIRED"),
    "LibraryLayout is absent"
  )
})


test_that("workflow 1 rejects simultaneous download sources before setup", {
  workflow.file = test_path("..", "..", "workflows", "workflow-1_preprocess.R")
  root = tempfile("preprocess-workflow-")
  dir.create(root)
  file.copy(workflow.file, file.path(root, "workflow-1_preprocess.R"))
  writeLines(c("dropbox.download = TRUE", "sra.download = TRUE"),
             file.path(root, "workflow-1_configuration-file.R"))

  expect_error(
    with_preprocess_test_directory(
      root,
      sys.source("workflow-1_preprocess.R", envir = new.env(parent = globalenv()))
    ),
    "one download source"
  )
  expect_false(dir.exists(file.path(root, "processed-reads")))
})


test_that("the contaminant identity denominator excludes soft clips", {
  root = tempfile("preprocess-identity-")
  dir.create(root)
  script.file = file.path(root, "filter.awk")
  sam.file = file.path(root, "reads.sam")
  counts.file = file.path(root, "counts.txt")
  stats.file = file.path(root, "stats.txt")
  clean.file = file.path(root, "clean.sam")
  PhyloProcessR:::.writeContaminantAwk(script.file)
  soft.sequence = paste(rep("A", 150), collapse = "")
  aligned.sequence = paste(rep("A", 50), collapse = "")
  writeLines(c(
    "@HD\tVN:1.6",
    paste("soft", 99, "ref", 1, 60, "100S50M", "=", 1, 0,
          soft.sequence, "*", "NM:i:10", sep = "\t"),
    paste("soft", 147, "ref", 1, 60, "100S50M", "=", 1, 0,
          soft.sequence, "*", "NM:i:10", sep = "\t"),
    paste("full", 99, "ref", 1, 60, "50M", "=", 1, 0,
          aligned.sequence, "*", "NM:i:5", sep = "\t"),
    paste("full", 147, "ref", 1, 60, "50M", "=", 1, 0,
          aligned.sequence, "*", "NM:i:5", sep = "\t")
  ), sam.file)

  command = paste0("awk -v MAXMM=0.1 -v COUNTS=", shQuote(counts.file),
                   " -v STATS=", shQuote(stats.file), " -f ", shQuote(script.file),
                   " ", shQuote(sam.file), " > ", shQuote(clean.file))
  expect_equal(system(command), 0)
  expect_identical(readLines(stats.file), "1\t1")
  expect_true(any(grepl("^soft", readLines(clean.file))))
  expect_false(any(grepl("^full", readLines(clean.file))))
})


test_that("contaminant indexes use checksums and require every BWA component", {
  root = tempfile("preprocess-index-")
  dir.create(root)
  reference.directory = file.path(root, "references")
  dir.create(reference.directory)
  reference.file = file.path(reference.directory, "test.fa")
  writeLines(c(">ref", "ACGT"), reference.file)

  fake.bwa = file.path(root, "bwa")
  writeLines(c(
    "#!/bin/sh",
    "shift",
    "if [ \"$1\" = '-p' ]; then prefix=\"$2\"; else prefix=\"$1\"; fi",
    "for suffix in amb ann bwt pac sa; do : > \"${prefix}.${suffix}\"; done"
  ), fake.bwa)
  Sys.chmod(fake.bwa, mode = "0755")

  with_preprocess_test_directory(root, {
    first.map = PhyloProcessR:::.buildContaminantIndex(
      reference.directory, shQuote(fake.bwa), quiet = TRUE
    )
    first.identity = attr(first.map, "reference.identity")
    writeLines(c(">ref", "TGCA"), reference.file)
    second.map = PhyloProcessR:::.buildContaminantIndex(
      reference.directory, shQuote(fake.bwa), quiet = TRUE
    )
    expect_false(identical(first.identity, attr(second.map, "reference.identity")))

    unlink(file.path(root, "ref-index", "reference.sa"))
    PhyloProcessR:::.buildContaminantIndex(
      reference.directory, shQuote(fake.bwa), quiet = TRUE
    )
    expect_true(file.exists(file.path(root, "ref-index", "reference.sa")))
  })
})


test_that("temporary reference indexes are removed when preprocessing exits", {
  root = tempfile("preprocess-temporary-indexes-")
  dir.create(root)
  input.directory = file.path(root, "input")
  reference.directory = file.path(root, "references")
  dir.create(input.directory)
  dir.create(reference.directory)
  target.file = file.path(root, "targets.fa")
  writeLines(c(">target", "ACGT"), target.file)

  fake.bwa = file.path(root, "bwa")
  writeLines(c(
    "#!/bin/sh",
    "shift",
    "prefix=\"$1\"",
    "for suffix in amb ann bwt pac sa; do : > \"${prefix}.${suffix}\"; done"
  ), fake.bwa)
  Sys.chmod(fake.bwa, mode = "0755")
  fake.samtools = file.path(root, "samtools")
  writeLines("#!/bin/sh", fake.samtools)
  Sys.chmod(fake.samtools, mode = "0755")

  with_preprocess_test_directory(root, {
    dir.create("ref-index")
    removeContamination(input.reads = input.directory,
                        output.directory = "decontaminated-reads",
                        decontamination.path = reference.directory,
                        bwa.path = fake.bwa, samtools.path = fake.samtools)
    expect_false(dir.exists("ref-index"))

    assessCaptureEfficiency(input.reads = input.directory,
                            output.directory = "sample-capture-assessment",
                            target.fasta = target.file,
                            bwa.path = fake.bwa, samtools.path = fake.samtools)
    expect_false(dir.exists("sample-capture-assessment/target-index"))
  })
})


test_that("zero contaminant counts retain typed report columns", {
  result = PhyloProcessR:::.readContaminantCounts(tempfile(),
                                                  data.frame(Contig = character()))
  expect_identical(names(result), c("Contaminant", "Accession", "Reads"))
  expect_equal(nrow(result), 0)
  expect_type(result$Reads, "double")
})


test_that("active contaminant references exclude cached extras", {
  root = tempfile("preprocess-references-")
  dir.create(root)
  list.file = file.path(root, "references.csv")
  writeLines("Genome,GenBank_Accession", list.file)
  local.fasta = file.path(root, "local.fa")
  writeLines(c(">local", "ACGT"), local.fasta)
  output.directory = file.path(root, "database")
  dir.create(output.directory)
  writeLines(c(">old", "AAAA"), file.path(output.directory, "old.fa"))

  result = createContaminantDB(decontamination.list = list.file,
                               output.directory = output.directory,
                               include.univec = FALSE,
                               include.fasta = local.fasta,
                               overwrite = FALSE)
  active = read.csv(file.path(output.directory, "active-references.csv"),
                    stringsAsFactors = FALSE)
  expect_identical(basename(result), "manually-included-data.fa")
  expect_identical(active$File, "manually-included-data.fa")
  expect_true(file.exists(file.path(output.directory, "old.fa")))
})
