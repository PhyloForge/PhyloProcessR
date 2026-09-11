make_copy_mafft <- function(directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  executable <- file.path(directory, "mafft")
  writeLines(c(
    "#!/bin/sh",
    "for arg in \"$@\"; do input=\"$arg\"; done",
    "cat \"$input\""
  ), executable)
  Sys.chmod(executable, "0755")
  executable
}

test_that("alignTargets uses exact locus IDs and diagnoses duplicate references", {
  root <- tempfile("workflow4-exact-")
  dir.create(root)
  tools <- file.path(root, "mock tools")
  make_copy_mafft(tools)
  input <- file.path(root, "input.fa")
  reference <- file.path(root, "reference.fa")
  output <- file.path(root, "alignments")

  sequences <- Biostrings::DNAStringSet(c(
    "gene1_|_a" = "ACGT", "gene1_|_b" = "ACGT",
    "xgene1_|_a" = "ACGT", "gene.1_|_a" = "ACGT",
    "geneX1_|_a" = "ACGT"))
  Biostrings::writeXStringSet(sequences, input)
  references <- Biostrings::DNAStringSet(c(
    gene1 = "ACGT", xgene1 = "ACGT", "gene.1" = "ACGT",
    geneX1 = "ACGT", "gene1_|_duplicate" = "ACGT"))
  Biostrings::writeXStringSet(references, reference)

  alignTargets(input, reference, output, min.taxa = 1, mafft.path = tools)
  log <- read.csv(file.path(output, "logs", "alignTargets_locus_summary.csv"))
  expect_equal(log$Status[log$Locus == "gene1"], "duplicate_reference")
  expect_true(file.exists(file.path(output, "xgene1.phy")))
  expect_equal(rownames(as.matrix(ape::read.dna(
    file.path(output, "xgene1.phy"), format = "sequential"))), "a")
})

test_that("alignTargets subsets are disjoint and partial overwrite is scoped", {
  root <- tempfile("workflow4-subsets-")
  dir.create(root)
  tools <- file.path(root, "tools")
  make_copy_mafft(tools)
  loci <- paste0("l", seq_len(8))
  sequences <- Biostrings::DNAStringSet(rep("ACGT", length(loci)))
  names(sequences) <- paste0(loci, "_|_sample")
  references <- Biostrings::DNAStringSet(rep("ACGT", length(loci)))
  names(references) <- loci
  input <- file.path(root, "input.fa")
  reference <- file.path(root, "reference.fa")
  output <- file.path(root, "output")
  Biostrings::writeXStringSet(sequences, input)
  Biostrings::writeXStringSet(references, reference)

  alignTargets(input, reference, output, min.taxa = 1, subset.end = 0.5,
               overwrite = TRUE, mafft.path = tools)
  expect_equal(sort(sub("\\.phy$", "", list.files(output, pattern = "\\.phy$"))),
               sort(loci[1:4]))
  sentinel <- file.path(output, "unrelated.phy")
  writeLines("sentinel", sentinel)
  alignTargets(input, reference, output, min.taxa = 1, subset.start = 0.5,
               overwrite = TRUE, mafft.path = tools)
  expect_true(file.exists(sentinel))
  expect_true(all(file.exists(file.path(output, paste0(loci[5:8], ".phy")))))

  empty.output <- file.path(root, "empty")
  expect_invisible(alignTargets(input, reference, empty.output, min.taxa = 1,
                                subset.start = 0.25, subset.end = 0.25,
                                mafft.path = tools))
  expect_length(list.files(empty.output, pattern = "\\.phy$"), 0)
  expect_error(alignTargets(input, reference, output, subset.start = 0.7,
                            subset.end = 0.2), "0 <= start <= end <= 1")
})

test_that("alignTargets includes extracted genome target sequences", {
  root <- tempfile("workflow4-genomes-")
  capture.file <- file.path(root, "capture.fa")
  target.file <- file.path(root, "targets.fa")
  genome.directory <- file.path(root, "genome-targets")
  genome.sample.directory <- file.path(genome.directory, "genome_sample")
  output.directory <- file.path(root, "alignments")
  tools <- file.path(root, "tools")
  dir.create(genome.sample.directory, recursive = TRUE)
  make_copy_mafft(tools)

  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c("locus1_|_capture_sample" = "ACGT")),
    capture.file
  )
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c("locus1_|_genome_sample" = "ACGA")),
    file.path(genome.sample.directory, "genome_sample_target-matches.fa")
  )
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(locus1 = "ACGT")),
    target.file
  )

  alignTargets(
    targets.to.align = capture.file,
    target.file = target.file,
    output.directory = output.directory,
    min.taxa = 1,
    removal.threshold = 1,
    mafft.path = tools,
    additional.sequence.directory = genome.directory
  )

  alignment <- ape::read.dna(
    file.path(output.directory, "locus1.phy"),
    format = "sequential"
  )
  expect_equal(
    sort(rownames(as.matrix(alignment))),
    c("capture_sample", "genome_sample")
  )
})

test_that("MAFFT failures cannot reuse partial output", {
  root <- tempfile("workflow4-mafft-failure-")
  dir.create(root)
  tools <- file.path(root, "tools")
  dir.create(tools)
  executable <- file.path(tools, "mafft")
  writeLines(c("#!/bin/sh", "printf '>partial\\nAAAA\\n'", "exit 7"), executable)
  Sys.chmod(executable, "0755")
  sequences <- Biostrings::DNAStringSet(c(a = "AAAA", b = "AAAA"))
  expect_error(runMafft(sequences, save.name = file.path(root, "result"),
                        mafft.path = tools), "MAFFT alignment.*failed")
  expect_false(file.exists(file.path(root, "result_align.fa")))
})

test_that("pairwiseDistanceTarget reports undefined comparisons as NA", {
  alignment <- Biostrings::DNAStringSet(c(reference = "ACGT----",
                                           sample = "----ACGT"))
  distance <- pairwiseDistanceTarget(alignment, "reference")
  expect_true(is.na(distance[["sample"]]))
})

test_that("alignTargets removes sequences with undefined reference distance", {
  root <- tempfile("workflow4-no-overlap-")
  dir.create(root)
  tools <- file.path(root, "tools")
  make_copy_mafft(tools)
  input <- file.path(root, "input.fa")
  reference <- file.path(root, "reference.fa")
  output <- file.path(root, "output")
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(c(
    "locus_|_overlap" = "ACGT----", "locus_|_disjoint" = "----ACGT")), input)
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(locus = "ACGT----")), reference)
  alignTargets(input, reference, output, min.taxa = 1, mafft.path = tools)
  result <- ape::read.dna(file.path(output, "locus.phy"), format = "sequential")
  expect_equal(rownames(as.matrix(result)), "overlap")
})

test_that("annotateTargets combines disjoint hit coverage for one contig", {
  root <- tempfile("workflow4-fragments-")
  input <- file.path(root, "input")
  output <- file.path(root, "output")
  tools <- file.path(root, "mock tools")
  dir.create(input, recursive = TRUE)
  dir.create(tools)
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(source = strrep("A", 200))),
    file.path(input, "sample.fa"))
  reference <- file.path(root, "targets.fa")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(locus = strrep("A", 500))), reference)

  writeLines(c("#!/bin/sh", "exit 0"), file.path(tools, "lastdb"))
  hits <- c(
    "locus\tcontig_1\t99\t80\t0\t0\t1\t80\t1\t80\t0\t80\t500\t200\t0",
    "locus\tcontig_1\t99\t80\t0\t0\t201\t280\t121\t200\t0\t80\t500\t200\t0")
  writeLines(c("#!/bin/sh", paste0("printf '%s\\n' ", shQuote(hits))),
             file.path(tools, "lastal"))
  Sys.chmod(file.path(tools, c("lastdb", "lastal")), "0755")

  old <- setwd(root)
  on.exit(setwd(old), add = TRUE)
  annotateTargets(input, reference, output.directory = output,
                  min.match.length = 50, min.match.coverage = 30,
                  last.path = tools, threads = 1,
                  overwrite = TRUE)
  result <- Biostrings::readDNAStringSet(file.path(output, "sample.fa"))
  expect_equal(length(result), 1)
  expect_equal(names(result), "locus_|_sample")
})

test_that("annotateTargets separates competitive copies from weak alternatives", {
  root <- tempfile("workflow4-paralogs-")
  input <- file.path(root, "input")
  tools <- file.path(root, "tools")
  dir.create(input, recursive = TRUE)
  dir.create(tools)
  contigs <- Biostrings::DNAStringSet(rep(strrep("A", 60), 6))
  names(contigs) <- paste0("source", seq_along(contigs))
  Biostrings::writeXStringSet(contigs, file.path(input, "sample.fa"))
  reference <- file.path(root, "targets.fa")
  targets <- Biostrings::DNAStringSet(rep(strrep("A", 100), 3))
  names(targets) <- c("close", "clear", "shortbest")
  Biostrings::writeXStringSet(targets, reference)

  writeLines(c("#!/bin/sh", "exit 0"), file.path(tools, "lastdb"))
  hits <- c(
    "close\tcontig_1\t99\t50\t0\t0\t1\t50\t1\t50\t0\t100\t100\t60\t0",
    "close\tcontig_2\t96\t45\t0\t0\t1\t45\t1\t45\t0\t90\t100\t60\t0",
    "clear\tcontig_3\t99\t60\t0\t0\t1\t60\t1\t60\t0\t100\t100\t60\t0",
    "clear\tcontig_4\t95\t30\t0\t0\t1\t30\t1\t30\t0\t30\t100\t60\t0",
    "shortbest\tcontig_5\t99\t30\t0\t0\t1\t30\t1\t30\t0\t50\t100\t60\t0",
    "shortbest\tcontig_6\t99\t20\t0\t0\t1\t20\t1\t20\t0\t45\t100\t60\t0")
  writeLines(c("#!/bin/sh", paste0("printf '%s\\n' ", shQuote(hits))),
             file.path(tools, "lastal"))
  Sys.chmod(file.path(tools, c("lastdb", "lastal")), "0755")

  old <- setwd(root)
  on.exit(setwd(old), add = TRUE)
  excluded <- file.path(root, "excluded")
  paralogs <- file.path(root, "paralogs")
  annotateTargets(
    input, reference, alignment.contig.name = file.path(root, "excluded-all"),
    output.directory = excluded, paralog.directory = paralogs,
    min.match.length = 20, min.match.coverage = 30,
    paralog.action = "exclude", last.path = tools, overwrite = TRUE)

  primary <- Biostrings::readDNAStringSet(file.path(excluded, "sample.fa"))
  expect_equal(sort(sub("_\\|_.*", "", names(primary))),
               c("clear", "shortbest"))
  saved <- Biostrings::readDNAStringSet(file.path(paralogs, "sample.fa"))
  expect_equal(length(saved), 4)
  evidence <- read.csv(file.path("logs/sample_logs",
                                 "sample_target-candidates.csv"))
  expect_true(evidence$Competitive_with_best[evidence$Target == "close" &
                                               evidence$Candidate_rank == 2])
  expect_false(evidence$Competitive_with_best[evidence$Target == "clear" &
                                                evidence$Candidate_rank == 2])
  below.floor <- evidence[evidence$Target == "shortbest" &
                            evidence$Target_coverage < 0.30, ]
  expect_false(below.floor$Passes_absolute_filters)
  best.short <- which(evidence$Target == "shortbest" &
                        !is.na(evidence$Candidate_rank) &
                        evidence$Candidate_rank == 1)
  expect_equal(evidence$Decision[best.short], "selected")

  best <- file.path(root, "best")
  annotateTargets(
    input, reference, alignment.contig.name = file.path(root, "best-all"),
    output.directory = best,
    paralog.directory = file.path(root, "best-paralogs"),
    min.match.length = 20, min.match.coverage = 30,
    paralog.action = "best", last.path = tools, overwrite = TRUE)
  best.primary <- Biostrings::readDNAStringSet(file.path(best, "sample.fa"))
  expect_equal(sort(sub("_\\|_.*", "", names(best.primary))),
               c("clear", "close", "shortbest"))
})
