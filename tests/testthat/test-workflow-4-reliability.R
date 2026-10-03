random_dna <- function(length) {
  paste(sample(c("A", "C", "G", "T"), length, replace = TRUE), collapse = "")
}

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

# Writes contigs, targets and a mock LAST that prints the given hits, then runs
# annotateTargets in root. Returns the primary output and the evidence log.
run_annotation <- function(root, contigs, targets, hits, ...) {
  input <- file.path(root, "input")
  tools <- file.path(root, "tools")
  dir.create(input, recursive = TRUE)
  dir.create(tools)
  Biostrings::writeXStringSet(contigs, file.path(input, "sample.fa"))
  reference <- file.path(root, "targets.fa")
  Biostrings::writeXStringSet(targets, reference)
  writeLines(c("#!/bin/sh", "exit 0"), file.path(tools, "lastdb"))
  writeLines(c("#!/bin/sh", paste0("printf '%s\\n' ", shQuote(hits))),
             file.path(tools, "lastal"))
  Sys.chmod(file.path(tools, c("lastdb", "lastal")), "0755")

  old <- setwd(root)
  on.exit(setwd(old), add = TRUE)
  annotateTargets(
    input, reference, alignment.contig.name = file.path(root, "all"),
    output.directory = file.path(root, "output"),
    paralog.directory = file.path(root, "paralogs"),
    min.match.length = 50, min.match.coverage = 30,
    paralog.action = "exclude", last.path = tools, overwrite = TRUE, ...)
  list(primary = Biostrings::readDNAStringSet(file.path(root, "output", "sample.fa")),
       evidence = read.csv(file.path(root, "logs/sample_logs",
                                     "sample_target-candidates.csv")))
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
  writeLines(c("1 4", "sample ACGT"), sentinel)
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
  # Different sequences, so that the two copies of "close" are not one copy
  set.seed(4)
  contigs <- Biostrings::DNAStringSet(vapply(1:6, function(i) random_dna(60),
                                             character(1)))
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

test_that("annotateTargets gives each contig only to its best target", {
  set.seed(1)
  root <- tempfile("workflow4-best-target-")
  contigs <- Biostrings::DNAStringSet(c(random_dna(200), random_dna(200)))
  names(contigs) <- c("geneA", "geneB")
  targets <- Biostrings::DNAStringSet(c(geneA = random_dna(200),
                                        geneB = random_dna(200)))
  # Both targets match both contigs. Under the score, coverage and identity
  # tests alone, each target would see a competitive second copy.
  hits <- c(
    "geneA\tcontig_1\t95\t180\t0\t0\t1\t180\t1\t180\t0\t300\t200\t200\t0",
    "geneA\tcontig_2\t92\t180\t0\t0\t1\t180\t1\t180\t0\t280\t200\t200\t0",
    "geneB\tcontig_1\t92\t180\t0\t0\t1\t180\t1\t180\t0\t280\t200\t200\t0",
    "geneB\tcontig_2\t95\t180\t0\t0\t1\t180\t1\t180\t0\t300\t200\t200\t0")
  result <- run_annotation(root, contigs, targets, hits)

  expect_equal(sort(sub("_\\|_.*", "", names(result$primary))), c("geneA", "geneB"))
  expect_equal(unname(as.character(result$primary[["geneA_|_sample"]])),
               unname(as.character(contigs[["geneA"]])))
  other <- result$evidence[result$evidence$Target == "geneA" &
                             result$evidence$Contig == "contig_2", ]
  expect_equal(other$Decision, "other_target_better")
  expect_equal(other$Contig_best_target, "geneB")
})

test_that("annotateTargets keeps a target whose competitor is the same copy", {
  set.seed(2)
  root <- tempfile("workflow4-same-copy-")
  # Two overlapping fragments of one copy, and two copies that differ by
  # about 12 percent (a real paralog pair).
  locus <- random_dna(300)
  paralog <- random_dna(250)
  diverged <- strsplit(paralog, "")[[1]]
  changed <- seq(5, 250, by = 8)
  diverged[changed] <- ifelse(diverged[changed] == "A", "C", "A")
  contigs <- Biostrings::DNAStringSet(c(
    substr(locus, 1, 250), substr(locus, 51, 300),
    paralog, paste(diverged, collapse = "")))
  names(contigs) <- paste0("source", 1:4)
  targets <- Biostrings::DNAStringSet(c(fragmented = random_dna(300),
                                        duplicated = random_dna(250)))
  hits <- c(
    "fragmented\tcontig_1\t95\t250\t0\t0\t1\t250\t1\t250\t0\t400\t300\t250\t0",
    "fragmented\tcontig_2\t95\t250\t0\t0\t51\t300\t1\t250\t0\t390\t300\t250\t0",
    "duplicated\tcontig_3\t95\t250\t0\t0\t1\t250\t1\t250\t0\t400\t250\t250\t0",
    "duplicated\tcontig_4\t94\t250\t0\t0\t1\t250\t1\t250\t0\t380\t250\t250\t0")
  result <- run_annotation(root, contigs, targets, hits)

  expect_equal(sub("_\\|_.*", "", names(result$primary)), "fragmented")
  fragment <- result$evidence[result$evidence$Target == "fragmented" &
                                result$evidence$Candidate_rank == 2, ]
  expect_equal(fragment$Decision, "same_copy_as_best")
  expect_equal(fragment$Same_copy_identity, 100)
  copy <- result$evidence[result$evidence$Target == "duplicated" &
                            result$evidence$Candidate_rank == 2, ]
  expect_equal(copy$Decision, "competitive_paralog")
  expect_lt(copy$Same_copy_identity, 98)
})

test_that("annotateTargets renames redundant targets to their representative", {
  set.seed(3)
  root <- tempfile("workflow4-redundant-")
  dir.create(root)
  contigs <- Biostrings::DNAStringSet(c(source1 = random_dna(200)))
  targets <- Biostrings::DNAStringSet(c(uce = random_dna(200), busco = random_dna(200)))
  # Two redundant probes match the same contig with nearly equal scores
  hits <- c(
    "uce\tcontig_1\t95\t180\t0\t0\t1\t180\t1\t180\t0\t300\t200\t200\t0",
    "busco\tcontig_1\t95\t180\t0\t0\t1\t180\t1\t180\t0\t310\t200\t200\t0")
  map <- file.path(root, "redundant-targets.csv")
  write.csv(data.frame(Target = c("busco", "uce"), Representative = "uce",
                       Group_size = 2, Samples_with_target = 5),
            map, row.names = FALSE)
  result <- run_annotation(root, contigs, targets, hits, redundant.targets = map)

  expect_equal(names(result$primary), "uce_|_sample")
  expect_equal(result$evidence$Original_target, "busco")
  expect_error(annotateTargets(root, map, same.copy.identity = 101),
               "same.copy.identity")
})

test_that("alignTargets renames redundant genome targets to their representative", {
  root <- tempfile("workflow4-genome-redundant-")
  genomes <- file.path(root, "genomes")
  dir.create(genomes, recursive = TRUE)
  tools <- file.path(root, "tools")
  make_copy_mafft(tools)

  samples <- Biostrings::DNAStringSet(rep(strrep("ACGT", 25), 5))
  names(samples) <- paste0("uce_|_sample", 1:5)
  to.align <- file.path(root, "to-align.fa")
  Biostrings::writeXStringSet(samples, to.align)
  reference <- file.path(root, "targets.fa")
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(c(
    uce = strrep("ACGT", 25), busco = strrep("ACGT", 25))), reference)

  # Both redundant probes extracted the same region from the genome
  genome <- Biostrings::DNAStringSet(c(strrep("ACGT", 25), strrep("ACGT", 24)))
  names(genome) <- c("uce_|_genomeA", "busco_|_genomeA")
  Biostrings::writeXStringSet(genome, file.path(genomes, "genomeA_target-matches.fa"))
  map <- file.path(root, "redundant-targets.csv")
  write.csv(data.frame(Target = c("busco", "uce"), Representative = "uce",
                       Group_size = 2, Samples_with_target = 5),
            map, row.names = FALSE)

  output <- file.path(root, "alignments")
  alignTargets(to.align, reference, output.directory = output, min.taxa = 4,
               additional.sequence.directory = genomes, redundant.targets = map,
               mafft.path = tools, overwrite = TRUE)
  expect_equal(list.files(output, pattern = "\\.phy$"), "uce.phy")
  aligned <- readLines(file.path(output, "uce.phy"))
  expect_equal(sum(grepl("^genomeA", aligned)), 1)
})
