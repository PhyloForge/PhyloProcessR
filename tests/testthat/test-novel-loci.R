test_that("novel contig collection keeps the longest contig and counts input contigs", {
  root = tempfile("novel collection ")
  dir.create(root)
  old.directory = setwd(root)
  on.exit(setwd(old.directory), add = TRUE)
  dir.create("contigs")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(chr1_0_20_contig_1 = "ACGTACGT",
                               chr1_0_20_contig_2 = "ACGTACGTACGT",
                               chr2_0_20_contig_1 = "AC")),
    "contigs/sample1.fa")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(chr1_0_20_contig_1 = "ACGTACGT")),
    "contigs/sample2.fa")

  result = collectNovelContigs("contigs", "nested/novel", min.contig.length = 4,
                               min.taxa = 2)
  output = Biostrings::readDNAStringSet("nested/novel_to-align.fa")
  expect_equal(length(output), 2)
  expect_equal(as.character(output[["chr1_0_20_|_sample1"]]), "ACGTACGTACGT")
  expect_equal(result$InputContigs, c(3L, 1L))
  expect_equal(result$LociPassingFilter, c(1L, 1L))

  result = collectNovelContigs("contigs", "nested/novel", min.contig.length = 100,
                               min.taxa = 2, overwrite = TRUE)
  expect_equal(file.size("nested/novel_to-align.fa"), 0)
  expect_equal(result$InputContigs, c(3L, 1L))
  expect_equal(result$LociPassingFilter, c(0L, 0L))
  expect_true(file.exists("logs/novel_locus_summary.csv"))
})


test_that("novel contig collection writes an empty result when occupancy fails", {
  root = tempfile("novel occupancy ")
  dir.create(root)
  old.directory = setwd(root)
  on.exit(setwd(old.directory), add = TRUE)
  dir.create("contigs")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(chr1_0_20_contig_1 = "ACGT")), "contigs/sample.fa")
  result = collectNovelContigs("contigs", "novel", min.contig.length = 1, min.taxa = 2)
  expect_equal(file.size("novel_to-align.fa"), 0)
  expect_equal(result$LociPassingFilter, 0L)
})


test_that("shared region selection clears old results when too few samples remain", {
  root = tempfile("shared empty ")
  dir.create(root)
  bed = file.path(root, "sample.bed")
  output = file.path(root, "shared.bed")
  writeLines("chr1\t0\t100", bed)
  writeLines("chr1\t0\t100\t4", output)
  PhyloProcessR:::.sharedCoveredRegions(bed, output, 2, 1, 0, "unused")
  expect_equal(file.size(output), 0)
})


test_that("shared regions require coverage on the same bases", {
  skip_if(Sys.which("bedtools") == "", "bedtools is not installed")
  root = tempfile("shared intervals ")
  dir.create(root)
  beds = file.path(root, c("a.bed", "b.bed", "c.bed"))
  writeLines("chr1\t0\t100", beds[1])
  writeLines("chr1\t90\t190", beds[2])
  writeLines("chr1\t180\t280", beds[3])
  output = file.path(root, "shared.bed")
  command = PhyloProcessR:::.toolCommand("bedtools")

  PhyloProcessR:::.sharedCoveredRegions(beds, output, 3, 1, 0, command)
  expect_equal(file.size(output), 0)
  PhyloProcessR:::.sharedCoveredRegions(beds, output, 2, 1, 0, command)
  expect_equal(readLines(output), c("chr1\t90\t100\t2", "chr1\t180\t190\t2"))
  PhyloProcessR:::.sharedCoveredRegions(beds, output, 2, 1, 80, command)
  expect_equal(readLines(output), "chr1\t90\t190\t2")
  PhyloProcessR:::.sharedCoveredRegions(beds[1], output, 1, 1, 0, command)
  expect_equal(readLines(output), "chr1\t0\t100\t1")
})


test_that("shared region selection reports a failed tool without replacing results", {
  root = tempfile("shared failure ")
  dir.create(root)
  beds = file.path(root, c("a.bed", "b.bed"))
  for (bed in beds) writeLines("chr1\t0\t100", bed)
  output = file.path(root, "shared.bed")
  writeLines("original", output)
  expect_error(
    PhyloProcessR:::.sharedCoveredRegions(beds, output, 2, 1, 0, "false"),
    "shared region selection.*failed")
  expect_equal(readLines(output), "original")
})


test_that("assembly reports sample failures and cannot reuse an overwritten FASTA", {
  root = tempfile("novel assembly ")
  dir.create(root)
  dir.create(file.path(root, "sample-bams"))
  output = file.path(root, "contigs")
  dir.create(output)
  writeLines("chr1\t0\t100\t1", file.path(root, "novel_regions.bed"))
  writeLines(c(">chr1_0_100", "ACGT"), file.path(root, "novel_targets.fa"))
  writeLines("placeholder", file.path(root, "sample-bams", "sample.bam"))
  old.fasta = file.path(output, "sample.fa")
  writeLines(c(">old", "ACGT"), old.fasta)

  testthat::local_mocked_bindings(
    .toolCommand = function(program, program.path = NULL) program,
    .runCommand = function(command, ...) {
      if (startsWith(command, "makeblastdb")) return(invisible(0L))
      stop("Read extraction failed")
    },
    .package = "PhyloProcessR")
  expect_error(assembleSharedRegions(root, output, overwrite = TRUE),
               "Assembly failed for one or more samples")
  expect_false(file.exists(old.fasta))
  expect_false(file.exists(file.path(output, ".sample.complete")))
})


test_that("discovery checks both indexes and passes all paired lanes to HISAT2", {
  root = tempfile("novel discovery ")
  dir.create(root)
  alignments = file.path(root, "alignments")
  reads = file.path(root, "reads")
  output = file.path(root, "output")
  dir.create(alignments)
  dir.create(reads)
  dir.create(output)
  sample.dir = file.path(reads, "sample")
  dir.create(sample.dir)
  writeLines("placeholder", file.path(alignments, "known.phy"))
  genome = file.path(root, "genome.fa")
  writeLines(c(">chr1", "ACGT"), genome)
  writeLines(c(">known", "ACGT"), file.path(output, "known_loci_consensus.fa"))
  for (i in 1:8) {
    writeLines("index", file.path(output, paste0("known_loci_index.", i, ".ht2l")))
  }
  for (lane in c("L001", "L002")) {
    for (mate in 1:2) {
      writeLines(c("@read", "ACGT", "+", "IIII"),
                 file.path(sample.dir, paste0("sample_", lane, "_R", mate, ".fastq")))
    }
  }
  commands = character()
  testthat::local_mocked_bindings(
    .toolCommand = function(program, program.path = NULL) program,
    .runCommand = function(command, ...) {
      commands <<- c(commands, command)
      if (startsWith(command, "hisat2-build")) return(invisible(0L))
      stop("Mapping failed")
    },
    .package = "PhyloProcessR")

  expect_error(discoverSharedRegions(alignments, read.directory = reads,
                                     genome.file = genome, output.directory = output),
               "Discovery failed for one or more samples")
  expect_length(commands, 2)
  expect_match(commands[1], "genome_index", fixed = TRUE)
  expect_false(grepl("known_loci_index", commands[1], fixed = TRUE))
  for (lane in c("L001", "L002")) {
    for (mate in 1:2) {
      expect_match(commands[2], paste0("sample_", lane, "_R", mate, ".fastq"),
                   fixed = TRUE)
    }
  }
  expect_false(file.exists(file.path(output, "novel_regions.bed")))
})
