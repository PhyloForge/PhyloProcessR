random.dna = function(length) {
  paste(sample(c("A", "C", "G", "T"), length, replace = TRUE), collapse = "")
}

test_that("LAST failures are distinct from successful searches with no hits", {
  root = tempfile("last-status-")
  dir.create(root)
  query = file.path(root, "query.fa")
  writeLines(c(">q", "AAAA"), query)

  expect_error(
    .lastSearch(query.file = query, db.prefix = file.path(root, "db"),
                out.file = file.path(root, "failed.txt"),
                lastal.command = "false", quiet = TRUE),
    "LAST search.*failed"
  )

  empty.last = file.path(root, "empty-last")
  writeLines(c("#!/bin/sh", "printf '# header\\n'"), empty.last)
  Sys.chmod(empty.last, "0755")
  out = file.path(root, "empty.txt")
  .lastSearch(query.file = query, db.prefix = file.path(root, "db"),
              out.file = out, lastal.command = shQuote(empty.last), quiet = TRUE)
  expect_true(file.exists(out))
  expect_equal(file.size(out), 0)
})

test_that("curation orders, orients, gaps, retains flanks, and keeps copies", {
  root = tempfile("curation-known-answer-")
  input = file.path(root, "input")
  output = file.path(root, "output")
  tools = file.path(root, "mock tools")
  dir.create(input, recursive = TRUE)
  dir.create(output)
  dir.create(tools)

  contigs = Biostrings::DNAStringSet(c(
    strrep("C", 60), strrep("A", 60),
    strrep("A", 120), strrep("C", 110),
    strrep("G", 60), strrep("T", 60),
    paste0(strrep("C", 30), strrep("A", 100), strrep("G", 30)),
    strrep("A", 300), strrep("A", 200)))
  names(contigs) = paste0("source", seq_along(contigs))
  Biostrings::writeXStringSet(contigs, file.path(input, "sample.fasta"))

  targets = Biostrings::DNAStringSet(rep(strrep("A", 200), 8))
  names(targets) = c("joined", "copies", "short", "flanked",
                     "left-target", "right-target", "overlap-one",
                     "overlap-two")
  target.file = file.path(root, "targets.fa")
  Biostrings::writeXStringSet(targets, target.file)

  cdhit = file.path(tools, "cd-hit-est")
  writeLines(c(
    "#!/bin/sh", "while [ $# -gt 0 ]; do",
    "case \"$1\" in -i) in_file=\"$2\"; shift 2;; -o) out_file=\"$2\"; shift 2;; *) shift;; esac",
    "done", "cp \"$in_file\" \"$out_file\""), cdhit)
  lastdb = file.path(tools, "lastdb")
  writeLines(c("#!/bin/sh", "exit 0"), lastdb)
  lastal = file.path(tools, "lastal")
  hits = c(
    "joined\tcontig_1\t99\t60\t0\t0\t141\t200\t60\t1\t0\t100\t200\t60\t0",
    "joined\tcontig_2\t99\t60\t0\t0\t1\t60\t1\t60\t0\t100\t200\t60\t0",
    "copies\tcontig_3\t99\t120\t0\t0\t1\t120\t1\t120\t0\t120\t200\t120\t0",
    "copies\tcontig_4\t98\t110\t0\t0\t1\t110\t1\t110\t0\t110\t200\t110\t0",
    "short\tcontig_5\t99\t60\t0\t0\t1\t60\t1\t60\t0\t90\t200\t60\t0",
    "short\tcontig_6\t98\t60\t0\t0\t1\t60\t1\t60\t0\t80\t200\t60\t0",
    "flanked\tcontig_7\t99\t100\t0\t0\t51\t150\t31\t130\t0\t100\t200\t160\t0",
    "left-target\tcontig_8\t99\t100\t0\t0\t1\t100\t31\t130\t0\t100\t200\t300\t0",
    "right-target\tcontig_8\t99\t100\t0\t0\t1\t100\t171\t270\t0\t100\t200\t300\t0",
    "overlap-one\tcontig_9\t99\t100\t0\t0\t1\t100\t31\t130\t0\t100\t200\t200\t0",
    "overlap-two\tcontig_9\t99\t100\t0\t0\t1\t100\t71\t170\t0\t100\t200\t200\t0")
  writeLines(c("#!/bin/sh", paste0("printf '%s\\n' ", shQuote(hits))), lastal)
  Sys.chmod(c(cdhit, lastdb, lastal), "0755")

  old = setwd(root)
  on.exit(setwd(old), add = TRUE)
  curateTargetContigs(
    assembly.directory = input, target.file = target.file,
    output.directory = output, min.match.percent = 60,
    min.match.length = 50, min.match.coverage = 50,
    search.method = "last", threads = 1, memory = 1,
    cdhit.path = tools, last.path = tools, overwrite = TRUE, quiet = TRUE)

  result = Biostrings::readDNAStringSet(file.path(output, "sample.fa"))
  expect_equal(unname(as.character(result["joined"])),
               paste0(strrep("A", 60), strrep("N", 80), strrep("G", 60)))
  expect_equal(sum(sub("_.*$", "", names(result)) == "copies"), 2)
  expect_false(any(sub("_.*$", "", names(result)) == "short"))
  expect_equal(Biostrings::width(result["flanked"]), 160)
  expect_equal(Biostrings::width(result["left-target"]), 150)
  expect_equal(Biostrings::width(result["right-target"]), 150)
  # The two targets match the same region with equal scores, so the name that
  # sorts first keeps the contig and the other target gets nothing.
  expect_equal(Biostrings::width(result["overlap-one"]), 200)
  expect_false(any(names(result) == "overlap-two"))

  # The match log keeps the matches of both targets for findRedundantTargets
  matches = read.csv(file.path(root, "logs/sample_logs/sample_curation-matches.csv"))
  expect_true("overlap-two" %in% matches$qName)
})

test_that("curation gives a shared region to the target with the higher score", {
  root = tempfile("curation-best-target-")
  input = file.path(root, "input")
  output = file.path(root, "output")
  tools = file.path(root, "tools")
  dir.create(input, recursive = TRUE)
  dir.create(output)
  dir.create(tools)

  # contig_1 is gene A and contig_2 is gene B. Both targets match both contigs,
  # as the targets of two genes in one family do.
  contigs = Biostrings::DNAStringSet(c(strrep("A", 150), strrep("C", 150)))
  names(contigs) = c("source1", "source2")
  Biostrings::writeXStringSet(contigs, file.path(input, "sample.fa"))
  targets = Biostrings::DNAStringSet(rep(strrep("A", 150), 2))
  names(targets) = c("geneA", "geneB")
  target.file = file.path(root, "targets.fa")
  Biostrings::writeXStringSet(targets, target.file)

  cdhit = file.path(tools, "cd-hit-est")
  writeLines(c(
    "#!/bin/sh", "while [ $# -gt 0 ]; do",
    "case \"$1\" in -i) in_file=\"$2\"; shift 2;; -o) out_file=\"$2\"; shift 2;; *) shift;; esac",
    "done", "cp \"$in_file\" \"$out_file\""), cdhit)
  writeLines(c("#!/bin/sh", "exit 0"), file.path(tools, "lastdb"))
  hits = c(
    "geneA\tcontig_1\t99\t150\t0\t0\t1\t150\t1\t150\t0\t280\t150\t150\t0",
    "geneA\tcontig_2\t85\t150\t0\t0\t1\t150\t1\t150\t0\t200\t150\t150\t0",
    "geneB\tcontig_1\t85\t150\t0\t0\t1\t150\t1\t150\t0\t200\t150\t150\t0",
    "geneB\tcontig_2\t99\t150\t0\t0\t1\t150\t1\t150\t0\t280\t150\t150\t0")
  writeLines(c("#!/bin/sh", paste0("printf '%s\\n' ", shQuote(hits))),
             file.path(tools, "lastal"))
  Sys.chmod(c(cdhit, file.path(tools, c("lastdb", "lastal"))), "0755")

  old = setwd(root)
  on.exit(setwd(old), add = TRUE)
  curateTargetContigs(
    assembly.directory = input, target.file = target.file,
    output.directory = output, min.match.percent = 60,
    min.match.length = 50, min.match.coverage = 50,
    search.method = "last", threads = 1, memory = 1,
    cdhit.path = tools, last.path = tools, overwrite = TRUE, quiet = TRUE)

  result = Biostrings::readDNAStringSet(file.path(output, "sample.fa"))
  expect_equal(sort(names(result)), c("geneA", "geneB"))
  expect_equal(unname(as.character(result["geneA"])), strrep("A", 150))
  expect_equal(unname(as.character(result["geneB"])), strrep("C", 150))
})

test_that("findRedundantTargets groups targets that share a best contig", {
  root = tempfile("redundant-targets-")
  logs = file.path(root, "logs")
  dir.create(logs, recursive = TRUE)

  # In each sample: A, B and C match the same region of contig_1, so they form
  # one group. D and E share a region in one sample only, which is below
  # min.samples. F and G match different regions of contig_4 (neighbouring
  # exons), so they are not grouped.
  for (sample in paste0("s", 1:3)) {
    log = data.frame(
      Sample = sample,
      qName = c("A", "B", "C", "D", "E", "F", "G"),
      tName = c("contig_1", "contig_1", "contig_1", "contig_2",
                if (sample == "s1") "contig_2" else "contig_3",
                "contig_4", "contig_4"),
      pident = 99, matches = 100, bitscore = c(100, 120, 90, 100, 100, 100, 100),
      evalue = 0, qLen = 100, tLen = 400,
      tStart = c(1, 11, 21, 1, 1, 1, 201), tEnd = c(100, 110, 120, 100, 100, 100, 300))
    write.csv(log, file.path(logs, paste0(sample, "_curation-matches.csv")),
              row.names = FALSE)
  }
  # One more sample has only target A, so A is the representative
  write.csv(data.frame(Sample = "s4", qName = "A", tName = "contig_1", pident = 99,
                       matches = 100, bitscore = 100, evalue = 0, qLen = 100,
                       tLen = 100, tStart = 1, tEnd = 100),
            file.path(logs, "s4_curation-matches.csv"), row.names = FALSE)

  output = file.path(root, "output")
  result = findRedundantTargets(log.directory = logs, output.directory = output,
                                min.samples = 3, min.shared.proportion = 0.5)
  expect_equal(sort(result$Target), c("A", "B", "C"))
  expect_true(all(result$Representative == "A"))
  pairs = read.csv(file.path(output, "shared-contig-pairs.csv"))
  expect_false(pairs$Grouped[pairs$Target_a == "D" & pairs$Target_b == "E"])
  expect_false(any(pairs$Target_a == "F" & pairs$Target_b == "G"))
  expect_error(findRedundantTargets(log.directory = logs, min.samples = 0),
               "min.samples")
})

test_that("curation similarity is configurable and validated", {
  root = tempfile("curation-similarity-")
  dir.create(root)
  target.file = file.path(root, "targets.fa")
  writeLines(c(">target", "AAAA"), target.file)
  expect_error(curateTargetContigs(
    assembly.directory = root, target.file = target.file,
    output.directory = file.path(root, "output"), similarity = 0.79),
    "similarity must")

  function.body = paste(deparse(body(curateTargetContigs)), collapse = "\n")
  expect_match(function.body, '" -c ", similarity', fixed = TRUE)
  expect_equal(formals(curateTargetContigs)$similarity, 0.98)
})

test_that("SPAdes restart accepts valid completed outputs", {
  root = tempfile("spades-restart-")
  reads = file.path(root, "reads", "sample")
  output = file.path(root, "work")
  assembly = file.path(root, "assembly")
  dir.create(reads, recursive = TRUE)
  dir.create(assembly, recursive = TRUE)
  writeLines(c("@r", "AAAA", "+", "IIII"), file.path(reads, "sample_R1.fastq"))
  writeLines(c("@r", "TTTT", "+", "IIII"), file.path(reads, "sample_R2.fastq"))
  writeLines(c(">done", "AAAA"), file.path(assembly, "sample.fa"))

  expect_invisible(assembleSpades(
    input.reads = dirname(reads), output.directory = output,
    assembly.directory = assembly, spades.path = file.path(root, "missing"),
    overwrite = FALSE))
})

test_that("read pairs require identical lane keys and ignore sidecars", {
  root = tempfile("lane-pairs-")
  dir.create(root)
  files = c("sample_lane1_R1.fastq", "sample_lane1_R2.fastq",
            "sample_lane10_R1.fastq", "sample_lane10_R2.fastq",
            "sample_lane1_R1.fastq.md5")
  file.create(file.path(root, files))
  pair = .pairSampleReads(root)
  expect_equal(basename(pair$read1), c("sample_lane1_R1.fastq",
                                       "sample_lane10_R1.fastq"))
  expect_equal(basename(pair$read2), c("sample_lane1_R2.fastq",
                                       "sample_lane10_R2.fastq"))

  unlink(file.path(root, "sample_lane10_R2.fastq"))
  file.create(file.path(root, "sample_lane2_R2.fastq"))
  expect_null(.pairSampleReads(root))
})

test_that("assembly preserves numeric sample names and separate lanes", {
  root = tempfile("numeric-read-names-")
  dir.create(root)
  reads = c("Test_dataset_1_L001_READ1.fastq.gz",
            "Test_dataset_1_L001_READ2.fastq.gz",
            "Test_dataset_1_L002_READ1.fastq.gz",
            "Test_dataset_1_L002_READ2.fastq.gz")
  paths = file.path(root, reads)
  file.create(paths)
  pair = .pairSampleReads(root)
  expect_identical(basename(pair$read1), reads[c(1, 3)])
  expect_identical(basename(pair$read2), reads[c(2, 4)])
  expect_identical(.readMateNumbers(paths), c(1L, 2L, 1L, 2L))
})

test_that("an early pipeline failure is reported", {
  expect_error(.runCommand("false | true", quiet = TRUE, task = "test pipeline"),
               "test pipeline.*failed")
})

test_that("curation keeps one piece when two pieces of a target are one copy", {
  set.seed(5)
  root = tempfile("curation-same-copy-")
  input = file.path(root, "input")
  output = file.path(root, "output")
  tools = file.path(root, "tools")
  dir.create(input, recursive = TRUE)
  dir.create(output)
  dir.create(tools)

  # contig_1 and contig_2 are overlapping fragments of one copy. contig_3 and
  # contig_4 are two copies that differ at every eighth base.
  locus = random.dna(300)
  paralog = random.dna(200)
  diverged = strsplit(paralog, "")[[1]]
  changed = seq(5, 200, by = 8)
  diverged[changed] = ifelse(diverged[changed] == "A", "C", "A")
  contigs = Biostrings::DNAStringSet(c(substr(locus, 1, 200), substr(locus, 101, 300),
                                       paralog, paste(diverged, collapse = "")))
  names(contigs) = paste0("source", 1:4)
  Biostrings::writeXStringSet(contigs, file.path(input, "sample.fa"))
  targets = Biostrings::DNAStringSet(c(fragmented = random.dna(300), copies = random.dna(200)))
  target.file = file.path(root, "targets.fa")
  Biostrings::writeXStringSet(targets, target.file)

  cdhit = file.path(tools, "cd-hit-est")
  writeLines(c(
    "#!/bin/sh", "while [ $# -gt 0 ]; do",
    "case \"$1\" in -i) in_file=\"$2\"; shift 2;; -o) out_file=\"$2\"; shift 2;; *) shift;; esac",
    "done", "cp \"$in_file\" \"$out_file\""), cdhit)
  writeLines(c("#!/bin/sh", "exit 0"), file.path(tools, "lastdb"))
  hits = c(
    "fragmented\tcontig_1\t99\t200\t0\t0\t1\t200\t1\t200\t0\t300\t300\t200\t0",
    "fragmented\tcontig_2\t99\t200\t0\t0\t101\t300\t1\t200\t0\t290\t300\t200\t0",
    "copies\tcontig_3\t99\t200\t0\t0\t1\t200\t1\t200\t0\t300\t200\t200\t0",
    "copies\tcontig_4\t95\t200\t0\t0\t1\t200\t1\t200\t0\t280\t200\t200\t0")
  writeLines(c("#!/bin/sh", paste0("printf '%s\\n' ", shQuote(hits))),
             file.path(tools, "lastal"))
  Sys.chmod(c(cdhit, file.path(tools, c("lastdb", "lastal"))), "0755")

  old = setwd(root)
  on.exit(setwd(old), add = TRUE)
  curateTargetContigs(
    assembly.directory = input, target.file = target.file,
    output.directory = output, min.match.percent = 60,
    min.match.length = 50, min.match.coverage = 30,
    search.method = "last", threads = 1, memory = 1,
    cdhit.path = tools, last.path = tools, overwrite = TRUE, quiet = TRUE)

  result = Biostrings::readDNAStringSet(file.path(output, "sample.fa"))
  expect_equal(sum(sub("_[0-9]+$", "", names(result)) == "fragmented"), 1)
  expect_equal(sum(sub("_[0-9]+$", "", names(result)) == "copies"), 2)
  expect_error(curateTargetContigs(
    assembly.directory = input, target.file = target.file,
    output.directory = output, same.copy.identity = 101), "same.copy.identity")
})
