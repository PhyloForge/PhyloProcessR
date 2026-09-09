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

test_that("curation orders, orients, gaps, and retains copies independently", {
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
    strrep("G", 60), strrep("T", 60)))
  names(contigs) = paste0("source", seq_along(contigs))
  Biostrings::writeXStringSet(contigs, file.path(input, "sample.fasta"))

  targets = Biostrings::DNAStringSet(rep(strrep("A", 200), 3))
  names(targets) = c("joined", "copies", "short")
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
    "short\tcontig_6\t98\t60\t0\t0\t1\t60\t1\t60\t0\t80\t200\t60\t0")
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
  expect_equal(formals(curateTargetContigs)$similarity, 0.9)
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

test_that("an early pipeline failure is reported", {
  expect_error(.runCommand("false | true", quiet = TRUE, task = "test pipeline"),
               "test pipeline.*failed")
})
