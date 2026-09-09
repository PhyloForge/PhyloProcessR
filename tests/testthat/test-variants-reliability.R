test_that("workflow 3 resource budgets never produce zero heaps", {
  x = .validateResources(threads = 8, memory = 0.01, samples = 2)
  expect_equal(x$workers, 2)
  expect_gt(x$heap.mb, 0)
  pipe = .validateResources(threads = 2, memory = 1, samples = 2, simultaneous.jvms = 2)
  expect_lte(pipe$heap.mb * pipe$workers * 2, 1024)
  expect_error(.validateResources(0, 1, 1), "positive integer")
})

test_that("lane discovery excludes derived merged BAMs", {
  root = tempfile(); dir.create(root)
  dir.create(file.path(root, "Lane_2")); dir.create(file.path(root, "Lane_4"))
  dir.create(file.path(root, "Lane_Merge")); dir.create(file.path(root, "Lane_notes"))
  expect_setequal(basename(.laneDirectories(root)), c("Lane_2", "Lane_4"))
})

test_that("depth rules use strict cutoffs and preserve ambiguity semantics", {
  seqs = Biostrings::DNAStringSet(c(locus = "ARYNNN"))
  depth = tempfile(fileext = ".tsv")
  writeLines(paste("locus", 1:6, c(0,1,2,9,10,20), sep="\t"), depth)
  report = tempfile(fileext = ".tsv")
  default = .filterDepthSequences(seqs, depth, .validateDepthSettings("site", 1, 1, NULL), report)
  expect_equal(unname(as.character(default)), "NRYNNN")
  at10 = .filterDepthSequences(seqs, depth, .validateDepthSettings("site", 10, 1, NULL), report)
  expect_length(at10, 0)
  row = data.table::fread(report)
  expect_equal(row$preexisting_N, 3)
  expect_equal(row$newly_masked, 3)
  expect_equal(row$final_N_proportion, 1)
})

test_that("public scalar defaults and invalid depth combinations are explicit", {
  expect_equal(formals(VCFtoContigs)$vcf.file, "SNP")
  expect_equal(formals(genotypeSamples)$custom.SNP.QD, 2)
  expect_error(.validateDepthSettings("site", -1, 1, NULL), "min.site.depth")
})
