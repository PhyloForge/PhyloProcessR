test_that("consensus methods return documented result types", {
  alignment <- Biostrings::DNAStringSet(c(
    taxon1 = "ACGT",
    taxon2 = "A-NN",
    taxon3 = "ACGT",
    taxon4 = "ACGT"
  ))

  majority <- makeConsensus(alignment, method = "majority", remove.gaps = FALSE)
  profile <- makeConsensus(alignment, method = "profile")
  threshold <- makeConsensus(
    alignment,
    method = "threshold",
    threshold = 0.8,
    remove.gaps = FALSE
  )

  expect_equal(unname(as.character(majority)), "ACGT")
  expect_true(is.matrix(profile))
  expect_equal(unname(as.character(threshold)), "ANNN")
})

test_that("replaceAlignmentCharacter respects both character arguments", {
  alignment <- Biostrings::DNAStringSet(c(taxon1 = "ACGA", taxon2 = "AAGA"))

  observed <- replaceAlignmentCharacter(alignment, "A", "T")

  expect_equal(unname(as.character(observed)), c("TCGT", "TTGT"))
  expect_equal(names(observed), names(alignment))
  expect_error(
    replaceAlignmentCharacter(alignment, c("A", "C"), "T"),
    "must each contain one character"
  )
})

test_that("FASTA and strict PHYLIP writers produce valid headers", {
  fasta_file <- tempfile(fileext = ".fa")
  writeFasta(
    sequences = list(c("A", "C", "G", "T")),
    names = "taxon1",
    file.out = fasta_file,
    nbchar = 2
  )
  expect_equal(readLines(fasta_file), c(">taxon1", "AC", "GT"))

  phylip_file <- tempfile(fileext = ".phy")
  alignment <- matrix(
    c("A", "C", "G", "T"),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("abcdefghijk", "mnopqrstuvw"), NULL)
  )
  writePhylip(alignment, phylip_file, strict = TRUE)
  phylip <- readLines(phylip_file)

  expect_equal(phylip[1], "2 2")
  expect_match(phylip[2], "^abcdefghij")
  expect_match(phylip[3], "^mnopqrstuv")
})
