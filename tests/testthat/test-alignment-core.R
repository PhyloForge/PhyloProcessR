test_that("alignment gap statistics include gaps and missing bases", {
  alignment <- Biostrings::DNAStringSet(c(
    taxon1 = "ACGT",
    taxon2 = "A-NN",
    taxon3 = "ACGT",
    taxon4 = "ACGT"
  ))

  observed <- countAlignmentGaps(alignment)

  expect_equal(unname(observed), c(3, 16, 18.75))
  expect_named(observed, c("gaps", "alignment", "percent_gaps"))
})

test_that("informativeSites preserves DNAbin dimensions", {
  alignment <- as.matrix(ape::as.DNAbin(list(
    taxon1 = c("a", "c", "g"),
    taxon2 = c("a", "c", "g"),
    taxon3 = c("g", "c", "t"),
    taxon4 = c("g", "t", "t")
  )))

  expect_equal(informativeSites(alignment), 2)
  expect_equal(informativeSites(alignment, count = FALSE), 0.667)
})

test_that("row and column trimming use the requested gap threshold", {
  alignment <- Biostrings::DNAStringSet(c(
    taxon1 = "ACGT",
    taxon2 = "A-NN",
    taxon3 = "ACGT",
    taxon4 = "ACGT"
  ))

  row_trimmed <- trimAlignmentRows(alignment, min.gap.percent = 50)
  column_trimmed <- trimAlignmentColumns(alignment, min.gap.percent = 25)

  expect_equal(names(row_trimmed), c("taxon1", "taxon3", "taxon4"))
  expect_equal(unname(as.character(column_trimmed)), rep("A", 4))
})

test_that("alignmentAssess applies all three boundaries", {
  alignment <- Biostrings::DNAStringSet(c(
    taxon1 = "ACGTAC",
    taxon2 = "ACGTAC",
    taxon3 = "ACGTAC",
    taxon4 = "ACGTAC"
  ))

  expect_true(alignmentAssess(alignment, 10, 3, 5))
  expect_false(alignmentAssess(alignment, 0, 3, 5))
  expect_false(alignmentAssess(alignment, 10, 4, 5))
  expect_false(alignmentAssess(alignment, 10, 3, 6))
})

test_that("pairwiseDistanceTarget supports documented and inherited target types", {
  alignment <- Biostrings::DNAStringSet(c(
    taxon1 = "ACGT",
    taxon2 = "ACGA"
  ))
  target <- structure("taxon1", class = c("target_name", "character"))
  reference <- Biostrings::DNAStringSet(c(reference = "ACGT"))

  observed_name <- pairwiseDistanceTarget(alignment, target)
  observed_sequence <- pairwiseDistanceTarget(alignment, reference)

  expect_equal(unname(observed_name), c(0, 0.25))
  expect_named(observed_name, names(alignment))
  expect_equal(unname(observed_sequence), c(0, 0.25, 0))
  expect_named(observed_sequence, c(names(alignment), "reference"))
  expect_error(
    pairwiseDistanceTarget(alignment, 1),
    "target must be a DNAStringSet or a character sequence name"
  )
})
