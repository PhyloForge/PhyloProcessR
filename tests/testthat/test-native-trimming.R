test_that("native HMM and entropy trimming preserve a clean alignment", {
  alignment <- Biostrings::DNAStringSet(c(
    taxon1 = "ATGAAATAA",
    taxon2 = "ATGAAATAA",
    taxon3 = "ATGAAATAA"
  ))

  hmm_trimmed <- trimSampleHMM(alignment)
  entropy_trimmed <- trimStatisticalColumns(
    alignment,
    method = "entropy",
    max.entropy = 0
  )

  expect_equal(hmm_trimmed, alignment)
  expect_equal(entropy_trimmed, alignment)
})

test_that("trimExonORF preserves an in-frame coding alignment", {
  alignment <- Biostrings::DNAStringSet(c(
    taxon1 = "ATGAAATAA",
    taxon2 = "ATGAAATAA",
    taxon3 = "ATGAAATAA"
  ))

  observed <- trimExonORF(alignment, locus.name = "exon_001")

  expect_equal(observed, alignment)
  expect_true(all(Biostrings::width(observed) %% 3 == 0))
})

test_that("native trimmers leave small alignments unchanged", {
  alignment <- Biostrings::DNAStringSet(c(taxon1 = "ACGT", taxon2 = "ACGT"))

  expect_identical(trimSampleHMM(alignment), alignment)
  expect_identical(trimExonORF(alignment), alignment)
})

test_that("trimExternal reports an uncovered alignment with a numeric sentinel", {
  alignment <- Biostrings::DNAStringSet(c(
    taxon1 = "----",
    taxon2 = "----",
    taxon3 = "----",
    taxon4 = "----"
  ))

  observed <- trimExternal(alignment, min.n.seq = 4)

  expect_true(is.numeric(observed))
  expect_equal(observed, 0)
})
