test_that("alignMACSE converts MACSE frameshift marks to gaps", {
  skip_on_os("windows")

  root <- tempfile()
  dir.create(root)
  input <- file.path(root, "input")
  output <- file.path(root, "output")
  bin <- file.path(root, "bin")
  dir.create(input)
  dir.create(bin)

  writeLines(c("3 9", "s1 ATGAACCCC", "s2 ATGAACCCC", "s3 ATGCCCCCC"),
             file.path(input, "locus1.phy"))

  # A stand-in for macse that writes a nucleotide alignment with "!" marks.
  fake.macse <- file.path(bin, "macse")
  writeLines(c(
    "#!/bin/sh",
    "while [ $# -gt 0 ]; do",
    "  case \"$1\" in",
    "    -out_NT) nt=\"$2\"; shift ;;",
    "    -out_AA) aa=\"$2\"; shift ;;",
    "  esac",
    "  shift",
    "done",
    "printf '>s1\\nATG!!-CCC\\n>s2\\nATGAA-CCC\\n>s3\\nATG---CCC\\n' > \"$nt\"",
    "printf '>s1\\nM\\n' > \"$aa\""
  ), fake.macse)
  Sys.chmod(fake.macse, "755")

  alignMACSE(alignment.folder = input, output.folder = output,
             macse.path = bin, threads = 1)

  result <- ape::read.dna(file.path(output, "locus1.phy"), format = "sequential",
                          as.character = TRUE)
  expect_equal(dim(result), c(3L, 9L))
  expect_equal(paste(result["s1", ], collapse = ""), "atg---ccc")
  # Temp files and logs are not left in the output folder.
  expect_identical(list.files(output, all.files = TRUE, no.. = TRUE), "locus1.phy")
})

test_that("alignMACSE aligns only markers with a gene in the metadata", {
  skip_on_os("windows")

  root <- tempfile()
  dir.create(root)
  input <- file.path(root, "input")
  output <- file.path(root, "output")
  bin <- file.path(root, "bin")
  dir.create(input)
  dir.create(bin)

  for (locus in c("exon1", "UCE_1")) {
    writeLines(c("2 6", "s1 ATGAAC", "s2 ATGAAC"),
               file.path(input, paste0(locus, ".phy")))
  }
  metadata <- file.path(root, "gene_metadata.txt")
  writeLines(c("marker\tgene", "exon1\tgeneA"), metadata)

  fake.macse <- file.path(bin, "macse")
  writeLines(c(
    "#!/bin/sh",
    "while [ $# -gt 0 ]; do",
    "  case \"$1\" in",
    "    -out_NT) nt=\"$2\"; shift ;;",
    "  esac",
    "  shift",
    "done",
    "printf '>s1\\nATGAAC\\n>s2\\nATGAAC\\n' > \"$nt\""
  ), fake.macse)
  Sys.chmod(fake.macse, "755")

  alignMACSE(alignment.folder = input, output.folder = output,
             feature.gene.names = metadata, macse.path = bin, threads = 1)

  expect_true(file.exists(file.path(output, "exon1.phy")))
  expect_false(file.exists(file.path(output, "UCE_1.phy")))
})

test_that("alignMACSE keeps only the log when MACSE fails", {
  skip_on_os("windows")

  root <- tempfile()
  dir.create(root)
  input <- file.path(root, "input")
  output <- file.path(root, "output")
  bin <- file.path(root, "bin")
  dir.create(input)
  dir.create(bin)

  writeLines(c("2 6", "s1 ATGAAC", "s2 ATGAAC"), file.path(input, "locus1.phy"))
  fake.macse <- file.path(bin, "macse")
  writeLines(c("#!/bin/sh", "echo 'MACSE error'", "exit 1"), fake.macse)
  Sys.chmod(fake.macse, "755")

  old_directory <- setwd(root)
  on.exit(setwd(old_directory), add = TRUE)
  expect_error(
    alignMACSE(alignment.folder = input, output.folder = output,
               macse.path = bin, threads = 1),
    "locus1"
  )
  expect_length(list.files(output, all.files = TRUE, no.. = TRUE), 0)
  expect_true(file.exists(file.path(root, "logs", "macse_logs", "locus1_macse.log")))
})
