with_test_working_directory <- function(path, code) {
  old_directory <- setwd(path)
  on.exit(setwd(old_directory), add = TRUE)
  force(code)
}

write_test_phylip <- function(path, sequences) {
  stopifnot(length(unique(nchar(sequences))) == 1L)
  writeLines(
    c(
      paste(length(sequences), nchar(sequences[[1]])),
      paste(names(sequences), unname(sequences))
    ),
    path
  )
}

find_calls_named <- function(expression, function_name) {
  matches <- list()
  if (is.call(expression) && identical(expression[[1]], as.name(function_name))) {
    matches <- list(expression)
  }
  if (is.recursive(expression)) {
    for (component in as.list(expression)) {
      matches <- c(matches, find_calls_named(component, function_name))
    }
  }
  matches
}

find_workflow_file <- function(filename, start_directory = getwd()) {
  directory <- normalizePath(start_directory, mustWork = TRUE)
  repeat {
    candidates <- c(
      file.path(directory, "workflows", filename),
      file.path(
        directory, "00_pkg_src", "PhyloProcessR", "workflows", filename
      )
    )
    existing <- candidates[file.exists(candidates)]
    if (length(existing) > 0) {
      return(existing[[1]])
    }
    parent <- dirname(directory)
    if (identical(parent, directory)) {
      return(character())
    }
    directory <- parent
  }
}

test_that("alignTargets creates and overwrites nested output directories", {
  root <- tempfile("align-targets-directory-")
  dir.create(root)

  with_test_working_directory(root, {
    sequences_file <- file.path(root, "sequences.fa")
    targets_file <- file.path(root, "targets.fa")
    output_directory <- file.path(root, "nested", "alignment", "output")

    Biostrings::writeXStringSet(
      Biostrings::DNAStringSet(c("locus1_|_taxon1" = "ACGT")),
      sequences_file
    )
    Biostrings::writeXStringSet(
      Biostrings::DNAStringSet(c(locus1 = "ACGT")),
      targets_file
    )

    alignTargets(
      targets.to.align = sequences_file,
      target.file = targets_file,
      output.directory = output_directory,
      min.taxa = 4,
      overwrite = FALSE
    )
    expect_true(dir.exists(output_directory))

    sentinel <- file.path(output_directory, "remove-on-overwrite.txt")
    writeLines("old output", sentinel)
    alignTargets(
      targets.to.align = sequences_file,
      target.file = targets_file,
      output.directory = output_directory,
      min.taxa = 4,
      overwrite = TRUE
    )

    expect_true(dir.exists(output_directory))
    expect_false(file.exists(sentinel))
  })
})

test_that("trimAlignmentTargets creates and overwrites nested output directories", {
  root <- tempfile("trim-targets-directory-")
  dir.create(root)

  with_test_working_directory(root, {
    input_directory <- file.path(root, "input")
    output_directory <- file.path(root, "nested", "target", "output")
    targets_file <- file.path(root, "targets.fa")
    dir.create(input_directory)

    write_test_phylip(
      file.path(input_directory, "locus1.phy"),
      c(taxon1 = "ACGT", taxon2 = "ACGT", taxon3 = "ACGT",
        taxon4 = "ACGT", taxon5 = "ACGT")
    )
    Biostrings::writeXStringSet(
      Biostrings::DNAStringSet(c(different_locus = "ACGT")),
      targets_file
    )

    trimAlignmentTargets(
      alignment.directory = input_directory,
      output.directory = output_directory,
      target.file = targets_file,
      threads = 1,
      overwrite = FALSE
    )
    expect_true(dir.exists(output_directory))

    sentinel <- file.path(output_directory, "remove-on-overwrite.txt")
    writeLines("old output", sentinel)
    trimAlignmentTargets(
      alignment.directory = input_directory,
      output.directory = output_directory,
      target.file = targets_file,
      threads = 1,
      overwrite = TRUE
    )

    expect_true(dir.exists(output_directory))
    expect_false(file.exists(sentinel))
  })
})

test_that("superTrimmer creates and overwrites nested output directories", {
  root <- tempfile("super-trimmer-directory-")
  dir.create(root)

  with_test_working_directory(root, {
    input_directory <- file.path(root, "input")
    output_directory <- file.path(root, "nested", "trimmed", "output")
    dir.create(input_directory)
    write_test_phylip(
      file.path(input_directory, "locus1.phy"),
      c(taxon1 = "ACGT", taxon2 = "ACGT", taxon3 = "ACGT", taxon4 = "ACGT")
    )

    run_trimmer <- function(overwrite) {
      superTrimmer(
        alignment.dir = input_directory,
        output.dir = output_directory,
        TrimAl = FALSE,
        trim.similarity = FALSE,
        trim.external = FALSE,
        trim.coverage = FALSE,
        trim.column = FALSE,
        alignment.assess = FALSE,
        min.taxa.alignment = 0,
        threads = 1,
        overwrite = overwrite
      )
    }

    run_trimmer(FALSE)
    expect_true(file.exists(file.path(output_directory, "locus1.phy")))

    sentinel <- file.path(output_directory, "remove-on-overwrite.txt")
    writeLines("old output", sentinel)
    run_trimmer(TRUE)

    expect_true(file.exists(file.path(output_directory, "locus1.phy")))
    expect_false(file.exists(sentinel))
  })
})

test_that("Workflow 5 forwards alignment assessment configuration", {
  workflow_file <- find_workflow_file("workflow-5_trimming.R")
  expect_gt(length(workflow_file), 0)
  if (length(workflow_file) == 0) {
    return(invisible(NULL))
  }

  expressions <- parse(workflow_file, keep.source = FALSE)
  trimmer_calls <- find_calls_named(expressions, "superTrimmer")

  expect_length(trimmer_calls, 4)
  for (trimmer_call in trimmer_calls) {
    arguments <- as.list(trimmer_call)[-1]
    expect_true("alignment.assess" %in% names(arguments))
    expect_true("max.alignment.gap.percent" %in% names(arguments))
    expect_identical(arguments[["alignment.assess"]], quote(alignment.assess))
    expect_identical(
      arguments[["max.alignment.gap.percent"]],
      quote(max.alignment.gap.percent)
    )
  }
})

test_that("workflow GitHub installation is explicit and disabled by default", {
  workflow_directory <- dirname(find_workflow_file("workflow-5_trimming.R"))
  workflow_scripts <- list.files(
    workflow_directory,
    pattern = "^workflow-.*[.]R$",
    full.names = TRUE
  )
  workflow_scripts <- workflow_scripts[
    !grepl("configuration-file[.]R$", workflow_scripts)
  ]
  configuration_files <- list.files(
    workflow_directory,
    pattern = "configuration-file[.]R$",
    full.names = TRUE
  )

  expect_length(workflow_scripts, 9)
  expect_length(configuration_files, 9)

  for (configuration_file in configuration_files) {
    configuration <- paste(readLines(configuration_file, warn = FALSE),
                           collapse = "\n")
    expect_match(
      configuration,
      "install[.]latest[.]github[[:space:]]*=[[:space:]]*FALSE"
    )
  }

  for (workflow_script in workflow_scripts) {
    workflow <- readLines(workflow_script, warn = FALSE)
    workflow_text <- paste(workflow, collapse = "\n")
    source_line <- grep("^source\\(", workflow)[[1]]
    option_line <- grep("get0\\(\"install[.]latest[.]github\"", workflow)[[1]]

    expect_lt(source_line, option_line)
    expect_match(
      workflow_text,
      'remotes::install_github\\("PhyloForge/PhyloProcessR"'
    )
    expect_false(grepl("chutter/PhyloProcessR", workflow_text, fixed = TRUE))
    expect_false(grepl("devtools::install_github", workflow_text, fixed = TRUE))
  }
})
