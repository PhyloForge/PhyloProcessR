source("workflow-4_configuration-file.R")
if (isTRUE(get0("install.latest.github", ifnotfound = FALSE))) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    stop("Install the remotes package to use install.latest.github = TRUE.")
  }
  remotes::install_github("PhyloForge/PhyloProcessR", upgrade = "never",
                          dependencies = FALSE)
}
library(PhyloProcessR)
setwd(working.directory)

##################################################################################################
##################################################################################################
## Annotation and paralog filtering
##################################################################################################

if (file.exists("data-analysis/contigs") == FALSE){ dir.create("data-analysis/contigs") }

if (heterozygote.filter == TRUE){
  #Remove contigs with too much heterozygosity
  filterHeterozygosity(
    iupac.directory = contig.directory,
    output.directory = "data-analysis/contigs/7_filtered-contigs",
    removed.directory = "data-analysis/contigs/6_removed-contigs",
    threshold = heterozygote.filter.threshold,
    min.length = heterozygote.min.length,
    threads = threads,
    memory = memory,
    overwrite = overwrite
  )

  input.contigs = "data-analysis/contigs/7_filtered-contigs"

} else { input.contigs = contig.directory }

if (annotate.targets == TRUE) {
  # annotates targets
  annotateTargets(
    assembly.directory = input.contigs,
    target.file = target.file,
    alignment.contig.name = paste0("data-analysis/", dataset.name),
    output.directory = "data-analysis/contigs/8_annotated-contigs",
    min.match.percent = min.match.percent,
    min.match.length = min.match.length,
    min.match.coverage = min.match.coverage,
    paralog.action = paralog.action,
    paralog.score.ratio = paralog.score.ratio,
    paralog.coverage.ratio = paralog.coverage.ratio,
    paralog.identity.delta = paralog.identity.delta,
    paralog.directory = "data-analysis/contigs/9_paralog-contigs",
    threads = threads,
    memory = memory,
    overwrite = overwrite,
    quiet = quiet,
    blast.path = blast.path,
    last.path = last.path,
    search.method = annotate.search.method
  )
}#end if

# Create alignments folder
dir.create("data-analysis/alignments", recursive = TRUE, showWarnings = FALSE)

if (align.targets == TRUE) {
  # Aligns target markers from annotation files
  alignTargets(
    targets.to.align = paste0("data-analysis/", dataset.name, "_to-align.fa"),
    target.file = target.file,
    additional.sequence.directory = if (include.genomes == TRUE) {
      genome.target.directory
    } else {
      NULL
    },
    output.directory = "data-analysis/alignments/untrimmed_all-markers",
    min.taxa = min.taxa.alignment,
    removal.threshold = removal.threshold,
    algorithm = alignment.algorithm,
    subset.start = subset.start,
    subset.end = subset.end,
    threads = threads,
    memory = memory,
    overwrite = overwrite,
    quiet = quiet,
    mafft.path = mafft.path
  )
}#end if

if (include.legacy == TRUE) {
  alignment.output.directory = "data-analysis/alignments/untrimmed_all-markers"

  readLegacyRenameTable = function(path) {
    if (!file.exists(path)) {
      stop("Legacy rename file not found: ", path)
    }

    extension = tolower(tools::file_ext(path))
    if (extension %in% c("xls", "xlsx")) {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("Install the readxl package to use an XLS or XLSX rename file.")
      }
      rename.data = as.data.frame(
        readxl::read_excel(path, col_names = TRUE),
        stringsAsFactors = FALSE
      )
    } else {
      separator = if (extension == "csv") {
        ","
      } else if (extension == "tsv") {
        "\t"
      } else if (extension == "txt") {
        first.line = readLines(path, n = 1, warn = FALSE)
        if (grepl("\t", first.line)) {
          "\t"
        } else if (grepl(",", first.line)) {
          ","
        } else {
          ""
        }
      } else {
        stop("legacy.rename.file must be a CSV, TSV, TXT, XLS, or XLSX file.")
      }
      rename.data = utils::read.table(
        path,
        header = TRUE,
        sep = separator,
        quote = "\"",
        comment.char = "",
        stringsAsFactors = FALSE,
        fill = TRUE
      )
    }

    if (ncol(rename.data) < 2 || nrow(rename.data) == 0) {
      stop(
        "Legacy rename file must contain at least two columns and one mapping."
      )
    }
    if (all(c("Legacy_Name", "SeqCap_Name") %in% names(rename.data))) {
      rename.data = rename.data[, c("Legacy_Name", "SeqCap_Name"), drop = FALSE]
    } else {
      rename.data = rename.data[, 1:2, drop = FALSE]
    }
    names(rename.data) = c("Legacy_Name", "SeqCap_Name")
    rename.data[] = lapply(rename.data, function(values) {
      trimws(as.character(values))
    })

    invalid = is.na(rename.data$Legacy_Name) |
      is.na(rename.data$SeqCap_Name) |
      rename.data$Legacy_Name == "" |
      rename.data$SeqCap_Name == ""
    if (any(invalid)) {
      stop("Legacy rename file contains a blank legacy or sequence-capture name.")
    }
    if (anyDuplicated(rename.data$Legacy_Name)) {
      stop("Each Legacy_Name must occur only once in the legacy rename file.")
    }
    rename.data
  }

  mergeLegacyRows = function(alignment, rename.data) {
    applied = character()
    missing.characters = c("-", "N", "n", "?")

    for (i in seq_len(nrow(rename.data))) {
      legacy.name = rename.data$Legacy_Name[i]
      seqcap.name = rename.data$SeqCap_Name[i]
      legacy.index = which(names(alignment) == legacy.name)
      if (length(legacy.index) == 0) next
      if (length(legacy.index) > 1) {
        stop("Legacy alignment contains more than one row named: ", legacy.name)
      }

      seqcap.index = which(names(alignment) == seqcap.name)
      if (length(seqcap.index) > 1) {
        stop("Legacy alignment contains more than one row named: ", seqcap.name)
      }

      if (length(seqcap.index) == 0 || legacy.name == seqcap.name) {
        names(alignment)[legacy.index] = seqcap.name
      } else {
        legacy.bases = strsplit(
          as.character(alignment[legacy.index]), ""
        )[[1]]
        seqcap.bases = strsplit(
          as.character(alignment[seqcap.index]), ""
        )[[1]]
        use.legacy = seqcap.bases %in% missing.characters &
          !legacy.bases %in% missing.characters
        seqcap.bases[use.legacy] = legacy.bases[use.legacy]
        merged = Biostrings::DNAStringSet(paste0(seqcap.bases, collapse = ""))
        names(merged) = seqcap.name
        alignment[seqcap.index] = merged
        alignment = alignment[-legacy.index]
      }
      applied = c(applied, legacy.name)
    }

    if (anyDuplicated(names(alignment))) {
      stop("A legacy rename creates duplicate sample names in an alignment.")
    }
    list(alignment = alignment, applied = applied)
  }

  if (!dir.exists(legacy.alignment.directory)) {
    stop("Workflow X3 legacy alignment directory not found: ",
         legacy.alignment.directory)
  }

  legacy.files = list.files(
    legacy.alignment.directory,
    pattern = "[.]phy$",
    full.names = TRUE
  )
  if (length(legacy.files) == 0) {
    stop("No PHYLIP alignments found in workflow X3 output: ",
         legacy.alignment.directory)
  }

  rename.data = NULL
  if (!is.null(legacy.rename.file)) {
    rename.data = readLegacyRenameTable(legacy.rename.file)
  }

  dir.create(alignment.output.directory, recursive = TRUE,
             showWarnings = FALSE)
  copied = file.copy(
    from = legacy.files,
    to = file.path(alignment.output.directory, basename(legacy.files)),
    overwrite = TRUE
  )
  if (any(!copied)) {
    stop("Could not copy all workflow X3 alignments into: ",
         alignment.output.directory)
  }

  if (!is.null(rename.data)) {
    renamed.samples = character()

    for (alignment.file in file.path(
      alignment.output.directory, basename(legacy.files)
    )) {
      alignment = Biostrings::DNAStringSet(
        Biostrings::readDNAMultipleAlignment(alignment.file, format = "phylip")
      )
      result = mergeLegacyRows(alignment, rename.data)
      if (length(result$applied) == 0) next

      alignment.matrix = as.matrix(ape::as.DNAbin(
        strsplit(as.character(result$alignment), "")
      ))
      PhyloProcessR::writePhylip(
        alignment = alignment.matrix,
        file = alignment.file,
        interleave = FALSE,
        strict = FALSE
      )
      renamed.samples = union(renamed.samples, result$applied)
    }

    unused.names = setdiff(rename.data$Legacy_Name, renamed.samples)
    if (length(unused.names) > 0) {
      warning(
        "Legacy rename entries not found in the imported alignments: ",
        paste(unused.names, collapse = ", ")
      )
    }
  }

  print(paste0(
    "Added ", length(legacy.files),
    " workflow X3 alignment(s) to ", alignment.output.directory, "."
  ))
}#end include.legacy
