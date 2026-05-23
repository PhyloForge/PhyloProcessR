#' @title convertNexusPartitions
#'
#' @description Splits a single concatenated NEXUS alignment into separate per-locus
#' alignment files using the \code{charset} partition definitions in the file's
#' \code{BEGIN SETS} block. For each charset, the corresponding columns are extracted
#' from the full matrix, samples with excessive missing data are optionally removed, and
#' the alignment is written as a separate file. Multi-range charsets (e.g.
#' \code{charset COI = 1-100 200-300}) are supported: the constituent ranges are
#' concatenated in order. Output files are named after the charset (e.g.
#' \code{16S.phy}).
#'
#' @param nexus.file full path to the input NEXUS file containing a concatenated
#' alignment and a \code{BEGIN SETS} block with \code{charset} definitions.
#'
#' @param output.directory path to the directory where per-locus alignment files
#' will be written.
#'
#' @param output.format format for the output alignment files. Accepted: "phylip" or
#' "fasta". Default "phylip".
#'
#' @param min.taxa.alignment minimum number of taxa required to write an alignment.
#' Loci with fewer taxa after missing-data filtering are skipped. Default 4.
#'
#' @param max.missing.percent maximum percent of missing or gap characters (\code{?},
#' \code{-}, \code{N}) allowed per sample per locus. Samples exceeding this threshold
#' are removed from that locus's alignment before writing. Default 100 (keep all
#' samples; set to e.g. 50 to remove samples with more than half missing data).
#'
#' @param overwrite logical. If TRUE, deletes and recreates the output directory.
#' Default FALSE.
#'
#' @param quiet logical. If TRUE, suppresses per-locus progress messages. Default FALSE.
#'
#' @return Writes one alignment file per charset partition to \code{output.directory}.
#' No value is returned to R.
#'
#' @export

convertNexusPartitions = function(nexus.file = NULL,
                                   output.directory = NULL,
                                   output.format = "phylip",
                                   min.taxa.alignment = 4,
                                   max.missing.percent = 100,
                                   overwrite = FALSE,
                                   quiet = FALSE) {

  # nexus.file = "/Users/chutter/Downloads/Glassfrog_matrix_2026.nex"
  # output.directory = "data-analysis/legacy-integration/legacy-alignments"
  # output.format = "phylip"
  # min.taxa.alignment = 4
  # max.missing.percent = 100
  # overwrite = TRUE
  # quiet = FALSE

  #Input checks
  if (is.null(nexus.file)) { print("nexus.file not provided."); return(NULL) }
  if (is.null(output.directory)) { print("output.directory not provided."); return(NULL) }
  if (!file.exists(nexus.file)) {
    print(paste0("nexus.file not found: ", nexus.file)); return(NULL)
  }

  #Overwrite
  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE) {
      system(paste0("rm -r ", output.directory))
    } else {
      stop("Overwrite = FALSE and output directory exists. Either change to TRUE or overwrite manually.")
    }
  }
  dir.create(output.directory)

  ##################################################################################################
  ## Step 1: Read the full concatenated alignment
  ##################################################################################################
  print("Reading NEXUS alignment...")
  nex.data = ape::read.nexus.data(nexus.file)
  # nex.data is a named list; each element is a character vector of individual bases
  n.taxa = length(nex.data)
  n.chars = length(nex.data[[1]])
  print(paste0("  ", n.taxa, " taxa, ", n.chars, " total characters."))

  ##################################################################################################
  ## Step 2: Parse charset partitions from the SETS block
  ##################################################################################################
  nex.lines = readLines(nexus.file)
  charset.lines = nex.lines[grepl("charset", nex.lines, ignore.case = TRUE)]

  if (length(charset.lines) == 0) {
    print("No charset definitions found in the NEXUS file. Verify the BEGIN SETS block exists.")
    return(NULL)
  }

  # Extract gene/locus names: everything between 'charset' and '='
  gene.names = trimws(gsub("(?i)^[[:space:]]*charset[[:space:]]+([^=]+)=.*$",
                           "\\1", charset.lines, perl = TRUE))

  # Extract all start-end ranges per charset (supports multi-range charsets)
  ranges.list = lapply(charset.lines, function(line) {
    raw = regmatches(line, gregexpr("[0-9]+[[:space:]]*-[[:space:]]*[0-9]+", line))[[1]]
    if (length(raw) == 0) { return(NULL) }
    do.call(rbind, lapply(raw, function(r) {
      parts = as.integer(strsplit(trimws(r), "[[:space:]]*-[[:space:]]*")[[1]])
      c(start = parts[1], end = parts[2])
    }))
  })

  # Drop any lines that failed to parse
  valid = !sapply(ranges.list, is.null) & nchar(gene.names) > 0
  gene.names  = gene.names[valid]
  ranges.list = ranges.list[valid]

  print(paste0("Found ", length(gene.names), " charset partitions."))

  ##################################################################################################
  ## Step 3: Extract, filter, and write each partition
  ##################################################################################################
  for (i in seq_along(gene.names)) {

    gene    = gene.names[i]
    ranges  = ranges.list[[i]]

    # Concatenate columns across all ranges for this charset
    gene.seqs = lapply(nex.data, function(seq) {
      unlist(lapply(seq_len(nrow(ranges)), function(r) {
        seq[ranges[r, "start"]:ranges[r, "end"]]
      }))
    })

    gene.len = sum(ranges[, "end"] - ranges[, "start"] + 1)

    # Remove samples exceeding max.missing.percent
    if (max.missing.percent < 100) {
      miss.prop = sapply(gene.seqs, function(seq) {
        mean(tolower(seq) %in% c("-", "?", "n"))
      })
      gene.seqs = gene.seqs[miss.prop <= (max.missing.percent / 100)]
    }

    # Replace '?' with 'n': '?' is valid NEXUS missing-data notation but is not a
    # recognised IUPAC code and will cause Biostrings to error downstream.
    gene.seqs = lapply(gene.seqs, function(seq) { seq[seq == "?"] = "n"; seq })

    # Skip if too few taxa
    if (length(gene.seqs) <= min.taxa.alignment) {
      if (quiet == FALSE) {
        print(paste0(gene, ": only ", length(gene.seqs),
                     " taxa after filtering. Skipping (min.taxa.alignment = ",
                     min.taxa.alignment, ")."))
      }
      next
    }

    # Convert list of character vectors → DNAbin matrix
    align.mat = as.matrix(ape::as.DNAbin(gene.seqs))

    # Write output
    if (output.format == "phylip") {
      PhyloProcessR::writePhylip(alignment = align.mat,
                                  file = paste0(output.directory, "/", gene, ".phy"),
                                  interleave = FALSE,
                                  strict = FALSE)
    } else {
      # fasta: convert DNAbin rows back to strings
      seq.strings = apply(as.character(align.mat), 1, paste0, collapse = "")
      out.seqs = Biostrings::DNAStringSet(seq.strings)
      Biostrings::writeXStringSet(out.seqs,
                                   filepath = paste0(output.directory, "/", gene, ".fa"))
    }

    if (quiet == FALSE) {
      print(paste0("Written: ", gene, " (", length(gene.seqs), " taxa, ", gene.len, " bp)"))
    }

    rm(gene.seqs, align.mat)
    gc()

  }#end gene loop

  print(paste0("Done. ", length(gene.names), " loci written to: ", output.directory))

}#end function

#END SCRIPT
