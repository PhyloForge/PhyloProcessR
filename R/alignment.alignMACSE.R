#' @title alignMACSE
#'
#' @description Refines and aligns exon sequences using MACSE to ensure proper reading frames and codon alignment.
#'
#' @param alignment.folder character string; path to the folder containing the alignments.
#' @param output.folder character string; path to the folder to save the MACSE-refined alignments.
#' @param alignment.format character string; the format of the input alignments (e.g. "phylip", "fasta").
#' @param output.format character string; output alignment format. Currently
#'   used to name the intended output format; PHYLIP output is supported.
#' @param feature.gene.names character string; optional gene metadata file with
#'   "marker" and "gene" columns. If given, only markers with a gene are aligned,
#'   so UCEs and other non-coding markers are skipped.
#' @param macse.path character string; system path to the directory containing the macse executable. If NULL, searches the system PATH.
#' @param genetic.code integer; the genetic code table to use (default: 1 for standard nuclear).
#' @param threads integer; number of threads to use.
#' @param memory integer; reserved for compatibility. This value does not set the JVM memory limit.
#' @param overwrite logical; if TRUE, overwrite existing output files.
#' @param quiet logical; if TRUE, suppress output messages.
#'
#' @return A folder of MACSE-refined codon alignments. The MACSE log of each
#'   failed alignment is saved in logs/macse_logs.
#'
#' @export

alignMACSE = function(alignment.folder = NULL,
                      output.folder = NULL,
                      alignment.format = "phylip",
                      output.format = "phylip",
                      feature.gene.names = NULL,
                      macse.path = NULL,
                      genetic.code = 1,
                      threads = 1,
                      memory = 4,
                      overwrite = FALSE,
                      quiet = TRUE) {

  if (is.null(alignment.folder)) { stop("An input alignment folder must be provided.") }
  if (is.null(output.folder)) { stop("An output folder must be provided.") }

  # Creates output directory if it doesn't exist
  if (dir.exists(output.folder) == FALSE) {
    dir.create(output.folder, recursive = TRUE)
  }

  macse.command = if (is.null(macse.path)) "macse" else
    file.path(macse.path, "macse")

  # Gets all alignment files
  if (alignment.format == "phylip") {
    align.files = list.files(alignment.folder, full.names = TRUE, pattern = "\\.phy$")
  } else {
    align.files = list.files(alignment.folder, full.names = TRUE, pattern = "\\.fa$|\\.fasta$")
  }

  if (length(align.files) == 0) { stop("No alignments found in the input folder.") }

  # Keeps only the coding markers that have a gene in the metadata
  if (!is.null(feature.gene.names)) {
    metadata = .readGeneMetadata(feature.gene.names)
    exon.markers = metadata$marker[!is.na(metadata$gene) & nzchar(as.character(metadata$gene))]
    align.files = align.files[.alignmentId(align.files) %in% exon.markers]
    if (length(align.files) == 0) {
      stop("No alignment names match a marker with a gene in the gene metadata.")
    }
  }

  # Sets up foreach loop
  
  cl = makeCluster(threads)
  registerDoParallel(cl)
  on.exit({
    stopCluster(cl)
    foreach::registerDoSEQ()
  }, add = TRUE)

  cat(paste0("Refining ", length(align.files), " alignments using MACSE...\n"))

  # Aligns one locus. MACSE input, output and log files go to a scratch folder
  # that is always deleted, so the output folder holds only the final alignments.
  align.bases = .alignmentId(align.files)
  alignOne = function(align.file, file.base) {
    final.file = if (output.format == "phylip") {
      file.path(output.folder, paste0(file.base, ".phy"))
    } else {
      file.path(output.folder, paste0(file.base, ".fa"))
    }
    if (file.exists(final.file) && file.info(final.file)$size > 0 && !overwrite) {
      return(list(status = "skipped", locus = file.base))
    }

    work.dir = tempfile(paste0("macse_", file.base, "_"))
    dir.create(work.dir)
    on.exit(unlink(work.dir, recursive = TRUE), add = TRUE)
    out.file = file.path(work.dir, "out_NT.fa")
    out.aa.file = file.path(work.dir, "out_AA.fa")
    log.file = file.path(work.dir, "macse.log")

    # Keeps the MACSE log in logs/macse_logs only when the alignment fails
    fail = function(message) {
      log.directory = file.path("logs", "macse_logs")
      dir.create(log.directory, recursive = TRUE, showWarnings = FALSE)
      saved.log = file.path(log.directory, paste0(file.base, "_macse.log"))
      if (file.exists(log.file)) {
        file.copy(log.file, saved.log, overwrite = TRUE)
      }
      cat(paste0("\nPhyloProcessR: ", message, "\n"), file = saved.log,
          append = file.exists(saved.log))
      list(status = "error", locus = file.base, message = message)
    }

    # Records a locus that cannot be aligned, without stopping the other loci.
    skip = function(message) {
      log.directory = file.path("logs", "macse_logs")
      dir.create(log.directory, recursive = TRUE, showWarnings = FALSE)
      writeLines(message, file.path(log.directory, paste0(file.base, "_macse.log")))
      list(status = "excluded", locus = file.base, message = message)
    }

    # MACSE cannot align a sequence that contains only gaps or ambiguous bases.
    # Target trimming can create these rows when a taxon has no target data.
    if (alignment.format == "phylip") {
      align = Biostrings::DNAStringSet(Biostrings::readDNAMultipleAlignment(
        align.file, format = "phylip"
      ))
    } else {
      align = Biostrings::readDNAStringSet(align.file, format = "fasta")
    }
    sequence.strings = as.character(align)
    usable = grepl("[ACGT]", sequence.strings, ignore.case = TRUE)
    align = align[usable]
    if (length(align) < 2) {
      return(skip(paste0(
        "Only ", length(align),
        " sequence(s) contain a definite nucleotide after missing-data rows were removed."
      )))
    }
    macse.input = file.path(work.dir, "input.fa")
    Biostrings::writeXStringSet(align, macse.input, format = "fasta")

    macse.cmd = paste0(
      shQuote(macse.command), " -prog alignSequences ",
      "-seq ", shQuote(macse.input), " ",
      "-gc_def ", genetic.code, " ",
      "-out_NT ", shQuote(out.file), " ",
      "-out_AA ", shQuote(out.aa.file)
    )
    status = system(paste0(macse.cmd, " > ", shQuote(log.file), " 2>&1"))
    if (status != 0) {
      return(fail(paste0("MACSE exited with status ", status, ".")))
    }
    if (!file.exists(out.file) || file.info(out.file)$size == 0) {
      return(fail("MACSE did not create a nucleotide alignment."))
    }

    # MACSE marks frameshifts with "!". ape drops that character, which gives
    # rows of different length, so change it to a gap.
    macse.lines = readLines(out.file)
    seq.lines = !startsWith(macse.lines, ">")
    macse.lines[seq.lines] = gsub("!", "-", macse.lines[seq.lines], fixed = TRUE)
    writeLines(macse.lines, out.file)

    if (output.format == "phylip") {
      align.macse = ape::read.FASTA(out.file, type = "DNA")
      align.mat = tryCatch(as.matrix(align.macse), error = function(e) NULL)
      if (is.null(align.mat)) {
        return(fail("MACSE sequences are not all the same length."))
      }
      rownames(align.mat) = labels(align.macse)
      .writePhylipAtomic(align.mat, final.file)
    } else {
      file.copy(out.file, final.file, overwrite = TRUE)
    }
    list(status = "success", locus = file.base)
  }

  results = foreach(i = seq_along(align.files),
                    .packages = c("Biostrings", "ape"),
                    .export = ".writePhylipAtomic") %dopar% {
    alignOne(align.files[i], align.bases[i])
  }

  cat("MACSE refinement complete.\n")
}
