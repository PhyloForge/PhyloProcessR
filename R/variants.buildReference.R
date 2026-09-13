#' Build one dataset-owned shared reference for joint genotyping
#'
#' @description Produces a single reference FASTA and its BWA, samtools, and GATK
#'   index sidecars in a dataset-owned directory, together with a build record.
#'   Every sample in a joint-genotyping run maps and is genotyped against this
#'   one reference, so the contig names, sequences, and coordinates stay the same
#'   for all samples. Three modes select the reference source:
#'   \describe{
#'     \item{consensus}{One majority consensus sequence per phylip alignment in
#'       `alignment.directory` (the default; the historical workflow X1 source).}
#'     \item{target}{The sequence-capture target markers FASTA in `target.file`.}
#'     \item{user}{A reference FASTA the user supplies in `reference.file`.}
#'   }
#'   The build record stores the mode and the md5 sum of every input, so a resumed
#'   run reuses an unchanged reference. When an input changed while a reference
#'   already exists, the function stops instead of silently mixing coordinate
#'   systems; rerun with `overwrite = TRUE` or use a new dataset directory.
#'
#' @param reference.path path to the reference FASTA to create, for example
#'   `data-analysis/<dataset>/reference/reference.fa`.
#' @param reference.mode one of `consensus`, `target`, or `user`.
#' @param alignment.directory phylip alignment directory for `consensus` mode.
#' @param target.file capture-target markers FASTA for `target` mode.
#' @param reference.file user-supplied reference FASTA for `user` mode.
#' @param samtools.path,bwa.path,gatk4.path tool directories, executable paths,
#'   or NULL to search the system PATH.
#' @param threads number of alignment workers used to build the consensus.
#' @param overwrite logical; if TRUE the reference directory is rebuilt.
#' @param quiet logical; if TRUE tool stdout/stderr is suppressed while logs are
#'   retained.
#'
#' @return invisibly the reference path.
#'
#' @export

buildReference = function(reference.path = NULL,
                          reference.mode = c("consensus", "target", "user"),
                          alignment.directory = NULL,
                          target.file = NULL,
                          reference.file = NULL,
                          samtools.path = NULL,
                          bwa.path = NULL,
                          gatk4.path = NULL,
                          threads = 1,
                          overwrite = FALSE,
                          quiet = TRUE) {

  reference.mode = match.arg(reference.mode)
  if (is.null(reference.path)) { stop("Please provide reference.path.") }

  # Resolves the source input for the selected mode
  if (reference.mode == "consensus") {
    if (is.null(alignment.directory) || !dir.exists(alignment.directory)) {
      stop("consensus mode needs an existing alignment.directory.")
    }
    source.files = list.files(alignment.directory, pattern = "\\.(phy|phylip)$",
                              full.names = TRUE)
    if (length(source.files) == 0) {
      stop("No phylip alignment files (.phy or .phylip) were found in ", alignment.directory)
    }
  } else if (reference.mode == "target") {
    if (is.null(target.file) || !file.exists(target.file)) {
      stop("target mode needs an existing target.file.")
    }
    source.files = target.file
  } else {
    if (is.null(reference.file) || !file.exists(reference.file)) {
      stop("user mode needs an existing reference.file.")
    }
    source.files = reference.file
  }

  reference.dir = dirname(reference.path)
  record.path = file.path(reference.dir, "reference-build.rds")
  sidecars = paste0(reference.path, c(".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"))
  dict.path = sub("\\.fa$", ".dict", reference.path)
  reference.outputs = c(reference.path, sidecars, dict.path)

  # Describes the current request from the mode and the md5 of every input file
  source.files = sort(source.files)
  current.record = list(mode = reference.mode,
                        inputs = tools::md5sum(source.files))

  # Resume: reuse an unchanged reference, refuse a changed one, rebuild on overwrite
  if (file.exists(record.path) && overwrite == FALSE) {
    previous.record = readRDS(record.path)
    same.inputs = identical(previous.record$mode, current.record$mode) &&
      identical(previous.record$inputs, current.record$inputs)
    if (same.inputs && all(file.exists(reference.outputs))) {
      return(invisible(reference.path))
    }
    if (!same.inputs) {
      stop("The reference inputs changed since the existing reference was built ",
           "in ", reference.dir, ".\nRerun with overwrite = TRUE to rebuild the ",
           "reference and all downstream results, or use a new dataset directory.")
    }
  }
  if (file.exists(reference.path) && overwrite == FALSE && !file.exists(record.path)) {
    stop("A reference already exists in ", reference.dir, " without a build record. ",
         "Rerun with overwrite = TRUE to rebuild it.")
  }

  # Builds the owned reference directory from scratch
  if (dir.exists(reference.dir)) { unlink(reference.dir, recursive = TRUE) }
  .ensureDirectory(reference.dir, "reference directory")

  if (reference.mode == "consensus") {
    # Locus IDs keep any internal periods; only the final extension is removed.
    locus.ids = sub("\\.(phy|phylip)$", "", basename(source.files))

    # Builds one majority consensus per alignment, recording any read failure.
    results = parallel::mclapply(seq_along(source.files), function(i) {
      tryCatch({
        align = Biostrings::DNAStringSet(
          Biostrings::readDNAMultipleAlignment(file = source.files[i], format = "phylip"))
        if (length(align) == 0) { stop("the alignment has no sequences") }
        list(sequence = as.character(makeConsensus(align)), error = NA_character_)
      }, error = function(e) list(sequence = NA_character_, error = conditionMessage(e)))
    }, mc.cores = threads)

    # Stops on any unreadable alignment instead of building a smaller reference
    failed = vapply(results, function(r) !is.na(r$error), logical(1))
    if (any(failed)) {
      stop("Could not read alignment(s): ",
           paste0(basename(source.files[failed]), " (",
                  vapply(results[failed], function(r) r$error, character(1)), ")", collapse = "; "))
    }

    sequences = vapply(results, function(r) r$sequence, character(1))

    # An empty or all-missing consensus has no informative reference sequence
    informative = grepl("[ACGT]", sequences, ignore.case = TRUE)
    if (any(nchar(sequences) == 0 | !informative)) {
      bad = locus.ids[nchar(sequences) == 0 | !informative]
      stop("These loci produced an empty or all-missing consensus: ", paste(bad, collapse = ", "))
    }

    # Locus names must be usable as file names and GATK interval names
    if (any(!nzchar(locus.ids)) || any(grepl("\\s|/", locus.ids))) {
      stop("Locus names must be non-empty and contain no spaces or slashes.")
    }
    if (anyDuplicated(locus.ids)) {
      stop("Duplicate locus names after removing the file extension: ",
           paste(unique(locus.ids[duplicated(locus.ids)]), collapse = ", "))
    }

    writeFasta(sequences = as.list(sequences), names = locus.ids,
               reference.path, nbchar = 1000000, as.string = TRUE)
  } else {
    if (!file.copy(source.files, reference.path, overwrite = TRUE)) {
      stop("Could not copy the reference FASTA from ", source.files)
    }
  }

  # Indexes the reference as one complete owned set
  bwa = .toolCommand("bwa", bwa.path)
  samtools = .toolCommand("samtools", samtools.path)
  gatk = .toolCommand("gatk", gatk4.path)
  log = file.path(reference.dir, "reference-build.stderr.log")
  .runCommand(paste(bwa, "index -a bwtsw", shQuote(reference.path)),
              quiet, "BWA reference indexing", stderr.log = log)
  .runCommand(paste(samtools, "faidx", shQuote(reference.path)),
              quiet, "samtools reference indexing", stderr.log = log)
  .runCommand(paste(gatk, "CreateSequenceDictionary --REFERENCE", shQuote(reference.path),
                    "--OUTPUT", shQuote(dict.path),
                    "--USE_JDK_DEFLATER true --USE_JDK_INFLATER true"),
              quiet, "reference dictionary", stderr.log = log)

  if (!all(file.exists(reference.outputs))) {
    stop("Reference indexing did not produce every required file in ", reference.dir)
  }

  # Records the reference identity last, after every product is present
  current.record$reference = tools::md5sum(reference.path)
  saveRDS(current.record, record.path)
  invisible(reference.path)
}#end function

# END SCRIPT
