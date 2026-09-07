#' @title createContaminantDB
#'
#' @description Builds a local directory of contaminant reference genomes to be
#'   used by removeContamination(). Genomes are downloaded from NCBI with the
#'   NCBI Datasets API for assembly accessions (GCA_ or GCF_) and with Entrez
#'   efetch for nucleotide accessions. The function can also add the NCBI UniVec
#'   vector sequences, extra GenBank accessions, or a custom FASTA file. Each
#'   genome is downloaded only once, so an interrupted run continues where it
#'   stopped. Cached references that are no longer requested are retained but
#'   omitted from active-references.csv and the resulting contaminant index.
#'
#' @param decontamination.list path to a CSV file with at least two columns:
#'   Genome (a short display name) and GenBank_Accession (the accession to
#'   download). A column named Accession is also accepted.
#'
#' @param output.directory path to the directory where contaminant genome files
#'   will be saved.
#'
#' @param include.univec logical; if TRUE, the NCBI UniVec vector/adaptor
#'   sequence database is downloaded and added to the contaminant directory.
#'
#' @param include.genbank character vector of additional GenBank accessions to
#'   download and include; NULL skips this step.
#'
#' @param include.fasta path to an existing FASTA file to copy directly into
#'   the contaminant directory; NULL skips this step.
#'
#' @param overwrite logical; if TRUE an existing output.directory is deleted
#'   and every genome is downloaded again. If FALSE the function keeps the
#'   genomes that are already present and downloads only the missing ones.
#'   Downloads are published only after their FASTA content is verified.
#'
#' @return invisibly returns the paths of the contaminant files; side effect is
#'   a populated output.directory containing compressed genome FASTA files
#'   ready for use as a decontamination reference.
#'
#' @export

createContaminantDB = function(decontamination.list = NULL,
                               output.directory = "contaminant-references",
                               include.univec = TRUE,
                               include.genbank = NULL,
                               include.fasta = NULL,
                               overwrite = FALSE) {

  #Quick checks
  if (is.null(decontamination.list) == TRUE){ stop("Please provide list of sequences and genbank numbers.") }
  if (file.exists(decontamination.list) == F){ stop("Input list not found.") }
  if (length(include.univec) != 1 || is.logical(include.univec) == FALSE || is.na(include.univec) ||
      length(overwrite) != 1 || is.logical(overwrite) == FALSE || is.na(overwrite)) {
    stop("include.univec and overwrite must be TRUE or FALSE.")
  }

  sample.data = read.csv(file = decontamination.list, stringsAsFactors = FALSE)

  #The configuration file and the older tables use different column names
  if ("GenBank_Accession" %in% names(sample.data) == FALSE){
    if ("Accession" %in% names(sample.data) == TRUE){
      sample.data$GenBank_Accession = sample.data$Accession
    } else {
      stop("decontamination.list must have a GenBank_Accession or Accession column.")
    }
  }
  if ("Genome" %in% names(sample.data) == FALSE){
    stop("decontamination.list must have a Genome column.")
  }

  sample.data$Genome = trimws(as.character(sample.data$Genome))
  if (any(is.na(sample.data$Genome)) || any(nchar(sample.data$Genome) == 0)) {
    stop("Genome values cannot be missing or blank.")
  }

  # UniVec is downloaded separately and does not need an accession value.
  sample.data = sample.data[!tolower(sample.data$Genome) %in% "univec", , drop = FALSE]
  sample.data$GenBank_Accession = trimws(as.character(sample.data$GenBank_Accession))
  if (any(is.na(sample.data$GenBank_Accession)) ||
      any(nchar(sample.data$GenBank_Accession) == 0)) {
    stop("Accession values cannot be missing or blank.")
  }
  if (any(grepl("^https?://", sample.data$GenBank_Accession))) {
    stop("URL rows are not supported in decontamination.list. Use accessions or include.fasta.")
  }
  if (is.null(include.genbank) == FALSE) {
    include.genbank = trimws(as.character(include.genbank))
    if (any(is.na(include.genbank)) || any(nchar(include.genbank) == 0)) {
      stop("include.genbank cannot contain missing or blank accessions.")
    }
  }
  all.accessions = c(sample.data$GenBank_Accession, include.genbank)
  if (length(all.accessions) > 0 &&
      any(grepl("^[A-Za-z0-9_.-]+$", all.accessions) == FALSE)) {
    stop("Contaminant accessions may contain only letters, numbers, periods, underscores, and hyphens.")
  }
  .checkFileOutsideOutput(decontamination.list, output.directory)
  if (is.null(include.fasta) == FALSE) {
    if (file.exists(include.fasta) == F){ stop("include.fasta file not found.") }
    .checkFileOutsideOutput(include.fasta, output.directory)
  }

  #Sets up the output directory after every source has been validated
  if (dir.exists(output.directory) == F){
    dir.create(output.directory, recursive = TRUE)
  } else {
    if (overwrite == TRUE){ .resetDirectory(output.directory) }
  }#end else

  active.files = character()
  failed.accessions = character()

  if (is.null(include.fasta) != TRUE){
    fasta.extension = if (grepl("\\.gz$", include.fasta, ignore.case = TRUE)) ".fa.gz" else ".fa"
    manual.path = file.path(output.directory, paste0("manually-included-data", fasta.extension))
    temp.path = tempfile(pattern = "manual-reference-", tmpdir = output.directory,
                         fileext = fasta.extension)
    on.exit(unlink(temp.path), add = TRUE)
    if (file.copy(include.fasta, temp.path, overwrite = TRUE) == FALSE) {
      stop("Could not copy include.fasta into the contaminant directory.")
    }
    .checkFastaFile(temp.path)
    .publishFiles(temp.path, manual.path)
    active.files = c(active.files, manual.path)
  }#end if

  if (include.univec == TRUE) {
    univec.path = file.path(output.directory, "UniVec.fa")
    if (.contaminantFilePresent(univec.path) == FALSE){
      message("Downloading UniVec")
      temp.path = tempfile(pattern = "UniVec-", tmpdir = output.directory, fileext = ".fa")
      download.error = tryCatch({
        utils::download.file("https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec",
                             destfile = temp.path, quiet = TRUE)
        .checkFastaFile(temp.path)
        .publishFiles(temp.path, univec.path)
        NULL
      }, error = function(e) {
        unlink(temp.path)
        e
      })
      if (is.null(download.error) == FALSE) {
        failed.accessions = c(failed.accessions, "UniVec")
      }
    }
    if (.contaminantFilePresent(univec.path) == TRUE) {
      active.files = c(active.files, univec.path)
    }
  }

  if (is.null(include.genbank) != TRUE) {
    for (i in seq_along(include.genbank)) {
      out.path = file.path(output.directory, paste0("include-genbank_", include.genbank[i], ".fna.gz"))
      if (.contaminantFilePresent(out.path) == FALSE) {
        message("Downloading include.genbank: ", include.genbank[i])
        download.error = tryCatch(
          .downloadAccession(include.genbank[i], out.path),
          error = function(e) e
        )
        if (inherits(download.error, "error")) {
          failed.accessions = c(failed.accessions, include.genbank[i])
        }
      }
      if (.contaminantFilePresent(out.path) == TRUE) { active.files = c(active.files, out.path) }
    }
  }

  if (nrow(sample.data) > 0) {
    for (i in 1:nrow(sample.data)) {
      accession = trimws(sample.data$GenBank_Accession[i])
      out.path = file.path(output.directory,
                           paste0(.sanitizeName(sample.data$Genome[i]), "-", accession, ".fna.gz"))
      # Skips a genome that is already downloaded, so an interrupted run continues
      if (.contaminantFilePresent(out.path) == FALSE) {
        message("Downloading ", sample.data$Genome[i], " (", accession, ")")
        download.error = tryCatch(
          .downloadAccession(accession, out.path),
          error = function(e) e
        )
        if (inherits(download.error, "error")) {
          failed.accessions = c(failed.accessions, accession)
        }
      }
      if (.contaminantFilePresent(out.path) == TRUE) { active.files = c(active.files, out.path) }
    }
  }

  if (length(failed.accessions) > 0) {
    stop("Could not create the required contaminant references: ",
         paste(unique(failed.accessions), collapse = ", "), ".")
  }
  active.files = unique(active.files)
  if (length(active.files) == 0){
    stop("No contaminant reference files could be created in ", output.directory, ".")
  }

  active.table = data.frame(File = basename(active.files), stringsAsFactors = FALSE)
  active.path = file.path(output.directory, "active-references.csv")
  temp.active = tempfile(pattern = "active-references-", tmpdir = output.directory,
                         fileext = ".csv")
  write.csv(active.table, temp.active, row.names = FALSE)
  .publishFiles(temp.active, active.path)

  print("Creation of contamination database was successful.")

  return(invisible(active.files))
} #end function


# Internal helper: reports whether a contaminant file is already downloaded and
# holds data.
.contaminantFilePresent = function(file.path = NULL) {

  if (file.exists(file.path) == FALSE) { return(FALSE) }
  if (file.info(file.path)$size == 0) { return(FALSE) }
  return(tryCatch({ .checkFastaFile(file.path); TRUE }, error = function(e) FALSE))
}#end .contaminantFilePresent


# Internal helper: checks that a file holds FASTA sequences. NCBI answers a bad
# request with an error page and an HTTP 200 status, so the content must be
# checked. Without this check the error page became part of the reference.
.checkFastaFile = function(file.path = NULL) {

  if (grepl("\\.gz$", file.path) == TRUE) {
    .runCommand(paste0("gzip -t ", shQuote(file.path)), quiet = TRUE,
                task = "compressed FASTA check")
    first.line = .runCommandOutput(paste0("gzip -cd ", shQuote(file.path), " | head -1"),
                                   task = "FASTA check")
  } else {
    first.line = .runCommandOutput(paste0("head -1 ", shQuote(file.path)), task = "FASTA check")
  }

  if (length(first.line) == 0 || substr(first.line[1], 1, 1) != ">") {
    stop("The downloaded file is not in FASTA format: ", file.path)
  }

  return(invisible(TRUE))
}#end .checkFastaFile


# Internal helper: downloads a genome by accession using the NCBI Datasets API
# (GCA/GCF) or Entrez efetch (nucleotide accessions like NC_*), and saves it as
# a gzipped FASTA. A partial file is deleted so the next run downloads it again.
.downloadAccession = function(accession = NULL,
                              out.path = NULL) {

  temp.out = tempfile(pattern = paste0(basename(out.path), "-"),
                      tmpdir = dirname(out.path), fileext = ".fna.gz")
  on.exit(unlink(temp.out), add = TRUE)

  is.assembly = grepl("^GC[AF]_", accession)

  if (is.assembly) {
    zip.file = tempfile(fileext = ".zip")
    tmp.dir = tempfile()
    on.exit({ unlink(zip.file); unlink(tmp.dir, recursive = TRUE) }, add = TRUE)

    url = paste0(
      "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/",
      accession,
      "/download?include_annotation_type=GENOME_FASTA"
    )
    utils::download.file(url, destfile = zip.file, quiet = TRUE, mode = "wb")
    dir.create(tmp.dir)
    utils::unzip(zip.file, exdir = tmp.dir)
    fna.file = list.files(tmp.dir, pattern = "\\.fna$", recursive = TRUE, full.names = TRUE)
    if (length(fna.file) == 0) stop(paste("No .fna found in NCBI zip for", accession))
    .runCommand(paste0("gzip -c ", shQuote(fna.file[1]), " > ", shQuote(temp.out)),
                quiet = TRUE, task = "genome compression", keep.stdout = TRUE)
  } else {
    # Nucleotide accession (NC_*, AY_*, etc.) via Entrez efetch
    fa.file = tempfile(fileext = ".fa")
    on.exit(unlink(fa.file), add = TRUE)

    url = paste0(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
      "?db=nucleotide&id=", accession, "&rettype=fasta&retmode=text"
    )
    utils::download.file(url, destfile = fa.file, quiet = TRUE)
    .checkFastaFile(fa.file)
    .runCommand(paste0("gzip -c ", shQuote(fa.file), " > ", shQuote(temp.out)),
                quiet = TRUE, task = "genome compression", keep.stdout = TRUE)
  }

  .checkFastaFile(temp.out)
  .publishFiles(temp.out, out.path)

  return(invisible(out.path))
}#end .downloadAccession
