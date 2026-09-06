#' @title createContaminantDB
#'
#' @description Builds a local directory of contaminant reference genomes to be
#'   used by removeContamination(). Genomes are downloaded from NCBI with the
#'   NCBI Datasets API for assembly accessions (GCA_ or GCF_) and with Entrez
#'   efetch for nucleotide accessions. The function can also add the NCBI UniVec
#'   vector sequences, extra GenBank accessions, or a custom FASTA file. Each
#'   genome is downloaded only once, so an interrupted run continues where it
#'   stopped.
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

  #Sets up the output directory
  if (dir.exists(output.directory) == F){
    dir.create(output.directory, recursive = TRUE)
  } else {
    if (overwrite == TRUE){ .resetDirectory(output.directory) }
  }#end else

  sample.data = read.csv(file = decontamination.list)

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

  #UniVec is downloaded separately, so its row is dropped here
  sample.data = sample.data[!tolower(sample.data$Genome) %in% "univec", , drop = FALSE]

  if (is.null(include.fasta) != TRUE){
    if (file.exists(include.fasta) == F){ stop("include.fasta file not found.") }
    file.copy(include.fasta, paste0(output.directory, "/manually-included-data.fa"), overwrite = TRUE)
  }#end if

  if (include.univec == TRUE) {
    univec.path = file.path(output.directory, "UniVec.fa")
    if (.contaminantFilePresent(univec.path) == FALSE){
      message("Downloading UniVec")
      tryCatch({
        utils::download.file("https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec",
                             destfile = univec.path, quiet = TRUE)
        .checkFastaFile(univec.path)
      }, error = function(e) {
        unlink(univec.path)
        warning("Could not download UniVec: ", conditionMessage(e))
      })
    }
  }

  if (is.null(include.genbank) != TRUE) {
    for (i in seq_along(include.genbank)) {
      out.path = file.path(output.directory, paste0("include-genbank_", include.genbank[i], ".fna.gz"))
      if (.contaminantFilePresent(out.path) == TRUE) { next }
      message("Downloading include.genbank: ", include.genbank[i])
      tryCatch(
        .downloadAccession(include.genbank[i], out.path),
        error = function(e) warning("Could not download ", include.genbank[i], ": ", conditionMessage(e))
      )
    }
  }

  if (nrow(sample.data) > 0) {
    for (i in 1:nrow(sample.data)) {
      accession = trimws(sample.data$GenBank_Accession[i])
      # Skip rows whose accession column contains a URL (handled elsewhere)
      if (grepl("^https?://", accession)) next
      out.path = file.path(output.directory,
                           paste0(.sanitizeName(sample.data$Genome[i]), "-", accession, ".fna.gz"))
      # Skips a genome that is already downloaded, so an interrupted run continues
      if (.contaminantFilePresent(out.path) == TRUE) { next }
      message("Downloading ", sample.data$Genome[i], " (", accession, ")")
      tryCatch(
        .downloadAccession(accession, out.path),
        error = function(e) warning("Could not download ", accession, ": ", conditionMessage(e))
      )
    }
  }

  final.files = list.files(output.directory, full.names = TRUE)
  final.files = final.files[grep("\\.fna\\.gz$|\\.fa\\.gz$|\\.fasta\\.gz$|\\.fna$|\\.fa$|\\.fasta$", final.files)]

  if (length(final.files) == 0){
    stop("No contaminant reference files could be created in ", output.directory, ".")
  }

  print("Creation of contamination database was successful.")

  return(invisible(final.files))
} #end function


# Internal helper: reports whether a contaminant file is already downloaded and
# holds data.
.contaminantFilePresent = function(file.path = NULL) {

  if (file.exists(file.path) == FALSE) { return(FALSE) }
  return(file.info(file.path)$size > 0)
}#end .contaminantFilePresent


# Internal helper: checks that a file holds FASTA sequences. NCBI answers a bad
# request with an error page and an HTTP 200 status, so the content must be
# checked. Without this check the error page became part of the reference.
.checkFastaFile = function(file.path = NULL) {

  if (grepl("\\.gz$", file.path) == TRUE) {
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

  on.exit({
    if (file.exists(out.path) == TRUE && file.info(out.path)$size == 0) { unlink(out.path) }
  }, add = TRUE)

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
    .runCommand(paste0("gzip -c ", shQuote(fna.file[1]), " > ", shQuote(out.path)),
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
    .runCommand(paste0("gzip -c ", shQuote(fa.file), " > ", shQuote(out.path)),
                quiet = TRUE, task = "genome compression", keep.stdout = TRUE)
  }

  .checkFastaFile(out.path)

  return(invisible(out.path))
}#end .downloadAccession
