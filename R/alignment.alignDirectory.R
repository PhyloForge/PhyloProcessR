#' @title alignDirectory
#'
#' @description Aligns each FASTA file found in a directory using MAFFT (local pair
#' algorithm by default) and saves the resulting alignments in phylip format. Loci with
#' fewer than \code{min.taxa} sequences are skipped. A subset of loci can be processed
#' by specifying fractional start and end positions. Optionally renames sequences by
#' stripping everything after a pipe-delimited sample tag and removes the MAFFT reverse-
#' complement prefix from sequence names.
#'
#' @param marker.directory path to the directory containing unaligned FASTA files, one per locus.
#'
#' @param input.format format of the input sequence files. Currently "fasta" is supported.
#'
#' @param output.directory path to the directory where phylip alignments will be saved.
#' Default "alignments".
#'
#' @param min.taxa minimum number of sequences required for a locus to be aligned. Loci
#' with fewer sequences are skipped. Default 4.
#'
#' @param subset.start fractional position (0 to 1) in the list of loci at which to begin
#' processing. Useful for splitting a run across multiple jobs. Default 0.
#'
#' @param subset.end fractional position (0 to 1) in the list of loci at which to stop
#' processing. Default 1 (process all loci).
#'
#' @param adjust.direction logical. If TRUE, MAFFT adjusts sequence direction before
#' aligning. Default TRUE.
#'
#' @param remove.reverse.tag logical. If TRUE, the leading \code{_R_} tag added by MAFFT
#' to reverse-complemented sequences is stripped from sequence names. Default TRUE.
#'
#' @param sample.rename controls renaming of sequences before alignment. "none" leaves
#' names unchanged; "space" trims everything after the first space; "_|_" trims everything
#' before and including the pipe delimiter. Default "none".
#'
#' @param threads number of threads passed to MAFFT. Default 1.
#'
#' @param memory not currently used; reserved for future parallelisation. Default 1.
#'
#' @param overwrite logical. If TRUE, previously completed alignments are overwritten;
#' if FALSE, they are skipped. Default FALSE.
#'
#' @param quiet logical. If TRUE, suppresses MAFFT screen output. Default TRUE.
#'
#' @param mafft.path path to the directory containing the MAFFT executable. If NULL, MAFFT
#' is expected to be on the system PATH.
#'
#' @return Writes phylip alignment files to \code{output.directory}. No value is returned to R.
#'
#' @export

alignDirectory = function(marker.directory = NULL,
                          input.format = "fasta",
                          output.directory = "alignments",
                          min.taxa = 4,
                          subset.start = 0,
                          subset.end = 1,
                          adjust.direction = TRUE,
                          remove.reverse.tag = TRUE,
                          sample.rename = c("none", "space", "_|_"),
                          threads = 1,
                          memory = 1,
                          overwrite = FALSE,
                          quiet = TRUE,
                          mafft.path = NULL) {
  #Debug setup
  # setwd("/Users/chutter/Dropbox/Elapid_probe_design/Sanger_venom_genes_from_NCBI")
  # marker.directory = "/Users/chutter/Dropbox/Elapid_probe_design/Sanger_venom_genes_from_NCBI/01_all_sanger_loci_clustered_(cd_hit_80)"
  # output.directory = "alignments"
  # input.format = "fasta"
  #
  # #Main settings
  # subset.start = 0
  # subset.end = 1
  # min.taxa = 2
  # threads = 1
  # memory = 2
  # overwrite = TRUE
  # quiet = TRUE
  # remove.reverse.tag = TRUE
  # adjust.direction = TRUE
  #
  #program paths
  #mafft.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin"

  #Same adds to bbmap path
  if (is.null(mafft.path) == FALSE){
    b.string = unlist(strsplit(mafft.path, ""))
    if (b.string[length(b.string)] != "/") {
      mafft.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { mafft.path = NULL }

  #Initial checks
  if (is.null(marker.directory) == T){ stop("A fasta file of targets is needed for alignment.") }

  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } else { dir.create(output.directory) }

  #Load sample file from probe matching step
  locus.names = list.files(marker.directory)

  #Checks for alignments already done and removes from the to-do list
  if (overwrite == FALSE){
    done = list.files(output.directory)
    locus.names = locus.names[!locus.names %in% gsub(".phy$", "", done)]
  }
  if (length(locus.names) == 0){ return(invisible(NULL)) }

  #Figures out start and end of subset
  sub.start = floor(subset.start * length(locus.names))
  if (sub.start == 0){ sub.start = 1}
  sub.end = floor(subset.end * length(locus.names))

  #Loops through each locus and writes each species to end of file
  for (i in sub.start:sub.end) {

    #Match probe names to contig names to acquire data
    match.data = Biostrings::readDNAStringSet(paste0(marker.directory, "/", locus.names[i]),format = "fasta")   # loads up fasta file

    ##############
    #STEP 1: Throw out loci if there are too few taxa
    ##############
    if (length(names(match.data)) < min.taxa){
      print(paste0(locus.names[i], " had too few taxa"))
      next
    }

    #Rename samples based on pipe
    if (sample.rename == "_|_"){
      names(match.data) = gsub(pattern = ".*_\\|_", replacement = "", x = names(match.data))
    }

    #Rename samples here deleting everything after the first space
    if (sample.rename == "space"){
      names(match.data) = gsub("\\ .*", "", names(match.data))
    }

    #Gets reference locus
    final.loci = match.data

    #Aligns and then reverses back to correction orientation
    alignment = runMafft(sequence.data = match.data,
                         save.name = paste0(output.directory, "/", locus.names[i]),
                         algorithm = "localpair",
                         adjust.direction = adjust.direction,
                         threads = threads,
                         cleanup.files = T,
                         quiet = quiet,
                         mafft.path = mafft.path)

    #Checks for failed mafft run
    if (length(alignment) == 0){
      print(paste0(locus.names[i], " did not successfully align."))
      next }

    #Removes the reverse name
    if (remove.reverse.tag == TRUE){
      names(alignment) = gsub(pattern = "^_R_", replacement = "", x = names(alignment))
    }#end if

    #Saves prelim exon file
    new.align = strsplit(as.character(alignment), "")
    aligned.set = as.matrix(ape::as.DNAbin(new.align))

    #readies for saving
    writePhylip(alignment = aligned.set,
                file=paste0(output.directory, "/", locus.names[i], ".phy"),
                interleave = F,
                strict = F)

    #Deletes old files
    print(paste0(locus.names[i], " alignment saved."))

  }# end big i loop

}# end function


#END SCRIPT

