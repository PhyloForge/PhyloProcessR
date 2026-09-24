#' @title copyUnmatchedMarkers
#'
#' @description Copies alignments whose marker name is absent from a gene-metadata file into
#' an existing output directory. \code{gatherUnlinked} only keeps markers listed in the gene
#' metadata, so markers that carry no gene assignment (for example newly discovered novel
#' markers) would be dropped from the unlinked set. This function adds those unmatched
#' markers back so every locus is represented once.
#'
#' @param alignment.directory path to the directory of per-marker alignment files to check.
#'
#' @param output.directory path to the directory that unmatched alignments are added to,
#' typically the unlinked set produced by \code{gatherUnlinked}.
#'
#' @param feature.gene.names path to a tab-delimited metadata file with at minimum columns
#' named \code{marker} and \code{gene}. Markers absent from its \code{marker} column are the
#' ones copied.
#'
#' @param overwrite logical. If TRUE, an existing alignment of the same name is replaced;
#' if FALSE, it is left in place. Default FALSE.
#'
#' @return Copies alignment files to \code{output.directory} and invisibly returns the
#' number of alignments added.
#'
#' @export

copyUnmatchedMarkers = function(alignment.directory = NULL,
                                output.directory = NULL,
                                feature.gene.names = NULL,
                                overwrite = FALSE
                                ) {

  # Parameter checks
  if (is.null(alignment.directory) == TRUE) {
    stop("Error: alignment.directory is required.")
  }
  if (is.null(output.directory) == TRUE) {
    stop("Error: output.directory is required.")
  }
  if (is.null(feature.gene.names) == TRUE) {
    stop("Error: a table associating each marker with a gene is needed.")
  }
  if (dir.exists(alignment.directory) == FALSE) {
    stop("alignment.directory not found. Please check the path.")
  }
  if (!dir.exists(output.directory)) dir.create(output.directory, recursive = TRUE)

  metadata = .readGeneMetadata(feature.gene.names)
  align.files = .alignmentFiles(alignment.directory, full.names = TRUE)
  unmatched = align.files[!.alignmentId(align.files) %in% metadata$marker]

  copied = .copyAlignmentsWithLogs(
    unmatched, output.directory, overwrite,
    file.path("logs", "unmatched_marker_logs", basename(output.directory)),
    "_unmatched_marker.log"
  )

  message(sum(copied), " unmatched marker(s) added individually.")
  invisible(sum(copied))

}#end function
