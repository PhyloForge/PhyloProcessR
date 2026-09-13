#' @title combineMarkerAlignments
#'
#' @description Combines two or more directories of per-marker alignments into a single
#' output directory. Every alignment file is copied under its own name, so the result holds
#' one alignment per marker across all sources. A marker that occurs in more than one input
#' directory is an error, because the sources would otherwise overwrite each other; this
#' keeps, for example, ordinary target-marker alignments and novel-marker alignments
#' distinct when they are merged for downstream trimming.
#'
#' @param alignment.directories character vector of paths to the directories of alignment
#' files to combine. Each directory holds one alignment file per marker.
#'
#' @param output.directory path to the directory where the combined alignment files will be
#' saved.
#'
#' @param overwrite logical. If TRUE, the output directory is removed and recreated before
#' copying; if FALSE, files are added to an existing directory. Default FALSE.
#'
#' @return Copies alignment files to \code{output.directory} and invisibly returns the path.
#'
#' @export

combineMarkerAlignments = function(alignment.directories = NULL,
                                   output.directory = NULL,
                                   overwrite = FALSE
                                   ) {

  # Parameter checks
  if (is.null(alignment.directories) == TRUE) {
    stop("Error: alignment.directories is required.")
  }
  if (is.null(output.directory) == TRUE) {
    stop("Error: output.directory is required.")
  }

  missing.dirs = alignment.directories[!dir.exists(alignment.directories)]
  if (length(missing.dirs) > 0) {
    stop("alignment.directories not found: ", paste(missing.dirs, collapse = ", "), ".")
  }

  # List the alignments in each source and reject markers that occur in more than
  # one source, which would otherwise overwrite each other in the output.
  source.files = lapply(alignment.directories, .alignmentFiles, full.names = TRUE)
  all.files = unlist(source.files)
  if (length(all.files) == 0) {
    stop("No alignment files were found in any of the alignment.directories.")
  }
  duplicate.ids = unique(.alignmentId(all.files)[duplicated(.alignmentId(all.files))])
  if (length(duplicate.ids) > 0) {
    stop("Marker(s) occur in more than one input directory: ",
         paste(utils::head(duplicate.ids, 5), collapse = ", "), ".")
  }

  # Checks output overwrite after validating all inputs.
  if (overwrite == TRUE) {
    if (dir.exists(output.directory) == TRUE) unlink(output.directory, recursive = TRUE)
    dir.create(output.directory, recursive = TRUE)
  } else {
    if (!dir.exists(output.directory)) dir.create(output.directory, recursive = TRUE)
  }#end overwrite if

  for (source in all.files) {
    dest = file.path(output.directory, basename(source))
    if (overwrite == FALSE && file.exists(dest)) { next }
    .copyAlignment(source, dest, overwrite = overwrite)
  }

  invisible(output.directory)

}#end function
