#' @title organizeReads
#'
#' @description Copies raw fastq files from a source directory into a
#'   standardised per-sample directory structure, renaming files to a
#'   consistent convention (SampleName_L00N_READ1/2.fastq.gz). Sample-to-file
#'   mapping is provided by a two-column CSV (File, Sample) such as that
#'   produced by makeFileRename(). Each row in the CSV represents one
#'   paired-end lane for a sample; multiple rows for the same sample generate
#'   multiple lanes.
#'
#' @param read.directory path to the directory containing the raw fastq or
#'   fastq.gz files to be organised.
#'
#' @param output.directory path to the directory where the reorganised reads
#'   will be saved. A sub-directory is created for each unique sample name.
#'
#' @param rename.file path to a CSV file with at least two columns: File
#'   (partial or full file name used to locate the reads) and Sample (the
#'   desired output sample name).
#'
#' @param link.reads logical; if TRUE the reads are hard-linked instead of
#'   copied, and a symbolic link is used when a hard link is not possible. This
#'   avoids a second copy of a read set that is often hundreds of gigabytes.
#'   FALSE copies the files. Plain .fastq and .fq inputs are compressed once
#'   because organized output always uses the .fastq.gz convention.
#'
#' @param overwrite logical; if TRUE the output directory is deleted and
#'   recreated. Completed lanes are skipped when overwrite is FALSE.
#'
#' @return invisibly; side effect is a populated output.directory with one
#'   sub-directory per sample containing renamed fastq.gz read pairs.
#'
#' @export

organizeReads = function(read.directory = NULL,
                        output.directory = "organized-reads",
                        rename.file = NULL,
                        link.reads = FALSE,
                        overwrite = FALSE) {

  #Quick checks
  options(stringsAsFactors = FALSE)
  if (is.null(read.directory) == TRUE){ stop("Please provide a directory of raw reads.") }
  if (file.exists(read.directory) == F){ stop("Input reads not found.") }
  if (is.null(rename.file) == TRUE){ stop("Please provide a table of file to sample name conversions.") }
  if (file.exists(rename.file) == F){ stop("Rename file not found.") }
  if (length(link.reads) != 1 || is.logical(link.reads) == FALSE || is.na(link.reads)) {
    stop("link.reads must be TRUE or FALSE.")
  }
  if (length(overwrite) != 1 || is.logical(overwrite) == FALSE || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE.")
  }

  sample.data = read.csv(rename.file, stringsAsFactors = FALSE)
  sample.data = .validateRenameTable(sample.data)
  if (nrow(sample.data) == 0){ return("no samples available to organize.") }
  .checkDirectoryOverlap(read.directory, output.directory)
  .checkFileOutsideOutput(rename.file, output.directory)

  #Sets directory and reads in
  if (dir.exists(output.directory) == F){
    dir.create(output.directory, recursive = TRUE)
  } else {
    if (overwrite == TRUE){ .resetDirectory(output.directory) }
  }#end else

  #Read in sample data and finds reads
  read.directory = sub("/+$", "", read.directory)
  reads = .listFastqFiles(read.directory)
  read.names = .relativePaths(reads, read.directory)
  if (dir.exists("logs/sample_logs") == FALSE) { dir.create("logs/sample_logs", recursive = TRUE) }

  sample.names = unique(sample.data$Sample)

  for (i in seq_along(sample.names)){

    temp.data = sample.data[sample.data$Sample %in% sample.names[i], ]

    for (j in 1:nrow(temp.data)) {
      #################################################
      ### Part A: prepare for loading and checks
      #################################################
      # Sets up the output paths first so a finished lane can be skipped
      out.path = paste0(output.directory, "/", temp.data$Sample[j])
      lane.tag = sprintf("L%03d", j)
      outread.1 = paste0(out.path, "/", temp.data$Sample[j], "_", lane.tag, "_READ1.fastq.gz")
      outread.2 = paste0(out.path, "/", temp.data$Sample[j], "_", lane.tag, "_READ2.fastq.gz")

      # Finds all files for this given sample
      sample.reads = .matchPrefix(reads, read.names, temp.data$File[j])
      if (length(sample.reads) == 0) {
        sample.reads = reads[grepl(temp.data$File[j], read.names, fixed = TRUE)]
      }
      # Checks the Sample column in case already renamed
      if (length(sample.reads) == 0) {
        sample.reads = .matchPrefix(reads, read.names, temp.data$Sample[j])
      }

      # Returns an error if reads are not found
      if (length(sample.reads) == 0) {
        stop(paste0(
          temp.data$Sample[j], " does not have any reads present for files ",
          temp.data$File[j], " from the input spreadsheet."
        ))
      } # end if statement

      sample.reads = .orderReadPair(sample.reads)
      metadata.file = file.path("logs/sample_logs", temp.data$Sample[j],
                                paste0(temp.data$Sample[j], "_", lane.tag,
                                       "_organization-metadata.csv"))
      metadata = .laneMetadata(sample.reads,
                               list(file.match = temp.data$File[j], sample = temp.data$Sample[j]))

      if (overwrite == FALSE &&
          .laneComplete(c(outread.1, outread.2), metadata.file = metadata.file,
                        metadata = metadata) == TRUE) { next }
      if (overwrite == FALSE && .metadataConflicts(metadata.file, metadata) == TRUE) {
        stop(temp.data$Sample[j], " ", lane.tag,
             " was organized from different input reads. Use overwrite = TRUE to replace it.")
      }

      #################################################
      ### Part B: Create directories and move files
      #################################################
      # Create sample directory
      if (file.exists(out.path) == FALSE) {
        dir.create(out.path, recursive = TRUE)
      }

      .linkOrCopyRead(sample.reads[1], outread.1, link.reads)
      .linkOrCopyRead(sample.reads[2], outread.2, link.reads)
      .writeLaneMetadata(metadata, metadata.file)
    } # end j loop

  }#end i loop

}#end function


# Internal helper: places one read file at the output path. A hard link costs no
# disk space. The function falls back to a symbolic link and then to a copy, for
# example when the input is on a different file system.
.linkOrCopyRead = function(source.file = NULL,
                           target.file = NULL,
                           link.reads = FALSE) {

  temp.file = tempfile(pattern = paste0(basename(target.file), "-"),
                       tmpdir = dirname(target.file), fileext = ".fastq.gz")
  on.exit(unlink(temp.file), add = TRUE)

  compressed = grepl("\\.gz$", source.file, ignore.case = TRUE)

  if (link.reads == TRUE && compressed == TRUE) {
    if (isTRUE(file.link(source.file, temp.file)) == FALSE &&
        isTRUE(file.symlink(normalizePath(source.file), temp.file)) == FALSE &&
        file.copy(source.file, temp.file, overwrite = TRUE) == FALSE) {
      stop("Could not place ", source.file, " at ", target.file, ".")
    }
  } else if (compressed == TRUE) {
    if (file.copy(source.file, temp.file, overwrite = TRUE) == FALSE) {
      stop("Could not place ", source.file, " at ", target.file, ".")
    }
  } else {
    .runCommand(paste0("gzip -c ", shQuote(source.file), " > ", shQuote(temp.file)),
                quiet = TRUE, task = "FASTQ compression", keep.stdout = TRUE)
  }

  .publishFiles(temp.file, target.file)

  return(invisible(TRUE))
}#end .linkOrCopyRead
