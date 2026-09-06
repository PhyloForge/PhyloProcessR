#' @title fastqStats
#'
#' @description Counts reads in fastq files and computes summary statistics for
#'   each sample lane, including per-read-file counts, total reads, read pairs,
#'   megabase pairs sequenced, and millions of reads. Compressed and
#'   uncompressed fastq files are both accepted. Results are written to a CSV
#'   file and returned as a data frame.
#'
#' @param read.directory path to the top-level directory containing sample
#'   read files or sample sub-directories.
#'
#' @param sub.directory optional sub-directory name within each sample folder
#'   to restrict which reads are counted (e.g. "cleaned-reads-snp"); NULL uses
#'   the top-level directory structure.
#'
#' @param output.name base name (without extension) for the output CSV file.
#'
#' @param read.length expected read length in base pairs, used to calculate
#'   megabase pairs (MegaBasePairs = read.length * read.pairs * 2 / 1e6).
#'
#' @param threads number of read files to count at the same time.
#'
#' @param mem amount of RAM in GB (currently reserved).
#'
#' @param overwrite logical; if TRUE an existing output CSV is deleted before
#'   writing.
#'
#' @return a data frame with columns Sample, Lane, Read1_Count, Read2_Count,
#'   Read3_Count, Total_Reads, Read_Pairs, Read_Length, MegaBasePairs, and
#'   Reads_Per_Million. Total_Reads counts the two paired mates and excludes
#'   Read3_Count, which holds merged reads or singletons. Reads_Per_Million is
#'   the total read count in millions, including Read3_Count, and is used as a
#'   scaling factor. Results are also written to output.name.csv.
#'
#' @export

fastqStats = function(read.directory = NULL,
                      sub.directory = NULL,
                      output.name = "fastq-stats",
                      read.length = 150,
                      threads = 1,
                      mem = 1,
                      overwrite = FALSE) {

  #Quick checks
  if (is.null(read.directory) == TRUE){ stop("Please provide input reads.") }
  if (file.exists(read.directory) == F){ stop("Input reads not found.") }

  #Removes a previous run of the same output file
  if (file.exists(paste0(output.name, ".csv")) == T){
    if (overwrite == TRUE){ unlink(paste0(output.name, ".csv")) }
  }#end if

  read.directory = sub("/+$", "", read.directory)
  reads = list.files(read.directory, recursive = T, full.names = T)
  reads = reads[grep("\\.fastq\\.gz$|\\.fq\\.gz$|\\.fastq$|\\.fq$", reads)]
  read.names = .relativePaths(reads, read.directory)

  if (is.null(sub.directory) != TRUE) {
    keep = grepl(paste0("/", sub.directory, "/"), reads, fixed = TRUE)
    reads = reads[keep]
    read.names = read.names[keep]
    sample.names = unique(gsub(paste0("/", sub.directory, "/.*"), "", read.names))
  } else {
    sample.names = .listSampleNames(read.directory)
  }

  if (length(sample.names) == 0){ return("no samples remain to analyze.") }

  #################################################
  ### Part A: build the list of lanes to count
  #################################################
  lane.list = list()

  for (i in seq_along(sample.names)) {

    sample.reads = .matchPrefix(reads, read.names, sample.names[i])

    #Returns an error if reads are not found
    if (length(sample.reads) == 0 ){
      warning(sample.names[i], " does not have any reads present. Skipping.")
      next
    } #end if statement

    lane.prefixes = .stripReadSuffix(sample.reads)

    for (j in seq_along(lane.prefixes)){

      lane.reads = sort(.matchPrefix(reads, reads, lane.prefixes[j]))

      #Returns an error if reads are not found
      if (length(lane.reads) == 0 ){
        warning(lane.prefixes[j], " does not have any reads present. Skipping.")
        next
      } #end if statement

      lane.list[[length(lane.list) + 1]] = list(sample = sample.names[i],
                                                lane = gsub(".*_", "", basename(lane.prefixes[j])),
                                                files = lane.reads)
    }#end lane j loop

  }#end sample i loop

  if (length(lane.list) == 0){ return("no samples remain to analyze.") }

  #################################################
  ### Part B: count the reads in parallel
  #################################################
  # Counting decompresses each file once, so the lanes are counted at the same
  # time on the requested number of threads.
  count.list = parallel::mclapply(lane.list, mc.cores = max(1, threads), FUN = function(lane.job) {

    read1.count = .countFastqReads(lane.job$files[1])
    read2.count = if (length(lane.job$files) >= 2) .countFastqReads(lane.job$files[2]) else 0
    read3.count = if (length(lane.job$files) >= 3) .countFastqReads(lane.job$files[3]) else 0

    scale.factor = (read1.count + read2.count + read3.count) / 1000000
    if (isTRUE(read1.count == read2.count)){ read.pairs = read1.count } else { read.pairs = NA_real_ }

    data.frame(Sample = lane.job$sample,
               Lane = lane.job$lane,
               Read1_Count = read1.count,
               Read2_Count = read2.count,
               Read3_Count = read3.count,
               Total_Reads = read1.count + read2.count,
               Read_Pairs = read.pairs,
               Read_Length = read.length,
               MegaBasePairs = (read.length * read.pairs * 2) / 1000000,
               Reads_Per_Million = scale.factor,
               stringsAsFactors = FALSE)
  })

  failed = vapply(count.list, function(x) inherits(x, "try-error"), logical(1))
  if (any(failed) == TRUE){
    stop("Read counting failed for ", sum(failed), " lane(s). First error: ",
         as.character(count.list[failed][[1]]))
  }

  summary.data = do.call(rbind, count.list)

  for (sample.name in unique(summary.data$Sample)) {
    print(paste0(sample.name, " Completed fastq counting!"))
  }

  write.csv(summary.data, file = paste0(output.name, ".csv"), row.names = FALSE)
  return(summary.data)
}#end function
