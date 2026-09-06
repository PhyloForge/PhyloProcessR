#' @title mergePairedEndReads
#'
#' @description Merges overlapping paired-end reads using fastp's --merge mode.
#'   Read pairs that overlap are combined into a single merged read (saved as
#'   READ3/merged), while non-overlapping pairs are kept as READ1 and READ2.
#'   Adapter trimming, quality filtering, and length filtering are disabled so
#'   that only the merging step is performed.
#'
#' @param input.reads path to a directory of processed paired-end reads in
#'   fastq.gz format, organised in per-sample sub-directories.
#'
#' @param output.directory path to the directory where merged read files will
#'   be saved (one sub-directory per sample).
#'
#' @param fastp.path system path to the directory that contains the fastp
#'   executable, or the full path to the executable; NULL searches the system
#'   PATH.
#'
#' @param threads number of CPU threads to pass to fastp.
#'
#' @param mem amount of RAM in GB (currently reserved).
#'
#' @param overwrite logical; if TRUE the output directory is deleted and
#'   recreated before processing. FALSE resumes and skips only the lanes that
#'   already have complete output.
#'
#' @param quiet logical; if TRUE fastp stdout and stderr are suppressed.
#'
#' @return invisibly returns the summary data frame; writes merged fastq.gz
#'   files to output.directory, the fastp HTML and JSON reports to
#'   logs/sample_logs, and a CSV summary to
#'   logs/mergePairedEndReads_summary.csv. The summary holds a mergedReads
#'   column with the number of merged single reads per lane.
#'
#' @export

mergePairedEndReads = function(input.reads = NULL,
                               output.directory = "read-processing/pe-merged-reads",
                               fastp.path = NULL,
                               threads = 1,
                               mem = 8,
                               overwrite = FALSE,
                               quiet = TRUE) {

  .runFastpStep(input.reads = input.reads,
                output.directory = output.directory,
                fastp.path = fastp.path,
                fastp.args = paste0("--merge",
                                    " --disable_adapter_trimming --disable_quality_filtering",
                                    " --disable_length_filtering --compression 6"),
                task = "merge-pe-reads",
                report.tag = "pe-merged_fastp",
                summary.csv = "logs/mergePairedEndReads_summary.csv",
                merge.reads = TRUE,
                threads = threads,
                mem = mem,
                overwrite = overwrite,
                quiet = quiet)

}#end function
