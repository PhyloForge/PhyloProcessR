#' @title removeDuplicateReads
#'
#' @description Removes PCR and optical duplicate read pairs from processed
#'   paired-end fastq files using fastp's --dedup mode. Adapter trimming,
#'   quality filtering, and length filtering are disabled so that only
#'   duplicate removal is performed.
#'
#' @param input.reads path to a directory of processed paired-end reads in
#'   fastq.gz format, organised in per-sample sub-directories.
#'
#' @param output.directory path to the directory where deduplicated reads will
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
#' @return invisibly returns the summary data frame; writes deduplicated
#'   fastq.gz files to output.directory, the fastp HTML and JSON reports to
#'   logs/sample_logs, and a CSV summary to
#'   logs/removeDuplicateReads_summary.csv.
#'
#' @export

removeDuplicateReads = function(input.reads = NULL,
                                output.directory = "deduped-reads",
                                fastp.path = NULL,
                                threads = 1,
                                mem = 1,
                                overwrite = FALSE,
                                quiet = FALSE) {

  .runFastpStep(input.reads = input.reads,
                output.directory = output.directory,
                fastp.path = fastp.path,
                fastp.args = paste0("--dedup --dup_calc_accuracy 6",
                                    " --disable_adapter_trimming --disable_quality_filtering",
                                    " --disable_length_filtering --compression 6"),
                task = "deduplicate-reads",
                report.tag = "dedup_fastp",
                summary.csv = "logs/removeDuplicateReads_summary.csv",
                threads = threads,
                mem = mem,
                overwrite = overwrite,
                quiet = quiet)

}#end function
