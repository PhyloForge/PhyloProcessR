#' @title qualityTrimReads
#'
#' @description Performs sliding-window quality trimming on paired-end fastq
#'   reads using fastp's --cut_front and --cut_tail options. Adapter trimming
#'   and length filtering are disabled so that only low-quality bases at the
#'   ends of reads are removed. Intended to be run after adapter removal and
#'   deduplication.
#'
#' @param input.reads path to a directory of processed paired-end reads in
#'   fastq.gz format, organised in per-sample sub-directories.
#'
#' @param output.directory path to the directory where quality-trimmed reads
#'   will be saved (one sub-directory per sample).
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
#' @return invisibly returns the summary data frame; writes quality-trimmed
#'   fastq.gz files to output.directory, the fastp HTML and JSON reports to
#'   logs/sample_logs, and a CSV summary to logs/qualityTrimReads_summary.csv.
#'
#' @export

qualityTrimReads = function(input.reads = "deduped-reads",
                            output.directory = "quality-trimmed-reads",
                            fastp.path = NULL,
                            threads = 1,
                            mem = 1,
                            overwrite = FALSE,
                            quiet = TRUE) {

  .runFastpStep(input.reads = input.reads,
                output.directory = output.directory,
                fastp.path = fastp.path,
                fastp.args = paste0("--cut_front --cut_tail",
                                    " --disable_adapter_trimming --disable_length_filtering",
                                    " --compression 6"),
                task = "quality-trim",
                report.tag = "quality-trim_fastp",
                summary.csv = "logs/qualityTrimReads_summary.csv",
                threads = threads,
                mem = mem,
                overwrite = overwrite,
                quiet = quiet)

}#end function
