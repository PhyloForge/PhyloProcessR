#' @title removeAdaptors
#'
#' @description Removes Illumina adapter sequences from raw paired-end reads
#'   using fastp with automatic adapter detection (--detect_adapter_for_pe).
#'   Also applies a minimum length filter (30 bp), low-complexity filtering,
#'   and poly-X tail trimming. Quality-based filtering is not applied here.
#'
#' @param input.reads path to a directory of raw paired-end reads in fastq.gz
#'   format. Samples may be in per-sample sub-directories or identified by a
#'   shared filename prefix.
#'
#' @param output.directory path to the directory where adapter-trimmed reads
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
#' @return invisibly returns the summary data frame; writes adapter-trimmed
#'   fastq.gz files to output.directory, the fastp HTML and JSON reports to
#'   logs/sample_logs, and a CSV summary to logs/removeAdaptors_summary.csv.
#'
#' @export

removeAdaptors = function(input.reads = NULL,
                          output.directory = NULL,
                          fastp.path = NULL,
                          threads = 1,
                          mem = 8,
                          overwrite = FALSE,
                          quiet = TRUE) {

  .runFastpStep(input.reads = input.reads,
                output.directory = output.directory,
                fastp.path = fastp.path,
                fastp.args = paste0("--length_required 30 --low_complexity_filter --complexity_threshold 30",
                                    " --trim_poly_x --detect_adapter_for_pe --compression 6"),
                task = "trim-adaptors+filter-complex",
                report.tag = "adapter-trim_fastp",
                summary.csv = "logs/removeAdaptors_summary.csv",
                threads = threads,
                mem = mem,
                overwrite = overwrite,
                quiet = quiet)

}#end function
