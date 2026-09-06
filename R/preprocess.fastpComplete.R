#' @title fastpComplete
#'
#' @description Runs a comprehensive fastp cleaning step on raw paired-end
#'   Illumina reads. In a single pass, performs adapter detection and trimming,
#'   length filtering (minimum 60 bp), low-complexity filtering, poly-X
#'   trimming, overlap-based base correction, and deduplication. Intended as an
#'   all-in-one alternative to running the individual preprocess steps
#'   separately.
#'
#' @param input.reads path to a directory of raw reads. Each sample may be
#'   organised in its own sub-directory or identified by a shared filename
#'   prefix containing "_L00".
#'
#' @param output.directory path to the directory where cleaned reads will be
#'   saved (one sub-directory per sample).
#'
#' @param fastp.path system path to the directory that contains the fastp
#'   executable, or the full path to the executable; NULL searches the system
#'   PATH.
#'
#' @param threads number of CPU threads to pass to fastp.
#'
#' @param mem amount of RAM in GB (currently unused by fastp directly, reserved
#'   for future use).
#'
#' @param overwrite logical; if TRUE the output directory is deleted and
#'   recreated before processing. FALSE resumes and skips only the lanes that
#'   already have complete output.
#'
#' @param quiet logical; if TRUE fastp stdout and stderr are suppressed.
#'
#' @return invisibly returns the summary data frame; writes cleaned fastq.gz
#'   files to output.directory, the fastp HTML and JSON reports to
#'   logs/sample_logs, and a CSV summary to logs/fastpComplete_summary.csv.
#'
#' @export

fastpComplete = function(input.reads = NULL,
                         output.directory = NULL,
                         fastp.path = NULL,
                         threads = 1,
                         mem = 8,
                         overwrite = FALSE,
                         quiet = TRUE) {

  .runFastpStep(input.reads = input.reads,
                output.directory = output.directory,
                fastp.path = fastp.path,
                fastp.args = paste0("--length_required 60 --low_complexity_filter --complexity_threshold 30",
                                    " --trim_poly_x --correction --detect_adapter_for_pe",
                                    " --dedup --dup_calc_accuracy 5 --compression 8"),
                task = "fastp-complete",
                report.tag = "fastp-complete",
                summary.csv = "logs/fastpComplete_summary.csv",
                threads = threads,
                mem = mem,
                overwrite = overwrite,
                quiet = quiet)

}#end function
