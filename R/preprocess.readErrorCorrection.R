#' @title readErrorCorrection
#'
#' @description Performs overlap-based read error correction on paired-end
#'   fastq reads using fastp's --correction flag. Only the correction step is
#'   applied; adapter trimming, quality filtering, and length filtering are all
#'   disabled. fastp corrects mismatches in the overlapping region of paired
#'   reads by majority vote.
#'
#' @param input.reads path to a directory of processed paired-end reads in
#'   fastq.gz format, organised in per-sample sub-directories.
#'
#' @param output.directory path to the directory where error-corrected reads
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
#' @return invisibly returns the summary data frame; writes error-corrected
#'   fastq.gz files to output.directory, the fastp HTML and JSON reports to
#'   logs/sample_logs, and a CSV summary to
#'   logs/readErrorCorrection_summary.csv.
#'
#' @export

readErrorCorrection = function(input.reads = "adaptor-removed-reads",
                               output.directory = "corrected-reads",
                               fastp.path = NULL,
                               threads = 1,
                               mem = 1,
                               overwrite = FALSE,
                               quiet = TRUE) {

  .runFastpStep(input.reads = input.reads,
                output.directory = output.directory,
                fastp.path = fastp.path,
                fastp.args = paste0("--correction",
                                    " --disable_adapter_trimming --disable_quality_filtering",
                                    " --disable_length_filtering --compression 6"),
                task = "error-correction",
                report.tag = "error-correct_fastp",
                summary.csv = "logs/readErrorCorrection_summary.csv",
                threads = threads,
                mem = mem,
                overwrite = overwrite,
                quiet = quiet)

}#end function
