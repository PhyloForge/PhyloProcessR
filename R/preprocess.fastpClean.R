#' @title fastpClean
#'
#' @description Cleans raw paired-end Illumina reads with a single pass of
#'   fastp. The fastp command is built from the logical arguments below, so one
#'   pass does every step that is set to TRUE. This is much faster than a
#'   separate pass for each step, because the reads are read and written once.
#'   Set an argument to FALSE to leave that step out of the command.
#'
#' @param input.reads path to a directory of raw reads. Each sample may be
#'   organised in its own sub-directory or identified by a shared filename
#'   prefix containing "_L00".
#'
#' @param output.directory path to the directory where cleaned reads will be
#'   saved (one sub-directory per sample).
#'
#' @param remove.adaptors logical; TRUE trims adapters and detects the adapter
#'   sequence from the read pairs (--detect_adapter_for_pe). FALSE disables
#'   adapter trimming.
#'
#' @param remove.duplicate.reads logical; TRUE removes PCR and optical
#'   duplicate read pairs (--dedup).
#'
#' @param error.correction logical; TRUE corrects mismatched bases in the
#'   overlap of a read pair by majority vote (--correction).
#'
#' @param quality.trim.reads logical; TRUE trims low-quality bases from both
#'   read ends with a sliding window (--cut_front --cut_tail). This shortens
#'   reads and can hurt assembly, so the default is FALSE.
#'
#' @param quality.filter logical; TRUE keeps the fastp quality filter, which
#'   discards a read when more than 40 percent of its bases are below Q15.
#'   FALSE disables the filter.
#'
#' @param low.complexity.filter logical; TRUE discards low-complexity reads
#'   below a complexity threshold of 30 percent.
#'
#' @param trim.poly.x logical; TRUE trims poly-X tails, for example the poly-G
#'   tails of two-colour Illumina chemistry.
#'
#' @param min.read.length minimum read length in bp. A shorter read is
#'   discarded. Set to 0 or NULL to disable the length filter.
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
#'   logs/sample_logs, and a CSV summary to logs/fastpClean_summary.csv.
#'
#' @export

fastpClean = function(input.reads = NULL,
                      output.directory = NULL,
                      remove.adaptors = TRUE,
                      remove.duplicate.reads = TRUE,
                      error.correction = TRUE,
                      quality.trim.reads = FALSE,
                      quality.filter = TRUE,
                      low.complexity.filter = TRUE,
                      trim.poly.x = TRUE,
                      min.read.length = 60,
                      fastp.path = NULL,
                      threads = 1,
                      mem = 8,
                      overwrite = FALSE,
                      quiet = TRUE) {

  #Quick check that the command would do something to the reads
  if (all(c(remove.adaptors, remove.duplicate.reads, error.correction,
            quality.trim.reads, quality.filter, low.complexity.filter,
            trim.poly.x) == FALSE) &&
      (is.null(min.read.length) == TRUE || min.read.length <= 0)){
    stop(paste0("Every fastpClean step is FALSE, so the reads would not change.",
                " Skip this step instead, or set a step to TRUE."))
  }

  # Builds the fastp command from the arguments above. fastp trims adapters,
  # filters on quality and filters on length by default, so a step that is
  # FALSE needs its --disable flag rather than the absence of a flag.
  fastp.args = c()

  if (remove.adaptors == TRUE){
    fastp.args = c(fastp.args, "--detect_adapter_for_pe")
  } else { fastp.args = c(fastp.args, "--disable_adapter_trimming") }

  if (quality.filter == FALSE){
    fastp.args = c(fastp.args, "--disable_quality_filtering")
  }

  if (quality.trim.reads == TRUE){
    fastp.args = c(fastp.args, "--cut_front --cut_tail")
  }

  if (is.null(min.read.length) == TRUE || min.read.length <= 0){
    fastp.args = c(fastp.args, "--disable_length_filtering")
  } else { fastp.args = c(fastp.args, paste0("--length_required ", min.read.length)) }

  if (low.complexity.filter == TRUE){
    fastp.args = c(fastp.args, "--low_complexity_filter --complexity_threshold 30")
  }

  if (trim.poly.x == TRUE){ fastp.args = c(fastp.args, "--trim_poly_x") }

  if (remove.duplicate.reads == TRUE){
    fastp.args = c(fastp.args, "--dedup --dup_calc_accuracy 5")
  }

  if (error.correction == TRUE){ fastp.args = c(fastp.args, "--correction") }

  fastp.args = c(fastp.args, "--compression 6")

  if (quiet == FALSE){
    print(paste0("fastp arguments: ", paste(fastp.args, collapse = " ")))
  }

  .runFastpStep(input.reads = input.reads,
                output.directory = output.directory,
                fastp.path = fastp.path,
                fastp.args = paste(fastp.args, collapse = " "),
                task = "fastp-clean",
                report.tag = "fastp-clean",
                summary.csv = "logs/fastpClean_summary.csv",
                threads = threads,
                mem = mem,
                overwrite = overwrite,
                quiet = quiet)

}#end function
