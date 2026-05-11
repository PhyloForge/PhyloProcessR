#' @title trimTAPER
#'
#' @description Wrapper function for running TAPER (via its Julia script correction.jl) on a single alignment. The alignment is written to a temporary fasta file, the TAPER Julia script is called to mask low-quality sequence regions (replacing them with N), and the corrected alignment is read back. All-gap columns introduced by TAPER are then removed with trimAlignmentColumns. If TAPER produces no output the original alignment is returned unchanged. TAPER and Julia must be installed and accessible.
#'
#' @param alignment a DNAStringSet containing the aligned sequences to clean
#'
#' @param TAPER.path path to the TAPER correction.jl Julia script; defaults to "correction.jl"
#'
#' @param julia.path path to the julia executable; defaults to "julia"
#'
#' @param quiet if TRUE, suppress TAPER screen output
#'
#' @param delete.temp if TRUE, delete the temporary fasta files created during the TAPER run
#'
#' @return a DNAStringSet of the TAPER-cleaned alignment with all-gap columns removed; returns the original alignment unchanged if TAPER produces no output
#'
#' @export

trimTAPER = function(alignment = NULL,
                     TAPER.path = "correction.jl",
                     julia.path = "julia",
                     quiet = FALSE,
                     delete.temp = TRUE) {

  # #Debug section
  #  alignment = non.align
  #  TAPER.path = "/usr/local/bin/correction.jl"
  #  julia.path = "/Applications/Julia-1.5.app/Contents/Resources/julia/bin/julia"
  #  quiet = FALSE
  #  delete.temp = TRUE

  if (length(alignment) <= 3){ return(alignment) }

  #Creates random name and saves it
  write.align = as.list(as.character(alignment))
  input.file = paste0("temp-", sample(1:1000000, 1), ".fa")
  writeFasta(sequences = write.align,
             names = names(write.align),
             file.out = input.file,
             nbchar = 1000000,
             as.string = T)

  #Runs the thing
  system(paste0(julia.path, "julia ", TAPER.path, "correction.jl -m N -a N ",
                input.file, " > out_", input.file),
         ignore.stdout = F, ignore.stderr = quiet)

  #Deletes input
  system(paste0("rm ", input.file))

  #Checks for file and loads
  if (file.exists(paste0("out_", input.file)) == TRUE){
    taper.align = Rsamtools::scanFa(Rsamtools::FaFile(paste0("out_", input.file)) )
  } else {
    if (delete.temp == TRUE){ system(paste0("out_", input.file)) }
    return(alignment)
  }#end file check

  if (delete.temp == TRUE){ system(paste0("rm ", "out_", input.file)) }

  #Removes the edge gaps
  out.align = trimAlignmentColumns(alignment = taper.align,
                                   min.gap.percent = 100)

  return(out.align)

}#end function
