#' @title trimTAPER
#'
#' @description wrapper function for running TAPER Must be installed.
#'
#' @param alignment alignment in DNAStringSet format
#'
#' @param TAPER.path Absolute path to TAPER if R cannot find it in your path
#'
#' @param julia.path Absolute path to julia if R cannot find it in your path
#'
#' @param quiet TRUE to supress TAPER screen output
#'
#' @param delete.temp TRUE to delete temporary files made by TAPIR
#'
#' @return returns a data.table with the raw summary statistics calculated for each alignment in the alignment.path. A csv file can optionally be saved by giving a file name to file.export
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
#'
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
