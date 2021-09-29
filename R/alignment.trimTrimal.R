#' @title trimTrimal
#'
#' @description wrapper function for the program TrimAl that corrects for NAs introduced from TrimAl removing samples from the alignment
#'
#' @param alignment path to a folder of sequence alignments in phylip format.
#'
#' @param trimal.path provide system absolute path to trimal if R cannot find it
#'
#' @param method trimming method, "auto" is automated. Other options coming soon.
#'
#' @param quiet TRUE to supress TrimAl screen output
#'
#' @return returns an alignment in DNAStringSet format
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

trimTrimal = function(alignment = NULL,
                      trimal.path = NULL,
                      method = "auto",
                      quiet = TRUE) {
  #Debug
   # alignment = non.align
   # quiet = FALSE
   # trimal.path = "/Users/chutter/conda/PhyloCap/bin"

  #Same adds to bbmap path
  if (is.null(trimal.path) == FALSE){
    b.string = unlist(strsplit(trimal.path, ""))
    if (b.string[length(b.string)] != "/") {
      trimal.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { trimal.path = "" }

  if (length(alignment) <= 3){ return(alignment) }

  #Finds probes that match to two or more contigs
  save.rownames = names(alignment)
  write.align = as.list(as.character(alignment))

  #Creates random name and saves it
  input.file = paste0("temp_", sample(1:1000000, 1), ".fa")
  writeFasta(sequences = write.align,
             names = names(write.align),
             file.out = input.file,
             nbchar = 1000000,
             as.string = T)

  #Runs trimal command with input file
  system(paste0(trimal.path, "trimal -in ", input.file, " -out tm-", input.file, " -automated1"),
         ignore.stdout = quiet, ignore.stderr = quiet)

  if (file.exists(paste0("tm-", input.file)) == F) {
    system(paste0("rm ", input.file))
    return(alignment)
  } else { system(paste0("mv tm-", input.file, " ", input.file)) }

  out.align = Rsamtools::scanFa(Rsamtools::FaFile(input.file))

  #Fixes any terrible NA names introduced by trimal
  new.names = c()
  for (j in 1:length(names(out.align))){
    new.names[j] = save.rownames[grep(pattern = paste0(names(out.align)[j], "$"), x = save.rownames)]
  }

  temp = names(out.align)[is.na(names(out.align)) == T]
  if (length(temp) > 0){ stop("there are NAs in the names") }
  names(out.align) = new.names
  system(paste0("rm ", input.file))
  return(out.align)

}#end function
