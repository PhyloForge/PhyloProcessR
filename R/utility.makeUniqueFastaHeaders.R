#' @title makeUniqueFastaHeaders
#'
#' @description Resolves duplicate FASTA sequence headers in a single file.
#'   Duplicate entries are detected and either renamed by appending a numeric
#'   suffix ("append-number") or replaced with sequential Sequence_N names
#'   ("rename"). The deduplicated sequences are saved to a new FASTA file using
#'   writeFasta().
#'
#' @param fasta.file path to the input FASTA file that may contain duplicate
#'   sequence headers.
#'
#' @param output.name path or file name for the output FASTA file with unique
#'   headers.
#'
#' @param type character; "append-number" to append a numeric suffix to
#'   duplicate header names, or "rename" to replace all headers with
#'   zero-padded sequential names (Sequence_0001, Sequence_0002, ...).
#'
#' @param overwrite logical; if TRUE an existing output file is deleted before
#'   writing.
#'
#' @param quiet logical; currently unused, reserved for future output control.
#'
#' @return invisibly; the deduplicated sequences are written to output.name.
#'
#' @export

makeUniqueFastaHeaders = function(fasta.file = NULL,
                                  output.name = "unique-headers",
                                  type = c("append-number", "rename"),
                                  overwrite = FALSE,
                                  quiet = TRUE) {


  #Read in basic genome info
  #fasta.directory = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Mantellidae_Subfamily/Mantella_Genome/maker-analysis/proteomes"
  #output.directory = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Mantellidae_Subfamily/Mantella_Genome/maker-analysis/shortened-headers"
  #number.characters = 70
  #fasta.file = "/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/Elapid_Probes/venom_loci_updated_Mar12_cdhit95.fa"
  #output.name = "/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/Elapid_Probes/venom_loci_updated_Mar12_cdhit95_unique.fa"
  #type = "append-number"

  if (is.null(fasta.file) == T){ stop("A directory of genome(s) is needed.") }

  #Overwrites
  if (file.exists(output.name) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm ", output.directory))
    } else { stop("file exists and overwrite = FALSE") }
  }

  old.file = Biostrings::readAAStringSet(fasta.file, format = "fasta")

  dup.entry = old.file[duplicated(names(old.file)) == T]
  keep.entry = old.file[duplicated(names(old.file)) == F]

  if (type == "append-number"){
    names(dup.entry) = paste0(names(dup.entry), "_", seq(1:length(dup.entry)))
  }

  save.entry = append(keep.entry, dup.entry)
  save.entry = save.entry[order(names(save.entry))]

  #Renames them all Sequence_X
  if (type == "rename"){
    names(save.entry) = paste0("Sequence_", stringr::str_pad(string = seq(1:length(save.entry)),
                                                             width = nchar(length(save.entry)),
                                                             side = "left",
                                                             pad = "0") )
  }

  #Creates random name and saves it
  write.file = as.list(as.character(save.entry))
  writeFasta(sequences = write.file,
             names = names(write.file),
             file.out = output.name,
             nbchar = 100000000,
             as.string = T)

}#end funtion


