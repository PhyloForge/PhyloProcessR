#' @title reformatAlignmentFolder
#'
#' @description Converts all alignments in a directory from one file format to another. Supported input formats are "phylip", "nexus", and "fasta". Supported output formats are "phylip", "nexus", and "fasta". Each alignment is read with the appropriate ape reader and written in the target format to the output directory.
#'
#' @param alignment.path path to a directory containing the input alignment files
#'
#' @param out.dir path to the output directory where reformatted alignment files will be saved
#'
#' @param overwrite if TRUE, delete and recreate the output directory before writing; if FALSE, return early if the output directory already exists
#'
#' @param from.format the format of the input alignments: "phylip", "nexus", or "fasta"
#'
#' @param to.format the desired output format: "phylip", "nexus", or "fasta"
#'
#' @return saves reformatted alignment files to out.dir; nothing is returned to R
#'
#' @export

#function
reformatAlignmentFolder = function(alignment.path = NULL,
                                   out.dir = NULL,
                                   overwrite = FALSE,
                                   from.format = "nexus",
                                   to.format = "phylip") {

  #Gets list of alignments from path
  if (overwrite == TRUE){
    if (dir.exists(out.dir) == TRUE){ system(paste0("rm -r ", out.dir)) }
    dir.create(out.dir)
  } else {
    if (dir.exists(out.dir) == TRUE){ return("Directory exists and overwrite == FALSE") }
    dir.create(out.dir)
  }#end overwrite if

  if (alignment.path == out.dir){ stop("Error: output directory cannot be the same as input directory.")}

  align.names = list.files(alignment.path)

  #Loops through each alignment to gather statistics
  for (x in 1:length(align.names)){
    #Reads in alignment

    if (from.format == "phylip"){
      save.name = gsub(".phy$", "", align.names[x])
      align = ape::read.dna(paste0(alignment.path, "/", align.names[x]),
                            format = "sequential")
    }

    if (from.format == "nexus"){
      save.name = gsub(".nexus$|.nex$", "", align.names[x])
      l.align = ape::read.nexus.data(paste0(alignment.path, "/", align.names[x]))
      align = ape::as.DNAbin(matrix(unlist(l.align), ncol = length(l.align[[1]]), byrow = TRUE))
      rownames(align) = names(l.align)
    }

    if (from.format == "fasta"){
      save.name = gsub(".fasta$|.fa$", "", align.names[x])
      align = ape::read.FASTA(paste0(alignment.path, "/", align.names[x]), type = "DNA")
    }

    if (to.format == "phylip"){
      align = as.matrix(align)
      rownames(align) = labels(align)
      writePhylip(align, file = paste0(out.dir, "/",  save.name, ".phy"))
    }

    if (to.format == "nexus"){
      ape::write.nexus.data(align,  file = paste0(out.dir, "/",  save.name, ".nexus"), interleaved = FALSE)
    }

    #Saves as a fasta file
    if (to.format == "fasta"){
      temp.align = as.character(as.list(align))
      temp.align2 = lapply(temp.align, FUN = function(x) paste(x, collapse = ""))
      align.out = Biostrings::DNAStringSet(unlist(temp.align2))
      write.loci = as.list(as.character(align.out))
      seqinr::write.fasta(sequences = write.loci, names = names(write.loci),
                          paste0(out.dir, "/",  save.name, ".fa"),
                          nbchar = 1000000, as.string = T)
    }#end fasta if

  }#end x loop

}#end function

