#' @title runSpades
#'
#' @description Function for running the program spades to assemble short read sequencing data
#'
#' @param read.paths path to a folder of sequence alignments in phylip format.
#'
#' @param full.path.spades contigs are added into existing alignment if algorithm is "add"
#'
#' @param mismatch.corrector algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param read.contigs TRUE applies the adjust sequence direction function of MAFFT
#'
#' @param save.file if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param save.name if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param threads number of computation processing threads
#'
#' @param mem amount of system memory to use
#'
#' @param resume TRUE to skip samples already completed
#'
#' @param overwrite TRUE to overwrite a folder of samples with output.dir
#'
#' @param quiet TRUE to supress screen output

#' @return an alignment of provided sequences in DNAStringSet format. Also can save alignment as a file with save.name
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

runSpades = function(read.paths = NULL,
                     full.path.spades = "spades.py",
                     mismatch.corrector = FALSE,
                     read.contigs = T,
                     save.file = T,
                     save.name = NULL,
                     threads = 1,
                     memory = 4,
                     overwrite = T,
                     quiet = T) {

  # #debug
  # full.path.spades = spades.path
  # read.paths = temp.read.path
  # supress.print = T
  # save.file = F
  # read.contigs = T
  # mismatch.corrector = F
  # overwrite = T

  if (file.exists(read.paths[1]) == FALSE){ stop("Read files not found.") }
  if (overwrite == T){
    if (dir.exists("spades") == TRUE){ system(paste0("rm -r spades")) }
  }#end

  if (mismatch.corrector == TRUE){ mismatch.string = "--careful " }
  if (mismatch.corrector == FALSE){ mismatch.string = "" }

  #Run SPADES on sample
  k = c(9,13,21,33,55,77,99,127)
  k.val = paste(k, collapse = ",")

  #Checks to see if one kmer failed or not
  while (file.exists("spades/contigs.fasta") == F){
    #stop("the while loop messed up K")
    #if (counter == 1){

    #Single end reads
    if (length(read.paths) == 1){
      system(paste0(full.path.spades, " --s1 ", read.paths[1],
                    " -o spades -k ",k.val," ", mismatch.string, "-t ", threads, " -m ", memory),
             ignore.stdout = quiet)
    }#end 2 reads

    if (length(read.paths)  == 2){
      system(paste0(full.path.spades, " --pe1-1 ", read.paths[1], " --pe1-2 ", read.paths[2],
                    " -o spades -k ",k.val," ", mismatch.string, "-t ", threads, " -m ", memory),
             ignore.stdout = quiet)
    }#end 2 reads

    if (length(read.paths)  == 3){
      system(paste0(full.path.spades, " --pe1-1 ", read.paths[1],
                    " --pe1-2 ", read.paths[2], " --merged ", read.paths[3],
                    " -o spades -k ",k.val, " ", mismatch.string, "-t ", threads, " -m ", memory),
             ignore.stdout = quiet)
    }#end 3 reads
    #subtract Ks until it works
    k = k[-length(k)]
    if (length(k) == 0) { break }
    k.val = paste(k, collapse = ",")
  }#end while

  #If the k-mers are all run out, therefore nothing can be assembled
  if (length(k) == 0) {
    print("k-mer values all used up, cannot assemble!")
    system("rm -r spades")
    contigs = Biostrings::DNAStringSet()
    return(contigs)
  }#end k

  if (read.contigs == T){
    if (file.exists("spades/contigs.fasta") == TRUE){
      contigs = Rsamtools::scanFa(Rsamtools::FaFile("spades/contigs.fasta"))
    } else {
      contigs = Rsamtools::scanFa(Rsamtools::FaFile("spades/scaffolds.fasta"))
    }#end else
  } #end if

  if (save.file == T){
    if (file.exists("spades/contigs.fasta") == TRUE){
      system(paste0("cp spades/contigs.fasta ", save.name, ".fa"))
    } else {
      system(paste0("cp spades/scaffolds.fasta ", save.name, ".fa"))
    }#end else
  }#end save file
  system("rm -r spades")

  if (length(contigs) == 0){
    print("No contigs were assembled.")
    return(contigs) }

  if (read.contigs == T) {return(contigs) }
  if (save.file == T) {return("Contigs were saved to file.") }

  return("Nothing was saved.")
  ##############################
}# end spades function
