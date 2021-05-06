#' @title setupCheck
#'
#' @description Function for removing contamination from other organisms from adaptor trimmed Illumina sequence data using BWA
#'
#' @param anaconda.bin path to a folder of adaptor trimmed reads in fastq format.
#'
#' @param fastp.path the new directory to save the adaptor trimmed sequences
#'
#' @param samtools.path directory of genomes contaminants to scan samples
#'
#' @param spades.path system path to samtools in case it can't be found
#'
#' @param bbmap.path system path to bwa in case it can't be found
#'
#' @param blast.path number of computation processing threads
#'
#' @param mafft.path amount of system memory to use
#'
#' @param iqtree.path TRUE to skip samples already completed
#'
#' @param trimAl.path TRUE to overwrite a folder of samples with output.dir
#'
#' @param julia.path TRUE to supress screen output
#'
#' @param taper.path TRUE to supress screen output
#'
#' @return a new directory of adaptor trimmed reads and a summary of the trimming in logs/
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

setupCheck = function(anaconda.bin = "PhyloCap/bin",
                      fastp.path = NULL,
                      samtools.path = NULL,
                      bwa.path = NULL,
                      spades.path = NULL,
                      bbmap.path = NULL,
                      blast.path = NULL,
                      mafft.path = NULL,
                      iqtree.path = NULL,
                      trimAl.path = NULL,
                      julia.path = NULL,
                      taper.path = NULL) {

  fastp.path = "/Users/chutter/conda/PhyloCap/bin"
  samtools.path = "/Users/chutter/conda/PhyloCap/bin"
  bwa.path = "/Users/chutter/conda/PhyloCap/bin"
  spades.path = "/Users/chutter/conda/PhyloCap/bin"
  bbmap.path = "/Users/chutter/conda/PhyloCap/bin"
  blast.path = "/Users/chutter/conda/PhyloCap/bin"
  mafft.path = "/Users/chutter/conda/PhyloCap/bin"
  iqtree.path = "/Users/chutter/conda/PhyloCap/bin"
  trimAl.path = "/Users/chutter/conda/PhyloCap/bin"
  taper.path = "/Users/chutter/conda/PhyloCap/bin"
  julia.path = "/Users/chutter/conda/PhyloCap/bin"

  #Same adds to bbmap path
  if (is.null(samtools.path) == FALSE){
    b.string = unlist(strsplit(samtools.path, ""))
    if (b.string[length(b.string)] != "/") {
      samtools.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { samtools.path = "" }

  #Same adds to bbmap path
  if (is.null(bwa.path) == FALSE){
    b.string = unlist(strsplit(bwa.path, ""))
    if (b.string[length(b.string)] != "/") {
      bwa.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { bwa.path = "" }



}#end function

