#' @title setupCheck
#'
#' @description Function for removing contamination from other organisms from adaptor trimmed Illumina sequence data using BWA
#'
#' @param anaconda.environment path to a folder of adaptor trimmed reads in fastq format.
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

setupCheck = function(anaconda.environment = "conda/PhyloCap",
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

  #anaconda.environment = "/Users/chutter/conda/PhyloCap"

  if (is.null(anaconda.environment) == FALSE) {
    fastp.path = paste0(anaconda.environment, "/bin")
    samtools.path = paste0(anaconda.environment, "/bin")
    bwa.path = paste0(anaconda.environment, "/bin")
    spades.path = paste0(anaconda.environment, "/bin")
    bbmap.path = paste0(anaconda.environment, "/bin")
    blast.path = paste0(anaconda.environment, "/bin")
    mafft.path = paste0(anaconda.environment, "/bin")
    iqtree.path = paste0(anaconda.environment, "/bin")
    trimAl.path = paste0(anaconda.environment, "/bin")
    taper.path = paste0(anaconda.environment, "/bin")
    julia.path = paste0(anaconda.environment, "/bin")
  }

  #Check paths from above
  pass = TRUE
  if (file.exists(paste0(fastp.path, "/fastp")) == TRUE){
    print("Fastp was found.")
  } else {
    pass = FALSE
    print("Fastp could not be found.")
  } #end else

  if (file.exists(paste0(samtools.path, "/samtools")) == TRUE){
    print("Samtools was found.")
  } else {
    pass = FALSE
    print("Samtools could not be found.")
  } #end else

  if (file.exists(paste0(bwa.path, "/bwa")) == TRUE){
    print("BWA was found.")
  } else {
    pass = FALSE
    print("BWA could not be found.")
  } #end else

  if (file.exists(paste0(spades.path, "/spades.py")) == TRUE){
    print("Spades was found.")
  } else {
    pass = FALSE
    print("Spades could not be found.")
  } #end else

  if (file.exists(paste0(bbmap.path, "/bbmap.sh")) == TRUE){
    print("BBMap was found.")
  } else {
    pass = FALSE
    print("BBMap could not be found.")
  } #end else

  if (file.exists(paste0(blast.path, "/blastn")) == TRUE){
    print("blastn was found.")
  } else {
    pass = FALSE
    print("blastn could not be found.")
  } #end else

  if (file.exists(paste0(blast.path, "/makeblastdb")) == TRUE){
    print("makeblastdb was found.")
  } else {
    pass = FALSE
    print("makeblastdb could not be found.")
  } #end else

  if (file.exists(paste0(mafft.path, "/mafft")) == TRUE){
    print("mafft was found.")
  } else {
    pass = FALSE
    print("mafft could not be found.")
  } #end else

  if (file.exists(paste0(iqtree.path, "/iqtree2")) == TRUE){
    print("iqtree was found.")
  } else {
    pass = FALSE
    print("iqtree could not be found.")
  } #end else

  if (file.exists(paste0(trimAl.path, "/trimal")) == TRUE){
    print("trimal was found.")
  } else {
    pass = FALSE
    print("trimal could not be found.")
  } #end else
#
#   if (file.exists(paste0(taper.path, "/taper")) == TRUE){
#     print("TAPER was found.")
#   } else {
#     pass = FALSE
#     print("TAPER could not be found.")
#   } #end else
#
#   if (file.exists(paste0(julia.path, "/julia")) == TRUE){
#     print("Julia was found.")
#   } else {
#     pass = FALSE
#     print("Julia could not be found.")
#   } #end else

  return(pass)

}#end function

