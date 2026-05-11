#' @title setupCheck
#'
#' @description Verifies that the required external bioinformatics programs are
#'   accessible on the system. Checks for fastp, samtools, bwa, spades.py,
#'   bbmap.sh, bbnorm.sh, blastn, makeblastdb, mafft, iqtree2, and trimal. If
#'   an anaconda.environment path is provided all tool paths are automatically
#'   set to that environment's bin/ directory before checking. Prints a FOUND or
#'   NOT FOUND message for each program and returns a logical indicating overall
#'   pass/fail.
#'
#' @param anaconda.environment path to a conda environment directory (e.g.
#'   "~/miniconda3/envs/PhyloCap"). If provided, all individual path arguments
#'   are overridden and set to anaconda.environment/bin. Set to NULL to supply
#'   paths individually.
#'
#' @param fastp.path system path to the directory containing fastp; NULL skips
#'   this check.
#'
#' @param samtools.path system path to the directory containing samtools; NULL
#'   skips this check.
#'
#' @param bwa.path system path to the directory containing bwa; NULL skips this
#'   check.
#'
#' @param spades.path system path to the directory containing spades.py; NULL
#'   skips this check.
#'
#' @param bbmap.path system path to the directory containing bbmap.sh; NULL
#'   skips this check.
#'
#' @param bbnorm.path system path to the directory containing bbnorm.sh; NULL
#'   skips this check.
#'
#' @param blast.path system path to the directory containing blastn and
#'   makeblastdb; NULL skips this check.
#'
#' @param mafft.path system path to the directory containing mafft; NULL skips
#'   this check.
#'
#' @param iqtree.path system path to the directory containing iqtree2; NULL
#'   skips this check.
#'
#' @param trimAl.path system path to the directory containing trimal; NULL skips
#'   this check.
#'
#' @param julia.path system path to the directory containing julia; NULL skips
#'   this check (currently not checked by the function body).
#'
#' @param taper.path system path to the directory containing taper; NULL skips
#'   this check (currently not checked by the function body).
#'
#' @return logical TRUE if all checked programs were found, FALSE otherwise.
#'
#' @export

setupCheck = function(anaconda.environment = "conda/PhyloCap",
                      fastp.path = NULL,
                      samtools.path = NULL,
                      bwa.path = NULL,
                      spades.path = NULL,
                      bbmap.path = NULL,
                      bbnorm.path = NULL,
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
    bbnorm.path = paste0(anaconda.environment, "/bin")
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
  #CHecks if it was even inputted
  if (is.null(fastp.path) == TRUE) { pass = TRUE }

  if (file.exists(paste0(samtools.path, "/samtools")) == TRUE){
    print("Samtools was found.")
  } else {
    pass = FALSE
    print("Samtools could not be found.")
  } #end else
  #CHecks if it was even inputted
  if (is.null(samtools.path) == TRUE) { pass = TRUE }

  if (file.exists(paste0(bwa.path, "/bwa")) == TRUE){
    print("BWA was found.")
  } else {
    pass = FALSE
    print("BWA could not be found.")
  } #end else
  #CHecks if it was even inputted
  if (is.null(bwa.path) == TRUE) { pass = TRUE }

  if (file.exists(paste0(spades.path, "/spades.py")) == TRUE){
    print("Spades was found.")
  } else {
    pass = FALSE
    print("Spades could not be found.")
  } #end else
  #CHecks if it was even inputted
  if (is.null(spades.path) == TRUE) { pass = TRUE }

  if (file.exists(paste0(bbmap.path, "/bbmap.sh")) == TRUE){
    print("BBMap was found.")
  } else {
    pass = FALSE
    print("BBMap could not be found.")
  } #end else
  #CHecks if it was even inputted
  if (is.null(bbmap.path) == TRUE) { pass = TRUE }

  if (file.exists(paste0(bbnorm.path, "/bbnorm.sh")) == TRUE){
    print("BBNorm was found.")
  } else {
    pass = FALSE
    print("BBNorm could not be found.")
  } #end else
  #CHecks if it was even inputted
  if (is.null(bbnorm.path) == TRUE) { pass = TRUE }

  if (file.exists(paste0(blast.path, "/blastn")) == TRUE){
    print("blastn was found.")
  } else {
    pass = FALSE
    print("blastn could not be found.")
  } #end else
  #CHecks if it was even inputted
  if (is.null(blast.path) == TRUE) { pass = TRUE }

  if (file.exists(paste0(blast.path, "/makeblastdb")) == TRUE){
    print("makeblastdb was found.")
  } else {
    pass = FALSE
    print("makeblastdb could not be found.")
  } #end else
  #CHecks if it was even inputted
  if (is.null(blast.path) == TRUE) { pass = TRUE }

  if (file.exists(paste0(mafft.path, "/mafft")) == TRUE){
    print("mafft was found.")
  } else {
    pass = FALSE
    print("mafft could not be found.")
  } #end else
  #CHecks if it was even inputted
  if (is.null(mafft.path) == TRUE) { pass = TRUE }

  if (file.exists(paste0(iqtree.path, "/iqtree2")) == TRUE){
    print("iqtree was found.")
  } else {
    pass = FALSE
    print("iqtree could not be found.")
  } #end else
  #CHecks if it was even inputted
  if (is.null(iqtree.path) == TRUE) { pass = TRUE }

  if (file.exists(paste0(trimAl.path, "/trimal")) == TRUE){
    print("trimal was found.")
  } else {
    pass = FALSE
    print("trimal could not be found.")
  } #end else
  #CHecks if it was even inputted
  if (is.null(trimAl.path) == TRUE) { pass = TRUE }

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

