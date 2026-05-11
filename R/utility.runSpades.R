#' @title runSpades
#'
#' @description Assembles short reads using SPAdes. Supports single-end,
#'   paired-end, and paired-end with merged reads (1, 2, or 3 read files).
#'   If assembly fails for the full set of k-mer values, the largest k-mer is
#'   progressively removed and SPAdes is re-run until assembly succeeds or all
#'   k-mer values are exhausted. The assembled scaffolds (preferred) or contigs
#'   are optionally read back into R and/or saved to a file.
#'
#' @param read.paths character vector of 1, 2, or 3 paths to the input fastq.gz
#'   read files (READ1, READ2, and optionally merged READ3).
#'
#' @param full.path.spades system path to the directory containing spades.py;
#'   NULL searches the system PATH.
#'
#' @param mismatch.corrector logical; if TRUE passes --careful to SPAdes to
#'   enable the mismatch correction module. Cannot be used with isolate = TRUE.
#'
#' @param isolate logical; if TRUE passes --isolate to SPAdes, which is
#'   optimised for high-coverage isolate data. Cannot be used with
#'   mismatch.corrector = TRUE.
#'
#' @param kmer.values integer vector of k-mer sizes to try; if the largest
#'   k-mer causes failure it is dropped and SPAdes is rerun with the remaining
#'   values.
#'
#' @param read.contigs logical; if TRUE the assembled scaffolds (or contigs if
#'   no scaffold file exists) are read into R and returned as a DNAStringSet.
#'
#' @param save.name base name (without extension) for an output FASTA file
#'   where the assembly is copied; NULL skips file saving.
#'
#' @param clean logical; if TRUE the spades/ working directory is deleted after
#'   assembly.
#'
#' @param threads number of CPU threads to pass to SPAdes.
#'
#' @param memory amount of RAM in GB to allocate to SPAdes.
#'
#' @param overwrite logical; if TRUE an existing spades/ directory is deleted
#'   before running.
#'
#' @param quiet logical; if TRUE SPAdes stdout is suppressed.
#'
#' @return if read.contigs is TRUE, a DNAStringSet of assembled sequences; if
#'   save.name is provided, a character string "Contigs were saved to file.";
#'   otherwise "Nothing was saved.". An empty DNAStringSet is returned if no
#'   k-mer values succeed.
#'
#' @export

runSpades = function(read.paths = NULL,
                     full.path.spades = NULL,
                     mismatch.corrector = TRUE,
                     isolate = FALSE,
                     kmer.values = c(33,55,77,99,127),
                     read.contigs = F,
                     save.name = NULL,
                     clean = FALSE,
                     threads = 1,
                     memory = 4,
                     overwrite = T,
                     quiet = T) {

  # #debug
  # full.path.spades = spades.path
  # read.paths = temp.read.path
  # #read.paths = paste0("/Volumes/LaCie/Mantellidae/Wakea_madinika_2001F54/", list.files("/Volumes/LaCie/Mantellidae/Wakea_madinika_2001F54"))
  # quiet = F
  # save.name = "iterative_temp/contigs"
  # clean = T
  # read.contigs = F
  # mismatch.corrector = F
  # isolate = T
  # overwrite = T
  # kmer.values = c(21,33,55,77,99,127)
  # memory = 1024
  # threads = 8

  #Same adds to bbmap path
  if (is.null(full.path.spades) == FALSE){
    b.string = unlist(strsplit(full.path.spades, ""))
    if (b.string[length(b.string)] != "/") {
      full.path.spades = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { full.path.spades = "" }

  if (file.exists(read.paths[1]) == FALSE){ stop("Read files not found.") }
  if (overwrite == T){
    if (dir.exists("spades") == TRUE){ system(paste0("rm -r spades")) }
  }#end


  if (isolate == TRUE && mismatch.corrector == TRUE) {stop("Both --isolate or --careful (mismatch corrector) can not be used together. Please choose only one.")}
  if (mismatch.corrector == FALSE && isolate == FALSE){ mismatch.string = "" }
  if (isolate == TRUE){ mismatch.string = "--isolate " }
  if (mismatch.corrector == TRUE){ mismatch.string = "--careful " }

  #Run SPADES on sample
  k = kmer.values
  k.val = paste(k, collapse = ",")

  #Checks to see if one kmer failed or not
  while (file.exists("spades/contigs.fasta") == F){
    #stop("the while loop messed up K")
    #if (counter == 1){

    #Single end reads
    if (length(read.paths) == 1){
      system(paste0(full.path.spades, "spades.py --s1 ", read.paths[1],
                    " -o spades -k ",k.val," ", mismatch.string, "-t ", threads, " -m ", memory),
             ignore.stdout = quiet)
    }#end 2 reads

    if (length(read.paths)  == 2){
      system(paste0(full.path.spades, "spades.py --pe1-1 ", read.paths[1], " --pe1-2 ", read.paths[2],
                    " -o spades -k ",k.val," ", mismatch.string, "-t ", threads, " -m ", memory),
             ignore.stdout = quiet)
    }#end 2 reads

    if (length(read.paths)  == 3){
      system(paste0(full.path.spades, "spades.py --pe1-1 ", read.paths[1],
                    " --pe1-2 ", read.paths[2], " --pe1-m ", read.paths[3],
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
    if (file.exists("spades/scaffolds.fasta") == TRUE){
      contigs = Biostrings::readDNAStringSet("spades/scaffolds.fasta")
    } else {
      contigs = Biostrings::readDNAStringSet("spades/contigs.fasta")
    }#end else

    if (length(contigs) == 0){
      print("No contigs were assembled.")
      return(contigs) }

  } #end if

  if (is.null(save.name) == FALSE){
    if (file.exists("spades/scaffolds.fasta") == TRUE){
      system(paste0("cp spades/scaffolds.fasta ", save.name, ".fa"))
    } else {
      system(paste0("cp spades/contigs.fasta ", save.name, ".fa"))
    }#end else
  }#end save file

  if (clean == TRUE){ system("rm -r spades") }

  if (read.contigs == T) {return(contigs) }
  if (is.null(save.name) == F) {return("Contigs were saved to file.") }

  return("Nothing was saved.")
  ##############################
}# end spades function
