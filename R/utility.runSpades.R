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
#'   read files (READ1, READ2, and optionally merged READ3). More than 3 paths
#'   raise an error.
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

  if (length(read.paths) == 0){ stop("No read files were supplied.") }
  if (length(read.paths) > 3){ stop("runSpades accepts 1, 2, or 3 read files.") }
  if (file.exists(read.paths[1]) == FALSE){ stop("Read files not found.") }
  if (overwrite == T){
    if (dir.exists("spades") == TRUE){ unlink("spades", recursive = TRUE) }
  }#end


  if (isolate == TRUE && mismatch.corrector == TRUE) {stop("Both --isolate or --careful (mismatch corrector) can not be used together. Please choose only one.")}
  if (mismatch.corrector == FALSE && isolate == FALSE){ mismatch.string = "" }
  if (isolate == TRUE){ mismatch.string = "--isolate " }
  if (mismatch.corrector == TRUE){ mismatch.string = "--careful " }

  #Run SPADES on sample
  k = kmer.values

  #Removes the largest k-mer after each failure. The k-mer vector is only reduced
  #when the previous attempt produced no contigs, so a run that succeeds on the
  #last k-mer value is kept.
  if (file.exists("spades/contigs.fasta") == F){
    repeat {
      k.val = paste(k, collapse = ",")

      #Single end reads
      if (length(read.paths) == 1){
        system(paste0(full.path.spades, "spades.py --s1 ", shQuote(read.paths[1]),
                      " -o spades -k ",k.val," ", mismatch.string, "-t ", threads, " -m ", memory),
               ignore.stdout = quiet)
      }#end 1 read

      if (length(read.paths)  == 2){
        system(paste0(full.path.spades, "spades.py --pe1-1 ", shQuote(read.paths[1]),
                      " --pe1-2 ", shQuote(read.paths[2]),
                      " -o spades -k ",k.val," ", mismatch.string, "-t ", threads, " -m ", memory),
               ignore.stdout = quiet)
      }#end 2 reads

      if (length(read.paths)  == 3){
        system(paste0(full.path.spades, "spades.py --pe1-1 ", shQuote(read.paths[1]),
                      " --pe1-2 ", shQuote(read.paths[2]), " --pe1-m ", shQuote(read.paths[3]),
                      " -o spades -k ",k.val, " ", mismatch.string, "-t ", threads, " -m ", memory),
               ignore.stdout = quiet)
      }#end 3 reads

      if (file.exists("spades/contigs.fasta") == TRUE) { break }
      #subtract Ks until it works
      k = k[-length(k)]
      if (length(k) == 0) { break }
    }#end repeat
  }#end assembly

  #If the k-mers are all run out, therefore nothing can be assembled
  if (file.exists("spades/contigs.fasta") == F) {
    print("k-mer values all used up, cannot assemble!")
    unlink("spades", recursive = TRUE)
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
      file.copy("spades/scaffolds.fasta", paste0(save.name, ".fa"), overwrite = TRUE)
    } else {
      file.copy("spades/contigs.fasta", paste0(save.name, ".fa"), overwrite = TRUE)
    }#end else
  }#end save file

  if (clean == TRUE){ unlink("spades", recursive = TRUE) }

  if (read.contigs == T) {return(contigs) }
  if (is.null(save.name) == F) {return("Contigs were saved to file.") }

  return("Nothing was saved.")
  ##############################
}# end spades function
