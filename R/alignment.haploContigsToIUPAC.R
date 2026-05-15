#' @title haploContigsToIUPAC
#'
#' @description Converts phased haplotype alignments (where each taxon is represented by
#' two sequences distinguished by a \code{_0} / \code{_1} suffix) into single-sequence
#' alignments using IUPAC ambiguity codes. For each taxon, the two haplotype sequences
#' are merged into a consensus using \code{makeConsensus} with the IUPAC method. The
#' resulting alignments are saved in phylip format. Processing is parallelised across
#' alignment files.
#'
#' @param alignment.directory path to the directory containing haplotype alignment files.
#'
#' @param alignment.format format of the input alignment files. Accepted values: "phylip"
#' or "fasta".
#'
#' @param output.directory path to the directory where IUPAC consensus alignment files
#' will be saved.
#'
#' @param output.format format for the output alignment files. Currently "phylip" is
#' supported.
#'
#' @param threads number of parallel threads to use. Default 1.
#'
#' @param memory total memory (in GB) to allocate across all threads. Default 1.
#'
#' @param overwrite logical. If TRUE, the output directory is deleted and recreated before
#' processing. If FALSE and the output directory already exists, the function stops with an
#' error. Default FALSE.
#'
#' @param quiet logical. If TRUE, suppresses per-file progress messages. Default FALSE.
#'
#' @return Writes IUPAC consensus alignment files in phylip format to
#' \code{output.directory}. No value is returned to R.
#'
#' @export

haploContigsToIUPAC = function(alignment.directory = NULL,
                               alignment.format = "phylip",
                               output.directory = NULL,
                               output.format = "phylip",
                               threads = 1,
                               memory = 1,
                               overwrite = FALSE,
                               quiet = FALSE) {



  # setwd("/Users/chutter/Dropbox/VCF Limnonectes/phyllofolia_SNP-Analysis/gphocs_input-phyllofolia")
  # alignment.directory = "Haplotype_contigs"
  # alignment.format = "phylip"
  # output.directory = "IUPAC_alignments"
  # output.format = "phylip"
  # threads = 6
  # memory = 12
  # overwrite = TRUE
  # quiet = TRUE


  require(PhyloCap)
  #Checks this
  if (alignment.directory == output.directory){ stop("You should not overwrite the original alignments.") }

  #Overwrite
  if (dir.exists(paste0(output.directory)) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(paste0(output.directory))
    } else {
      stop("Overwrite = FALSE and directory exists. Either change to TRUE or overwrite manually.")
    }
  } else {
    dir.create(paste0(output.directory))
  }#end else

  #Gathers alignments
  align.files = list.files(alignment.directory)

  if (length(align.files) == 0) { stop("alignment files could not be found.") }

  mem.cl = floor(memory/threads)

  #Loops through each locus and does operations on them
  parallel::mclapply(seq_along(align.files), function(i) {
  tryCatch({

    ##############
    #STEP 2: Runs MAFFT to add
    ##############
    #Load in alignments
    if (alignment.format == "phylip"){
      old.align = Biostrings::readAAMultipleAlignment(file = paste0(alignment.directory, "/", align.files[i]), format = "phylip")
      old.align = Biostrings::DNAStringSet(old.align)
      save.name = gsub(".phy$", "", align.files[i])
      save.name = gsub(".phylip$", "", save.name)
    }#end phylip

    if (alignment.format == "fasta"){
      old.align = Biostrings::readDNAStringSet(paste0(alignment.directory, "/", align.files[i]) )
      save.name = gsub(".fa$", "", align.files[i])
      save.name = gsub(".fasta$", "", save.name)
    }#end phylip

    taxa.names = unique(gsub("_0|_1", "", names(old.align)))

    final.align = Biostrings::DNAStringSet()
    for (j in 1:length(taxa.names)){

      temp.align = old.align[gsub("_0|_1", "", names(old.align)) %in% taxa.names[j]]
      n.align = Biostrings::DNAStringSet(gsub("-", "N", as.character(temp.align)))

      #makes a consensus using ambig codes
      new.align = makeConsensus(alignment = n.align,
                                method = "IUPAC",
                                type = "DNA")

      names(new.align) = taxa.names[j]
      final.align = append(final.align, new.align)
    }#end j loop

    #save file
    #Saves them
    write.temp = strsplit(as.character(final.align), "")
    aligned.set = as.matrix(ape::as.DNAbin(write.temp) )

    #readies for saving
    writePhylip(alignment = aligned.set,
                file=paste0(output.directory, "/", save.name, ".phy"),
                interleave = F,
                strict = F)

    print(paste0("Finished ", save.name, " making IUPAC alignment!"))

  }, error = function(e) {
    warning(align.files[i], " failed: ", conditionMessage(e))
  })
  }, mc.cores = threads) #end loop

}


#########################
###### END SCRIPT
#########################


