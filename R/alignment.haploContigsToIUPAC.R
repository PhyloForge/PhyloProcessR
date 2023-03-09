#' @title haploContigsToIUPAC
#'
#' @description Function for batch trimming a folder of alignments, with the various trimming functions available to select from
#'
#' @param alignment.dir path to a folder of sequence alignments in phylip format.
#'
#' @param alignment.format available input alignment formats: fasta or phylip
#'
#' @param output.dir contigs are added into existing alignment if algorithm is "add"
#'
#' @param output.format available output formats: phylip
#'
#' @param threads path to a folder of sequence alignments in phylip format.
#'
#' @param memory give a save name if you wnat to save the summary to file.
#'
#' @param overwrite TRUE to supress mafft screen output
#'
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
  require(foreach)
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

  #Sets up multiprocessing
  cl = parallel::makeCluster(threads, outfile = "")
  doParallel::registerDoParallel(cl)
  mem.cl = floor(memory/threads)

  #Loops through each locus and does operations on them
  foreach(i=1:length(align.files), .packages = c("PhyloCap", "foreach", "Biostrings", "ape", "stringr")) %dopar% {
    #Loops through each locus and does operations on the
    #for (i in 1:length(align.files)) {

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

  }#end loop

  parallel::stopCluster(cl)

}


#########################
###### END SCRIPT
#########################


