#' @title alignMatchingTargets
#'
#' @description Function for gather summary statistics on your alignment. Can be used for filtering or summarizing data.
#'
#' @param matched.targets.file path to a folder of sequence alignments in phylip format.
#'
#' @param subset.start contigs are added into existing alignment if algorithm is "add"
#'
#' @param subset.end contigs are added into existing alignment if algorithm is "add
#'
#' @param output.dir contigs are added into existing alignment if algorithm is "add"
#'
#' @param algorithm algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param adjust.direction TRUE applies the adjust sequence direction function of MAFFT
#'
#' @param remove.reverse.tag if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param threads path to a folder of sequence alignments in phylip format.
#'
#' @param cleanup.files give a save name if you wnat to save the summary to file.
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

alignMatchingTargets = function(matched.targets.file = NULL,
                                subset.start = "all",
                                subset.end = "all",
                                output.dir = NULL,
                                algorithm = "localpair",
                                adjust.direction = TRUE,
                                remove.reverse.tag = TRUE,
                                threads = 1,
                                overwrite = FALSE,
                                resume = TRUE,
                                quiet = TRUE) {

  #Debug
  # matched.targets.file = "/Volumes/Rodents/Murinae/Data_Processing/alignment_fastas/Emily-Data_iupac-no-trim-noncoding_contigs.fa"
  # output.dir = "full-dataset/coding-untrimmed"
  # dir.create("full-dataset")
  # remove.reverse.tag = TRUE
  # subset.start = "all"
  # subset.end = "all"
  # overwrite = FALSE
  # resume = TRUE
  # quiet = TRUE
  # adjust.direction = TRUE
  # algorithm = "localpair"

  #Checks and creates output directory
  if (dir.exists(output.dir) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.dir))
      dir.create(output.dir)
    }
  } else { dir.create(output.dir) }

  #Load sample file from probe matching step
  all.data = Biostrings::readDNAStringSet(matched.targets.file)
  locus.names = unique(gsub("_\\|_.*", "", names(all.data)))

  #Skips files done already if resume = TRUE
  if (resume == TRUE){
    done.files = list.files(output.dir)
    locus.names = locus.names[!locus.names %in% gsub(".phy$", "", done.files)]
  }

  #Sets and checks bounds
  if (subset.end == "all") {
    subset.start  = 1
    subset.end = length(locus.names)
  }#end all
  #if the length is longer than the number of loci
  if (subset.end > length(locus.names)){ subset.end = length(locus.names) }

  #Loops through each locus and writes each species to end of file
  for (i in subset.start:subset.end) {

    #Match probe names to contig names to acquire data
    match.data = all.data[grep(pattern = paste0(locus.names[i], "_"), x = names(all.data))]

    ##############
    #STEP 1: Throw out loci if there are too few taxa
    ##############
    if (length(names(match.data)) <= 2){
      print(paste0(locus.names[i], " had too few taxa."))
      next
    }

    ##############
    #STEP 2: Sets up fasta for aligning
    ##############
    names(match.data) = gsub(pattern = ".*_\\|_", replacement = "", x = names(match.data))
    #Gets reference locus
    final.loci = match.data

    ##############
    #STEP 3: Runs MAFFT to align
    ##############
    #Aligns and then reverses back to correction orientation
    alignment = runMafft(sequence.data = final.loci,
                         save.name = locus.names[i],
                         threads = threads,
                         algorithm = algorithm,
                         adjust.direction = adjust.direction,
                         cleanup.files = FALSE)

    #Checks for failed mafft run
    if (length(alignment) == 0){
      print(paste0(locus.names[i], " did not successfully align."))
      next }

    #Removes the reverse name
    if (remove.reverse.tag == TRUE){
      names(alignment) = gsub(pattern = "^_R_", replacement = "", x = names(alignment))
    }#end if

    #If no alignment assessing is done, saves
    write.temp = strsplit(as.character(alignment), "")
    aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
    #readies for saving
    writePhylip(aligned.set, file= paste0(output.dir, "/", locus.names[i], ".phy"), interleave = F)

    #Deletes old files
    system(paste0("rm ", locus.names[i], "_align.fa "))
    print(paste0(locus.names[i], " alignment saved."))

  }# end big i loop

}#end function
