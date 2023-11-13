#' @title trimAlignmentTargets
#'
#' @description Function to trim a set of alignments to a provided target. The function operates across a directory of alignments that correspond to a single fasta file of capture targets
#'
#' @param alignment.directory path to a folder of sequence alignments
#'
#' @param alignment.format available input alignment formats: fasta or phylip
#'
#' @param output.directory new alignment directory where the trimmed output files are saved
#'
#' @param output.format available output formats: phylip
#'
#' @param target.file path to the fasta file with the target sequences. These should be the entire marker, not the probe.
#'
#' @param target.direction TRUE ensures output alignments are the same direction as the targets
#'
#' @param min.alignment.length minimum alignment length to save in bp (default: 100)
#'
#' @param min.taxa.alignment mininum number of taxa to save alignment (default: 4)
#'
#' @param threads number of CPU threads / processes
#'
#' @param memory memory in GB
#'
#' @param overwrite TRUE to overwrite output files with the same name
#'
#' @param mafft.path system path to the mafft program
#'
#' @return a new output directory with the trimmed alignments
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

trimAlignmentTargets = function(alignment.directory = NULL,
                                alignment.format = "phylip",
                                output.directory = NULL,
                                output.format = "phylip",
                                target.file = NULL,
                                target.direction = TRUE,
                                min.alignment.length = 100,
                                min.taxa.alignment = 4,
                                threads = 1,
                                memory = 1,
                                overwrite = FALSE,
                                mafft.path = NULL) {

  # alignment.directory = "data-analysis/alignments/untrimmed_all-markers"
  # alignment.format = "phylip"
  # target.file = target.file
  # output.directory = "data-analysis/alignments/untrimmed_no-flanks"
  # output.format = "phylip"
  # min.alignment.length = 100
  # min.taxa.alignment = min.taxa.alignment
  # threads = threads
  # memory = memory
  # overwrite = overwrite
  # mafft.path = mafft.path
  # target.direction = TRUE

  require(foreach)

  #Same adds to bbmap path
  if (is.null(mafft.path) == FALSE){
    b.string = unlist(strsplit(mafft.path, ""))
    if (b.string[length(b.string)] != "/") {
      mafft.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { mafft.path = "" }

  if (alignment.directory == output.directory){ stop("You should not overwrite the original alignments.") }

 # if (dir.exists(output.dir) == FALSE) { dir.create(output.dir) }

  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
     }
  } else { dir.create(output.directory) }

  #Gathers alignments
  align.files = list.files(alignment.directory)
  target.loci = Biostrings::readDNAStringSet(file = target.file, format = "fasta")

  if (length(align.files) == 0) { stop("alignment files could not be found.") }

  #Skips files done already if resume = TRUE
  if (overwrite == FALSE){
    done.files = list.files(output.directory)
    align.files = align.files[!gsub("\\..*", "", align.files) %in% gsub("\\..*", "", done.files)]
  }

  if (length(align.files) == 0) { return("All alignments have already been completed and overwrite = FALSE.") }

  #Sets up multiprocessing
  cl = parallel::makeCluster(threads, outfile = "")
  doParallel::registerDoParallel(cl)
  mem.cl = floor(memory/threads)

  #Loops through each locus and does operations on them
  foreach::foreach(i=1:length(align.files), .packages = c("PhyloProcessR", "foreach", "Biostrings", "ape", "stringr")) %dopar% {
  #Loops through each locus and does operations on them
  for (i in 1:length(align.files)){
    #Load in alignments
    if (alignment.format == "phylip"){
      align = Biostrings::readAAMultipleAlignment(file = paste0(alignment.directory, "/", align.files[i]), format = "phylip")

      #  align = Biostrings::readDNAStringSet(file = paste0(alignment.dir, "/", align.files[i]), format = "phylip")
      #  align = readLines(paste0(alignment.dir, "/", align.files[i]))[-1]
      #  align = gsub(".*\\ ", "", align)
      #  char.count = nchar(align)

      align = Biostrings::DNAStringSet(align)
      save.name = gsub(".phy$", "", align.files[i])
      save.name = gsub(".phylip$", "", save.name)
    }#end phylip

    if (alignment.format == "fasta"){
      align = Biostrings::readDNAStringSet(paste0(alignment.directory, "/", align.files[i]) )
      save.name = gsub(".fa$", "", align.files[i])
      save.name = gsub(".fasta$", "", save.name)
    }#end phylip

    #Loads in a pulls out relevant target file
    target.seq = target.loci[names(target.loci) %in% save.name]
    if (length(target.seq) == 0) { 
      print("ALIGNMENT NOT FOUND IN TARGET MARKERS.") 
      next
    }
    names(target.seq) = "Reference_Locus"

    #Checks for correct target sequence amount
    if (length(target.seq) == 0){ stop("Locus not found in target file.")}
    if (length(target.seq) >= 2){ stop("Duplicate loci found in target file.")}

    ##############
    #STEP 2: Runs MAFFT to add
    ##############
    #Aligns and then reverses back to correction orientation
    alignment <- PhyloProcessR::runMafft(
      sequence.data = align,
      add.contigs = target.seq,
      save.name = paste0(output.directory, "/", align.files[i]),
      algorithm = "add",
      adjust.direction = TRUE,
      threads = 1,
      cleanup.files = T,
      quiet = TRUE,
      mafft.path = mafft.path
    )

    #Checks if you want to keep to target direction or not
    if (target.direction == TRUE){
      #Aligns and then reverses back to correction orientation
      reversed = names(alignment)[grep(pattern = "_R_", names(alignment))]
      if (length(reversed[grep(pattern = "Reference_Locus", reversed)]) == 1){ alignment = Biostrings::reverseComplement(alignment) }
      names(alignment) = gsub(pattern = "_R_", replacement = "", x = names(alignment))
    } else {
      #Regular fixes
      names(alignment) = gsub(pattern = "_R_", replacement = "", x = names(alignment))
    }#end if

    # ##############
    # #STEP 3: Removes exon from the intron part
    # ##############
    # #Removes the edge gaps
    ref.aligned = as.character(alignment['Reference_Locus'])
    not.gaps = stringr::str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
    ref.start = min(not.gaps)
    ref.finish = max(not.gaps)

    #Finds weird gaps to fix
    temp.gaps = as.numeric(1)
    for (k in 1:length(not.gaps)-1){ temp.gaps = append(temp.gaps, not.gaps[k+1]-not.gaps[k]) }
    temp.gaps = temp.gaps-1
    names(temp.gaps) = not.gaps
    bad.gaps = which(temp.gaps >= 30)
    front.gaps = bad.gaps[bad.gaps <= length(not.gaps) *0.10]
    end.gaps = bad.gaps[bad.gaps >= length(not.gaps) *0.90]

    #Fix big gaps if there are any
    if (length(front.gaps) != 0){
      temp.change = (max(as.numeric(names(front.gaps))-ref.start))-(max(front.gaps)-1)
      ref.start = ref.start+temp.change
    }#end gap if

    #Fix big gaps if there are any
    if (length(end.gaps) != 0){
      add.bp = length(temp.gaps)-min(end.gaps)
      #add.bp<-(ref.finish-min(as.numeric(names(end.gaps))))
      min.gaps = temp.gaps[min(end.gaps)]
      temp.change = as.numeric(names(min.gaps))-as.numeric(min.gaps)
      ref.finish = temp.change+add.bp
    }#end gap if

    target.region = Biostrings::subseq(alignment, ref.start, ref.finish)
    target.region = target.region[names(target.region) != "Reference_Locus"]

    ##############
    #STEP 5: Cleanup and save
    ##############
    #removes loci with too few taxa
    skip.file = FALSE
    if (length(names(target.region)) < as.numeric(min.taxa.alignment)){
      print(paste0(align.files[i], " deleted. Too few taxa after trimming.") )
      skip.file = TRUE
    }

    #removes too short loci
    if (Biostrings::width(target.region)[1] < as.numeric(min.alignment.length)){
      print(paste(align.files[i], " deleted. Trimmed alignment length below threshold.") )
      skip.file = TRUE
    }

    if (skip.file != TRUE){
      #string splitting
      write.temp = strsplit(as.character(target.region), "")
      aligned.set = as.matrix(ape::as.DNAbin(write.temp) )

      #readies for saving
      PhyloProcessR::writePhylip(
        alignment = aligned.set,
        file = paste0(output.directory, "/", align.files[i]),
        interleave = F,
        strict = F
      )

      #Deletes old files
      print(paste0(align.files[i], " alignment saved."))
    } else { print(paste0(align.files[i], " alignment discarded. Not enough data to save.")) }

    rm()
    gc()

  } #end i loop

  parallel::stopCluster(cl)

} #end function
