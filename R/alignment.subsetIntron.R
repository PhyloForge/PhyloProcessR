#' @title makeIntronAlignments
#'
#' @description Function for batch trimming a folder of alignments, with the various trimming functions available to select from
#'
#' @param alignment.directory path to a folder of sequence alignments
#'
#' @param alignment.format available input alignment formats: fasta or phylip
#'
#' @param output.directory new alignment directory where the trimmed output files are saved
#'
#' @param output.format available output formats: phylip
#'
#' @param HmmCleaner algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param HmmCleaner.path TRUE applies the adjust sequence direction function of MAFFT
#'
#' @param TrimAl if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param TrimAl.path path to a folder of sequence alignments in phylip format.
#'
#' @param trim.external give a save name if you wnat to save the summary to file.
#'
#' @param min.external.percent TRUE to supress mafft screen output
#'
#' @param trim.coverage path to a folder of sequence alignments in phylip format.
#'
#' @param min.coverage.percent contigs are added into existing alignment if algorithm is "add"
#'
#' @param trim.column algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param min.column.gap.percent TRUE applies the adjust sequence direction function of MAFFT
#'
#' @param alignment.assess if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param min.sample.bp path to a folder of sequence alignments in phylip format.
#'
#' @param min.alignment.length give a save name if you wnat to save the summary to file.
#'
#' @param min.taxa.alignment TRUE to supress mafft screen output
#'
#' @param min.gap.percent if a file name is provided, save.name will be used to save aligment to file as a fasta
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

makeIntronAlignments = function(alignment.directory = NULL,
                                alignment.format = "phylip",
                                output.directory = NULL,
                                output.format = "phylip",
                                reference.type = c("target", "alignment"),
                                reference.path = NULL,
                                target.direction = TRUE,
                                concatenate.intron.flanks = TRUE,
                                threads = 1,
                                memory = 1,
                                overwrite = FALSE,
                                mafft.path = NULL) {

#
#   alignment.directory = "alignments/untrimmed_all-markers"
#   alignment.format = "phylip"
#   output.directory = "alignments/untrimmed_introns"
#   output.format = "phylip"
#   reference.type = "target"
#   reference.path = target.file
#   target.direction = TRUE
#   concatenate.intron.flanks = TRUE
#   threads = threads
#   memory = memory
#   overwrite = overwrite
#   mafft.path = mafft.path


  #Same adds to bbmap path
  if (is.null(mafft.path) == FALSE){
    b.string = unlist(strsplit(mafft.path, ""))
    if (b.string[length(b.string)] != "/") {
      mafft.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { mafft.path = "" }

  if (alignment.directory == output.directory){ stop("You should not overwrite the original alignments.") }

  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } else { dir.create(output.directory) }

  #Gathers alignments
  align.files = list.files(alignment.directory)

  if (reference.type == "target"){
    target.loci = Biostrings::readDNAStringSet(file = reference.path, format = "fasta")
  }#end if

  if (reference.type == "alignment"){
    ref.align = list.files(reference.path)
  }#end if

  if (length(align.files) == 0) { stop("alignment files could not be found.") }

  #Skips files done already if resume = TRUE
  if (overwrite == FALSE){
    done.files = list.files(output.directory)
    align.files = align.files[!gsub("\\..*", "", align.files) %in% gsub("\\..*", "", done.files)]
  }

  if (length(align.files) == 0) { stop("All alignments have already been completed and overwrite = FALSE.") }

  #Sets up multiprocessing
  cl = parallel::makeCluster(threads, outfile = "")
  doParallel::registerDoParallel(cl)
  mem.cl = floor(memory/threads)

  #Loops through each locus and does operations on them
  foreach::foreach(i=1:length(align.files), .packages = c("PhyloCap", "foreach", "Biostrings", "ape", "stringr")) %dopar% {
  #Loops through each locus and does operations on them
  #for (i in 1:length(align.files)) {
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

    if (reference.type == "target"){
      #Loads in a pulls out relevant target file
      target.seq = target.loci[names(target.loci) %in% save.name]
      names(target.seq) = "Reference_Locus"

    }#end target if

    #If using the alignments
    if (reference.type == "alignment"){

      if (file.exists(paste0(reference.path, "/", align.files[i])) == FALSE){ return(NULL) }

      ref.align = Biostrings::readAAMultipleAlignment(file = paste0(reference.path, "/", align.files[i]), format = "phylip")
      ref.align = Biostrings::DNAStringSet(ref.align)
      #Gets consensus seq for trimming more
      con.seq = makeConsensus(ref.align, method = "majority")
      #Removes the edge gaps
      target.seq = Biostrings::DNAStringSet(gsub("\\+|-", "", as.character(con.seq)))
      names(target.seq) = "Reference_Locus"
    }#end alignment if


    #Checks for correct target sequence amount
    if (length(target.seq) == 0){ return(NULL) }
    if (length(target.seq) >= 2){ return(NULL) }
    if (Biostrings::width(target.seq) <= 10) { return(NULL) }

    ##############
    #STEP 2: Runs MAFFT to add
    ##############
    #Aligns and then reverses back to correction orientation
    alignment = runMafft(sequence.data = align,
                         add.contigs = target.seq,
                         save.name = paste0(output.directory, "/", align.files[i]),
                         algorithm = "add",
                         adjust.direction = TRUE,
                         threads = 1,
                         cleanup.files = T,
                         quiet = TRUE,
                         mafft.path = mafft.path)

    #Checks for failed mafft run
    if (length(alignment) == 0){ return(NULL) }

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

    #Cuts out the intron pieces
    intron.left = Biostrings::subseq(alignment, 1, ref.start-1)
    intron.right = Biostrings::subseq(alignment, ref.finish+1, Biostrings::width(alignment))
    save.names  = names(alignment)

    #saves intron flanks separately
    if (concatenate.intron.flanks == FALSE) {
      #Remove gap only alignments
      gap.align = strsplit(as.character(intron.left), "")
      gap.count = unlist(lapply(gap.align, function(x) length(x[x != "-"]) ) )
      gap.rem = gap.count[gap.count <= 10]
      left.intron = intron.left[!names(intron.left) %in% names(gap.rem)]

      #Remove gap only alignments
      gap.align = strsplit(as.character(intron.right), "")
      gap.count = unlist(lapply(gap.align, function(x) length(x[x != "-"]) ) )
      gap.rem = gap.count[gap.count <= 10]
      right.intron = intron.right[!names(intron.right) %in% names(gap.rem)]

      #Saves them
      write.temp = strsplit(as.character(left.intron), "")
      aligned.set = as.matrix(ape::as.DNAbin(write.temp) )

      save.name = gsub(".phy$", "", align.files[i])

      #readies for saving
      writePhylip(alignment = aligned.set,
                  file=paste0(output.directory, "/", save.name, "_1.phy"),
                  interleave = F,
                  strict = F)

      #Saves them
      write.temp = strsplit(as.character(right.intron), "")
      aligned.set = as.matrix(ape::as.DNAbin(write.temp) )

      #readies for saving
      writePhylip(alignment = aligned.set,
                  file=paste0(output.directory, "/", save.name, "_2.phy"),
                  interleave = F,
                  strict = F)

      #Deletes old files
      print(paste0(align.files[i], " alignment saved."))

    } else {

      #Merges the alignments
      intron.align = Biostrings::DNAStringSet(paste0(as.character(intron.left), as.character(intron.right)))
      names(intron.align) = save.names
      intron.align = intron.align[names(intron.align) != "Reference_Locus"]

      #Remove gap only alignments
      gap.align = strsplit(as.character(intron.align), "")
      gap.count = unlist(lapply(gap.align, function(x) length(x[x != "-"]) ) )
      gap.rem = gap.count[gap.count <= 10]
      rem.align = intron.align[!names(intron.align) %in% names(gap.rem)]

      ##############
      #STEP 5: Cleanup and save
      ##############

      if (length(rem.align) != 0 && length(Biostrings::width(rem.align)) != 0){
        #string splitting
        write.temp = strsplit(as.character(rem.align), "")
        aligned.set = as.matrix(ape::as.DNAbin(write.temp) )

        #readies for saving
        writePhylip(alignment = aligned.set,
                    file=paste0(output.directory, "/", align.files[i]),
                    interleave = F,
                    strict = F)

        #Deletes old files
        print(paste0(align.files[i], " alignment saved."))
      } else { print(paste0(align.files[i], " alignment not saved. Not enough data."))  }
    }#end

    rm()
    gc()

  }#end i loop

  #close multithread
  parallel::stopCluster(cl)

} #end function
