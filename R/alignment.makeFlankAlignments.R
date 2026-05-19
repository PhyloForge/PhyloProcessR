#' @title makeFlankAlignments
#'
#' @description Extracts intron flank regions from a directory of whole-marker alignments by aligning a reference sequence (from a target fasta or a reference alignment directory) to each alignment with MAFFT, identifying the coordinates of the reference, and cutting out the flanking intronic sequence on either side. Flanks can be saved as two separate files or concatenated into a single alignment. Runs in parallel using foreach/doParallel.
#'
#' @param alignment.directory path to a folder of input sequence alignments
#'
#' @param alignment.format format of the input alignments; currently "phylip" or "fasta"
#'
#' @param output.directory path to the directory where flank alignments will be saved
#'
#' @param reference.type whether to use a "target" fasta file or an "alignment" directory as the reference for locating the exon coordinates
#'
#' @param reference.path path to the reference target fasta file (when reference.type = "target") or a directory of reference alignments (when reference.type = "alignment")
#'
#' @param target.direction if TRUE, output alignments are oriented to match the direction of the reference sequence
#'
#' @param concatenate.intron.flanks if TRUE, the left and right flanking regions are concatenated into a single alignment file; if FALSE they are saved as separate files with "_1" and "_2" suffixes
#'
#' @param threads number of parallel processing threads
#'
#' @param memory total memory in GB to allocate across threads
#'
#' @param overwrite if TRUE, overwrite existing output files; if FALSE, skip alignments already present in the output directory
#'
#' @param mafft.path system path to the MAFFT executable directory; NULL to use the system PATH
#'
#' @return saves flank alignments to output.directory in phylip format; nothing is returned to R
#'
#' @export

makeFlankAlignments = function(alignment.directory = NULL,
                                alignment.format = "phylip",
                                output.directory = NULL,
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
#   reference.type = "target"
#   reference.path = target.file
#   target.direction = TRUE
#   concatenate.intron.flanks = TRUE
#   threads = threads
#   memory = memory
#   overwrite = overwrite
#   mafft.path = mafft.path

  reference.type = match.arg(reference.type)

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

  if (length(align.files) == 0) { return("All alignments have already been completed and overwrite = FALSE.") }

  mem.cl = floor(memory/threads)

  #Loops through each locus and does operations on them
  parallel::mclapply(seq_along(align.files), function(i) {
  tryCatch({
    #Load in alignments
    if (alignment.format == "phylip"){
      align = Biostrings::DNAStringSet(Biostrings::readDNAMultipleAlignment(
        file = paste0(alignment.directory, "/", align.files[i]), format = "phylip"
      ))
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
      if (length(target.seq) == 0) {
        print(paste0(save.name, ": no matching reference found in target file — skipping."))
        return(NULL)
      }
      names(target.seq) = "Reference_Locus"

    }#end target if

    #If using the alignments
    if (reference.type == "alignment"){

      if (file.exists(paste0(reference.path, "/", align.files[i])) == FALSE){ return(NULL) }

      ref.align = Biostrings::readDNAMultipleAlignment(file = paste0(reference.path, "/", align.files[i]), format = "phylip")
      ref.align = Biostrings::DNAStringSet(ref.align)
      #Gets consensus seq for trimming more
      con.seq = PhyloProcessR::makeConsensus(ref.align, method = "majority")
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
    # Aligns and then reverses back to correction orientation
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

    #Checks for failed mafft run
    if (length(alignment) == 0){ return(NULL) }

    #Checks if you want to keep to target direction or not
    if (target.direction == TRUE){
      #Aligns and then reverses back to correction orientation
      reversed = names(alignment)[grep(pattern = "^_R_", names(alignment))]
      if (length(reversed[grep(pattern = "Reference_Locus", reversed)]) == 1){ alignment = Biostrings::reverseComplement(alignment) }
      names(alignment) = gsub(pattern = "^_R_", replacement = "", x = names(alignment))
    } else {
      #Regular fixes
      names(alignment) = gsub(pattern = "^_R_", replacement = "", x = names(alignment))
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
      PhyloProcessR::writePhylip(alignment = aligned.set,
                  file=paste0(output.directory, "/", save.name, "_1.phy"),
                  interleave = F,
                  strict = F)

      #Saves them
      write.temp = strsplit(as.character(right.intron), "")
      aligned.set = as.matrix(ape::as.DNAbin(write.temp) )

      #readies for saving
      PhyloProcessR::writePhylip(alignment = aligned.set,
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

        # readies for saving
        PhyloProcessR::writePhylip(
          alignment = aligned.set,
          file = paste0(output.directory, "/", align.files[i]),
          interleave = F,
          strict = F
        )

        #Deletes old files
        print(paste0(align.files[i], " alignment saved."))
      } else { print(paste0(align.files[i], " alignment not saved. Not enough data."))  }
    }#end

    rm(align, alignment, intron.align)
    gc()

  }, error = function(e) {
    warning(align.files[i], " failed: ", conditionMessage(e))
  })
  }, mc.cores = threads) #end i loop

} #end function
