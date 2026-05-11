#' @title mergeDuplicates
#'
#' @description Processes a directory of alignments to handle samples that appear more than once (duplicate sequence names). When method is "merge-sample", duplicate sequences for the same sample are merged column by column: identical characters are kept, gap/N/? characters are resolved to the non-gap character when possible, and samples with irresolvable conflicts between their duplicate sequences are dropped. When method is "delete-shortest", the shorter duplicate is removed (not yet fully implemented). Runs in parallel using foreach/doParallel.
#'
#' @param alignment.directory path to a folder of input sequence alignments
#'
#' @param alignment.format format of the input alignments; "phylip" or "fasta"
#'
#' @param output.directory path to the directory where merged alignments will be saved
#'
#' @param output.format format for the output alignments; currently "phylip"
#'
#' @param method how to handle duplicate sample names: "merge-sample" to merge duplicate sequences column by column, or "delete-shortest" to remove the shorter duplicate
#'
#' @param threads number of parallel processing threads
#'
#' @param memory total memory in GB to allocate across threads
#'
#' @param overwrite if TRUE, delete and recreate the output directory; if FALSE, stop if the directory already exists
#'
#' @param quiet if TRUE, suppress progress messages
#'
#' @param mafft.path system path to the MAFFT executable directory; NULL to use the system PATH
#'
#' @return saves merged alignments to output.directory in phylip format; nothing is returned to R
#'
#' @export

mergeDuplicates = function(alignment.directory = NULL,
                           alignment.format = "phylip",
                           output.directory = NULL,
                           output.format = "phylip",
                           method = c("delete-shortest", "merge-sample"),
                           threads = 1,
                           memory = 1,
                           overwrite = FALSE,
                           quiet = FALSE,
                           mafft.path = NULL) {

  # setwd("/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/New-Guinea_Frogs")
  # alignment.directory = "/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/New-Guinea_Frogs/alignments/trimmed_all-markers"
  # alignment.format = "phylip"
  # output.directory = "alignments/rm-dup-trimmed_all-markers"
  # output.format = "phylip"
  # method = "merge-sample"
  # threads = 4
  # memory = 8
  # overwrite = TRUE
  # quiet = TRUE
  # mafft.path = "/usr/local/bin"

  # *** combine non-overlapping seqs
  # *** check for duplicates and remove
  # *** include loci not in sequence capture too
  # *** trim to legacy

  #Same adds to bbmap path
  if (is.null(mafft.path) == FALSE){
    b.string = unlist(strsplit(mafft.path, ""))
    if (b.string[length(b.string)] != "/") {
      mafft.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { mafft.path = "" }

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
  new.align = c()

  #Sets up multiprocessing
  cl = parallel::makeCluster(threads, outfile = "")
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
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

    ###################################################################################
    if (method == "delete-shortest"){


    }

    ###################################################################################
    # Duplication changes and such
    if (method == "merge-sample"){

      dup.names = names(old.align)[duplicated(names(old.align)) == TRUE]

      if (length(dup.names) == 0){
        #save file
        #Saves them
        write.temp = strsplit(as.character(old.align), "")
        aligned.set = as.matrix(ape::as.DNAbin(write.temp) )

        #readies for saving
        writePhylip(alignment = aligned.set,
                    file=paste0(output.directory, "/", save.name, ".phy"),
                    interleave = F,
                    strict = F)

        print(paste0("Finished ", save.name, " duplicate merging successful!"))
      }

      if (length(dup.names) != 0){

        nodup.align = old.align[!names(old.align) %in% dup.names]

        #Goes through each duplicate
        save.seqs = Biostrings::DNAStringSet()
        no.del = 0
        for (j in 1:length(dup.names)){

          dup.sample = old.align[names(old.align) %in% dup.names[j]]
          #Converts alignment to matrix of characters to be used
          new.align = strsplit(as.character(dup.sample), "")
          align.in = matrix(unlist(new.align), ncol = length(new.align[[1]]), byrow = T)

          save.string = as.character()
          skip.taxa = FALSE
          for (k in 1:ncol(align.in)){

            temp.col = align.in[,k]

            if (temp.col[1] == temp.col[2]){ save.char = temp.col[1] }
            if (temp.col[1] != temp.col[2]){
              temp.col = temp.col[temp.col != "-"]
              temp.col = temp.col[temp.col != "N"]
              temp.col = temp.col[temp.col != "?"]

              if (length(temp.col) == 0){ save.char = "-" }
              if (length(temp.col) == 1){ save.char = temp.col }
              if (length(temp.col) == 2){
                no.del = no.del + 1
                skip.taxa = TRUE
                break
              }

              # temp.col = temp.col[temp.col != "Y"]
              # temp.col = temp.col[temp.col != "B"]
              # temp.col = temp.col[temp.col != "D"]
              # temp.col = temp.col[temp.col != "H"]
              # temp.col = temp.col[temp.col != "V"]
              # temp.col = temp.col[temp.col != "K"]
              # temp.col = temp.col[temp.col != "S"]
              # temp.col = temp.col[temp.col != "W"]
              # temp.col = temp.col[temp.col != "R"]
              # temp.col = temp.col[temp.col != "M"]

            }#end upper if

            save.string = append(save.string, save.char)

          }#end k loop

          if (skip.taxa == FALSE){
            out.consensus = Biostrings::DNAStringSet(paste0(save.string, collapse = ""))
            names(out.consensus) = dup.names[j]
            save.seqs = append(save.seqs, out.consensus)
          }#end false

        }#end j

        new.align = append(save.seqs, nodup.align)

        #Saves them
        write.temp = strsplit(as.character(new.align), "")
        aligned.set = as.matrix(ape::as.DNAbin(write.temp) )

        #readies for saving
        writePhylip(alignment = aligned.set,
                    file=paste0(output.directory, "/", save.name, ".phy"),
                    interleave = F,
                    strict = F)

        print(paste0("Finished ", save.name, " duplicate merging successfully!"))
        if (no.del != 0){
          print(paste0("***Removed ", no.del, " non-matching duplicates from ", save.name, "."))
        }


      }#end if to do the thing

    }#end combine.same.sample if

    ####################################################################################

  }#end loop

  parallel::stopCluster(cl)


}#end fuction

#END SCRIPT
