#' @title mergeDuplicates
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

mergeDuplicates = function(alignment.directory = NULL,
                           alignment.format = "phylip",
                           output.directory = NULL,
                           output.format = "phylip",
                           method = c("delete-shortest", "merge-sample"),
                           include.all.together = FALSE,
                           threads = 1,
                           memory = 1,
                           overwrite = FALSE,
                           quiet = FALSE,
                           mafft.path = NULL) {

  # setwd("/Volumes/Armored/Brygomantis")
  # library(foreach)
  # alignment.directory = "data-processing/alignments-3/trimmed_all-markers"
  # alignment.format = "phylip"
  # output.directory = "data-processing/alignments-3/rm-dup-trimmed_all-markers"
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
        for (j in 1:length(dup.names)){

          dup.sample = old.align[names(old.align) %in% dup.names[j]]
          #Converts alignment to matrix of characters to be used
          new.align = strsplit(as.character(dup.sample), "")
          align.in = matrix(unlist(new.align), ncol = length(new.align[[1]]), byrow = T)

          save.string = as.character()
          for (k in 1:ncol(align.in)){

            temp.col = align.in[,k]

            if (temp.col[1] == temp.col[2]){ save.char = temp.col[1] }
            if (temp.col[1] != temp.col[2]){
              temp.col = temp.col[temp.col != "-"]
              temp.col = temp.col[temp.col != "N"]
              temp.col = temp.col[temp.col != "?"]

              if (length(temp.col) == 0){ save.char = "-" }
              if (length(temp.col) == 1){ save.char = temp.col }
              if (length(temp.col) == 2){ save.char = temp.col[1] }

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

          out.consensus = Biostrings::DNAStringSet(paste0(save.string, collapse = ""))
          names(out.consensus) = dup.names[j]
          save.seqs = append(save.seqs, out.consensus)

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

      }#end if to do the thing

    }#end combine.same.sample if

    ####################################################################################

  }#end loop

  parallel::stopCluster(cl)


}#end fuction

#END SCRIPT
