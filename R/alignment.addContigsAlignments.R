#' @title addTaxaAlignments
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

addContigsAlignments = function(alignment.directory = NULL,
                                alignment.format = "phylip",
                                output.directory = NULL,
                                output.format = "phylip",
                                sample.contigs = NULL,
                                copy.all = TRUE,
                                duplicate.handling = "merge",
                                threads = 1,
                                memory = 1,
                                overwrite = FALSE,
                                quiet = FALSE,
                                mafft.path = NULL,
                                blast.path = NULL) {

  # library(foreach)
  # library(PhyloCap)
  # setwd("/Volumes/LaCie/Mantellidae_Subfamily/Transcriptomes_new")
  # alignment.directory = "Alignments/original"
  # alignment.format = "phylip"
  # output.directory = "Alignments/new_alignments"
  # output.format = "phylip"
  # sample.contigs = "contigs"
  # duplicate.handling = "merge"
  # threads = 8
  # memory = 24
  # overwrite = FALSE
  # quiet = FALSE
  # copy.all = TRUE
  # mafft.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin"
  # blast.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin"
  # bbmap.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin"

  #Adds slash to path
  if (is.null(mafft.path) == FALSE){
    b.string = unlist(strsplit(mafft.path, ""))
    if (b.string[length(b.string)] != "/") {
      mafft.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { mafft.path = "" }

  if (is.null(blast.path) == FALSE){
    b.string = unlist(strsplit(blast.path, ""))
    if (b.string[length(b.string)] != "/") {
      blast.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { blast.path = "" }

  if (is.null(bbmap.path) == FALSE){
    b.string = unlist(strsplit(bbmap.path, ""))
    if (b.string[length(b.string)] != "/") {
      bbmap.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { bbmap.path = "" }

  if (alignment.directory == output.directory){ stop("You should not overwrite the original alignments.") }

  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } else { dir.create(output.directory) }

  #Gathers alignments
  align.files = list.files(alignment.directory)
  if (length(align.files) == 0) { stop("alignment files could not be found.") }

  #Make consensus sequences from alignments
  align.consensus = Biostrings::DNAStringSet()
  for (i in 1:length(align.files)){

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

    con.seq = makeConsensus(alignment = align)
    names(con.seq) = save.name
    align.consensus = append(align.consensus, con.seq)

  }

  #Finds probes that match to two or more contigs
  final.loci = as.list(as.character(align.consensus))
  writeFasta(sequences = final.loci, names = names(final.loci),
             paste0("consensus_alignments.fa"), nbchar = 1000000, as.string = T)

  #Go thorugh each contig and blast to to consensus sequences
  contig.files = list.files(sample.contigs)

  #Uses the match targets function
  matchTargets(assembly.directory = sample.contigs,
               target.file = "consensus_alignments.fa",
               alignment.contig.name = "add_contigs",
               output.directory = "match-targets",
               min.match.percent = 50,
               min.match.length = 40,
               min.match.coverage = 35,
               trim.target = FALSE,
               threads = threads,
               memory = memory,
               overwrite = TRUE,
               resume = TRUE,
               quiet = TRUE,
               blast.path = blast.path,
               bbmap.path = bbmap.path)

  # Go through each alignment and add the sample in

  ### Reads in the additional stuff
  add.taxa = Biostrings::readDNAStringSet("add_contigs_to-align.fa")

  #Sets up multiprocessing
  cl = parallel::makeCluster(threads, outfile = "")
  doParallel::registerDoParallel(cl)
  mem.cl = floor(memory/threads)

  #Loops through each locus and does operations on them
  foreach(i=1:length(align.files), .packages = c("PhyloCap", "foreach", "Biostrings", "ape", "stringr")) %dopar% {
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

      ##############
      #STEP 1: Sets up new for adding into alignment
      ##############
      #Match probe names to contig names to acquire data
      add.seqs = add.taxa[gsub("_\\|_.*", "", names(add.taxa)) == save.name]
      #add.seqs = add.taxa[grep(pattern = paste0(name.locus, "_"), x = names(add.taxa))]
      names(add.seqs) = gsub(pattern = ".*_\\|_", replacement = "", x = names(add.seqs))

      if (length(add.seqs) == 0){

        if (copy.all == TRUE){
          #Saves prelim exon file
          new.align = strsplit(as.character(align), "")
          aligned.set = as.matrix(ape::as.DNAbin(new.align))

          #readies for saving
          writePhylip(alignment = aligned.set,
                      file=paste0(output.directory, "/", align.files[i]),
                      interleave = F,
                      strict = F)
        }#end

        return(NULL)
      }

      ##############
      #STEP 2: load up and align additional samples to old alignment
      ##############
      #Align
      alignment = runMafft(sequence.data = align,
                           add.contigs = add.seqs,
                           adjust.direction = T,
                           algorithm = "add",
                           save.name = save.name,
                           cleanup.files = TRUE,
                           quiet = TRUE,
                           mafft.path = mafft.path)

      reversed = names(alignment)[grep(pattern = "_R_", names(alignment))]
      if (length(reversed[grep(pattern = "Reference_Locus", reversed)]) == 1){
        alignment = reverseComplement(alignment)
      }

      names(alignment) = gsub(pattern = "_R_", replacement = "", x = names(alignment))

      dup.taxa = names(alignment)[duplicated(names(alignment))]

      if (length(dup.taxa) != 0){

        #To merge duplicates in case extra data was added in or deletes the shortest
        if (duplicate.handling == "merge"){

          nodup.align = alignment[!names(alignment) %in% dup.taxa]

          #Goes through each duplicate
          save.seqs = Biostrings::DNAStringSet()
          for (j in 1:length(dup.taxa)){

            dup.sample = alignment[names(alignment) %in% dup.taxa[j]]
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

              }#end upper if

              save.string = append(save.string, save.char)

            }#end k loop

            out.consensus = Biostrings::DNAStringSet(paste0(save.string, collapse = ""))
            names(out.consensus) = dup.taxa[j]
            save.seqs = append(save.seqs, out.consensus)

          }#end j loop

          dup.alignment = append(nodup.align, save.seqs)

        }# if statement merge


        if (duplicate.handling == "longest"){

          nodup.align = alignment[!names(alignment) %in% dup.taxa]

          #Goes through each duplicate
          save.seqs = Biostrings::DNAStringSet()
          for (j in 1:length(dup.taxa)){

            dup.sample = alignment[names(alignment) %in% dup.taxa[j]]
            #Converts alignment to matrix of characters to be used
            new.align = strsplit(as.character(dup.sample), "")
            align.in = matrix(unlist(new.align), ncol = length(new.align[[1]]), byrow = T)

            gap.align = lapply(new.align, function(x) gsub("\\?|N|n", "-", x) )
            base.count = unlist(lapply(gap.align, function(x) length(x[x != "-"]) ) )

            taxa.keep = base.count[base.count == max(base.count)][1]
            out.sample = dup.sample[names(dup.sample) %in% names(taxa.keep)][1]
            names(out.sample) = dup.taxa[j]
            save.seqs = append(save.seqs, out.sample)

          }#end j loop

          dup.alignment = append(nodup.align, save.seqs)

        }

        if (duplicate.handling == "original"){

          nodup.align = alignment[!names(alignment) %in% dup.taxa]

          #Goes through each duplicate
          save.seqs = Biostrings::DNAStringSet()
          for (j in 1:length(dup.taxa)){

            dup.sample = alignment[names(alignment) %in% dup.taxa[j]]
            out.sample = dup.sample[1]
            names(out.sample) = dup.taxa[j]
            save.seqs = append(save.seqs, out.sample)

          }#end j loop

          dup.alignment = append(nodup.align, save.seqs)


        }
      }#end if

      if (length(dup.taxa) == 0){ dup.alignment = alignment }

      #Saves prelim exon file
      new.align = strsplit(as.character(dup.alignment), "")
      aligned.set = as.matrix(ape::as.DNAbin(new.align))

      #readies for saving
      writePhylip(alignment = aligned.set,
                  file=paste0(output.directory, "/", align.files[i]),
                  interleave = F,
                  strict = F)

      #Deletes old files
     # print(paste0(align.files[i], " alignment saved."))

    }#end i loop

    #close multithread
    parallel::stopCluster(cl)

} #end function
