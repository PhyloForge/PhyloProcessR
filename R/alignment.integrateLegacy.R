#' @title integrateLegacy
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

integrateLegacy = function(alignment.directory = NULL,
                           alignment.format = "phylip",
                           output.directory = NULL,
                           output.format = "phylip",
                           legacy.directory = NULL,
                           legacy.format = "phylip",
                           target.markers = NULL,
                           combine.same.sample = TRUE,
                           include.uncaptured.legacy = FALSE,
                           include.all.together = FALSE,
                           threads = 1,
                           memory = 1,
                           overwrite = FALSE,
                           quiet = FALSE,
                           mafft.path = NULL,
                           blast.path = NULL) {

#
#   alignment.directory = "/Volumes/LaCie/mitocap_2/Alignments/untrimmed-alignments"
#   alignment.format = "phylip"
#   output.directory = "/Volumes/LaCie//mitocap_2/untrimmed_mt_legacy"
#   output.format = "phylip"
#   legacy.directory = "/Volumes/LaCie//mitocap_2/genbank-legacy"
#   legacy.format = "phylip"
#   target.markers = "/Volumes/LaCie/mitocap_2/reference/refMarkers.fa"
#   combine.same.sample = FALSE
#   include.uncaptured.legacy = FALSE
#   include.all.together = TRUE
#   threads = 10
#   memory = 50
#   overwrite = TRUE
#   quiet = FALSE
#   mafft.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin"
#   blast.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin"

  # setwd("/Volumes/Armored/Boophis_Tree")
  # alignment.directory = "alignments/untrimmed_all-markers"
  # alignment.format = "phylip"
  # output.directory = "alignments/untrimmed_legacy"
  # output.format = "phylip"
  # legacy.directory = "alignments/legacy-genbank-nuc"
  # legacy.format = "phylip"
  # target.markers = "Hyloidea_All-Markers_Apr21-2019.fa"
  # combine.same.sample = TRUE
  # include.uncaptured.legacy = TRUE
  # include.all.together = TRUE
  # threads = 1
  # memory = 4
  # overwrite = TRUE
  # quiet = TRUE
  # mafft.path = "/usr/local/bin"
  # blast.path = "/Users/chutter/miniconda3/bin"

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

  #Same adds to bbmap path
  if (is.null(blast.path) == FALSE){
    b.string = unlist(strsplit(blast.path, ""))
    if (b.string[length(b.string)] != "/") {
      blast.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { blast.path = "" }

  #Checks this
  if (alignment.directory == output.directory){ stop("You should not overwrite the original alignments.") }

  #Overwrite
  if (dir.exists(paste0(output.directory, "-all")) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory, "-all"))
      dir.create(paste0(output.directory, "-all"))
      system(paste0("rm -r ", output.directory, "-only"))
      dir.create(paste0(output.directory, "-only"))
    } else {
      stop("Overwrite = FALSE and directory exists. Either change to TRUE or overwrite manually.")
    }
  } else {
    dir.create(paste0(output.directory, "-all"))
    dir.create(paste0(output.directory, "-only"))
  }#end else

  #Gathers alignments
  align.files = list.files(alignment.directory)
  legacy.files = list.files(legacy.directory)

  if (length(align.files) == 0) { stop("alignment files could not be found.") }
  if (length(legacy.files) == 0) { stop("alignment legacy files could not be found.") }

  #headers for the blast db
  headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
              "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

  #Make blast database for the probe loci
  system(paste0(blast.path, "makeblastdb -in ", target.markers,
                " -parse_seqids -dbtype nucl -out target_nucl-blast_db"), ignore.stdout = quiet)

  ### Reads in the additional stuff
  #add.taxa = Biostrings::readDNAStringSet(sample.markers)
  bait.loci = Biostrings::readDNAStringSet(target.markers)  # loads up fasta file

  #Sets up multiprocessing
  #cl = parallel::makeCluster(threads, outfile = "")
  #doParallel::registerDoParallel(cl)
  #mem.cl = floor(memory/threads)

  #Loops through each locus and does operations on them
  #foreach(i=1:length(legacy.files), .packages = c("PhyloCap", "foreach", "Biostrings", "ape", "stringr")) %dopar% {
  #Loops through each locus and does operations on the
  delete.old = c()
  for (i in 1:length(legacy.files)) {

    ##############
    #STEP 0: Load in legacy alignment
    ##############
    #Load in alignments
    if (legacy.format == "phylip"){
      align = Biostrings::readAAMultipleAlignment(file = paste0(legacy.directory, "/", legacy.files[i]), format = "phylip")
      align = Biostrings::DNAStringSet(align)
      save.name = gsub(".phy$", "", legacy.files[i])
      save.name = gsub(".phylip$", "", save.name)
    }#end phylip

    if (legacy.format == "fasta"){
      align = Biostrings::readDNAStringSet(paste0(legacy.directory, "/", legacy.files[i]) )
      save.name = gsub(".fa$", "", legacy.files[i])
      save.name = gsub(".fasta$", "", save.name)
    }#end phylip

    # ##############
    # #STEP 1: Sets up new for adding into alignment
    # ##############
    # #Match probe names to contig names to acquire data
    # name.locus = gsub(".phy$", "", locus.names[i])
    # add.seqs = add.taxa[grep(pattern = paste0(name.locus, "_"), x = names(add.taxa))]
    # names(add.seqs) = gsub(pattern = ".*_\\|_", replacement = "", x = names(add.seqs))
    #
    # if (length(add.seqs) == 0){ return(NULL) }
    #
    # #Gets reference locus
    # ref.locus = bait.loci[grep(pattern = paste0(name.locus, "$"), x = names(bait.loci))]
    # names(ref.locus) = paste("Reference_Locus")
    # add.seqs = append(add.seqs, ref.locus)

    ##############
    #STEP 1: Blast to targets
    ##############

    #Create consensus to blast
    con.seq = makeConsensus(alignment = align,
                  method = c("majority"),
                  warn.non.IUPAC = FALSE,
                  remove.gaps = TRUE,
                  type = c("DNA"))

    #Finds probes that match to two or more contigs
    final.loci = as.list(as.character(con.seq))
    writeFasta(sequences = final.loci, names = names(final.loci),
               paste0(save.name, "_consensus.fa"), nbchar = 1000000, as.string = T)

    #Matches samples to loci
    system(paste0(blast.path, "blastn -task dc-megablast -db target_nucl-blast_db -evalue 0.001",
                  " -query ", save.name, "_consensus.fa -out ", save.name, "_target-blast-match.txt",
                  " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
                  " -num_threads ", threads))

    #Loads in match data
    match.data = data.table::fread(paste0(save.name, "_target-blast-match.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
    data.table::setnames(match.data, headers)

    if (nrow(match.data) == 0) {
      print(paste0(sample, " had no matches. Skipping"))
      next }

    if (nrow(match.data) >= 2){
      match.data = match.data[match.data$bitscore == max(match.data$bitscore),]
    }

    if (nrow(match.data) >= 2){
      stop("too many matches...")
    }

    found.align = align.files[grep(match.data$tName, align.files)]

    if (include.uncaptured.legacy == FALSE){
      if (length(found.align) == 0){
        print(paste0(save.name, " not found in the alignments. Moving to next."))
        next
      }
    }#end if

    if (include.uncaptured.legacy == TRUE){
      if (length(found.align) == 0){
        #Saves them
        write.temp = strsplit(as.character(align), "")
        aligned.set = as.matrix(ape::as.DNAbin(write.temp) )

        #readies for saving
        writePhylip(alignment = aligned.set,
                    file=paste0(output.directory, "-only/", save.name, ".phy"),
                    interleave = F,
                    strict = F)

        print(paste0("No sequence capture alignment found. Included uncaptured ", save.name, "  successfully!"))
        next
      }
    }#end if

    ##############
    #STEP 2: Runs MAFFT to add
    ##############
    #Load in alignments
    if (alignment.format == "phylip"){
      old.align = Biostrings::readAAMultipleAlignment(file = paste0(alignment.directory, "/", found.align), format = "phylip")
      old.align = Biostrings::DNAStringSet(old.align)
      found.name = gsub(".phy$", "", found.align)
      found.name = gsub(".phylip$", "", found.name)
    }#end phylip

    if (alignment.format == "fasta"){
      old.align = Biostrings::readDNAStringSet(paste0(alignment.directory, "/", found.align) )
      found.name = gsub(".fa$", "", found.align)
      found.name = gsub(".fasta$", "", found.name)
    }#end phylip

    #Aligns and then reverses back to correction orientation
    combo.align = runMafft(sequence.data = old.align,
                           add.contigs = align,
                           save.name = paste0(output.directory, "-only/", found.name),
                           algorithm = "add",
                           adjust.direction = TRUE,
                           threads = 1,
                           cleanup.files = T,
                           quiet = quiet,
                           mafft.path = mafft.path)

    #Checks for failed mafft run
    if (length(combo.align) == 0){ return(NULL) }
    #Aligns and then reverses back to correction orientation
    names(combo.align) = gsub(pattern = "^_R_", replacement = "", x = names(combo.align))

    # Duplication changes and such
    if (combine.same.sample == TRUE){

      dup.names = names(combo.align)[duplicated(names(combo.align)) == TRUE]

      if (length(dup.names) != 0){
        nodup.align = combo.align[!names(combo.align) %in% dup.names]

        #Goes through each duplicate
        save.seqs = Biostrings::DNAStringSet()
        for (j in 1:length(dup.names)){

          dup.sample = combo.align[names(combo.align) %in% dup.names[j]]
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

        combo.align = append(save.seqs, nodup.align)
      }#end if to do the thing

    }#end combine.same.sample if

    #Saves them
    write.temp = strsplit(as.character(combo.align), "")
    aligned.set = as.matrix(ape::as.DNAbin(write.temp) )

    #readies for saving
    writePhylip(alignment = aligned.set,
                file=paste0(output.directory, "-only/", save.name, ".phy"),
                interleave = F,
                strict = F)

    print(paste0("Finished ", save.name, " legacy integration successfully!"))
    system(paste0("rm ", save.name, "*"))


  }#end loop

  system(paste0("rm *nucl-blast_db*"))

  ####################################################################################
  if (include.all.together == TRUE){
    ### Copies all files to new folder and then yeah
    for (i in 1:length(align.files)){

      system(paste0("cp ", alignment.directory, "/", align.files[i],
                    " ", output.directory, "-all/", align.files[i]))

    }#end i

    #Delete old
    delete.old = delete.old[duplicated(delete.old) != T]

    if (is.null(delete.old) != TRUE){
      del.string = paste0(output.directory, "-all/", delete.old, collapse = " ")
      system(paste0("rm ", del.string))
    }

    #Add new
    new.files = list.files(paste0(output.directory, "-only"))

    for (i in 1:length(new.files)){
      system(paste0("cp ", output.directory, "-only/", new.files[i], " ",
                    output.directory, "-all/", new.files[i]))
    }#end i loop
  }#end if

  ####################################################################################

}#end fuction

#END SCRIPT
