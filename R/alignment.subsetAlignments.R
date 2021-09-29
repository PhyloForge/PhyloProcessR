#' @title makeAlignmentSubset
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

makeAlignmentSubset = function(alignment.directory = NULL,
                               alignment.format = "phylip",
                               output.directory = NULL,
                               output.format = "phylip",
                               subset.reference = c("fasta", "grep", "blast"),
                               subset.fasta.file = NULL,
                               subset.grep.string = NULL,
                               subset.blast.targets = NULL,
                               blast.path = NULL,
                               threads = 1,
                               memory = 1,
                               overwrite = FALSE) {

  # alignment.directory = "alignments/untrimmed_all-markers"
  # alignment.format = "phylip"
  # output.directory = "alignments/subset_UCEs"
  # output.format = "phylip"
  # subset.reference = "fasta"
  # subset.fasta.file = "/Users/chutter/Dropbox/Research/0_Github/FrogCap_Pipeline/Source_Files/Hutter_uce5k_loci.fa"
  # subset.grep.string = "uce"
  # subset.blast.targets = "/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/Microhylidae_SeqCap/New_Work_2021/Master_Ranoidea_All-Markers_Apr21-2019.fa"

  #Checks the program
  if (is.null(blast.path) == FALSE){
    b.string = unlist(strsplit(blast.path, ""))
    if (b.string[length(b.string)] != "/") {
      blast.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { blast.path = NULL }

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

  #Uses the reference names to isolate
  if (subset.reference == "fasta"){
    fasta.loci = Biostrings::readDNAStringSet(file = subset.fasta.file, format = "fasta")
    ref.names = gsub(".phy$", "", names(fasta.loci) )
    subset.files = align.files[gsub(".phy$", "", align.files) %in% ref.names]

    if (length(subset.files) == 0){ stop("Reference fasta names could not be matched.")}
  }#end if

  #Greps for the file name
  if (subset.reference == "grep"){
    subset.files = align.files[grep(subset.grep.string, align.files)]
    if (length(subset.files) == 0){ stop("Grep name could not be matched.")}
  }#end if

  #uses blast if the names don't match
  if (subset.reference == "blast"){

    #headers for the blast db
    headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
                "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

    #Make blast database for the probe loci
    system(paste0(blast.path, "makeblastdb -in ", subset.fasta.file,
                  " -parse_seqids -dbtype nucl -out subset_nucl-blast_db"), ignore.stdout = quiet)

    #Matches samples to loci
    system(paste0(blast.path, "blastn -task dc-megablast -db subset_nucl-blast_db -evalue 0.001",
                  " -query ", subset.blast.targets, " -out subset-blast-match.txt",
                  " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
                  " -num_threads ", threads))

    #Loads in match data
    match.data = data.table::fread(paste0("subset-blast-match.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
    data.table::setnames(match.data, headers)

    filt.data = match.data[match.data$matches >= ( 0.5 * match.data$tLen),]

    if (nrow(match.data) == 0) {
      stop(paste0("no subset matches to targets. Skipping"))  }

    #Sorting: exon name, contig name, bitscore higher first, evalue
    data.table::setorder(filt.data, qName, tName, -pident, -bitscore, evalue)

    #########################################################################
    #Part B: Multiple sample contigs (tName) matching to one target (qName)
    #########################################################################
    #Pulls out
    contig.names = unique(filt.data[duplicated(filt.data$tName) == T,]$tName)

    save.match = c()
    for (j in 1:length(contig.names)){

      sub.match = filt.data[filt.data$tName %in% contig.names[j],]
      sub.match = sub.match[sub.match$bitscore == max(sub.match$bitscore),][1]
      save.match = rbind(save.match, sub.match)

    }# end j

    good.data = filt.data[!filt.data$tName %in% contig.names,]
    good.data = rbind(good.data, save.match)

    #########################################################################
    #Pulls out
    contig.names = unique(good.data[duplicated(good.data$qName) == T,]$qName)

    save.match = c()
    for (j in 1:length(contig.names)){

      sub.match = good.data[good.data$qName %in% contig.names[j],]
      sub.match = sub.match[sub.match$bitscore == max(sub.match$bitscore),][1]
      save.match = rbind(save.match, sub.match)

    }# end j

    final.data = good.data[!good.data$qName %in% contig.names,]
    final.data = rbind(final.data, save.match)

    #Gets final set of alignments
    subset.files = align.files[gsub(".phy$", "", align.files) %in% final.data$qName]

  }#end blast if

  #save subset files separately
  for (i in 1:length(subset.files)){
    system(paste0("cp ", alignment.directory, "/", subset.files[i]," ",
                  output.directory, "/", subset.files[i]))
  }#end loop

} #end function
