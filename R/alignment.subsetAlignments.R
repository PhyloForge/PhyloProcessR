#' @title makeAlignmentSubset
#'
#' @description Selects a subset of alignments from alignment.directory and copies them to output.directory. The subset can be defined in three ways: by matching alignment file names to names in a reference fasta file ("fasta"), by grepping for a pattern in alignment file names ("grep"), or by using BLAST to match alignment names against a set of target sequences ("blast"). This allows flexible subsetting of large alignment directories by locus type or identity.
#'
#' @param alignment.directory path to the directory of input alignment files
#'
#' @param alignment.format format of the input alignments; "phylip" or "fasta"
#'
#' @param output.directory path to the directory where the subset alignment files will be copied
#'
#' @param output.format format for output alignments (currently informational; files are copied directly)
#'
#' @param subset.reference method used to define the subset: "fasta" matches alignment names to sequence names in a fasta file, "grep" matches alignment file names to a pattern string, "blast" uses BLAST similarity to a set of target sequences
#'
#' @param subset.fasta.file path to a fasta file whose sequence names are used to select alignments (required when subset.reference is "fasta")
#'
#' @param subset.grep.string a regular expression pattern passed to grep to select alignment files by name (required when subset.reference is "grep")
#'
#' @param subset.blast.targets path to a fasta file of query sequences used as BLAST queries against a database built from subset.fasta.file (required when subset.reference is "blast")
#'
#' @param blast.path system path to the directory containing BLAST executables; NULL to use the system PATH
#'
#' @param threads number of CPU threads to pass to BLAST
#'
#' @param memory total memory in GB (currently reserved for future use)
#'
#' @param overwrite if TRUE, overwrite an existing output directory; if FALSE, keep existing files
#'
#' @return copies the selected alignment files to output.directory; nothing is returned to R
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
