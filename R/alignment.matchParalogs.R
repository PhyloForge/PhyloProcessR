#' @title matchParalogs
#'
#' @description For each sample assembly in assembly.directory, uses BLAST (dc-megablast) to identify contigs that match a provided set of paralog target sequences. Matching contigs are extracted and saved per sample as a fasta file of candidate paralog sequences, enabling downstream paralog detection and filtering.
#'
#' @param assembly.directory path to a directory of per-sample assembly fasta files (one file per sample)
#'
#' @param target.file path to a fasta file of paralog target sequences to match against each assembly
#'
#' @param alignment.contig.name label prefix used for naming the output contig alignment files
#'
#' @param output.directory path to the directory where per-sample paralog output files will be saved
#'
#' @param min.match.percent minimum BLAST percent identity required to retain a contig match
#'
#' @param min.match.length minimum number of matching base pairs required to retain a BLAST hit
#'
#' @param min.match.coverage minimum percent of the query length that must be covered by a BLAST hit
#'
#' @param trim.target if TRUE, trim target sequences before matching (currently unused placeholder)
#'
#' @param threads number of CPU threads to pass to BLAST
#'
#' @param memory total memory in GB (currently unused, reserved for future use)
#'
#' @param overwrite if TRUE, redo and overwrite samples already completed; if FALSE, skip them
#'
#' @param resume if TRUE, skip samples whose output files already exist
#'
#' @param quiet if TRUE, suppress BLAST screen output
#'
#' @param blast.path system path to the directory containing BLAST executables; NULL to use the system PATH
#'
#' @param bbmap.path system path to the directory containing BBMap executables; NULL to use the system PATH
#'
#' @return saves per-sample fasta files of matching paralog contigs to output.directory; nothing is returned to R
#'
#' @export

matchParalogs = function(assembly.directory = NULL,
                         target.file = NULL,
                         alignment.contig.name = NULL,
                         output.directory = "match-targets",
                         min.match.percent = 50,
                         min.match.length = 40,
                         min.match.coverage = 50,
                         trim.target = FALSE,
                         threads = 1,
                         memory = 1,
                         overwrite = FALSE,
                         resume = TRUE,
                         quiet = TRUE,
                         blast.path = NULL,
                         bbmap.path = NULL) {
#

  # comb.loci = list.files("/Volumes/Rodents/Murinae/paralogs", full.names = T)
  #
  # save.loci = DNAStringSet()
  # for (j in 1:length(comb.loci)){
  #
  #   temp.loci = readDNAStringSet(comb.loci[j], format = "fasta")
  #   names(temp.loci) = gsub("_sequence.fa", "", gsub(".*\\/", "", comb.loci[j]))
  #   save.loci = append(save.loci, temp.loci)
  #
  # }
  #
  # final.loci = as.list(as.character(save.loci))
  #
  # PhyloCap::writeFasta(sequences = final.loci, names = names(final.loci),
  #                      "/Volumes/Rodents/Murinae/paralogs/mouse_paralogs.fa", nbchar = 1000000, as.string = T)
  #


#   #Debug setup
  # setwd("/Volumes/Rodents/Murinae/Data_Processing") #Your main project directory
  # assembly.directory<-"/Volumes/Rodents/Murinae/Data_Processing/iupacAssemblies"
  # target.file<-"/Volumes/Rodents/Murinae/paralogs/mouse_paralogs.fa"
  # output.directory = "paralogFinder"
  # alignment.contig.name = "paralog_murinae"
  # #
  # # #Main settings
  #  threads = 4
  #  memory = 8
  #  trim.target = FALSE
  #  overwrite = FALSE
  #  resume = TRUE
  #  quiet = TRUE
  # #
  # # #tweak settings (make some statements to check these)
  #  min.match.percent = 60
  #  min.match.length = 50
  #  min.match.coverage = 50
  # #
  # # #program paths
  #  blast.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin"
  #  bbmap.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin"

  #Add the slash character to path
  if (is.null(blast.path) == FALSE){
    b.string = unlist(strsplit(blast.path, ""))
    if (b.string[length(b.string)] != "/") {
      blast.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { blast.path = "" }

  #Same adds to bbmap path
  if (is.null(bbmap.path) == FALSE){
    b.string = unlist(strsplit(bbmap.path, ""))
    if (b.string[length(b.string)] != "/") {
      bbmap.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { bbmap.path = "" }

  #Initial checks
  if (assembly.directory == output.directory){ stop("You should not overwrite the original contigs.") }
  if (is.null(target.file) == T){ stop("A fasta file of targets to match to assembly contigs is needed.") }

  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } else { dir.create(output.directory) }

  #Gets contig file names
  file.names = list.files(assembly.directory)

  #headers for the blast db
  headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
             "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

  #Matching and processing for each sample
  for (i in 1:length(file.names)) {

    #Sets up working directories for each species
    sample = gsub(pattern = ".fa$", replacement = "", x = file.names[i])
    species.dir = paste0(output.directory, "/", sample)

    #Creates species directory if none exists
    if (file.exists(species.dir) == F){ dir.create(species.dir) }

    #Checks if this has been done already
    if (overwrite == FALSE){
      if (file.exists(paste0(species.dir, "/", sample, "_matching-contigs.fa")) == T){
        print(paste0(sample, " already finished, skipping. Set overwrite to T if you want to overwrite."))
        next
      }
    }#end

    #########################################################################
    #Part A: Blasting
    #########################################################################

    #Reads in contigs
    contigs = Biostrings::readDNAStringSet(paste0(assembly.directory, "/", file.names[i]), format = "fasta")
    names(contigs) = paste0("contig_", stringr::str_pad(seq(1:length(contigs)), 6, pad = "0"))

    #Finds probes that match to two or more contigs
    final.loci = as.list(as.character(contigs))
    PhyloCap::writeFasta(sequences = final.loci, names = names(final.loci),
               paste0(species.dir, "/", sample, "_renamed-contigs.fa"), nbchar = 1000000, as.string = T)

    #Make blast database for the probe loci
    system(paste0(blast.path, "makeblastdb -in ", species.dir, "/", sample, "_renamed-contigs.fa",
                  " -parse_seqids -dbtype nucl -out ", species.dir, "/", sample, "_nucl-blast_db"), ignore.stdout = quiet)

    #Matches samples to loci
    system(paste0(blast.path, "blastn -task dc-megablast -db ", species.dir, "/", sample, "_nucl-blast_db -evalue 0.001",
                  " -query ", target.file, " -out ", species.dir, "/", sample, "_target-blast-match.txt",
                  " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
                  " -num_threads ", threads))

    #Need to load in transcriptome for each species and take the matching transcripts to the database
    system(paste0("rm ", species.dir, "/*nucl-blast_db*"))

    #Loads in match data
    match.data = data.table::fread(paste0(species.dir, "/", sample, "_target-blast-match.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
    data.table::setnames(match.data, headers)

    #Matches need to be greater than 12
    filt.data = match.data[match.data$matches > min.match.length,]
    #Percent identitiy must match 50% or greater
    filt.data = filt.data[filt.data$pident >= min.match.percent,]

    if (nrow(filt.data) == 0) {
      print(paste0(sample, " had no matches. Skipping"))
      next }

    #Sorting: exon name, contig name, bitscore higher first, evalue
    data.table::setorder(filt.data, qName, tName, -pident, -bitscore, evalue)

    #Make sure the hit is greater than 50% of the reference length
    #filt.data = filt.data[filt.data$matches >= ( (min.match.coverage/100) * filt.data$qLen),]

    #Reads in contigs
    contigs = Biostrings::readDNAStringSet(paste0(species.dir, "/", sample, "_renamed-contigs.fa"), format = "fasta")
    names(contigs) = gsub(" .*", "", names(contigs))

    paralog.loci = unique(filt.data$qName)

    all.paralogs = Biostrings::DNAStringSet()
    for (j in 1:length(paralog.loci)){

        temp.para = filt.data[filt.data$qName %in% paralog.loci[j],]
        para.contigs = contigs[names(contigs) %in% temp.para$tName]
        names(para.contigs) = paste0(sample, "_|_", paralog.loci[j], "-p", seq(1:length(para.contigs)))
        all.paralogs = append(all.paralogs, para.contigs)

    }

    #Finds probes that match to two or more contigs
    final.loci = as.list(as.character(all.paralogs))
    PhyloCap::writeFasta(sequences = final.loci, names = names(final.loci),
               paste0(species.dir, "/", sample, "_matching-paralogs.fa"), nbchar = 1000000, as.string = T)

    print(paste0(sample, " paralog matching complete. ", length(final.loci), " targets found!"))

    system(paste0("rm ", species.dir, "/", sample, "_target-blast-match.txt"))
    system(paste0("rm ", species.dir, "/", sample, "_renamed-contigs.fa"))

  }# end i loop
}



#### END SCRIPT
