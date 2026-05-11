#' @title genomeStats
#'
#' @description Computes assembly quality statistics for every genome FASTA file
#'   in a directory. Reported metrics include: total assembly size in gigabases,
#'   number of contigs/scaffolds, mean and median contig length, maximum and
#'   minimum contig length, N50, L50, N90, L90, and contig counts above
#'   size thresholds of 1 KB, 10 KB, 100 KB, 500 KB, 1 MB, 10 MB, and 100 MB.
#'   All results are written to a single CSV file.
#'
#' @param genome.directory path to a directory of genome FASTA files, one file
#'   per assembly.
#'
#' @param output.directory path to the directory where the output CSV will be
#'   written. Default: \code{"alignments"}.
#'
#' @param threads number of CPU threads (not currently used). Default: \code{1}.
#'
#' @param memory RAM in GB (not currently used). Default: \code{1}.
#'
#' @param overwrite logical; if \code{TRUE} the output directory is deleted and
#'   recreated. Default: \code{FALSE}.
#'
#' @param resume logical; if \code{TRUE} already-completed samples are skipped.
#'   Cannot be \code{TRUE} when \code{overwrite} is \code{TRUE}. Default:
#'   \code{TRUE}.
#'
#' @param quiet logical; reserved for future use. Default: \code{TRUE}.
#'
#' @return Invisibly returns nothing. Writes
#'   \code{<output.directory>/Spades_assembly_stats.csv} with one row per
#'   genome assembly.
#'
#' @export

genomeStats = function(genome.directory = NULL,
                       output.directory = "alignments",
                       threads = 1,
                       memory = 1,
                       overwrite = FALSE,
                       resume = TRUE,
                       quiet = TRUE) {


  #Read in basic genome info
  setwd("/Volumes/Rodents/Shrew_Genome")
  genome.directory = "/Volumes/Rodents/Shrew_Genome/Draft_assemblies"
  output.directory = "/Volumes/Rodents/Shrew_Genome/genomeStats"
  threads = 4
  memory = 4

  if (is.null(genome.directory) == T){ stop("A directory of genome(s) is needed.") }

  #So I don't accidentally delete everything while testing resume
  if (resume == TRUE & overwrite == TRUE){
    overwrite = FALSE
    stop("Error: resume = T and overwrite = T, cannot resume if you are going to delete everything!")
  }

  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } else { dir.create(output.directory) }

  genome.files = list.files(genome.directory)

  #Stats table prepare
  header.data = c("Assembly", "Total_bp_gb", "N_Contigs", "Mean_bp","Median_bp", "Max_contig_bp", "Min_contig_bp",
                  "N50", "L50","N90", "L90",
                  "N_Contigs_1KB", "N_Contigs_10KB", "N_Contigs_100KB", "N_Contigs_500KB", "N_Contigs_1MB", "N_Contigs_10MB", "N_Contigs_100MB")
  collect.data = data.table::data.table(matrix(as.numeric(0), nrow = length(genome.files), ncol = length(header.data)))
  data.table::setnames(collect.data, header.data)
  collect.data[, Assembly:=as.character(Assembly)]

  for (i in 1:length(genome.files)){

    genome = Biostrings::readDNAStringSet(paste0(genome.directory, "/", genome.files[i]), format = "fasta")

    data.table::set(collect.data, i = as.integer(i), j = match("Assembly", header.data), value = genome.files[i])
    data.table::set(collect.data, i = as.integer(i), j = match("Total_bp_gb", header.data), value = sum(Biostrings::width(genome))/1000000000)
    data.table::set(collect.data, i = as.integer(i), j = match("N_Contigs", header.data), value = length(Biostrings::width(genome)) )
    data.table::set(collect.data, i = as.integer(i), j = match("Mean_bp", header.data), value = mean(Biostrings::width(genome)))
    data.table::set(collect.data, i = as.integer(i), j = match("Median_bp", header.data), value = median(Biostrings::width(genome)))
    data.table::set(collect.data, i = as.integer(i), j = match("Max_contig_bp", header.data), value = max(Biostrings::width(genome)))
    data.table::set(collect.data, i = as.integer(i), j = match("Min_contig_bp", header.data), value = min(Biostrings::width(genome)))

    ### N50: N50 scaffold/contig length is calculated by summing lengths of scaffolds/contigs from the longest to the shortest and determining at what point you reach 50% of the total assembly size.
    ######## The length of the scaffold/contig at that point is the N50 length.
    total.size = sum(Biostrings::width(genome))
    sort.genome = genome[order(Biostrings::width(genome), decreasing = T)]
    index = 0
    current.size = 0
    while (current.size <= total.size/2){
      index = index+1
      current.size = current.size + Biostrings::width(sort.genome[index])
    }

    N50.value = Biostrings::width(sort.genome[index])
    L50.value = index

    ##N90
    index = 0
    current.size = 0
    while (current.size <= total.size*0.9){
      index = index+1
      current.size = current.size + Biostrings::width(sort.genome[index])
    }

    N90.value = Biostrings::width(sort.genome[index])
    L90.value = index

    data.table::set(collect.data, i = as.integer(i), j = match("N50", header.data), value = N50.value)
    data.table::set(collect.data, i = as.integer(i), j = match("L50", header.data), value = L50.value)
    data.table::set(collect.data, i = as.integer(i), j = match("N90", header.data), value = N90.value)
    data.table::set(collect.data, i = as.integer(i), j = match("L90", header.data), value = L90.value)

    data.table::set(collect.data, i = as.integer(i), j = match("N_Contigs_1KB", header.data), value = length(Biostrings::width(genome)[Biostrings::width(genome) > 1000]) )
    data.table::set(collect.data, i = as.integer(i), j = match("N_Contigs_10KB", header.data), value = length(Biostrings::width(genome)[Biostrings::width(genome) > 10000]) )
    data.table::set(collect.data, i = as.integer(i), j = match("N_Contigs_100KB", header.data), value = length(Biostrings::width(genome)[Biostrings::width(genome) > 100000]) )
    data.table::set(collect.data, i = as.integer(i), j = match("N_Contigs_500KB", header.data), value = length(Biostrings::width(genome)[Biostrings::width(genome) > 500000]) )
    data.table::set(collect.data, i = as.integer(i), j = match("N_Contigs_1MB", header.data), value = length(Biostrings::width(genome)[Biostrings::width(genome) > 1000000]) )
    data.table::set(collect.data, i = as.integer(i), j = match("N_Contigs_10MB", header.data), value = length(Biostrings::width(genome)[Biostrings::width(genome) > 10000000]) )
    data.table::set(collect.data, i = as.integer(i), j = match("N_Contigs_100MB", header.data), value = length(Biostrings::width(genome)[Biostrings::width(genome) > 100000000]) )

    #### DEPTH OF COVERAGE





  }# end i loop

write.csv(collect.data, file = paste0(output.directory, "/Spades_assembly_stats.csv"), row.names = F)

}#end funtion



#
# genome.string = paste0(genome.dir, "/",genome.files, collapse = " ")
#
# system(paste0("quast/quast.py ", genome.string,
#               " -r Mus_Genome/Mus_GRCm38.p6_genomic.fna.gz",
#               " -g Mus_Genome/Mus_GRCm38.p6_genomic.gff.gz",
#               " --gene-finding --conserved-genes-finding --eukaryote --large --threads ", threads))
#
# #
#
# ##########################################################################################################
# # Match to targets and get stuff
# ##########################################################################################################
#
#
#
# work.dir = "/Volumes/Rodents/Spades-assembly/K99"
# setwd(work.dir)
#
# #Make blast database for the probe loci
# system(paste("makeblastdb -in final_contigs.fasta -parse_seqids -dbtype nucl ",
#              " -out ", work.dir, "/blast_db", sep = ""))
#
# #Matches samples to loci
# system(paste0("blastn -task dc-megablast -db ", work.dir, "/blast_db",
#              " -query nucleotide-mouse_exons.fa -out genome_match.txt",
#              " -outfmt 6 -num_threads 6"))
#
#
# #headers for the blast db
# headers<-c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
#            "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore")
#
# #Matching
# match.data = fread("genome_match.txt", sep = "\t", header = F, stringsAsFactors = FALSE)
# setnames(match.data, headers)
# match.data = match.data[match.data$matches >= 30,]
# match.data = match.data[match.data$evalue <= 0.01,]
#
# length(unique(match.data$qName))/length(exon.data)
#
# exon.data = scanFa(FaFile("nucleotide-mouse_exons.fa"))
#
#
#
#
