#' @title paralogStats
#'
#' @description Function for batch trimming a folder of alignments, with the various trimming functions available to select from
#'
#' @param assembly.directory path to a folder of sequence alignments in phylip format.
#'
#' @param target.file available input alignment formats: fasta or phylip
#'
#' @param alignment.contig.name contigs are added into existing alignment if algorithm is "add"
#'
#' @param output.directory available output formats: phylip
#'
#' @param min.match.percent algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param min.match.length TRUE applies the adjust sequence direction function of MAFFT
#'
#' @param min.match.coverage if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param threads path to a folder of sequence alignments in phylip format.
#'
#' @param memory give a save name if you wnat to save the summary to file.
#'
#' @param trim.target TRUE to supress mafft screen output
#'
#' @param overwrite path to a folder of sequence alignments in phylip format.
#'
#' @param resume contigs are added into existing alignment if algorithm is "add"
#'
#' @param quiet algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param blast.path algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param bbmap.path algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
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

geneFamilyStats = function(assembly.directory = NULL,
                           target.file = NULL,
                           gene.family.file = NULL,
                           output.directory = "gene-family-stats",
                           min.match.percent = 50,
                           min.match.length = 60,
                           min.match.coverage = 50,
                           max.match.copy = 100,
                           threads = 1,
                           memory = 1,
                           overwrite = FALSE,
                           quiet = TRUE,
                           blast.path = NULL) {


  #Debug setup
  # assembly.directory = "/Users/chutter/Dropbox/SharewithCarl/renamed_contigs"
  # output.directory = "/Users/chutter/Dropbox/SharewithCarl/gene-family-stats"
  # target.file = "/Users/chutter/Dropbox/SharewithCarl/venom_loci_updated_April17-2024.fa"
  # gene.family.file = "/Users/chutter/Dropbox/SharewithCarl/family_file.csv"
  # blast.path = "/Users/chutter/Bioinformatics/anaconda3/envs/mitocap/bin/"
  #
  # quiet = TRUE
  # overwrite = FALSE
  # threads = 6
  # memory = 6
  #
  # # #tweak settings (make some statements to check these)
  #  min.match.percent = 50
  #  min.match.length = 50
  #  min.match.coverage = 50
  #  max.match.copy = 100

  #Add the slash character to path
  if (is.null(blast.path) == FALSE){
    b.string = unlist(strsplit(blast.path, ""))
    if (b.string[length(b.string)] != "/") {
      blast.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { blast.path = "" }

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

  #Read in gene family file
  gene.data = read.csv(gene.family.file)
  gene.family.names = unique(gene.data$gene_family)

  #headers for the blast db
  headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
             "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

  #Stats table prepare
  header.data = c("sample", "total_nucleotides", "total_megabases", "total_contigs",
                  "unique_paralogs", "all_paralogs", "mean_copy", "median_copy", "min_copy", "max_copy",
                  "mean_contig_length","median_contig_length", "max_contig_length", "min_contig_length")
  collect.data = data.table::data.table(matrix(as.numeric(0), nrow = length(file.names), ncol = length(header.data)))
  data.table::setnames(collect.data, header.data)
  collect.data[, sample:=as.character(sample)]

  #Stats table prepare
  header.family = c("sample", gene.family.names)
  collect.data.family = data.table::data.table(matrix(as.numeric(0), nrow = length(file.names), ncol = length(header.family)))
  data.table::setnames(collect.data.family, header.family)
  collect.data.family[, sample:=as.character(sample)]

  combined.gene.family = c()

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
    filt.data = filt.data[(filt.data$matches / filt.data$qLen) * 100 >= min.match.coverage,]

    if (nrow(filt.data) == 0) {
      print(paste0(sample, " had no matches. Skipping"))
      next }

    #Sorting: exon name, contig name, bitscore higher first, evalue
    data.table::setorder(filt.data, qName, tName, -pident, -bitscore, evalue)

    #Make sure the hit is greater than 50% of the reference length
    filt.data = filt.data[filt.data$matches >= ( (min.match.coverage/100) * filt.data$qLen),]

    write.table(filt.data, file = paste0(species.dir, "/", sample, "_gene-family_filtered.txt"), sep = "\t", row.names = F, quote = F)

    ######################################################
    #### Filtering out too many paralogs
    ######################################################

    paralog.loci = unique(filt.data$qName)

    all.paralogs = Biostrings::DNAStringSet()
    para.count = c()
    remove.data = c()
    for (j in 1:length(paralog.loci)){

      temp.para = filt.data[filt.data$qName %in% paralog.loci[j],]

      if (nrow(temp.para) > max.match.copy){
        remove.data = append(remove.data, paralog.loci[j])
        next
      }

      para.contigs = contigs[names(contigs) %in% temp.para$tName]
      if (length(para.contigs) <= 1){ next}
      names(para.contigs) = paste0(sample, "_|_", paralog.loci[j], "-p", seq(1:length(para.contigs)))
      all.paralogs = append(all.paralogs, para.contigs)
      para.count = append(para.count, length(para.contigs))
    }

    if (is.null(para.count) == TRUE){ para.count = 0 }

    filt.data = filt.data[!filt.data$qName %in% remove.data,]

    ######################################################
    #### Individual sample section
    ######################################################

    #Stats table prepare
    header.gene = c("sample", "gene_family", "number_copies", "total_targets", "percent_total_targets", "total_bases",
                    "mean_copy_length","median_copy_length",
                    "max_copy_length", "min_copy_length")
    collect.data.gene = data.table::data.table(matrix(as.numeric(0), nrow = length(gene.family.names), ncol = length(header.gene)))
    data.table::setnames(collect.data.gene, header.gene)
    collect.data.gene[, gene_family:=as.character(gene.family.names)]
    collect.data.gene[, sample:=as.character(sample)]

    for (j in 1:length(gene.family.names)){

      temp.gene = gene.data[gene.data$gene_family %in% gene.family.names[j],]

      temp.data = filt.data[filt.data$qName %in% temp.gene$target_name,]

      data.table::set(collect.data.gene, i = match(gene.family.names[j], gene.family.names), j = match("sample", header.gene), value =  gsub(".fa$|.fasta$", "", file.names[i]) )

      if (nrow(temp.data) == 0) { next }

      data.table::set(collect.data.gene, i = match(gene.family.names[j], gene.family.names), j = match("number_copies", header.gene), value = nrow(temp.data) )
      data.table::set(collect.data.gene, i = match(gene.family.names[j], gene.family.names), j = match("total_bases", header.gene), value = sum(temp.data$matches) )
      data.table::set(collect.data.gene, i = match(gene.family.names[j], gene.family.names), j = match("total_targets", header.gene), value = length(unique(temp.data$qName)) )
      data.table::set(collect.data.gene, i = match(gene.family.names[j], gene.family.names), j = match("percent_total_targets", header.gene), value = length(unique(temp.data$qName))/length(unique(temp.gene$target_name)) )
      data.table::set(collect.data.gene, i = match(gene.family.names[j], gene.family.names), j = match("mean_copy_length", header.gene), value = mean(temp.data$matches) )
      data.table::set(collect.data.gene, i = match(gene.family.names[j], gene.family.names), j = match("median_copy_length", header.gene), value = median(temp.data$matches) )
      data.table::set(collect.data.gene, i = match(gene.family.names[j], gene.family.names), j = match("max_copy_length", header.gene), value = max(temp.data$matches) )
      data.table::set(collect.data.gene, i = match(gene.family.names[j], gene.family.names), j = match("min_copy_length", header.gene), value = min(temp.data$matches) )

      data.table::set(collect.data.family, i = as.integer(i), j = match("sample", header.family), value = gsub(".fa$|.fasta$", "", file.names[i]) )
      data.table::set(collect.data.family, i = as.integer(i), j = match(gene.family.names[j], header.family), value = nrow(temp.data) )

    } #end j loop

    write.table(collect.data.gene, file = paste0(species.dir, "/", sample, "_gene-family_summary.txt"),
                sep = "\t", row.names = F, quote = F)

    combined.gene.family = rbind(combined.gene.family, collect.data.gene)

    ######################################################
    #### Sample Summary section
    ######################################################

    #Reads in contigs
    contigs = Biostrings::readDNAStringSet(paste0(species.dir, "/", sample, "_renamed-contigs.fa"), format = "fasta")
    names(contigs) = gsub(" .*", "", names(contigs))

    data.table::set(collect.data, i = as.integer(i), j = match("sample", header.data), value = gsub(".fa$|.fasta$", "", file.names[i]))
    data.table::set(collect.data, i = as.integer(i), j = match("total_nucleotides", header.data), value = sum(Biostrings::width(all.paralogs)))
    data.table::set(collect.data, i = as.integer(i), j = match("total_megabases", header.data), value = sum(Biostrings::width(all.paralogs))/1000000)
    data.table::set(collect.data, i = as.integer(i), j = match("total_contigs", header.data), value = length(Biostrings::width(all.paralogs)) )

    data.table::set(collect.data, i = as.integer(i), j = match("unique_paralogs", header.data), value = length(paralog.loci) )
    data.table::set(collect.data, i = as.integer(i), j = match("all_paralogs", header.data), value = length(para.count) )
    data.table::set(collect.data, i = as.integer(i), j = match("mean_copy", header.data), value = mean(para.count) )
    data.table::set(collect.data, i = as.integer(i), j = match("median_copy", header.data), value = median(para.count) )
    data.table::set(collect.data, i = as.integer(i), j = match("min_copy", header.data), value = min(para.count) )
    data.table::set(collect.data, i = as.integer(i), j = match("max_copy", header.data), value = max(para.count) )

    data.table::set(collect.data, i = as.integer(i), j = match("mean_contig_length", header.data), value = mean(Biostrings::width(all.paralogs)))
    data.table::set(collect.data, i = as.integer(i), j = match("median_contig_length", header.data), value = median(Biostrings::width(all.paralogs)))
    data.table::set(collect.data, i = as.integer(i), j = match("max_contig_length", header.data), value = max(Biostrings::width(all.paralogs)))
    data.table::set(collect.data, i = as.integer(i), j = match("min_contig_length", header.data), value = min(Biostrings::width(all.paralogs)))

    print(paste0(sample, " gene family stats completed. ", length(all.paralogs), " paralogs found!"))

    system(paste0("rm ", species.dir, "/", sample, "_renamed-contigs.fa"))

  }# end i loop

  write.table(collect.data, file = paste0(output.directory, "/sample_summary.txt"), sep = "\t", row.names = F, quote = F)
  write.table(collect.data.family, file = paste0(output.directory, "/gene-family_summary_short.txt"), sep = "\t", row.names = F, quote = F)
  write.table(combined.gene.family, file = paste0(output.directory, "/gene-family_summary_long.txt"), sep = "\t", row.names = F, quote = F)

}



#### END SCRIPT
