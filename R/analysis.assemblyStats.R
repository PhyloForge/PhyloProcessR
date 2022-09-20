#' @title assemblyStats
#'
#' @description Function for batch trimming a folder of alignments, with the various trimming functions available to select from
#'
#' @param genome.directory path to a folder of sequence alignments in phylip format.
#'
#' @param output.directory available input alignment formats: fasta or phylip
#'
#' @param threads contigs are added into existing alignment if algorithm is "add"
#'
#' @param threads path to a folder of sequence alignments in phylip format.
#'
#' @param memory give a save name if you wnat to save the summary to file.
#'
#' @param overwrite TRUE to supress mafft screen output
#'
#' @param resume TRUE to supress mafft screen output
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

analysis.assemblyStats = function(assembly.directory = NULL,
                                  output.directory = "assembly-stats",
                                  overwrite = FALSE) {


  # assembly.directory = "/Volumes/LaCie/VenomCap/data-analysis/draft-assemblies"
  # output.directory = "/Volumes/LaCie/VenomCap/data-analysis/assembly-stats"
  #
  # samtools.path = "/Users/chutter/Bioinformatics/anaconda3/envs/mitocap/bin/"
  # bwa.path = "/Users/chutter/Bioinformatics/anaconda3/envs/mitocap/bin/"
  # picard.path = "/Users/chutter/Bioinformatics/anaconda3/envs/mitocap/bin/"
  #
  # quiet = TRUE
  # overwrite = FALSE
  # threads = 6
  # memory = 6

  if (is.null(assembly.directory) == T){ stop("A directory of genome(s) is needed.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (file.exists(output.directory) == F){ dir.create(output.directory) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  }#end else

  assembly.files = list.files(assembly.directory)
  sample.names = gsub(".fa|.fasta", "", assembly.files)

  #Stats table prepare
  header.data = c("sample", "total_nucleotides", "total_megabases", "total_contigs",
                  "mean_contig_length","median_contig_length", "max_contig_length", "min_contig_length")
  collect.data = data.table::data.table(matrix(as.numeric(0), nrow = length(assembly.files), ncol = length(header.data)))
  data.table::setnames(collect.data, header.data)
  collect.data[, sample:=as.character(sample)]

  for (i in 1:length(assembly.files)){

    assembly.data = Biostrings::readDNAStringSet(paste0(assembly.directory, "/", assembly.files[i]), format = "fasta")

    data.table::set(collect.data, i = as.integer(i), j = match("sample", header.data), value = gsub(".fa|.fasta", "", assembly.files[i]))
    data.table::set(collect.data, i = as.integer(i), j = match("total_nucleotides", header.data), value = sum(Biostrings::width(assembly.data)))
    data.table::set(collect.data, i = as.integer(i), j = match("total_megabases", header.data), value = sum(Biostrings::width(assembly.data))/1000000)
    data.table::set(collect.data, i = as.integer(i), j = match("total_contigs", header.data), value = length(Biostrings::width(assembly.data)) )
    data.table::set(collect.data, i = as.integer(i), j = match("mean_contig_length", header.data), value = mean(Biostrings::width(assembly.data)))
    data.table::set(collect.data, i = as.integer(i), j = match("median_contig_length", header.data), value = median(Biostrings::width(assembly.data)))
    data.table::set(collect.data, i = as.integer(i), j = match("max_contig_length", header.data), value = max(Biostrings::width(assembly.data)))
    data.table::set(collect.data, i = as.integer(i), j = match("min_contig_length", header.data), value = min(Biostrings::width(assembly.data)))

  }# end i loop

  write.table(collect.data, file = paste0(output.directory, "/sample-assembly-summary.txt"), sep = "\t", row.names = F)

}#end funtion


#END SCRIPT

