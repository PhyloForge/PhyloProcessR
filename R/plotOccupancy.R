#' @title plotOccupancy
#'
#' @description Function for batch trimming a folder of alignments, with the various trimming functions available to select from
#'
#' @param species.tree path to a folder of sequence alignments in phylip format.
#'
#' @param file.directory available input alignment formats: fasta or phylip
#'
#' @param type contigs are added into existing alignment if algorithm is "add"
#'
#' @param out.name available output formats: phylip
#'
#' @param sample.order algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param custom.order TRUE applies the adjust sequence direction function of MAFFT
#'
#' @param exclude.taxa if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param proportion path to a folder of sequence alignments in phylip format.
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

plotOccupancy = function(species.tree = NULL,
                         file.directory = NULL,
                         type = c("alignment", "tree"),
                         out.name = "species-occupancy",
                         sample.order = c("value", "alphabetical", "custom"),
                         save.width = 10,
                         save.height = 8,
                         custom.order = NULL,
                         exclude.taxa = NULL,
                         proportion = TRUE,
                         overwrite = TRUE) {


  #Debug
  # sample.order = "custom"
  # custom.order = c("Boophis_tephraeomystax_CRH1675", "Tsingymantis_antitra_FGZC1128",
  #                  "Aglyptodactylus_securifer_CRH1644", "Laliostoma_labrosum_ZCMV5617",
  #                  "Blommersia_grandisonae_CRH792", "Boehmantis_microtympanum_FGZC132",
  #                  "Gephyromantis_redimitus_CRH1628", "Guibemantis_depressiceps_CRH535",
  #                  "Mantella_baroni_CRH1027", "Mantidactylus_melanopleura_CRH1998",
  #                  "Spinomantis_aglavei_JJW2354", "Wakea_madinika_2001F54")
  #
  # proportion = TRUE
  # exclude.taxa = c("Kalophrynus_pleurostigma_CDS6115", "Amolops_ricketti_KUFS65", "Platymantis_corrugatus_RMB15045", "Rhacophorus_bipunctatus_221351")
  # out.name = "sequence-capture_occupancy"
  # work.dir = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Mantellidae_Subfamily/Statistics"
  # file.directory = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Mantellidae_Subfamily/Trees/Sequence_capture/Alignments/all-markers_untrimmed"
  # species.tree = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Mantellidae_Subfamily/Statistics/species_tree.tre"
  # setwd(work.dir)

  #Load in loci and sample names
  spp.tree = ape::read.tree(species.tree)
  sample.names = spp.tree$tip.label

  #Get locus names
  locus.names = unique(list.files(file.directory))
  save.names = gsub(".phy$", "", locus.names)

  #Gathers number of samples from each alignment for each sample
  header.data = c("Sample", save.names)
  sc.data = data.table::data.table(matrix(as.integer(0), nrow = length(sample.names), ncol = length(header.data)))
  data.table::setnames(sc.data, header.data)
  sc.data[, Sample:=as.character(sample.names)]

  #Loop through each locus and make a big matrix
  for (j in 1:length(locus.names)){
    #Reads in files
    align = Biostrings::readAAMultipleAlignment(file = paste0(file.directory, "/", locus.names[j]), format = "phylip")

    #Use taxa remove
    tax.names = names(Biostrings::DNAStringSet(align))
    data.table::set(sc.data, i =  match(tax.names, sample.names), j = as.integer(j+1), value = as.integer(1) )
  }#end j loop

  #Sets up final table
  if (proportion == TRUE){
    sample.total = data.table::data.table(Sample = sc.data$Sample,
                                          Value = rowSums(sc.data[,-1])/max(rowSums(sc.data[,-1])) )
  } else {
    sample.total = data.table::data.table(Sample = sc.data$Sample,
                                          Value = rowSums(sc.data[,-1]))
  }

  sample.total = sample.total[!sample.total$Sample %in% exclude.taxa]

  if (sample.order == "value"){
    sample.total = sample.total[order(sample.total$Value)]
  }
  if (sample.order == "alphabetical"){
    sample.total = sample.total[order(sample.total$Sample)]
  }
  if (sample.order == "custom"){
    sample.total = sample.total[pmatch(custom.order, sample.total$Sample),]
  }

  sample.total$Sample = factor(sample.total$Sample, levels=unique(sample.total$Sample))

  #Sets up PDF to save plot
  pdf(paste0(out.name, "_occupancy-plot.pdf"), width=save.width, height=save.height, compress=F)

  #Does the GGplot stuff
  p1 = ggplot(sample.total, aes(x = Sample, y = Value, fill = Sample)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size=12,angle = 45, hjust = 1, colour = "black"),
          axis.text.y = element_text(size=12,angle = 0, hjust = 1, colour = "black"))

  print(p1)

  dev.off()

  return(p1)

}#end function
