#' @title plotOccupancy
#'
#' @description Creates a bar plot showing the number (or proportion) of
#'   alignment files in which each sample is present. Samples are read from the
#'   tip labels of a species tree and their presence/absence is scored across all
#'   phylip alignment files in \code{file.directory}. The bar plot is saved as a
#'   PDF and the ggplot2 object is returned invisibly.
#'
#' @param species.tree path to a Newick tree file whose tip labels are the
#'   expected sample names.
#'
#' @param file.directory path to a directory of phylip alignment files (.phy).
#'   Each file is scanned for the presence of each sample.
#'
#' @param type character; type of input files. Accepted values are
#'   \code{"alignment"} and \code{"tree"} (currently only alignment is
#'   implemented). Default: \code{c("alignment", "tree")}.
#'
#' @param out.name base name (without extension) used for the output PDF file.
#'   Default: \code{"species-occupancy"}.
#'
#' @param sample.order character; how to order samples on the x-axis. One of
#'   \code{"value"} (ascending occupancy), \code{"alphabetical"}, or
#'   \code{"custom"} (uses \code{custom.order}). Default:
#'   \code{c("value", "alphabetical", "custom")}.
#'
#' @param save.width width of the saved PDF in inches. Default: \code{10}.
#'
#' @param save.height height of the saved PDF in inches. Default: \code{8}.
#'
#' @param custom.order character vector giving the desired sample order when
#'   \code{sample.order = "custom"}. Default: \code{NULL}.
#'
#' @param exclude.taxa character vector of sample names to exclude from the
#'   plot. Default: \code{NULL}.
#'
#' @param proportion logical; if \code{TRUE} values are expressed as a
#'   proportion of the maximum occupancy rather than raw counts. Default:
#'   \code{TRUE}.
#'
#' @param overwrite logical; reserved for future use. Default: \code{TRUE}.
#'
#' @return Returns the ggplot2 plot object invisibly. Also saves a PDF named
#'   \code{<out.name>_occupancy-plot.pdf} in the working directory.
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
