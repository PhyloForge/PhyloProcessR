#' @title hyphyOmega
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

hyphyParser = function(results.directory = NULL,
                      hyphy.analysis = c("omega", "BUSTED"),
                      output.name = "BUSTED",
                      threads = 1,
                      memory = 1,
                      overwrite = FALSE,
                      quiet = TRUE,
                      hyphy.path = NULL) {


  #Read in basic genome info
  # library(PhyloCap)
  # setwd("/Volumes/Rodents/Australian_Rodents/Data_Processing")
  # tree.directory= "/Volumes/Rodents/Australian_Rodents/Data_Processing/Trees/Ausfull/genes_trimmed_trees"
  # alignment.directory = "/Volumes/Rodents/Australian_Rodents/Data_Processing/Alignments/Ausfull/coding_trimmed/nt"
  # metadata.file = "/Volumes/Rodents/Australian_Rodents/Data_Processing/Mus-selected-sequences_metadata_final.csv"
  # dataset.name = "Omega"
  # threads = 4
  # memory = 4
  # resume = T
  # overwrite = F
  # hyphy.path = "/usr/local/bin"
  # mg94.path = "/Users/chutter/hyphy-analyses/FitMG94"

  #Directoires
  setwd("/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Rodent_Mitochondrial")
  results.directory = "/Volumes/Rodents/Australian_Rodents/Data_Processing/hyphy/busted_all"
  output.name = "busted_all"
  hyphy.analysis = "BUSTED"
  threads = 4
  memory = 8
  resume = T
  overwrite = F
  quiet = T
  hyphy.path = "/usr/local/bin/"


  #Same adds to bbmap path
  if (is.null(hyphy.path) == FALSE){
    b.string = unlist(strsplit(hyphy.path, ""))
    if (b.string[length(b.string)] != "/") {
      hyphy.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { hyphy.path = "" }

  if (is.null(results.directory) == T){ stop("A directory of results is needed.") }

  results.files = list.dirs(results.directory, full.names = F)
  results.files = results.files[results.files != ""]

  #################################################################################
  ################################# BUSTED ############################
  ######################################################
  if (hyphy.analysis == "BUSTED"){

    #Stats table prepare
    #header.data = c("sample", "constrained", "unconstrained", "gtr_model", "mg94xrev")
    #collect.data = data.table::data.table(matrix(as.numeric(0), nrow = length(results.files), ncol = length(header.data)))
    #data.table::setnames(collect.data, header.data)
    #collect.data[, sample:=as.character(sample)]

    #Run analyses through all results output folders
    all.data = data.frame()
    for (i in 1:length(results.files)){

      json.data = jsonlite::fromJSON(paste0(results.directory, "/", results.files[i], "/BUSTED-results.json"))
      branch.att = json.data$`branch attributes`
      branch.att = unlist(branch.att[[1]])

      #Gather stats
      con.stats = branch.att[grep("\\.constrained", names(branch.att) )]
      uncon.stats = branch.att[grep("\\.unconstrained", names(branch.att) )]
      nucl.stats = branch.att[grep("\\.Nucleotide", names(branch.att) )]
      omega.stats = branch.att[grep("\\.MG94xREV", names(branch.att) )]

      #Make table from stats
      con.data = data.frame(Sample = gsub("\\.constrained", "", names(con.stats)),
                 null_constrained = as.numeric(con.stats))

      uncon.data = data.frame(Sample = gsub("\\.unconstrained", "", names(uncon.stats)),
                 unconstrained = as.numeric(uncon.stats))

      nucl.data = data.frame(Sample = gsub("\\.Nucleotide GTR", "", names(nucl.stats)),
                              gtr_model = as.numeric(nucl.stats))

      omega.data = data.frame(Sample = gsub("\\.MG94xREV.*", "", names(omega.stats)),
                              mg94xrev = as.numeric(omega.stats))

      if (nrow(con.data) != 0) { a.data = merge(con.data, uncon.data, by = "Sample") } else {
        a.data = uncon.data
        a.data$constrained = NA
      }
      b.data = merge(nucl.data, omega.data, by = "Sample")
      w.data = merge(a.data, b.data, by = "Sample")

      save.data = cbind(Locus = results.files[i], w.data)
      all.data = rbind(all.data, save.data)
    } #end i loop

  write.csv(all.data, file = paste0(output.name, "_BUSTED-stats.csv"), row.names = F, quote = F)
  print(paste0("Finished ", hyphy.analysis, " data summary!"))

  }#end busted if
  #################################################################################



}#end function


