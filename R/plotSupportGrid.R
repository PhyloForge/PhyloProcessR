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

plotSupportGrid = function(species.tree = NULL,
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


  #Directory of trees
  tree.dir = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Mantellidae_Subfamily/Tree_Grid"

  #Outgroups
  outgroups = c("Kalophrynus_pleurostigma_CDS6115", "Amolops_ricketti_KUFS65",
                "Platymantis_corrugatus_RMB15045", "Rhacophorus_bipunctatus_221351")

  # Mantellidae
  taxa.set = list()
  taxa.set[[1]] = c("Boophis_tephraeomystax_CRH1675", "Tsingymantis_antitra_FGZC1128",
                   "Aglyptodactylus_securifer_CRH1644", "Laliostoma_labrosum_ZCMV5617",
                   "Blommersia_grandisonae_CRH792", "Boehmantis_microtympanum_FGZC132",
                   "Gephyromantis_redimitus_CRH1628", "Guibemantis_depressiceps_CRH535",
                   "Mantella_baroni_CRH1027", "Mantidactylus_melanopleura_CRH1998",
                   "Spinomantis_aglavei_JJW2354", "Wakea_madinika_2001F54")

  # Laliostominane
  taxa.set[[2]] = c("Aglyptodactylus_securifer_CRH1644", "Laliostoma_labrosum_ZCMV5617")

  #Tsingymantinae
  taxa.set[[3]] = c("Boophis_tephraeomystax_CRH1675", "Tsingymantis_antitra_FGZC1128",
                    "Blommersia_grandisonae_CRH792", "Boehmantis_microtympanum_FGZC132",
                    "Gephyromantis_redimitus_CRH1628", "Guibemantis_depressiceps_CRH535",
                    "Mantella_baroni_CRH1027", "Mantidactylus_melanopleura_CRH1998",
                    "Spinomantis_aglavei_JJW2354", "Wakea_madinika_2001F54")

  #Boophinae
  taxa.set[[4]] = c("Boophis_tephraeomystax_CRH1675",
                    "Blommersia_grandisonae_CRH792", "Boehmantis_microtympanum_FGZC132",
                    "Gephyromantis_redimitus_CRH1628", "Guibemantis_depressiceps_CRH535",
                    "Mantella_baroni_CRH1027", "Mantidactylus_melanopleura_CRH1998",
                    "Spinomantis_aglavei_JJW2354", "Wakea_madinika_2001F54")

  #Mantellinae
  taxa.set[[5]] = c("Blommersia_grandisonae_CRH792", "Boehmantis_microtympanum_FGZC132",
                    "Gephyromantis_redimitus_CRH1628", "Guibemantis_depressiceps_CRH535",
                    "Mantella_baroni_CRH1027", "Mantidactylus_melanopleura_CRH1998",
                    "Spinomantis_aglavei_JJW2354", "Wakea_madinika_2001F54")

  names(taxa.set) = c("Mantellidae", "Laliostominae", "Tsingymantis", "Boophinae", "Mantellinae")

  #####
  #### Start function here

  #Gets the new tree fiels it made
  header.data = c("File_name", "Dataset", "Analysis", names(taxa.set))

  #Gets the new tree fiels it made
  tree.files = list.files(tree.dir, full.names = F, recursive = T)
  tree.files = tree.files[tree.files != ""]

  ### Transcriptome database
  tree.files = tree.files[grep("SeqCap", tree.files)]

  collect.data = data.table(matrix(as.numeric(0), nrow = length(tree.files), ncol = length(header.data)))
  setnames(collect.data, header.data)
  collect.data[, File_name:=as.character(File_name)]
  collect.data[, Analysis:=as.character(Analysis)]
  collect.data[, Dataset:=as.character(Dataset)]

  supp.grid = data.frame()
  for (i in 1:length(tree.files)){

    temp.tree = c()
    if (length(grep("Jackknife", tree.files[i])) == 1){
      temp.tree = read.tree(paste0(tree.dir, "/", tree.files[i]))
      dataset.name = gsub("_jackknife.tre$", "", tree.files[i])
      dataset.name = gsub(".*/", "", dataset.name)
      analysis.name = "Jackknife"

    }

    if (length(grep("SVD", tree.files[i])) == 1){
      temp.tree = read.nexus(paste0(tree.dir, "/", tree.files[i]))
      dataset.name = gsub("_svd.tre$", "", tree.files[i])
      dataset.name = gsub(".*/", "", dataset.name)
      analysis.name = "SVD"
    }

    if (length(grep("Astral", tree.files[i])) == 1) {
      temp.tree = read.tree(paste0(tree.dir, "/", tree.files[i]))
      dataset.name = gsub("_astral.tre$", "", tree.files[i])
      dataset.name = gsub(".*/", "", dataset.name)
      analysis.name = "Astral"

    }#end if

    if (length(grep("Concat", tree.files[i])) == 1) {
      temp.tree = read.tree(paste0(tree.dir, "/", tree.files[i]))
      dataset.name = gsub("_concat_.*$", "", tree.files[i])
      dataset.name = gsub(".*/", "", dataset.name)
      analysis.name = "Concat"

    }#end if

    temp.outgroups = outgroups[outgroups %in% temp.tree$tip.label]
    tree.file = root(temp.tree, temp.outgroups, resolve.root = T)

    #Make a support value table
    if (analysis.name == "Astral"){
      node.no = as.numeric(length(tree.file$tip.label)+1:tree.file$Nnode)

      #Extract node numbers and associated branch lengths
      edge.node = tree.file$edge
      edge.node = data.frame(edge = rep(1:nrow(edge.node)), node1 = edge.node[,1], node2 = edge.node[,2])
      edge.node = edge.node[edge.node$node1 > length(tree.file$tip.label),]
      edge.node = edge.node[edge.node$node2 > length(tree.file$tip.label),]

      #Obtain support
      tree.file$edge.length[is.na(tree.file$edge.length) == T] = 1

      node.vals = stringr::str_split(pattern = ";", tree.file$node.label)
      node.vals = node.vals[-1]
      node.data = as.data.frame(do.call(rbind, node.vals))
      node.data$V1 = gsub("\\[", "", node.data$V1)
      colnames(node.data) = c("q1", "q2", "q3", "f1", 'f2', "f3", "pp1", "pp2", "pp3",
                              "QC", "EN")
      #Adds in node number
      node.data = cbind(node = (length(tree.file$tip.label)+2):(length(tree.file$tip.label)+tree.file$Nnode), node.data)

      #Node data
      node.data$pp1 = as.numeric(gsub("pp1=", "", node.data$pp1))
      node.data$q1 = as.numeric(gsub("q1=", "", node.data$q1))

      node.data = data.frame(Node = node.data$node, Support = node.data$pp1*100)
      node.data$Support[is.na(node.data$Support) == T] = 100
    }#end Astral if

    #Make a support value table
    if (analysis.name == "Concat"){
      node.no = as.numeric(length(tree.file$tip.label)+2:tree.file$Nnode)

      #Extract node numbers and associated branch lengths
      edge.node = tree.file$edge
      edge.node = data.frame(edge = rep(1:nrow(edge.node)), node1 = edge.node[,1], node2 = edge.node[,2])
      edge.node = edge.node[edge.node$node1 > length(tree.file$tip.label),]
      edge.node = edge.node[edge.node$node2 > length(tree.file$tip.label),]

      #Obtain support
      tree.file$edge.length[is.na(tree.file$edge.length) == T] = 1

      node.vals = tree.file$node.label
      node.vals = node.vals[node.vals != ""]
      node.vals = node.vals[-1]

      #Adds in node number
      node.data = data.frame(Node = as.numeric((length(tree.file$tip.label)+3):(length(tree.file$tip.label)+tree.file$Nnode)),
                             Support = as.numeric(node.vals) )
    }#end Concat if

    #Make a support value table
    if (analysis.name == "Jackknife"){
      node.no = as.numeric(length(tree.file$tip.label)+2:tree.file$Nnode)

      #Extract node numbers and associated branch lengths
      edge.node = tree.file$edge
      edge.node = data.frame(edge = rep(1:nrow(edge.node)), node1 = edge.node[,1], node2 = edge.node[,2])
      edge.node = edge.node[edge.node$node1 > length(tree.file$tip.label),]
      edge.node = edge.node[edge.node$node2 > length(tree.file$tip.label),]

      #Obtain support
      tree.file$edge.length[is.na(tree.file$edge.length) == T] = 1

      node.vals = tree.file$node.label
      node.vals = node.vals[node.vals != ""]
      node.vals = node.vals[-1]
      node.vals = as.numeric(node.vals) * 100

      #Adds in node number
      node.data = data.frame(Node = as.numeric((length(tree.file$tip.label)+2):(length(tree.file$tip.label)+tree.file$Nnode)),
                             Support = as.numeric(node.vals) )
    }#end Concat if

    #Make a support value table
    if (analysis.name == "SVD"){
      node.no = as.numeric(length(tree.file$tip.label)+2:tree.file$Nnode)

      #Extract node numbers and associated branch lengths
      edge.node = tree.file$edge
      edge.node = data.frame(edge = rep(1:nrow(edge.node)), node1 = edge.node[,1], node2 = edge.node[,2])
      edge.node = edge.node[edge.node$node1 > length(tree.file$tip.label),]
      edge.node = edge.node[edge.node$node2 > length(tree.file$tip.label),]

      #Obtain support
      tree.file$edge.length[is.na(tree.file$edge.length) == T] = 1

      node.vals = tree.file$node.label
      node.vals = node.vals[node.vals != ""]
      node.vals = node.vals[-1]

      #Adds in node number
      node.data = data.frame(Node = as.numeric((length(tree.file$tip.label)+2):(length(tree.file$tip.label)+tree.file$Nnode)),
                             Support = as.numeric(node.vals) )
    }#end Concat if

    #Saves the data
    set(collect.data, i = as.integer(i), j = match("File_name", header.data), value = tree.files[i] )
    set(collect.data, i = as.integer(i), j = match("Analysis", header.data), value = analysis.name )
    set(collect.data, i = as.integer(i), j = match("Dataset", header.data), value = dataset.name )

    for (j in 1:length(taxa.set)){
      #Monophyly
      mono.p = is.monophyletic(tree.file, taxa.set[[j]])
      if (mono.p == T){
        node.sup = node.data[node.data$Node == getMRCA(tree.file, taxa.set[[j]]),]$Support
      } else { node.sup = 0 }

      set(collect.data, i = as.integer(i),
          j = match(names(taxa.set)[j], header.data), value = node.sup )
    }#end j loop

  }#end i loop

  supp.grid = collect.data[,4:10]
  supp.grid[supp.grid == 0] = 1
  supp.grid[supp.grid >= 95] = 2
  supp.grid[supp.grid >= 70] = 3
  supp.grid[supp.grid >= 50] = 4
  supp.grid[supp.grid >= 5] = 5

  supp.grid = as.matrix(supp.grid)
  rownames(supp.grid) = paste0(collect.data$Analysis, "_", collect.data$Dataset)

  row.colors = collect.data$Analysis
  row.colors[row.colors == "SVD"] = "Pink"
  row.colors[row.colors == "Jackknife"] = "Purple"
  row.colors[row.colors == "Astral"] = "Blue"
  row.colors[row.colors == "Concat"] = "Green"

  supp.grid = supp.grid[seq(dim(supp.grid)[1],1),]

  heatmap(x = supp.grid, Colv = NA, Rowv = NA, scale = "none",
          col = c("black", "blue", "orange", "yellow", "red"), RowSideColors = row.colors)


}#end function
