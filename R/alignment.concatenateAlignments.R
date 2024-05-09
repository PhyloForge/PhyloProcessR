#' @title concatenateAlignments
#'
#' @description Function for concatenating a large number of alignments
#'
#' @param alignment.folder folder that contains aligmnents to be concatenated
#'
#' @param output.name output file name
#'
#' @param partition.file TRUE to save a partition file
#'
#' @param output.format output file format. "all" to save all three types phylip, nexus, and fasta
#'
#' @param partition.format partition file format. "all" to save both raxml and table format
#'
#' @param overwrite TRUE to overwrite file. Default FALSE.
#'
#' @return saves to file concatenated alignments and partition files delimiting the coordinates of each indidividual marker
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
#'
#' @export

concatenateAlignments = function(alignment.folder = NULL,
                                 output.name = NULL,
                                 partition.file = TRUE,
                                 output.format = c("phylip", "nexus", "fasta", "all"),
                                 partition.format = c("raxml", "table", "all"),
                                 overwrite = FALSE) {

  #Debug
#  alignment.folder = paste0(output.directory, "/alignments/")
# output.name = "test-dataset"
#  overwrite = FALSE
# partition.file = TRUE
#  output.format = "fasta"
#  partition.format = "all"

  require(data.table)

  #Parameter checks
  if(is.null(alignment.folder) == TRUE){ stop("Error: a folder of alignments is needed.") }
  if(is.null(output.name) == TRUE){ stop("Error: an output file name is needed.") }
  if(length(output.format) != 1){ stop("Error: please select a single output format or 'all' to save files in all formats.") }
  if(length(partition.format) != 1){ stop("Error: please select a single partition format or 'all' to save files in all formats.") }

  #Check if files exist or not
  if (dir.exists(alignment.folder) == F){
    return(paste0("Directory of alignments could not be found. Exiting."))
  }#end file check

  #Checks output overwrite
  if (overwrite == TRUE){
    if (file.exists(paste0(output.name, ".phy")) == TRUE){ system(paste0("rm ", output.name, ".phy")) }
    if (file.exists(paste0(output.name, ".nex")) == TRUE){ system(paste0("rm ", output.name, ".nex")) }
    if (file.exists(paste0(output.name, ".fa")) == TRUE){ system(paste0("rm ", output.name, ".fa")) }
  } else {
    if (file.exists(paste0(output.name, ".phy")) == TRUE){ return("FIle exists and overwrite == FALSE") }
    if (file.exists(paste0(output.name, ".nex")) == TRUE){ return("FIle exists and overwrite == FALSE") }
    if (file.exists(paste0(output.name, ".fa")) == TRUE){ return("FIle exists and overwrite == FALSE") }
  }#end overwrite if

  #Gets list of alignments
  align.files = list.files(alignment.folder, full.names = T)
  align.list = lapply(align.files, function (x) data.table::fread(x, header = T))
  sample.names = unique(unlist(lapply(align.list, function (x) x[,1])))
  sample.names = sample.names[order(sample.names)]

  #sets up data collection for proportion
  feat.headers = c("Marker", "Length", "Start", "End")
  feature.data = data.table::data.table(matrix(as.numeric(0),
                                   nrow = length(align.files),
                                   ncol = length(feat.headers)))
  data.table::setnames(feature.data, feat.headers)
  feature.data[, Marker:=as.character(Marker)]

  #Gets the names and number of columns
  names(align.list) = gsub(".*/", "", align.files)
  n.cols = length(align.list)

  #For the concatenated alignment
  concat.headers = c("Sample", names(align.list))
  concat.data = data.table::data.table(matrix(as.character("n"), nrow = length(sample.names),
                            ncol = length(concat.headers)))
  data.table::setnames(concat.data, concat.headers)
  concat.data[, Sample:=as.character(sample.names)]

  # Add columns with generated name with standardized length
  start = 1
  for (x in 2:ncol(concat.data)) {

    align = data.table::data.table(align.list[[x-1]])
    align.len = as.numeric(colnames(align)[2])
    colnames(align) = c("Sample", "Sequence")

    #This fills in missing samples
    current.samples = unlist(c(align[,1]))
    names(current.samples) = NULL
    miss.samples = sample.names[!sample.names %in% current.samples]

    #Fill in missing samples
    if (length(miss.samples) != 0){
      #Adds Ns in for each missing sample
      for (j in 1:length(miss.samples)){
        temp.seq = paste0(rep("n", align.len), collapse = "")
        save.seq = data.table::data.table(Sample = miss.samples[j],
                              Sequence = temp.seq)
        align = rbind(align, save.seq)
      }#end j loop
    } #end if

    align = align[order(align$Sample),]

    data.table::set(concat.data, i = match(concat.data$Sample, align$Sample),
        j = as.integer(x), value = align$Sequence)

    align.name = colnames(concat.data)[x]

    if (start == 1){ end = start + align.len - 1 } else { end = start + align.len - 1 }

    #Saves location data
    data.table::set(feature.data, i = as.integer(x-1), j = match("Marker", feat.headers), value =  align.name)
    data.table::set(feature.data, i = as.integer(x-1), j = match("Length", feat.headers), value =  align.len)
    data.table::set(feature.data, i = as.integer(x-1), j = match("Start", feat.headers), value =  start)
    data.table::set(feature.data, i = as.integer(x-1), j = match("End", feat.headers), value =  end)
    start = end + 1
  }#end x loop

  ############
  ### Turn data into a phylip format to then save differently
  ############

  #### SAVE AS FASTA
  ###################
  if (output.format == "fasta" | output.format == "fa" | output.format == "fas" | output.format == "all"){
    #Fasta easy
    fileConn = file(paste0(output.name, ".fa"), open = "a")
    #Saves each line
    for (x in 1:nrow(concat.data)){
      writeLines(paste0(">", concat.data$Sample[x]), con = fileConn, sep = "\n")
      sequence.data = paste0(concat.data[x,2:ncol(concat.data)], collapse = "")
      writeLines(sequence.data, con = fileConn, sep = "\n")
    }#end x loop
    close(fileConn)
    print(paste0(output.name, ".fa concatenated alignment written in fasta format."))
  }#End fasta if

  #### SAVE AS NEXUS
  ###################
  if (output.format == "nexus" | output.format == "nex" | output.format == "all"){

    #Gets header information
    ntax = nrow(concat.data)
    nchar = sum(feature.data$Length)
    nex.header = paste0("dimensions ntax=", ntax, " nchar=", nchar, ";")
    for.header = paste0("format datatype=dna missing=n gap=-;")

    #Prep sample names to all have the same length padded with spaces
    name.length = max(nchar(sample.names)) + 1
    for (x in 1:length(sample.names)){
      temp.name = concat.data$Sample[x]
      space.pad = name.length - nchar(temp.name)
      space.add = paste(rep(" ", space.pad), collapse = "")
      new.name = paste0(temp.name, space.add, collapse = "")
      concat.data$Sample[x] = new.name
    }#end x

    #Start saving file
    fileConn = file(paste0(output.name, ".nex"), open = "a")
    writeLines("#NEXUS", con = fileConn, sep = "\n")
    writeLines("begin data;", con = fileConn, sep = "\n")
    writeLines(paste0("\t", nex.header), con = fileConn, sep = "\n")
    writeLines(paste0("\t", for.header), con = fileConn, sep = "\n")
    writeLines("matrix", con = fileConn, sep = "\n")

    #Saves each line
    for (x in 1:nrow(concat.data)){
      sequence.data = paste0(concat.data[x,], collapse = "")
      writeLines(sequence.data, con = fileConn, sep = "\n")
    }#end x loop

    #Closes file
    writeLines(";", con = fileConn, sep = "\n")
    writeLines("end;", con = fileConn, sep = "\n")
    close(fileConn)
    #Done
    print(paste0(output.name, ".nex concatenated alignment written in nexus format."))

  }#End nexus if

  #### SAVE AS PHYLIP
  ###################
  if (output.format == "phylip" | output.format == "phy" | output.format == "all"){
    #Now have to reformat
    ntax = nrow(concat.data)
    nchar = sum(feature.data$Length)
    phy.header = paste(ntax, nchar)

    #Prep sample names to all have the same length padded with spaces
    name.length = max(nchar(sample.names)) + 4
    for (x in 1:length(sample.names)){
      temp.name = concat.data$Sample[x]
      space.pad = name.length - nchar(temp.name)
      space.add = paste(rep(" ", space.pad), collapse = "")
      new.name = paste0(temp.name, space.add, collapse = "")
      concat.data$Sample[x] = new.name
    }#end x

    fileConn = file(paste0(output.name, ".phy"), open = "a")
    writeLines(phy.header, con = fileConn, sep = "\n")
    #Saves each line
    for (x in 1:nrow(concat.data)){
      sequence.data = paste0(concat.data[x,], collapse = "")
      writeLines(sequence.data, con = fileConn, sep = "\n")
    }#end x loop
    close(fileConn)
    print(paste0(output.name, ".phy concatenated alignment written in phylip format."))
  }#End phylip if

  ####################################################################
  ##### Saves partition files
  ####################################################################
  if (partition.file == TRUE){

    #Saves in a simple table format
    if (length(partition.format[partition.format == "table"]) == 1){
      #Checks to overwrite
      if (file.exists(paste0(output.name, "_partitions_table.txt")) == TRUE){
        unlink(paste0(output.name, "_partitions_table.txt"))}

      write.table(feature.data, file = paste0( output.name, "_partitions_table.txt"),
                  sep = "\t", row.names = FALSE, quote = FALSE)

      print(paste0(output.name, "_partitions_table.txt written in table format."))
    }#ends table if

    #Save in raxml format
    if (length(partition.format[partition.format == "raxml"]) == 1){
      #Checks to overwrite
      if (file.exists(paste0(output.name, "_partitions_raxml.txt")) == TRUE){
        unlink(paste0(output.name, "_partitions_raxml.txt"))}

      fileConn = file(paste0(output.name, "_partitions_raxml.txt"), open = "a")

      for (x in 1:nrow(feature.data)){

        writeLines(paste0("DNA, ", feature.data$Marker[x], " = ",
                          feature.data$Start[x], "-", feature.data$End[x]),
                   con = fileConn, sep = "\n")

      }#end x loop

      close(fileConn)

      print(paste0(output.name, "_partitions_raxml.txt written in raxml format."))

    }#end raxml if
  }#end partition if

}#end function

