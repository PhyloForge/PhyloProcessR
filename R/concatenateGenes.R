#' @title concatenateGenes
#'
#' @description Function for concatenating a large number of alignments
#'
#' @param alignment.folder folder that contains aligmnents to be concatenated
#'
#' @param output.folder output file name
#'
#' @param exon.gene.names output file name
#'
#' @param input.format input file format. Save three types: phylip, nexus, and fasta
#'
#' @param output.format output file format. Save three types: phylip, nexus, and fasta
#'
#' @param remove.reverse TRUE to remove "_R_" placed before reversed sequences in some alignments. Default FALSE.
#'
#' @param overwrite TRUE to overwrite file. Default FALSE.
#'
#' @param threads TRUE to overwrite file. Default FALSE.
#'
#' @param memory TRUE to overwrite file. Default FALSE.
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

concatenateGenes = function(alignment.folder = NULL,
                            output.folder = NULL,
                            feature.gene.names = NULL,
                            input.format = c("phylip", "nexus", "fasta"),
                            output.format = c("phylip", "nexus", "fasta"),
                            remove.reverse = FALSE,
                            overwrite = FALSE,
                            resume = TRUE,
                            threads = 1,
                            memory = 1) {

  #Debug
  # library(foreach)
  # setwd("/Volumes/Rodents/Murinae/Trimming")
  # alignment.folder = "/Volumes/Rodents/Murinae/Trimming"
  # output.folder = "01_emily-subset-genes"
  # overwrite = TRUE
  # input.format = "fasta"
  # output.format = "phylip"
  # feature.gene.names = "/Volumes/Rodents/Murinae/Selected_Transcripts/Mus_gene_metadata.csv"
  # remove.reverse = TRUE
  # threads = 6
  # memory = 6

  #Parameter checks
  if(is.null(alignment.folder) == TRUE){ stop("Error: a folder of alignments is needed.") }
  if(is.null(output.folder) == TRUE){ stop("Error: an output file name is needed.") }
  if(length(output.format) != 1){ stop("Error: please select a single output format or 'all' to save files in all formats.") }
  if(is.null(feature.gene.names) == TRUE){ stop("Error: a table associating each exon with a gene is needed.") }

  #Check if files exist or not
  if (dir.exists(alignment.folder) == F){
    return(paste0("Directory of alignments could not be found. Exiting."))
  }#end file check

  #Checks output overwrite
  if (overwrite == TRUE){
    if (file.exists(paste0(output.folder)) == TRUE){ system(paste0("rm -r ", output.folder)) }
    dir.create(output.folder)
  } else {
    dir.create(output.folder)
  }#end overwrite if

  #Gets list of alignments
  align.files = list.files(alignment.folder, full.names = F, recursive = T)
  exon.data = data.table::fread(file = feature.gene.names, header = T)
  gene.names = unique(exon.data$gene_id)

  length(align.files)

  #Skips files done already if resume = TRUE
  if (resume == TRUE){
    done.files = list.files(output.folder)
    align.files = align.files[!gsub("\\..*", "", align.files) %in% gsub("\\..*", "", done.files)]
  }

  length(done.files)
  length(align.files)
  # #Sets up multiprocessing
   #cl = parallel::makeCluster(threads, outfile = "")
   #doParallel::registerDoParallel(cl)
   #mem.cl = floor(memory/threads)

  # foreach::foreach(i=1:length(gene.names), .combine = rbind, .packages = c("PHYLOCAP", "foreach", "Biostrings","Rsamtools", "ape", "stringr", "data.table")) %dopar% {
  for (i in 1:length(gene.names)){
    #Find exon data for this gene
    gene.data = exon.data[exon.data$gene_id %in% gene.names[i],]
    #Match to the files to obtain
    gene.exons = gene.data$target_name

    temp.dir = paste0(output.folder, "/temp-", gene.names[i])
    if (dir.exists(temp.dir) == TRUE){
      system(paste0("rm -r ", temp.dir))
      dir.create(temp.dir)
    } else { dir.create(temp.dir) }

    save.names = c()
    save.length = 0
    for (y in 1:length(gene.exons)){
      #Finds alignment name
      temp.exon = align.files[grep(gene.exons[y], align.files)]

      if (length(temp.exon) >= 2){
        temp.exon = temp.exon[gsub("-mafft.fa", "", temp.exon) %in% gene.exons[y]]
      }

      if (length(temp.exon) == 0) {  next }

      #Load in alignments
      if (input.format == "phylip"){
        exon.align = Biostrings::readAAMultipleAlignment(file = paste0(alignment.folder, "/", temp.exon), format = "phylip")
        exon.align = Biostrings::DNAStringSet(exon.align)
      }#end phylip

      if (input.format == "fasta"){
        exon.align = Biostrings::readDNAStringSet(paste0(alignment.folder, "/", temp.exon) )
      }#end phylip

      #Removes the reverse from mafft if present
      if (remove.reverse == TRUE){ names(exon.align) = gsub("^_R_", "", names(exon.align) ) }
      #Gathers sample names
      exon.save = paste0(temp.dir, "/", gene.names[i], "-", y, ".phy")
      save.names = append(save.names, exon.save)
      #Saves alignment as phylip
      write.align = alignmentConversion(input.alignment = exon.align, end.format = "matrix")
      writePhylip(write.align, file = exon.save, interleave = F )
      save.length = save.length + Biostrings::width(exon.align)[1]
    }#end y loop

    if (is.null(save.names) == FALSE) {

      #Loads in the data
      align.list = lapply(save.names, function (x) data.table::fread(x, header = T))
      sample.names = unique(unlist(lapply(align.list, function (x) x[,1])))
      sample.names = sample.names[order(sample.names)]

      #Gets the names and number of columns
      names(align.list) = gsub(".*/", "", save.names)
      n.cols = length(align.list)

      #For the concatenated alignment
      concat.headers = c("Sample", names(align.list))
      concat.data = data.table::data.table(matrix(as.character("n"), nrow = length(sample.names),
                                                  ncol = length(concat.headers)))
      data.table::setnames(concat.data, concat.headers)
      concat.data[, Sample:=as.character(sample.names)]

      # Add columns with generated name with standardized length
      start = 1
      for (z in 2:ncol(concat.data)) {

        align = data.table::data.table(align.list[[z-1]])
        align.len = as.numeric(colnames(align)[2])
        colnames(align) = c("Sample", "Sequence")

        #This fills in missing samples
        current.samples = unlist(c(align[,1]))
        names(current.samples) = NULL
        miss.samples = sample.names[!sample.names %in% current.samples]

        #Fill in missing samples
        if (length(miss.samples) != 0){
          #Adds Ns in for each missing sample
          for (y in 1:length(miss.samples)){
            temp.seq = paste0(rep("n", align.len), collapse = "")
            save.seq = data.table::data.table(Sample = miss.samples[y],
                                              Sequence = temp.seq)
            align = rbind(align, save.seq)
          }#end j loop
        } #end if

        align = align[order(align$Sample),]

        data.table::set(concat.data, i = match(concat.data$Sample, align$Sample),
                        j = as.integer(z), value = align$Sequence)

        align.name = colnames(concat.data)[z]

        if (start == 1){ end = start + align.len - 1 } else { end = start + align.len - 1 }

        start = end + 1
      }#end x loop

      ###########
      ## Turn data into a phylip format to then save differently
      ###########
      # if (remove.reverse == TRUE){ names(align) = gsub("^_R_", "", names(align) ) }
      #
      # align.list = lapply(gene.exons, function (x) data.table::fread(x, header = T))
      # sample.names = unique(unlist(lapply(align.list, function (x) x[,1])))
      # sample.names = sample.names[order(sample.names)]

      # sample.names = c()
      # for (j in 1:length(gene.exons)){
      #   #Reads in files
      #   temp.exon = align.files[grep(gene.exons[j], align.files)]
      #   align = Biostrings::readDNAStringSet(paste0(alignment.folder, "/", temp.exon) )
      #   sample.names = append(sample.names, names(align))
      # }#end j
      #
      # sample.names = unique(gsub("^_R_", "", sample.names) )
      # concat.data = Biostrings::DNAStringSet()
      # for (j in 1:length(gene.exons)){
      #   #Reads in files
      #   temp.exon = align.files[grep(gene.exons[j], align.files)]
      #   align = Biostrings::readDNAStringSet(paste0(alignment.folder, "/", temp.exon) )
      #   if (remove.reverse == TRUE){ names(align) = gsub("^_R_", "", names(align) ) }
      #
      #   #Adds blanks
      #   add.taxa = sample.names[!sample.names %in% names(align)]
      #   blank.align = Biostrings::DNAStringSet()
      #   if (length(add.taxa) != 0){
      #     for (y in 1:length(add.taxa)){
      #       blank.align = append(blank.align,
      #                            Biostrings::DNAStringSet(paste0(rep("-", max(Biostrings::width(align))), collapse = "")) )
      #     }
      #     names(blank.align) = add.taxa
      #   }#end rem seqs if
      #
      #   #Saves the slices and cats
      #   new.align = append(align, blank.align)
      #   new.align = new.align[order(names(new.align))]
      #   save.names = names(new.align)
      #   concat.data = Biostrings::DNAStringSet(paste0(as.character(concat.data), as.character(new.align)))
      #   names(concat.data) = save.names
      # }#end j loop

      output.name = paste0(output.folder, "/", gene.names[i])
      # write.align = alignmentConversion(input.alignment = concat.data, end.format = "matrix")
      # writePhylip(write.align, file = paste0(output.name, ".phy"), interleave = F )

      #### SAVE AS FASTA
      ###################
      if (output.format == "fasta" | output.format == "fa" | output.format == "fas"){
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
      if (output.format == "nexus" | output.format == "nex"){

        #Gets header information
        ntax = nrow(concat.data)
        nchar = save.length
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
      if (output.format == "phylip" | output.format == "phy"){
        #Now have to reformat
        ntax = nrow(concat.data)
        nchar = save.length
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

      system(paste0("rm -r ", temp.dir))

    }#end if way above

   }#end i loop

  parallel::stopCluster(cl)

  system(paste0("rm -r ", output.folder, "/temp*"))

}#end function

