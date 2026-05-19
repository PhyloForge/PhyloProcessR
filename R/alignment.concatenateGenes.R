#' @title concatenateGenes
#'
#' @description Concatenates exon-level alignments into gene-level alignments using a
#' metadata table that maps each exon (marker) to a gene. For each gene, all matching
#' exon alignments are joined in order, with missing taxa filled in with N characters.
#' Genes represented by fewer than \code{minimum.exons} present exon alignments are
#' skipped. Output alignments are written in the chosen format and processing is
#' parallelised across genes.
#'
#' @param alignment.folder path to the directory containing per-exon alignment files.
#'
#' @param output.folder path to the directory where concatenated gene alignment files
#' will be saved.
#'
#' @param feature.gene.names path to a tab-delimited metadata file with at minimum columns
#' named \code{marker} and \code{gene}, mapping each exon alignment name to a gene name.
#'
#' @param find.exon.names logical. If TRUE, exon-to-gene associations are inferred
#' automatically rather than read from \code{feature.gene.names}. Default FALSE.
#'
#' @param input.format format of the input exon alignment files. Accepted values: "phylip",
#' "nexus", or "fasta".
#'
#' @param output.format format for the output gene alignment files. Accepted values:
#' "phylip", "nexus", or "fasta".
#'
#' @param minimum.exons minimum number of exon alignments that must be present for a gene
#' to be concatenated and saved. Genes with fewer exons are skipped. Default 2.
#'
#' @param remove.reverse logical. If TRUE, the leading \code{_R_} tag added by MAFFT to
#' reverse-complemented sequences is stripped from sequence names before processing.
#' Default FALSE.
#'
#' @param remove.duplicates logical. If TRUE, exon alignments containing duplicate sample
#' names are skipped rather than causing an error. Default FALSE.
#'
#' @param overwrite logical. If TRUE, previously completed gene alignments are overwritten;
#' if FALSE, they are skipped. Default FALSE.
#'
#' @param threads number of parallel threads to use. Default 1.
#'
#' @param memory total memory (in GB) to allocate across all threads. Default 1.
#'
#' @return Writes concatenated gene alignment files to \code{output.folder}. No value is
#' returned to R.
#'
#' @export

concatenateGenes = function(alignment.folder = NULL,
                            output.folder = NULL,
                            feature.gene.names = NULL,
                            find.exon.names = FALSE,
                            input.format = c("phylip", "nexus", "fasta"),
                            output.format = c("phylip", "nexus", "fasta"),
                            minimum.exons = 2,
                            remove.reverse = FALSE,
                            remove.duplicates = FALSE,
                            overwrite = FALSE,
                            threads = 1,
                            memory = 1) {

  # #Debug
  # work.dir = "/Volumes/LaCie/data-analysis/alignments"
  # setwd(work.dir)
  # alignment.folder = "untrimmed_all-markers"
  # output.folder = "untrimmed_genes"
  # input.format = "phylip"
  # output.format = "phylip"
  # remove.reverse = FALSE
  # overwrite = FALSE
  # feature.gene.names = "/Volumes/LaCie/data-analysis/gene_metadata.txt"
  # minimum.exons = 2
  # threads = 8
  # memory = 24

  input.format  = match.arg(input.format)
  output.format = match.arg(output.format)

  #Parameter checks
  if(is.null(alignment.folder) == TRUE){ stop("Error: a folder of alignments is needed.") }
  if(is.null(output.folder) == TRUE){ stop("Error: an output file name is needed.") }
  if(is.null(feature.gene.names) == TRUE && find.exon.names == FALSE){ stop("Error: a table associating each exon with a gene is needed.") }

  #Check if files exist or not
  if (dir.exists(alignment.folder) == F){
    return(paste0("Directory of alignments could not be found. Exiting."))
  }#end file check

  #Checks output overwrite
  if (overwrite == TRUE){
    if (file.exists(paste0(output.folder)) == TRUE){ system(paste0("rm -r ", output.folder)) }
    dir.create(output.folder)
  } else {
    if (!dir.exists(output.folder)) { dir.create(output.folder) }
  }#end overwrite if

  #Gets list of alignments
  align.files = list.files(alignment.folder, full.names = FALSE, recursive = TRUE)
  exon.data = data.table::fread(file = feature.gene.names, header = TRUE)
  gene.names = unique(exon.data$gene)
  gene.names = gene.names[is.na(gene.names) != TRUE]

  #Skips files done already if resume = TRUE
  if (overwrite == FALSE){
    done.files = list.files(output.folder)
    done.genes = gsub("\\..*$", "", done.files)
    gene.names = gene.names[!gene.names %in% done.genes]
  }

  mem.cl = floor(memory/threads)

  parallel::mclapply(seq_along(gene.names), function(i) {
  tryCatch({
    #Find exon data for this gene
    gene.data = exon.data[exon.data$gene %in% gene.names[i],]
    #Match to the files to obtain
    gene.exons = gene.data$marker

    temp.dir = paste0(output.folder, "/temp-", gene.names[i])
    if (dir.exists(temp.dir) == TRUE){
      system(paste0("rm -r ", temp.dir))
      dir.create(temp.dir)
    } else { dir.create(temp.dir) }

    save.names = c()
    save.length = 0
    exon.count = 0
    for (y in 1:length(gene.exons)){
      # Finds alignment name
      temp.exon = align.files[gsub("\\..*", "", align.files) %in% gene.exons[y]]

      if (length(temp.exon) >= 2){
        stop("more than one exon matches to spreadsheet entry.")
      }

      if (length(temp.exon) == 0) {  next }

      #Load in alignments
      if (input.format == "phylip"){
        exon.align = Biostrings::readDNAMultipleAlignment(file = paste0(alignment.folder, "/", temp.exon), format = "phylip")
        exon.align = Biostrings::DNAStringSet(exon.align)
      }#end phylip

      if (input.format == "fasta"){
        exon.align = Biostrings::readDNAStringSet(paste0(alignment.folder, "/", temp.exon) )
      }#end phylip

      #Removes the reverse from mafft if present
      if (remove.reverse == TRUE){ names(exon.align) = gsub("^_R_", "", names(exon.align) ) }

      #Checks for duplicates
      dup.names <- exon.align[duplicated(names(exon.align)), ]
      if (length(dup.names) != 0) {
        if (remove.duplicates == TRUE) {
          next
        } else {
          stop("DUPLICATE names found. Please check alignments.")
        }#end else
      }#end if

      #Gathers sample names
      exon.save = paste0(temp.dir, "/", gene.names[i], "-", y, ".phy")
      save.names = append(save.names, exon.save)
      #Saves alignment as phylip
      write.align = PhyloProcessR::alignmentConversion(input.alignment = exon.align, end.format = "matrix")
      PhyloProcessR::writePhylip(write.align, file = exon.save, interleave = F )
      save.length = save.length + Biostrings::width(exon.align)[1]
      exon.count = exon.count + 1
    }#end y loop

    if (exon.count < minimum.exons){
      print(paste0("Below minimum exon count of ", minimum.exons, ". Skipped gene."))
      system(paste0("rm -r ", temp.dir))
      return(NULL)
    }

    if (length(save.names) > 0) {

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

        end = start + align.len - 1

        start = end + 1
      }#end x loop

      output.name = paste0(output.folder, "/", gene.names[i])

      #### SAVE AS FASTA
      ###################
      if (output.format == "fasta" || output.format == "fa" || output.format == "fas"){
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
      if (output.format == "nexus" || output.format == "nex"){

        #Gets header information
        ntax = nrow(concat.data)
        n.char = save.length
        nex.header = paste0("dimensions ntax=", ntax, " nchar=", n.char, ";")
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
      if (output.format == "phylip" || output.format == "phy"){
        #Now have to reformat
        ntax = nrow(concat.data)
        n.char = save.length
        phy.header = paste0(" ", ntax, " ", n.char)

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
    }#end if way above

    system(paste0("rm -r ", temp.dir))

    rm(align.list, concat.data, exon.align)
    gc()

  }, error = function(e) {
    warning(gene.names[i], " failed: ", conditionMessage(e))
  })
  }, mc.cores = threads) #end i loop

}#end function
