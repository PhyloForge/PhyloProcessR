#' @title batchTrimAlignments
#'
#' @description Function for batch trimming a folder of alignments, with the various trimming functions available to select from
#'
#' @param alignment.dir path to a folder of sequence alignments in phylip format.
#'
#' @param alignment.format available input alignment formats: fasta or phylip
#'
#' @param output.dir contigs are added into existing alignment if algorithm is "add"
#'
#' @param output.format available output formats: phylip
#'
#' @param HmmCleaner algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param HmmCleaner.path TRUE applies the adjust sequence direction function of MAFFT
#'
#' @param TrimAl if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param TrimAl.path path to a folder of sequence alignments in phylip format.
#'
#' @param trim.external give a save name if you wnat to save the summary to file.
#'
#' @param min.external.percent TRUE to supress mafft screen output
#'
#' @param trim.coverage path to a folder of sequence alignments in phylip format.
#'
#' @param min.coverage.percent contigs are added into existing alignment if algorithm is "add"
#'
#' @param trim.column algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param min.column.gap.percent TRUE applies the adjust sequence direction function of MAFFT
#'
#' @param alignment.assess if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param min.sample.bp path to a folder of sequence alignments in phylip format.
#'
#' @param min.align.length give a save name if you wnat to save the summary to file.
#'
#' @param min.taxa.count TRUE to supress mafft screen output
#'
#' @param min.gap.percent if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param threads path to a folder of sequence alignments in phylip format.
#'
#' @param mem give a save name if you wnat to save the summary to file.
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

batchTrimAlignments = function(alignment.dir = NULL,
                               alignment.format = "phylip",
                               output.dir = NULL,
                               output.format = "phylip",
                               PreQual = TRUE,
                               PreQual.path = "prequal",
                               HmmCleaner = FALSE,
                               HmmCleaner.path = "HmmCleaner.pl",
                               TrimAl = TRUE,
                               TrimAl.path = "trimal",
                               trim.external = TRUE,
                               min.external.percent = 50,
                               trim.coverage = TRUE,
                               min.coverage.percent = 50,
                               trim.column = TRUE,
                               min.column.gap.percent = 100,
                               alignment.assess = TRUE,
                               min.sample.bp = 0,
                               min.align.length = 0,
                               min.taxa.count = 0,
                               min.gap.percent = 0,
                               threads = 1,
                               mem = 8,
                               overwrite = FALSE,
                               resume = TRUE) {


  # work.dir = "/Volumes/Rodents/Murinae/Trimming"
  # align.dir = "/Volumes/Rodents/Murinae/Trimming/01_emily-subset-mafft"
  # setwd(work.dir)
  # alignment.dir = align.dir
  # alignment.format = "fasta"
  # output.format = "phylip"
  # output.dir = paste0(work.dir, "/01_emily_trimmed")
  # TrimAl = TRUE
  # PreQual = TRUE
  # HmmCleaner = TRUE
  # trim.column = TRUE
  # alignment.assess = TRUE
  # trim.external = TRUE
  # trim.coverage = TRUE
  # min.coverage.percent = 30
  # min.external.percent = 50
  # min.column.gap.percent = 100
  # overwrite = FALSE
  # min.align.length = 60
  # min.taxa.count = 3
  # min.gap.percent = 50
  # min.sample.bp = 40
  # threads = 4
  # mem = 8
  # resume = TRUE

  if (alignment.dir == output.dir){ stop("You should not overwrite the original alignments.") }

  if (dir.exists(output.dir) == FALSE) { dir.create(output.dir) }

  #So I don't accidentally delete everything while testing resume
  if (resume == TRUE & overwrite == TRUE){
    overwrite = FALSE
    stop("Error: resume = T and overwrite = T, cannot resume if you are going to delete everything!")
  }

  if (dir.exists(output.dir) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.dir))
      dir.create(output.dir)
     }
  }#end dir exists

  #Gathers alignments
  align.files = list.files(alignment.dir)
  #Skips files done already if resume = TRUE
  if (resume == TRUE){
    done.files = list.files(output.dir)
    align.files = align.files[!gsub("\\..*", "", align.files) %in% gsub("\\..*", "", done.files)]
  }

  #Data to collect
  header.data = c("Alignment", "Pass", "startSamples", "hmmSamples", "trimalSamples",
                  "edgeSamples", "covSamples", "columnSamples",
                  "startLength", "hmmLength", "trimalLength",
                  "edgeLength", "covLength", "columnLength",
                  "startBasepairs", "hmmBasepairs", "trimalBasepairs",
                  "edgeBasepairs", "covBasepairs", "columnBasepairs",
                  "startGaps", "hmmGaps", "trimalGaps",
                  "edgeGaps", "covGaps", "columnGaps",
                  "startPerGaps", "hmmPerGaps", "trimalPerGaps",
                  "edgePerGaps", "covPerGaps", "columnPerGaps")

  save.data = data.table::data.table(matrix(as.double(0), nrow = length(align.files), ncol = length(header.data)))
  data.table::setnames(save.data, header.data)
  save.data[, Alignment:=as.character(Alignment)]
  save.data[, Pass:=as.logical(Pass)]

  #Sets up multiprocessing
  #cl = parallel::makeCluster(threads, outfile = "")
  #doParallel::registerDoParallel(cl)
  #mem.cl = floor(mem/threads)

  #Loops through each locus and does operations on them
  #out.data = foreach(i=1:length(align.files), .combine = rbind, .packages = c("PHYLOCAP", "foreach", "Biostrings","Rsamtools", "ape", "stringr", "data.table")) %dopar% {
  for (i in 1:length(align.files)){
    print(paste0(align.files[i], " Started..."))

     #Load in alignments
    if (alignment.format == "phylip"){
      align = Biostrings::readAAMultipleAlignment(file = paste0(alignment.dir, "/", align.files[i]), format = "phylip")
      align = Biostrings::DNAStringSet(align)
      save.name = gsub(".phy$", "", align.files[i])
      save.name = gsub(".phylip$", "", save.name)
    }#end phylip

    if (alignment.format == "fasta"){
      align = Rsamtools::scanFa(Rsamtools::FaFile(paste0(alignment.dir, "/", align.files[i])))   # loads up fasta file
      save.name = gsub(".fa$", "", align.files[i])
      save.name = gsub(".fasta$", "", save.name)
    }#end phylip

    temp.data = save.data[1,]

    # Runs the functions
    #######
    #Step 1: Strip Ns
    non.align = replaceAlignmentCharacter(alignment = align,
                                          char.find = "N",
                                          char.replace = "-")

    #Summarize all this, no functoin
    data.table::set(temp.data, i = as.integer(1), j = match("Alignment", header.data), value = save.name)
    data.table::set(temp.data, i = as.integer(1), j = match("startSamples", header.data), value = length(non.align))
    data.table::set(temp.data, i = as.integer(1), j = match("startLength", header.data), value = Biostrings::width(non.align)[1] )
    gap.count = countAlignmentGaps(non.align)
    data.table::set(temp.data, i = as.integer(1), j = match("startBasepairs", header.data), value = gap.count[2] - gap.count[1])
    data.table::set(temp.data, i = as.integer(1), j = match("startGaps", header.data), value = gap.count[1])
    data.table::set(temp.data, i = as.integer(1), j = match("startPerGaps", header.data), value = gap.count[3])

    #Step 2. Runs hmmCleaner
    # if (PreQual == TRUE){
    #
    #   preql.align = trimPreQual(alignment = non.align,
    #                             specificity = TRUE,
    #                             large = FALSE,
    #                             hmmcleaner.path = "HmmCleaner.pl",
    #                             quiet = T,
    #                             delete.temp = T)
    #   non.align = preql.align
    #   #Saves the data
    #   data.table::set(temp.data, i = as.integer(1), j = match("hmmSamples", header.data), value = length(hmmcl.align))
    #   data.table::set(temp.data, i = as.integer(1), j = match("hmmLength", header.data), value = Biostrings::width(hmmcl.align)[1])
    #   gap.count = countAlignmentGaps(non.align)
    #   data.table::set(temp.data, i = as.integer(1), j = match("hmmBasepairs", header.data), value = gap.count[2] - gap.count[1])
    #   data.table::set(temp.data, i = as.integer(1), j = match("hmmGaps", header.data), value = gap.count[1])
    #   data.table::set(temp.data, i = as.integer(1), j = match("hmmPerGaps", header.data), value = gap.count[3])
    # }#end if


    #Step 2. Runs hmmCleaner
    if (HmmCleaner == TRUE){

      hmmcl.align = trimHmmCleaner(alignment = non.align,
                                   specificity = TRUE,
                                   large = FALSE,
                                   hmmcleaner.path = "HmmCleaner.pl",
                                   quiet = T,
                                   delete.temp = T)
      non.align = hmmcl.align
      #Saves the data
      data.table::set(temp.data, i = as.integer(1), j = match("hmmSamples", header.data), value = length(hmmcl.align))
      data.table::set(temp.data, i = as.integer(1), j = match("hmmLength", header.data), value = Biostrings::width(hmmcl.align)[1])
      gap.count = countAlignmentGaps(non.align)
      data.table::set(temp.data, i = as.integer(1), j = match("hmmBasepairs", header.data), value = gap.count[2] - gap.count[1])
      data.table::set(temp.data, i = as.integer(1), j = match("hmmGaps", header.data), value = gap.count[1])
      data.table::set(temp.data, i = as.integer(1), j = match("hmmPerGaps", header.data), value = gap.count[3])
    }#end if

    #Step 3. Trimal trimming
    if (TrimAl == TRUE){

      trimal.align = trimTrimal(alignment = non.align,
                                quiet = TRUE)
      non.align = trimal.align
      #Saves stat data
      data.table::set(temp.data, i = as.integer(1), j = match("trimalSamples", header.data), value = length(trimal.align))
      data.table::set(temp.data, i = as.integer(1), j = match("trimalLength", header.data), value = Biostrings::width(trimal.align)[1])
      gap.count = countAlignmentGaps(non.align)
      data.table::set(temp.data, i = as.integer(1), j = match("trimalBasepairs", header.data), value = gap.count[2] - gap.count[1])
      data.table::set(temp.data, i = as.integer(1), j = match("trimalGaps", header.data), value = gap.count[1])
      data.table::set(temp.data, i = as.integer(1), j = match("trimalPerGaps", header.data), value = gap.count[3])
    }#end if

    # Step 4. Edge trimming
    if (trim.external == TRUE){
      #external trimming function
      edge.align = trimExternal(alignment = non.align,
                                min.n.seq = ceiling(length(non.align) * (min.external.percent/100)),
                                codon.trim = F)
      non.align = edge.align
      #Saves stat data
      data.table::set(temp.data, i = as.integer(1), j = match("edgeSamples", header.data), value = length(edge.align))
      data.table::set(temp.data, i = as.integer(1), j = match("edgeLength", header.data), value = Biostrings::width(edge.align)[1])
      gap.count = countAlignmentGaps(non.align)
      data.table::set(temp.data, i = as.integer(1), j = match("edgeBasepairs", header.data), value = gap.count[2] - gap.count[1])
      data.table::set(temp.data, i = as.integer(1), j = match("edgeGaps", header.data), value = gap.count[1])
      data.table::set(temp.data, i = as.integer(1), j = match("edgePerGaps", header.data), value = gap.count[3])
    }#end trim external

    #Step 5. Evaluate and cut out each sample
    if (trim.coverage == TRUE){
      #sample coverage function
      cov.align = trimSampleCoverage(alignment = non.align,
                                     min.coverage.percent = min.coverage.percent,
                                     min.sample.bp = min.sample.bp)
      non.align = cov.align
      #Saves stat data
      data.table::set(temp.data, i = as.integer(1), j = match("covSamples", header.data), value = length(cov.align))
      data.table::set(temp.data, i = as.integer(1), j = match("covLength", header.data), value = Biostrings::width(cov.align)[1])
      gap.count = countAlignmentGaps(non.align)
      data.table::set(temp.data, i = as.integer(1), j = match("covBasepairs", header.data), value = gap.count[2] - gap.count[1])
      data.table::set(temp.data, i = as.integer(1), j = match("covGaps", header.data), value = gap.count[1])
      data.table::set(temp.data, i = as.integer(1), j = match("covPerGaps", header.data), value = gap.count[3])
    }#end trim.external

    if (trim.column == TRUE){
      #Trim alignment colums
      col.align = trimAlignmentColumns(alignment = non.align,
                                       min.gap.percent = min.column.gap.percent)
      non.align = col.align
      #Saves stat data
      data.table::set(temp.data, i = as.integer(1), j = match("columnSamples", header.data), value = length(col.align))
      data.table::set(temp.data, i = as.integer(1), j = match("columnLength", header.data), value = Biostrings::width(col.align)[1])
      gap.count = countAlignmentGaps(non.align)
      data.table::set(temp.data, i = as.integer(1), j = match("columnBasepairs", header.data), value = gap.count[2] - gap.count[1])
      data.table::set(temp.data, i = as.integer(1), j = match("columnGaps", header.data), value = gap.count[1])
      data.table::set(temp.data, i = as.integer(1), j = match("columnPerGaps", header.data), value = gap.count[3])
    }#end trim column.

    if (alignment.assess == TRUE) {
      #Assesses the alignment returning TRUE for pass and FALSE for fail
      test.result = alignmentAssess(alignment = non.align,
                                    min.gap.percent = min.gap.percent,
                                    min.taxa.count = min.taxa.count,
                                    min.align.length = min.align.length)

      data.table::set(temp.data, i = as.integer(1), j = match("Pass", header.data), value = test.result)

      if (test.result == FALSE){
        print(paste0(align.files[i], " Failed and was discarded."))
      } else {
        write.temp = strsplit(as.character(non.align), "")
        aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
        #readies for saving
        writePhylip(aligned.set, file= paste0(output.dir, "/", save.name, ".phy"), interleave = F)
      }#end else test result
    } else {
      #If no alignment assessing is done, saves
      write.temp = strsplit(as.character(non.align), "")
      aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
      #readies for saving
      writePhylip(aligned.set, file= paste0(output.dir, "/", save.name, ".phy"), interleave = F)
    }#end else

    data.frame(temp.data)
    print(paste0(align.files[i], " Completed."))

  }#end i loop

  #parallel::stopCluster(cl)

  #Print and save summary table
  write.csv(out.data, file = paste0(output.dir, "_trimming-summary.csv"), row.names = F)

  #Saves log file of things
  if (file.exists(paste0(output.dir, ".log")) == TRUE){ system(paste0("rm ", output.dir, ".log")) }
  fileConn = file(paste0(output.dir, ".log"), open = "w")
  writeLines(paste0("Log file for ", output.dir), fileConn)
  writeLines(paste0("\n"), fileConn)
  writeLines(paste0("Overall trimming summary:"), fileConn)
  writeLines(paste0("------------------------------------------------------------------"), fileConn)
  writeLines(paste0(""), fileConn)
  writeLines(paste0("Starting alignments: ", length(align.files)), fileConn)
  writeLines(paste0("Trimmed alignments: ", length(out.data$Pass[out.data$Pass == TRUE])), fileConn)
  writeLines(paste0("Discarded alignments: ", length(out.data$Pass[out.data$Pass == FALSE])), fileConn)
  writeLines(paste0(""), fileConn)
  writeLines(paste0("Mean samples removed per alignment: ",
                    mean(out.data$startSamples - out.data$columnSamples)), fileConn)
  writeLines(paste0("Mean alignment length trimmed per alignment: ",
                    mean(out.data$startLength - out.data$columnLength)), fileConn)
  writeLines(paste0("Mean basepairs trimmed per alignment: ",
                    mean(out.data$startBasepairs - out.data$columnBasepairs)), fileConn)
  writeLines(paste0("Mean gaps trimmed per alignment: ",
                    mean(out.data$startGaps - out.data$columnGaps)), fileConn)
  writeLines(paste0("Mean gap percent trimmed per alignment: ",
                    mean(out.data$startPerGaps - out.data$columnPerGaps)), fileConn)
  writeLines(paste0(""), fileConn)
  writeLines(paste0(""), fileConn)
  writeLines(paste0("Individual trimming step summary:"), fileConn)
  writeLines(paste0("------------------------------------------------------------------"), fileConn)
  writeLines(paste0(""), fileConn)
  writeLines(paste0("Starting alignments:"), fileConn)
  writeLines(paste0("Mean samples: ",
                    mean(out.data$startSamples)), fileConn)
  writeLines(paste0("Mean alignment length: ",
                    mean(out.data$startLength)), fileConn)
  writeLines(paste0("Mean basepairs: ",
                    mean(out.data$startBasepairs)), fileConn)
  writeLines(paste0("Mean gaps: ",
                    mean(out.data$startGaps)), fileConn)
  writeLines(paste0("Mean gap percentage: ",
                    mean(out.data$startPerGaps)), fileConn)
  writeLines(paste0(""), fileConn)
  writeLines(paste0("hmmCleaner:"), fileConn)
  writeLines(paste0("Mean samples removed: ",
                    mean(out.data$startSamples - out.data$hmmSamples)), fileConn)
  writeLines(paste0("Mean alignment length reduction: ",
                    mean(out.data$startLength - out.data$hmmLength)), fileConn)
  writeLines(paste0("Mean basepairs trimmed: ",
                    mean(out.data$startBasepairs - out.data$hmmBasepairs)), fileConn)
  writeLines(paste0("Mean gap change: ",
                    mean(out.data$startGaps - out.data$hmmGaps)), fileConn)
  writeLines(paste0("Mean gap percent change: ",
                    mean(out.data$startPerGaps - out.data$hmmPerGaps)), fileConn)
  writeLines(paste0(""), fileConn)
  writeLines(paste0("Trimal:"), fileConn)
  writeLines(paste0("Mean samples removed: ",
                    mean(out.data$hmmSamples - out.data$trimalSamples)), fileConn)
  writeLines(paste0("Mean alignment length reduction: ",
                    mean(out.data$hmmLength - out.data$trimalLength)), fileConn)
  writeLines(paste0("Mean basepairs trimmed: ",
                    mean(out.data$hmmBasepairs - out.data$trimalBasepairs)), fileConn)
  writeLines(paste0("Mean gap change: ",
                    mean(out.data$hmmGaps - out.data$trimalGaps)), fileConn)
  writeLines(paste0("Mean gap percent change: ",
                    mean(out.data$hmmPerGaps - out.data$trimalPerGaps)), fileConn)
  writeLines(paste0(""), fileConn)
  writeLines(paste0("External Trimming:"), fileConn)
  writeLines(paste0("Mean samples removed: ",
                    mean(out.data$trimalSamples - out.data$edgeSamples)), fileConn)
  writeLines(paste0("Mean alignment length reduction: ",
                    mean(out.data$trimalLength - out.data$edgeLength)), fileConn)
  writeLines(paste0("Mean basepairs trimmed: ",
                    mean(out.data$trimalBasepairs - out.data$edgeBasepairs)), fileConn)
  writeLines(paste0("Mean gap change: ",
                    mean(out.data$trimalGaps - out.data$edgeGaps)), fileConn)
  writeLines(paste0("Mean gap percent change: ",
                    mean(out.data$trimalPerGaps - out.data$edgePerGaps)), fileConn)
  writeLines(paste0(""), fileConn)
  writeLines(paste0("Sample Coverage Trimming:"), fileConn)
  writeLines(paste0("Mean samples removed: ",
                    mean(out.data$edgeSamples - out.data$covSamples)), fileConn)
  writeLines(paste0("Mean alignment length reduction: ",
                    mean(out.data$edgeLength - out.data$covLength)), fileConn)
  writeLines(paste0("Mean basepairs trimmed: ",
                    mean(out.data$edgeBasepairs - out.data$covBasepairs)), fileConn)
  writeLines(paste0("Mean gap change: ",
                    mean(out.data$edgeGaps - out.data$covGaps)), fileConn)
  writeLines(paste0("Mean gap percent change: ",
                    mean(out.data$edgePerGaps - out.data$covPerGaps)), fileConn)
  writeLines(paste0(""), fileConn)
  writeLines(paste0("Column Coverage Trimming:"), fileConn)
  writeLines(paste0("Mean samples removed: ",
                    mean(out.data$covSamples - out.data$columnSamples)), fileConn)
  writeLines(paste0("Mean alignment length reduction: ",
                    mean(out.data$covLength - out.data$columnLength)), fileConn)
  writeLines(paste0("Mean basepairs trimmed: ",
                    mean(out.data$covBasepairs - out.data$columnBasepairs)), fileConn)
  writeLines(paste0("Mean gap change: ",
                    mean(out.data$covGaps - out.data$columnGaps)), fileConn)
  writeLines(paste0("Mean gap percent change: ",
                    mean(out.data$covPerGaps - out.data$columnPerGaps)), fileConn)
  close(fileConn)

} #end function
