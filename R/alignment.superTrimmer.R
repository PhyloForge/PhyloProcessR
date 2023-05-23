#' @title superTrimmer
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
#' @param min.alignment.length give a save name if you wnat to save the summary to file.
#'
#' @param min.taxa.alignment TRUE to supress mafft screen output
#'
#' @param min.gap.percent if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param threads path to a folder of sequence alignments in phylip format.
#'
#' @param memory give a save name if you wnat to save the summary to file.
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

superTrimmer = function(alignment.dir = NULL,
                        alignment.format = "phylip",
                        output.dir = NULL,
                        output.format = "phylip",
                        TAPER = FALSE,
                        TAPER.path = NULL,
                        julia.path = NULL,
                        TrimAl = FALSE,
                        TrimAl.path = NULL,
                        trim.external = TRUE,
                        min.external.percent = 50,
                        trim.coverage = TRUE,
                        min.coverage.percent = 50,
                        trim.column = TRUE,
                        min.column.gap.percent = 100,
                        convert.ambiguous.sites = FALSE,
                        alignment.assess = TRUE,
                        min.coverage.bp = 0,
                        min.alignment.length = 0,
                        min.taxa.alignment = 0,
                        max.alignment.gap.percent = 0,
                        threads = 1,
                        memory = 1,
                        overwrite = FALSE,
                        resume = TRUE) {

  # #devtools::install_github("chutter/PHYLOCAP", upgrade = "never")
  # library(PhyloCap)
  # library(foreach)
  #
  #
  # #Save directory
  # #work.dir = "/Volumes/Rodents/Murinae/Trimming"
  # #align.dir = "/Volumes/Rodents/Murinae/Trimming/genes-untrimmed"
  # #feat.gene.names = "/Volumes/Rodents/Murinae/Selected_Transcripts/Mus_gene_metadata.csv"
  # #work.dir = "/home/c111h652/scratch/Rodents/Trimming"
  #
  # work.dir = "/home/c111h652/scratch/Rodents/Trimming"
  # align.dir = "/home/c111h652/scratch/Rodents/Trimming/Emily/genes_untrimmed"
  #
  # work.dir = "/Volumes/Rodents/Murinae/Trimming"
  # align.dir =  "/Volumes/Rodents/Murinae/Trimming/Emily/genes_untrimmed"
  # out.name = "Emily"
  #
  #
  # setwd(work.dir)
  # library(foreach)
  # alignment.dir = "new_alignments"
  # alignment.format = "phylip"
  # output.dir = "trimmed_alignments"
  # output.format = "phylip"
  # overwrite = FALSE
  # TAPER = FALSE
  # TAPER.path = NULL
  # julia.path = NULL
  # TrimAl = TRUE
  # TrimAl.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin"
  # trim.column = TRUE
  # convert.ambiguous.sites = TRUE
  # alignment.assess = F
  # trim.external = TRUE
  # trim.coverage = TRUE
  # min.coverage.percent = 60
  # min.external.percent = 60
  # min.column.gap.percent = 60
  # min.alignment.length = 200
  # min.taxa.alignment = 4
  # min.coverage.bp = 60
  # max.alignment.gap.percent = 50
  # threads = 8
  # memory = 24

  if (is.null(TAPER.path) == FALSE){
    b.string = unlist(strsplit(TAPER.path, ""))
    if (b.string[length(b.string)] != "/") {
      TAPER.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { TAPER.path = NULL }

  if (is.null(julia.path) == FALSE){
    b.string = unlist(strsplit(julia.path, ""))
    if (b.string[length(b.string)] != "/") {
      julia.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { julia.path = NULL }

  if (is.null(TrimAl.path) == FALSE){
    b.string = unlist(strsplit(TrimAl.path, ""))
    if (b.string[length(b.string)] != "/") {
      TrimAl.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { TrimAl.path = NULL }

  if (alignment.dir == output.dir){ stop("You should not overwrite the original alignments.") }

  if (dir.exists(output.dir) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.dir))
      dir.create(output.dir)
     }
  } else { dir.create(output.dir) }

  #Gathers alignments
  align.files = list.files(alignment.dir)

  if (length(align.files) == 0) { stop("alignment files could not be found.") }

  #Skips files done already if resume = TRUE
  if (overwrite == FALSE){
    done.files = list.files(output.dir)
    align.files = align.files[!gsub("\\..*", "", align.files) %in% gsub("\\..*", "", done.files)]
  }

  if (length(align.files) == 0) { return("All alignments have already been completed and overwrite = FALSE.") }

  #Data to collect
  header.data = c("Alignment", "Pass", "startSamples", "tapirSamples", "trimalSamples",
                  "edgeSamples", "columnSamples", "covSamples",
                  "startLength", "tapirLength", "trimalLength",
                  "edgeLength", "columnLength", "covLength",
                  "startBasepairs", "tapirBasepairs", "trimalBasepairs",
                  "edgeBasepairs", "columnBasepairs", "covBasepairs",
                  "startGaps", "tapirGaps", "trimalGaps",
                  "edgeGaps", "columnGaps", "covGaps",
                  "startPerGaps", "tapirPerGaps", "trimalPerGaps",
                  "edgePerGaps", "columnPerGaps", "covPerGaps")

  save.data = data.table::data.table(matrix(as.double(0), nrow = length(align.files), ncol = length(header.data)))
  data.table::setnames(save.data, header.data)
  save.data[, Alignment:=as.character(Alignment)]
  save.data[, Pass:=as.logical(Pass)]

  #Sets up multiprocessing
  cl = parallel::makeCluster(threads, outfile = "")
  doParallel::registerDoParallel(cl)
  mem.cl = floor(memory/threads)

  #Loops through each locus and does operations on them
  out.data = foreach::foreach(i=1:length(align.files), .combine = rbind, .packages = c("PhyloProcessR", "foreach", "Biostrings","Rsamtools", "ape", "stringr", "data.table")) %dopar% {
  #for (i in 1:length(align.files)){
    print(paste0(align.files[i], " Starting..."))

     #Load in alignments
    if (alignment.format == "phylip"){
      align = Biostrings::readAAMultipleAlignment(file = paste0(alignment.dir, "/", align.files[i]), format = "phylip")

    #  align = Biostrings::readDNAStringSet(file = paste0(alignment.dir, "/", align.files[i]), format = "phylip")
    #  align = readLines(paste0(alignment.dir, "/", align.files[i]))[-1]
    #  align = gsub(".*\\ ", "", align)
    #  char.count = nchar(align)

      align = Biostrings::DNAStringSet(align)
      save.name = gsub(".phy$", "", align.files[i])
      save.name = gsub(".phylip$", "", save.name)
    }#end phylip

    if (alignment.format == "fasta"){
      align = Biostrings::readDNAStringSet(paste0(alignment.dir, "/", align.files[i]) )
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

    #Convert ambiguous sites
    if (convert.ambiguous.sites == TRUE){
      non.align = convertAmbiguousConsensus(alignment = non.align)
    }#end phylip


    #Summarize all this, no functoin
    data.table::set(temp.data, i = as.integer(1), j = match("Alignment", header.data), value = save.name)
    data.table::set(temp.data, i = as.integer(1), j = match("startSamples", header.data), value = length(non.align))
    data.table::set(temp.data, i = as.integer(1), j = match("startLength", header.data), value = Biostrings::width(non.align)[1] )
    gap.count = countAlignmentGaps(non.align)
    data.table::set(temp.data, i = as.integer(1), j = match("startBasepairs", header.data), value = gap.count[2] - gap.count[1])
    data.table::set(temp.data, i = as.integer(1), j = match("startGaps", header.data), value = gap.count[1])
    data.table::set(temp.data, i = as.integer(1), j = match("startPerGaps", header.data), value = gap.count[3])

    # Run the TAPER
    if (TAPER == TRUE && length(non.align) != 0){
      #Runs the TAPER
      taper.align = trimTAPER(alignment = non.align,
                              TAPER.path = TAPER.path,
                              julia.path = julia.path,
                              quiet = F,
                              delete.temp = T)
      non.align = taper.align
      #Saves the data
      data.table::set(temp.data, i = as.integer(1), j = match("tapirSamples", header.data), value = length(taper.align))
      data.table::set(temp.data, i = as.integer(1), j = match("tapirLength", header.data), value = Biostrings::width(taper.align)[1])
      gap.count = countAlignmentGaps(non.align)
      data.table::set(temp.data, i = as.integer(1), j = match("tapirBasepairs", header.data), value = gap.count[2] - gap.count[1])
      data.table::set(temp.data, i = as.integer(1), j = match("tapirGaps", header.data), value = gap.count[1])
      data.table::set(temp.data, i = as.integer(1), j = match("tapirPerGaps", header.data), value = gap.count[3])
    }#end if

    #Step 3. Trimal trimming
    if (TrimAl == TRUE && length(non.align) != 0){

      trimal.align = trimTrimal(alignment = non.align,
                                trimal.path = TrimAl.path,
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
    if (trim.external == TRUE && length(non.align) != 0){
      #external trimming function
      edge.align = trimExternal(alignment = non.align,
                                min.n.seq = ceiling(length(non.align) * (min.external.percent/100)),
                                codon.trim = F)

      if (class(edge.align) == "numeric") { edge.align = Biostrings::DNAStringSet() }
      non.align = edge.align

      #Saves stat data
      data.table::set(temp.data, i = as.integer(1), j = match("edgeSamples", header.data), value = length(edge.align))
      data.table::set(temp.data, i = as.integer(1), j = match("edgeLength", header.data), value = Biostrings::width(edge.align)[1])
      gap.count = countAlignmentGaps(non.align)
      data.table::set(temp.data, i = as.integer(1), j = match("edgeBasepairs", header.data), value = gap.count[2] - gap.count[1])
      data.table::set(temp.data, i = as.integer(1), j = match("edgeGaps", header.data), value = gap.count[1])
      data.table::set(temp.data, i = as.integer(1), j = match("edgePerGaps", header.data), value = gap.count[3])
    }#end trim external

    if (trim.column == TRUE && length(non.align) != 0){
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

    #Step 5. Evaluate and cut out each sample
    if (trim.coverage == TRUE && length(non.align) != 0){
      #sample coverage function
      cov.align = trimSampleCoverage(alignment = non.align,
                                     min.coverage.percent = min.coverage.percent,
                                     min.coverage.bp = min.coverage.bp,
                                     relative.width = "sample")

      non.align = cov.align
      #Saves stat data
      data.table::set(temp.data, i = as.integer(1), j = match("covSamples", header.data), value = length(cov.align))
      data.table::set(temp.data, i = as.integer(1), j = match("covLength", header.data), value = Biostrings::width(cov.align)[1])
      gap.count = countAlignmentGaps(non.align)
      data.table::set(temp.data, i = as.integer(1), j = match("covBasepairs", header.data), value = gap.count[2] - gap.count[1])
      data.table::set(temp.data, i = as.integer(1), j = match("covGaps", header.data), value = gap.count[1])
      data.table::set(temp.data, i = as.integer(1), j = match("covPerGaps", header.data), value = gap.count[3])
    }#end trim.external



    #Step 6
    if (length(non.align) != 0){
      if (alignment.assess == TRUE) {
        #Assesses the alignment returning TRUE for pass and FALSE for fail
        test.result = alignmentAssess(alignment = non.align,
                                      max.alignment.gap.percent = max.alignment.gap.percent,
                                      min.taxa.alignment = min.taxa.alignment,
                                      min.alignment.length = min.alignment.length)

        data.table::set(temp.data, i = as.integer(1), j = match("Pass", header.data), value = test.result)

        if (test.result == FALSE){
          print(paste0(align.files[i], " failed filtering and was discarded."))
        } else {
          print(paste0(align.files[i], " passed filters and was saved to file."))
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
    }#outer if

    print(paste0(align.files[i], " Completed."))
    print(data.frame(temp.data))

    #rm()
    #gc()

  }#end i loop

  parallel::stopCluster(cl)

  #Print and save summary table
  write.csv(out.data, file = "alignment-trimming_summary.csv", row.names = F)
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
