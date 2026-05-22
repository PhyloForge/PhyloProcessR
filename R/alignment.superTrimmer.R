#' @title superTrimmer
#'
#' @description Applies a configurable pipeline of trimming steps to every alignment in a directory and saves the trimmed alignments to an output directory. Steps are applied in this order: (1) replace N with gap; (2) optionally convert ambiguous sites; (3) optionally run TrimAl (automated mode); (4) optionally trim external (edge) columns; (5) optionally trim high-gap columns; (6) optionally remove low-coverage samples. After trimming, alignments can optionally be assessed against minimum length and minimum taxa thresholds before saving. A per-alignment summary CSV and a log file are written to disk. Runs in parallel using foreach/doParallel.
#'
#' @param alignment.dir path to the directory of input alignment files
#'
#' @param alignment.format format of the input alignments; "phylip" or "fasta"
#'
#' @param output.dir path to the directory where trimmed alignments will be saved
#'
#' @param TrimAl if TRUE, run TrimAl (automated1 mode) as a trimming step
#'
#' @param TrimAl.path system path to the directory containing the trimal executable; NULL to use the system PATH
#'
#' @param trim.external if TRUE, trim poorly covered columns from the external (5' and 3') ends of the alignment
#'
#' @param min.external.percent minimum percentage of sequences that must have data at a column for it to be retained during external trimming
#'
#' @param trim.coverage if TRUE, remove samples whose sequence coverage falls below the minimum thresholds
#'
#' @param min.coverage.percent minimum percentage of the alignment (or the longest sample) that each sample must cover to be retained
#'
#' @param trim.column if TRUE, remove alignment columns where the gap percentage meets or exceeds min.column.gap.percent
#'
#' @param min.column.gap.percent gap percentage threshold (0-100) at or above which a column is removed
#'
#' @param convert.ambiguous.sites if TRUE, convert ambiguous (IUPAC) bases to gaps before trimming
#'
#' @param alignment.assess if TRUE, apply final minimum length and minimum taxa filters and record pass/fail status in the summary table
#'
#' @param min.coverage.bp minimum number of non-gap base pairs a sample must have to be retained during sample coverage trimming
#'
#' @param min.alignment.length minimum alignment length in bp required to save an alignment
#'
#' @param min.taxa.alignment minimum number of sequences required to save an alignment
#'
#' @param max.alignment.gap.percent maximum overall gap percentage allowed for an alignment to pass the assessment step
#'
#' @param threads number of parallel processing threads
#'
#' @param memory total memory in GB to allocate across threads
#'
#' @param overwrite if TRUE, delete and recreate the output directory; if FALSE, keep existing trimmed files and skip them
#'
#' @return saves trimmed alignment files to output.dir, writes a summary CSV (alignment-trimming_summary.csv) and a log file; nothing is returned to R
#'
#' @export

superTrimmer = function(alignment.dir = NULL,
                        alignment.format = "phylip",
                        output.dir = NULL,
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
                        overwrite = FALSE) {


  # TrimAl.path = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  # alignment.dir = "/Users/chutter/Downloads/untrimmed_all-markers_new"
  # alignment.format = "phylip"
  # output.dir = "/Users/chutter/Downloads/trimmed_all-markers_new"
  # overwrite = FALSE
  # TrimAl = TRUE
  # trim.column = TRUE
  # convert.ambiguous.sites = FALSE
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
  # threads = 6
  # memory = 24
#
  # alignment.dir = "data-analysis/alignments/untrimmed_only-flanks-unlinked"
  # alignment.format = "phylip"
  # output.dir = "data-analysis/alignments/trimmed_only-flanks-unlinked"
  # overwrite = overwrite
  # TrimAl = run.TrimAl
  # TrimAl.path = trimAl.path
  # trim.column = trim.column
  # convert.ambiguous.sites = convert.ambiguous.sites
  # alignment.assess = FALSE
  # trim.external = trim.external
  # trim.coverage = trim.coverage
  # min.coverage.percent = min.coverage.percent
  # min.external.percent = min.external.percent
  # min.column.gap.percent = min.column.gap.percent
  # min.alignment.length = min.alignment.length
  # min.taxa.alignment = min.taxa.alignment
  # min.coverage.bp = min.coverage.bp
  # threads = threads
  # memory = memory

  if (is.null(TrimAl.path) == FALSE){
    b.string = unlist(strsplit(TrimAl.path, ""))
    if (b.string[length(b.string)] != "/") {
      TrimAl.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { TrimAl.path = "" }

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
  header.data = c("Alignment", "Pass", "startSamples", "trimalSamples",
                  "edgeSamples", "columnSamples", "covSamples",
                  "startLength", "trimalLength",
                  "edgeLength", "columnLength", "covLength",
                  "startBasepairs", "trimalBasepairs",
                  "edgeBasepairs", "columnBasepairs", "covBasepairs",
                  "startGaps",  "trimalGaps",
                  "edgeGaps", "columnGaps", "covGaps",
                  "startPerGaps", "trimalPerGaps",
                  "edgePerGaps", "columnPerGaps", "covPerGaps")

  save.data = data.table::data.table(matrix(as.double(0), nrow = 1L, ncol = length(header.data)))
  data.table::setnames(save.data, header.data)
  save.data[, Alignment:=as.character(Alignment)]
  save.data[, Pass:=as.logical(Pass)]

  mem.cl = floor(memory/threads)

  #Loops through each locus and does operations on them
  out.data = do.call(rbind, parallel::mclapply(seq_along(align.files), function(i) {
  tryCatch({
    print(paste0(align.files[i], " Starting..."))

     #Load in alignments
    if (alignment.format == "phylip"){
      align = Biostrings::DNAStringSet(Biostrings::readDNAMultipleAlignment(
        file = paste0(alignment.dir, "/", align.files[i]), format = "phylip"
      ))
      save.name = gsub(".phy$", "", align.files[i])
      save.name = gsub(".phylip$", "", save.name)
    }#end phylip

    if (alignment.format == "fasta"){
      align = Biostrings::readDNAStringSet(paste0(alignment.dir, "/", align.files[i]) )
      save.name = gsub(".fa$", "", align.files[i])
      save.name = gsub(".fasta$", "", save.name)
    }#end phylip

    if (length(align) <= min.taxa.alignment) {
      print("too few taxa in alignment, skipping")
      #next
      return(NULL)
    }

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

      # Remove gap-only columns created by sample removal
      if (length(non.align) != 0) {
        non.align = trimAlignmentColumns(alignment = non.align, min.gap.percent = 100)
      }
    }#end trim.external

    if (length(non.align) <= min.taxa.alignment) {
      print("too few taxa in alignment, skipping")
      #next
      return(NULL)
    }

    if (length(non.align) == 0 || unique(Biostrings::width(non.align)) < min.alignment.length) {
      print("alignment below minimum length, skipping")
      #next
      return(NULL)
    }

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
          PhyloProcessR::writePhylip(aligned.set, file= paste0(output.dir, "/", save.name, ".phy"), interleave = F)
        }#end else test result
      } else {
        #If no alignment assessing is done, saves
        write.temp = strsplit(as.character(non.align), "")
        aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
        #readies for saving
        PhyloProcessR::writePhylip(aligned.set, file= paste0(output.dir, "/", save.name, ".phy"), interleave = F)
      }#end else
    }#outer if

    print(paste0(align.files[i], " Completed."))
    print(data.frame(temp.data))

    rm(align, non.align)
    gc()
    temp.data  # explicit return value for do.call(rbind, ...)

  }, error = function(e) {
    warning(align.files[i], " failed: ", conditionMessage(e))
    NULL
  })
  }, mc.cores = threads)) #end i loop

  if (is.null(out.data) ==  TRUE){ return("No alignments were trimmed.") }

  #Print and save summary table and log to the logs/ directory
  dir.create("logs", recursive = TRUE, showWarnings = FALSE)
  log.base = paste0("logs/", basename(output.dir))
  write.csv(out.data, file = paste0(log.base, "_trimming_summary.csv"), row.names = F)

  #Saves log file of things
  if (file.exists(paste0(log.base, ".log")) == TRUE){ system(paste0("rm ", log.base, ".log")) }
  fileConn = file(paste0(log.base, ".log"), open = "w")
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
                    mean(out.data$startSamples - out.data$covSamples)), fileConn)
  writeLines(paste0("Mean alignment length trimmed per alignment: ",
                    mean(out.data$startLength - out.data$covLength)), fileConn)
  writeLines(paste0("Mean basepairs trimmed per alignment: ",
                    mean(out.data$startBasepairs - out.data$covBasepairs)), fileConn)
  writeLines(paste0("Mean gaps trimmed per alignment: ",
                    mean(out.data$startGaps - out.data$covGaps)), fileConn)
  writeLines(paste0("Mean gap percent trimmed per alignment: ",
                    mean(out.data$startPerGaps - out.data$covPerGaps)), fileConn)
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
  writeLines(paste0("Trimal:"), fileConn)
  writeLines(paste0("Mean samples removed: ",
                    mean(out.data$startSamples - out.data$trimalSamples)), fileConn)
  writeLines(paste0("Mean alignment length reduction: ",
                    mean(out.data$startLength - out.data$trimalLength)), fileConn)
  writeLines(paste0("Mean basepairs trimmed: ",
                    mean(out.data$startBasepairs - out.data$trimalBasepairs)), fileConn)
  writeLines(paste0("Mean gap change: ",
                    mean(out.data$startGaps - out.data$trimalGaps)), fileConn)
  writeLines(paste0("Mean gap percent change: ",
                    mean(out.data$startPerGaps - out.data$trimalPerGaps)), fileConn)
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
  writeLines(paste0("Column Coverage Trimming:"), fileConn)
  writeLines(paste0("Mean samples removed: ",
                    mean(out.data$edgeSamples - out.data$columnSamples)), fileConn)
  writeLines(paste0("Mean alignment length reduction: ",
                    mean(out.data$edgeLength - out.data$columnLength)), fileConn)
  writeLines(paste0("Mean basepairs trimmed: ",
                    mean(out.data$edgeBasepairs - out.data$columnBasepairs)), fileConn)
  writeLines(paste0("Mean gap change: ",
                    mean(out.data$edgeGaps - out.data$columnGaps)), fileConn)
  writeLines(paste0("Mean gap percent change: ",
                    mean(out.data$edgePerGaps - out.data$columnPerGaps)), fileConn)
  writeLines(paste0(""), fileConn)
  writeLines(paste0("Sample Coverage Trimming:"), fileConn)
  writeLines(paste0("Mean samples removed: ",
                    mean(out.data$columnSamples - out.data$covSamples)), fileConn)
  writeLines(paste0("Mean alignment length reduction: ",
                    mean(out.data$columnLength - out.data$covLength)), fileConn)
  writeLines(paste0("Mean basepairs trimmed: ",
                    mean(out.data$columnBasepairs - out.data$covBasepairs)), fileConn)
  writeLines(paste0("Mean gap change: ",
                    mean(out.data$columnGaps - out.data$covGaps)), fileConn)
  writeLines(paste0("Mean gap percent change: ",
                    mean(out.data$columnPerGaps - out.data$covPerGaps)), fileConn)
  close(fileConn)

} #end function
