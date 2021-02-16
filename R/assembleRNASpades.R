#' @title assembleRNASpades
#'
#' @description Function for running the program spades to assemble short read sequencing data
#'
#' @param read.directory directory of processed reads
#'
#' @param output.directory save name for the output directory
#'
#' @param full.path.spades contigs are added into existing alignment if algorithm is "add"
#'
#' @param kmer.values if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param threads number of computation processing threads
#'
#' @param mem amount of system memory to use
#'
#' @param resume TRUE to skip samples already completed
#'
#' @param overwrite TRUE to overwrite a folder of samples with output.dir
#'
#' @param quiet TRUE to supress screen output

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

assembleRNASpades = function(input.reads = NULL,
                             output.directory = "processed-reads/rnaspades-assembly",
                             assembly.directory = "draft-transcripts",
                             spades.path = "spades.py",
                             kmer.values = c(21,33,55,77,99,127),
                             threads = 1,
                             memory = 4,
                             overwrite = FALSE,
                             resume = TRUE,
                             save.corrected.reads = FALSE,
                             quiet = TRUE) {

  # #debug
  # setwd("/Users/chutter/Dropbox/Research/0_Github/Test-dataset")
  # input.reads = "/Users/chutter/Dropbox/Research/0_Github/Test-dataset/processed-reads/pe-merged-reads"
  # spades.path = "/Users/chutter/miniconda3/bin/rnaspades.py"
  # output.directory = "processed-reads/rnaspades-assembly"
  # assembly.directory = "draft-transcripts"
  # kmer.values = c(21,33,55,77,99,127)
  # threads = 1
  # memory = 4
  # overwrite = FALSE
  # quiet = FALSE
  # resume = TRUE

  #Quick checks
  options(stringsAsFactors = FALSE)
  if (is.null(input.reads) == TRUE){ stop("Please provide input reads.") }
  if (file.exists(input.reads) == F){ stop("Input reads not found.") }
  if (is.null(output.directory) == TRUE){ stop("Please provide an output directory.") }
  if (is.null(assembly.directory) == TRUE){ stop("Please provide an contig save directory.") }

  #Sets directory and reads
  if (dir.exists(output.directory) == F){ dir.create(output.directory) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  }#end else

  #Sets directory and reads
  if (dir.exists(assembly.directory) == F){ dir.create(assembly.directory) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", assembly.directory))
      dir.create(assembly.directory)
    }
  }#end else

  #Creates output directory
  if (dir.exists("logs") == F){ dir.create("logs") }

  #Sets up the reads
  files = list.files(path = input.reads, full.names = T, recursive = T)
  reads = files[grep(pattern = "fastq|fq|clustS", x = files)]

  samples = gsub(paste0(input.reads, "/"), "", reads)
  samples = unique(gsub("/.*", "", samples))

  #Skips samples already finished
  if (resume == TRUE){
    done.names = list.files(assembly.directory)
    samples = samples[!samples %in% gsub(".fa$", "", done.names)]
  } else { samples = samples }

  if (length(samples) == 0){ stop("No samples to run or incorrect directory.") }

  #Header data for features and whatnot
  for (i in 1:length(samples)){

    sample.reads = reads[grep(samples[i], reads)]

    #Run SPADES on sample
    k.val = paste(kmer.values, collapse = ",")

    #Creates assembly reads folder if not present
    save.assem = paste0(output.directory, "/", samples[i])
    dir.create(save.assem)

    #Sorts reads
    sample.names = unique(gsub("_READ.*", "", sample.reads))

    #Creates a spades character string to run different reads and library configurations
    final.read.string = c()
    for (j in 1:length(sample.names)) {
      #Gets the sample reads
      lib.reads = sample.reads[grep(paste0(sample.names[j], "_"), sample.reads)]
      #Checks for different read lengths
      if (length(lib.reads) == 1){ read.string = paste0("--s1 ", lib.reads[1], " ") }
      if (length(lib.reads) >= 2){ read.string = paste0("--pe", j,"-1 ", lib.reads[1],
                                                        " --pe", j, "-2 ", lib.reads[2], " ") }
      if (length(lib.reads) == 3){ read.string = paste0(read.string, "--pe", j , "-m ", lib.reads[3], " ") }
      final.read.string = paste0(final.read.string, read.string)
    }#end j loop

    #Runs spades command
    system(paste0(spades.path, " --rna ", final.read.string,
                  "-o ", save.assem, " -k ", k.val, " ",
                  "-t ", threads, " -m ", memory),
           ignore.stdout = quiet)

    #Crashes function if spades failed, also copies new assemblies to assembly.directory
    if (file.exists(paste0(save.assem, "/transcripts.fasta")) == TRUE ){
      system(paste0("cp ", save.assem, "/transcripts.fasta ",  assembly.directory,
                    "/", samples[i], ".fa"))
    } else { stop(paste0("spades error for ", samples[i], ", check spades.log file in rnaspades-assembly folder.")) }

    if (save.corrected.reads == FALSE) {
      #copies to contigs folder
      system(paste0("rm -r ", save.assem, "/corrected"))
    } # end if

    system(paste0("rm -r ", save.assem, "/tmp"))
    print(paste0(samples[i], " Completed rnaSPAdes asssembly!"))

  }#end sample loop

}#end function

