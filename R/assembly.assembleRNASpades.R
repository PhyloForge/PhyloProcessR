#' @title assembleRNASpades
#'
#' @description Runs rnaSPAdes (\code{spades.py --rna}) on a directory of
#'   processed reads to produce de novo transcriptome assemblies. Each sample
#'   subdirectory under \code{input.reads} is assembled independently with
#'   support for multi-lane and paired/single-end read configurations. The
#'   resulting \code{transcripts.fasta} file from each sample is copied to
#'   \code{assembly.directory} as \code{<sample>.fa}. Samples for which a
#'   \code{.fa} file already exists in \code{assembly.directory} are skipped
#'   unless \code{overwrite = TRUE}.
#'
#' @param input.reads path to a directory of processed reads. Each sample must
#'   occupy its own subdirectory containing FASTQ files whose names encode read
#'   number (e.g. \code{_READ1}, \code{_READ2}).
#'
#' @param output.directory path to the directory where per-sample rnaSPAdes
#'   working directories will be written. Default:
#'   \code{"processed-reads/rnaspades-assembly"}.
#'
#' @param assembly.directory path to the directory where the final transcript
#'   FASTA files (.fa) are copied after assembly. Default:
#'   \code{"draft-transcripts"}.
#'
#' @param spades.path path to the directory containing \code{spades.py}. If
#'   \code{NULL} expected on the system PATH. Default: \code{NULL}.
#'
#' @param kmer.values integer vector of k-mer sizes passed to rnaSPAdes with
#'   \code{-k}. Default: \code{c(21, 33, 55, 77, 99, 127)}.
#'
#' @param threads number of CPU threads passed to rnaSPAdes with \code{-t}.
#'   Default: \code{1}.
#'
#' @param memory RAM in GB passed to rnaSPAdes with \code{-m}. Default:
#'   \code{4}.
#'
#' @param overwrite logical; if \code{TRUE} existing output and assembly
#'   directories are deleted and recreated and all samples are rerun. Default:
#'   \code{FALSE}.
#'
#' @param save.corrected.reads logical; if \code{FALSE} (default) the
#'   \code{corrected/} subdirectory produced by rnaSPAdes is deleted after
#'   assembly to save disk space. Default: \code{FALSE}.
#'
#' @param quiet logical; if \code{TRUE} rnaSPAdes screen output is suppressed.
#'   Default: \code{TRUE}.
#'
#' @return Invisibly returns nothing. Assembled transcripts for each sample are
#'   saved as \code{<assembly.directory>/<sample>.fa}.
#'
#' @export

assembleRNASpades = function(input.reads = NULL,
                             output.directory = "processed-reads/rnaspades-assembly",
                             assembly.directory = "draft-transcripts",
                             spades.path = NULL,
                             kmer.values = c(21,33,55,77,99,127),
                             threads = 1,
                             memory = 4,
                             overwrite = FALSE,
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


  #Same adds to bbmap path
  if (is.null(spades.path) == FALSE){
    b.string = unlist(strsplit(spades.path, ""))
    if (b.string[length(b.string)] != "/") {
      spades.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { spades.path = "" }

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
  if (dir.exists("logs/sample_logs") == F){ dir.create("logs/sample_logs", recursive = TRUE) }

  #Sets up the reads
  files = list.files(path = input.reads, full.names = T, recursive = T)
  reads = files[grep(pattern = "fastq|fq|clustS", x = files)]

  samples = gsub(paste0(input.reads, "/"), "", reads)
  samples = unique(gsub("/.*", "", samples))

  #Skips samples already finished
  if (overwrite == FALSE){
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
    system(paste0(spades.path, "spades.py --rna ", final.read.string,
                  "-o ", save.assem, " -k ", k.val, " ",
                  "-t ", threads, " -m ", memory),
           ignore.stdout = quiet)

    #Crashes function if spades failed, also copies new assemblies to assembly.directory
    if (file.exists(paste0(save.assem, "/transcripts.fasta")) == TRUE ){
      system(paste0("cp ", save.assem, "/transcripts.fasta ",  assembly.directory,
                    "/", samples[i], ".fa"))
    } else { print(paste0("spades error for ", samples[i], ", check spades.log file in rnaspades-assembly folder.")) }

    if (save.corrected.reads == FALSE) {
      #copies to contigs folder
      system(paste0("rm -r ", save.assem, "/corrected"))
    } # end if

    system(paste0("rm -r ", save.assem, "/tmp"))
    print(paste0(samples[i], " Completed rnaSPAdes asssembly!"))

  }#end sample loop

}#end function

