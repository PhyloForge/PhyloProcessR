#' @title assembleSpades
#'
#' @description Runs SPAdes (\code{spades.py}) on a directory of processed reads
#'   to produce de novo genome assemblies. Each sample subdirectory under
#'   \code{input.reads} is assembled independently with support for multi-lane
#'   and mixed paired/single-end/merged read configurations. The
#'   \code{scaffolds.fasta} output from each sample is copied to
#'   \code{assembly.directory} as \code{<sample>.fa}. Samples for which a
#'   \code{.fa} file already exists in \code{assembly.directory} are skipped
#'   unless \code{overwrite = TRUE}. The \code{--careful} and \code{--isolate}
#'   SPAdes modes cannot be used simultaneously.
#'
#' @param input.reads path to a directory of processed reads. Each sample must
#'   occupy its own subdirectory containing FASTQ files whose names encode read
#'   number and pair identity.
#'
#' @param output.directory path to the directory where per-sample SPAdes working
#'   directories will be written. Default:
#'   \code{"processed-reads/spades-assembly"}.
#'
#' @param assembly.directory path to the directory where the final scaffold
#'   FASTA files (.fa) are copied after assembly. Default:
#'   \code{"draft-assemblies"}.
#'
#' @param spades.path path to the directory containing \code{spades.py}. If
#'   \code{NULL} expected on the system PATH. Default: \code{NULL}.
#'
#' @param mismatch.corrector logical; if \code{TRUE} passes \code{--careful} to
#'   SPAdes to reduce mismatches and indels in the assembly. Cannot be \code{TRUE}
#'   when \code{isolate = TRUE}. Default: \code{TRUE}.
#'
#' @param isolate logical; if \code{TRUE} passes \code{--isolate} to SPAdes,
#'   recommended for highly covered isolate genomes. Cannot be \code{TRUE} when
#'   \code{mismatch.corrector = TRUE}. Default: \code{FALSE}.
#'
#' @param kmer.values integer vector of k-mer sizes passed to SPAdes with
#'   \code{-k}. Default: \code{c(33, 55, 77, 99, 127)}.
#'
#' @param threads number of CPU threads passed to SPAdes with \code{-t}.
#'   Default: \code{1}.
#'
#' @param memory RAM in GB passed to SPAdes with \code{-m}. Default: \code{4}.
#'
#' @param overwrite logical; if \code{TRUE} existing output and assembly
#'   directories are deleted and recreated and all samples are rerun. Default:
#'   \code{FALSE}.
#'
#' @param save.corrected.reads logical; if \code{FALSE} (default) the
#'   \code{corrected/} subdirectory produced by SPAdes is deleted after assembly
#'   to save disk space. Default: \code{FALSE}.
#'
#' @param quiet logical; if \code{TRUE} SPAdes screen output is suppressed.
#'   Default: \code{TRUE}.
#'
#' @return Invisibly returns nothing. Assembled scaffolds for each sample are
#'   saved as \code{<assembly.directory>/<sample>.fa}.
#'
#' @export

assembleSpades = function(input.reads = NULL,
                          output.directory = "processed-reads/spades-assembly",
                          assembly.directory = "draft-assemblies",
                          spades.path = NULL,
                          mismatch.corrector = TRUE,
                          isolate = FALSE,
                          kmer.values = c(33,55,77,99,127),
                          threads = 1,
                          memory = 4,
                          overwrite = FALSE,
                          save.corrected.reads = FALSE,
                          quiet = TRUE) {

  # #debug
  # setwd("/Users/chutter/Dropbox/Research/0_Github/Test-dataset")
  # input.reads = "/Users/chutter/Dropbox/Research/0_Github/Test-dataset/processed-reads/pe-merged-reads"
  # spades.path = "/Users/chutter/miniconda3/bin/spades.py"
  # output.directory = "processed-reads/spades-assembly"
  # assembly.directory = "draft-assemblies"
  # mismatch.corrector = F
  # kmer.values = c(21,33,55,77,99,127)
  # threads = 1
  # memory = 4
  # overwrite = FALSE
  # quiet = TRUE
  # resume = TRUE

  # Same adds to bbmap path
  if (is.null(spades.path) == FALSE) {
    b.string = unlist(strsplit(spades.path, ""))
    if (b.string[length(b.string)] != "/") {
      spades.path = paste0(append(b.string, "/"), collapse = "")
    } # end if
  } else {
    spades.path = ""
  }

  # Quick checks
  if (is.null(input.reads) == TRUE) {
    stop("Please provide input reads.")
  }
  if (file.exists(input.reads) == F) {
    stop("Input reads not found.")
  }
  if (is.null(output.directory) == TRUE) {
    stop("Please provide an output directory.")
  }
  if (is.null(assembly.directory) == TRUE) {
    stop("Please provide an contig save directory.")
  }

  # Sets directory and reads
  if (dir.exists(output.directory) == F) {
    dir.create(output.directory)
  } else {
    if (overwrite == TRUE) {
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } # end else

  # Sets directory and reads
  if (dir.exists(assembly.directory) == F) {
    dir.create(assembly.directory)
  } else {
    if (overwrite == TRUE) {
      system(paste0("rm -r ", assembly.directory))
      dir.create(assembly.directory)
    }
  } # end else

  # Creates output directory
  if (dir.exists("logs") == F) {
    dir.create("logs")
  }

  if (isolate == TRUE && mismatch.corrector == TRUE) {
    stop("Both --isolate or --careful (mismatch corrector) can not be used together. Please choose only one.")
  }

  if (mismatch.corrector == FALSE && isolate == FALSE) {
    mismatch.string = ""
  }
  
  if (isolate == TRUE) {
    mismatch.string = "--isolate "
  }
  
  if (mismatch.corrector == TRUE) {
    mismatch.string = "--careful "
  }

  # Sets up the reads
  files <- list.files(path = input.reads, full.names = T, recursive = T)
  reads <- files[grep(pattern = "fastq|fq|clustS", x = files)]

  samples <- gsub(paste0(input.reads, "/"), "", reads)
  samples <- unique(gsub("/.*", "", samples))

  # Skips samples already finished
  if (overwrite == FALSE) {
    done.names <- list.files(assembly.directory)
    samples <- samples[!samples %in% gsub(".fa$", "", done.names)]
  } else {
    samples <- samples
  }

  if (length(samples) == 0) {
    stop("No samples to run or incorrect directory.")
  }
  #Header data for features and whatnot
  for (i in seq_along(samples)){

    sample.reads = reads[grep(pattern = paste0(samples[i], "_"), x = reads)]

    #Checks the Sample column in case already renamed
    if (length(sample.reads) == 0){ sample.reads = reads[grep(pattern = samples[i], x = reads)] }

    #Returns an error if reads are not found
    if (length(sample.reads) == 0 ){
      stop(sample.names[i], " does not have any reads present. Sample folder may be empty. ")
    } #end if statement

    #Run SPADES on sample
    k.val = paste(kmer.values, collapse = ",")

    #Creates assembly reads folder if not present
    save.assem = paste0(output.directory, "/", samples[i])
    dir.create(save.assem)

    #Sorts reads
    sample.lanes = unique(gsub("_1.f.*|_2.f.*|_3.f.*|-1.f.*|-2.f.*|-3.f.*|_R1_.*|_R2_.*|_R3_.*|_READ1_.*|_READ2_.*|_READ3_.*|_R1.f.*|_R2.f.*|_R3.f.*|-R1.f.*|-R2.f.*|-R3.f.*|_READ1.f.*|_READ2.f.*|_READ3.f.*|-READ1.f.*|-READ2.f.*|-READ3.f.*|_singleton.*|-singleton.*|READ-singleton.*|READ_singleton.*|_READ-singleton.*|-READ_singleton.*|-READ-singleton.*|_READ_singleton.*", "", sample.reads))

    #Creates a spades character string to run different reads and library configurations
    final.read.string = c()
    for (j in 1:length(sample.lanes)) {
      #Gets the sample reads
      lib.reads = sample.reads[grep(paste0(sample.lanes[j]), sample.reads)]

      #Concatenate together
      lib.read1 = lib.reads[grep("_1.f.*|-1.f.*|_R1_.*|-R1_.*|_R1-.*|-R1-.*|READ1.*|_R1.fast.*|-R1.fast.*", lib.reads)]
      lib.read2 = lib.reads[grep("_2.f.*|-2.f.*|_R2_.*|-R2_.*|_R2-.*|-R2-.*|READ2.*|_R2.fast.*|-R2.fast.*", lib.reads)]
      lib.read3 = lib.reads[grep("_3.f.*|-3.f.*|_R3_.*|-R3_.*|_R3-.*|-R3-.*|READ3.*|_R3.fast.*|-R3.fast.*|_READ3.fast.*|-READ3.fast.*|_singleton.*|-singleton.*|READ-singleton.*|READ_singleton.*|_READ-singleton.*|-READ_singleton.*|-READ-singleton.*|_READ_singleton.*", lib.reads)]

      #Checks for different read lengths
      if (length(lib.read1) == 1){ read.string = paste0("--s1 ", lib.read1, " ") }
      if (length(lib.read2) == 1){ read.string = paste0("--pe", j,"-1 ", lib.read1,
                                                        " --pe", j, "-2 ", lib.read2, " ") }
      if (length(lib.read3) == 1){ read.string = paste0(read.string, "--pe", j , "-m ", lib.read3, " ") }
      final.read.string = paste0(final.read.string, read.string)
    }#end j loop

#     pe.read1 = sample.reads[grep("_READ1", sample.reads)]
#     if (length(pe.read1) != 0){ pe.read1.string = paste0("--pe", rep(1:length(pe.read1)), "-1 ", pe.read1, collapse = " ") }
#     pe.read2 = sample.reads[grep("_READ2", sample.reads)]
#     if (length(pe.read2) != 0){ pe.read2.string = paste0("--pe", rep(1:length(pe.read2)), "-2 ", pe.read2, collapse = " ") }
#     mg.read3 = sample.reads[grep("_READ3", sample.reads)]
#     if (length(mg.read3) != 0){ mg.read3.string = paste0("--pe", rep(1:length(mg.read3)), "-m ", mg.read3, collapse = " ") }

    #Runs spades command
    system(paste0(spades.path, "spades.py ", final.read.string,
                  "-o ", save.assem, " -k ", k.val, " ", mismatch.string,
                  "-t ", threads, " -m ", memory),
           ignore.stdout = quiet)

    #Crashes function if spades failed, also copies new assemblies to assembly.directory
    if (file.exists(paste0(save.assem, "/scaffolds.fasta")) == TRUE ){
      system(paste0("cp ", save.assem, "/scaffolds.fasta ",  assembly.directory,
                    "/", samples[i], ".fa"))
    } else { stop(paste0("spades error for ", samples[i], ", check spades.log file in spades-assembly folder.")) }

    if (save.corrected.reads == FALSE) {
      #copies to contigs folder
      system(paste0("rm -r ", save.assem, "/corrected"))
    } # end if

    system(paste0("rm -r ", save.assem, "/tmp"))
    print(paste0(samples[i], " Completed Spades asssembly!"))

  }#end sample loop

}#end function

