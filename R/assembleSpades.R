#' @title assembleSpades
#'
#' @description Function for running the program spades to assemble short read sequencing data
#'
#' @param read.directory directory of processed reads
#'
#' @param output.directory save name for the output directory
#'
#' @param full.path.spades contigs are added into existing alignment if algorithm is "add"
#'
#' @param mismatch.corrector algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
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

assembleSpades = function(input.reads = NULL,
                          output.directory = "processed-reads/spades-assembly",
                          assembly.directory = "draft-assemblies",
                          spades.path = "spades.py",
                          mismatch.corrector = TRUE,
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

  #mismatch corrector
  if (mismatch.corrector == TRUE){ mismatch.string = "--careful " }
  if (mismatch.corrector == FALSE){ mismatch.string = " " }

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

    sample.reads = reads[grep(pattern = paste0(samples[i], "_"), x = reads)]

    #Checks the Sample column in case already renamed
    if (length(sample.reads) == 0){ sample.reads = reads[grep(pattern = samples[i], x = reads)] }

    #sample.reads = unique(gsub("_1.f.*|_2.f.*|_3.f.*|_R1_.*|_R2_.*|_R3_.*|_READ1_.*|_READ2_.*|_READ3_.*|_R1.fast.*|_R2.fast.*|_R3.fast.*|_READ1.fast.*|_READ2.fast.*|_READ3.fast.*", "", sample.reads))

    #Returns an error if reads are not found
    if (length(sample.reads) == 0 ){
      stop(sample.names[i], " does not have any reads present for files ")
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

      lib.read1 = lib.reads[grep("_1.f.*|-1.f.*|_R1_.*|_READ1_.*|_R1.fast.*|-R1.fast.*|_READ1.fast.*|-READ1.fast.*", sample.reads)]
      lib.read2 = lib.reads[grep("_2.f.*|-2.f.*|_R2_.*|_READ2_.*|_R2.fast.*|-R2.fast.*|_READ2.fast.*|-READ2.fast.*", sample.reads)]
      lib.read3 = lib.reads[grep("_3.f.*|-3.f.*|_R3_.*|_READ3_.*|_R3.fast.*|-R3.fast.*|_READ3.fast.*|-READ3.fast.*|_singleton.*|-singleton.*|READ-singleton.*|READ_singleton.*|_READ-singleton.*|-READ_singleton.*|-READ-singleton.*|_READ_singleton.*", sample.reads)]

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
    system(paste0(spades.path, " ", final.read.string,
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

