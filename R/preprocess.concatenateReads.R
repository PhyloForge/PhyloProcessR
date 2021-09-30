#' @title concatenateReads
#'
#' @description Function for running the program spades to assemble short read sequencing data
#'
#' @param input.reads directory of processed reads
#'
#' @param output.directory save name for the output directory
#'
#' @param assembly.directory save name for the output directory
#'
#' @param spades.path contigs are added into existing alignment if algorithm is "add"
#'
#' @param fastqsplitter.path contigs are added into existing alignment if algorithm is "add"
#'
#' @param number.chunks contigs are added into existing alignment if algorithm is "add"
#'
#' @param mismatch.corrector algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param kmer.values if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param threads number of computation processing threads
#'
#' @param memory amount of system memory to use
#'
#' @param overwrite TRUE to overwrite a folder of samples with output.dir
#'
#' @param save.corrected.reads TRUE to overwrite a folder of samples with output.dir
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

#### NOT DONE
concatenateReads = function(input.reads = NULL,
                            output.directory = "processed-reads/spades-assembly",
                            assembly.directory = "draft-assemblies",
                            spades.path = "spades.py",
                            fastqsplitter.path = "fastqsplitter",
                            cap3.path = "cap3",
                            number.chunks = 1,
                            mismatch.corrector = FALSE,
                            kmer.values = c(21,33,55,77,99,127),
                            threads = 1,
                            memory = 4,
                            overwrite = FALSE,
                            save.all.files = FALSE,
                            quiet = TRUE) {

  # #debug
  # setwd("/Users/chutter/Dropbox/Research/0_Github/Test-dataset")
  # input.reads = "/Users/chutter/Dropbox/Research/0_Github/Test-dataset/processed-reads/pe-merged-reads/Test_Sample_No1"
  # spades.path = "/Users/chutter/miniconda3/bin/spades.py"
  # fastqsplitter.path = "/Users/chutter/miniconda3/bin/fastqsplitter"
  # output.directory = "processed-reads/spades-assembly"
  # assembly.directory = "draft-assemblies"
  # mismatch.corrector = F
  # kmer.values = c(21,33,55,77,99,127)
  # number.chunks = 4
  # threads = 1
  # memory = 4
  # overwrite = FALSE
  # quiet = FALSE
  # resume = TRUE
  # save.all.files = FALSE

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
  if (mismatch.corrector == TRUE){ mismatch.string = " --careful " }
  if (mismatch.corrector == FALSE){ mismatch.string = " " }

  #Sets up the reads
  files = list.files(path = input.reads, full.names = T, recursive = T)
  reads = files[grep(pattern = "fastq|fq|clustS", x = files)]
  samples = gsub(paste0(input.reads, "/"), "", reads)
  samples = unique(gsub("/.*", "", samples))
  sample.name = unique(gsub("_L00.*", "", samples))

  #Creates temp directory
  dir.create(paste0(output.directory, "/chunked-reads/"))
  sample.chunk = paste0(output.directory, "/chunked-reads/", sample.name)

  # #Skips samples already finished
  # if (resume == TRUE){
  #   done.names = list.files(assembly.directory)
  #   samples = samples[!samples %in% gsub(".fa$", "", done.names)]
  # } else { samples = samples }

  if (length(samples) == 0){ stop("No samples to run or incorrect directory.") }

  if (dir.exists(sample.chunk) == FALSE) {

    #Concatenate together
    dir.create(sample.chunk)
    read1.reads = reads[grep("_READ1", reads)]
    read2.reads = reads[grep("_READ2", reads)]
    read3.reads = reads[grep("_READ3", reads)]

    system(paste0("cat ", paste0(read1.reads, collapse = " "), " > ", sample.chunk,
                  "/", sample.name, "_ALL_READ1.fastq.gz"))
    system(paste0("cat ", paste0(read2.reads, collapse = " "), " > ", sample.chunk,
                  "/", sample.name, "_ALL_READ2.fastq.gz"))
    system(paste0("cat ", paste0(read3.reads, collapse = " "), " > ", sample.chunk,
                  "/", sample.name, "_ALL_READ3.fastq.gz"))

    ### Create chunks
    chunk.no = sprintf("%03d", rep(1:number.chunks))
    chunk.strings = paste0(sample.chunk, "/", sample.name, "_CHNK", chunk.no, "_READ1.fastq.gz")
    chunk.strings = paste0("-o ", chunk.strings, collapse = " ")
    system(paste0(fastqsplitter.path, " -i ", sample.chunk, "/", sample.name, "_ALL_READ1.fastq.gz",
                  " ", chunk.strings, " -c 4"))

    chunk.strings = paste0(sample.chunk, "/", sample.name, "_CHNK", chunk.no, "_READ2.fastq.gz")
    chunk.strings = paste0("-o ", chunk.strings, collapse = " ")
    system(paste0(fastqsplitter.path, " -i ", sample.chunk, "/", sample.name, "_ALL_READ2.fastq.gz",
                  " ", chunk.strings, " -c 4"))

    chunk.strings = paste0(sample.chunk, "/", sample.name, "_CHNK", chunk.no, "_READ3.fastq.gz")
    chunk.strings = paste0("-o ", chunk.strings, collapse = " ")
    system(paste0(fastqsplitter.path, " -i ", sample.chunk, "/", sample.name, "_ALL_READ3.fastq.gz",
                  " ", chunk.strings, " -c 4"))

    system(paste0("rm ", sample.chunk, "/", sample.name, "_ALL_*"))

  } #end if

  #Loops through each chunk
  chunk.no = sprintf("%03d", rep(1:number.chunks))
  chunk.names = paste0("CHNK", chunk.no)
  chunk.reads = list.files(sample.chunk, full.names = T)

  #Header data for features and whatnot
  for (i in 1:length(chunk.names)){

    sample.reads = chunk.reads[grep(paste0("_", chunk.names[i], "_"), chunk.reads)]

    #Run SPADES on sample
    k.val = paste(kmer.values, collapse = ",")
    save.assem = paste0(sample.chunk, "/", chunk.names[i])

    #Creates assembly reads folder if not present
    if (dir.exists(save.assem) == TRUE) {
      #Run spades
      system(paste0(spades.path, " -o ", save.assem,
                    " --continue"),
             ignore.stdout = quiet)
    } else {

      #Create new directory
      dir.create(save.assem)

      #Sorts reads
      sample.names = unique(gsub("_READ.*", "", sample.reads))

      #Creates a spades character string to run different reads and library configurations
      read.string = c()
      #Checks for different read lengths
      if (length(sample.reads) == 1){ read.string = paste0("--s1 ", sample.reads[1], " ") }
      if (length(sample.reads) >= 2){ read.string = paste0("--pe1-1 ", sample.reads[1],
                                                           " --pe1-2 ", sample.reads[2], " ") }
      if (length(sample.reads) == 3){ read.string = paste0(read.string, "--pe1-m ", sample.reads[3], " ") }

      if (file.exists(paste0(assembly.directory, "/", chunk.names[i-1], "_temp-contigs.fa")) == TRUE) {
        untrust.string = paste0("--untrusted-contigs ", assembly.directory, "/", chunk.names[i-1], "_temp-contigs.fa")
      } else { untrust.string = "" }

      #Run spades
      system(paste0(spades.path, " ", read.string,
                    "-o ", save.assem, " -k ", k.val,  mismatch.string, untrust.string,
                    " -t ", threads, " -m ", memory),
             ignore.stdout = quiet)

    }#end else

    #Crashes function if spades failed, also copies new assemblies to assembly.directory
    if (file.exists(paste0(save.assem, "/contigs.fasta")) == TRUE ){
      #cop new contig
      system(paste0("cp ", save.assem, "/contigs.fasta ",  assembly.directory,
                    "/", chunk.names[i], "_temp-contigs.fa"))
      #concatenate old
      if (file.exists(paste0(assembly.directory, "/all-chunks_temp-contigs.fa")) == TRUE ) {
        system(paste0("cat ", assembly.directory, "/all-chunks_temp-contigs.fa ",
                      assembly.directory, "/", chunk.names[i], "_temp-contigs.fa > ",
                      assembly.directory, "/temp_all-chunks_temp-contigs.fa"))
        system(paste0("rm ", assembly.directory, "/all-chunks_temp-contigs.fa"))

        #Trys to combine the contigs
        input.contigs = paste0(assembly.directory, "/temp_all-chunks_temp-contigs.fa")
        runCap3(contigs = input.contigs,
                output.name = paste0(assembly.directory, "/all-chunks_temp-contigs.fa"),
                read.R = FALSE,
                cap3.path = cap3.path)
      } else {
        system(paste0("mv ", assembly.directory, "/", chunk.names[i], "_temp-contigs.fa ",
                      assembly.directory, "/all-chunks_temp-contigs.fa"))
      }#end else
    } #end if

    if (save.all.files == FALSE) {
      #copies to contigs folder
      system(paste0("rm -rf ", save.assem, ""))
    } # end if

    print(paste0(chunk.names[i], " Completed Spades asssembly!"))

  }#end sample loop

}#end function

