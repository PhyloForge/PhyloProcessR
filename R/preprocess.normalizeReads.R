#' @title normalizeReads
#'
#' @description Function for removing contamination from other organisms from adaptor trimmed Illumina sequence data using BWA
#'
#' @param input.reads path to a folder of adaptor trimmed reads in fastq format.
#'
#' @param output.directory the new directory to save the adaptor trimmed sequences
#'
#' @param decontamination.path directory of genomes contaminants to scan samples
#'
#' @param samtools.path system path to samtools in case it can't be found
#'
#' @param bwa.path system path to bwa in case it can't be found
#'
#' @param threads number of computation processing threads
#'
#' @param memory amount of system memory to use
#'
#' @param resume TRUE to skip samples already completed
#'
#' @param overwrite TRUE to overwrite a folder of samples with output.dir
#'
#' @param quiet TRUE to supress screen output
#'
#' @return a new directory of adaptor trimmed reads and a summary of the trimming in logs/
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

normalizeReads = function(input.reads = NULL,
                          output.directory = "normalized-reads",
                          orna.path = NULL,
                          threads = 1,
                          memory = 1,
                          overwrite = FALSE,
                          quiet = TRUE) {

  # #Debug
  # setwd("/Volumes/LaCie/Mantellidae_Genome")
  # input.reads = "adaptor-removed-reads"
  # output.directory = "normalized-reads"
  # orna.path = "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"
  # threads = 4
  # memory = 8
  # resume = TRUE
  # overwrite = FALSE
  # quiet = TRUE

  #Same adds to bbmap path
  if (is.null(orna.path) == FALSE){
    b.string = unlist(strsplit(orna.path, ""))
    if (b.string[length(b.string)] != "/") {
      orna.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { orna.path = "" }

  #Quick checks
  if (is.null(input.reads) == TRUE){ stop("Please provide input reads.") }
  if (dir.exists(input.reads) == F){ stop("Input reads not found.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.directory) == F){ dir.create(output.directory) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  }#end else

  #Creates output directory
  if (dir.exists("logs") == F){ dir.create("logs") }

  #Read in sample data **** sample is run twice?!
  reads = list.files(input.reads, recursive = T, full.names = T)
  sample.names = list.files(input.reads, recursive = F, full.names = F)

  #Resumes file download
  if (overwrite == FALSE){
    done.files = list.files(output.directory)
    sample.names = sample.names[!sample.names %in% done.files]
  }

  if (length(sample.names) == 0){ return("no samples remain to analyze.") }

  #Creates the summary log
  summary.data =  data.frame(Sample = as.character(),
                             Lane = as.character(),
                             Task = as.character(),
                             Program = as.character(),
                             startPairs = as.numeric(),
                             removePairs = as.numeric())

  # #Sets up multiprocessing
  cl = parallel::makeCluster(threads, outfile = "")
  doParallel::registerDoParallel(cl)
  mem.cl = floor(memory/threads)

  foreach::foreach(i=1:length(sample.names), .packages = c("PhyloCap", "foreach", "Biostrings","Rsamtools", "ape", "stringr", "data.table")) %dopar% {
    #Runs through each sample
    #for (i in 1:length(sample.names)) {
    #################################################
    ### Part A: prepare for loading and checks
    #################################################
    sample.reads = reads[grep(pattern = paste0(sample.names[i], "_"), x = reads)]

    #Checks the Sample column in case already renamed
    if (length(sample.reads) == 0){ sample.reads = reads[grep(pattern = sample.names[i], x = reads)] }

    sample.reads = unique(gsub("_1.f.*|_2.f.*|_3.f.*|-1.f.*|-2.f.*|-3.f.*|_R1_.*|_R2_.*|_R3_.*|_READ1_.*|_READ2_.*|_READ3_.*|_R1.f.*|_R2.f.*|_R3.f.*|-R1.f.*|-R2.f.*|-R3.f.*|_READ1.f.*|_READ2.f.*|_READ3.f.*|-READ1.f.*|-READ2.f.*|-READ3.f.*|_singleton.*|-singleton.*|READ-singleton.*|READ_singleton.*|_READ-singleton.*|-READ_singleton.*|-READ-singleton.*|_READ_singleton.*", "", sample.reads))

    #Returns an error if reads are not found
    if (length(sample.reads) == 0 ){
      stop(sample.names[i], " does not have any reads present for files ")
    } #end if statement

    #CReates new directory
    out.path = paste0(output.directory, "/", sample.names[i])
    report.path = paste0("logs/", sample.names[i])
    if (file.exists(out.path) == FALSE) { dir.create(out.path) }
    if (file.exists(report.path) == FALSE) { dir.create(report.path) }

    for (j in 1:length(sample.reads)){
      #Gets reads for the lane / separately sequenced set
      lane.reads = reads[grep(pattern = paste0(sample.reads[j], "_"), x = reads)]

      #Checks the Sample column in case already renamed
      if (length(lane.reads) == 0){ lane.reads = reads[grep(pattern = sample.reads[j], x = reads)] }
      #Returns an error if reads are not found
      if (length(lane.reads) == 0 ){
        stop(sample.reads[j], " does not have any reads present for files ")
      } #end if statement

      lane.save = gsub(".*/", "", lane.reads)
      lane.name = gsub(".*/", "", sample.reads[j])

      #################################################
      ### Part C: Runs fastp
      #################################################
      #sets up output reads
      outreads = paste0(out.path, "/", lane.name)

      #Runs fastp: only does adapter trimming, no quality stuff
      system(paste0(orna.path, "ORNA -type fastq",
                    " -pair1 ",lane.reads[1], " -pair2 ", lane.reads[2],
                    " -output ", outreads),
             ignore.stderr = quiet, ignore.stdout = quiet)

      system(paste0("gzip ", outreads, "_1.fq"))
      system(paste0("mv ", outreads, "_1.fq.gz", " ", outreads, "_READ1.fastq.gz"))
      system(paste0("gzip ", outreads, "_2.fq"))
      system(paste0("mv ", outreads, "_2.fq.gz", " ", outreads, "_READ2.fastq.gz"))

      #Gathers stats on initial data
      start.reads = as.numeric(system(paste0("zcat < ", lane.reads[1], " | echo $((`wc -l`/4))"), intern = T))
      end.reads = as.numeric(system(paste0("zcat < ", outreads, "_READ1.fastq.gz | echo $((`wc -l`/4))"), intern = T))

      temp.remove = data.frame(Sample = sample.names[i],
                               Lane = gsub(".*_", "", lane.name),
                               Task = "normalize",
                               Program = "ORNA",
                               startPairs = start.reads,
                               removePairs = start.reads-end.reads,
                               endPairs = end.reads)

      summary.data = rbind(summary.data, temp.remove)

      system(paste0("rm ", lane.name, "_READ2.fastq.h5"))

    }#end sample j loop

    print(paste0(sample.names[i], " Completed read normalization!"))

  }#end sample i loop

  parallel::stopCluster(cl)


  write.csv(summary.data, file = paste0("logs/normalizeReads_summary.csv"), row.names = FALSE)
}

