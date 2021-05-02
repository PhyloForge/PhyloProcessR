#' @title mergePairedEndReads
#'
#' @description Function for merging paired-end reads from processed Illumina sequence data using the program fastp
#'
#' @param input.reads path to a folder of processed reads in fastq format.
#'
#' @param output.dir the new directory to save the adaptor trimmed sequences
#'
#' @param fastp.path system path to fastp in case it can't be found
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

mergePairedEndReads = function(input.reads = NULL,
                               output.directory = "read-processing/pe-merged-reads",
                               fastp.path = NULL,
                               threads = 1,
                               mem = 8,
                               resume = TRUE,
                               overwrite = TRUE,
                               quiet = TRUE) {

  #Debegging
  # setwd("/Users/chutter/Dropbox/Research/0_Github/Test-dataset")
  # input.reads = "/Users/chutter/Dropbox/Research/0_Github/Test-dataset/read-processing/decontaminated-reads"
  # fastp.path = "/Users/chutter/miniconda3/bin/fastp"
  # output.directory = "read-processing/pe-merged-reads"
  # mode = "directory"
  # threads = 4
  # mem = 8
  # resume = FALSE
  # overwrite = FALSE
  # quiet = TRUE

  #Same adds to bbmap path
  if (is.null(fastp.path) == FALSE){
    b.string = unlist(strsplit(fastp.path, ""))
    if (b.string[length(b.string)] != "/") {
      fastp.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { fastp.path = "" }


  #Quick checks
  options(stringsAsFactors = FALSE)
  if (is.null(input.reads) == TRUE){ stop("Please provide raw reads.") }
  if (file.exists(input.reads) == F){ stop("Input reads not found.") }
  if (is.null(file.rename) == TRUE){ stop("Please provide a table of file to sample name conversions.") }

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
  if (resume == TRUE){
    done.files = list.files(output.directory)
    sample.names = sample.names[!sample.names %in% done.files]
  }

  if (length(sample.names) == 0){ stop("no samples remain to analyze.") }

  #Creates the summary log
  summary.data =  data.frame(Sample = as.character(),
                             Lane = as.character(),
                             Task = as.character(),
                             Program = as.character(),
                             startPairs = as.numeric(),
                             removePairs = as.numeric())

  for (i in 1:length(sample.names)) {
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

      lane.reads = reads[grep(pattern = paste0(sample.reads[j], "_"), x = reads)]

      #Checks the Sample column in case already renamed
      if (length(lane.reads) == 0){ lane.reads = reads[grep(pattern = paste0(sample.reads[j], x = reads))] }
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
      outreads[1] = paste0(out.path, "/", lane.name, "_READ1.fastq.gz")
      outreads[2] = paste0(out.path, "/", lane.name, "_READ2.fastq.gz")
      outread.m = paste0(out.path, "/", lane.name, "_READ3.fastq.gz")

      #Runs fastp: only does adapter trimming, no quality stuff
      system(paste0(fastp.path, "fastp --merge --disable_adapter_trimming --disable_quality_filtering",
                    " --in1 ", lane.reads[1], " --in2 ", lane.reads[2],
                    " --out1 ", outreads[1], " --out2 ", outreads[2], " --merged_out ", outread.m,
                    " --html pe-merged_fastp.html --json pe-merged_fastp.json",
                    " --report_title ", sample.reads[j]," --thread ", threads),
             ignore.stderr = quiet, ignore.stdout = quiet)

      system(paste0("cp pe-merged_fastp.html ",
                    report.path, "/", lane.name, "_pe-merged_fastp.html"))
      system(paste0("rm pe-merged_fastp*"))

      #Gathers stats on initial data
      start.reads = as.numeric(system(paste0("zcat < ", lane.reads[1], " | echo $((`wc -l`/4))"), intern = T))
      end.reads = as.numeric(system(paste0("zcat < ", outreads[1], " | echo $((`wc -l`/4))"), intern = T))

      temp.remove = data.frame(Sample = sample.names[i],
                               Lane = gsub(".*_", "", lane.name),
                               Task = "merge-pe-reads",
                               Program = "fastp",
                               startPairs = start.reads,
                               removePairs = start.reads-end.reads,
                               endPairs = end.reads)

      summary.data = rbind(summary.data, temp.remove)

      print(paste0(lane.name, " Completed paired-end read merging!"))
    }#end sample j loop

    print(paste0(sample.names[i], " Completed paired-end read merging!"))
  }#end i loop

  write.csv(summary.data, file = paste0("logs/mergePairedEndReads_summary.csv"), row.names = FALSE)

}#end function

