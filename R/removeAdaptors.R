#' @title removeAdaptors
#'
#' @description Function for removing adaptor sequences from raw Illumina sequence data using the program fastp
#'
#' @param raw.reads path to a folder of raw reads in fastq format.
#'
#' @param output.dir the new directory to save the adaptor trimmed sequences
#'
#' @param mode "Sample" to run on a single sample or "Directory" to run on a directory of samples
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

removeAdaptors = function(raw.reads = NULL,
                          output.dir = NULL,
                          mode = c("sample", "directory"),
                          fastp.path = "fastp",
                          threads = 1,
                          mem = 8,
                          resume = TRUE,
                          overwrite = TRUE,
                          quiet = TRUE) {

  #Debegging
  # raw.reads = "/Users/chutter/Dropbox/Mammals/Philippine_Shrews_Raw_Data"
  # fastp.path = "/Users/chutter/miniconda3/bin/fastp"
  # output.dir = "adaptor-removed-reads"
  # mode = "directory"
  # threads = 4
  # mem = 8
  # resume = FALSE
  # overwrite = TRUE
  # quiet = TRUE

  #Quick checks
  options(stringsAsFactors = FALSE)
  if (is.null(raw.reads) == TRUE){ stop("Please provide raw reads.") }
  if (is.null(file.rename) == TRUE){ stop("Please provide a table of file to sample name conversions.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.dir) == F){ dir.create(output.dir) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.dir))
      dir.create(output.dir)
    }
  }#end else

  #Creates output directory
  if (dir.exists("logs") == F){ dir.create("logs") }

  #Read in sample data
  reads = list.files(raw.reads, recursive = T, full.names = T)
  sample.names = gsub(".*/", "", reads)
  sample.names = gsub("_R1_.*|_R2_.*|_READ1_.*|_READ2_.*", "", sample.names)
  sample.data = data.frame(File = sample.names, Sample = sample.names)

  #Creates the summary log
  summary.data =  data.frame(Sample = as.character(),
                             rawReads =as.character(),
                             Task = as.character(),
                             Program = as.character(),
                             startPairs = as.numeric(),
                             removePairs = as.numeric())

  for (i in 1:nrow(sample.data)) {
    #################################################
    ### Part A: prepare for loading and checks
    #################################################
    #Finds all files for this given sample and turns into a cat string
    sample.reads = reads[grep(pattern = paste0(sample.data$File[i], "_"), x = reads)]
    #Checks the Sample column in case already renamed
    if (length(sample.reads) == 0){ sample.reads = reads[grep(pattern = paste0(sample.data$Sample[i], "_"), x = reads)] }
    #Returns an error if reads are not found
    if (length(sample.reads) == 0 ){
      stop(sample.data$Sample[i], " does not have any reads present for files ",
           sample.data$File[i], " from the input spreadsheet. ")
    } #end if statement

    #################################################
    ### Part B: Create directories and move files
    #################################################
    #Create sample directory
    out.path = paste0(output.dir, "/", sample.data$Sample[i])
    dir.create(out.path)
    #Creates sample report directory
    report.path = paste0("logs/", sample.data$Sample[i])
    dir.create(report.path)

    #################################################
    ### Part C: Runs fastp
    #################################################
    #sets up output reads
    outread.1 = paste0(out.path, "/", sample.data$Sample[i], "_READ1_L001.fastq.gz")
    outread.2 = paste0(out.path, "/", sample.data$Sample[i], "_READ2_L001.fastq.gz")

    #Runs fastp: only does adapter trimming, no quality stuff
    system(paste0(fastp.path, " --in1 ",sample.reads[1], " --in2 ", sample.reads[2],
                  " --out1 ", outread.1, " --out2 ", outread.2,
                  " --length_required 30 --low_complexity_filter --complexity_threshold 30",
                  " --html adapter-trim_fastp.html --json adapter-trim_fastp.json",
                  " --report_title ", sample.data$Sample[i]," --thread ", threads),
           ignore.stderr = quiet, ignore.stdout = quiet)

    system(paste0("cp adapter-trim_fastp.html ",
                  report.path, "/adapter-trim_fastp.html"))
    system(paste0("rm adapter-trim_fastp*"))

    #Gathers stats on initial data
    start.reads = as.numeric(system(paste0("zcat < ", sample.reads[1], " | echo $((`wc -l`/4))"), intern = T))
    end.reads = as.numeric(system(paste0("zcat < ", outread.1, " | echo $((`wc -l`/4))"), intern = T))

    temp.remove = data.frame(Sample = sample.data$Sample[i],
                             rawReads = sample.data$File[i],
                             Task = "trim-adaptors+filter-complex",
                             Program = "fastp",
                             startPairs = start.reads,
                             removePairs = start.reads-end.reads,
                             endPairs = end.reads)

    summary.data = rbind(summary.data, temp.remove)

    print(paste0(sample.data$Sample[i], " Completed adaptor removal!"))
  }#end sample i loop

  write.csv(summary.data, file = paste0("logs/fastp_summary.csv"), row.names = FALSE)

}#end function
