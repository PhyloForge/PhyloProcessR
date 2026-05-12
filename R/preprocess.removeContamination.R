#' @title removeContamination
#'
#' @description Removes reads that map to a set of contaminant genomes (e.g.
#'   human, mouse, or vector sequences) using BWA-MEM and samtools. A combined
#'   BWA index is built from all FASTA files in decontamination.path, reads are
#'   mapped against it, and only the unmapped (non-contaminant) read pairs are
#'   retained and saved. Per-contaminant read counts are written to the logs
#'   directory.
#'
#' @param input.reads path to a directory of adapter-trimmed paired-end reads
#'   in fastq.gz format, organised in per-sample sub-directories.
#'
#' @param output.directory path to the directory where decontaminated reads
#'   will be saved (one sub-directory per sample).
#'
#' @param decontamination.path path to a directory of contaminant reference
#'   genome FASTA files (e.g. produced by createContaminantDB()).
#'
#' @param map.match numeric; minimum mapping identity threshold (currently
#'   passed for reference but filtering is handled by samtools flags rather than
#'   this value directly).
#'
#' @param samtools.path system path to the directory containing the samtools
#'   executable; NULL searches the system PATH.
#'
#' @param bwa.path system path to the directory containing the bwa executable;
#'   NULL searches the system PATH.
#'
#' @param threads number of CPU threads for BWA and samtools.
#'
#' @param mem amount of RAM in GB (currently reserved).
#'
#' @param overwrite logical; if TRUE the output directory is deleted and
#'   recreated before processing.
#'
#' @param overwrite.reference logical; if TRUE an existing BWA index in
#'   ref-index/ is deleted and rebuilt.
#'
#' @param quiet logical; if TRUE BWA and samtools stdout/stderr are suppressed.
#'
#' @return invisibly; writes decontaminated fastq.gz files to output.directory,
#'   per-sample contamination counts to logs/, and a CSV summary to
#'   logs/removeContamination_summary.csv.
#'
#' @export

removeContamination = function(input.reads = "cleaned-reads",
                               output.directory = "decontaminated-reads",
                               decontamination.path = NULL,
                               map.match = 0.90,
                               samtools.path = NULL,
                               bwa.path = NULL,
                               threads = 1,
                               mem = 1,
                               overwrite = FALSE,
                               overwrite.reference = FALSE,
                               quiet = TRUE) {

  #Debug
  # setwd("/Users/chutter/Dropbox/Research/0_Github/Test-dataset")
  # read.directory = "read-processing/adaptor-removed-reads"
  # output.directory = "read-processing/decontaminated-reads"
  # decontamination.path = "/Users/chutter/Dropbox/Research/0_Github/Contamination_Genomes"
  # samtools.path = "/Users/chutter/miniconda3/bin/samtools"
  # bwa.path = "/usr/local/bin/bwa"
  # mode = "directory"
  # threads = 4
  # mem = 8
  # resume = TRUE
  # overwrite = FALSE
  # quiet = TRUE
  # map.match = 0.99

  #Same adds to bbmap path
  if (is.null(samtools.path) == FALSE){
    b.string = unlist(strsplit(samtools.path, ""))
    if (b.string[length(b.string)] != "/") {
      samtools.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { samtools.path = "" }

  #Same adds to bbmap path
  if (is.null(bwa.path) == FALSE){
    b.string = unlist(strsplit(bwa.path, ""))
    if (b.string[length(b.string)] != "/") {
      bwa.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { bwa.path = "" }

  #Quick checks
  options(stringsAsFactors = FALSE)
  if (is.null(input.reads) == TRUE){ stop("Please provide input reads.") }
  if (file.exists(input.reads) == F){ stop("Input reads not found.") }
  if (is.null(decontamination.path) == TRUE){ stop("Please provide decontamination genomes / sequences.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.directory) == F){ dir.create(output.directory) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  }#end else

  #Creates output directory
  if (dir.exists("logs/sample_logs") == F){ dir.create("logs/sample_logs", recursive = TRUE) }

  #Read in sample data **** sample is run twice?!
  reads = list.files(input.reads, recursive = T, full.names = T)
  sample.names = list.files(input.reads, recursive = F, full.names = F)

  #Resumes file download
  if (overwrite == FALSE){
    done.files = list.files(output.directory)
    sample.names = sample.names[!sample.names %in% done.files]
  }

  if (length(sample.names) == 0){ return("no samples remain to analyze.") }

  #Create combined and indexed reference
  if (dir.exists("ref-index") == TRUE && overwrite.reference == TRUE){ system(paste0("rm -rf ref-index")) }

  if (dir.exists("ref-index") == FALSE){
    dir.create("ref-index")
    reference.list = list.files(decontamination.path, full.names = TRUE)

    #Build contig-to-genome mapping before concatenating
    contig.map = data.frame(Contig = character(), Genome = character(), Accession = character())
    for (ref.file in reference.list) {
      file.base = gsub("\\.fna\\.gz$|\\.fa\\.gz$|\\.fa$", "", basename(ref.file))
      if (grepl("-", file.base)) {
        parts = strsplit(file.base, "-")[[1]]
        genome.name = parts[1]
        accession = paste(parts[-1], collapse = "-")
      } else {
        genome.name = file.base
        accession = file.base
      }
      headers = system(paste0("zcat -f ", shQuote(ref.file), " | grep '^>'"), intern = TRUE)
      contig.names = gsub("^>([^ ]+).*", "\\1", headers)
      contig.map = rbind(contig.map, data.frame(Contig = contig.names, Genome = genome.name, Accession = accession))
    }
    write.csv(contig.map, file = "ref-index/contig_mapping.csv", row.names = FALSE)

    system(paste0("cat ", paste(shQuote(reference.list), collapse = " "),
                  " | zcat -f > ref-index/reference.fa"))
    system(paste0(bwa.path, "bwa index -p ref-index/reference ref-index/reference.fa"),
           ignore.stderr = quiet, ignore.stdout = quiet)
  }#end dir exists false

  #Creates the summary log
  summary.data =  data.frame(Sample = as.character(),
                             Lane = as.character(),
                             Task = as.character(),
                             Program = as.character(),
                             startPairs = as.numeric(),
                             removePairs = as.numeric())

  contam.summary = data.frame(Sample = as.character(),
                              Lane = as.character(),
                              Contaminant = as.character(),
                              Accession = as.character(),
                              ReadPairs = as.numeric())

  contig.map = read.csv("ref-index/contig_mapping.csv")

  #Runs through each sample
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
    report.path = paste0("logs/sample_logs/", sample.names[i])
    if (file.exists(out.path) == FALSE) { dir.create(out.path) }
    if (file.exists(report.path) == FALSE) { dir.create(report.path) }

    for (j in 1:length(sample.reads)){

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
      outreads[1] = paste0(out.path, "/", lane.name, "_READ1.fastq.gz")
      outreads[2] = paste0(out.path, "/", lane.name, "_READ2.fastq.gz")

      #BWA mapping
      system(paste0(bwa.path, "bwa mem -M",
                    " -t ", threads, " ref-index/reference ",
                    lane.reads[1], " ", lane.reads[2], " | ", samtools.path ,"samtools sort -@ ", threads,
                    " -o ", out.path, "/decontam-all.bam"),  ignore.stderr = quiet, ignore.stdout = quiet)

      system(paste0(samtools.path, "samtools index ", out.path, "/decontam-all.bam"),  ignore.stderr = quiet, ignore.stdout = quiet)

      #Convert map.match identity threshold to max allowed mismatch rate for samtools -e filter
      max.mismatch = 1 - map.match

      #Extract mapped reads (contam reads), extracts both mapped
      system(paste0(samtools.path, "samtools view -b -f 1 -F 12 -e '[NM]/qlen<=", max.mismatch, "' ",
                    out.path, "/decontam-all.bam | ",
                    samtools.path, "samtools sort -@ ", threads, " -o ", out.path, "/decontam-mapped-sort-1.bam"),
             ignore.stderr = quiet, ignore.stdout = quiet)

      #Extracts R2 mapped
      system(paste0(samtools.path, "samtools view -b -f 4 -F 264 -e '[NM]/qlen<=", max.mismatch, "' ",
                    out.path, "/decontam-all.bam | ",
                    samtools.path, "samtools sort -@ ", threads, " -o ", out.path, "/decontam-mapped-sort-2.bam"),
             ignore.stderr = quiet, ignore.stdout = quiet)

      #Extracts R3 mapped
      system(paste0(samtools.path, "samtools view -b -f 8 -F 260 -e '[NM]/qlen<=", max.mismatch, "' ",
                    out.path, "/decontam-all.bam | ",
                    samtools.path, "samtools sort -@ ", threads, " -o ", out.path, "/decontam-mapped-sort-3.bam"),
             ignore.stderr = quiet, ignore.stdout = quiet)

      system(paste0(samtools.path, "samtools merge -f ",
                    out.path, "/decontam-mapped-sort.bam ",
                    out.path, "/decontam-mapped-sort-1.bam ",
                    out.path, "/decontam-mapped-sort-2.bam ",
                    out.path, "/decontam-mapped-sort-3.bam"),
             ignore.stderr = quiet, ignore.stdout = quiet)

      system(paste0(samtools.path, "samtools index ", out.path, "/decontam-mapped-sort.bam"),
             ignore.stderr = quiet, ignore.stdout = quiet)

      #Extract unmapped reads (good reads): extracts unmapped for R2 and R1
      system(paste0(samtools.path, "samtools view -b -f 12 -F 256 ", out.path, "/decontam-all.bam | ",
                    samtools.path, "samtools sort -@ ", threads, " -n -o ", out.path, "/decontam-unmapped-sort.bam"),
             ignore.stderr = quiet, ignore.stdout = quiet)

      #Saves as fastq
      system(paste0(samtools.path, "samtools fastq -@ ", threads, " ",
                    out.path, "/decontam-unmapped-sort.bam -1 ", outreads[1], " -2 ", outreads[2]),
             ignore.stderr = quiet, ignore.stdout = quiet)

      #Checks stats
      #system(paste0(samtools.path, " flagstat ", out.path, "/decontam-all.bam"))
      system(paste0(samtools.path, "samtools view -c ", out.path, "/decontam-all.bam"),
             ignore.stderr = quiet, ignore.stdout = quiet)

      #system(paste0(samtools.path, " flagstat ", out.path, "/decontam-unmapped-sort.bam"))
      system(paste0(samtools.path, "samtools view -c ", out.path, "/decontam-unmapped-sort.bam"),
             ignore.stderr = quiet, ignore.stdout = quiet)

      #system(paste0(samtools.path, " flagstat ", out.path, "/decontam-mapped-sort.bam"))
      system(paste0(samtools.path, "samtools view -c ", out.path, "/decontam-mapped-sort.bam"),
             ignore.stderr = quiet, ignore.stdout = quiet)

      #Gets match statistics
      table.dat = (system(paste0(samtools.path, "samtools idxstats ",
                                 out.path, "/decontam-mapped-sort.bam | cut -f 1,3"), intern = T))
      table.dat = data.frame(Organism = gsub("\t.*", "", table.dat),
                             Count = as.numeric(gsub(".*\t", "", table.dat)) )
      table.dat = merge(table.dat, contig.map, by.x = "Organism", by.y = "Contig", all.x = TRUE)
      table.dat$Genome[is.na(table.dat$Genome)] = table.dat$Organism[is.na(table.dat$Genome)]
      table.dat$Accession[is.na(table.dat$Accession)] = table.dat$Organism[is.na(table.dat$Accession)]

      contam.data = aggregate(table.dat$Count, FUN = sum, by = list(table.dat$Genome, table.dat$Accession))
      colnames(contam.data) = c("Contaminant", "Accession", "ReadPairs")
      contam.data = contam.data[contam.data$Contaminant != "*",]

      system(paste0("rm ", out.path, "/decontam-*"))

      #Gathers stats on initial data
      start.reads = as.numeric(system(paste0("zcat < ", lane.reads[1], " | echo $((`wc -l`/4))"), intern = T))
      end.reads = as.numeric(system(paste0("zcat < ", outreads[1], " | echo $((`wc -l`/4))"), intern = T))

      temp.remove = data.frame(Sample = sample.names[i],
                               Lane = gsub(".*_", "", lane.name),
                               Task = "decontamination",
                               Program = "bwa",
                               startPairs = start.reads,
                               removePairs = start.reads-end.reads,
                               endPairs = end.reads)

      summary.data = rbind(summary.data, temp.remove)

      lane.contam = contam.data
      lane.contam$Sample = sample.names[i]
      lane.contam$Lane = gsub(".*_", "", lane.name)
      contam.summary = rbind(contam.summary, lane.contam[, c("Sample", "Lane", "Contaminant", "Accession", "ReadPairs")])

      write.csv(contam.data, file = paste0(report.path, "/contamination-read-counts.csv"), row.names = FALSE)

    }#end sample j loop

    print(paste0(sample.names[i], " Completed decontamination removal!"))

  }#end sample i loop

  write.csv(summary.data, file = "logs/removeContamination_summary.csv", row.names = FALSE)
  write.csv(contam.summary, file = "logs/removeContamination_contaminants.csv", row.names = FALSE)
}

