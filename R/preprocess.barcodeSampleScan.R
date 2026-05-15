#' @title barcodeSampleScan
#'
#' @description Assembles a target barcode marker (e.g. 16S rRNA, COI) from
#'   cleaned reads for each sample, then identifies the best-matching species
#'   via BLAST against a user-supplied reference database. Useful for verifying
#'   sample identity and flagging mislabelled or cross-contaminated libraries
#'   before committing to full assembly.
#'
#'   For each sample: (1) reads are mapped to the barcode reference with BWA to
#'   fish out on-target reads; (2) those reads are extracted and assembled with
#'   SPAdes; (3) assembled contigs are queried against the identification
#'   database with dc-megablast. A cross-sample summary is appended to
#'   logs/barcodeSampleScan_summary.csv after every sample so the file grows
#'   safely across successive single-sample runs.
#'
#' @param input.reads path to a directory of cleaned reads. Each sample must
#'   occupy its own sub-directory or be identified by a shared filename prefix.
#'
#' @param output.directory path to the directory where per-sample assembly and
#'   BLAST results will be saved. Default: \code{"barcode-assessment"}.
#'
#' @param barcode.fasta path to a FASTA file of target barcode reference
#'   sequence(s) (e.g. a single 16S or COI representative). Used to recruit
#'   matching reads before assembly.
#'
#' @param database.fasta path to a FASTA file of named barcode sequences used
#'   as the BLAST identification database. Species names should be in the FASTA
#'   headers. If \code{NULL} (the default), BLAST is run remotely against the
#'   NCBI \code{nt} database — no local database is required, but an internet
#'   connection is needed and queries will be slower due to NCBI rate limits.
#'
#' @param hits.per.sample number of top BLAST hits to keep per sample via
#'   \code{-max_target_seqs}. Default: \code{5}.
#'
#' @param min.mapping.reads minimum number of reads that must map to the
#'   barcode reference for assembly to proceed. Samples below this threshold
#'   are recorded as low-coverage and skipped. Default: \code{10}.
#'
#' @param bwa.path system path to the directory containing the \code{bwa}
#'   executable; NULL searches the system PATH.
#'
#' @param samtools.path system path to the directory containing
#'   \code{samtools}; NULL searches the system PATH.
#'
#' @param spades.path system path to the directory containing
#'   \code{spades.py}; NULL searches the system PATH.
#'
#' @param blast.path system path to the directory containing \code{blastn} and
#'   \code{makeblastdb}; NULL searches the system PATH.
#'
#' @param threads number of CPU threads to pass to BWA, samtools, SPAdes, and
#'   BLAST. Default: \code{1}.
#'
#' @param mem amount of RAM in GB passed to SPAdes with \code{-m}.
#'   Default: \code{8}.
#'
#' @param overwrite logical; if TRUE the output directory is deleted and
#'   recreated before processing. Default: \code{FALSE}.
#'
#' @param quiet logical; if TRUE tool screen output is suppressed.
#'   Default: \code{TRUE}.
#'
#' @return invisibly; per-sample contig FASTAs and BLAST result tables are
#'   written to output.directory, and a cross-sample summary is appended to
#'   logs/barcodeSampleScan_summary.csv.
#'
#' @export

barcodeSampleScan = function(input.reads = NULL,
                             output.directory = "barcode-assessment",
                             barcode.fasta = NULL,
                             database.fasta = NULL,
                             hits.per.sample = 5,
                             min.mapping.reads = 10,
                             bwa.path = NULL,
                             samtools.path = NULL,
                             spades.path = NULL,
                             blast.path = NULL,
                             threads = 1,
                             mem = 8,
                             overwrite = FALSE,
                             quiet = TRUE) {

  # #Debug
  # setwd("/Users/chutter/Dropbox/Research/0_Github/Test-dataset")
  # input.reads = "processed-reads/cleaned-reads"
  # output.directory = "barcode-assessment"
  # barcode.fasta = "/path/to/16S_reference.fa"
  # database.fasta = "/path/to/16S_database.fa"
  # bwa.path = "/Users/chutter/miniconda3/bin"
  # samtools.path = "/Users/chutter/miniconda3/bin"
  # spades.path = "/Users/chutter/miniconda3/bin"
  # blast.path = "/Users/chutter/miniconda3/bin"
  # threads = 4
  # mem = 8
  # overwrite = TRUE
  # quiet = TRUE

  # Adds trailing slash to tool paths
  if (is.null(bwa.path) == FALSE) {
    b.string = unlist(strsplit(bwa.path, ""))
    if (b.string[length(b.string)] != "/") { bwa.path = paste0(append(b.string, "/"), collapse = "") }
  } else { bwa.path = "" }

  if (is.null(samtools.path) == FALSE) {
    b.string = unlist(strsplit(samtools.path, ""))
    if (b.string[length(b.string)] != "/") { samtools.path = paste0(append(b.string, "/"), collapse = "") }
  } else { samtools.path = "" }

  if (is.null(spades.path) == FALSE) {
    b.string = unlist(strsplit(spades.path, ""))
    if (b.string[length(b.string)] != "/") { spades.path = paste0(append(b.string, "/"), collapse = "") }
  } else { spades.path = "" }

  if (is.null(blast.path) == FALSE) {
    b.string = unlist(strsplit(blast.path, ""))
    if (b.string[length(b.string)] != "/") { blast.path = paste0(append(b.string, "/"), collapse = "") }
  } else { blast.path = "" }

  # Quick checks
  if (is.null(input.reads) == TRUE) { stop("Please provide input reads.") }
  if (is.null(barcode.fasta) == TRUE) { stop("Please provide a barcode reference FASTA file.") }
  if (file.exists(barcode.fasta) == FALSE) { stop("Barcode reference FASTA not found.") }

  # Determine BLAST mode: remote NCBI nt or local database
  use.remote.blast = is.null(database.fasta)
  if (!use.remote.blast && !file.exists(database.fasta)) {
    stop("Barcode database FASTA not found: ", database.fasta)
  }

  # Sets up output directory
  if (dir.exists(output.directory) == FALSE) {
    dir.create(output.directory)
  } else {
    if (overwrite == TRUE) {
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  }#end else

  # Creates log directories
  if (dir.exists("logs/sample_logs") == FALSE) { dir.create("logs/sample_logs", recursive = TRUE) }

  #################################################
  ### One-time index and database setup
  #################################################

  # BWA index of barcode reference — built once, reused across all samples
  index.path = paste0(output.directory, "/barcode-index")
  if (dir.exists(index.path) == FALSE) { dir.create(index.path) }
  if (!file.exists(paste0(index.path, "/barcode.fa.bwt"))) {
    system(paste0("cp ", barcode.fasta, " ", index.path, "/barcode.fa"))
    system(paste0(bwa.path, "bwa index ", index.path, "/barcode.fa"),
           ignore.stdout = quiet, ignore.stderr = quiet)
  }

  # BLAST database for identification — built once from a local FASTA, or
  # left unset when querying NCBI nt remotely.
  db.path = paste0(output.directory, "/blast-db")
  if (!use.remote.blast) {
    if (dir.exists(db.path) == FALSE) {
      dir.create(db.path)
      system(paste0("cp ", database.fasta, " ", db.path, "/database.fa"))
      system(paste0(blast.path, "makeblastdb -in ", db.path, "/database.fa",
                    " -parse_seqids -dbtype nucl -out ", db.path, "/barcode_db"),
             ignore.stdout = quiet, ignore.stderr = quiet)
    }
  } else {
    cat(" Barcode BLAST mode: remote NCBI nt (internet connection required)\n")
  }

  # Read in sample data
  reads = list.files(input.reads, recursive = TRUE, full.names = TRUE)
  sample.names = list.dirs(input.reads, recursive = FALSE, full.names = FALSE)

  if (length(sample.names) == 0) {
    sample.names = list.files(input.reads, recursive = FALSE, full.names = FALSE)
    sample.names = unique(gsub("_L00.*|_R[12][._].*|_READ[123][._].*|\\.fastq.*|\\.fq.*", "", sample.names))
    sample.names = sample.names[nchar(sample.names) > 0]
  }

  # Resumes: skip samples already done (ignore index and db subdirectories)
  if (overwrite == FALSE) {
    done.files = list.files(output.directory)
    done.files = done.files[!done.files %in% c("barcode-index", "blast-db")]
    sample.names = sample.names[!sample.names %in% done.files]
  }

  if (length(sample.names) == 0) { return("No samples remain to analyze.") }

  blast.headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
                    "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

  for (i in 1:length(sample.names)) {
    #################################################
    ### Part A: prepare reads
    #################################################
    sample.reads = reads[grep(pattern = paste0(sample.names[i], "_"), x = reads)]
    if (length(sample.reads) == 0) { sample.reads = reads[grep(pattern = sample.names[i], x = reads)] }

    sample.reads = unique(gsub("_1.f.*|_2.f.*|_3.f.*|-1.f.*|-2.f.*|-3.f.*|_R1_.*|_R2_.*|_R3_.*|_READ1_.*|_READ2_.*|_READ3_.*|_R1.f.*|_R2.f.*|_R3.f.*|-R1.f.*|-R2.f.*|-R3.f.*|_READ1.f.*|_READ2.f.*|_READ3.f.*|-READ1.f.*|-READ2.f.*|-READ3.f.*|_singleton.*|-singleton.*|READ-singleton.*|READ_singleton.*|_READ-singleton.*|-READ_singleton.*|-READ-singleton.*|_READ_singleton.*", "", sample.reads))

    if (length(sample.reads) == 0) {
      warning(sample.names[i], " does not have any reads present. Skipping.")
      next
    }#end if

    # Check for empty or near-empty files
    raw.reads = reads[grep(pattern = sample.names[i], x = reads)]
    file.sizes = file.info(raw.reads)$size
    if (any(is.na(file.sizes)) || max(file.sizes, na.rm = TRUE) < 1000) {
      warning(sample.names[i], " has empty input read files. Skipping.")
      next
    }

    out.path = paste0(output.directory, "/", sample.names[i])
    if (file.exists(out.path) == FALSE) { dir.create(out.path) }

    # Collect read1/read2 across all lanes
    all.read1 = c()
    all.read2 = c()
    for (j in 1:length(sample.reads)) {
      lane.reads = reads[grep(pattern = paste0(sample.reads[j], "_"), x = reads)]
      if (length(lane.reads) == 0) { lane.reads = reads[grep(pattern = sample.reads[j], x = reads)] }

      r1 = lane.reads[grep("_1.f.*|-1.f.*|_R1_.*|-R1_.*|_R1-.*|-R1-.*|READ1.*|_R1.fast.*|-R1.fast.*", lane.reads)]
      r2 = lane.reads[grep("_2.f.*|-2.f.*|_R2_.*|-R2_.*|_R2-.*|-R2-.*|READ2.*|_R2.fast.*|-R2.fast.*", lane.reads)]
      if (length(r1) == 1 && length(r2) == 1) {
        all.read1 = c(all.read1, r1)
        all.read2 = c(all.read2, r2)
      }
    }#end j loop

    if (length(all.read1) == 0 || length(all.read2) == 0) {
      warning(sample.names[i], " read pairs could not be identified. Skipping.")
      next
    }

    # Concatenate multi-lane reads into single files for mapping
    if (length(all.read1) > 1) {
      r1.in = paste0(out.path, "/combined_READ1.fastq.gz")
      r2.in = paste0(out.path, "/combined_READ2.fastq.gz")
      system(paste0("cat ", paste(all.read1, collapse = " "), " > ", r1.in))
      system(paste0("cat ", paste(all.read2, collapse = " "), " > ", r2.in))
    } else {
      r1.in = all.read1
      r2.in = all.read2
    }

    #################################################
    ### Part B: map reads to barcode reference
    #################################################
    bam.file = paste0(out.path, "/barcode_mapped.bam")
    system(paste0(bwa.path, "bwa mem -M -t ", threads, " ",
                  index.path, "/barcode.fa ",
                  r1.in, " ", r2.in,
                  " | ", samtools.path, "samtools sort -@", threads, " -O BAM",
                  " -o ", bam.file, " -"),
           ignore.stdout = quiet, ignore.stderr = quiet)

    system(paste0(samtools.path, "samtools index ", bam.file),
           ignore.stdout = quiet, ignore.stderr = quiet)

    # Count mapped reads
    n.mapped = as.integer(trimws(
      system(paste0(samtools.path, "samtools view -c -F 4 ", bam.file), intern = TRUE)
    ))

    if (is.na(n.mapped) || n.mapped < min.mapping.reads) {
      warning(sample.names[i], " had too few reads map to barcode reference (",
              n.mapped, " < ", min.mapping.reads, "). Skipping assembly.")
      writeLines(paste0("Too few reads mapped to barcode reference: ", n.mapped,
                        " (minimum: ", min.mapping.reads, ")"),
                 paste0("logs/sample_logs/", sample.names[i], "_barcode-lowmap.txt"))
      system(paste0("rm ", bam.file, " ", bam.file, ".bai"))
      if (length(all.read1) > 1) { unlink(c(r1.in, r2.in)) }

      # Record in summary as low-coverage
      temp.row = data.frame(Sample = sample.names[i], MappedReads = n.mapped,
                            ContigLength = NA_integer_, BestMatch = "Low-coverage",
                            Pident = NA_real_, AlignLength = NA_integer_,
                            Evalue = NA_real_, Bitscore = NA_real_,
                            stringsAsFactors = FALSE)
      out.csv = "logs/barcodeSampleScan_summary.csv"
      if (file.exists(out.csv)) {
        existing = read.csv(out.csv, stringsAsFactors = FALSE)
        existing = existing[!existing$Sample %in% temp.row$Sample, ]
        temp.row = rbind(existing, temp.row)
      }
      write.csv(temp.row, file = out.csv, row.names = FALSE)
      next
    }

    #################################################
    ### Part C: extract mapped reads and assemble
    #################################################
    # Sort by name for paired FASTQ extraction
    bam.namesort = paste0(out.path, "/barcode_namesort.bam")
    system(paste0(samtools.path, "samtools sort -n -@", threads,
                  " -o ", bam.namesort, " ", bam.file),
           ignore.stdout = quiet, ignore.stderr = quiet)

    fq1 = paste0(out.path, "/barcode_READ1.fastq")
    fq2 = paste0(out.path, "/barcode_READ2.fastq")
    system(paste0(samtools.path, "samtools fastq -@ ", threads,
                  " -1 ", fq1, " -2 ", fq2,
                  " -0 /dev/null -s /dev/null -n -F 4 ", bam.namesort),
           ignore.stdout = quiet, ignore.stderr = quiet)

    system(paste0("rm ", bam.file, " ", bam.file, ".bai ", bam.namesort))
    if (length(all.read1) > 1) { unlink(c(r1.in, r2.in)) }

    # SPAdes assembly of on-target reads
    spades.dir = paste0(out.path, "/spades")
    system(paste0(spades.path, "spades.py --careful",
                  " -1 ", fq1, " -2 ", fq2,
                  " -o ", spades.dir,
                  " -t ", threads, " -m ", mem),
           ignore.stdout = quiet, ignore.stderr = quiet)

    system(paste0("rm -f ", fq1, " ", fq2))

    contigs.file = paste0(spades.dir, "/scaffolds.fasta")
    if (file.exists(contigs.file) == FALSE) {
      warning(sample.names[i], " SPAdes assembly failed. Skipping BLAST step.")
      system(paste0("rm -rf ", spades.dir))
      temp.row = data.frame(Sample = sample.names[i], MappedReads = n.mapped,
                            ContigLength = NA_integer_, BestMatch = "SPAdes-failed",
                            Pident = NA_real_, AlignLength = NA_integer_,
                            Evalue = NA_real_, Bitscore = NA_real_,
                            stringsAsFactors = FALSE)
      out.csv = "logs/barcodeSampleScan_summary.csv"
      if (file.exists(out.csv)) {
        existing = read.csv(out.csv, stringsAsFactors = FALSE)
        existing = existing[!existing$Sample %in% temp.row$Sample, ]
        temp.row = rbind(existing, temp.row)
      }
      write.csv(temp.row, file = out.csv, row.names = FALSE)
      next
    }

    # Save contigs and clean up SPAdes working directory
    out.contigs = paste0(out.path, "/", sample.names[i], "_barcode-contigs.fa")
    system(paste0("cp ", contigs.file, " ", out.contigs))
    system(paste0("rm -rf ", spades.dir))

    #################################################
    ### Part D: BLAST contigs against database
    #################################################
    blast.out = paste0(out.path, "/", sample.names[i], "_blast-results.txt")

    if (use.remote.blast) {
      # Remote NCBI nt query — -num_threads is ignored remotely so omitted
      cat(" BLASTing", sample.names[i], "against NCBI nt (remote)...\n")
      system(paste0(blast.path, "blastn -task dc-megablast -remote",
                    " -db nt",
                    " -query ", out.contigs,
                    " -out ", blast.out,
                    " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\"",
                    " -max_target_seqs ", hits.per.sample),
             ignore.stdout = quiet, ignore.stderr = quiet)
    } else {
      # Local database query
      system(paste0(blast.path, "blastn -task dc-megablast",
                    " -db ", db.path, "/barcode_db",
                    " -query ", out.contigs,
                    " -out ", blast.out,
                    " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\"",
                    " -max_target_seqs ", hits.per.sample,
                    " -num_threads ", threads),
             ignore.stdout = quiet, ignore.stderr = quiet)
    }

    if (!file.exists(blast.out) || file.info(blast.out)$size == 0) {
      warning(sample.names[i], " no BLAST hits found against barcode database.")
      temp.row = data.frame(Sample = sample.names[i], MappedReads = n.mapped,
                            ContigLength = NA_integer_, BestMatch = "No-hit",
                            Pident = NA_real_, AlignLength = NA_integer_,
                            Evalue = NA_real_, Bitscore = NA_real_,
                            stringsAsFactors = FALSE)
    } else {
      match.data = read.table(blast.out, sep = "\t", header = FALSE,
                              stringsAsFactors = FALSE)
      colnames(match.data) = blast.headers

      # Best hit by bitscore
      best = match.data[which.max(match.data$bitscore), ]
      temp.row = data.frame(Sample = sample.names[i],
                            MappedReads = n.mapped,
                            ContigLength = best$qLen,
                            BestMatch = best$tName,
                            Pident = round(best$pident, 2),
                            AlignLength = best$matches,
                            Evalue = best$evalue,
                            Bitscore = best$bitscore,
                            stringsAsFactors = FALSE)
    }

    #################################################
    ### Part E: append to growing summary CSV
    #################################################
    out.csv = "logs/barcodeSampleScan_summary.csv"
    if (file.exists(out.csv)) {
      existing = read.csv(out.csv, stringsAsFactors = FALSE)
      existing = existing[!existing$Sample %in% temp.row$Sample, ]
      temp.row = rbind(existing, temp.row)
    }
    write.csv(temp.row, file = out.csv, row.names = FALSE)

    print(paste0(sample.names[i], " barcode scan complete!"))
  }#end i loop

}#end function
