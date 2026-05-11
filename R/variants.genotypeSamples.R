#' @title genotypeSamples
#'
#' @description Performs per-sample genotyping and variant filtering using
#'   GATK4. For each sample, GenotypeGVCFs converts a haplotype caller GVCF
#'   to a genotyped VCF, SelectVariants separates SNPs and indels, and
#'   VariantFiltration applies user-specified hard filters. Final filtered VCFs
#'   for SNPs, indels, and their combination are saved per sample. Samples are
#'   processed in parallel.
#'
#' @param mapping.directory path to the directory containing per-sample
#'   reference index sub-directories.
#'
#' @param haplotype.caller.directory path to the directory of per-sample
#'   haplotype caller GVCF files (output of haplotypeCaller()).
#'
#' @param output.directory path where per-sample genotype VCF sub-directories
#'   will be saved.
#'
#' @param use.base.recalibration logical; if TRUE the BQSR-recalibrated GVCF
#'   (gatk4-bqsr-haplotype-caller.g.vcf.gz) is used as input; if FALSE the
#'   standard GVCF (gatk4-haplotype-caller.g.vcf.gz) is used.
#'
#' @param custom.SNP.QD numeric QD filter threshold for SNPs (variant is
#'   filtered if QD < this value).
#' @param custom.SNP.QUAL numeric QUAL filter threshold for SNPs.
#' @param custom.SNP.SOR numeric SOR filter threshold for SNPs (filtered if
#'   SOR > this value).
#' @param custom.SNP.FS numeric FS filter threshold for SNPs (filtered if
#'   FS > this value).
#' @param custom.SNP.MQ numeric MQ filter threshold for SNPs (filtered if
#'   MQ < this value).
#' @param custom.SNP.MQRankSum numeric MQRankSum filter threshold for SNPs
#'   (filtered if MQRankSum < this value).
#' @param custom.SNP.ReadPosRankSum numeric ReadPosRankSum filter threshold for
#'   SNPs (filtered if ReadPosRankSum < this value).
#' @param custom.INDEL.QD numeric QD filter threshold for indels.
#' @param custom.INDEL.QUAL numeric QUAL filter threshold for indels.
#' @param custom.INDEL.FS numeric FS filter threshold for indels.
#' @param custom.INDEL.ReadPosRankSum numeric ReadPosRankSum filter threshold
#'   for indels.
#'
#' @param gatk4.path system path to the directory containing the gatk
#'   executable; NULL searches the system PATH.
#'
#' @param temp.directory path to a temporary directory for GATK JVM temp files;
#'   NULL uses the current working directory.
#'
#' @param threads number of parallel samples to process simultaneously.
#'
#' @param memory total RAM in GB to allocate as the JVM heap (-Xmx).
#'
#' @param overwrite logical; if TRUE the output directory is deleted and
#'   recreated before processing.
#'
#' @param quiet logical; currently unused.
#'
#' @return invisibly; writes filtered VCF files to per-sample sub-directories
#'   in output.directory.
#'
#' @export

genotypeSamples = function(mapping.directory = "sample-mapping",
                            haplotype.caller.directory = "haplotype-caller",
                            output.directory = "sample-genotypes",
                            use.base.recalibration = FALSE,
                            custom.SNP.QD =  NULL,
                            custom.SNP.QUAL =  NULL,
                            custom.SNP.SOR =  NULL,
                            custom.SNP.FS =  NULL,
                            custom.SNP.MQ =  NULL,
                            custom.SNP.MQRankSum =  NULL,
                            custom.SNP.ReadPosRankSum =  NULL,
                            custom.INDEL.QD =  NULL,
                            custom.INDEL.QUAL =  NULL,
                            custom.INDEL.FS =  NULL,
                            custom.INDEL.ReadPosRankSum =  NULL,
                            gatk4.path = NULL,
                            temp.directory = NULL,
                            threads = 1,
                            memory = 1,
                            overwrite = TRUE,
                            quiet = TRUE) {

 #Debugging
  #Home directoroies
  # library(PhyloProcessR)
  # setwd("/Volumes/LaCie/Mantellidae")

  # custom.SNP.QD <- 2
  # custom.SNP.QUAL <- 30
  # custom.SNP.SOR <- 3
  # custom.SNP.FS <- 60
  # custom.SNP.MQ <- 40
  # custom.SNP.MQRankSum <- -12.5
  # custom.SNP.ReadPosRankSum <- -8
  # custom.INDEL.QD <- 2
  # custom.INDEL.QUAL <- 30
  # custom.INDEL.FS <- 60
  # custom.INDEL.ReadPosRankSum <- -8


  # dataset.name = "variant-calling"
  # mapping.directory = paste0("data-analysis/", dataset.name, "/sample-mapping")
  # haplotype.caller.directory = paste0("data-analysis/", dataset.name, "/haplotype-caller")
  # output.directory = paste0("data-analysis/", dataset.name, "/sample-genotypes")

  # output.directory <- "sample-genotypes"
  # gatk4.path <- "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  # use.base.recalibration = FALSE
  # filtering.thresholds = "high"
  # threads <- 4
  # memory <- 8
  # quiet <- FALSE
  # overwrite <- TRUE

  # Same adds to bbmap path
  require(foreach)

  # Same adds to bbmap path
  if (is.null(gatk4.path) == FALSE) {
    b.string <- unlist(strsplit(gatk4.path, ""))
    if (b.string[length(b.string)] != "/") {
      gatk4.path <- paste0(append(b.string, "/"), collapse = "")
    } # end if
  } else {
    gatk4.path <- ""
  }

  #Quick checks
  if (is.null(haplotype.caller.directory) == TRUE) {
    stop("Please provide the haplotype caller directory.")
  }
  if (file.exists(haplotype.caller.directory) == F) {
    stop("Haplotype caller directory not found.")
  }

  # Quick checks
  if (is.null(mapping.directory) == TRUE) {
    stop("Please provide the sample mapping directory.")
  }
  if (file.exists(mapping.directory) == FALSE) {
    stop("Sample mapping directory not found.")
  }

  if (is.null(temp.directory) == TRUE){
    temp.directory = getwd()
  }


  #Creates output directory
  if (dir.exists("logs") == F){ dir.create("logs") }

  # Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.directory) == FALSE) {
    dir.create(output.directory)
  } else {
    if (overwrite == TRUE) {
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } # end else

  # Get multifile databases together
  sample.names <- list.dirs(haplotype.caller.directory, recursive = F, full.names = F)

  # Resumes file download
  if (overwrite == FALSE) {
    done.files <- list.files(output.directory, full.names = T, recursive = T)
    done.files <- done.files[grep("gatk4-final-genotypes.vcf$", done.files)]
    done.names <- gsub("/gatk4-final-genotypes.vcf$", "", done.files)
    done.names <- gsub(".*\\/", "", done.names)
    sample.names <- sample.names[!sample.names %in% done.names]
  }

  if (length(sample.names) == 0) {
    return("no samples available to analyze.")
  }

  if (length(sample.names) == 0){ return("no samples available to analyze.") }


  SNP.QD.string <- paste0(" -filter \"QD<", format(custom.SNP.QD, nsmall = 1), "\" --filter-name \"QD\"")
  SNP.QUAL.string <- paste0(" -filter \"QUAL<", format(custom.SNP.QUAL, nsmall = 1), "\" --filter-name \"QUAL\"")
  SNP.SOR.string <- paste0(" -filter \"SOR>", format(custom.SNP.SOR, nsmall = 1), "\" --filter-name \"SOR\"")
  SNP.FS.string <- paste0(" -filter \"FS>", format(custom.SNP.FS, nsmall = 1), "\" --filter-name \"FS\"")
  SNP.MQ.string <- paste0(" -filter \"MQ<", format(custom.SNP.MQ, nsmall = 1), "\" --filter-name \"MQ\"")
  SNP.MQRankSum.string <- paste0(" -filter \"MQRankSum<", format(custom.SNP.MQRankSum, nsmall = 1), "\" --filter-name \"MQRankSum\"")
  SNP.ReadPosRankSum.string <- paste0(" -filter \"ReadPosRankSum<", format(custom.SNP.ReadPosRankSum, nsmall = 1), "\" --filter-name \"ReadPosRankSum\"")

  IN.QD.string <- paste0(" -filter \"QD<", format(custom.INDEL.QD, nsmall = 1), "\" --filter-name \"QD\"")
  IN.QUAL.string <- paste0(" -filter \"QUAL<", format(custom.INDEL.QUAL, nsmall = 1), "\" --filter-name \"QUAL\"")
  IN.FS.string <- paste0(" -filter \"FS>", format(custom.INDEL.FS, nsmall = 1), "\" --filter-name \"FS\"")
  IN.ReadPosRankSum.string <- paste0(" -filter \"ReadPosRankSum<", format(custom.INDEL.ReadPosRankSum, nsmall = 1), "\" --filter-name \"ReadPosRankSum\"")

  ############################################################################################
  ########### Step 1 #########################################################################
  ##### Start up loop for each sample
  ############################################################################################

  # Sets up multiprocessing
  cl <- snow::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  mem.cl <- floor(memory / threads)

  #Loops through each locus and does operations on them
  foreach::foreach(i = seq_along(sample.names), .packages = c("foreach")) %dopar% {
    # Loops through each locus and does operations on them
    # for (i in 1:length(loci.names)){

    dir.create(paste0(output.directory, "/", sample.names[i]))

    reference.path <- paste0(mapping.directory, "/", sample.names[i], "/index/reference.fa")

    if (use.base.recalibration == TRUE) {
      haplocaller.file <- paste0(haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-haplotype-caller.g.vcf.gz")
      if (file.exists(haplocaller.file) == FALSE) {
        stop("the haplotype caller file for use.base.recalibration does not exist.")
      }
    }

    if (use.base.recalibration == FALSE) {
      haplocaller.file <- paste0(haplotype.caller.directory, "/", sample.names[i], "/gatk4-haplotype-caller.g.vcf.gz")
      if (file.exists(haplocaller.file) == FALSE) {
        stop("the haplotype caller file does not exist.")
      }
    }

    # Genotype haplotype caller results
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " GenotypeGVCFs -R ", reference.path,
      " -V ", haplocaller.file,
      " --use-new-qual-calculator true -O ", output.directory, "/", sample.names[i], "/gatk4-unfiltered-genotypes.vcf"
    ))

    # Selects only the SNPs from the VCF
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " SelectVariants",
      " -V ", output.directory, "/", sample.names[i], "/gatk4-unfiltered-genotypes.vcf",
      " -O ", output.directory, "/", sample.names[i], "/gatk4-unfiltered-snps.vcf",
      " --select-type SNP"
    ))

    #Selects only the indels from the VCF
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " SelectVariants",
      " -V ", output.directory, "/", sample.names[i], "/gatk4-unfiltered-genotypes.vcf",
      " -O ", output.directory, "/", sample.names[i], "/gatk4-unfiltered-indels.vcf",
      " --select-type INDEL"
    ))

    ###########################################################################################
    ######### Filtering begineth
    ###########################################################################################

    # Custom filtering
    #########################
    # Applies filters to SNPs
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " VariantFiltration -R ", reference.path,
      " -V ", output.directory, "/", sample.names[i], "/gatk4-unfiltered-snps.vcf",
      " -O ", output.directory, "/", sample.names[i], "/gatk4-filtered-snps.vcf",
      SNP.QD.string,
      SNP.QUAL.string,
      SNP.SOR.string,
      SNP.FS.string,
      SNP.MQ.string,
      SNP.MQRankSum.string,
      SNP.ReadPosRankSum.string
    ))

    # Applies filters to indels
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " VariantFiltration -R ", reference.path,
      " -V ", output.directory, "/", sample.names[i], "/gatk4-unfiltered-indels.vcf",
      " -O ", output.directory, "/", sample.names[i], "/gatk4-filtered-indels.vcf",
      IN.QD.string,
      IN.QUAL.string,
      IN.FS.string,
      IN.ReadPosRankSum.string
    ))

    ###########################################################################################
    ######### Filtering finished
    ###########################################################################################

    # Combine them into a single VCF
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " SortVcf",
      " -I ", output.directory, "/", sample.names[i], "/gatk4-filtered-snps.vcf",
      " -I ", output.directory, "/", sample.names[i], "/gatk4-filtered-indels.vcf",
      " -O ", output.directory, "/", sample.names[i], "/gatk4-filtered-genotypes.vcf"
    ))

    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " SelectVariants",
      " -V ", output.directory, "/", sample.names[i], "/gatk4-filtered-genotypes.vcf",
      " -O ", output.directory, "/", sample.names[i], "/gatk4-final-genotypes.vcf",
      " --exclude-filtered TRUE"
    ))

    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " SelectVariants",
      " -V ", output.directory, "/", sample.names[i], "/gatk4-filtered-snps.vcf",
      " -O ", output.directory, "/", sample.names[i], "/gatk4-final-snps.vcf",
      " --exclude-filtered TRUE"
    ))

    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " SelectVariants",
      " -V ", output.directory, "/", sample.names[i], "/gatk4-filtered-indels.vcf",
      " -O ", output.directory, "/", sample.names[i], "/gatk4-final-indels.vcf",
      " --exclude-filtered TRUE"
    ))

    print(paste0(sample.names[i], " completed GATK4 sample genotyping!"))

  }#end i loop

  snow::stopCluster(cl)

}#end function

# END SCRIPT
