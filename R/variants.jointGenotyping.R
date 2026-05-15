#' @title jointGenotyping
#'
#' @description Performs joint genotyping across all samples at each locus
#'   using GATK4. For each locus, all per-sample GVCFs are combined into a
#'   GenomicsDB, then GenotypeGVCFs genotypes the combined data. SNPs and
#'   indels are separated, hard-filtered with user-specified thresholds, merged
#'   back, and written to per-locus VCF files in subdirectories. Loci are
#'   processed in parallel. The function optionally removes unfiltered and/or
#'   SNP/indel-specific VCF directories at the end.
#'
#' @param haplotype.caller.directory path to the directory of per-sample
#'   haplotype caller GVCF files (output of haplotypeCaller()).
#'
#' @param output.directory path where the genotype database and per-locus VCF
#'   subdirectories (filtered-all, filtered-snps, filtered-indels, and their
#'   unfiltered counterparts) will be created.
#'
#' @param use.base.recalibration logical; if TRUE the BQSR-recalibrated GVCFs
#'   are used as input.
#'
#' @param save.unfiltered logical; if FALSE unfiltered VCF subdirectories are
#'   deleted after filtering.
#' @param save.SNPs logical; if FALSE SNP-specific VCF subdirectories are
#'   deleted.
#' @param save.indels logical; if FALSE indel-specific VCF subdirectories are
#'   deleted.
#' @param save.combined logical; if FALSE the combined (all variants) VCF
#'   subdirectories are deleted.
#'
#' @param custom.SNP.QD,custom.SNP.QUAL,custom.SNP.SOR,custom.SNP.FS,custom.SNP.MQ,custom.SNP.MQRankSum,custom.SNP.ReadPosRankSum
#'   numeric hard-filter thresholds for SNPs (see GATK VariantFiltration
#'   documentation for field definitions).
#' @param custom.INDEL.QD,custom.INDEL.QUAL,custom.INDEL.FS,custom.INDEL.ReadPosRankSum
#'   numeric hard-filter thresholds for indels.
#'
#' @param gatk4.path system path to the directory containing the gatk
#'   executable; NULL searches the system PATH.
#'
#' @param temp.directory path to a GATK JVM temp directory; NULL uses the
#'   current working directory.
#'
#' @param threads number of loci to process in parallel.
#'
#' @param memory total RAM in GB to allocate as the JVM heap (-Xmx).
#'
#' @param overwrite logical; if TRUE the output directory is deleted and
#'   recreated before processing.
#'
#' @param quiet logical; currently unused.
#'
#' @return invisibly; writes per-locus VCF files to output.directory
#'   subdirectories.
#'
#' @export

jointGenotyping = function(haplotype.caller.directory = "haplotype-caller",
                          output.directory = "genotype-database",
                          use.base.recalibration = FALSE,
                          save.unfiltered = TRUE,
                          save.SNPs = TRUE,
                          save.indels = TRUE,
                          save.combined = TRUE,
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
                          overwrite = FALSE,
                          quiet = TRUE) {

 #Debugging
  # #Home directoroies
  # library(PhyloProcessR)
  # setwd("/Volumes/LaCie/Anax")

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
  # save.unfiltered = TRUE
  # save.SNPs = TRUE
  # save.indels = TRUE
  # save.combined = TRUE

  # dataset.name = "joint-genotyping"
  # haplotype.caller.directory = paste0("data-analysis/", dataset.name, "/haplotype-caller")
  # output.directory = paste0("data-analysis/", dataset.name, "/genotype-database")

  # gatk4.path <- "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  # use.base.recalibration = FALSE
  # threads <- 4
  # memory <- 8
  # quiet <- FALSE
  # overwrite <- FALSE

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

  if (is.null(temp.directory) == TRUE){
    temp.directory = getwd()
  }


  #Creates output directory
  if (dir.exists("logs/sample_logs") == F){ dir.create("logs/sample_logs", recursive = TRUE) }

  # Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.directory) == FALSE) {
    dir.create(output.directory)
    dir.create(paste0(output.directory, "/unfiltered-all"))
    dir.create(paste0(output.directory, "/unfiltered-snps"))
    dir.create(paste0(output.directory, "/unfiltered-indels"))
    dir.create(paste0(output.directory, "/filtered-all"))
    dir.create(paste0(output.directory, "/filtered-snps"))
    dir.create(paste0(output.directory, "/filtered-indels"))

  } else {
    if (overwrite == TRUE) {
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
      dir.create(paste0(output.directory, "/unfiltered-all"))
      dir.create(paste0(output.directory, "/unfiltered-snps"))
      dir.create(paste0(output.directory, "/unfiltered-indels"))
      dir.create(paste0(output.directory, "/filtered-all"))
      dir.create(paste0(output.directory, "/filtered-snps"))
      dir.create(paste0(output.directory, "/filtered-indels"))
    }
  } # end else

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

  # Get loci names I guess
  reference.path = "index/reference.fa"
  reference.seq = Biostrings::readDNAStringSet((paste0("index/reference.fa")))
  loci.names = names(reference.seq)

  # Resumes file download
  if (overwrite == FALSE) {
    done.files <- list.files(paste0(output.directory, "/filtered-all"), full.names = T, recursive = T)
    done.names <- gsub(".vcf$", "", done.files)
    done.names <- gsub(".*\\/", "", done.names)
    loci.names <- loci.names[!loci.names %in% done.names]
  }

  # Get multifile databases together
  sample.names <- list.files(haplotype.caller.directory, recursive = TRUE, full.names = TRUE)

  if (length(sample.names) == 0) {
    return("no samples available to analyze.")
  }

  if (use.base.recalibration == TRUE) {
    sample.names = sample.names[grep("gatk4-bqsr-haplotype-caller.g.vcf.gz$", sample.names)]
    if (length(sample.names) == 0) {
      stop("the haplotype caller file for use.base.recalibration does not exist.")
    }
  }

  if (use.base.recalibration == FALSE) {
    sample.names = sample.names[grep("gatk4-haplotype-caller.g.vcf.gz$", sample.names)]
    if (length(sample.names) == 0) {
      stop("the haplotype caller file does not exist.")
    }
  }

  #Gathers input file names together
  sample.datasets = paste0("-V ", sample.names)
  vcf.files = paste0(sample.datasets, collapse = " ")

  mem.cl = floor(memory / threads)

  # Loops through each locus and does operations on them
  parallel::mclapply(seq_along(loci.names), function(i) {
  tryCatch({

    # Genotype haplotype caller results
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", mem.cl, "G\"",
      " GenomicsDBImport ", vcf.files,
      " --genomicsdb-workspace-path ", output.directory, "/", loci.names[i],
      " --intervals ", loci.names[i]
    ))

    # Genotype haplotype caller results
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", mem.cl, "G\"",
      " GenotypeGVCFs -R ", reference.path,
      " -V gendb://", output.directory, "/", loci.names[i],
      " --use-new-qual-calculator true",
      " -O ", output.directory, "/unfiltered-all/", loci.names[i], ".vcf"
    ))

    # Selects only the SNPs from the VCF
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", mem.cl, "G\"",
      " SelectVariants",
      " -V ", output.directory, "/unfiltered-all/", loci.names[i], ".vcf",
      " -O ", output.directory, "/unfiltered-snps/", loci.names[i], ".vcf",
      " --select-type SNP"
    ))

    # Selects only the indels from the VCF
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", mem.cl, "G\"",
      " SelectVariants",
      " -V ", output.directory, "/unfiltered-all/", loci.names[i], ".vcf",
      " -O ", output.directory, "/unfiltered-indels/", loci.names[i], ".vcf",
      " --select-type INDEL"
    ))

    ###########################################################################################
    ######### Filtering begineth
    ###########################################################################################

    # Custom filtering
    #########################
    # Applies filters to SNPs
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", mem.cl, "G\"",
      " VariantFiltration -R ", reference.path,
      " -V ", output.directory, "/unfiltered-snps/", loci.names[i], ".vcf",
      " -O ", output.directory, "/unfiltered-snps/", loci.names[i], "_filter.vcf",
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
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", mem.cl, "G\"",
      " VariantFiltration -R ", reference.path,
      " -V ", output.directory, "/unfiltered-indels/", loci.names[i], ".vcf",
      " -O ", output.directory, "/unfiltered-indels/", loci.names[i], "_filter.vcf",
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
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", mem.cl, "G\"",
      " SortVcf",
      " -I ", output.directory, "/unfiltered-snps/", loci.names[i], "_filter.vcf",
      " -I ", output.directory, "/unfiltered-indels/", loci.names[i], "_filter.vcf",
      " -O ", output.directory, "/unfiltered-all/", loci.names[i], "_filter.vcf"
    ))

    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", mem.cl, "G\"",
      " SelectVariants",
      " -V ", output.directory, "/unfiltered-all/", loci.names[i], "_filter.vcf",
      " -O ", output.directory, "/filtered-all/", loci.names[i], ".vcf",
      " --exclude-filtered TRUE"
    ))

    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", mem.cl, "G\"",
      " SelectVariants",
      " -V ", output.directory, "/unfiltered-snps/", loci.names[i], "_filter.vcf",
      " -O ", output.directory, "/filtered-snps/", loci.names[i], ".vcf",
      " --exclude-filtered TRUE"
    ))

    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", mem.cl, "G\"",
      " SelectVariants",
      " -V ", output.directory, "/unfiltered-indels/", loci.names[i], "_filter.vcf",
      " -O ", output.directory, "/filtered-indels/", loci.names[i], ".vcf",
      " --exclude-filtered TRUE"
    ))

    print(paste0(loci.names[i], " completed GATK4 joint sample genotyping!"))
    system(paste0("rm -r ", output.directory, "/", loci.names[i]))

  }, error = function(e) {
    warning(loci.names[i], " failed: ", conditionMessage(e))
  })
  }, mc.cores = threads) #end i loop

  #Datasets to save
  #Save SNPs
  if (save.SNPs == FALSE) {
    system(paste0("rm -r ", output.directory, "/unfiltered-snps ", output.directory, "/filtered-snps"))
  }

  #Save indels
  if (save.indels == FALSE) {
    system(paste0("rm -r ", output.directory, "/unfiltered-indels ", output.directory, "/filtered-indels"))
  }

  #Save all
  if (save.combined == FALSE) {
    system(paste0("rm -r ", output.directory, "/unfiltered-all ", output.directory, "/filtered-all"))
  }

  #Save unfiltered
  if (save.unfiltered == FALSE) {
    system(paste0("rm -r ", output.directory, "/unfiltered*"))
  }

}#end function

# END SCRIPT
