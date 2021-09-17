#' @title variants.sequenceAlignment
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

variants.sequenceAlignment = function(database.directory = NULL,
                                      reference.path = "ref-index",
                                      output.directory = "variant-calling/alignments",
                                      invariant.sites = TRUE,
                                      ambiguity.codes = FALSE,
                                      haplotype.split = FALSE,
                                      save.all.types = FALSE,
                                      threads = 1,
                                      memory = 1,
                                      resume = TRUE,
                                      overwrite = TRUE,
                                      quiet = TRUE) {

  #Debugging
  #Home directoroies
  # setwd("/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/Microhylidae_SeqCap/New_Work_2021/Extract_Haplotypes")
  # database.directory = "variant-calling/variant-database"
  # reference.path = "ref-index"
  # #subreference.name = "rag1"
  # output.directory = "variant-calling/alignments"
  # threads = 4 #number of threads, 8-10 is a recommended amount
  # memory = 8
  # gatk4.path = "/Users/chutter/miniconda3/bin"
  # resume = FALSE
  # quiet = FALSE
  # overwrite = TRUE
  # invariant.sites = TRUE
  # ambiguity.codes = TRUE
  # haplotype.split = FALSE
  # save.all.types = TRUE


  #Same adds to bbmap path
  require(doParallel)
  if (is.null(gatk4.path) == FALSE){
    b.string = unlist(strsplit(gatk4.path, ""))
    if (b.string[length(b.string)] != "/") {
      gatk4.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { gatk4.path = "" }

  #Quick checks
  if (is.null(database.directory) == TRUE){ stop("Please provide the directory to the variant database.") }
  if (file.exists(database.directory) == F){ stop("Database directory not found.") }

  #Creates output directory
  if (dir.exists("logs") == F){ dir.create("logs") }
  if (dir.exists(output.directory) == F){ dir.create(output.directory) }

  #Creates various directories
  dir.create(paste0(output.directory, "/contigs-haplotype"))
  dir.create(paste0(output.directory, "/contigs-consensus"))
  dir.create(paste0(output.directory, "/contigs-ambiguity"))
  dir.create(paste0(output.directory, "/variants-haplotype"))
  dir.create(paste0(output.directory, "/variants-consensus"))
  dir.create(paste0(output.directory, "/variants-ambiguity"))

  #Get loci names I guess
  probe.loci = Biostrings::readDNAStringSet(paste0(reference.path, "/reference.fa"))
  locus.names = names(probe.loci)

  #Get loci names I guess
  loci.files = list.files(database.directory, full.names = TRUE)
  loci.files = loci.files[grep(".vcf$", loci.files)]
  file.names = gsub(".vcf$", "", loci.files)
  file.names = gsub(".*\\/", "", file.names)
  locus.names = loci.names[locus.names %in% file.names]

  if (length(locus.names) == 0){ return("no loci remain to analyze.") }

  # #Resumes file download
  # if (resume == TRUE){
  #   done.files = list.files(output.directory, full.names = T, recursive = T)
  #   done.files = done.files[grep("gatk4-haplotype-caller.g.vcf.gz$", done.files)]
  #   done.names = gsub("/gatk4-haplotype-caller.g.vcf.gz$", "", done.files)
  #   done.names = gsub(".*\\/", "", done.names)
  #   sample.names = sample.names[!sample.names %in% done.names]
  # }
  #
  # #Resumes file download
  # if (overwrite == FALSE){
  #   done.files = list.files(output.directory, full.names = T, recursive = T)
  #   done.files = done.files[grep("gatk4-haplotype-caller.g.vcf.gz$", done.files)]
  #   done.names = gsub("/gatk4-haplotype-caller.g.vcf.gz$", "", done.files)
  #   done.names = gsub(".*\\/", "", done.names)
  #   sample.names = sample.names[!sample.names %in% done.names]
  # }

  ############################################################################################
  ########### Step 1 #########################################################################
  ##### Start up loop for each sample
  ############################################################################################

  #Sets up multiprocessing
  cl = makeCluster(threads)
  registerDoParallel(cl)
  mem.cl = floor(memory/threads)

  #Loops through each locus and does operations on them
  foreach(i=1:length(locus.names), .packages = c("foreach", "vcfR", "PhyloCap", "ape", "Biostrings")) %dopar% {

    #Loads in alignment data
    sub.file = loci.files[grep(locus.names[i], loci.files)]
    ref.locus = probe.loci[grep(locus.names[i], names(probe.loci))]

    #Loads in VCF data
    vcf.locus = tryCatch(expr = {vcfR::read.vcfR(sub.file, verbose = FALSE) },
                         error = function(x) { x<-data.frame(); return(x) })

    #Catches exceptions or bad files
    if (nrow(vcf.locus) == 0){
      print(paste0(vcf.files[i], " Could not be read. File is empty or corrupted."))
      next
    }

    #Also catches other types of bad files that load in
    if(nrow(vcf.locus@fix) == 0){
      print(paste0(locus.names[i], " Had no initial variants. Skip."))
      next
    }

    #Extracts genotypes
    #geno.types = extract.gt(vcf.locus)
    #Gts the data
    chrom = vcfR::create.chromR(name='Supercontig', vcf=vcf.locus)
    dp.filt = quantile(chrom@var.info$DP, c(0.01, 0.99))

    #Creates a mask
    chrom = vcfR::masker(chrom, min_QUAL=1, min_DP=dp.filt[1], max_DP=dp.filt[2])

    #makes a new contig with ambiguities
    temp.ref = strsplit(as.character(ref.locus), "")
    temp.ref = lapply(temp.ref, tolower)
    temp.ref = as.matrix(ape::as.DNAbin(temp.ref))

    #####################################################
    ###### output 1: contigs-ambiguity
    ######
    if (haplotype.split == FALSE && invariant.sites == TRUE && ambiguity.codes == TRUE | save.all.types == TRUE){

      new.align = vcfR::vcfR2DNAbin(chrom,
                                    consensus = T, #<= set to TRUE
                                    extract.haps = F,
                                    extract.indels = T,
                                    asterisk_as_del = T,
                                    unphased_as_NA = F,
                                    ref.seq = temp.ref,
                                    start.pos = 1,
                                    verbose = quiet)

      #Fixes messed up NAs
      temp.align = as.character(new.align)
      temp.align[is.na(temp.align) == T] = "n"
      temp.align = as.list(data.frame(t(temp.align)))
      temp.align2 = lapply(temp.align, FUN = function(x) paste(x, collapse = ""))
      align.out = Biostrings::DNAStringSet(unlist(temp.align2))

      if (length(align.out) == 0) {
        print(paste0(locus.names[i], " Had no variants. Skip."))
        next
      }

      #If no alignment assessing is done, saves
      write.temp = strsplit(as.character(align.out), "")
      aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
      #readies for saving
      writePhylip(aligned.set, file= paste0(output.directory, "/contigs-ambiguity/", locus.names[i], ".phy"), interleave = F)
    }#end if

    #####################################################
    ###### output 2: contigs-haplocontigs
    ######
    if (haplotype.split == TRUE && invariant.sites == TRUE | save.all.types == TRUE){

      new.align = vcfR::vcfR2DNAbin(chrom,
                                    consensus = F, #<- set to FALSE
                                    extract.haps = T, #<- set to TRUE
                                    extract.indels = T,
                                    asterisk_as_del = T,
                                    unphased_as_NA = F,
                                    ref.seq = temp.ref,
                                    start.pos = 1,
                                    verbose = quiet)

      #Fixes messed up NAs
      temp.align = as.character(new.align)
      temp.align[is.na(temp.align) == T] = "n"
      temp.align = as.list(data.frame(t(temp.align)))
      temp.align2 = lapply(temp.align, FUN = function(x) paste(x, collapse = ""))
      align.out = Biostrings::DNAStringSet(unlist(temp.align2))

      if (length(align.out) == 0) {
        print(paste0(locus.names[i], " Had no variants. Skip."))
        next
      }

      #If no alignment assessing is done, saves
      write.temp = strsplit(as.character(align.out), "")
      aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
      #readies for saving
      writePhylip(aligned.set, file= paste0(output.directory, "/contigs-haplotype/", locus.names[i], ".phy"), interleave = F)
    }#end if


    #####################################################
    ###### output 3: contigs-consensus
    ######
    if (haplotype.split == FALSE && invariant.sites == TRUE && ambiguity.codes == FALSE | save.all.types == TRUE){

      new.align = vcfR::vcfR2DNAbin(chrom,
                                    consensus = T, # <-- changed here
                                    extract.haps = F,
                                    extract.indels = T,
                                    asterisk_as_del = T,
                                    unphased_as_NA = T,
                                    ref.seq = temp.ref,
                                    start.pos = 1,
                                    verbose = quiet)

      #Fixes messed up NAs
      temp.align = as.character(new.align)
      temp.align[is.na(temp.align) == T] = "n"
      temp.align = as.list(data.frame(t(temp.align)))
      temp.align2 = lapply(temp.align, FUN = function(x) paste(x, collapse = ""))
      align.out = Biostrings::DNAStringSet(unlist(temp.align2))
      non.align = convertAmbiguousConsensus(alignment = align.out)

      if (length(non.align) == 0) {
        print(paste0(locus.names[i], " Had no variants. Skip."))
        next
      }

      #If no alignment assessing is done, saves
      write.temp = strsplit(as.character(non.align), "")
      aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
      #readies for saving
      writePhylip(aligned.set, file= paste0(output.directory, "/contigs-consensus/", locus.names[i], ".phy"), interleave = F)
    }#end if

    #####################################################
    #####################################################
    #####################################################
    #####################################################
    ###### output 4: variants-ambiguity
    ######
    if (haplotype.split == FALSE && invariant.sites == FALSE && ambiguity.codes == TRUE | save.all.types == TRUE){

      new.align = vcfR::vcfR2DNAbin(chrom,
                                    consensus = T, #<= set to TRUE
                                    extract.haps = F,
                                    extract.indels = T,
                                    asterisk_as_del = T,
                                    unphased_as_NA = F,
                                    start.pos = 1,
                                    verbose = quiet)

      #Fixes messed up NAs
      temp.align = as.character(new.align)
      temp.align[is.na(temp.align) == T] = "n"
      temp.align = as.list(data.frame(t(temp.align)))
      temp.align2 = lapply(temp.align, FUN = function(x) paste(x, collapse = ""))
      align.out = Biostrings::DNAStringSet(unlist(temp.align2))

      if (length(align.out) == 0) {
        print(paste0(locus.names[i], " Had no variants. Skip."))
        next
      }

      #If no alignment assessing is done, saves
      write.temp = strsplit(as.character(align.out), "")
      aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
      #readies for saving
      writePhylip(aligned.set, file= paste0(output.directory, "/variants-ambiguity/", locus.names[i], ".phy"), interleave = F)
    }#end if

    #####################################################
    ###### output 5: variants-haplocontigs
    ######
    if (haplotype.split == TRUE && invariant.sites == FALSE | save.all.types == TRUE){

      new.align = vcfR::vcfR2DNAbin(chrom,
                                    consensus = F, #<- set to FALSE
                                    extract.haps = T, #<- set to TRUE
                                    extract.indels = T,
                                    asterisk_as_del = T,
                                    unphased_as_NA = F,
                                    start.pos = 1,
                                    verbose = quiet)

      #Fixes messed up NAs
      temp.align = as.character(new.align)
      temp.align[is.na(temp.align) == T] = "n"
      temp.align = as.list(data.frame(t(temp.align)))
      temp.align2 = lapply(temp.align, FUN = function(x) paste(x, collapse = ""))
      align.out = Biostrings::DNAStringSet(unlist(temp.align2))

      if (length(align.out) == 0) {
        print(paste0(locus.names[i], " Had no variants. Skip."))
        next
      }

      #If no alignment assessing is done, saves
      write.temp = strsplit(as.character(align.out), "")
      aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
      #readies for saving
      writePhylip(aligned.set, file= paste0(output.directory, "/variants-haplotype/", locus.names[i], ".phy"), interleave = F)
    }#end if

    #####################################################
    ###### output 6: variants-consensus
    ######
    if (haplotype.split == FALSE && invariant.sites == FALSE && ambiguity.codes == FALSE | save.all.types == TRUE){

      new.align = vcfR::vcfR2DNAbin(chrom,
                                    consensus = T, # <-- changed here
                                    extract.haps = F,
                                    extract.indels = T,
                                    asterisk_as_del = T,
                                    unphased_as_NA = T,
                                    start.pos = 1,
                                    verbose = quiet)

      #Fixes messed up NAs
      temp.align = as.character(new.align)
      temp.align[is.na(temp.align) == T] = "n"
      temp.align = as.list(data.frame(t(temp.align)))
      temp.align2 = lapply(temp.align, FUN = function(x) paste(x, collapse = ""))
      align.out = Biostrings::DNAStringSet(unlist(temp.align2))
      non.align = convertAmbiguousConsensus(alignment = align.out)

      if (length(non.align) == 0) {
        print(paste0(locus.names[i], " Had no variants. Skip."))
        next
      }

      #If no alignment assessing is done, saves
      write.temp = strsplit(as.character(non.align), "")
      aligned.set = as.matrix(ape::as.DNAbin(write.temp) )

      #readies for saving
      writePhylip(aligned.set, file= paste0(output.directory, "/variants-consensus/", locus.names[i], ".phy"), interleave = F)
    }#end if

  }#end i loop

  stopCluster(cl)

}#end function




# #END SCRIPT

