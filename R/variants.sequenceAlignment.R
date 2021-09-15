#' @title variants.haplotypeCallerGATK4
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

# variants.haplotypeCallerGATK4 = function(bam.directory = NULL,
#                                          output.directory = "variant-discovery",
#                                          reference.path = "ref-index",
#                                          samtools.path = NULL,
#                                          bwa.path = NULL,
#                                          picard.path = NULL,
#                                          gatk4.path = NULL,
#                                          threads = 1,
#                                          memory = 1,
#                                          resume = TRUE,
#                                          overwrite = TRUE,
#                                          quiet = TRUE) {
#
#   #Debugging
#   #Home directoroies
#   library(doParallel)
#   work.dir = "/Volumes/Armored/Test/variant-calling" #Your main project directory
#   dir.create(work.dir)
#   setwd(work.dir)
#
#   bam.directory = "/Volumes/Armored/Test/variant-calling/variant-discovery"
#   reference.path = "ref-index"
#   subreference.name = "rag1"
#   output.directory = "variant-discovery"
#   auto.readgroup = T #Keep it T unless it crashes.
#   threads = 4 #number of threads, 8-10 is a recommended amount
#   memory = 8
#   samtools.path = "/Users/chutter/miniconda3/bin"
#   bwa.path = "/usr/local/bin"
#   picard.path = "/Users/chutter/miniconda3/bin"
#   gatk4.path = "/Users/chutter/miniconda3/bin"
#   resume = FALSE
#   quiet = FALSE
#   overwrite = TRUE
#
#   require(doParallel)
#
#
# ###############################################################################
# ######## DO NOT EDIT BELOW THIS POINT  ########################################
# ###############################################################################
#
# ############################################################################################
# ########### Step 1 #########################################################################
# ##### Merge all these samfiles and stuff
# ############################################################################################
#
# #Sets directory and reads in sample file
# setwd(work.dir)
#
# #Checks for missing files if completed
# raw.dirs = list.dirs(path = paste0(proc.dir, "/."), full.names = F, recursive = F)
# merge.dirs = list.dirs(path = paste0(proc.dir, "_MERGE/."), full.names = F, recursive = F)
# merge.names = gsub("_MERGE$", "", merge.dirs)
#
# #Gets all files
# all.files = list.files(proc.dir, full.names = F, recursive = T)
#
# #Sets up multiprocessing
# cl = makeCluster(threads)
# registerDoParallel(cl)
# mem.val = floor(mem.val/threads)
#
# #Loops through each locus and does operations on them
# foreach(i=1:24, .packages = c("foreach", "vcfR", "seqinr","data.table", "Rsamtools")) %dopar% {
#
# #Loops through each sample to combine the VCFs
# #for (i in 1:length(merge.names)){
#
#   #Goes back to home
#   setwd(paste0(proc.dir))
#
#   #Prepares and finds files for merging
#   sample.files = all.files[grep(merge.names[i], all.files)]
#   snp.files = sample.files[grep("snp_ready.bam", sample.files)]
#   snp.file.paths = paste0(proc.dir, "/", snp.files)
#
#   #Prepares and finds files for merging
#   ref.files = sample.files[grep("reference", sample.files)]
#   ref.files.all = ref.files[grep(gsub("/.*","", ref.files)[1], ref.files)]
#   ref.file.paths = paste0(proc.dir, "/", ref.files.all)
#
#   #Goes to directory
#   sample.merge = merge.dirs[grep(merge.names[i], merge.dirs)]
#   setwd(paste0(proc.dir, "_MERGE/", sample.merge) )
#
#   #Creates the new snp-analysis folder
#   dir.create("haplocontigs")
#   setwd(paste0("haplocontigs"))
#
#   if (file.exists("varscan_vars.vcf.gz") == T){ next }
#   del.files = list.files(".")
#   unlink(del.files)
#
#   #Copy files from other directories
#   for (j in 1:length(snp.file.paths)){
#     system(paste0("cp ", snp.file.paths[j], " ", as.character(getwd()), "/snp_ready_", j, ".bam" ) )
#   } #J loop end
#
#   #Next combine .bam files together!
#   snp.merge = list.files(".")
#
#   system(paste0(picard.path,
#                 " MergeSamFiles ", paste0("I=", snp.merge, collapse = " "),
#                 " O=merged_snp.bam use_jdk_deflater=true"))
#
#   #Sort by coordinate for input into MarkDuplicates  use_jdk_deflater=true
#   system(paste0(picard.path,
#                 " SortSam INPUT=merged_snp.bam OUTPUT=cleaned_final_sort.bam",
#                 " CREATE_INDEX=true SORT_ORDER=coordinate use_jdk_deflater=true"))
#
#   #Marks duplicate reads  use_jdk_deflater=true
#   system(paste0(picard.path,
#                 " MarkDuplicates INPUT=cleaned_final_sort.bam OUTPUT=cleaned_final_md.bam",
#                 " CREATE_INDEX=true METRICS_FILE=duplicate_metrics.txt use_jdk_deflater=true"))
#
#   #### Merges all the reference files
#   system(paste0("cp ", paste0(ref.file.paths, collapse = " "), " ",
#          as.character(getwd())) )
#
#   #Sorts and stuff  use_jdk_deflater=true
#   system(paste0(picard.path,
#                 " SetNmAndUqTags INPUT=cleaned_final_md.bam OUTPUT=snp_ready.bam CREATE_INDEX=true R=reference.fa use_jdk_deflater=true"))
#
#   #Delete old files to make more space  use_jdk_deflater=true
#   system("rm snp_ready_* merged_snp.bam cleaned_final* cleaned_final_sort* cleaned_final_md*")
#
#   # #Starts to finally look for Haplotypes!  use_jdk_deflater=true
#   # system(paste0("/usr/local/bin/gatk4/gatk --java-options '-Xmx", mem.val, "G' HaplotypeCaller",
#   #                " -R reference.fa -O gatk_output.g.vcf.gz -I snp_ready.bam",
#   #                " -ERC GVCF --max-alternate-alleles 3 -bamout snp_call.bam"))
#   #
#   # #### Varscan alternative
#   system(paste0("samtools mpileup -B -A -f reference.fa snp_ready.bam |",
#                 " java -jar /usr/local/bin/varscan.jar mpileup2cns --min-coverage 4",
#                 " -output-vcf 1 >varscan_output.vcf"))
#
#   system(paste0("gzip varscan_vars.vcf"))
#
#   ###############################################################################
#   ######################Step 1 pull out consensus loci  #########################
#   ###############################################################################
#
#   #Gets sample data
#   vcf.file = paste0(sample.datasets[i], "/haplocontigs/varscan_vars.vcf.gz")
#   sample.name = gsub(".*/", "", sample.datasets[i])
#   sample.name = gsub("_MERGE", "", sample.name)
#
#   #Get loci names
#   ref.file = paste0(sample.datasets[i], "/haplocontigs/reference.fa")
#   sample.exons = scanFa(file = ref.file)
#   loci.names = names(sample.exons)
#
#   #Loads in VCF data
#   vcf.locus = tryCatch(expr = {read.vcfR(vcf.file, verbose = T) },
#                        error = function(x) { x<-data.frame(); return(x) })
#
#   #Loops through each locus and does operations on them
#   save.seq = list()
#   #foreach(j=1:length(loci.names), .packages = c("foreach")) %dopar% {
#   vcf.chroms = data.table(vcf.locus@fix)
#
#   for (j in 1:length(loci.names)){
#
#     #temp.ref = sample.exons[names(sample.exons) %in% loci.names[j]]
#     #chrom = create.chromR(name=loci.names[j], vcf=vcf.locus, verbose = F, seq = temp.ref)
#
#     locus.data = vcf.chroms[vcf.chroms$CHROM == loci.names[j],]
#     locus.temp = paste0(locus.data$REF, collapse = "")
#
#     #Converts and saves sequence
#     final.locus = as.list(as.character(locus.temp))
#     names(final.locus) = paste0(loci.names[j], "_|_", sample.name)
#     save.seq = append(save.seq, final.locus)
#
#   }#end j loop
#
#   # writes the data
#   seqinr::write.fasta(sequences = save.seq, names = names(save.seq),
#               paste(sample.name, "_consensus-haplocontigs.fa", sep = ""), nbchar = 1000000, as.string = T)
#
# } # end j loop
#
# #END SCRIPT

