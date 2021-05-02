#' @title runCap3
#'
#' @description Function for running the program spades to assemble short read sequencing data
#'
#' @param contigs path to a folder of sequence alignments in phylip format.
#'
#' @param output.name contigs are added into existing alignment if algorithm is "add"
#'
#' @param cap3.path contigs are added into existing alignment if algorithm is "add"
#'
#' @param z algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param o TRUE applies the adjust sequence direction function of MAFFT
#'
#' @param e if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param s if a file name is provided, save.name will be used to save aligment to file as a fasta
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

runCap3 = function(contigs = input.contigs,
                   output.name = NULL,
                   cap3.path = NULL,
                   read.R = FALSE,
                   a = 20,
                   b = 20,
                   c = 12,
                   d = 200,
                   e = 30,
                   f = 20,
                   g = 6,
                   h = 20,
                   i = 40,
                   j = 80,
                   k = 1,
                   m = 2,
                   n = -5,
                   o = 40,
                   p = 90,
                   r = 1,
                   s = 900,
                   t = 300,
                   u = 3,
                   v = 2,
                   y = 100,
                   z = 3) {

  # If the file of reads is named 'xyz', then
  # the file of quality values must be named 'xyz.qual',
  # and the file of constraints named 'xyz.con'.
  # Options (default values):
  # -a  N  specify band expansion size N > 10 (20)
  # -b  N  specify base quality cutoff for differences N > 15 (20)
  # -c  N  specify base quality cutoff for clipping N > 5 (12)
  # -d  N  specify max qscore sum at differences N > 20 (200)
  # -e  N  specify clearance between no. of diff N > 10 (30)
  # -f  N  specify max gap length in any overlap N > 1 (20)
  # -g  N  specify gap penalty factor N > 0 (6)
  # -h  N  specify max overhang percent length N > 2 (20)
  # -i  N  specify segment pair score cutoff N > 20 (40)
  # -j  N  specify chain score cutoff N > 30 (80)
  # -k  N  specify end clipping flag N >= 0 (1)
  # -m  N  specify match score factor N > 0 (2)
  # -n  N  specify mismatch score factor N < 0 (-5)
  # -o  N  specify overlap length cutoff > 15 (40)
  # -p  N  specify overlap percent identity cutoff N > 65 (90)
  # -r  N  specify reverse orientation value N >= 0 (1)
  # -s  N  specify overlap similarity score cutoff N > 250 (900)
  # -t  N  specify max number of word matches N > 30 (300)
  # -u  N  specify min number of constraints for correction N > 0 (3)
  # -v  N  specify min number of constraints for linking N > 0 (2)
  # -w  N  specify file name for clipping information (none)
  # -x  N  specify prefix string for output file names (cap)
  # -y  N  specify clipping range N > 5 (100)
  # -z  N  specify min no. of good reads at clip pos N > 0 (3)

  #Writes contigs for cap3

  #Same adds to bbmap path
  if (is.null(cap3.path) == FALSE){
    b.string = unlist(strsplit(cap3.path, ""))
    if (b.string[length(b.string)] != "/") {
      cap3.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { cap3.path = "" }


  if (class(contigs) != "character") {
    write.loci = as.list(as.character(contigs))
    writeFasta(sequences = write.loci, names = names(write.loci),
               "input_contigs.fa", nbchar = 1000000, as.string = T)
    contig.file = "input_contigs.fa"
  } else { contig.file = contigs }

  #runs cap3 to merge similar contigs (pull only clustered contigs out?)
  system(paste0(cap3.path, "cap3 ", contig.file, " -z ", z, " -o ", o, " -e ", e, " -s ", s, " > ",
                "input_contigs.fa.cap.txt"))

  #Get cap3 files and deletes
  if (is.null(output.name) == F){
    system(paste0("cat ", contig.file, ".cap.contigs ",
                  contig.file, ".cap.singlets > ", output.name))
    system(paste0("rm ", contig.file, "*"))
  }#end output

  #Reads in results files
  if (read.R == TRUE){
    temp.assembled = Rsamtools::scanFa(Rsamtools::FaFile(paste0("input_contigs.fa.cap.contigs")))
    temp.singlets = Rsamtools::scanFa(Rsamtools::FaFile(paste0("input_contigs.fa.cap.singlets")))
    final.save = append(temp.assembled, temp.singlets)
    system(paste("rm input_contigs.fa*"))
    return(final.save)
  }#end if

}#end function
