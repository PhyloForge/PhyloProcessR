#' @title runCap3
#'
#' @description Assembles or merges a set of DNA contigs using CAP3. Input can
#'   be either a path to an existing FASTA file or a DNAStringSet object (which
#'   is written to a temporary file). Assembled contigs and optionally singlet
#'   sequences are concatenated into a single output file. The function also
#'   supports reading the assembled contigs back into R as a DNAStringSet.
#'
#' @param contigs either a character string giving the path to an input FASTA
#'   file, or a DNAStringSet of sequences to assemble.
#'
#' @param output.name output file name for the merged contigs (assembled +
#'   singlets). If NULL, no output file is saved.
#'
#' @param cap3.path system path to the directory containing the cap3 executable;
#'   NULL searches the system PATH.
#'
#' @param read.R logical; if TRUE the assembled contigs are read back into R and
#'   returned as a DNAStringSet.
#'
#' @param include.singlets logical; if TRUE and read.R is TRUE, unassembled
#'   singlet sequences are appended to the returned DNAStringSet.
#'
#' @param a,b,c,d,e,f,g,h,i,j,k,m,n,o,p,r,s,t,u,v,y,z numeric CAP3 algorithm
#'   parameters passed directly to the command line. See the CAP3 manual for
#'   definitions (defaults match CAP3 built-in defaults). Key parameters: o is
#'   the minimum overlap length (default 40), p is the minimum percent identity
#'   (default 90), s is the overlap similarity score cutoff (default 900), and
#'   z is the minimum number of good reads at a clipping position (default 3).
#'
#' @return if read.R is TRUE, a DNAStringSet of assembled (and optionally
#'   singlet) sequences; otherwise the function returns invisibly after writing
#'   output files.
#'
#' @export

runCap3 = function(contigs = input.contigs,
                   output.name = NULL,
                   cap3.path = NULL,
                   read.R = FALSE,
                   include.singlets = FALSE,
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
  # contigs = spades.contigs
  #  cap3.path = "cap3"
  # z = 3
  #  o = 40
  # e = 30
  #  s = 900
  #
  if (is.null(cap3.path) == FALSE){
    b.string = unlist(strsplit(cap3.path, ""))
    if (b.string[length(b.string)] != "/") {
      cap3.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { cap3.path = "" }

  if (class(contigs) != "character") {
    write.loci = as.list(as.character(contigs))
    PhyloCap::writeFasta(sequences = write.loci, names = names(write.loci),
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
    temp.assembled = Biostrings::readDNAStringSet(paste0("input_contigs.fa.cap.contigs"))
    if (include.singlets == TRUE){
      temp.singlets = Biostrings::readDNAStringSet(paste0("input_contigs.fa.cap.singlets"))
      final.save = append(temp.assembled, temp.singlets)
    } else{ final.save = temp.assembled }

      system(paste("rm input_contigs.fa*"))
    return(final.save)
  }#end if

  system(paste0("rm ", contig.file, "*"))

}#end function
