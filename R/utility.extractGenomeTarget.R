#' @title extractGenomeTarget
#'
#' @description Extracts genomic sequences from one or more genome assemblies
#'   that correspond to a set of target loci. Targets can be provided either as
#'   a FASTA file (matched via dc-megablast) or as a BED coordinate file.
#'   Matching genomic regions are extracted with optional flanking sequence and
#'   saved as per-genome FASTA files. Optionally also writes BED and match-table
#'   CSV outputs.
#'
#' @param genome.path path to a single genome FASTA file or to a directory of
#'   genome files to process.  May be \code{NULL} when \code{genome.accessions}
#'   is supplied.
#'
#' @param genome.accessions optional path to a CSV file whose rows each
#'   identify one NCBI genome assembly to download, extract targets from, and
#'   then delete.  Use \strong{GCA accession numbers} (GenBank assemblies) —
#'   they cover every genome, including those without a RefSeq (GCF) record.
#'   Required column: \code{Accession} (or the name given in
#'   \code{accession.column}).  Optional column: \code{Name} (or
#'   \code{name.column}) for a human-readable sample label; the accession is
#'   used when absent.  Downloads are performed with the NCBI \code{datasets}
#'   CLI tool one genome at a time and the FASTA is deleted immediately after
#'   extraction to minimise disk usage.
#'
#' @param accession.column name of the column in \code{genome.accessions} that
#'   holds the GCA/GCF accession strings.  Default \code{"Accession"}.
#'
#' @param name.column optional name of the column in \code{genome.accessions}
#'   that holds a human-readable sample label.  When \code{NULL} (default) the
#'   accession string is used as the sample name.
#'
#' @param genome.taxon optional taxonomic name or NCBI taxon ID (e.g.
#'   \code{"Centrolenidae"} or \code{"8404"}) to search NCBI automatically for
#'   matching genome assemblies.  All assemblies at the requested
#'   \code{assembly.level}(s) are discovered, their accessions resolved, and
#'   each is downloaded, processed, and deleted in sequence — no user-provided
#'   table required.  Can be combined with \code{genome.path} and/or
#'   \code{genome.accessions}; results are pooled.
#'
#' @param assembly.level character vector of NCBI assembly level(s) to include
#'   when \code{genome.taxon} is supplied.  One or more of \code{"contig"},
#'   \code{"scaffold"}, \code{"chromosome"}, \code{"complete"}.  Default
#'   \code{c("scaffold","chromosome","complete")} excludes raw contig-level
#'   assemblies, which are usually too fragmented for reliable target
#'   extraction.
#'
#' @param max.assemblies optional integer; maximum number of assemblies to
#'   download when using \code{genome.taxon}.  \code{NULL} (default) downloads
#'   all that match the query.  Useful for large clades or testing.
#'
#' @param ncbi.api.key optional NCBI API key string.  Without a key NCBI
#'   allows 3 requests/second; with one it allows 10/second.  Obtain a free
#'   key at \url{https://www.ncbi.nlm.nih.gov/account/}.  Only needed when
#'   processing large numbers of accessions where rate-limiting would otherwise
#'   slow the FTP-URL lookup step.
#'
#' @param input.file path to the target sequence file; either a multi-sequence
#'   FASTA (when input.type = "fasta") or a BED coordinate file (when
#'   input.type = "bed").
#'
#' @param input.type character; "fasta" to BLAST targets against each genome,
#'   or "bed" to extract coordinates directly from a BED file.
#'
#' @param output.name name of the output directory where results are saved.
#'
#' @param output.bed logical; if TRUE a BED file of matched coordinates is
#'   written for each genome.
#'
#' @param bed.headers logical; if TRUE the BED file is treated as having a
#'   header row (only relevant when input.type = "bed").
#'
#' @param name.bed.names logical; reserved for naming BED entries from the
#'   target name column.
#'
#' @param output.table logical; if TRUE a CSV table of BLAST match statistics
#'   is written for each genome.
#'
#' @param match.by.chr logical; reserved for filtering matches to chromosome-
#'   level sequences only.
#'
#' @param duplicate.matches character; how to handle cases where multiple
#'   targets match the same genomic contig: "best" keeps the highest-scoring
#'   match, "all" keeps all matches, "none" discards duplicates.
#'
#' @param merge.matches logical; if TRUE, multiple BLAST hits from different
#'   parts of the same target to the same contig are merged into a single
#'   coordinate span before extraction.
#'
#' @param minimum.match.length integer minimum number of aligned bases required
#'   to retain a BLAST match.
#'
#' @param minimum.match.identity numeric minimum percent identity (0-100) for a
#'   BLAST match to be retained.
#'
#' @param minimum.match.coverage numeric minimum proportion of the target length
#'   that must be covered by the BLAST match.
#'
#' @param add.flanks integer; base pairs to add upstream and downstream of each
#'   matched region when extracting sequences.  Default 1000 bp, which is
#'   sufficient to capture the full sequenced fragment for most Illumina
#'   paired-end libraries (typical insert size 300-600 bp) with room to spare.
#'   Increase to 1000-5000 bp for UCE analyses where flanking sequence is the
#'   primary phylogenetic signal.  Flanks can always be trimmed back later
#'   during the alignment trimming step.
#'
#' @param genome.search.string optional character string to filter files in a
#'   genome directory by name (e.g. "_genomic.fna.gz").
#'
#' @param threads number of CPU threads to pass to blastn.
#'
#' @param memory amount of RAM in GB (currently reserved).
#'
#' @param max.retries integer; number of times to retry a failed NCBI FTP URL
#'   lookup or genome download before giving up and moving to the next genome.
#'   Default 3.  Between each attempt the function waits \code{retry.delay}
#'   seconds.
#'
#' @param retry.delay numeric; seconds to wait between retry attempts.
#'   Default 60.  Increase for flaky connections or heavily loaded NCBI FTP
#'   servers.
#'
#' @param overwrite logical; if TRUE existing output directories are deleted and
#'   recreated.
#'
#' @param quiet logical; if TRUE BLAST stdout is suppressed.
#'
#' @param blast.path system path to the directory containing blastn and
#'   makeblastdb; NULL searches the system PATH.
#'
#' @return invisibly; extracted FASTA sequences are written to per-genome
#'   sub-directories inside output.name.
#'
#' @export

extractGenomeTarget = function(genome.path = NULL,
                               genome.accessions = NULL,
                               accession.column = "Accession",
                               name.column = NULL,
                               genome.taxon = NULL,
                               assembly.level = c("scaffold", "chromosome", "complete"),
                               max.assemblies = NULL,
                               ncbi.api.key = NULL,
                               max.retries = 3,
                               retry.delay = 60,
                               input.file = NULL,
                               input.type = c("fasta", "bed"),
                               output.name = NULL,
                               output.bed = TRUE,
                               bed.headers = FALSE,
                               name.bed.names = TRUE,
                               output.table = TRUE,
                               match.by.chr = FALSE,
                               duplicate.matches = c("best", "all", "none"),
                               merge.matches = FALSE,
                               minimum.match.length = 100,
                               minimum.match.identity = 0.75,
                               minimum.match.coverage = 0.75,
                               add.flanks = 1000,
                               genome.search.string = NULL,
                               threads = 1,
                               memory = 1,
                               overwrite = FALSE,
                               quiet = FALSE,
                               blast.path = NULL) {

  #todo: make single function for one sample
  #todo: blast type
  #todo: names

  # setwd("/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/Microhylidae_SeqCap")
  # genome.path = "/Volumes/LaCie/Reference_Genomes/Amphibians/Frogs/Gastrophryne_carolinensis_GCA-027917425/GCA_027917425.1_aGasCar1.pri_genomic.fna.gz"
  # input.file = "COGEDA_81_marker.fasta"
  # #input.file = "/Volumes/Armored/Anolis_UCE/Anolis_genome/Anolis_carolinensis_final-table.bed"
  # output.name = "amphibians"
  # input.type = "fasta"
  # match.by.chr = TRUE
  # output.bed = TRUE
  # bed.headers = TRUE
  # output.table = TRUE
  # duplicate.matches = "best"
  # minimum.match.length = 50
  # minimum.match.identity = 0.50
  # minimum.match.coverage = 0.5
  # add.flanks = 500
  # threads = 8
  # memory = 24
  # overwrite = TRUE
  # quiet = FALSE
  # merge.matches = TRUE
  # blast.path = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  # genome.search.string = "_genomic.fna.gz"
#
  # setwd("/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/Elapid_Probes/2_sea_snake_genomes")
  # genome.path = "/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/Elapid_Probes/Gene_based/genomes/GCA_019472885_1_HCur_v2_genomic.fa"
  # input.file = "Hcur_coordinates.txt"
  # blast.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin/"
  #
  # output.name = "Hcur"
  # input.type = "bed"
  # output.bed = TRUE
  # bed.headers = FALSE
  # output.table = FALSE
  # match.by.chr = FALSE
  # duplicate.matches = "none"
  # merge.matches = FALSE
  # minimum.match.length = 100
  # minimum.match.identity = 0.75
  # add.flanks = 30
  # genome.search.string = "_genomic.fna.gz"
  # threads = 8
  # memory = 16
  # overwrite = FALSE
  # quiet = FALSE

  #Same adds to blast path
  if (is.null(blast.path) == FALSE){
    b.string = unlist(strsplit(blast.path, ""))
    if (b.string[length(b.string)] != "/") {
      blast.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { blast.path = "" }

  # Validate choice parameters
  input.type        = match.arg(input.type)
  duplicate.matches = match.arg(duplicate.matches)

  #Initial checks
  if (is.null(genome.path) && is.null(genome.accessions) && is.null(genome.taxon)) {
    stop("Error: supply at least one of genome.path, genome.accessions, or genome.taxon.")
  }

  # Validate assembly.level values
  valid.levels = c("contig", "scaffold", "chromosome", "complete")
  bad.levels = assembly.level[!assembly.level %in% valid.levels]
  if (length(bad.levels) > 0) {
    stop(paste0("Error: invalid assembly.level value(s): ",
                paste(bad.levels, collapse=", "),
                ". Must be one or more of: ", paste(valid.levels, collapse=", ")))
  }
  if (!is.null(output.name) && !is.null(genome.path) && output.name == genome.path){
    stop("You should not overwrite the original genome file.")
  }
  if (output.name == input.file){ stop("You should not overwrite the original target file.") }
  if (is.null(input.file) == T){ stop("A fasta file of targets to match to assembly contigs is needed.") }

  if (dir.exists(output.name) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.name))
      dir.create(output.name)
    }
  } else { dir.create(output.name) }

  ##############################################################################
  # Build a unified list of genome sources.
  # Each entry: list(type, file, path, sample, accession)
  #   type="file"      — local genome file (existing behaviour)
  #   type="accession" — download from NCBI, extract, delete
  ##############################################################################
  genome.sources = list()

  # --- Local files ---
  if (!is.null(genome.path)) {
    if (dir.exists(genome.path)) {
      file.names   = list.files(genome.path, recursive = TRUE)
      genome.files = file.names[grep("\\.fai$", file.names, invert = TRUE)]
      if (!is.null(genome.search.string)) {
        genome.files = genome.files[grep(genome.search.string, genome.files)]
      }
    } else {
      genome.files = basename(genome.path)
      genome.path  = dirname(genome.path)
    }
    for (f in genome.files) {
      sn = gsub(pattern = "\\.fna\\.gz$|\\.fa\\.gz$|\\.fna$|\\.fa$|\\.gz$",
                replacement = "", x = basename(f))
      genome.sources = append(genome.sources,
                              list(list(type="file", file=f,
                                        genome.path=genome.path, sample=sn)))
    }
  }

  # --- NCBI accession downloads ---
  if (!is.null(genome.accessions)) {
    if (!file.exists(genome.accessions)) {
      stop(paste0("Error: genome.accessions file not found: ", genome.accessions))
    }
    acc.table = read.csv(genome.accessions, stringsAsFactors = FALSE)
    if (!accession.column %in% colnames(acc.table)) {
      stop(paste0("Error: column '", accession.column, "' not found in genome.accessions CSV."))
    }
    accessions = acc.table[[accession.column]]
    sample.nms = if (!is.null(name.column) && name.column %in% colnames(acc.table))
      acc.table[[name.column]]
    else
      accessions
    for (k in seq_along(accessions)) {
      genome.sources = append(genome.sources,
                              list(list(type="accession",
                                        accession = accessions[k],
                                        sample    = sample.nms[k])))
    }
  }

  # --- NCBI taxon search (pure R REST API — no CLI tool required) ---
  if (!is.null(genome.taxon)) {
    print(paste0("Querying NCBI for '", genome.taxon, "' assemblies (levels: ",
                 paste(assembly.level, collapse=", "), ")..."))

    key.str = if (!is.null(ncbi.api.key)) paste0("&api_key=", ncbi.api.key) else ""
    # Map user-facing level names to the NCBI v2 API filter values:
    # "complete" must be sent as "complete_genome" in the query string.
    api.levels   = ifelse(assembly.level == "complete", "complete_genome", assembly.level)
    level.params = paste(paste0("filters.assembly_level=", api.levels), collapse="&")
    base.url     = paste0("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/",
                          utils::URLencode(genome.taxon, reserved=TRUE),
                          "/dataset_report?", level.params, "&page_size=1000", key.str)

    # Page through results (NCBI returns max 1000 per page).
    # Use url() + readLines to fetch and jsonlite::fromJSON() to parse,
    # which is more robust than passing the URL string directly to fromJSON.
    ncbi.get.json = function(query.url) {
      tryCatch({
        con  = url(query.url, open="rb")
        on.exit(close(con), add=TRUE)
        txt  = paste(readLines(con, warn=FALSE), collapse="")
        jsonlite::fromJSON(txt, simplifyVector=TRUE)
      }, error=function(e) NULL)
    }

    all.reports = list()
    next.token  = NULL
    repeat {
      page.url = if (is.null(next.token)) base.url
                 else paste0(base.url, "&page_token=",
                             utils::URLencode(next.token, reserved=TRUE))
      resp = ncbi.get.json(page.url)
      if (is.null(resp) || is.null(resp$reports)) break
      all.reports = c(all.reports, list(resp$reports))
      next.token  = resp$next_page_token
      if (is.null(next.token) || identical(next.token, "")) break
      Sys.sleep(if (!is.null(ncbi.api.key)) 0.1 else 0.35)
    }

    if (length(all.reports) == 0) {
      warning(paste0("No assemblies found for taxon: '", genome.taxon, "'"))
    } else {
      reports = do.call(rbind, all.reports)
      taxon.entries = Filter(Negate(is.null), lapply(seq_len(nrow(reports)), function(k) {
        tryCatch({
          acc = reports$accession[k]
          nm  = reports$organism$organism_name[k]
          if (is.null(nm) || is.na(nm)) nm = acc
          nm.clean  = sub("_$", "", gsub("_+", "_", gsub("[^A-Za-z0-9]", "_", nm)))
          acc.clean = gsub("\\.", "_", acc)
          list(accession=acc, sample=paste0(nm.clean, "_", acc.clean))
        }, error=function(e) NULL)
      }))

      if (length(taxon.entries) == 0) {
        warning(paste0("Could not parse assembly records for taxon: '", genome.taxon, "'"))
      } else {
        if (!is.null(max.assemblies) && length(taxon.entries) > max.assemblies) {
          print(paste0("Found ", length(taxon.entries), " assemblies; limiting to ", max.assemblies))
          taxon.entries = taxon.entries[seq_len(max.assemblies)]
        } else {
          print(paste0("Found ", length(taxon.entries), " assemblies for '", genome.taxon, "'"))
        }
        for (entry in taxon.entries) {
          genome.sources = append(genome.sources,
                                  list(list(type="accession",
                                            accession=entry$accession,
                                            sample=entry$sample)))
        }
      }
    }
  }

  if (length(genome.sources) == 0) { stop("Error: no genome sources found.") }

  ##############################################################################
  # Helper closure: resolve NCBI FTP download URL for a genome accession.
  # Uses Entrez esearch (accession → UID) then esummary (UID → FTP path).
  # Returns the HTTPS URL of the _genomic.fna.gz file, or NULL on failure.
  ##############################################################################
  ncbi.genome.url = function(acc) {
    rate.sleep = if (!is.null(ncbi.api.key)) 0.1 else 0.4
    key.str    = if (!is.null(ncbi.api.key)) paste0("&api_key=", ncbi.api.key) else ""

    # Robust JSON fetch helper (reuses the closure defined above if in scope,
    # but redefined here so the helper is self-contained)
    get.json = function(u) {
      tryCatch({
        con = url(u, open="rb")
        on.exit(close(con), add=TRUE)
        jsonlite::fromJSON(paste(readLines(con, warn=FALSE), collapse=""),
                           simplifyVector=TRUE)
      }, error=function(e) NULL)
    }

    Sys.sleep(rate.sleep)
    # Square brackets in the Entrez field qualifier must be percent-encoded
    # so they are not misinterpreted as URL syntax by url() / curl.
    search.url = paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?",
                        "db=assembly&term=", utils::URLencode(acc, reserved=TRUE),
                        "%5BAssembly%20Accession%5D&retmode=json", key.str)
    sr = get.json(search.url)
    if (is.null(sr) || length(sr$esearchresult$idlist) == 0) return(NULL)
    uid = sr$esearchresult$idlist[1]

    Sys.sleep(rate.sleep)
    sum.url = paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?",
                     "db=assembly&id=", uid, "&retmode=json", key.str)
    smr = get.json(sum.url)
    if (is.null(smr)) return(NULL)

    ftp = tryCatch(smr$result[[uid]]$ftppath_genbank, error=function(e) NULL)
    if (is.null(ftp) || ftp %in% c("", "na", "NA")) return(NULL)

    ftp.https = sub("^ftp://", "https://", ftp)
    aname     = basename(ftp.https)
    paste0(ftp.https, "/", aname, "_genomic.fna.gz")
  }


  # if (match.by.chr == TRUE){
  #   chrom.samples = file.names[grep("chromosomes", file.names)]
  #   chrom.names = gsub("/chromosomes.*", "",  chrom.samples)
  #   chrom.names = unique(gsub(".*/", "", chrom.names))
  #   #remove and add from original files
  #   genome.files = genome.files[grep(paste(chrom.names, collapse = "|"), genome.files, invert = T)]
  #   genome.files = append(genome.files, chrom.samples)
  # }#end if


  if (input.type == "fasta"){

    #headers for the blast db
    headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
                "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

    ### Reads in the additional stuff
    target.markers = Biostrings::readDNAStringSet(input.file)  # loads up fasta file

    #Make blast database for the probe loci
    system(paste0(blast.path, "makeblastdb -in ", input.file,
                  " -parse_seqids -dbtype nucl -out ", output.name, "/target_nucl-blast_db"), ignore.stdout = quiet)
  }#end fasta


  #Loops through each genome source (local files and/or NCBI downloads)
  for (i in seq_along(genome.sources)) {
    src = genome.sources[[i]]

    sample      = src$sample
    species.dir = paste0(output.name, "/", sample)
    downloaded  = FALSE   # tracks whether we need to clean up a downloaded genome

    #Creates species directory if none exists
    if (!dir.exists(species.dir)) { dir.create(species.dir) }

    ##########################################################################
    # NCBI download preamble — only for accession sources
    # Uses pure R: Entrez API to resolve FTP URL, download.file() to fetch.
    # No CLI tools required.
    ##########################################################################
    if (src$type == "accession") {
      acc          = src$accession
      sentinel     = file.path(species.dir, paste0(sample, "_COMPLETED"))
      matches.file = file.path(species.dir, paste0(sample, "_target-matches.fa"))

      # ── Already-done check ─────────────────────────────────────────────────
      # A COMPLETED sentinel file is written at the very end of a successful
      # run.  The output FASTA is also checked as a fallback for runs that
      # finished before the sentinel was introduced.
      if (!overwrite &&
          (file.exists(sentinel) || file.exists(matches.file))) {
        print(paste0("[", i, "/", length(genome.sources), "] ",
                     sample, " already completed — skipping."))
        next
      }

      # ── Step 1: resolve FTP URL via Entrez (with retries) ─────────────────
      print(paste0("[", i, "/", length(genome.sources), "] ",
                   "Resolving download URL for ", acc, "..."))
      genome.url = NULL
      for (attempt in seq_len(max.retries)) {
        genome.url = ncbi.genome.url(acc)
        if (!is.null(genome.url)) break
        if (attempt < max.retries) {
          warning(paste0(acc, ": URL lookup attempt ", attempt, " failed — ",
                         "retrying in ", retry.delay, "s..."))
          Sys.sleep(retry.delay)
        }
      }
      if (is.null(genome.url)) {
        warning(paste0(acc, ": could not resolve genome FTP URL after ",
                       max.retries, " attempt(s). Skipping."))
        next
      }

      # ── Step 2: download gzipped genome FASTA (with retries) ──────────────
      genome.gz = file.path(species.dir, paste0(sample, "_genome.fna.gz"))
      print(paste0("  Downloading ", basename(genome.url), "..."))
      dl.ok = FALSE
      for (attempt in seq_len(max.retries)) {
        if (file.exists(genome.gz)) file.remove(genome.gz)
        dl.ret = tryCatch(
          download.file(genome.url, destfile=genome.gz,
                        method="curl", quiet=quiet, mode="wb"),
          error=function(e) 1L
        )
        if (dl.ret == 0 && file.exists(genome.gz) && file.size(genome.gz) > 0) {
          dl.ok = TRUE; break
        }
        if (attempt < max.retries) {
          warning(paste0(acc, ": download attempt ", attempt, " failed — ",
                         "retrying in ", retry.delay, "s..."))
          Sys.sleep(retry.delay)
        }
      }
      if (!dl.ok) {
        warning(paste0(acc, ": download failed after ", max.retries,
                       " attempt(s). Skipping."))
        if (file.exists(genome.gz)) file.remove(genome.gz)
        next
      }

      # ── Step 3: decompress — delete .gz immediately to free disk space ─────
      species.genome.path = file.path(species.dir, paste0(sample, "_genome.fa"))
      system(paste0("gzip -dc \"", genome.gz, "\" > \"", species.genome.path, "\""))
      file.remove(genome.gz)

      downloaded = TRUE
      zipped.up  = integer(0)
    }

    if (input.type == "fasta"){

      #Checks if this has been done already (file sources only; accession already checked above)
      if (src$type == "file") {
        if (overwrite == FALSE){
          if (file.exists(paste0(species.dir, "/", sample, "_target-matches.fa")) == T){ next }
        }#end

        #######################################################################
        #Part A: Resolve genome path for local files
        #######################################################################
        zipped.up = grep("\\.gz$", src$file)

        if (length(zipped.up) == 1){
          system(paste0("gzip -dc \"", file.path(src$genome.path, src$file), "\" > \"",
                        file.path(species.dir, paste0(sample, "_genome.fa")), "\""))
          species.genome.path = file.path(species.dir, paste0(sample, "_genome.fa"))
        } else {
          species.genome.path = file.path(src$genome.path, src$file)
        }#end else
      }#end file source

      #Matches samples to loci
      system(paste0(blast.path, "blastn -task dc-megablast -db ", output.name, "/target_nucl-blast_db -evalue 0.001",
                    " -query ", species.genome.path,
                    " -out ", species.dir, "/", sample, "_target-blast-match.txt",
                    " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
                    " -num_threads ", threads))

      #Matches samples to proteins
      # system(paste0("tblastn -task tblastn-fast -db ", sample, "_nucl-blast_db -evalue 0.001 -seg no",
      #               " -query ", prot.file, " -out ", sample, "_prot-match.txt",
      #               " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
      #               " -num_threads ", threads))
      #  #Matches need to be greater than 12
      #filt.data = match.data[match.data$matches > 12,]
      #Percent identitiy must match 50% or greater
      #filt.data = filt.data[filt.data$pident >= 50,]

      #Loads in match data
      match.data = data.table::fread(paste0(species.dir, "/", sample, "_target-blast-match.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)

      #Guard: empty BLAST output
      if (nrow(match.data) == 0) {
        print(paste0(sample, " had no matches. Skipping"))
        next
      }
      data.table::setnames(match.data, headers)

      #Matches need to be >= minimum.match.length bp
      filt.data = match.data[match.data$matches >= minimum.match.length,]
      #Remove identical rows
      filt.data = filt.data[duplicated(filt.data) != T,]
      #Percent identity: minimum.match.identity is 0-1 fraction; pident is 0-100
      filt.data = filt.data[filt.data$pident >= minimum.match.identity * 100,]

      #Make sure the hit is greater than 50% of the reference length
      filt.data = filt.data[filt.data$matches >= ( (minimum.match.coverage) * filt.data$tLen),]

      #Skips if no matches
      if (nrow(filt.data) == 0) {
        print(paste0(sample, " had no matches. Skipping"))
        next }

      #Sorting: exon name, contig name, bitscore higher first, evalue
      data.table::setorder(filt.data, qName, tName, -pident, -bitscore, evalue)

      #########################################################################
      #Part B: Merge tName that match to multiple same contig qName
      #########################################################################

      if (merge.matches == TRUE){

        merge.names = filt.data$tName[duplicated(filt.data$tName)]
        merge.match = filt.data[filt.data$tName %in% merge.names,]
        merge.data = merge.match[order(merge.match$tName)]

        #Loops through each potential duplicate
        merge.loci = unique(merge.data$tName)

        if (length(merge.loci) != 0){
          save.merge = c()
          del.merge = c()
          for (j in 1:length(merge.loci)){
            #pulls out data that matches to multiple contigs
            sub.data = merge.data[merge.data$tName %in% merge.loci[j],]
            sub.data = sub.data[order(sub.data$tStart)]

            dup.merge = unique(sub.data$qName)
            del.merge = append(del.merge, merge.loci[j])

            for (k in 1:length(dup.merge)){
              #merge duplicates
              m.data = sub.data[sub.data$qName %in% dup.merge[k],]
              m.data = m.data[order(m.data$tStart)]

              if (nrow(m.data) == 1){ save.merge = rbind(save.merge, m.data) }

              if (nrow(m.data) != 1){
                #Merge into one record
                m.data$tStart[1] = min(c(m.data$tEnd, m.data$tStart))
                m.data$tEnd[1] = max(c(m.data$tEnd, m.data$tStart))
                m.data$qStart[1] = min(c(m.data$qEnd, m.data$qStart))
                m.data$qEnd[1] = max(c(m.data$qEnd, m.data$qStart))
                m.data$bitscore[1] = sum(m.data$bitscore)
                m.data$matches[1] = sum(m.data$matches)
                m.data$misMatches[1] = sum(m.data$misMatches)
                m.data$gapopen[1] = sum(m.data$gapopen)

                save.merge = rbind(save.merge, m.data[1,])
              }#end if

            }#end k loop

          }#end j

          #Filters to new dataset
          filt.data = filt.data[!filt.data$tName %in% del.merge,]
          filt.data = rbind(filt.data, save.merge)
        }#end if
      }#end merge.matches if

      #########################################################################
      #Part C: Multiple targets (qName) matching to one sample contig (tName)
      #########################################################################

      dup.names = filt.data$tName[duplicated(filt.data$tName)]
      dup.match = filt.data[filt.data$tName %in% dup.names,]
      dup.data = dup.match[order(dup.match$tName)]

      #Loops through each potential duplicate
      dup.loci = unique(dup.data$tName)

      if (length(dup.loci) != 0){
        save.dup = c()
        for (j in 1:length(dup.loci)){
          #pulls out data that matches to multiple contigs
          sub.data = dup.data[dup.data$tName %in% dup.loci[j],]
          sub.data = sub.data[order(sub.data$tStart)]

          #Here decide what to do with the duplicate matches based on selection
          if (duplicate.matches == "none"){ next }

          #Saves all the decent matches
          if (duplicate.matches == "all"){
            save.temp = sub.data
          }

          #Here find the best match for this set of duplicates
          if (duplicate.matches == "best"){
            save.temp = sub.data[sub.data$bitscore == max(sub.data$bitscore),][1,]
          }

          save.dup = rbind(save.dup, save.temp)
        } #end j loop

        temp.filt = filt.data[!filt.data$tName %in% dup.names,]
        if (duplicate.matches == "none") {
          filt.data = temp.filt
        } else {
          filt.data = rbind(temp.filt, save.dup)
        }

      }#end if

      final.table = filt.data[order(filt.data$tName)]

      #Extracts the genomic data using the final.table coordinates
      Rsamtools::indexFa(species.genome.path)
      fa = Rsamtools::FaFile(species.genome.path)

      #adds genome columns — initialise to NA so skipped rows are filtered below
      final.table[, gStart:=NA_real_]
      final.table[, gEnd:=NA_real_]
      header.data = colnames(final.table)

      for (j in 1:nrow(final.table)){
        #subsets data
        sub.match = final.table[j,]

        #Gets starts and stops
        gen.start = min(c(sub.match$qStart, sub.match$qEnd)) - add.flanks
        gen.end  = max(c(sub.match$qEnd, sub.match$qStart)) + add.flanks
        if (gen.start <= 0){ gen.start = as.numeric(1) }
        if (gen.end > sub.match$qLen) { gen.end = sub.match$qLen }
        # Skip if the extracted span is implausibly large (> 100x target length)
        # Leave gStart/gEnd as NA so this row is removed by the filter below
        if (gen.end - gen.start >= max(sub.match$tLen) * 100){ next }

        ## Save in table here
        data.table::set(final.table, i =  as.integer(j),
                        j = match("gStart", header.data), value = gen.start )

        data.table::set(final.table, i =  as.integer(j),
                        j = match("gEnd", header.data), value = gen.end )

      }#end j loop

      #removes rows that were skipped (gStart is still NA)
      final.table = final.table[!is.na(final.table$gStart),]

      #gets ranges of stuff and obtains sequences
      gr.ranges = data.frame(seqnames = final.table$qName,
                             start = final.table$gStart,
                             end = final.table$gEnd)
      fin.ranges = GenomicRanges::makeGRangesFromDataFrame(gr.ranges)

      #gets ranges of stuff and obtains sequences
      target.seqs = BSgenome::getSeq(fa, fin.ranges)
      names(target.seqs) = paste0(final.table$tName, "_|_", sample)

      #Writes the table
      if (output.table == TRUE){
        write.csv(final.table, file = paste0(species.dir, "/", sample, "_final-table.csv"),
                  row.names = FALSE)
      }

      #Writes the bed
      if (output.bed == TRUE){
        bam.table = data.frame(chr = final.table$qName, start = final.table$gStart, stop = final.table$gEnd)
        write.table(bam.table, file = paste0(species.dir, "/", sample, "_final-table.bed"),
                    row.names = FALSE,
                    quote = F,
                    sep = "\t")
      }#end if


    }#end input type fasta

    if (input.type == "bed"){
      # BED mode only makes sense for local files; skip accession sources
      if (src$type == "accession") {
        warning(paste0(acc, ": input.type='bed' not supported for NCBI downloads. Skipping."))
        # Clean up any partially downloaded genome file
        genome.fa = file.path(species.dir, paste0(sample, "_genome.fa"))
        if (file.exists(genome.fa)) file.remove(genome.fa)
        next
      }

      zipped.up = grep("\\.gz$", src$file)

      if (length(zipped.up) == 1){
        system(paste0("gzip -dc \"", file.path(src$genome.path, src$file), "\" > \"",
                      file.path(species.dir, paste0(sample, "_genome.fa")), "\""))
        species.genome.path = file.path(species.dir, paste0(sample, "_genome.fa"))
      } else {
        species.genome.path = file.path(src$genome.path, src$file)
      }#end else

      #Extracts the genomic data using the final.table coordinates
      Rsamtools::indexFa(species.genome.path)
      fa = Rsamtools::FaFile(species.genome.path)

      #gets ranges of stuff and obtains sequences
      final.table = read.table(input.file, header = bed.headers)

      if (bed.headers == TRUE){
        gr.ranges = data.frame(seqnames = final.table$chrom,
                               start = as.integer(final.table$start),
                               end = as.integer(final.table$stop) )
        fin.ranges = GenomicRanges::makeGRangesFromDataFrame(gr.ranges)

        #gets ranges of stuff and obtains sequences
        target.seqs = BSgenome::getSeq(fa, fin.ranges)
        names(target.seqs) = paste0(final.table$name, "_|_", sample)

      } else {
        gr.ranges = data.frame(seqnames = final.table[,1],
                               start = as.integer(final.table[,2]),
                               end = as.integer(final.table[,3]) )
        fin.ranges = GenomicRanges::makeGRangesFromDataFrame(gr.ranges)

        #gets ranges of stuff and obtains sequences
        target.seqs = BSgenome::getSeq(fa, fin.ranges)
        names(target.seqs) = paste0(gr.ranges$seqnames, "_", gr.ranges$start, "-",gr.ranges$end, "_|_", sample)
      }#end else

    }#end bed

    #Writes the final loci
    final.loci = as.list(as.character(target.seqs))
    PhyloProcessR::writeFasta(sequences = final.loci, names = names(final.loci),
               paste0(species.dir, "/", sample, "_target-matches.fa"),
               nbchar = 1000000, as.string = T)

    print(paste0("[", i, "/", length(genome.sources), "] ",
                 sample, " complete — ", length(final.loci), " targets extracted."))

    # Write a sentinel file so re-runs can skip this genome without re-checking
    # every output file.  Includes a brief summary for auditing.
    if (src$type == "accession") {
      writeLines(
        c(paste0("accession: ", src$accession),
          paste0("sample:    ", sample),
          paste0("targets:   ", length(final.loci)),
          paste0("completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))),
        con = file.path(species.dir, paste0(sample, "_COMPLETED"))
      )
    }

    # Clean up temporary genome files
    if (downloaded) {
      # Accession source: the .gz was already deleted after decompression.
      # Remove the decompressed genome FASTA and its index — output files are kept.
      genome.fa = file.path(species.dir, paste0(sample, "_genome.fa"))
      fai.path  = paste0(genome.fa, ".fai")
      if (file.exists(genome.fa)) file.remove(genome.fa)
      if (file.exists(fai.path))  file.remove(fai.path)
    } else if (length(zipped.up) == 1) {
      # Local gzipped file: delete the decompressed copy
      system(paste0("rm \"", file.path(species.dir, paste0(sample, "_genome.fa")), "\""))
      fai.path = file.path(species.dir, paste0(sample, "_genome.fa.fai"))
      if (file.exists(fai.path)) { system(paste0("rm \"", fai.path, "\"")) }
    }

  }# end i loop

} #End function


#END SCRIPT
