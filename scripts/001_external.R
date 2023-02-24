#### Functions which call external software from R
# Brendan Furneaux 2022

find_vsearch <- function() {
  vsearch <- Sys.getenv("VSEARCH")
  if (nchar(vsearch) == 0 || !file.exists(vsearch)) {
    vsearch <- Sys.which("vsearch")
  }
  if (nchar(vsearch) == 0 || !file.exists(vsearch)) {
    stop("cannot find vsearch")
  }
  vsearch
}

find_cutadapt <- function() {
  cutadapt <- Sys.getenv("CUTADAPT")
  if (nchar(cutadapt) == 0 || !file.exists(cutadapt)) {
    cutadapt <- Sys.which("cutadapt")
  }
  if (nchar(cutadapt) == 0 || !file.exists(cutadapt)) {
    stop("cannot find cutadapt")
  }
  cutadapt
}

vsearch_usearch_global <- function(query, ref, threshold, global = TRUE, ncpu = local_cpus()) {
  tquery <- tempfile("query", fileext = ".fasta")
  on.exit(unlink(c(tquery), force = TRUE))
  write_sequence(query, tquery)
  if (is.character(ref) && length(ref) == 1 && file.exists(ref)) {
    tref <- ref
  } else {
  tref <- tempfile("ref", fileext = ".fasta")
  on.exit(unlink(c(tref), force = TRUE), add = TRUE)
  write_sequence(ref, tref)
  }
  assertthat::assert_that(assertthat::is.flag(global))
  gap <- if (global) "1" else "1I/0E"
  uc = system(
    paste(
      find_vsearch(),
      "--usearch_global", tquery,
      "--db", tref,
      "--id", threshold,
      "--uc", "-",
      "--maxaccepts", "100",
      "--top_hits_only",
      "--threads", ncpu,
      "--gapopen", gap,
      "--gapext", gap,
      "--match", "1",
      "--mismatch", "-1",
      "| awk '$1==\"H\" {print $9,$10}'"
    ),
    intern = TRUE
  )
  stopifnot(attr(uc, "status") == 0)
  if (length(uc) > 0) {
    readr::read_delim(
      I(uc),
      col_names = c("seq_id", "cluster"),
      delim = " ",
      col_types = "cc"
    )
  } else {
    tibble::tibble(seq_id = character(), cluster = character())
  }
}

chimera_callback <- function(x, pos) {
  dplyr::mutate(
    x,
    qcov_length = qend - qstart,
    tcov_length = tend - tstart,
    qcov_pct = length/qcov_length,
    tcov_pct = length/tcov_length
  ) %>%
    dplyr::filter(
      (qcov_pct <= 0.8 & tcov_pct <= 0.8) |
        (qcov_pct <= 0.85 & tcov_pct <= 0.85) |
        (length < 200 & qend >= 200)
    )
}

vsearch_usearch_global_blast6out <- function(
  query,
  ref,
  threshold,
  strand = c("plus", "both"),
  ncpu = local_cpus(),
  vsearch = find_vsearch()
) {
  strand <- match.arg(strand)
  tquery <- tempfile("query", fileext = ".fasta")
  on.exit(unlink(c(tquery), force = TRUE))
  write_sequence(query, tquery)
  if (is.character(ref) && length(ref) == 1 && file.exists(ref)) {
    tref <- ref
  } else {
    tref <- tempfile("ref", fileext = ".fasta")
    on.exit(unlink(c(tref), force = TRUE), add = TRUE)
    write_sequence(ref, tref)
  }
  blast6out = pipe(
    paste(
      find_vsearch(),
      "--usearch_global", tquery,
      "--db", tref,
      "--id", threshold,
      "--blast6out", "-",
      "--top_hits_only",
      "--threads", ncpu
    )
  )
  readr::read_tsv_chunked(
    blast6out,
    callback = readr::DataFrameCallback$new(chimera_callback),
    col_names = c("query", "target", "id", "length", "nmismatch", "ngap",
                  "qstart", "qend", "tstart", "tend", "e", "score"),
    col_types = "ccniiiiiiinn"
  )
}

vsearch_uchime_ref <- function(query, ref, ncpu = local_cpus()) {
  tquery <- tempfile("query", fileext = ".fasta")
  on.exit(unlink(c(tquery), force = TRUE))
  write_sequence(query, tquery)
  if (is.character(ref) && length(ref) == 1 && file.exists(ref)) {
    tref <- ref
  } else {
    tref <- tempfile("ref", fileext = ".fasta")
    on.exit(unlink(c(tref), force = TRUE), add = TRUE)
    write_sequence(ref, tref)
  }
  tchimeras <- tempfile("chimeras", fileext = ".fasta")
  on.exit(unlink(tchimeras), TRUE)
  vs <- system2(
    find_vsearch(),
    c(
      "--uchime_ref", tquery,
      "--db", tref,
      "--chimeras", tchimeras,
      "--threads", ncpu
    )
  )
  stopifnot(vs == 0L)
  Biostrings::readDNAStringSet(tchimeras) %>%
    as.character() %>%
    tibble::enframe(name = "seq_id", value = "seq")
}

vsearch_usearch_global_closed_ref <- function(query, ref, threshold, ...) {
  out <- tibble::tibble(seq_id = character(0), cluster = character(0))
  while(sequence_size(query) > 0 && sequence_size(ref) > 0) {
    result <- vsearch_usearch_global(query, ref, threshold, ...)
    if (nrow(out) > 0) {
      result <- dplyr::left_join(
        result,
        out,
        by = c("cluster" = "seq_id"),
        suffix = c(".orig", "")
      ) %>%
        dplyr::select(seq_id, cluster)
    }
    out <- dplyr::bind_rows(out, result)
    ref <- select_sequence(query, result$seq_id)
    query <- select_sequence(query, result$seq_id, negate = TRUE)
  }
  out
}

# build a usearch database (UDB) file using USEARCH

build_udb <- function(infile, outfile, type = c("usearch", "sintax", "ublast"),
                      usearch = Sys.which("usearch")) {
  type <- match.arg(type)
  command <- paste0("-makeudb_", type)
  args <- c(
    command, infile,
    "-output", outfile
  )
  result <- system2(usearch, args)
  stopifnot(result == 0)
  outfile
}

build_filtered_udb <- function(
  infile,
  outfile,
  type = c("usearch", "sintax", "ublast"),
  blacklist,
  usearch = Sys.which("usearch")
) {
  # make sure we have valid arguments
  type <- match.arg(type)
  command <- paste0("-makeudb_", type)
  stopifnot(system2(usearch, "--version")==0)
  
  # make a temp file and a temp fifo
  blf <- tempfile(fileext = ".txt")
  tf <- tempfile(fileext = ".fasta")
  on.exit(unlink(c(tf, blf), force = TRUE))
  writeLines(blacklist, blf)
  stopifnot(system2("mkfifo", tf) == 0)
  
  # first usearch call removes the blacklisted sequences
  system2(
    usearch,
    args = c(
      "--fastx_getseqs", infile,
      "--labels", blf,
      "--label_substr_match",
      "--notmatched", tf
    ),
    wait = FALSE
  )
  # second usearch call creates the udb file
  result = system2(
    usearch,
    args = c(
      command, tf,
      "--output", outfile
    )
  )
  stopifnot(result == 0)
  outfile
}

blastclust_usearch <- function(
  seq,
  threshold,
  seq_id = names(seq),
  which = TRUE,
  ncpu = local_cpus(),
  usearch = Sys.which("usearch")
) {
  is_file <- length(seq) == 1 && file.exists(seq)
  if (isFALSE(is_file)) {
    stopifnot(length(seq) == length(seq_id))
    if (length(seq) == 0) {
      return(character(0))
    } else if (length(seq) == 1) {
      return(paste0(seq_id, " "))
    }
  }
  
  hits <- tempfile("hits")
  on.exit(unlink(hits, force = TRUE), TRUE)
  usearch_hitlist(
    seq,
    threshold = threshold/100,
    seq_id = seq_id,
    which = which,
    ncpu = ncpu,
    hits = hits,
    usearch = usearch
  )
  blastclust_reclust(
    hits,
    threshold = threshold,
    ncpu = ncpu
  )
}

vsearch_cluster_smallmem <- function(seq, threshold = 1, ncpu = local_cpus()) {
  tout <- tempfile("data", fileext = ".fasta")
  tin <- tempfile("data", fileext = ".uc")
  on.exit(unlink(tout, force = TRUE))
  write_sequence(seq, tout)
  uc = system(
    paste(
      find_vsearch(),
      "--cluster_smallmem", tout,
      "--usersort",
      "--id", threshold,
      "--uc -",
      "--threads", ncpu,
      "| awk '$1==\"H\" {print $9,$10}'"
    ),
    intern = TRUE
  )
  stopifnot(attr(uc, "status") == 0)
  if (length(uc) > 0) {
    readr::read_delim(
      I(uc),
      col_names = c("query", "hit"),
      delim = " ",
      col_types = "cc"
    )
  } else {
    tibble::tibble(query = character(), hit = character())
  }
}

collapseNoMismatch_vsearch <- function(seqtab, ncpu = local_cpus()) {
  seqs <- colnames(seqtab)
  names(seqs) <- seq_along(seqs)
  matches <- vsearch_cluster_smallmem(seqs, ncpu = ncpu)
  if (nrow(matches) > 0) {
    matches$query <- as.integer(matches$query)
    matches$hit <- as.integer(matches$hit)
    for (i in unique(matches$hit)) {
      seqtab[,i] <- seqtab[,i] +
        as.integer(rowSums(seqtab[,matches$query[matches$hit == i], drop = FALSE]))
    }
    seqtab <- seqtab[,-matches$query]
  }
  return(seqtab)
}

cutadapt_paired_filter_trim <- function(
  file_R1, file_R2,
  primer_R1, primer_R2,
  trim_R1, trim_R2,
  max_err = NULL, min_overlap = NULL,
  action = NULL,
  discard_untrimmed = FALSE,
  max_n = NULL,
  max_ee = NULL,
  min_length = NULL, max_length = NULL,
  truncQ_R1 = NULL, truncQ_R2 = NULL,
  cut_R1 = NULL, cut_R2 = NULL,
  ncpu = local_cpus(),
  cutadapt = find_cutadapt(),
  ...
) {
  args <- c(
    "-g", primer_R1,
    "-G", primer_R2,
    "-o", trim_R1,
    "-p", trim_R2
  )
  if (!is.null(max_err)) {
    assertthat::assert_that(assertthat::is.number(max_err))
    args <- c(args, "-e", max_err)
  }
  if (!is.null(min_overlap)) {
    assertthat::assert_that(assertthat::is.number(min_overlap))
    args <- c(args, "-O", min_overlap)
  }
  if (!is.null(action)) {
    args <- c(args, paste0("--action=", action))
  }
  if (isTRUE(discard_untrimmed)) {
    args <- c(args, "--discard-untrimmed")
  }
  if (!is.null(max_n)) {
    assertthat::assert_that(assertthat::is.number(max_n))
    args <- c(args, "--max-n", max_n)
  }
  if (!is.null(max_ee)) {
    assertthat::assert_that(assertthat::is.number(max_ee))
    args <- c(args, "--max-ee", max_ee)
  }
  if (!is.null(min_length)) {
    assertthat::assert_that(assertthat::is.count(min_length))
    args <- c(args, "-m", min_length)
  }
  if (!is.null(max_length)) {
    assertthat::assert_that(assertthat::is.count(max_length))
    args <- c(args, "-M", max_length)
  }
  if (!is.null(truncQ_R1)) {
    assertthat::assert_that(
      rlang::is_integerish(truncQ_R1),
      length(truncQ_R1) >= 1,
      length(truncQ_R1) <= 2
    )
    args <- c(args, "-q", paste(truncQ_R1, collapse = ","))
  }
  if (!is.null(truncQ_R2)) {
    assertthat::assert_that(
      rlang::is_integerish(truncQ_R2),
      length(truncQ_R2) >= 1,
      length(truncQ_R2) <= 2
    )
    args <- c(args, "-Q", paste(truncQ_R2, collapse = ","))
  }
  if (!is.null(cut_R1)) {
    assertthat::assert_that(
      rlang::is_integerish(cut_R1),
      length(cut_R1) >= 1,
      length(cut_R1) <= 2
    )
    args <- c(args, "-u", paste(cut_R1, collapse = ","))
  }
  if (!is.null(cut_R2)) {
    assertthat::assert_that(
      rlang::is_integerish(cut_R2),
      length(cut_R2) >= 1,
      length(cut_R2) <= 2
    )
    args <- c(args, "-U", paste(cut_R2, collapse = ","))
  }
  if (!is.null(ncpu)) {
    assertthat::assert_that(assertthat::is.count(ncpu))
    args <- c(args, "-j", ncpu)
  }
  args <- c(args, file_R1, file_R2)
  out <- system2(
    cutadapt,
    args = shQuote(args),
  )
  stopifnot(out == 0)
  c(trim_R1, trim_R2)
}

cutadapt_filter_trim <- function(
  file,
  primer,
  trim,
  max_err = NULL, min_overlap = NULL,
  action = NULL,
  discard_untrimmed = FALSE,
  max_n = NULL,
  max_ee = NULL,
  min_length = NULL, max_length = NULL,
  truncQ = NULL,
  cut = NULL,
  ncpu = local_cpus(),
  cutadapt = find_cutadapt(),
  ...
) {
  args <- c(
    "-g", primer,
    "-o", trim
  )
  if (!is.null(max_err)) {
    assertthat::assert_that(assertthat::is.number(max_err))
    args <- c(args, "-e", max_err)
  }
  if (!is.null(min_overlap)) {
    assertthat::assert_that(assertthat::is.number(min_overlap))
    args <- c(args, "-O", min_overlap)
  }
  if (!is.null(action)) {
    args <- c(args, paste0("--action=", action))
  }
  if (isTRUE(discard_untrimmed)) {
    args <- c(args, "--discard-untrimmed")
  }
  if (!is.null(max_n)) {
    assertthat::assert_that(assertthat::is.number(max_n))
    args <- c(args, "--max-n", max_n)
  }
  if (!is.null(max_ee)) {
    assertthat::assert_that(assertthat::is.number(max_ee))
    args <- c(args, "--max-ee", max_ee)
  }
  if (!is.null(min_length)) {
    assertthat::assert_that(assertthat::is.count(min_length))
    args <- c(args, "-m", min_length)
  }
  if (!is.null(max_length)) {
    assertthat::assert_that(assertthat::is.count(max_length))
    args <- c(args, "-M", max_length)
  }
  if (!is.null(truncQ)) {
    assertthat::assert_that(
      rlang::is_integerish(truncQ),
      length(truncQ) >= 1,
      length(truncQ) <= 2
    )
    args <- c(args, "-q", paste(truncQ, collapse = ","))
  }
  if (!is.null(cut)) {
    assertthat::assert_that(
      rlang::is_integerish(cut),
      length(cut) >= 1,
      length(cut) <= 2
    )
    args <- c(args, "-u", paste(cut, collapse = ","))
  }
  if (!is.null(ncpu)) {
    assertthat::assert_that(assertthat::is.count(ncpu))
    args <- c(args, "-j", ncpu)
  }
  args <- c(args, file)
  out <- system2(
    cutadapt,
    args = shQuote(args),
  )
  stopifnot(out == 0)
  trim
}

trim_primer <- function(seqs, primer, ...) {
  tempseqs <- tempfile(fileext = ".fasta")
  write_sequence(seqs, tempseqs)
  temptrimmed <- tempfile(fileext = ".fasta")
  on.exit(unlink(c(tempseqs, temptrimmed), force = TRUE))
  cutadapt_filter_trim(
    file = tempseqs,
    primer = primer,
    trim = temptrimmed,
    ...
  )
  Biostrings::readDNAStringSet(temptrimmed) %>%
    as.character() %>%
    tibble::enframe(name = "seq_id", value = "seq")
}

trim_seqtable <- function(seqtable, primer, ...) {
  tempseqs <- tempfile(fileext = ".fasta")
  colnames(seqtable) %>%
    Biostrings::DNAStringSet() %>%
    Biostrings::writeXStringSet(tempseqs)
  temptrimmed <- tempfile(fileext = ".fasta")
  on.exit(unlink(c(tempseqs, temptrimmed), force = TRUE))
  cutadapt_filter_trim(
    file = tempseqs,
    primer = primer,
    trim = temptrimmed,
    ...
  )
  Biostrings::readDNAStringSet(temptrimmed) %>%
    as.character() %>%
    magrittr::set_colnames(seqtable, .)
}

run_protax <- function(seqs, outdir, modeldir, ncpu = local_cpus()) {
  if (dir.exists(outdir)) unlink(outdir, recursive = TRUE)
  dir.create(outdir)
  write_sequence(seqs, file.path(outdir, "all.fa"))
  status <- system2(
    "scripts/runprotax",
    c(outdir, modeldir, ncpu)
  )
  stopifnot(status == 0)
  list.files(outdir, full.names = TRUE)
}

#' Cluster ASVs using blastclust
#'
#' @param seq (`character` vector, file name, or `Biostrings::DNAStringSet`)
#' sequences to cluster
#' @param threshold (`numeric` scalar) percentage similarity threshold for
#' clustering. A number between 3 and 100
#' @param seq_id (`character`) names for the sequences (Default: `names(seq)`)
#' @param which (`logical`, `integer`, or `character`) subset indices
#' indicating which sequences in `seq` should be clustered. Has no effect if
#' `seq` is a file name. (Default: `TRUE`, i.e. all the sequences)
#' @param ncpu (`integer` scalar) number of CPUs to use
#' @param outfile (`character` string giving a valid filename) filename to write
#' clustering output, if desired.
#' @param hits (`character` string giving a valid filename) filename to
#' write the hit table, if desired. For use in `blastclust_reclust()`.
#'
#' @return a `character` vector, where each line gives the sequence names
#' included in one cluster
#' @export
#'
#' @examples
blastclust <- function(seq, threshold, seq_id = names(seq), which = TRUE,
                       ncpu = local_cpus(), outfile = NULL, hits = NULL) {
  UseMethod("blastclust", seq)
}

blastclust.character <- function(seq, threshold, seq_id = names(seq),
                                 which = TRUE, ncpu = local_cpus(),
                                 outfile = NULL, hits = NULL) {
  if (length(seq) == 1 && file.exists(seq)) {
    if (!missing(seq_id))
      warning("'seq_id' has no effect when 'seq' is a file.")
    if (!missing(which))
      warning("'which' has no effect when 'seq' is a file.")
    blastclust_filename(seq, threshold, ncpu)
  }
  seq <- Biostrings::DNAStringSet(seq)
  blastclust.DNAStringSet(seq, threshold, seq_id, which, ncpu, outfile, hits)
}

blastclust.DNAStringSet <- function(seq, threshold, seq_id = names(seq),
                                    which = TRUE,
                                    ncpu = local_cpus(),
                                    outfile = NULL, hits = NULL) {
  # rename the sequences if necessary
  if (!isTRUE(all.equal(names(seq), seq_id))) names(seq) <- seq_id
  seq <- seq[which]
  # shortcut if only one sequence
  if (length(seq) == 1) return(names(seq))
  tf <- tempfile(pattern = "clust", fileext = ".fasta")
  Biostrings::writeXStringSet(seq, tf)
  on.exit(unlink(tf))
  blastclust_filename(tf, threshold, ncpu, outfile, hits)
}

blastclust_filename <- function(seq_file, threshold, ncpu, outfile = NULL, hits = NULL) {
  if (is.null(outfile)) {
    outfile <- tempfile(pattern = "out")
    on.exit(unlink(outfile), TRUE)
  }
  stopifnot(
    system2(
      "blastclust",
      c(
        "-i", seq_file, # input file
        "-S", threshold, # similarity threshold
        "-L", "0.95", # overlap threshold
        "-e", "F", # parse sequence names (FALSE)
        "-b", "F", # require coverage on both neighbors (false)
        "-W", "16", # word (kmer) size 16
        "-a", ncpu, # number of threads
        "-o", outfile, # output file
        if (!is.null(hits)) c("-s", hits) else NULL,
        "-p", "F" # sequences are proteins (FALSE)
      )
    ) == 0
  )
  readLines(outfile)
}


#' Recluster sequences already clustered with `blastclust()`
#'
#' @param hits (`character` filename) hit table generated by a call to
#' `blastclust()`
#' @param threshold (`numeric` scalar) clustering threshold
#' @param preclusters (`character` filename or vector) clustering results from
#' previous call to `blastclust()` or `blastclust_reclust()` at a lower
#' threshold (optional).
#' @param ncpu (`integer` scalar) number of threads to use
#' @param outfile (`character` filename) file to save output (e.g., to use for
#' `preclusters` in future calls to `blastclust_reclust()`)
#'
#' @return a `character` vector, where each line gives the sequence names
#' included in one cluster
#' @export
#'
#' @examples
blastclust_reclust <- function(hits, threshold, preclusters = NULL,
                               ncpu = local_cpus(), outfile = NULL) {
  if (is.list(preclusters) || length(preclusters) > 1L && (all(file.exists(preclusters)) || all(endsWith(preclusters, " ")))) {
    out <- lapply(preclusters, blastclust_reclust, hits = hits,
                  threshold = threshold, ncpu = ncpu, outfile = NULL)
    out <- unlist(out)
    if (!is.null(outfile)) {
      writeLines(out, outfile)
    }
    return(out)
  }
  if (!is.null(preclusters)) {
    stopifnot(is.character(preclusters))
    if (length(preclusters) > 1 || !file.exists(preclusters)) {
      preclusters <- trimws(preclusters)
      if (!grepl(" ", preclusters, fixed = TRUE)) {
        # shortcut if there is only one sequence
        if (!is.null(outfile)) {
          writelines(preclusters, outfile, sep = " ")
        }
        return(paste0(preclusters, " "))
      }
      tf <- tempfile("preclusters")
      writeLines(preclusters, tf, sep = " ")
      on.exit(unlink(tf))
      preclusters = tf
    }
  }
  if (is.null(outfile)) {
    outfile <- tempfile(pattern = "out")
    on.exit(unlink(outfile, force = TRUE), TRUE)
  }
  stopifnot(
    system2(
      "blastclust",
      c(
        "-r", hits, # input file
        if (!is.null(preclusters)) c("-l", preclusters),
        "-S", threshold, # similarity threshold
        "-L", "0.5", # overlap threshold
        "-e", "F", # parse sequence names (FALSE)
        "-b", "F", # require coverage on both neighbors (false)
        "-W", "16", # word (kmer) size 16
        "-a", "8", # number of threads
        "-p", "F", # sequences are proteins (FALSE)
        "-o", outfile # output file
      )
    ) == 0
  )
  readLines(outfile)
}

usearch_hitlist <- function(seq, threshold, seq_id = names(seq), which = TRUE,
                            ncpu = local_cpus(), hits = NULL,
                            usearch = Sys.which("usearch"),
                            timeout = c(600, 60),
                            buffer = 1000) {
  UseMethod("usearch_hitlist", seq)
}

usearch_hitlist.character <- function(seq, threshold, seq_id = names(seq),
                                      which = TRUE, ncpu = local_cpus(),
                                      hits = NULL, usearch = Sys.which("usearch"),
                                      timeout = c(600, 60),
                                      buffer = 1000) {
  if (length(seq) == 1 && file.exists(seq)) {
    if (!missing(seq_id))
      warning("'seq_id' has no effect when 'seq' is a file.")
    if (!missing(which))
      warning("'which' has no effect when 'seq' is a file.")
    index <- Biostrings::fasta.seqlengths(seq)
    do_usearch_hitlist(
      seq,
      seq_len = index,
      seq_id = names(index),
      threshold = threshold,
      ncpu = ncpu,
      hits = hits,
      usearch = usearch,
      timeout = timeout,
      buffer = buffer
    )
  } else {
    mycall <- match.call()
    mycall[[1]] <- usearch_hitlist.DNAStringSet
    mycall$seq <- quote(Biostrings::DNAStringSet(seq))
    eval(mycall)
  }
}

usearch_hitlist.DNAStringSet <- function(seq, threshold, seq_id = names(seq),
                                         which = TRUE,
                                         ncpu = local_cpus(),
                                         hits = NULL,
                                         usearch = Sys.which("usearch"),
                                         timeout = c(600, 60),
                                         buffer = 1000) {
  # rename the sequences if necessary
  if (!isTRUE(all.equal(names(seq), seq_id))) names(seq) <- seq_id
  seq <- seq[which]
  # shortcut if only one sequence
  if (length(seq) == 1) return(names(seq))
  tf <- tempfile(pattern = "clust", fileext = ".fasta")
  Biostrings::writeXStringSet(seq, tf)
  on.exit(unlink(tf))
  do_usearch_hitlist(tf, seq_len = Biostrings::nchar(seq), seq_id = names(seq),
                     threshold = threshold, ncpu = ncpu, hits = hits,
                     usearch = usearch, timeout = timeout, buffer = buffer)
}

do_usearch_hitlist <- function(seq_file, seq_len, seq_id, threshold, ncpu, hits,
                               usearch = Sys.which("usearch"),
                               timeout = c(600, 60),
                               buffer = 1000) {
  assertthat::assert_that(
    is.numeric(timeout),
    length(timeout) >= 1,
    assertthat::is.count(buffer)
  )
  if (length(timeout) == 1) timeout <- rep(timeout, 2) 
  if (!methods::is(hits, "connection")) {
    hits <- file(hits, open = "wb")
  }
  if (!isOpen(hits)) open(hits, "wb")
  
  on.exit(close(hits), TRUE)
  # list type 0 (names)
  # list size (characters in names list)
  writeBin(c(0L, sum(nchar(seq_id)) + length(seq_id)), hits, size = 4L)
  # list of names
  writeChar(paste0(seq_id, " ", collapse = ""), hits, eos = NULL)
  # sequence lengths
  writeBin(seq_len, hits, size = 4L)
  # map from name to (C-style) index
  seqidx <- seq_along(seq_id) - 1L
  names(seqidx) <- seq_id
  if (is.null(names(seq_len))) names(seq_len) <- seq_id
  # usearch wants to write to a file, we want to avoid the file system
  # so give it a fifo
  fifoname <- tempfile("fifo")
  stopifnot(system2("mkfifo", fifoname) == 0)
  on.exit(unlink(fifoname), TRUE)
  f <- fifo(fifoname)
  system2(
    usearch,
    c(
      "-calc_distmx", seq_file, # input file
      "-tabbedout", fifoname, # output fifo
      "-maxdist", 1-threshold, # similarity threshold
      "-termdist", min(1, 1.5*(1-threshold)), # threshold for udist
      "-lopen", "1", # gap opening
      "-lext", "1", # gap extend
      # "-pattern", "111010010111", # pattern gives better result than kmers maybe?
      "-threads", ncpu
    ),
    wait = FALSE
  )
  open(f, mode = "r")
  on.exit(close(f), TRUE)
  tryCatch(
    R.utils::withTimeout(
      d <- readLines(f, n = buffer),
      timeout = timeout[1] # 10 min for kmer indexing + 1000 lines
    ),
    TimeoutException = function(e) stop("Timed out waiting for USEARCH results")
  )
  while (length(d) > 0) {
    d <- strsplit(d, "\t", fixed = TRUE)
    d <- unlist(d)
    d <- matrix(d, ncol = 3, byrow = TRUE)
    d <- as.data.frame(d)
    seq1 <- matrix(writeBin(seqidx[d$V1], raw(), size = 4), nrow = 4)
    seq2 <- matrix(writeBin(seqidx[d$V2], raw(), size = 4), nrow = 4)
    hsp1 <- seq_len[d$V1]
    hsp2 <- seq_len[d$V2]
    ident <- 1 - as.numeric(d$V3)
    score <- matrix(writeBin(ident * pmin(hsp1, hsp2), raw(), size = 8), nrow = 8)
    hsp1 <- matrix(writeBin(hsp1, raw(), size = 4), nrow = 4)
    hsp2 <- matrix(writeBin(hsp2, raw(), size = 4), nrow = 4)
    ident <- matrix(writeBin(100*ident, raw(), size = 8), nrow = 8)
    d <- c(rbind(seq1, seq2, hsp1, hsp2, score, ident))
    writeBin(d, hits)
    tryCatch(
      R.utils::withTimeout(
        d <- readLines(f, n = buffer),
        timeout = timeout[2] # 1 min for 1000 lines
      ),
      TimeoutException = function(e) stop("Timed out waiting for USEARCH results")
    )
  }
}

#' Do single-linkage clustering at a series of increasing similarity thresholds.
#'
#' @param seqs (`character` vector, filename, or `Biostrings::DNAStringSet`)
#' sequences to cluster
#' @param threshold (`numeric`) increasing percentage similarity thresholds for
#' clustering. Number between 3 and 100.
#' @param seqnames (`character` vector) names for the sequences.  If they are
#' already named, this will replace the names.  Has no effect if `seqs` is a
#' filename.
#' @param threshold_name (`character` vector) names for the thresholds. If they
#' are already named, this will replace the names.
#' @param which
#' @param ncpu
#' @param ...
#'
#' @return `list` of `character` vectors giving clustering results (as in )
#' @export
#'
#' @examples
blastclust_repeat <- function(seq, threshold, seq_id = names(seq),
                              which = TRUE,
                              threshold_name = names(threshold),
                              ncpu = local_cpus(),
                              hitlist_method = c("blastclust", "usearch"),
                              usearch = Sys.which("usearch"),
                              ...) {
  stopifnot(is.null(threshold_name) || length(threshold) == length(threshold_name))
  hitlist_method <- match.arg(hitlist_method, several.ok = FALSE)
  is_file <- length(seq) == 1 && file.exists(seq)
  # if we only have one sequence, nothing to do.
  if (isFALSE(is_file)) {
    stopifnot(length(seq) == length(seq_id))
    s <- set_names(seq_id, seq_id)[which]
  }
  if (isFALSE(is_file) && length(s) <= 1) {
    out <- rep(list(names(s)), length(threshold))
  } else if (length(threshold) == 1) {
    # if there is only one threshold, then no need to cache the hit table
    out <- list(
      blastclust(
        seq = seq,
        threshold = threshold,
        seq_id = seq_id,
        which = which,
        ncpu = ncpu
      )
    )
  } else {
    # normal case: multiple sequences, multiple thresholds
    out <- list()
    outfiles <- tempfile("out")
    on.exit(unlink(outfiles, force = TRUE), TRUE)
    hits <- tempfile("hits")
    on.exit(unlink(hits, force = TRUE))
    if (hitlist_method == "blastclust") {
      unpruned_hits <- tempfile("unpruned")
      on.exit(unlink(unpruned_hits, force = TRUE), TRUE)
      out <- list(
        blastclust(seq,
                   threshold = threshold[1],
                   seq_id = seq_id,
                   which = which,
                   ncpu = ncpu,
                   outfile = outfiles[1],
                   hits = unpruned_hits
        )
      )
      prune_hitlist(unpruned_hits, hits)
      unlink(unpruned_hits)
    } else if (hitlist_method == "usearch") {
      cat(sprintf("%s Generating hitlist for %d queries.\n", Sys.time(), length(seq[which])))
      usearch_hitlist(
        seq,
        threshold = threshold[1]/100,
        seq_id = seq_id,
        which = which,
        ncpu = ncpu,
        hits = hits,
        usearch = usearch
      )
      out <- list(
        blastclust_reclust(
          hits,
          threshold = threshold[1],
          ncpu = ncpu,
          outfile = outfiles[1]
        )
      )
    } else {
      stop("unknown hitlist method: ", hitlist_method)
    }
    for (i in 2:length(threshold)) {
      outfiles[i] <- tempfile("out")
      out[[i]] <-
        blastclust_reclust(
          hits = hits,
          threshold = threshold[i],
          preclusters = outfiles[i-1],
          ncpu = ncpu,
          outfile = outfiles[i]
        )
    }
  }
  if (!is.null(threshold_name)) names(out) <- threshold_name
  out
}
