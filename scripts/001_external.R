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

vsearch_usearch_global <- function(query, ref, threshold, ncpu = local_cpus()) {
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
      "--gapopen", "1",
      "--gapext", "1",
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
      col_names = c("ASV", "cluster"),
      delim = " ",
      col_types = "cc"
    )
  } else {
    tibble::tibble(ASV = character(), cluster = character())
  }
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
        rowSums(seqtab[,matches$query[matches$hit == i], drop = FALSE])
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

run_protax <- function(seqs, outdir, ncpu = local_cpus()) {
  if (dir.exists(outdir)) unlink(outdir, recursive = TRUE)
  dir.create(outdir)
  write_sequence(seqs, file.path(outdir, "all.fa"))
  status <- system2(
    "scripts/runprotax",
    c(outdir, ncpu)
  )
  stopifnot(status == 0)
  list.files(outdir, full.names = TRUE)
}

#' Cluster ASVs using blastclust
#'
#' @param seqs (`character` vector, file name, or `Biostrings::DNAStringSet`)
#' sequences to cluster
#' @param threshold (`numeric` scalar) percentage similarity threshold for
#' clustering. A number between 3 and 100
#' @param seqnames (`character`) names for the sequences (Default: `names(seqs)`)
#' @param which (`logical`, `integer`, or `character`) subset indices
#' indicating which sequences in `seqs` should be clustered. Has no effect if
#' `seqs` is a file name. (Default: `TRUE`, i.e. all the sequences)
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
blastclust <- function(seqs, threshold, seqnames = names(seqs), which = TRUE,
                       ncpu = local_cpus(), outfile = NULL, hits = NULL) {
  UseMethod("blastclust", seqs)
}

blastclust.character <- function(seqs, threshold, seqnames = names(seqs),
                                 which = TRUE, ncpu = local_cpus(),
                                 outfile = NULL, hits = NULL) {
  if (length(seqs) == 1 && file.exists(seqs)) {
    if (!missing(seqnames))
      warning("'seqnames' has no effect when 'seqs' is a file.")
    if (!missing(which))
      warning("'which' has no effect when 'seqs' is a file.")
    blastclust_filename(seqs, threshold, ncpu)
  }
  seqs <- Biostrings::DNAStringSet(seqs)
  blastclust.DNAStringSet(seqs, threshold, seqnames, which, ncpu, outfile, hits)
}

blastclust.DNAStringSet <- function(seqs, threshold, seqnames = names(seqs),
                                    which = TRUE,
                                    ncpu = local_cpus(),
                                    outfile = NULL, hits = NULL) {
  # rename the sequences if necessary
  if (!isTRUE(all.equal(names(seqs), seqnames))) names(seqs) <- seqnames
  seqs <- seqs[which]
  # shortcut if only one sequence
  if (length(seqs) == 1) return(names(seqs))
  tf <- tempfile(pattern = "clust", fileext = ".fasta")
  Biostrings::writeXStringSet(seqs, tf)
  on.exit(unlink(tf))
  blastclust_filename(tf, threshold, ncpu, outfile, hits)
}

blastclust_filename <- function(seqs, threshold, ncpu, outfile = NULL, hits = NULL) {
  if (is.null(outfile)) {
    outfile <- tempfile(pattern = "out")
    on.exit(unlink(outfile), TRUE)
  }
  stopifnot(
    system2(
      "blastclust",
      c(
        "-i", seqs, # input file
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
        if (!is.null(preclusters)) c("-l", preclusters) else NULL,
        "-S", threshold, # similarity threshold
        "-L", "0.95", # overlap threshold
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

usearch_hitlist <- function(seqs, threshold, seqnames = names(seqs), which = TRUE,
                            ncpu = local_cpus(), hits = NULL, usearch = Sys.which("usearch")) {
  UseMethod("usearch_hitlist", seqs)
}

usearch_hitlist.character <- function(seqs, threshold, seqnames = names(seqs),
                                      which = TRUE, ncpu = local_cpus(), hits = NULL, usearch = Sys.which("usearch")) {
  if (length(seqs) == 1 && file.exists(seqs)) {
    if (!missing(seqnames))
      warning("'seqnames' has no effect when 'seqs' is a file.")
    if (!missing(which))
      warning("'which' has no effect when 'seqs' is a file.")
    index <- Biostrings::fasta.seqlengths(seqs)
    do_usearch_hitlist(
      seqs,
      seqlen = index,
      names = names(index),
      threshold = threshold,
      ncpu = ncpu,
      hits = hits,
      usearch = usearch
    )
  } else {
    seqs <- Biostrings::DNAStringSet(seqs)
    usearch_hitlist.DNAStringSet(seqs, threshold, seqnames, which, ncpu,
                                 hits, usearch = usearch)
  }
}

usearch_hitlist.DNAStringSet <- function(seqs, threshold, seqnames = names(seqs),
                                         which = TRUE,
                                         ncpu = local_cpus(),
                                         hits = NULL,
                                         usearch = Sys.which("usearch")) {
  # rename the sequences if necessary
  if (!isTRUE(all.equal(names(seqs), seqnames))) names(seqs) <- seqnames
  seqs <- seqs[which]
  # shortcut if only one sequence
  if (length(seqs) == 1) return(names(seqs))
  tf <- tempfile(pattern = "clust", fileext = ".fasta")
  Biostrings::writeXStringSet(seqs, tf)
  on.exit(unlink(tf))
  do_usearch_hitlist(tf, seqlen = Biostrings::nchar(seqs), names = names(seqs),
                     threshold = threshold, ncpu = ncpu, hits = hits,
                     usearch = usearch)
}

do_usearch_hitlist <- function(seqs, seqlen, names, threshold, ncpu, hits,
                               usearch = Sys.which("usearch")) {
  if (!methods::is(hits, "connection")) {
    hits <- file(hits, open = "wb")
  }
  if (!isOpen(hits)) open(hits, "wb")
  
  on.exit(close(hits), TRUE)
  # list type 0 (names)
  # list size (characters in names list)
  writeBin(c(0L, sum(nchar(names)) + length(names)), hits, size = 4L)
  # list of names
  writeChar(paste0(names, " ", collapse = ""), hits, eos = NULL)
  # sequence lengths
  writeBin(seqlen, hits, size = 4L)
  seqidx <- seq_along(names) - 1L
  names(seqidx) <- names
  if (is.null(names(seqlen))) names(seqlen) <- names
  fifoname <- tempfile("fifo")
  stopifnot(system2("mkfifo", fifoname) == 0)
  on.exit(unlink(fifoname), TRUE)
  f <- fifo(fifoname)
  system2(
    usearch,
    c(
      "-calc_distmx", seqs, # input file
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
  d <- readLines(f, n = 100000)
  while (length(d) > 0) {
    d <- strsplit(d, "\t", fixed = TRUE)
    d <- unlist(d)
    d <- matrix(d, ncol = 3, byrow = TRUE)
    d <- as.data.frame(d)
    seq1 <- matrix(writeBin(seqidx[d$V1], raw(), size = 4), nrow = 4)
    seq2 <- matrix(writeBin(seqidx[d$V2], raw(), size = 4), nrow = 4)
    hsp1 <- seqlen[d$V1]
    hsp2 <- seqlen[d$V2]
    ident <- 1 - as.numeric(d$V3)
    score <- matrix(writeBin(ident * pmin(hsp1, hsp2), raw(), size = 8), nrow = 8)
    hsp1 <- matrix(writeBin(hsp1, raw(), size = 4), nrow = 4)
    hsp2 <- matrix(writeBin(hsp2, raw(), size = 4), nrow = 4)
    ident <- matrix(writeBin(100*ident, raw(), size = 8), nrow = 8)
    d <- c(rbind(seq1, seq2, hsp1, hsp2, score, ident))
    writeBin(d, hits)
    d <- readLines(f, n = 100000)
  }
}

do_usearch_hitlist2 <- function(seqs, seqlen, names, threshold, ncpu, hits,
                               usearch = Sys.which("usearch")) {
  if (!methods::is(hits, "connection")) {
    hits <- file(hits, open = "wb")
  }
  if (!isOpen(hits)) open(hits, "wb")
  
  on.exit(close(hits), TRUE)
  # list type 0 (names)
  # list size (characters in names list)
  writeBin(c(0L, sum(nchar(names)) + length(names)), hits, size = 4L)
  # list of names
  writeChar(paste0(names, " ", collapse = ""), hits, eos = NULL)
  # sequence lengths
  writeBin(seqlen, hits, size = 4L)
  seqidx <- seq_along(names) - 1L
  names(seqidx) <- names
  if (is.null(names(seqlen))) names(seqlen) <- names
  fifoname <- tempfile("fifo")
  stopifnot(system2("mkfifo", fifoname) == 0)
  on.exit(unlink(fifoname), TRUE)
  f <- fifo(fifoname)
  system2(
    usearch,
    c(
      "-calc_distmx", seqs, # input file
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
  single_linkage(fifoname)
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
blastclust_repeat <- function(seqs, threshold, seqnames = names(seqs),
                              which = TRUE,
                              threshold_name = names(threshold),
                              ncpu = local_cpus(),
                              hitlist_method = c("blastclust", "usearch"),
                              usearch = Sys.which("usearch"),
                              ...) {
  stopifnot(is.null(threshold_name) || length(threshold) == length(threshold_name))
  hitlist_method <- match.arg(hitlist_method, several.ok = FALSE)
  is_file <- length(seqs) == 1 && file.exists(seqs)
  # if we only have one sequence, nothing to do.
  if (isFALSE(is_file)) {
    stopifnot(length(seqs) == length(seqnames))
    s <- set_names(seqnames, seqnames)[which]
  }
  if (isFALSE(is_file) && length(s) <= 1) {
    out <- rep(list(names(s)), length(threshold))
  } else if (length(threshold) == 1) {
    # if there is only one threshold, then no need to cache the hit table
    out <- list(
      blastclust(
        seqs = seqs,
        threshold = threshold,
        seqnames = seqnames,
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
        blastclust(seqs,
                   threshold = threshold[1],
                   seqnames = seqnames,
                   which = which,
                   ncpu = ncpu,
                   outfile = outfiles[1],
                   hits = unpruned_hits
        )
      )
      prune_hitlist(unpruned_hits, hits)
      unlink(unpruned_hits)
    } else if (hitlist_method == "usearch") {
      cat(sprintf("%s Generating hitlist for %d queries.\n", Sys.time(), length(seqs[which])))
      usearch_hitlist(
        seqs,
        threshold = threshold[1]/100,
        seqnames = seqnames,
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

Rcpp::sourceCpp("src/single_linkage.cpp")