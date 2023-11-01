#### Functions which call external software from R
# Brendan Furneaux 2022

# try to find the vsearch executable
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

# try to find the cutadapt executable
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

# try to find the hmmer executable
find_hmmer <- function() {
  hmmalign <- Sys.getenv("hmmalign")
  if (nchar(hmmalign) == 0 || !file.exists(hmmalign)) {
    hmmalign <- Sys.which("hmmalign")
  }
  if (nchar(hmmalign) == 0 || !file.exists(hmmalign)) {
    stop("cannot find hmmer")
  }
  hmmalign
}

#' "usearch_global" function of vsearch
#'
#' @param query (`data.frame`, `Biostrings::DNAStringSet`, or `character` vector) query
#' sequences
#' @param ref (`data.frame`, `Biostrings::DNAStringSet`, `character` vector, or
#' file name) reference sequences
#' @param threshold (`numeric` scalar) identity threshold, in range 0.0-1.0
#' @param global (`logical` flag) if `TRUE`, end gaps and internal gaps are
#' penalized equally.  Otherwise end gaps are not penalized.
#' @param ncpu (`integer` count) number of threads to use
#'
#' @return `tibble::tibble` with columns `seq_id` and `clust`, where `seq_id` is
#' the name of a sequence from `query`, and `clust` is the closest match to that
#' sequence in `ref`
#' @export
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
  map <- tibble::tibble(
    seq_id_in = seq_len(ncol(seqtab)),
    seq_id_out = seq_len(ncol(seqtab))
  )
  if (nrow(matches) > 0) {
    matches$query <- as.integer(matches$query)
    matches$hit <- as.integer(matches$hit)
    matches <- matches[order(matches$query),]
    for (i in unique(matches$hit)) {
      seqtab[,i] <- seqtab[,i] +
        as.integer(rowSums(seqtab[,matches$query[matches$hit == i], drop = FALSE]))
    }
    seqtab <- seqtab[,-matches$query]
    map$seq_id_out[matches$query] <- matches$hit
    map$seq_id_out = map$seq_id_out - findInterval(map$seq_id_out, matches$query)
  }
  attr(seqtab, "map") <- map
  return(seqtab)
}

cutadapt_paired_options <- function(
    max_err = 0.2,
    min_overlap = 10L,
    action = "retain",
    discard_untrimmed = TRUE,
    max_n = 0L,
    max_ee = NULL,
    min_length = NULL, max_length = NULL,
    truncQ_R1 = 2, truncQ_R2 = 2,
    cut_R1 = NULL, cut_R2 = NULL
) {
  checkmate::assert_number(max_err, lower = 0, null.ok = TRUE)
  checkmate::assert_count(min_overlap, positive = TRUE, null.ok = TRUE)
  checkmate::assert_choice(
    action,
    c("trim", "retain", "mask", "lowercase", "none"),
    null.ok = TRUE
  )
  checkmate::assert_flag(discard_untrimmed)
  checkmate::assert_count(max_n, null.ok = TRUE)
  checkmate::assert_number(max_ee, lower = 0, finite = TRUE, null.ok = TRUE)
  checkmate::assert_count(min_length, positive = TRUE, null.ok = TRUE)
  checkmate::assert_integerish(truncQ_R1, min.len = 1, max.len = 2, null.ok = TRUE)
  checkmate::assert_integerish(truncQ_R2, min.len = 1, max.len = 2, null.ok = TRUE)
  checkmate::assert_integerish(cut_R1, min.len = 1, max.len = 2, null.ok = TRUE)
  checkmate::assert_integerish(cut_R2, min.len = 1, max.len = 2, null.ok = TRUE)
  structure(
    list(
      max_err = max_err,
      min_overlap = min_overlap,
      action = action,
      discard_untrimmed = discard_untrimmed,
      max_n = max_n,
      max_ee = max_ee,
      min_length = min_length,
      truncQ_R1 = truncQ_R1,
      truncQ_R2 = truncQ_R2,
      cut_R1 = cut_R1,
      cut_R2 = cut_R2
    ),
    class = "cutadapt_paired_options"
  )
}

cutadapt_paired_filter_trim <- function(
  file_R1, file_R2,
  primer_R1, primer_R2,
  trim_R1, trim_R2,
  options = cutadapt_paired_options(),
  ncpu = local_cpus(),
  cutadapt = find_cutadapt(),
  ...
) {
  checkmate::assert_class(options, "cutadapt_paired_options")
  args <- c(
    "--quiet",
    "-g", primer_R1,
    "-G", primer_R2,
    "-o", trim_R1,
    "-p", trim_R2
  )
  if (!is.null(options$max_err)) {
    args <- c(args, "-e", options$max_err)
  }
  if (!is.null(options$min_overlap)) {
    args <- c(args, "-O", options$min_overlap)
  }
  if (!is.null(options$action)) {
    args <- c(args, paste0("--action=", options$action))
  }
  if (isTRUE(options$discard_untrimmed)) {
    args <- c(args, "--discard-untrimmed")
  }
  if (!is.null(options$max_n)) {
    args <- c(args, "--max-n", options$max_n)
  }
  if (!is.null(options$max_ee)) {
    args <- c(args, "--max-ee", options$max_ee)
  }
  if (!is.null(options$min_length)) {
    args <- c(args, "-m", options$min_length)
  }
  if (!is.null(options$max_length)) {
    args <- c(args, "-M", options$max_length)
  }
  if (!is.null(options$truncQ_R1)) {
    args <- c(args, "-q", paste(round(options$truncQ_R1), collapse = ","))
  }
  if (!is.null(options$truncQ_R2)) {
    args <- c(args, "-Q", paste(round(options$truncQ_R2), collapse = ","))
  }
  if (!is.null(options$cut_R1)) {
    args <- c(args, "-u", paste(round(options$cut_R1), collapse = ","))
  }
  if (!is.null(options$cut_R2)) {
    args <- c(args, "-U", paste(round(options$cut_R2), collapse = ","))
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

cutadapt_options <- function(
    max_err = 0.2,
    min_overlap = 10L,
    action = "retain",
    discard_untrimmed = TRUE,
    max_n = 0L,
    max_ee = NULL,
    min_length = NULL, max_length = NULL,
    truncQ = NULL,
    cut = NULL
) {
  checkmate::assert_number(max_err, lower = 0, null.ok = TRUE)
  checkmate::assert_count(min_overlap, positive = TRUE, null.ok = TRUE)
  checkmate::assert_choice(
    action,
    c("trim", "retain", "mask", "lowercase", "none"),
    null.ok = TRUE
  )
  checkmate::assert_flag(discard_untrimmed)
  checkmate::assert_count(max_n, null.ok = TRUE)
  checkmate::assert_number(max_ee, lower = 0, finite = TRUE, null.ok = TRUE)
  checkmate::assert_count(min_length, positive = TRUE, null.ok = TRUE)
  checkmate::assert_integerish(truncQ, min.len = 1, max.len = 2, null.ok = TRUE)
  checkmate::assert_integerish(cut, min.len = 1, max.len = 2, null.ok = TRUE)
  structure(
    list(
      max_err = max_err,
      min_overlap = min_overlap,
      action = action,
      discard_untrimmed = discard_untrimmed,
      max_n = max_n,
      max_ee = max_ee,
      min_length = min_length,
      truncQ = truncQ,
      cut = cut
    ),
    class = "cutadapt_options"
  )
}

cutadapt_filter_trim <- function(
  file,
  primer,
  trim,
  options = cutadapt_options(),
  ncpu = local_cpus(),
  cutadapt = find_cutadapt(),
  ...
) {
  args <- c(
    "-g", primer,
    "-o", trim
  )
  if (!is.null(max_err)) {
    args <- c(args, "-e", max_err)
  }
  if (!is.null(min_overlap)) {
    args <- c(args, "-O", round(min_overlap))
  }
  if (!is.null(action)) {
    args <- c(args, paste0("--action=", action))
  }
  if (isTRUE(discard_untrimmed)) {
    args <- c(args, "--discard-untrimmed")
  }
  if (!is.null(max_n)) {
    args <- c(args, "--max-n", max_n)
  }
  if (!is.null(max_ee)) {
    args <- c(args, "--max-ee", max_ee)
  }
  if (!is.null(min_length)) {
    args <- c(args, "-m", min_length)
  }
  if (!is.null(max_length)) {
    args <- c(args, "-M", max_length)
  }
  if (!is.null(truncQ)) {
    args <- c(args, "-q", paste(round(truncQ), collapse = ","))
  }
  if (!is.null(cut)) {
    args <- c(args, "-u", paste(round(cut), collapse = ","))
  }
  if (!is.null(ncpu)) {
    checkmate::assert_count(ncpu, positive = TRUE)
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

fastq_names <- function(fq) {
  if (!file.exists(fq)) return(character())
  if (endsWith(fq, ".gz")) {
    system(paste("zcat", fq, "| awk 'NR%4==1{print substr($1, 2)}'"), intern = TRUE)
  } else {
    system(paste(" awk 'NR%4==1{print substr($1, 2)}'", fq), intern = TRUE)
  }
}
