#### Functions which call external software from R
# Brendan Furneaux 2022


#### VSEARCH ####

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
      "| awk '$1==\"H\" {print $9,$10,$4}'"
    ),
    intern = TRUE
  )
  stopifnot(attr(uc, "status") == 0)
  if (length(uc) > 0) {
    readr::read_delim(
      I(uc),
      col_names = c("seq_id", "cluster", "dist"),
      delim = " ",
      col_types = "ccn"
    ) |>
      dplyr::mutate(dist = 1 - dist/100)
  } else {
    tibble::tibble(seq_id = character(), cluster = character(), dist = numeric())
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
    result <- vsearch_usearch_global(query, ref, threshold, ...)[,c("seq_id", "cluster")]
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

#### Cutadapt ####

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

#### ProtaxFungi ####

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
