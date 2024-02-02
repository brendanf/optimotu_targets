#### Functions which call external software from R
# Brendan Furneaux 2022

find_executable <- function(executable) {
  checkmate::assert_character(executable)
  out <- Sys.getenv(executable)
  if (nchar(out) == 0 || !file.exists(out)) {
    out <- Sys.getenv(toupper(executable))
  }
  if (nchar(out) == 0 || !file.exists(out)) {
    out <- Sys.which(executable)
  }
  if (nchar(out) == 0 || !file.exists(out)) {
    out <- list.files(path = "bin", pattern = executable, recursive = TRUE, full.names = TRUE)
  }
  checkmate::assert_file_exists(out, access = "x", .var.name = executable)
  out
}

# try to find the vsearch executable
find_vsearch <- function() {
  vsearch <- Sys.getenv("VSEARCH")
  if (nchar(vsearch) == 0 || !file.exists(vsearch)) {
    vsearch <- Sys.which("vsearch")
  }
  if (nchar(vsearch) == 0 || !file.exists(vsearch)) {
    vsearch <- list.files(path = "bin", pattern = "vsearch", recursive = TRUE, full.names = TRUE)
  }
  checkmate::assert_file_exists(vsearch, access = "x")
  vsearch
}

# try to find the cutadapt executable
find_cutadapt <- function() {
  cutadapt <- Sys.getenv("CUTADAPT")
  if (nchar(cutadapt) == 0 || !file.exists(cutadapt)) {
    cutadapt <- Sys.which("cutadapt")
  }
  if (nchar(cutadapt) == 0 || !file.exists(cutadapt)) {
    cutadapt <- list.files(path = "bin", pattern = "cutadapt", recursive = TRUE, full.names = TRUE)
  }
  checkmate::assert_file_exists(cutadapt, access = "x")
  cutadapt
}

# try to find the hmmer executable
find_hmmalign <- function() {
  find_executable("hmmalign")
}

# try to find the hmmsearch executable
find_hmmsearch <- function() {
  find_executable("hmmsearch")
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

vsearch_uchime_ref <- function(query, ref, ncpu = local_cpus(), id_only = FALSE, id_is_int = FALSE) {
  checkmate::assert_flag(id_is_int)
  if (checkmate::test_file_exists(query, "r")) {
    tquery <- query
  } else {
    tquery <- tempfile("query", fileext = ".fasta")
    on.exit(unlink(c(tquery), force = TRUE))
    write_sequence(query, tquery)
  }
  if (checkmate::test_file_exists(ref, "r")) {
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
    args = c(
      "--uchime_ref", tquery,
      "--db", tref,
      "--chimeras", tchimeras,
      "--threads", ncpu
    )
  )
  stopifnot(vs == 0L)
  if (id_only) {
    out <- names(Biostrings::fasta.seqlengths(tchimeras))
    if (id_is_int) {
      as.integer(out)
    } else {
      out
    }
  } else {
    out <- Biostrings::readDNAStringSet(tchimeras) %>%
    as.character() %>%
    tibble::enframe(name = "seq_id", value = "seq")
    if (id_is_int) {
      dplyr::transmute(out, seq_idx = as.integer(seq_id), seq)
    } else {
      out
    }
  }
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

collapseNoMismatch_vsearch <- function(seqtab, ..., ncpu = local_cpus()) {
  UseMethod("collapseNoMismatch_vsearch", seqtab)
}

collapseNoMismatch_vsearch.matrix <- function(seqtab, ..., ncpu = local_cpus()) {
  seqs <- colnames(seqtab)
  names(seqs) <- seq_along(seqs)
  matches <- vsearch_cluster_smallmem(seqs, ncpu = ncpu)
  map <- tibble::tibble(
    seq_idx_in = seq_len(ncol(seqtab)),
    seq_idx_out = seq_len(ncol(seqtab))
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
    map$seq_ixd_out[matches$query] <- matches$hit
    map$seq_idx_out = map$seq_idx_out - findInterval(map$seq_idx_out, matches$query)
  }
  attr(seqtab, "map") <- map
  return(seqtab)
}

collapseNoMismatch_vsearch.data.frame <- function(seqtab, seqs = NULL, abund_col = "nread", ..., ncpu = local_cpus()) {
  if (is.null(seqs)) {
    checkmate::assert_names(
      names(seqtab),
      must.include = c("seq", abund_col),
      disjunct.from = "seq_idx"
    )
    seqs <- unique(seqtab$seqs)
  } else {
    checkmate::assert_names(
      names(seqtab),
      must.include = c("seq_idx", abund_col),
      disjunct.from = "seq"
    )
  }
  names(seqs) <- as.character(seq_along(seqs))
  matches <- vsearch_cluster_smallmem(seqs, ncpu = ncpu)
  map <- tibble::tibble(
    seq_idx_in = seq_along(seqs),
    seq_idx_out = seq_along(seqs)
  )
  if (nrow(matches) > 0) {
    matches$query <- as.integer(matches$query)
    matches$hit <- as.integer(matches$hit)
    matches <- matches[order(matches$query),]
    for (i in unique(matches$hit)) {
      if (seq_idx %in% names(seqtab)) {
        seqtab$seq_idx <- ifelse(seqtab$seq_idx == matches$query[i], matches$hit[i], seqtab$seq_idx)

      } else {
        seqtab$seq <- ifelse(seqtab$seq == seqs[matches$query[i]], seqs[matches$hit[i]], seqtab$seq)
      }
    }
    seqtab <- dplyr::summarize(
      seqtab,
      dplyr::across(all_of(abund_col), sum),
      .by = c(sample, any_of("seq_idx", "seq"))
    )
    map$seq_idx_out[matches$query] <- matches$hit
    map$seq_idx_out = map$seq_idx_out - findInterval(map$seq_idx_out, matches$query)
  }
  attr(seqtab, "map") <- map
  seqtab
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

hmmalign <- function(seqs, hmm, outfile, outformat = "A2M", compress = endsWith(outfile, ".gz")) {
  checkmate::assert_string(hmm)
  checkmate::assert_file_exists(hmm, access = "r")
  ensure_directory(outfile)
  checkmate::assert_path_for_output(outfile, overwrite = TRUE)
  checkmate::assert_choice(outformat, c("A2M", "a2m", "afa", "AFA"))
  checkmate::assert_flag(compress)
  exec <- find_hmmalign()
  checkmate::assert_file_exists(exec, access = "x")
  if (checkmate::test_file_exists(seqs, "r")) {
    tseqs <- seqs
    n <- length(seqs)
  } else if (checkmate::test_list(seqs, types = c("character", "XStringSet", "data.frame"))) {
    n <- length(seqs)
    tseqs <- replicate(n, withr::local_tempfile(fileext = ".fasta"))
    purrr::pwalk(list(seq = seqs, fname = tseqs), write_sequence)
  } else {
    checkmate::assert_multi_class(seqs, c("data.frame", "character", "XStringSet"))
    n <- 1
    tseqs <- withr::local_tempfile(fileext = ".fasta")
    write_sequence(seqs, tseqs)
  }
  checkmate::assert(
    length(outfile) == 1,
    length(outfile) == n
  )
  if (length(outfile) == 1 && n > 1) {
    tout <- replicate(n, withr::local_tempfile(fileext = ".fasta"))
  } else {
    tout <- outfile
  }
  if (compress) {
    if (identical(tout, outfile)) {
      tout <- replicate(n, withr::local_tempfile(fileext = ".fasta"))
    }
  }
  mout <- replicate(n, withr::local_tempfile(fileext = ".fasta"))
  for (i in seq_len(n)) {
    processx::run("mkfifo", mout[i])
  }

  args <- data.frame(
    "--outformat", outformat,
    "--trim",
    "-o", mout,
    hmm,
    tseqs
  )
  args <- as.matrix(args)
  hmmer <- vector("list", n)
  deline <- vector("list", n)
  for (i in seq_len(n)) {
    hmmer[[i]] <- processx::process$new(
      command = exec,
      args = args[i,],
      supervise = TRUE
    )
    deline[[i]] <- processx::process$new(
      command = "awk",
      args = 'BEGIN{ORS=""};NR>1&&/^>/{print "\\n"};{print};/^>/{print "\\n"};END{print "\\n"}',
      stdin = mout[i],
      stdout = tout[i],
      supervise = TRUE
    )
  }
  hmmer_return <- integer()
  for (i in seq_len(n)) {
    hmmer[[i]]$wait()
    hmmer_return <- union(hmmer[[i]]$get_exit_status(), hmmer_return)
  }
  stopifnot(identical(hmmer_return, 0L))
  for (i in seq_len(n)) {
    deline[[i]]$wait()
  }

  if (compress && length(outfile) == n) {
    gzip <- vector("list", n)
    for (i in seq_len(n)) {
      gzip[[i]] <- processx::process$new(
        command = "gzip",
        error_on_status = TRUE,
        args = c("-c", tout[i]),
        stdout = outfile[i]
      )
    }
    gzip_return <- integer()
    for (i in seq_len(n)) {
      gzip_return <- union(gzip[[i]]$wait()$status, gzip_return)
    }
    stopifnot(identical(gzip_return, 0L))
  } else if (length(outfile) < n) {
    fastx_combine(tout, outfile)
  }
  outfile
}

read_domtblout <- function(file) {
  tibble::tibble(
    text = readLines(file),
    is_widths = grepl("^#[- ]+$", text),
    part = cumsum(is_widths)
  ) |>
    dplyr::filter(is_widths | !startsWith(text, "#")) |>
    dplyr::group_split(part, .keep = FALSE) |>
    purrr::discard(\(x) nrow(x) == 1) |>
    purrr::map_dfr(
      \(x) {
        paste(x$text, collapse = "\n") |>
          readr::read_fwf(
            col_positions =  stringr::str_locate_all(x$text[1], "#?-+")[[1]] |>
              tibble::as_tibble() |>
              tibble::add_column(
                col_names = c("seq_name", "seq_accno", "seq_length", "hmm_name", "hmm_accno", "hmm_length", "Evalue", "full_score", "full_bias", "hit_num", "total_hits", "c_Evalue", "i_Evalue", "hit_score", "hit_bias", "hmm_from", "hmm_to", "seq_from", "seq_to", "env_from", "env_to", "acc", "description")
              ) |>
              do.call(readr::fwf_positions, args = _),
            skip = 1,
            col_types = "cciccidddiiddddiiiiiidc"
          )
      }
    )

}

hmmsearch <- function(seqs, hmm) {
  checkmate::assert_string(hmm)
  checkmate::assert_file_exists(hmm, access = "r")
  exec <- find_hmmsearch()
  checkmate::assert_file_exists(exec, access = "x")
  if (checkmate::test_file_exists(seqs, "r")) {
    tseqs <- seqs
    n <- length(seqs)
  } else if (checkmate::test_list(seqs, types = c("character", "XStringSet", "data.frame"))) {
    n <- length(seqs)
    tseqs <- replicate(n, withr::local_tempfile(fileext = ".fasta"))
    purrr::pwalk(list(seq = seqs, fname = tseqs), write_sequence)
  } else {
    checkmate::assert_multi_class(seqs, c("data.frame", "character", "XStringSet"))
    n <- 1
    tseqs <- withr::local_tempfile(fileext = ".fasta")
    write_sequence(seqs, tseqs)
  }
  outfile <- replicate(n, tempfile(fileext = ".hmmout"))
  args <- data.frame(
    "--noali",
    "--notextw",
    "--domtblout", outfile,
    hmm,
    tseqs
  )
  args <- as.matrix(args)

  hmmer <- vector("list", n)
  for (i in seq_len(n)) {
    hmmer[[i]] <- processx::process$new(
      command = exec,
      args = args[i,],
      supervise = TRUE
    )
  }
  hmmer_return <- integer()
  for (i in seq_len(n)) {
    hmmer[[i]]$wait()
    hmmer_return <- union(hmmer[[i]]$get_exit_status(), hmmer_return)
    stopifnot(identical(hmmer_return, 0L))
  }
  purrr::map_dfr(
    outfile,
    read_domtblout
  )
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

# vectorized on aln_seqs
# attempts to run them ALL in parallel, be careful!
run_protax_animal <- function(protax, aln_seqs, modeldir, min_p = 0.1, strip_inserts = TRUE) {
  checkmate::assert_file_exists(protax, access = "x")
  checkmate::assert_file_exists(aln_seqs, access = "r")
  checkmate::assert_directory_exists(modeldir)
  priors <- file.path(modeldir, "taxonomy.priors")
  checkmate::assert_file_exists(priors, access = "r")
  refs <- file.path(modeldir, "refs.aln")
  checkmate::assert_file_exists(refs, access = "r")
  rseqs <- file.path(modeldir, "model.rseqs.numeric")
  checkmate::assert_file_exists(rseqs, access = "r")
  pars <- file.path(modeldir, "model.pars")
  checkmate::assert_file_exists(pars, access = "r")
  scs <- file.path(modeldir, "model.scs")
  checkmate::assert_file_exists(scs, access = "r")
  checkmate::check_number(min_p, lower = 0, upper = 1, finite = TRUE)
  checkmate::assert_flag(strip_inserts)
  args <- c(priors, refs, rseqs, pars, scs, as.character(min_p))

  is_gz <-endsWith(aln_seqs, ".gz")
  stopifnot(all(is_gz) | all(!is_gz))
  is_gz <- all(is_gz)

  n <- length(aln_seqs)
  if (strip_inserts | is_gz) {
    prepipe <- vector("list", n)
    protax_in <- replicate(n, withr::local_tempfile(fileext = ".fasta"))
    pipecommand <- if (is_gz) "zcat" else "cat"
    pipecommand <- paste(pipecommand, aln_seqs)
    if (strip_inserts) pipecommand <- paste(pipecommand, "| tr -d 'acgt'")
    pipecommand <- paste(pipecommand, ">", protax_in)
  } else {
    protax_in <- aln_seqs
  }

  protax_process <- vector("list", n)
  outfiles <- replicate(n, withr::local_tempfile())
  for (i in seq_len(n)) {
    if (strip_inserts | is_gz) {
      system(pipecommand[i])
    }
    protax_process[[i]] <- processx::process$new(
      command = protax,
      args = c(args, protax_in[i]),
      stdout = outfiles[i]
    )
  }
  protax_exit_status = 0L
  output <- vector("list", n)
  for (i in seq_len(n)) {
    protax_process[[i]]$wait()
    protax_exit_status <- max(protax_exit_status, protax_process[[i]]$get_exit_status())
    stopifnot(protax_exit_status == 0L)
    output[[i]] <- readLines(outfiles[i])
  }
  unlist(output)
}

parse_protaxAnimal_output <- function(x) {
  out <-
    tibble::tibble(output = x) |>
    tidyr::separate(
      output,
      into = c("seq_id", "assignment"),
      extra = "merge",
      fill = "right",
      convert = TRUE
    ) |>
    dplyr::mutate(
      assignment = gsub("([^ ]+) ([0-9.]+)", "\\1\x1f\\2", assignment)
    ) |>
    tidyr::separate_longer_delim(assignment, " ") |>
    tidyr::separate(
      assignment,
      into = c("taxonomy", "prob"),
      sep = "\x1f",
      convert = TRUE
    ) |>
    dplyr::mutate(
      rank = factor(stringr::str_count(taxonomy, ",") + 1, labels = TAXRANKS)
    ) |>
    tidyr::replace_na(list(rank = TAXRANKS[1], prob = 1)) |>
    tidyr::extract(taxonomy, into = c("parent_taxonomy", "taxon"), regex = "(?:(.+),)?([^,]+)$") |>
    dplyr::mutate(
      parent_taxonomy = dplyr::na_if(parent_taxonomy, ""),
      taxon = dplyr::na_if(taxon, "unk")
    ) |>
    dplyr::select(seq_id, rank, parent_taxonomy, taxon, prob) |>
    dplyr::arrange(seq_id, rank)

  if (is.integer(out$seq_id)) out <- dplyr::rename(out, seq_idx = seq_id)
  out
}

parse_protaxA_info <- function(info) {
  if (length(info) == 1 && file.exists(info)) info <- readLines(info)
  trimws(info) |>
    tibble::tibble(
      text = _,
      entry_id = cumsum(text == "") + 1L
    ) |>
    dplyr::filter(text != "") |>
    dplyr::mutate(seq_id = text[1], .by = entry_id) |>
    tidyr::separate_wider_delim(
      text,
      delim = " ",
      names = c("rank", "taxonomy", "prob", NA, "best_id", "best_dist", NA,
                "second_id", "second_dist"),
      too_few = "align_start"
    )  |>
    dplyr::filter(!is.na(taxonomy) | dplyr::n() == 1, .by = entry_id) |>
    tidyr::extract(
      taxonomy,
      into = c("parent_taxonomy", "taxon"),
      regex = "(?:(.+),)?([^,]+)$"
    ) |>
    dplyr::mutate(
      parent_taxonomy = dplyr::na_if(parent_taxonomy, ""),
      rank = int2rank(ifelse(is.na(parent_taxonomy), 1L, as.integer(rank))),
      taxon = dplyr::na_if(taxon, "unk"),
      prob = as.numeric(prob),
      best_id = as.integer(best_id),
      best_dist = as.numeric(best_dist),
      second_id = as.integer(second_id),
      second_dist = as.numeric(second_dist)
    ) |>
    dplyr::select(seq_id, rank, parent_taxonomy, taxon, prob, best_id,
                  best_dist, second_id, second_dist)
}

fastq_names <- function(fq) {
  if (!file.exists(fq)) return(character())
  if (endsWith(fq, ".gz")) {
    system(paste("zcat", fq, "| awk 'NR%4==1{print substr($1, 2)}'"), intern = TRUE)
  } else {
    system(paste(" awk 'NR%4==1{print substr($1, 2)}'", fq), intern = TRUE)
  }
}


#' Quickly retrieve sequences from a gzipped file using an index
#'
#' @param file (`character` filename) file to create an index for
#'
#' @return file name of the created index
#' @rdname fastx_gz
fastx_gz_index <- function(file) {
  index <- sprintf("%s.fqi", file)
  args <- c(
    "index",
    sprintf("-f=%s", file),
    sprintf("-i=%s", index),
    "-w"
  )
  out = system2("bin/fastqindex_0.9.0b", args)
  stopifnot(out == 0)
  checkmate::assert_file_exists(index, "r")
  index
}

#' @param infile (`character` filename) gzipped fasta or fastq file
#' @param index (`character` filename) index file for `infile`
#' @param i (`integer` vector) indices to extract
#' @param outfile (`character` filename) file to write the extracted sequences to
#' @param renumber (`logical` flag) if `TRUE`, replace the sequence names with
#'   integers, starting at 0.
#' @param append (`logical` flag) if `TRUE`, append to `outfile` if it already
#'   exists, rather than overwriting.
#'
#' @return filename of the output file
#' @rdname fastx_gz
fastx_gz_extract <- function(infile, index, i, outfile, renumber = FALSE, append = FALSE, hash = NULL) {
  checkmate::assert_file_exists(infile, "r")
  checkmate::assert_file_exists(index, "r")
  checkmate::assert_integerish(i, lower = 1)
  checkmate::assert_string(outfile)
  checkmate::assert_flag(renumber)
  checkmate::assert_flag(append)
  if (file.exists(outfile) && !append) unlink(outfile)
  ensure_directory(outfile)
  if (!file.exists(outfile)) file.create(outfile)
  start <- which(i != dplyr::lag(i, 1, -1) + 1L)
  end <- c(start[-1] - 1L, length(i))
  is_fastq <- endsWith(infile, "fastq.gz") || endsWith(infile, "fq.gz")
  command <- sprintf(
    "bin/fastqindex_0.9.0b extract -s=%i -n=%i -e=%i -f=%s -i=%s | tail -n+8",
    i[start] - 1L,
    end - start + 1L,
    if (is_fastq) 4 else 2,
    infile,
    index
  )
  if (renumber) {
    command = paste(
      command,
      sprintf(
        "| awk -v n=%i 'NR%%%i==1{print \"%c\" n; n++; next}; {print}'",
        cumsum(dplyr::lag(end - start + 1L, 1L, 0L)),
        if (is_fastq) 4 else 2,
        if (is_fastq) "@" else ">"
      )
    )
  }
  if (endsWith(outfile, ".gz")) {
    command = paste(command, "| gzip -c -")
  }
  command = paste(command, ">>", outfile)
  result <- vapply(command, system, 0L)
  stopifnot(all(result == 0))
  outfile
}

fastx_gz_multi_extract <- function(infile, index, ilist, outfiles, renumber = FALSE, append = FALSE) {
  checkmate::assert_file_exists(infile, "r")
  checkmate::assert_file_exists(index, "r")
  checkmate::assert_list(i, types = "integerish")
  checkmate::assert_path_for_output(outfiles)
  stopifnot(length(ilist) == length(outfiles))
  checkmate::assert_flag(renumber)
  checkmate::assert_flag(append)
  for (of in outfiles) {
    if (file.exists(of) && !append) unlink(of)
    ensure_directory(of)
    if (!file.exists(of)) file.create(of)
  }
  start <- which(i != dplyr::lag(i, 1, -1) + 1L)
  end <- c(start[-1] - 1L, length(i))
  is_fastq <- endsWith(infile, "fastq.gz") || endsWith(infile, "fq.gz")
  command <- sprintf(
    "bin/fastqindex_0.9.0b extract -s=%i -n=%i -e=%i -f=%s -i=%s | tail -n+8",
    i[start] - 1L,
    end - start + 1L,
    if (is_fastq) 4 else 2,
    infile,
    index
  )
  if (renumber) {
    command = paste(
      command,
      sprintf(
        "| awk -v n=%i 'NR%%%i==1{print \"%c\" n; n++; next}; {print}'",
        cumsum(dplyr::lag(end - start + 1L, 1L, 0L)),
        if (is_fastq) 4 else 2,
        if (is_fastq) "@" else ">"
      )
    )
  }
  stopifnot(all(endsWith(outfiles, ".gz")) || all(!endsWith(outfiles, ".gz")))
  if (all(endsWith(outfiles, ".gz"))) {
    command = paste(command, "| gzip -c -")
  }
  command = paste(command, ">>", outfiles)
  result <- vapply(command, system, 0L)
  stopifnot(all(result == 0))
  outfiles
}


#' @param start (`integer` scalar) one-based index to start hashing
#' @param n (`integer` scalar) number of sequences to hash
#'
#' @return (`character`) md5 hash
fastx_gz_hash <- function(infile, index, start, n) {
  checkmate::assert_file_exists(infile, "r")
  checkmate::assert_file_exists(index, "r")
  checkmate::assert_integerish(start, lower = 1)
  checkmate::assert_integerish(n, lower = 1)
  is_fastq <- endsWith(infile, "fastq.gz") || endsWith(infile, "fq.gz")
  command <- sprintf(
    "bin/fastqindex_0.9.0b extract -s=%i -n=%i -e=%i -f=%s -i=%s | tail -n+8 | md5sum",
    start - 1L,
    n,
    if (is_fastq) 4 else 2,
    infile,
    index
  )
  result <- system(command, intern = TRUE)
  stopifnot(attr(result, "status") == 0)
  c(strtrim(result, 32))
}

# splits a (possibly gzipped) fastx file into n subfiles
# never loads anything into memory, if subfiles are compressed it happens in
# parallel
# this would be more portable with a custom c function (I think?)
fastx_split <- function(infile, n, outroot = tempfile(), compress = FALSE) {
  checkmate::assert_string(infile)
  checkmate::assert_file(infile, access = "r")
  checkmate::assert_int(n, lower = 1, upper = 64)
  checkmate::assert_path_for_output(outroot)
  checkmate::assert_flag(compress)

  is_fastq <- grepl("\\.f(ast)?q)(\\.gz)?$", infile)
  is_gz <- endsWith(infile, ".gz")

  suffix <- if(is_fastq) ".fastq" else ".fasta"
  if (compress) suffix <- paste0(suffix, ".gz")
  command <- if (is_gz) "zcat" else "cat"
  command <- paste(command, infile, "| paste -d'\u1f' - -")
  if (is_fastq) command <- paste(command, "- -")

  digits <- floor(log10(n)) + 1
  command <- paste(
    command,
    " | split",
    sprintf("-nr/%i", n),
    "--numeric-suffixes=1",
    "-a", digits,
    "--additional-suffix", suffix,
    "--filter='tr \"\u1f\" \"\\n\""
    )
  if (compress) command <- paste(command, "| gzip -c -")
  command <- paste(
    command, ">$FILE'",
    "-", outroot
  )
  command <- paste("bash -c", shQuote(command))
  result <- system(command)
  stopifnot(result == 0)
  outfiles <- sprintf(paste0("%s%0", digits, "i%s"), outroot, seq_len(n), suffix)
  checkmate::assert_file(outfiles, "r")
  outfiles
}

# rejoins some split fastx files
# of course this is most useful if something happened to them in between.
# If no lines are missing or reordered in the infiles, then the outfile
# will end up in the same order as the original file which was split.
# Otherwise this will almost certainly not happen.
fastx_combine <- function(infiles, outfile) {
  checkmate::assert_file(infiles, "r")
  checkmate::assert_path_for_output(outfile, overwrite = TRUE)
  is_fastq <- endsWith(infiles, ".fastq") | endsWith(infiles, ".fastq.gz")
  stopifnot(all(is_fastq) | all(!is_fastq))
  is_fastq <- all(is_fastq)

  is_gz <-endsWith(infiles, ".gz")
  stopifnot(all(is_gz) | all(!is_gz))
  is_gz <- all(is_gz)

  compress <- endsWith(outfile, ".gz")

  command <- 'paste -d"\\n"'
  subcommand <- if (is_gz) "<( zcat" else "<( cat"
  subcommand <- paste(subcommand, infiles, "| paste -d'\u1f' - -")
  if (is_fastq) subcommand <- paste(subcommand, "- -")
  subcommand <- paste(subcommand, ")", collapse = " ")
  command <- paste(command, subcommand, "| tr '\u1f' \"\\n\" | sed -n /^$/!p")
  if (compress) command <- paste(command, "| gzip -c -")
  command <- paste(command, ">", outfile)
  command <- paste("bash -c", shQuote(command))
  result <- system(command)
  stopifnot(result == 0)
  checkmate::assert_file(outfile, "r")
  outfile
}

epa_ng <- function(
    ref_msa,
    tree,
    query,
    outdir = tempfile(),
    model,
    ncpu,
    exec = find_executable("epa-ng"),
    redo = TRUE,
    strip_inserts = FALSE
) {
  checkmate::assert_flag(strip_inserts)

  ref_msa_file <- tempfile("reference", fileext = ".fasta")
  if (methods::is(ref_msa, "XStringSet")) {
    assertthat::assert_that(length(unique(Biostrings::width(ref_msa))) == 1)
  } else if (methods::is(ref_msa, "MultipleAlignment")) {
    ref_msa <- methods::as(ref_msa, "XStringSet")
  } else if (is.character(ref_msa)) {
    if (length(ref_msa) == 1 && file.exists(ref_msa)) {
      ref_msa_file <- ref_msa
    } else {
      assertthat::assert_that(length(unique(nchar(ref_msa))) == 1)
      if (has_alphabet(ref_msa, Biostrings::DNA_ALPHABET)) {
        ref_msa <- Biostrings::DNAStringSet(ref_msa)
      } else if (has_alphabet(ref_msa, Biostrings::RNA_ALPHABET)) {
        ref_msa <- Biostrings::RNAStringSet(ref_msa)
      } else if (has_alphabet(ref_msa, Biostrings::AA_ALPHABET)) {
        ref_msa <- Biostrings::AAStringSet(ref_msa)
      } else {
        stop("Unknown alphabet in reference alignment.")
      }
    }
  } else {
    stop("'ref_msa' should be an XStringSet, MultipleAlignment, character vector",
         "of aligned sequences, or filename.")
  }
  if (methods::is(ref_msa, "XStringSet")) {
    Biostrings::writeXStringSet(ref_msa, ref_msa_file)
    on.exit(file.remove(ref_msa_file))
  }

  query_file <- tempfile("query", fileext = ".fasta")
  if (methods::is(query, "XStringSet")) {
    assertthat::assert_that(length(unique(Biostrings::width(query))) == 1)
  } else if (methods::is(query, "MultipleAlignment")) {
    query <- methods::as(query, "XStringSet")
  } else if (is.character(query)) {
    if (length(query) == 1 && file.exists(query)) {
      query_file <- query
    } else {
      assertthat::assert_that(length(unique(nchar(query))) == 1)
      if (has_alphabet(query, Biostrings::DNA_ALPHABET)) {
        query <- Biostrings::DNAStringSet(query)
      } else if (has_alphabet(query, Biostrings::RNA_ALPHABET)) {
        query <- Biostrings::RNAStringSet(query)
      } else if (has_alphabet(query, Biostrings::AA_ALPHABET)) {
        query <- Biostrings::AAStringSet(query)
      } else {
        stop("Unknown alphabet in query alignment.")
      }
    }
  } else {
    stop("'query' should be an XStringSet, MultipleAlignment, character vector",
         "of aligned sequences, or filename.")
  }
  if (methods::is(query, "XStringSet")) {
    Biostrings::writeXStringSet(query, query_file)
  }

  is_gz <- endsWith(query_file, ".gz")
  if (is_gz | strip_inserts) {
    pre_file <- query_file
    query_file <- withr::local_tempfile(fileext=".fasta")
    precommand <- paste(if (is_gz) "zcat" else "cat", pre_file)
    if (strip_inserts) precommand <- paste(precommand, "| tr -d 'acgt'")
    precommand <- paste(precommand, ">", query_file)
    cat("running precommand:", precommand, "\n")
    pre_status <- system(precommand)
    stopifnot(pre_status == 0)
  }

  if (methods::is(tree, "phylo")) {
    tree_file <- tempfile("tree", tempdir(), ".tree")
    ape::write.tree(tree, tree_file)
    on.exit(file.remove(tree_file))
  } else if (is.character(tree) && file.exists(tree)) {
    tree_file <- tree
  } else {
    stop("'tree' should be a phylo object or a file name.")
  }

  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }

  # if (length(model) == 1 && file.exists(model)) {
  #   model_file <- model
  # } else {
  #   model_file <- tempfile("info")
  #   on.exit(unlink(model_file))
  #   writeLines(model, model_file)
  # }

  args <- c(
    "--ref-msa", ref_msa_file,
    "--query", query_file,
    "--tree", tree_file,
    "--outdir", outdir,
    "--model", model
  )

  if (isTRUE(redo)) {
    args <- c(args, "--redo")
  }

  if (!missing(ncpu) && !is.null(ncpu)) {
    assertthat::assert_that(assertthat::is.count(ncpu))
    args <- c(args, "--threads", ncpu)
  }

  processx::run(exec, args = args, echo_cmd = TRUE, stdout = "")
  outfile <- file.path(outdir, "epa_result.jplace")
  checkmate::assert_file_exists(outfile, "r")
  outfile
}

is_jplace <- function(jplace) {
  is.list(jplace) &&
    setequal(
      names(jplace),
      c("tree", "placements", "metadata", "version", "fields")
    )
}

parse_iqtree_model <- function(model) {
  checkmate::assert_character(model)
  if (checkmate::test_file_exists(model, access = "r")) {
    model = readLines(model)
  }
  model <- model[cumsum(grepl("FINALIZING TREE SEARCH", model)) == 1]
  rates <- model[grepl("Rate parameters:", model)]
  outmodel <- ""
  if (length(rates) == 1) {
    rates <- regmatches(rates, regexec("A-C: *([0-9.]+) +A-G: *([0-9.]+) +A-T: *([0-9.]+) +C-G: *([0-9.]+) +C-T: *([0-9.]+) +G-T: *([0-9.]+)", rates))[[1]][-1]
    rates <- as.numeric(rates)
    checkmate::assert_numeric(rates, lower = 0, finite = TRUE, len = 6, any.missing = FALSE)
    outmodel <- paste0(outmodel, "GTR{", paste0(rates, collapse = "/"), "}")
  }
  base_freq <- model[grepl("Base frequencies:", model)]
  if (length(base_freq) == 1) {
    base_freq <- regmatches(
      base_freq,
      regexec(
        pattern = "A: *([0-9.]+) +C: *([0-9.]+) +G: *([0-9.]+) +T: *([0-9.]+)",
        text = base_freq
      )
    )[[1]][-1]
    base_freq <- as.numeric(base_freq)
    checkmate::assert_numeric(base_freq, lower = 0, upper = 1, len = 4, any.missing = FALSE)
    outmodel <- paste0(outmodel, "+FU{", paste0(base_freq, collapse = "/"), "}")
  }
  outmodel
}

gappa_assign <- function(jplace, taxonomy, outgroup, ranks, ncpu = NULL,
                        allow_file_overwriting = FALSE,
                        distribution_ratio = NULL,
                        verbose = FALSE, id_is_int = FALSE) {
  gappa <- find_executable("gappa")
  args <- c("examine", "assign", "--per-query-results")

  if (checkmate::test_file_exists(jplace, access = "r")) {
    jplace_file <- jplace
  } else if (is_jplace(jplace)) {
    jplace_file <- withr::local_tempfile(fileext = ".jplace")
    jsonlite::write_json(jplace, jplace_file, auto_unbox = TRUE)
  } else {
    stop("'jplace' should be a filename or a list containing jplace data")
  }
  args <- c(args, "--jplace-path", jplace_file)

  if (checkmate::test_file_exists(taxonomy, access = "r")) {
    taxonomy_file <- taxonomy
  } else {
    checkmate::assert_data_frame(taxonomy, min.cols = 2, max.cols = 2)
    taxonomy_file <- withr::local_tempfile(fileext = ".txt")
    readr::write_tsv(taxonomy, taxonomy_file, col_names = FALSE)
  }
  args <- c(args, "--taxon-file", taxonomy_file)

  if (checkmate::test_file_exists(outgroup, access = "r")) {
    outgroup_file <- outgroup
  } else {
    checkmate::assert_character(outgroup, any.missing = FALSE)
    outgroup_file <- withr::local_tempfile(fileext = ".txt")
    writeLines(outgroup, outgroup_file)
  }
  args <- c(args, "--root-outgroup", outgroup_file)

  checkmate::assert_character(ranks, any.missing = FALSE)
  ranks <- paste(ranks, collapse = "|")
  args <- c(args, "--ranks-string", ranks)

  checkmate::assert_flag(allow_file_overwriting, null.ok = TRUE)
  if (isTRUE(allow_file_overwriting)) args <- c(args, "--allow-file-overwriting")

  checkmate::assert_flag(verbose, null.ok = TRUE)
  if (isTRUE(verbose)) args <- c(args, "--verbose")

  checkmate::assert_flag(id_is_int)

  checkmate::assert_count(ncpu, null.ok = TRUE)
  if (!is.null(ncpu)) args <- c(args, "--threads", ncpu)

  if (!is.null(distribution_ratio)) {
    checkmate::assert(
      checkmate::check_number(distribution_ratio, lower = 0, upper = 1),
      checkmate::check_number(distribution_ratio, lower = -1, upper = -1)
    )
    args <- c(args, "--distribution-ratio", distribution_ratio)
  }

  out_file <- withr::local_tempfile(fileext = "_per_query.tsv")
  checkmate::assert_path_for_output(out_file)
  args <- c(
    args,
    "--out-dir", dirname(out_file),
    "--file-prefix", sub("per_query.tsv$", "", basename(out_file))
  )

  processx::run(gappa, args, error_on_status = TRUE, echo_cmd = TRUE, stdout = "")
  parse_gappa_per_query(out_file, ranks, id_is_int)
}

parse_gappa_per_query <- function(per_query, ranks, id_is_int = FALSE) {
  checkmate::assert_file(per_query, "r")
  checkmate::assert_character(ranks, any.missing = FALSE)
  checkmate::assert_flag(id_is_int)
  pq <- readr::read_tsv(per_query, col_types = if (id_is_int) "innnnc" else "cnnnnc") |>
    dplyr::mutate(
      rank = factor(stringr::str_count(taxopath, ";") + 1, levels = seq_along(ranks), labels = ranks, ordered = TRUE),
      taxopath = chartr(";", ",", taxopath)
    )
  pq <- dplyr::bind_rows(
    dplyr::filter(pq, fract > 0, rank != dplyr::last(ranks)) |>
    dplyr::transmute(
      name,
      rank = factor(ranks[as.integer(rank) + 1L], levels = ranks, ordered = TRUE),
      parent_taxon = taxopath,
      taxon = NA_character_,
      prob = fract
    ) |>
      dplyr::filter(prob > 0),
    tidyr::extract(pq, taxopath, c("parent_taxon", "taxon"), , regex = "(?:(.+),)?([^,]+)$") |>
      dplyr::transmute(
        name,
        rank,
        parent_taxon = dplyr::na_if(parent_taxon, ""),
        taxon,
        prob = afract
      )
  ) |>
    dplyr::arrange(name, rank)

  names(pq)[1] <- if (id_is_int) "seq_idx" else "seq_id"
  pq
}
