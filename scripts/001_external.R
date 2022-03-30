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
  list.files(outdir)
}
