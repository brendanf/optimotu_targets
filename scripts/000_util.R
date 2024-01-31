#### Utility functions ####

# are we running slurm?
is_slurm <- function() nchar(Sys.getenv("SLURM_JOB_ID")) > 0 || nchar(Sys.which("sbatch")) > 0
is_local <- function() !is_slurm()

# are we running snakemake?
is_snakemake <- function() !interactive() && exists("snakemake")

# how many cpus do we have on the local machine?
# if we're not running on the cluster, leave one cpu free.
local_cpus <- function() {
  if (is_snakemake()) {
    snakemake@threads
  } else if (is_slurm()) {
    out <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
    # if slurm doesn't know how many cores, that means we're probably on the
    # login node, so we should only use 1 core.
    if (!assertthat::is.count(out)) out <- 1L
    out
  } else {
    max(parallel::detectCores() - 1L, 1L)
  }
}

#### convenience function for writing a file and returning its name ####

ensure_directory <- function(file) {
  d <- dirname(file)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

write_and_return_file <- function(x, file, ...) {
  UseMethod("write_and_return_file")
}

write_and_return_file.XStringSet <- function(x, file, ...) {
  ensure_directory(file)
  Biostrings::writeXStringSet(x, file, ...)
  file
}

write_and_return_file.data.frame <- function(x, file, type = c("rds", "tsv"), ...) {
  ensure_directory(file)
  type = match.arg(type)
  switch(
    type,
    rds = saveRDS(x, file, ...),
    tsv = readr::write_tsv(x, file, ...),
    stop("Unknown file type: ", type)
  )
  file
}

write_and_return_file.character <- function(x, file, ...) {
  ensure_directory(file)
  writeLines(x, file, ...)
  file
}

write_and_return_file.ggplot <- function(x, file, ...) {
  ensure_directory(file)
  ggplot2::ggsave(file, plot = x, ...)
  file
}

write_and_return_file.default <- function(x, file, ...) {
  ensure_directory(file)
  saveRDS(x, file, ...)
  file
}

# helper functions to work with sequences sets that may be XStringSet,
# named character, or data.frame

# guess the column name in a data frame which refers to the sequence ID
find_name_col <- function(d) {
  stopifnot(is.data.frame(d))
  if ("seq_id" %in% names(d)) return("seq_id")
  if ("name" %in% names(d)) return("name")
  if ("ASV" %in% names(d)) return("ASV")
  if ("OTU" %in% names(d)) return("OTU")
  if (ncol(d) == 2) {
    if ("seq" %in% names(d)) return(names(d)[names(d) != "seq"])
    if ("sequence" %in% names(d)) return(names(d)[names(d) != "sequence"])
    return(names(d)[1])
  }
  stop("unable to determine sequence name column:", names(d))
}

# guess the column name in a data frame which refers to the sequence
find_seq_col <- function(d) {
  stopifnot(is.data.frame(d))
  if ("seq" %in% names(d)) return("seq")
  if ("sequence" %in% names(d)) return("sequence")
  if (ncol(d) == 2) {
    if ("seq_id" %in% names(d)) return(names(d)[names(d) != "seq_id"])
    if ("name" %in% names(d)) return(names(d)[names(d) != "name"])
    if ("ASV" %in% names(d)) return(names(d)[names(d) != "ASV"])
    if ("OTU" %in% names(d)) return(names(d)[names(d) != "OTU"])
    return(names(d)[2])
  }
  stop("unable to determine sequence column:", names(d))
}

# write the sequence as a fasta file
write_sequence <- function(seq, fname, ...) {
  UseMethod("write_sequence", seq)
}

write_sequence.data.frame <- function(seq, fname, seq_col = find_seq_col(seq),
                                      name_col = find_name_col(seq), ...) {
  dplyr::select(seq, !!name_col, !!seq_col) %>%
    tibble::deframe() %>%
    Biostrings::DNAStringSet() %>%
    write_and_return_file(fname, ...)
}

write_sequence.character <- function(seq, fname, ...) {
  Biostrings::DNAStringSet(seq) %>%
    write_and_return_file(fname, ...)
}

write_sequence.XStringSet <- function(seq, fname, ...) {
  write_and_return_file(seq, fname, ...)
}

# select elements from the sequence set.
select_sequence <- function(seq, which, negate = FALSE, ...) {
  UseMethod("select_sequence", seq)
}

select_sequence.data.frame <- function(seq, which, name_col = find_name_col(seq), negate = FALSE, ...) {
  if (isTRUE(negate)) {
    if (is.integer(which)  ||
        (is.numeric(which) && all(which == round(which)))) return(seq[-which,])
    if (is.logical(which)) return(seq[!which,])
    if (is.character(which)) return(dplyr::filter(seq, !.data[[name_col]] %in% which))
    stop("'which' should be integer, logical, or character")
  } else if (isFALSE(negate)) {
    if (is.integer(which)  ||
        (is.numeric(which) && all(which == round(which)))) return(seq[which,])
    if (is.logical(which)) return(seq[which,])
    if (is.character(which)) return(dplyr::filter(seq, .data[[name_col]] %in% which))
    stop("'which' should be integer, logical, or character")
  }
  stop("'negate' must be TRUE or FALSE")
}

select_sequence.default <- function(seq, which, negate = FALSE, ...) {
  if (isTRUE(negate)) {
    if (is.integer(which) || (is.numeric(which) && all(which == round(which)))) return(seq[-which])
    if (is.logical(which)) return(seq[!which])
    if (is.character(which)) return(seq[setdiff(names(seq), which)])
    stop("'which' should be integer, logical, or character")
  } else if (isFALSE(negate)) {
    return(seq[which])
  }
  stop("'negate' must be true or false")
}

# get the number of sequences in the sequence set
sequence_size <- function(seq, ...) {
  UseMethod("sequence_size", seq)
}

sequence_size.XStringSet <- function(seq, ...) {
  length(seq)
}

sequence_size.character <- function(seq, ...) {
  if (all(file.exists(seq))) {
    if (all(grepl("fq|fastq", seq))) {
      return(
        lapply(seq, Biostrings::fastq.seqlengths) %>%
        vapply(length, 1L)
      )
    } else if (all(grepl("fas?|fasta", seq))) {
      return(
        lapply(seq, Biostrings::fasta.seqlengths) %>%
        vapply(length, 1L)
      )
    }
  }
  sequence_size.default(seq, ...)
}

sequence_size.default <- function(seq, ...) {
  vctrs::vec_size(seq)
}

# convert a character to an ordered factor of taxonomic ranks
TAXRANKS <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
rank2factor <- function(x) {
  factor(x, levels = rev(TAXRANKS), ordered = TRUE)
}

int2rank <- function(x) {
  rank2factor(TAXRANKS[x])
}

superranks <- function(x, ranks = TAXRANKS) {
  ranks[rank2factor(ranks) > x]
}

subranks <- function(x, ranks = TAXRANKS) {
  ranks[rank2factor(ranks) < x]
}

# force a string to be ASCII

ascii_clean <- function(s) {
  gsub(
    s,
    pattern = "[^A-Za-z0-9[:punct:]]",
    replacement = "x",
    useBytes = TRUE,
    perl = TRUE
  )
}

# Get all the target names defined in a plan

get_target_names <- function(plan) {
  if (methods::is(plan, "tar_target")) {
    plan$settings$name
  } else {
    unname(unlist(lapply(plan, get_target_names)))
  }
}

# generate hash codes from sequences
seqhash <- digest::getVDigest("spookyhash")

# generate names like "ASV0001", "ASV0002", ...

make_seq_names <- function(n, prefix) {
  sprintf(
    sprintf("%s%%0%dd", prefix, floor(log10(n)) + 1),
    seq.int(n)
  )
}

name_seqs <- function(seq, prefix, ...) {
  UseMethod("name_seqs", seq)
}

name_seqs.XStringSet <- function(seq, prefix, ...) {
  names(seq) <- make_seq_names(length(seq), prefix)
  seq
}

name_seqs.character <- function(seq, prefix, ...) {
  names(seq) <- make_seq_names(length(seq), prefix)
  seq
}

name_seqs.data.frame <- function(seq, prefix, id_col = prefix, ...) {
  seq[[id_col]] <- make_seq_names(nrow(seq), prefix)
  seq
}

name_seqs.matrix <- function(seq, prefix, ...) {
  colnames(seq) <- make_seq_names(ncol(seq), prefix)
  seq
}

#' Convert an object into a long (i.e. sparse) sequence occurrence table
#'
#' @param x (`data.frame` as returned by `dada2::mergePairs`, or integer matrix as returned by `dada2::makeSequenceTable`, or a list of one of these.)
#' @param rc (logical flag) if TRUE, sequences in `x` will be reverse complemented.
#'
#' @return a `data.frame` with columns `sample`, `seq`, and `nread`

make_long_sequence_table <- function(x, rc = FALSE) {
  UseMethod("make_long_sequence_table", x)
}

make_long_sequence_table.data.frame <- function(x, rc = FALSE) {
  checkmate::assert_data_frame(x, col.names = "named")
  checkmate::assert_names(names(x), must.include = c("sequence", "abundance"))
  checkmate::assert_flag(rc)
  if ("accept" %in% names(x)) {
    checkmate::assert_logical(x$accept)
    x <- x[x$accept,]
  }
  out <- x[c("sequence", "abundance")]
  names(out) <- c("seq", "nread")
  if (isTRUE(rc)) out$seq <- dada2::rc(out$seq)
  out
}

make_long_sequence_table.matrix <- function(x, rc = FALSE) {
  checkmate::assert_integerish(x)
  checkmate::assert_flag(rc)
  if (isTRUE(rc)) colnames(x) <- dada2::rc(colnames(x))
  if (typeof(x) != "integer") mode(x) <- "integer"
  x[x==0L] <- NA_integer_
  as.data.frame(x) |>
    tibble::rownames_to_column("sample") |>
    tidyr::pivot_longer(-1, names_to = "seq", values_to = "nread", values_drop_na = TRUE)
}

make_long_sequence_table.list <- function(x, rc = FALSE) {
  out <- if (checkmate::test_list(x, types = "data.frame")) {
    checkmate::assert_named(x)
    purrr::map_dfr(x, make_long_sequence_table.data.frame, rc = rc, .id = "sample")
  } else if (checkmate::test_list(x, types = "matrix")) {
    purrr::map_dfr(x, make_long_sequence_table.matrix, rc = rc)
  } else {
    stop("cannot determine entry type in make_long_sequence_table.list")
  }
  dplyr::summarize(out, nread = sum(nread), .by = c(sample, seq))
}

#' Convert an object into a long (i.e. sparse) sequence occurrence table where
#' sequences are stored as integer indices to a master list
#'
#' @param x (`data.frame` as returned by `dada2::mergePairs`, or integer matrix as returned by `dada2::makeSequenceTable`, or a list of one of these.)
#' @param seqs (`character` vector) master list of sequences
#' @param rc (logical flag) if TRUE, sequences in `x` will be reverse complemented.
#'
#' @return a `data.frame` with columns `sample`, `seq`, and `nread`

make_mapped_sequence_table <- function(x, seqs, rc = FALSE) {
  UseMethod("make_mapped_sequence_table", x)
}

make_mapped_sequence_table.data.frame <- function(x, seqs, rc = FALSE) {
  checkmate::assert_data_frame(x, col.names = "named")
  checkmate::assert_names(names(x), must.include = c("sequence", "abundance"))
  checkmate::assert_flag(rc)
  if ("accept" %in% names(x)) {
    checkmate::assert_logical(x$accept)
    x <- x[x$accept,]
  }
  if (checkmate::test_file_exists(seqs, "r")) {
    seqs <- Biostrings::readDNAStringSet(seqs)
  }
  out <- x[c("sequence", "abundance")]
  names(out) <- c("seq_idx", "nread")
  if (isTRUE(rc)) {
    out$seq_idx <- BiocGenerics::match(dada2::rc(out$seq_idx), seqs)
  } else {
    out$seq_idx <- BiocGenerics::match(out$seq_idx, seqs)
  }
  out
}

make_mapped_sequence_table.matrix <- function(x, seqs, rc = FALSE) {
  checkmate::assert_integerish(x)
  checkmate::assert_flag(rc)
  if (isTRUE(rc)) colnames(x) <- dada2::rc(colnames(x))
  if (checkmate::test_file_exists(seqs, "r")) {
    seqs <- Biostrings::readDNAStringSet(seqs)
  }
  colnames(x) <- BiocGenerics::match(colnames(x), seqs)
  if (typeof(x) != "integer") mode(x) <- "integer"
  x[x==0L] <- NA_integer_
  as.data.frame(x) |>
    tibble::rownames_to_column("sample") |>
    tidyr::pivot_longer(
      -1,
      names_to = "seq_idx",
      names_transform = as.integer,
      values_to = "nread",
      values_drop_na = TRUE
    )
}

make_mapped_sequence_table.list <- function(x, seqs, rc = FALSE) {
  if (checkmate::test_file_exists(seqs, "r")) {
    seqs <- Biostrings::readDNAStringSet(seqs)
  }
  out <- if (checkmate::test_list(x, types = "data.frame")) {
    checkmate::assert_named(x)
    purrr::map_dfr(x, make_mapped_sequence_table.data.frame, seqs = seqs, rc = rc, .id = "sample")
  } else if (checkmate::test_list(x, types = "matrix")) {
    purrr::map_dfr(x, make_mapped_sequence_table.matrix, seqs = seqs, rc = rc)
  } else {
    stop("cannot determine entry type in make_long_sequence_table.list")
  }
  dplyr::summarize(out, nread = sum(nread), .by = c(sample, seq_idx))
}

unnest_yaml_list <- function(x) {
  checkmate::assert_list(x)
  if (
    is.null(names(x)) &&
    checkmate::check_list(x, types = "list") &&
    all(vapply(x, length, 1L) == 1)
  ) {
    do.call(c, x)
  } else {
    x
  }
}
