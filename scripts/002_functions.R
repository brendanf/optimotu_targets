bimera_denovo_table <- function(
  seqtab,
  minFoldParentOverAbundance = 1.5,
  minParentAbundance = 2,
  allowOneOff = FALSE,
  minOneOffParentDistance = 4,
  maxShift = 16,
  multithread = FALSE
) {
  if (isTRUE(multithread)) {
    RcppParallel::setThreadOptions(numThreads = "auto")
  } else if (isFALSE(multithread)) {
    RcppParallel::setThreadOptions(numThreads = 1)
  } else {
    assertthat::assert_that(
      assertthat::is.count(multithread),
      msg = "argument 'multithread' must be TRUE, FALSE, or a positive integer."
    )
    RcppParallel::setThreadOptions(numThreads = multithread)
  }
  dada2:::C_table_bimera2(
    mat = seqtab,
    seqs = colnames(seqtab),
    min_fold = minFoldParentOverAbundance,
    min_abund = minParentAbundance,
    allow_one_off = allowOneOff,
    min_one_off_par_dist = minOneOffParentDistance,
    match = dada2::getDadaOpt("MATCH"),
    mismatch = dada2::getDadaOpt("MISMATCH"),
    gap_p = dada2::getDadaOpt("GAP_PENALTY"),
    max_shift = maxShift
  ) %>%
    tibble::as_tibble() %>%
    tibble::add_column(seq = colnames(seqtab))
}

combine_bimera_denovo_tables <- function(
  bimdf,
  minSampleFraction = 0.9,
  ignoreNNegatives = 1L,
  verbose = FALSE
) {
  bimdf <- dplyr::group_by(bimdf, seq) %>%
    dplyr::summarize(dplyr::across(everything(), sum), .groups = "drop")
  ## This snippet modified from DADA2
  is.bim <- function(nflag, nsam, minFrac, ignoreN) {
    nflag >= nsam || (nflag > 0 && nflag >= (nsam - ignoreN) * 
                        minFrac)
  }
  bims.out <- mapply(is.bim, bimdf$nflag, bimdf$nsam, minFrac = minSampleFraction, 
                     ignoreN = ignoreNNegatives)
  names(bims.out) <- bimdf$seq
  if (verbose) 
    message("Identified ", sum(bims.out), " bimeras out of ", 
            length(bims.out), " input sequences.")
  ## end snippet from DADA2
  return(bims.out)
}

remove_bimera_denovo_tables <- function(
  seqtabs,
  bimdf,
  minSampleFraction = 0.9,
  ignoreNNegatives = 1L,
  verbose = FALSE
) {
  bims.out <- combine_bimera_denovo_tables(
    bimdf,
    minSampleFraction = minSampleFraction,
    ignoreNNegatives = ignoreNNegatives,
    verbose = verbose
  )
  remove_chimeras <- function(seqtab, ischim) {
    seqtab[,!ischim[colnames(seqtab)], drop = FALSE]
  }
  seqtabs <- lapply(seqtabs, remove_chimeras, ischim = bims.out)
  if (length(seqtabs) > 1) {
    dada2::mergeSequenceTables(tables = seqtabs)
  } else {
    seqtabs[[1]]
  }
}

sort_seq_table <- function(seqtable) {
  colorder <- order(
    -colSums(seqtable > 0), # prevalence, highest to lowest
    -colSums(seqtable), # abundance, highest to lowest
    -apply(seqtable, 2, var), # variance, highest to lowest
    colnames(seqtable) # sequence, alphabetical
  )
  if (is.null(attr(seqtable, "map"))) {
    seqtable[order(rownames(seqtable)), colorder]
  } else {
    structure(
      seqtable[order(rownames(seqtable)), colorder],
      map = dplyr::mutate(attr(seqtable, "map"), seq_id_out = order(colorder)[seq_id_out])
    )
  }
}

summarize_by_rank <- function(rank, superrank, data) {
  rank_sym <- as.symbol(rank)
  superrank_sym <- as.symbol(superrank)
  dplyr::filter(
    data,
    !startsWith(!!superrank_sym, "dummy_"),
    !startsWith(!!rank_sym, "dummy_"),
    !is.na(!!rank_sym)
  ) %>%
    dplyr::group_by(!!superrank_sym) %>%
    dplyr::summarize(
      superrank = superrank,
      rank = rank,
      n_taxa = dplyr::n_distinct(!!rank_sym),
      n_seq = dplyr::n_distinct(seq_id),
      seq_id = list(seq_id),
      true_taxa = list(as.integer(factor(!!rank_sym)))
    ) %>%
    dplyr::rename(supertaxon = !!superrank)
}

#' Calculate clustering thresholds for each taxon, falling back to its ancestor
#' taxa as necessary
#'
#' @param rank (`character` string) the rank within which clustering will be
#' performed
#' @param conf_level (`character` string) as in `fmeasure_optima`
#' @param taxon_table (`data.frame`) taxonomy table; column "seq_id" gives the
#' sequence ID, and columns "kingdom" to "species" give the taxonomy at each
#' rank.
#' @param fmeasure_optima (`data.frame`) optimum clustering thresholds within
#' various taxa; column "rank" gives the rank which is approximated by
#' clustering; "superrank" gives the rank of the taxon within which the
#' clustering threshold was optimized; "supertaxon" gives that taxon name;
#' "conf_level" gives a string description of the confidence level threshold for
#' taxonomic assignments; "threshold" gives the optimum clustering threshold;
#' "f_measure" gives the F measure at the optimum threshold.
#' @param default (`character string`) default taxon to define threshold to use
#' when taxonomy is unknown. default: "Fungi"
#'
#'
#' @return (named `numeric`, where names are taxa and values are clustering
#' thresholds)
#' unknown sequences, grouped by taxonomy at the parent rank. Sequences where
#' the parent rank is also unknown are in the item named `"_NA_"`.
calc_taxon_thresholds <- function(rank, conf_level, taxon_table,
                                  fmeasure_optima, default = "Fungi") {
  rank_name <- rlang::sym(rank)
  dplyr::select(taxon_table, kingdom:!!rank_name) %>%
    dplyr::filter(!is.na(!!rank_name)) %>%
    unique() %>%
    purrr::reduce(
      c(superranks(rank), rank),
      function(thresholds, r) {
        dplyr::left_join(
          thresholds,
          dplyr::filter(
            fmeasure_optima,
            rank == subranks(!!rank)[1],
            superrank == r,
            conf_level == !!conf_level
          ) %>%
            dplyr::select(
              !!r := supertaxon,
              !!paste0("threshold_", r) := threshold
            ),
          by = r
        )
      },
      .init = .
    ) %>%
    dplyr::transmute(
      !!rank_name := !!rank_name,
      threshold = {.} %>%
        dplyr::select(dplyr::starts_with("threshold")) %>%
        rev() %>%
        do.call(dplyr::coalesce, .)
    ) %>%
    tibble::deframe() %>%
    c("_NA_" = dplyr::filter(
      fmeasure_optima,
      rank == subranks(!!rank)[1],
      supertaxon == default,
      conf_level == !!conf_level
    )$threshold)
}

threshold_as_dist <- function(thresholds) {
  if (any(thresholds > 50)) {
    1 - 0.01 * thresholds
  } else if (any(thresholds > 1)) {
    0.01 * thresholds
  } else if (any(thresholds > 0.5) ) {
    1 - thresholds
  } else {
    thresholds
  }
}

#' Calculate clustering thresholds for each taxon, falling back to its ancestor
#' taxa as necessary
#'
#' @param rank (`character` string) the rank within which clustering will be
#' performed
#' @param conf_level (`character` string) as in `fmeasure_optima`
#' @param taxon_table (`data.frame`) taxonomy table; column "seq_id" gives the
#' sequence ID, and columns "kingdom" to "species" give the taxonomy at each
#' rank.
#' @param fmeasure_optima (`data.frame`) optimum clustering thresholds within
#' various taxa; column "rank" gives the rank which is approximated by
#' clustering; "superrank" gives the rank of the taxon within which the
#' clustering threshold was optimized; "supertaxon" gives that taxon name;
#' "conf_level" gives a string description of the confidence level threshold for
#' taxonomic assignments; "threshold" gives the optimum clustering threshold;
#' "f_measure" gives the F measure at the optimum threshold.
#' @param default (`character string`) default taxon to define threshold to use
#' when taxonomy is unknown. default: "Fungi"
#'
#'
#' @return (named `list` of `double` vectors)
calc_subtaxon_thresholds <- function(rank, conf_level, taxon_table,
                                  fmeasure_optima, default = "Fungi") {
  rank_name <- rlang::sym(rank)
  dplyr::select(taxon_table, kingdom:!!rank_name) %>%
    tidyr::crossing(subrank = subranks(rank)) %>%
    dplyr::filter(!is.na(!!rank_name)) %>%
    unique() %>%
    purrr::reduce(
      c(superranks(rank), rank),
      function(thresholds, r) {
        dplyr::left_join(
          thresholds,
          dplyr::filter(
            fmeasure_optima,
            rank %in% subranks(!!rank),
            superrank == r,
            conf_level == !!conf_level
          ) %>%
            dplyr::select(
              subrank = rank,
              !!r := supertaxon,
              !!paste0("threshold_", r) := threshold
            ),
          by = c("subrank", r)
        )
      },
      .init = .
    ) %>%
    dplyr::transmute(
      subrank = rank2factor(subrank),
      !!rank_name := !!rank_name,
      threshold = {.} %>%
        dplyr::select(dplyr::starts_with("threshold")) %>%
        rev() %>%
        do.call(dplyr::coalesce, .) %>%
        threshold_as_dist()
    ) %>%
    dplyr::arrange(subrank) %>%
    split(.[[rank]]) %>%
    lapply(dplyr::select, !any_of(rank)) %>%
    lapply(tibble::deframe) %>%
    lapply(cummax) %>%
    c(
      "_NA_" = dplyr::filter(
        fmeasure_optima,
        rank %in% subranks(!!rank),
        supertaxon == default,
        conf_level == !!conf_level
      ) %>% dplyr::transmute(rank = rank2factor(rank), threshold_as_dist(threshold)) %>%
        dplyr::arrange(rank) %>%
        tibble::deframe() %>%
        cummax() %>%
        list()
    )
}

parse_protax_nameprob <- function(nameprob) {
    set_names(nameprob, basename(nameprob)) |>
    lapply(readLines) |>
    tibble::enframe() |>
    tidyr::extract(
      name,
      into = "rank",
      regex = "query(\\d+)\\.nameprob",
      convert = TRUE
    ) |>
    tidyr::unchop(value) |>
    dplyr::mutate(
      rank = rank2factor(TAXRANKS[rank]),
      value = gsub("([^\t]+)\t([0-9.]+)", "\\1:\\2", value) %>%
        gsub("(:[0-9.]+)\t", "\\1;", .)
    ) |>
    tidyr::separate(value, into = c("seq_id", "nameprob"), sep = "\t", fill = "right") |>
    tidyr::separate_rows(nameprob, sep = ";") |>
    tidyr::separate(nameprob, into = c("name", "prob"), sep = ":", convert = TRUE) |>
    tidyr::extract(name, into = c("parent_taxonomy", "taxon"), regex = "(.+),([^,]+)$") |>
    dplyr::mutate(
      taxon = dplyr::na_if(taxon, "unk"),
      prob = ifelse(is.na(taxon), 0, prob)
    ) |>
    dplyr::arrange(seq_id, rank, dplyr::desc(prob))
}

read_sfile <- function(file) {
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
          col_positions =  stringr::str_locate_all(x$text[1], "-+")[[1]] |>
            tibble::as_tibble() |>
            tibble::add_column(
              col_names = c("idx", "seq_id", "match_len", "cm_from", "cm_to",
                            "trunc", "bit_sc", "avg_pp", "time_band_calc",
                            "time_alignment", "time_total", "mem_mb")
            ) |>
            do.call(readr::fwf_positions, args = _),
          skip = 1,
          col_types = "iciiicdddddd"
        )
      }
    )

}
