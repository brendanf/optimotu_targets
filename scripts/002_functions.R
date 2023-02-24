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
  dada2::mergeSequenceTables(tables = seqtabs)
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

# combine tip classifications to build a full PROTAX taxonomy

build_taxonomy <- function(...) {
  tax <- tibble::tibble(
    classification = union(...),
    rank = ifelse(
      classification == "root",
      0L,
      stringr::str_count(classification, stringr::fixed(",")) + 1L),
    parent = ifelse(
      rank <= 1,
      "root",
      sub(",[^,]+$", "", classification)
    )
  ) %>%
    split(.$rank)
  tax[[1]]$taxon_id = 0L
  tax[[1]]$prior = 1
  tax[[1]]$parent_id = 0L
  tax[[2]]$taxon_id = 1L:2L
  tax[[2]]$prior = c(0.99, 0.01)

  for (r in length(tax):3L) {
    if (r == length(tax)) {
      tax[[r]]$prior = 0.99/nrow(tax[[r]])
    } else {
      tax[[r]]$prior <- NULL
      tax[[r]] <- dplyr::left_join(
        tax[[r]],
        dplyr::group_by(tax[[r+1]], parent) %>%
          dplyr::summarise(prior = sum(prior)) %>%
          dplyr::rename(classification = parent),
        by = "classification"
      )
    }
  }
  for (r in 2L:length(tax)) {
    tax[[r]]$taxon_id <- seq_len(nrow(tax[[r]])) + max(tax[[r-1L]]$taxon_id)
    tax[[r]]$parent_id <- NULL
    tax[[r]] <- dplyr::left_join(
      tax[[r]],
      dplyr::select(
        tax[[r-1]],
        parent = classification,
        parent_id = taxon_id
      ),
      by = "parent"
    )
  }
  dplyr::bind_rows(tax) %>%
    dplyr::select(taxon_id, parent_id, rank, classification, prior)
}

# Format a classification for use as a Sintax reference DB
sintax_format <- function(s) {
  s <- sub(",", ";p:", s, fixed = TRUE)
  s <- sub(",", ";c:", s, fixed = TRUE)
  s <- sub(",", ";o:", s, fixed = TRUE)
  s <- sub(",", ";f:", s, fixed = TRUE)
  s <- sub(",", ";g:", s, fixed = TRUE)
  s <- sub(",", ";s:", s, fixed = TRUE)
  s <- chartr(";", ",", s)
  paste0("tax=d:", s, ";")
}

# Truncate classification(s) at a given rank
truncate_taxonomy <- function(s, rank) {
  regex <- paste0("(^([^,]+,){", rank-1L, "}[^,]+).*")
  out <- gsub(regex, "\\1", s)
  out[!grepl(regex, s)] <- NA_character_
  out
}

# Find OTUs which contain sequences with _any_ probability
# of being a target taxon
find_target_taxa <- function(target_taxa, asv_all_tax_prob, asv_taxonomy, otu_taxonomy) {
  asv_otu_key <-
    dplyr::inner_join(
      dplyr::select(asv_taxonomy, asv_seq_id = seq_id, species),
      dplyr::select(otu_taxonomy, seq_id, species),
      by = "species"
    ) %>%
    dplyr::select(-species)
  otu_long_taxonomy <- tidyr::pivot_longer(
    otu_taxonomy,
    kingdom:species,
    names_to = "rank",
    values_to = "otu_taxon",
    names_transform = list(rank = rank2factor)
  ) %>%
    dplyr::select(seq_id, rank, otu_taxon)
  dplyr::select(asv_all_tax_prob, asv_seq_id = seq_id, rank, protax_taxon = taxon, protax_prob = prob) %>%
    dplyr::inner_join(asv_otu_key, by = "asv_seq_id") %>%
    dplyr::group_by(seq_id) %>%
    dplyr::filter(any(protax_taxon %in% target_taxa)) %>%
    dplyr::left_join(otu_long_taxonomy, by = c("seq_id", "rank")) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(seq_id, asv_seq_id, rank) %>%
    dplyr::select(seq_id, asv_seq_id, otu_taxon, rank, everything())
}
