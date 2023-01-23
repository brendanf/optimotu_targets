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
  min_sample_fraction = 0.9,
  ignoreNNegatives = 1L,
  verbose = FALSE
) {
  bimdf <- dplyr::group_by(bimdf, seq) %>%
    dplyr::summarize(dplyr::across(everything(), sum), .groups = "drop")
  is.bim <- function(nflag, nsam, minFrac, ignoreN) {
    nflag >= nsam || (nflag > 0 && nflag >= (nsam - ignoreN) * 
                        minFrac)
  }
  bims.out <- mapply(is.bim, bimdf$nflag, bimdf$nsam, minFrac = minSampleFraction, 
                     ignoreN = ignoreNNegatives)
  names(bims.out) <- sqs
  if (verbose) 
    message("Identified ", sum(bims.out), " bimeras out of ", 
            length(bims.out), " input sequences.")
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
#' @param taxon_table (`data.frame`) taxonomy table; column "ASV" gives the
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
