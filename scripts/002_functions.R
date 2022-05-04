# faster than length(intersect(x, y)), but assumes x = unique(x) and y = unique(y)
# this could be even faster since we know they are sorted, but there is no fast R
# function I am aware of.
intersect_length <- function(x, y) sum(!is.na(match(x, y)))

# compiled c++ version gives ~4x speedup
Rcpp::sourceCpp("src/intersect.cpp")

inner_fmeasure <- function(cj, kpartition, nk) {
  nc <- length(cj)
  nc * max(purrr::map_int(kpartition, intersect_length_, cj)/(nc + nk))
}

# Calculate F-measure for delimitation by clustering
f_measure <- function(data, c, k) {
  nseq <- nrow(data)
  cpartition <- split(seq_along(data[[c]]), data[[c]])
  kpartition <- split(seq_along(data[[k]]), data[[k]])
  nk <- purrr::map_int(kpartition, length)
  2 / nseq * sum(
    purrr::map_dbl(
      cpartition,
      inner_fmeasure,
      kpartition,
      nk
    )
  )
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
