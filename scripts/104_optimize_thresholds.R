library(tarchetypes)

threshold_plan <- list(
  tar_fst_tbl(
    threshold_meta,
    tibble::tibble(
      rank_int = seq_along(TAXRANKS)[-1],
      rank = TAXRANKS[rank_int]
    ),
    deployment = "main"
  ),
  tar_file_fast(
    protax_reftax_file,
    file.path(protax_modeldir, sprintf("ref.tax%d", threshold_meta$rank_int)),
    pattern = map(threshold_meta),
    deployment = "main"
  ),
  
  tar_fst_tbl(
    reftax,
    dplyr::filter(asv_all_tax_prob, prob > 0.5) |>
      tidyr::pivot_wider(
        values_from = "taxon",
        names_from = "rank",
        id_cols = "seq_id"
      ) |>
      dplyr::select(seq_id, dplyr::any_of(TAXRANKS)),
    deployment = "main"
  ),
  tar_target(
    testset_select,
    purrr::map_dfr(
      superranks(threshold_meta$rank),
      summarize_by_rank,
      data = reftax,
      rank = threshold_meta$rank
    ) %>%
      dplyr::filter(n_taxa >= 5 | superrank == "kingdom", n_seq >= 10),
    tidy_eval = FALSE,
    pattern = map(threshold_meta),
    deployment = "main"
  ),
  tar_target(
    testset_rowwise,
    testset_select,
    deployment = "main"
  ),
  tar_target(
    threshold_testset,
    optimotu::seq_cluster_usearch(
      seq = asv_seq,
      threshold_config = optimotu::threshold_uniform(
        from = 0,
        to = 0.4,
        by = 0.001,
        thresh_names = as.character(1000 - 0:400)
      ),
      clust_config = optimotu::clust_tree(),
      parallel_config = optimotu::parallel_concurrent(local_cpus()%/%2.5),
      which = testset_select$seq_id,
      usearch = "bin/usearch",
      usearch_ncpu = local_cpus()
    ),
    iteration = "list",
    deployment = "worker"
  ),
  tar_target(
    threshold_ntaxa,
    lapply(threshold_testset, apply, 1, dplyr::n_distinct) |>
      lapply(tibble::enframe, name = "threshold", value = "ntaxa") |>
      tibble::add_column(testset_rowwise, ntaxa = _) |>
      dplyr::select(supertaxon, superrank, rank, ntaxa) |>
      tidyr::unnest(ntaxa),
    deployment = "worker"
  ),
  tar_fst_tbl(
    cluster_metrics,
    purrr::map_dfr(
      seq_along(threshold_testset),
        function(i) {
          purrr::map_dfc(
            .x = list(
              optimotu::confusion_matrix,
              optimotu::adjusted_mutual_information,
              mFM = optimotu::fmeasure
            ),
            .f = purrr::exec,
            k = threshold_testset[[i]],
            c = testset_rowwise$true_taxa[[i]],
            local_cpus()
          ) |>
            tibble::remove_rownames() %>%
            dplyr::mutate(
              MCC = optimotu::matthews_correlation_coefficient(.),
              RI = optimotu::rand_index(.),
              ARI = optimotu::adjusted_rand_index(.),
              FMI = optimotu::fowlkes_mallow_index(.),
              threshold = (1000 - 0:400)/10,
              rank = testset_rowwise$rank[i],
              superrank = testset_rowwise$superrank[i],
              supertaxon = testset_rowwise$supertaxon[i]
            )
        }
      ),
      deployment = "worker"
    ),
  tar_fst_tbl(
    optima,
    cluster_metrics %>%
      tidyr::pivot_longer(
        -c(threshold, rank, superrank, supertaxon),
        names_to = "metric", values_to = "score"
      ) |>
      dplyr::filter(!(metric %in% c("TP", "FP", "FN", "TN"))) |>
      dplyr::group_by(rank, superrank, supertaxon, metric) |>
      dplyr::arrange(dplyr::desc(score)) |>
      dplyr::summarize(
        threshold = threshold[which.max(score)],
        score = max(score),
        .groups = "drop"
      ),
    deployment = "main"
  ),
  tar_file_fast(
      optima_file,
      write_and_return_file(
        optima, file.path("data", sprintf("GSSP_thresholds.tsv", refset_name)),
        "tsv"),
      deployment = "main"
  ),
  tar_fst_tbl(
    fmeasure_optima,
    dplyr::filter(optima, metric == "mFM")
  )
)

