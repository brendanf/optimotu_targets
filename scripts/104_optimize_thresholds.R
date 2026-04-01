threshold_plan <- list(
  tar_fst_tbl(
    threshold_meta,
    tibble::tibble(
      rank_int = seq_along(TAXRANKS)[-1],
      rank = TAXRANKS[rank_int]
    ),
    deployment = "main"
  ),
  tar_fst_tbl(
    threshold_reftax,
    dplyr::filter(asv_all_tax_prob, prob > 0.5, taxon != "unk") |>
      tidyr::pivot_wider(
        values_from = "taxon",
        names_from = "rank",
        id_cols = "seq_id"
      ) |>
      optimotu::clean_taxonomy(ranks = !!optimotu.pipeline::unknown_ranks()),
    deployment = "main"
  ),
  tar_target(
    threshold_optima,
    optimotu::optimize_thresholds(
      taxonomy = threshold_reftax,
      refseq = optimotu.pipeline::select_sequence(asv_seq, threshold_reftax$seq_id),
      ranks = !!optimotu.pipeline::unknown_ranks(),
      dist_config = !!(
        if (optimotu.pipeline::cluster_dist_config()$method == "usearch") {
          substitute(
            update(dc, usearch_ncpu = optimotu.pipeline::local_cpus()),
            list(dc = optimotu.pipeline::cluster_dist_config())
          )
        } else {
          optimotu.pipeline::cluster_dist_config()
        }
      ),
      threshold_config = optimotu::threshold_uniform(0.0, 0.4, 0.001),
      parallel_config = !!(
        if (optimotu.pipeline::cluster_dist_config()$method == "usearch") {
          quote(optimotu::parallel_concurrent(2))
        } else {
          quote(optimotu::parallel_concurrent(optimotu.pipeline::local_cpus()))
        }
      )
    ),
    tar_resources(crew = tar_resources_crew(controller = "wide"))
  ),
  tar_file(
      optima_file,
      write_and_return_file(
        optima, file.path("output", "GSSP_thresholds.tsv"),
        "tsv"),
      deployment = "main"
  ),
  tar_fst_tbl(
    cluster_optima,
    dplyr::filter(threshold_optima, metric == "FM")
  )
)

