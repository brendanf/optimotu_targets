# test the performance of different clustering thresholds for replicating
# "true" taxonomic groups.

if (isFALSE(optimotu.pipeline::do_optimize_thresholds())) {
  cluster_optima <- suppressWarnings(
    readr::read_tsv(
      optimotu.pipeline::cluster_thresholds(),
      col_types = list(
        rank = "c",
        superrank = "c",
        supertaxon = "c",
        conf_level = "c",
        threshold = "n",
        measure = "c",
        .default = "-"
      )
    ),
    "vroom_mismatched_column_name"
  )

  if ("conf_level" %in% names(cluster_optima)) {
    cluster_optima <- dplyr::filter(cluster_optima, conf_level == "plausible")
    cluster_optima$conf_level <- NULL
  }
}
