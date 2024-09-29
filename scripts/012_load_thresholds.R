# test the performance of different clustering thresholds for replicating
# "true" taxonomic groups.

if (isFALSE(optimize_thresholds)) {
  fmeasure_optima <- suppressWarnings(readr::read_tsv(
    threshold_file,
    col_types = list(
      rank = "c",
      superrank = "c",
      supertaxon = "c",
      conf_level = "c",
      threshold = "n",
      .default = "-"
    )
  ),
  "vroom_mismatched_column_name"
  )

  if ("conf_level" %in% names(fmeasure_optima)) {
    fmeasure_optima <- dplyr::filter(fmeasure_optima, conf_level == "plausible")
    fmeasure_optima$conf_level <- NULL
  }
}
