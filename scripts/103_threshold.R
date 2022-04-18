# test the performance of different clustering thresholds for replicating
# "true" taxonomic groups.

threshold_test_plan <- list(
  tar_fst_tbl(
    fmeasure_optima,
    fmeasures %>%
      dplyr::group_by(rank, superrank, supertaxon, conf_level) %>%
      dplyr::arrange(dplyr::desc(f_measure)) %>%
      dplyr::summarize(
        threshold = threshold[which.max(f_measure)],
        f_measure = max(f_measure),
        .groups = "drop"
      )
  )
)
