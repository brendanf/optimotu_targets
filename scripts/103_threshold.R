# test the performance of different clustering thresholds for replicating
# "true" taxonomic groups.

threshold_test_plan <- list(
  tar_file(
    ref_db_file,
    "protaxFungi/addedmodel/sintaxits2train.fa",
    deployment = "main"
  ),
  tar_fst_tbl(
    ref_db,
    Biostrings::readDNAStringSet(ref_db_file) %>%
      as.character() %>%
      tibble::enframe(value = "seq") %>%
      tidyr::extract(name, c("seq_id", "tax"), regex = "^([^;]+);tax=d:(.+)") %>%
      tidyr::separate(tax, into = TAXRANKS, sep = ",[pcofgs]:", fill = "right")
  ),
  tar_fst_tbl(
    ref_tax,
    tidyr::pivot_longer(
      ref_db,
      cols = kingdom:species,
      names_to = "rank",
      values_to = "taxon",
      values_drop_na = TRUE
    )
  ),
  
  tar_map(
    values = list(
      rank = TAXRANKS[-1]
    ),
    tar_target(
      testset_select,
        purrr::map_dfr(
          superranks(rank),
          function(r, data) dplyr::group_by(data, dplyr::across(r)) %>%
            dplyr::summarize(
              superrank = r,
              n_taxa = dplyr::n_distinct(.data[[rank]]),
              n_seq = dplyr::n_distinct(seq_id),
              seq_id = list(seq_id)
            ) %>%
            dplyr::rename(supertaxon = !!r),
          data = dplyr::filter(ref_tax, .data[["rank"]] == rank)
        ) %>%
        dplyr::filter(n_taxa >= 10 | superrank == "kingdom", n_seq >= 50),
      tidy_eval = FALSE
    ),
    tar_target(
      threshold_testset,
      blastclust_repeat(
        seqs = ref_db$seq,
        seqnames = ref_db$seq_id,
        threshold = seq(60, 99.9, 0.1),
        threshold_name = seq(600, 999, 1),
        which = unlist(testset_select$seq_id),
        hitlist_method = "usearch",
        usearch = "bin/usearch"
      ) %>%
        purrr::map(trimws),
      pattern = map(testset_select),
      iteration = "list"
    ),
    tar_fst_tbl(
      fmeasures,
      purrr::map(
        threshold_testset,
        tibble::enframe,
        name = "cluster",
        value = "seq_id"
      ) %>%
        purrr::map_dfr(
          tidyr::separate_rows,
          "seq_id",
          sep = " ",
          .id = "threshold"
        ) %>%
        tidyr::pivot_wider(
          names_from = threshold,
          values_from = cluster
        ) %>%
        dplyr::right_join(
          dplyr::select(ref_tax, seq_id, rank),
          .,
          by = "seq_id"
        ) %>%
        purrr::map_dbl(
          set_names(names(.), names(.))[-1:-2],
          f_measure,
          data = .,
          c = rank
        ) %>%
        tibble::enframe(name = "threshold", value = "f_measure") %>%
        dplyr::mutate(
          rank = rank,
          superrank = testset_select$superrank,
          supertaxon = testset_select$supertaxon,
          threshold = as.numeric(threshold)/10,
          conf_level = conf_level
        ),
      pattern = map(testset_select, threshold_testset)
    )
  )
)

threshold_test_plan <- c(
  threshold_test_plan,
  tar_combine(
    fmeasures,
    threshold_test_plan[[4]][startsWith(names(threshold_test_plan[[4]]), "fmeasures")]
  ),
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
