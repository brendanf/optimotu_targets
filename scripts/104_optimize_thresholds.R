library(tarchetypes)

if (isTRUE(optimize_thresholds)) {

  threshold_plan <- list(

    #### threshold_meta ####
    # tibble:
    #  rank: ranks which we will cluster at
    #  rank_int: integer representation of the ranks
    tar_fst_tbl(
      threshold_meta,
      tibble::tibble(
        rank = UNKNOWN_RANKS,
        rank_int = as.integer(rank2factor(rank))
      ),
      deployment = "main"
    ),
    if (identical(optimize_train_data, "self")) {
      #### "self" ####
      list(
        ##### reftax #####
        tar_fst_tbl(
          reftax,
          dplyr::filter(
            asv_all_tax_prob,
            prob > optimize_min_conf,
            taxon != "unk"
          ) |>
            tidyr::pivot_wider(
              values_from = "taxon",
              names_from = "rank",
              id_cols = "seq_id"
            ) |>
            dplyr::bind_cols(as.list(`names<-`(KNOWN_TAXA, KNOWN_RANKS))) |>
            dplyr::select(seq_id, dplyr::any_of(TAX_RANKS)),
          deployment = "main"
        )
      )
      } else {
        list(
          if (identical(optimize_train_data, "reference")) {
            #### "reference" ####
            ###### reftax_file ######
            tar_file(
              reftax_file,
              !!(if (protax_aligned) {
                file.path(protax_dir, "refs.aln")
              } else if (protax_unaligned) {
                quote(file.path(protax_model, "sintaxits2.fasta"))
              } else if (do_sintax) {
                quote(sintax_ref)
              }),
              deployment = "main"
            )
          } else {
            #### explicit file name ####
            ##### reftax_file #####
            tar_file(
              reftax_file,
              optimize_train_data,
              deployment = "main"
            )
          },
          #### reftax ####
          tar_fst_tbl(
            reftax,
            parse_fasta_taxonomy(reftax_file) |>
              dplyr::bind_cols(as.list(`names<-`(KNOWN_TAXA, KNOWN_RANKS))) |>
              dplyr::select(seq_id, dplyr::any_of(TAX_RANKS)),
            deployment = "main"
          )
        )
      },

    ##### testset_select #####
    tar_target(
      testset_select,
      purrr::map_dfr(
        superranks(threshold_meta$rank),
        summarize_by_rank,
        data = reftax,
        rank = threshold_meta$rank
      ) |>
        dplyr::filter(n_taxa >= 5 | superrank == INGROUP_RANK, n_seq >= 10),
      tidy_eval = FALSE,
      pattern = map(threshold_meta), # per rank
      deployment = "main"
    ),

    ##### testset_rowwise #####
    tar_target(
      testset_rowwise,
      testset_select,
      deployment = "main"
    ),

    ##### threshold_testset #####
    tar_target(
      threshold_testset,
      optimotu::seq_cluster(
        seq = !!(
          if (identical(optimize_train_data, "self")) {
            quote(asv_seq)
          } else {
            quote(reftax_file)
          }
        ),
        dist_config = dist_config,
        threshold_config = optimotu::threshold_uniform(
          from = 0,
          to = optimize_max_dist,
          by = optimize_dist_step,
          thresh_names = as.character(seq(0, optimize_max_dist, optimize_dist_step))
        ),
        clust_config = optimotu::clust_slink(),
        parallel_config = optimotu::parallel_merge(local_cpus()),
        which = testset_select$seq_id
      ),
      iteration = "list",
      resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
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
            tibble::remove_rownames() |>
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
      resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
    ),
    tar_fst_tbl(
      optima,
      cluster_metrics |>
        tidyr::pivot_longer(
          -c(threshold, rank, superrank, supertaxon),
          names_to = "metric", values_to = "score"
        ) |>
        dplyr::filter(!(metric %in% c("TP", "FP", "FN", "TN", "EMI"))) |>
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
        optima, file.path("output", "GSSP_thresholds.tsv"),
        "tsv"),
      deployment = "main"
    ),
    tar_fst_tbl(
      fmeasure_optima,
      dplyr::filter(optima, metric == "mFM"),
      deployment = "main"
    )
  )

  optimotu_plan <- c(optimotu_plan, threshold_plan)
}

