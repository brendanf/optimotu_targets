if (isTRUE(optimotu.pipeline::do_optimize_thresholds())) {
  threshold_plan <- list(
    tar_fst_tbl(
      threshold_meta,
      tibble::tibble(
        rank_int = seq_along(TAXRANKS)[-1],
        rank = TAXRANKS[rank_int]
      ),
      deployment = "main"
    ),
    #### threshold_reftax ####
    # tibble: taxonomy to use for threshold optimization
    #   seq_id: IDs of sequences in threshold_refseq
    #   {rank}: identification of each seq at each rank.
    threshold_reftax = if (optimotu.pipeline::do_optimize_thresholds_self()) {
      tar_fst_tbl(
        threshold_reftax,
        dplyr::filter(
          asv_all_tax_prob,
          prob >= !!optimotu.pipeline::optimize_min_conf(),
          taxon != "unk"
        ) |>
          tidyr::pivot_wider(
            values_from = "taxon",
            names_from = "rank",
            id_cols = "seq_id"
          ) |>
          optimotu::clean_taxonomy(
            ranks = !!optimotu.pipeline::unknown_ranks()
          ),
        deployment = "main"
      )
    } else if (optimotu.pipeline::do_optimize_thresholds_reference()) {
      if (optimotu.pipeline::protax_aligned()) {
        # protaxA does not have a reference taxonomy file in the format we
        # need, so we need to assemble it from `taxonomy.priors`` and
        # `model.rseqs.numeric`
        tar_fst_tbl(
          threshold_reftax,
          dplyr::left_join(
            # mapping from tax_id to all reference sequences in that taxon
            readr::read_tsv(
              file.path(protax_dir, "model.rseqs.numeric"),
              col_names = c("tax_id", "seq_id"),
              col_types = "i-c"
            ) |>
              tidyr::separate_longer_delim(seq_id, " ") |>
              dplyr::mutate(seq_id = as.integer(seq_id) + 1L),
            # mapping from tax_id to taxonomy
            readr::read_tsv(
              file.path(protax_dir, "taxonomy.priors"),
              col_names = c("tax_id", "parent_id", "rank", "taxon"),
              col_types = c(
                tax_id = "i",
                parent_id = "-",
                rank = "i",
                taxon = "c",
                .default = "-"
              ),
              skip = 1
            ) |>
              dplyr::mutate(
                taxon = gsub(".+,", "", taxon),
                rank = optimotu.pipeline::int2rankfactor(
                  rank,
                  !!optimotu.pipeline::unknown_ranks()
                )
              ) |>
              dplyr::select(tax_id, taxon, rank),
            by = "tax_id"
          ) |>
            tidyr::pivot_wider(
              names_from = rank,
              values_from = taxon,
              id_cols = seq_id
            ) |>
            dplyr::arrange(seq_id) |>
            optimotu::clean_taxonomy(
              ranks = !!optimotu.pipeline::unknown_ranks()
            ),
          deployment = "main"
        )
      } else {
        tar_fst_tbl(
          threshold_reftax,
          optimotu.pipeline::parse_reference_taxonomy(
            !!if (optimotu.pipeline::protax_unaligned()) {
              quote(file.path(protax_model, "seqid2tax"))
            } else if (optimotu.pipeline::do_sintax()) {
              quote(sintax_ref_file)
            } else if (optimotu.pipeline::do_bayesant) {
              if (is.null(optimotu.pipeline::bayesant_ref())) {
                stop(
                  "No taxonomic reference file available for threshold optimization.\n",
                  "For \"train_data: 'reference'\" you must supply\n",
                  "  taxonomy:\n",
                  "    bayesant:\n",
                  "      reference:\n",
                  "(instead of or in addition to \"model:\")\n",
                  "in 'pipeline_options.yaml'"
                )
              } else {
                quote(bayesant_ref_file)
              }
            } else if (optimotu.pipeline::do_epa()) {
              quote(epa_taxonomy_file)
            } else {
              stop("Cannot determine which taxonomic classifier is in use.")
            },
            ranks = !!optimotu.pipeline::unknown_ranks()
          ) |>
            optimotu::clean_taxonomy(
              ranks = !!optimotu.pipeline::unknown_ranks()
            ),
          deployment = "main"
        )
      }
    } else {
      tar_fst_tbl(
        threshold_reftax,
        optimotu.pipeline::tar_fst_tbl(
          threshold_reftax,
          !!if (optimotu.pipeline::do_protax()) {
            if (optimotu.pipeline::protax_aligned()) {
              file.path(optimotu.pipeline::protax_location(), "refs.aln")
            } else {
              file.path(
                optimotu.pipeline::protax_location(),
                "addedmodel",
                "sintaxits2.fasta"
              )
            }
          }
        )
      )
    },
    tar_target(
      threshold_optima,
      optimotu::optimize_thresholds(
        taxonomy = threshold_reftax,
        refseq = optimotu.pipeline::select_sequence(
          asv_seq,
          threshold_reftax$seq_id
        ),
        ranks = !!optimotu.pipeline::unknown_ranks(),
        dist_config = !!(if (
          optimotu.pipeline::cluster_dist_config()$method == "usearch"
        ) {
          substitute(
            update(dc, usearch_ncpu = optimotu.pipeline::local_cpus()),
            list(dc = optimotu.pipeline::cluster_dist_config())
          )
        } else {
          optimotu.pipeline::cluster_dist_config()
        }),
        threshold_config = optimotu::threshold_uniform(0.0, 0.4, 0.001),
        parallel_config = !!(if (
          optimotu.pipeline::cluster_dist_config()$method == "usearch"
        ) {
          quote(optimotu::parallel_concurrent(2))
        } else {
          quote(optimotu::parallel_concurrent(optimotu.pipeline::local_cpus()))
        })
      ),
      tar_resources(crew = tar_resources_crew(controller = "wide"))
    ),
    tar_file(
      optima_file,
      write_and_return_file(
        optima,
        file.path("output", "GSSP_thresholds.tsv"),
        "tsv"
      ),
      deployment = "main"
    ),
    tar_fst_tbl(
      cluster_optima,
      dplyr::filter(threshold_optima, metric == "FM")
    )
  )
}
