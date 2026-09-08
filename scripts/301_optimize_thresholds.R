#### Threshold optimization ####

if (isTRUE(optimotu.pipeline::do_optimize_thresholds())) {
  unknown_ranks <- optimotu.pipeline::unknown_ranks()
  threshold_measure <- optimotu.pipeline::cluster_measure()
  if (is.null(threshold_measure)) {
    threshold_measure <- "FM"
  }
  threshold_memory_budget_mb <- optimotu.pipeline::cluster_memory_budget_mb()

  # Pointers so threshold_optima does not name targets that are absent in
  # other train_data modes.
  if (isTRUE(optimotu.pipeline::do_optimize_thresholds_self())) {
    .threshold_seq_file <- quote(asv_seq)
    .threshold_seq_index <- quote(asv_seq_index)
  } else {
    .threshold_seq_file <- quote(threshold_refseq_file)
    .threshold_seq_index <- quote(threshold_refseq_index)
  }

  # Taxonomy source: branch at plan definition time so !! unquoting in one
  # mode cannot evaluate stop() paths from another mode.
  if (isTRUE(optimotu.pipeline::do_optimize_thresholds_self())) {
    threshold_reftax_target <- tar_fst_tbl(
      threshold_reftax,
      asv_all_tax_prob |>
        dplyr::filter(
          rank %in% !!unknown_ranks,
          prob >= !!optimotu.pipeline::cluster_min_conf(),
          taxon != "unk"
        ) |>
        dplyr::summarize(
          prob = max(prob),
          taxon = taxon[which.max(prob)],
          .by = c(seq_id, rank)
        ) |>
        tidyr::pivot_wider(
          names_from = rank,
          values_from = taxon,
          names_expand = TRUE
        ) |>
        dplyr::select(seq_id, dplyr::any_of(!!unknown_ranks)) |>
        optimotu::clean_taxonomy(ranks = !!unknown_ranks),
      deployment = "main"
    )
  } else if (
    isTRUE(optimotu.pipeline::do_optimize_thresholds_reference()) &&
      isTRUE(optimotu.pipeline::protax_aligned())
  ) {
    # protaxA has no single reference taxonomy file in our format;
    # assemble it from taxonomy.priors and model.rseqs.numeric.
    threshold_reftax_target <- tar_fst_tbl(
      threshold_reftax,
      dplyr::left_join(
        readr::read_tsv(
          file.path(protax_dir, "model.rseqs.numeric"),
          col_names = c("tax_id", "seq_id"),
          col_types = "i-c"
        ) |>
          tidyr::separate_longer_delim(seq_id, " ") |>
          dplyr::mutate(seq_id = as.integer(seq_id) + 1L),
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
              !!unknown_ranks
            )
          ) |>
          dplyr::select(tax_id, rank, taxon),
        by = "tax_id"
      ) |>
        tidyr::pivot_wider(
          id_cols = seq_id,
          names_from = rank,
          values_from = taxon
        ) |>
        dplyr::arrange(seq_id) |>
        dplyr::select(seq_id, dplyr::any_of(!!unknown_ranks)) |>
        optimotu::clean_taxonomy(ranks = !!unknown_ranks),
      deployment = "main"
    )
  } else {
    tax_src <- if (isTRUE(optimotu.pipeline::do_optimize_thresholds_file())) {
      quote(threshold_train_file)
    } else if (isTRUE(optimotu.pipeline::protax_unaligned())) {
      quote(file.path(protax_model, "seqid2tax"))
    } else if (isTRUE(optimotu.pipeline::do_sintax())) {
      quote(sintax_ref_file)
    } else if (isTRUE(optimotu.pipeline::do_bayesant())) {
      if (is.null(optimotu.pipeline::bayesant_ref())) {
        stop(
          "No taxonomic reference file available for threshold ",
          "optimization.\n",
          "For \"train_data: 'reference'\" you must supply\n",
          "  taxonomy:\n",
          "    bayesant:\n",
          "      reference:\n",
          "(instead of or in addition to \"model:\")\n",
          "in 'pipeline_options.yaml'"
        )
      }
      quote(bayesant_ref_file)
    } else if (isTRUE(optimotu.pipeline::do_epa())) {
      quote(epa_taxonomy_file)
    } else {
      stop("Cannot determine which taxonomic classifier is in use.")
    }
    threshold_reftax_target <- tar_fst_tbl(
      threshold_reftax,
      optimotu.pipeline::parse_reference_taxonomy(
        !!tax_src,
        ranks = !!unknown_ranks
      ) |>
        optimotu::clean_taxonomy(ranks = !!unknown_ranks),
      deployment = "main"
    )
  }

  if (!isTRUE(optimotu.pipeline::do_optimize_thresholds_self())) {
    seq_src <- if (isTRUE(optimotu.pipeline::do_optimize_thresholds_file())) {
      quote(threshold_train_file)
    } else if (isTRUE(optimotu.pipeline::protax_aligned())) {
      quote(file.path(protax_dir, "refs.aln"))
    } else if (isTRUE(optimotu.pipeline::protax_unaligned())) {
      quote(file.path(protax_model, "its2.fa"))
    } else if (isTRUE(optimotu.pipeline::do_sintax())) {
      quote(sintax_ref_file)
    } else if (isTRUE(optimotu.pipeline::do_bayesant())) {
      quote(bayesant_ref_file)
    } else if (isTRUE(optimotu.pipeline::do_epa())) {
      quote(epa_ref_file)
    } else {
      stop("Cannot determine reference sequence source.")
    }
    # Annotated FASTA (SINTAX/BayesANT/custom file): pass seq_names into
    # optimize_thresholds(). Protax/EPA already use bare seq_id headers.
    threshold_refseq_file_target <- tar_file(
      threshold_refseq_file,
      !!seq_src,
      deployment = "main"
    )
    threshold_refseq_index_target <- tar_target(
      threshold_refseq_index,
      fastqindexr::create_index(threshold_refseq_file),
      deployment = "main"
    )
  }

  # Plan-time flag: annotated headers need fasta_header_seq_ids() + seq_names.
  use_annotated_headers <- !isTRUE(
    optimotu.pipeline::do_optimize_thresholds_self()
  ) &&
    (isTRUE(optimotu.pipeline::do_optimize_thresholds_file()) ||
      isTRUE(optimotu.pipeline::do_sintax()) ||
      isTRUE(optimotu.pipeline::do_bayesant()))

  threshold_plan <- c(
    if (isTRUE(optimotu.pipeline::do_optimize_thresholds_file())) {
      list(
        #### threshold_train_file ####
        # character: annotated reference FASTA used as train_data
        threshold_train_file = tar_file(
          threshold_train_file,
          !!optimotu.pipeline::optimize_thresholds_file(),
          deployment = "main"
        )
      )
    },
    list(
      #### threshold_reftax ####
      # tibble: taxonomy used for threshold optimization
      #   seq_id: IDs of sequences in the training FASTA / ASV set
      #   {rank}: identification of each seq at each unknown rank
      threshold_reftax = threshold_reftax_target
    ),
    if (!isTRUE(optimotu.pipeline::do_optimize_thresholds_self())) {
      list(
        #### threshold_refseq_file ####
        # character: training FASTA (original reference; may have annotated
        # headers — seq_id matching is handled in threshold_optima)
        threshold_refseq_file = threshold_refseq_file_target,
        #### threshold_refseq_index ####
        # fastqindexr_index object for threshold_refseq_file
        threshold_refseq_index = threshold_refseq_index_target
      )
    },
    list(
      #### threshold_optima ####
      # tibble: optimum thresholds for each rank / supertaxon / measure
      threshold_optima = tar_target(
        threshold_optima,
        {
          if (isTRUE(!!use_annotated_headers)) {
            file_ids <- optimotu.pipeline::fasta_header_seq_ids(
              names(Biostrings::fasta.seqlengths(!!.threshold_seq_file))
            )
            seq_idx <- match(threshold_reftax$seq_id, file_ids)
          } else {
            seq_ids <- names(Biostrings::fasta.seqlengths(
              !!.threshold_seq_file
            ))
            seq_idx <- match(threshold_reftax$seq_id, seq_ids)
          }
          if (anyNA(seq_idx)) {
            stop(
              "Some threshold_reftax$seq_id values are missing from the ",
              "training sequence file."
            )
          }
          optimotu::optimize_thresholds(
            taxonomy = threshold_reftax,
            refseq = !!.threshold_seq_index,
            seq_file = !!.threshold_seq_file,
            seq_idx = seq_idx,
            seq_names = if (isTRUE(!!use_annotated_headers)) {
              threshold_reftax$seq_id
            },
            ranks = !!unknown_ranks,
            dist_config = !!(if (
              optimotu.pipeline::cluster_dist_config()$method == "usearch"
            ) {
              substitute(
                update(dc, usearch_ncpu = optimotu.pipeline::local_cpus()),
                list(dc = optimotu.pipeline::cluster_dist_config()$call)
              )
            } else {
              optimotu.pipeline::cluster_dist_config()$call
            }),
            threshold_config = optimotu::threshold_uniform(
              0.0,
              !!optimotu.pipeline::cluster_dist_max(),
              !!optimotu.pipeline::cluster_dist_step()
            ),
            parallel_config = !!(if (
              optimotu.pipeline::cluster_dist_config()$method == "usearch"
            ) {
              quote(optimotu::parallel_concurrent(2))
            } else {
              quote(
                optimotu::parallel_concurrent(optimotu.pipeline::local_cpus())
              )
            }),
            min_taxa = !!as.integer(optimotu.pipeline::cluster_min_taxa()),
            min_refseq = !!as.integer(optimotu.pipeline::cluster_min_refseq()),
            measures = !!threshold_measure,
            clustering_memory_budget_mb = !!threshold_memory_budget_mb
          )
        },
        resources = tar_resources(
          crew = tar_resources_crew(controller = "wide")
        )
      ),

      #### cluster_optima ####
      # tibble: measure-filtered optima consumed by scripts/302_cluster.R
      cluster_optima = tar_fst_tbl(
        cluster_optima,
        {
          out <- threshold_optima
          if ("measure" %in% names(out)) {
            out <- dplyr::filter(out, measure == !!threshold_measure)
          } else if ("metric" %in% names(out)) {
            out <- dplyr::filter(out, metric == !!threshold_measure)
          }
          out
        },
        deployment = "main"
      )
    ),
    if (!is.null(optimotu.pipeline::cluster_thresholds())) {
      list(
        #### optima_file ####
        # character: optional TSV of optimized thresholds for reuse
        optima_file = tar_file(
          optima_file,
          optimotu.pipeline::write_and_return_file(
            threshold_optima,
            !!optimotu.pipeline::cluster_thresholds(),
            "tsv"
          ),
          deployment = "main"
        )
      )
    }
  )

  optimotu_plan <- c(optimotu_plan, threshold_plan)
}
