#### Threshold optimization ####

if (isTRUE(optimotu.pipeline::do_optimize_thresholds())) {
  unknown_ranks <- optimotu.pipeline::unknown_ranks()
  threshold_measure <- optimotu.pipeline::cluster_measure()
  if (is.null(threshold_measure)) {
    threshold_measure <- "FM"
  }
  threshold_memory_budget_mb <- optimotu.pipeline::cluster_memory_budget_mb()
  use_hamming <- identical(
    optimotu.pipeline::cluster_dist_config()$method,
    "hamming"
  )
  train_known_aligned <- isTRUE(
    optimotu.pipeline::do_optimize_thresholds_reference()
  ) &&
    (isTRUE(optimotu.pipeline::protax_aligned()) ||
      isTRUE(optimotu.pipeline::bayesant_aligned()) ||
      isTRUE(optimotu.pipeline::do_epa()))

  # Plan-time filesystem path of the train FASTA when we can resolve it
  # without a targets-only directory (Protax model dirs stay as targets).
  threshold_train_path <- if (
    isTRUE(optimotu.pipeline::do_optimize_thresholds_file())
  ) {
    optimotu.pipeline::optimize_thresholds_file()
  } else if (
    isTRUE(optimotu.pipeline::do_optimize_thresholds_reference()) &&
      !train_known_aligned
  ) {
    if (isTRUE(optimotu.pipeline::do_sintax())) {
      optimotu.pipeline::sintax_ref()
    } else if (isTRUE(optimotu.pipeline::do_bayesant())) {
      optimotu.pipeline::bayesant_ref()
    } else {
      NULL
    }
  } else {
    NULL
  }

  file_hamming_ready <- isTRUE(use_hamming) &&
    isTRUE(optimotu.pipeline::do_optimize_thresholds_file()) &&
    !is.null(threshold_train_path) &&
    isTRUE(optimotu.pipeline::fasta_is_hamming_ready(threshold_train_path))

  .reuse_seq_file <- NULL
  .reuse_seq_index <- NULL
  .reuse_need_index <- FALSE
  if (
    isTRUE(use_hamming) &&
      !isTRUE(optimotu.pipeline::do_optimize_thresholds_self()) &&
      !train_known_aligned &&
      !file_hamming_ready &&
      !is.null(threshold_train_path)
  ) {
    if (
      isTRUE(optimotu.pipeline::do_outgroup()) &&
        isTRUE(optimotu.pipeline::do_model_align()) &&
        optimotu.pipeline::is_same_path(
          threshold_train_path,
          optimotu.pipeline::outgroup_reference()
        )
    ) {
      .reuse_seq_file <- quote(outgroup_aligned)
      .reuse_seq_index <- quote(outgroup_aligned_index)
    } else if (
      isTRUE(optimotu.pipeline::do_supp_asv()) &&
        isTRUE(optimotu.pipeline::do_model_align())
    ) {
      set_names <- optimotu.pipeline::supp_asv_set_names()
      set_ids <- make.unique(make.names(set_names))
      for (i in seq_along(set_names)) {
        if (
          optimotu.pipeline::is_same_path(
            threshold_train_path,
            optimotu.pipeline::supp_asv_sequences(set_names[[i]])
          )
        ) {
          .reuse_seq_file <- rlang::sym(
            sprintf("supp_asv_aligned_seq_%s", set_ids[[i]])
          )
          .reuse_need_index <- TRUE
          break
        }
      }
    }
  }

  need_batch_align <- isTRUE(use_hamming) &&
    !isTRUE(optimotu.pipeline::do_optimize_thresholds_self()) &&
    !train_known_aligned &&
    !file_hamming_ready &&
    is.null(.reuse_seq_file)

  # Pointers so threshold_optima does not name targets that are absent in
  # other train_data modes.
  if (isTRUE(optimotu.pipeline::do_optimize_thresholds_self())) {
    if (isTRUE(use_hamming)) {
      .threshold_seq_file <- final_asv_seq
      .threshold_seq_index <- final_asv_seq_index
    } else {
      .threshold_seq_file <- final_asv_unaln_seq
      .threshold_seq_index <- final_asv_unaln_seq_index
    }
  } else if (!is.null(.reuse_seq_file)) {
    .threshold_seq_file <- .reuse_seq_file
    if (isTRUE(.reuse_need_index)) {
      .threshold_seq_index <- quote(threshold_refseq_hamming_index)
    } else {
      .threshold_seq_index <- .reuse_seq_index
    }
  } else if (isTRUE(need_batch_align)) {
    .threshold_seq_file <- quote(threshold_refseq_aligned)
    .threshold_seq_index <- quote(threshold_refseq_aligned_index)
  } else {
    .threshold_seq_file <- quote(threshold_refseq_file)
    .threshold_seq_index <- quote(threshold_refseq_index)
  }

  # Taxonomy source: branch at plan definition time so !! unquoting in one
  # mode cannot evaluate stop() paths from another mode.
  if (isTRUE(optimotu.pipeline::do_optimize_thresholds_self())) {
    if (isTRUE(optimotu.pipeline::do_supp_asv())) {
      threshold_reftax_target <- tar_fst_tbl(
        threshold_reftax,
        dplyr::bind_rows(
          asv_all_tax_prob,
          !!optimotu.pipeline::tar_map_bind_rows(
            supp_tax_file_plan,
            "supp_tax_prob"
          ),
          !!optimotu.pipeline::tar_map_bind_rows(
            supp_tax_header_plan,
            "supp_tax_prob"
          ),
          !!optimotu.pipeline::tar_map_bind_rows(
            supp_tax_classify_plan,
            "supp_all_tax_prob"
          )
        ) |>
          dplyr::filter(
            rank %in% !!unknown_ranks,
            prob >= !!optimotu.pipeline::cluster_min_conf(),
            taxon != "unk"
          ) |>
          dplyr::summarize(
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
    } else {
      threshold_reftax_target <- tar_fst_tbl(
        threshold_reftax,
        asv_all_tax_prob |>
          dplyr::filter(
            rank %in% !!unknown_ranks,
            prob >= !!optimotu.pipeline::cluster_min_conf(),
            taxon != "unk"
          ) |>
          dplyr::summarize(
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
    }
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
    if (isTRUE(.reuse_need_index)) {
      list(
        #### threshold_refseq_hamming_index ####
        # fastqindexr_index for reused supplemental aligned batches
        threshold_refseq_hamming_index = tar_target(
          threshold_refseq_hamming_index,
          fastqindexr::create_index(!!.reuse_seq_file),
          deployment = "main"
        )
      )
    },
    if (isTRUE(need_batch_align)) {
      list(
        #### threshold_refseq_batch ####
        # tibble: Hamming train-alignment batches (max_batchsize)
        #   batch, batch_id, from, to
        threshold_refseq_batch = tar_fst_tbl(
          threshold_refseq_batch,
          {
            n_seq <- optimotu.pipeline::sequence_size(threshold_refseq_file)
            n_seq <- sum(n_seq)
            n_batch <- as.integer(ceiling(
              n_seq / !!optimotu.pipeline::max_batchsize()
            ))
            batchsize <- as.integer(ceiling(n_seq / n_batch))
            tibble::tibble(
              batch = seq_len(n_batch),
              batch_id = optimotu.pipeline::make_seq_names(n_batch, "thref"),
              from = (batch - 1L) * batchsize + 1L,
              to = dplyr::lead(from, 1L, default = n_seq + 1L) - 1L
            )
          },
          resources = tar_resources(
            crew = tar_resources_crew(controller = "thin")
          )
        ),
        #### threshold_refseq_aligned ####
        # character: one batch of train sequences aligned to the amplicon
        # model (same helper as outgroup_aligned)
        threshold_refseq_aligned = tar_file(
          threshold_refseq_aligned,
          optimotu.pipeline::align_fasta_batch_to_amplicon_model(
            seqs = threshold_refseq_index,
            files = threshold_refseq_file,
            seq_idx = seq(
              threshold_refseq_batch$from,
              threshold_refseq_batch$to
            ),
            model = amplicon_model_file,
            model_type = !!optimotu.pipeline::amplicon_model_type(),
            outfile = file.path(
              !!optimotu.pipeline::aligned_path(),
              sprintf(
                "%s.fasta.gz",
                threshold_refseq_batch$batch_id
              )
            ),
            ncpu = optimotu.pipeline::local_cpus(),
            hmmalign = !!optimotu.pipeline::find_hmmalign()
          ),
          pattern = map(threshold_refseq_batch),
          resources = tar_resources(
            crew = tar_resources_crew(controller = "wide")
          )
        ),
        #### threshold_refseq_aligned_index ####
        # fastqindexr_index object for threshold_refseq_aligned batches
        threshold_refseq_aligned_index = tar_target(
          threshold_refseq_aligned_index,
          fastqindexr::create_index(threshold_refseq_aligned),
          resources = tar_resources(
            crew = tar_resources_crew(controller = "thin")
          )
        )
      )
    },
    list(
      #### threshold_optima ####
      # tibble: optimum thresholds for each rank / supertaxon / measure
      threshold_optima = tar_target(
        threshold_optima,
        {
          file_ids <- optimotu.pipeline::fasta_names(
            !!.threshold_seq_file
          )
          if (isTRUE(!!use_annotated_headers)) {
            file_ids <- optimotu.pipeline::fasta_header_seq_ids(file_ids)
          }
          seq_idx <- match(threshold_reftax$seq_id, file_ids)
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
            clust_config = !!optimotu.pipeline::cluster_clust_config(),
            parallel_config = !!(if (
              optimotu.pipeline::cluster_dist_config()$method == "usearch"
            ) {
              quote(optimotu::parallel_concurrent(2))
            } else {
              optimotu.pipeline::cluster_parallel_config(
                optimotu.pipeline::local_cpus()
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
