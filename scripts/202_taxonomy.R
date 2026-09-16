# Taxonomy assignment
#
# Dispatches to the primary classifier (and optional reclassify fallback).
# When taxonomy.reclassify is enabled, primary_all_tax_prob runs on all ASVs,
# fallback_all_tax_prob runs only on sequences that fail the confidence gate
# or the hard-coded best_hit scope screen, and all_tax_prob combines them.
# Otherwise all_tax_prob is the primary classifier output as before.

library(tarchetypes)

if (optimotu.pipeline::do_epa()) {
  epa_path <- "data/intermediate/epa"
  if (!dir.exists(epa_path)) {
    dir.create(epa_path, recursive = TRUE)
  }
}

.primary_name <- optimotu.pipeline::primary_classifier()
.do_reclassify <- optimotu.pipeline::do_reclassify()
.reclassify_name <- optimotu.pipeline::reclassify_classifier()
.do_best_hit <- optimotu.pipeline::do_best_hit()
.tax_prob_primary <- if (.do_reclassify) {
  "primary_all_tax_prob"
} else {
  "all_tax_prob"
}

.reclassify_ids_pattern <- local({
  p <- quote(map(primary_all_tax_prob, seqbatch))
  if (.do_best_hit) {
    p[[length(p) + 1L]] <- quote(best_hit_taxon)
  }
  p
})

.reclassify_ids_command <- rlang::expr({
  optimotu.pipeline::ids_needing_reclassify(
    tax_prob = primary_all_tax_prob,
    rank = !!optimotu.pipeline::reclassify_rank(),
    min_prob = !!optimotu.pipeline::reclassify_min_prob(),
    all_ids = seqbatch$seq_idx,
    id_col = "seq_idx",
    screen = if (!!.do_best_hit) {
      dplyr::semi_join(best_hit_taxon, seqbatch, by = "seq_idx")
    } else {
      NULL
    },
    scope_ranks = !!optimotu.pipeline::classifier_scope_ranks("primary"),
    scope_taxa = !!optimotu.pipeline::classifier_scope_taxa("primary")
  )
})

#### classifier setup targets (any configured slot) ####
taxonomy_setup <- c(
  if (optimotu.pipeline::do_protax()) {
    list(
      #### protax_dir ####
      # character : directory name
      protax_dir = tar_file(
        protax_dir,
        !!optimotu.pipeline::protax_location(),
        deployment = "main"
      )
    )
  },
  if (optimotu.pipeline::do_protax() && optimotu.pipeline::protax_unaligned()) {
    list(
      #### protax_script ####
      protax_script = tar_file(
        protax_script,
        file.path(script_dir, "runprotax"),
        deployment = "main"
      )
    )
  },
  if (optimotu.pipeline::do_sintax()) {
    list(
      #### sintax_ref_file ####
      sintax_ref_file = tar_file(
        sintax_ref_file,
        !!optimotu.pipeline::sintax_ref(),
        deployment = "main"
      )
    )
  },
  if (optimotu.pipeline::do_bayesant()) {
    list(
      if (!is.null(optimotu.pipeline::bayesant_ref())) {
        list(
          #### bayesant_ref_file ####
          bayesant_ref_file = tar_file(
            bayesant_ref_file,
            !!optimotu.pipeline::bayesant_ref(),
            deployment = "main"
          )
        )
      },
      bayesant_model = if (is.null(optimotu.pipeline::bayesant_model())) {
        tar_target(
          bayesant_model,
          BayesANT::read.BayesANT.data(
            fasta.file = bayesant_ref_file,
            rank = !!length(optimotu.pipeline::assigned_ranks(
              optimotu.pipeline::classifier_scope_ranks("primary")
            )),
            rank_names = !!optimotu.pipeline::assigned_ranks(
              optimotu.pipeline::classifier_scope_ranks("primary")
            )
          ) |>
            BayesANT::BayesANT(
              typeseq = !!(if (optimotu.pipeline::bayesant_aligned()) {
                "aligned"
              } else {
                " not aligned"
              })
            ),
          resources = tar_resources(
            crew = tar_resources_crew(controller = "wide")
          )
        )
      } else if (file.exists(optimotu.pipeline::bayesant_model())) {
        tar_file(
          bayesant_model,
          !!optimotu.pipeline::bayesant_model(),
          deployment = "main"
        )
      } else {
        tar_file(
          bayesant_model,
          BayesANT::BayesANT(
            BayesANT::read.BayesANT.data(
              fasta.file = bayesant_ref_file,
              rank = !!length(optimotu.pipeline::assigned_ranks(
                optimotu.pipeline::classifier_scope_ranks("primary")
              )),
              rank_names = !!optimotu.pipeline::assigned_ranks(
                optimotu.pipeline::classifier_scope_ranks("primary")
              )
            ),
            typeseq = !!(if (optimotu.pipeline::bayesant_aligned()) {
              "aligned"
            } else {
              " not aligned"
            })
          ) |>
            optimotu.pipeline::write_and_return_file(
              file = !!optimotu.pipeline::bayesant_model()
            ),
          resources = tar_resources(
            crew = tar_resources_crew(controller = "wide")
          )
        )
      }
    )
  },
  if (optimotu.pipeline::do_epa()) {
    list(
      #### epa_ref_file ####
      tar_file(
        epa_ref_file,
        !!optimotu.pipeline::epa_ref(),
        deployment = "main"
      ),
      #### epa_taxonomy_file ####
      tar_file(
        epa_taxonomy_file,
        !!optimotu.pipeline::epa_taxonomy(),
        deployment = "main"
      ),
      #### epa_tree_file ####
      tar_file(
        epa_tree_file,
        !!optimotu.pipeline::epa_tree(),
        deployment = "main"
      ),
      if (file.exists(optimotu.pipeline::epa_params())) {
        tar_file(
          epa_params,
          !!optimotu.pipeline::epa_params(),
          deployment = "main"
        )
      } else {
        tar_target(
          epa_params,
          optimotu.pipeline::epa_params(),
          deployment = "main"
        )
      },
      if (file.exists(optimotu.pipeline::epa_outgroup())) {
        tar_file(
          epa_outgroup,
          !!optimotu.pipeline::epa_outgroup(),
          deployment = "main"
        )
      } else {
        tar_target(
          epa_outgroup,
          optimotu.pipeline::epa_outgroup(),
          deployment = "main"
        )
      }
    )
  }
)

#### primary classifier ####
primary_taxonomy <- switch(
  .primary_name,
  protax = if (optimotu.pipeline::protax_aligned()) {
    list(
      tar_target_raw(
        .tax_prob_primary,
        rlang::expr({
          optimotu.pipeline::run_protax_animal(
            aln_seqs = seq_model_align,
            modeldir = protax_dir,
            id_is_int = TRUE,
            min_p = 0.02,
            info = TRUE,
            options = c("-m", "300"),
            ncpu = optimotu.pipeline::local_cpus(),
            executable = !!optimotu.pipeline::find_executable("classify_info")
          ) |>
            dplyr::transmute(
              seq_idx,
              rank = optimotu.pipeline::rank2factor(
                (!!optimotu.pipeline::assigned_ranks(
                  optimotu.pipeline::classifier_scope_ranks("primary")
                ))[rank],
                !!optimotu.pipeline::tax_ranks()
              ),
              parent_taxonomy = paste(
                paste(
                  !!optimotu.pipeline::classifier_scope_taxa("primary"),
                  collapse = ","
                ),
                taxonomy,
                sep = ","
              ) |>
                sub(",[^,]+$", "", x = _),
              taxon = sub(".*,", "", taxonomy),
              prob,
              best_id,
              best_dist,
              second_id,
              second_dist
            )
        }),
        pattern = quote(map(seq_model_align)),
        format = "fst_tbl",
        resources = tar_resources(
          crew = tar_resources_crew(controller = "wide")
        )
      )
    )
  } else {
    list(
      protax = tar_file(
        protax,
        withr::with_tempfile(
          "tempout",
          fileext = ".fasta",
          {
            protax_dir
            protax_script
            optimotu.pipeline::run_protax(
              seqs = optimotu.pipeline::fastx_gz_extract(
                infile = !!seq_all_trim,
                index = !!seq_index_file,
                i = seqbatch$seq_idx,
                outfile = tempout,
                hash = seqbatch_hash
              ),
              outdir = file.path(
                !!optimotu.pipeline::protax_path(),
                tar_name()
              ),
              modeldir = protax_model,
              script = protax_script
            )
          }
        ),
        pattern = map(seqbatch, seqbatch_hash),
        iteration = "list",
        resources = tar_resources(
          crew = tar_resources_crew(controller = "wide")
        )
      ),
      tar_target_raw(
        .tax_prob_primary,
        quote(
          grep("query\\d.nameprob", protax, value = TRUE) |>
            optimotu.pipeline::parse_protax_nameprob(id_is_int = TRUE)
        ),
        pattern = quote(map(protax)),
        format = "fst_tbl",
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      )
    )
  },
  sintax = list(
    tar_target_raw(
      .tax_prob_primary,
      rlang::expr({
        optimotu.pipeline::sintax(
          query = optimotu.pipeline::fastx_gz_extract(
            infile = !!seq_all_trim,
            index = !!seq_index_file,
            i = seqbatch$seq_idx,
            outfile = withr::local_tempfile(fileext = ".fasta"),
            hash = seqbatch_hash
          ),
          ref = sintax_ref_file,
          ncpu = local_cpus(),
          id_is_int = TRUE,
          vsearch = !!optimotu.pipeline::find_vsearch()
        )
      }),
      pattern = quote(map(seqbatch, seqbatch_hash)),
      format = "fst_tbl",
      resources = tar_resources(
        crew = tar_resources_crew(controller = "wide")
      )
    )
  ),
  bayesant = list(
    tar_target_raw(
      .tax_prob_primary,
      if (optimotu.pipeline::bayesant_aligned()) {
        rlang::expr({
          optimotu.pipeline::bayesant(
            query = seq_model_align,
            model = !!optimotu.pipeline::read_bayesant_model(),
            ncpu = local_cpus(),
            id_is_int = TRUE
          )
        })
      } else {
        rlang::expr({
          optimotu.pipeline::bayesant(
            query = !!seq_index_file,
            file = seq_all_trim_file,
            seq_idx = seqbatch$seq_idx,
            model = !!optimotu.pipeline::read_bayesant_model(),
            ncpu = local_cpus(),
            id_is_int = TRUE,
            hash = seqbatch_hash
          )
        })
      },
      pattern = if (optimotu.pipeline::bayesant_aligned()) {
        quote(map(seq_model_align))
      } else {
        quote(map(seqbatch, seqbatch_hash))
      },
      format = "fst_tbl",
      resources = tar_resources(
        crew = tar_resources_crew(controller = "wide")
      )
    )
  ),
  epa = list(
    epa_ng = tar_file(
      epa_ng,
      optimotu.pipeline::epa_ng(
        ref_msa = epa_ref_file,
        tree = epa_tree_file,
        query = seq_model_align,
        outdir = file.path(epa_path, tar_name()),
        model = epa_params,
        strip_inserts = TRUE,
        epa_ng = !!optimotu.pipeline::find_epa_ng()
      ),
      pattern = map(seq_model_align),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "wide")
      )
    ),
    tar_target_raw(
      .tax_prob_primary,
      rlang::expr({
        optimotu.pipeline::gappa_assign(
          jplace = epa_ng,
          taxonomy = epa_taxonomy_file,
          outgroup = epa_outgroup,
          ranks = !!optimotu.pipeline::assigned_ranks(
            optimotu.pipeline::classifier_scope_ranks("primary")
          ),
          ncpu = optimotu.pipeline::local_cpus(),
          id_is_int = TRUE,
          gappa = !!optimotu.pipeline::find_gappa()
        ) |>
          dplyr::rename(parent_taxonomy = parent_taxon)
      }),
      pattern = quote(map(epa_ng)),
      format = "fst_tbl",
      resources = tar_resources(
        crew = tar_resources_crew(controller = "wide")
      )
    )
  ),
  stop("No taxonomy assignment method selected")
)

#### fallback classifier (optional cascade) ####
fallback_taxonomy <- if (!.do_reclassify) {
  list()
} else {
  reclassify_ids_target <- tar_target_raw(
    "reclassify_ids",
    .reclassify_ids_command,
    pattern = .reclassify_ids_pattern,
    iteration = "list",
    resources = tar_resources(
      crew = tar_resources_crew(controller = "thin")
    )
  )
  combine_all_tax_prob <- tar_fst_tbl(
    all_tax_prob,
    optimotu.pipeline::combine_tax_prob(
      primary = primary_all_tax_prob,
      fallback = fallback_all_tax_prob,
      reclassify_ids = reclassify_ids,
      id_col = "seq_idx",
      primary_name = !!.primary_name,
      fallback_name = !!.reclassify_name
    ),
    pattern = map(
      primary_all_tax_prob,
      fallback_all_tax_prob,
      reclassify_ids
    ),
    resources = tar_resources(
      crew = tar_resources_crew(controller = "thin")
    )
  )
  switch(
    .reclassify_name,
    epa = list(
      reclassify_ids_target,
      #### epa_ng_fallback ####
      epa_ng_fallback = tar_file(
        epa_ng_fallback,
        {
          need <- reclassify_ids
          query <- if (length(need) == 0L) {
            character()
          } else {
            aln <- Biostrings::readBStringSet(seq_model_align)
            aln[as.character(need)]
          }
          optimotu.pipeline::epa_ng(
            ref_msa = epa_ref_file,
            tree = epa_tree_file,
            query = query,
            outdir = file.path(epa_path, paste0(tar_name(), "_fallback")),
            model = epa_params,
            strip_inserts = TRUE,
            epa_ng = !!optimotu.pipeline::find_epa_ng()
          )
        },
        pattern = map(seq_model_align, reclassify_ids),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "wide")
        )
      ),
      #### fallback_all_tax_prob ####
      fallback_all_tax_prob = tar_fst_tbl(
        fallback_all_tax_prob,
        optimotu.pipeline::gappa_assign(
          jplace = epa_ng_fallback,
          taxonomy = epa_taxonomy_file,
          outgroup = epa_outgroup,
          ranks = !!optimotu.pipeline::assigned_ranks(
            optimotu.pipeline::classifier_scope_ranks("reclassify")
          ),
          ncpu = optimotu.pipeline::local_cpus(),
          id_is_int = TRUE,
          gappa = !!optimotu.pipeline::find_gappa()
        ) |>
          dplyr::rename(parent_taxonomy = parent_taxon),
        pattern = map(epa_ng_fallback),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "wide")
        )
      ),
      combine_all_tax_prob
    ),
    sintax = list(
      reclassify_ids_target,
      fallback_all_tax_prob = tar_fst_tbl(
        fallback_all_tax_prob,
        {
          need <- reclassify_ids
          if (length(need) == 0L) {
            tibble::tibble(
              seq_idx = integer(),
              rank = optimotu.pipeline::rank2factor(
                character(),
                !!optimotu.pipeline::tax_ranks()
              ),
              parent_taxonomy = character(),
              taxon = character(),
              prob = numeric()
            )
          } else {
            optimotu.pipeline::sintax(
              query = optimotu.pipeline::fastx_gz_extract(
                infile = !!seq_all_trim,
                index = !!seq_index_file,
                i = need,
                outfile = withr::local_tempfile(fileext = ".fasta"),
                hash = seqbatch_hash
              ),
              ref = sintax_ref_file,
              ncpu = local_cpus(),
              id_is_int = TRUE,
              vsearch = !!optimotu.pipeline::find_vsearch()
            )
          }
        },
        pattern = map(seqbatch, seqbatch_hash, reclassify_ids),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "wide")
        )
      ),
      combine_all_tax_prob
    ),
    protax = if (optimotu.pipeline::protax_aligned()) {
      list(
        reclassify_ids_target,
        fallback_all_tax_prob = tar_fst_tbl(
          fallback_all_tax_prob,
          {
            need <- reclassify_ids
            if (length(need) == 0L) {
              tibble::tibble(
                seq_idx = integer(),
                rank = optimotu.pipeline::rank2factor(
                  character(),
                  !!optimotu.pipeline::tax_ranks()
                ),
                parent_taxonomy = character(),
                taxon = character(),
                prob = numeric(),
                best_id = character(),
                best_dist = numeric(),
                second_id = character(),
                second_dist = numeric()
              )
            } else {
              aln <- Biostrings::readBStringSet(seq_model_align)
              query_file <- withr::local_tempfile(fileext = ".fasta")
              Biostrings::writeXStringSet(aln[as.character(need)], query_file)
              optimotu.pipeline::run_protax_animal(
                aln_seqs = query_file,
                modeldir = protax_dir,
                id_is_int = TRUE,
                min_p = 0.02,
                info = TRUE,
                options = c("-m", "300"),
                ncpu = optimotu.pipeline::local_cpus(),
                executable = !!optimotu.pipeline::find_executable(
                  "classify_info"
                )
              ) |>
                dplyr::transmute(
                  seq_idx,
                  rank = optimotu.pipeline::rank2factor(
                    (!!optimotu.pipeline::assigned_ranks(
                      optimotu.pipeline::classifier_scope_ranks("reclassify")
                    ))[rank],
                    !!optimotu.pipeline::tax_ranks()
                  ),
                  parent_taxonomy = paste(
                    paste(
                      !!optimotu.pipeline::classifier_scope_taxa("reclassify"),
                      collapse = ","
                    ),
                    taxonomy,
                    sep = ","
                  ) |>
                    sub(",[^,]+$", "", x = _),
                  taxon = sub(".*,", "", taxonomy),
                  prob,
                  best_id,
                  best_dist,
                  second_id,
                  second_dist
                )
            }
          },
          pattern = map(seq_model_align, reclassify_ids),
          resources = tar_resources(
            crew = tar_resources_crew(controller = "wide")
          )
        ),
        combine_all_tax_prob
      )
    } else {
      stop(
        "Unaligned Protax as a reclassify fallback is not implemented yet; ",
        "use aligned Protax, SINTAX, or EPA-ng."
      )
    },
    bayesant = stop(
      "BayesANT as a reclassify fallback is not implemented yet; ",
      "use Protax, SINTAX, or EPA-ng."
    ),
    stop("Unknown reclassify classifier: ", .reclassify_name)
  )
}

#### post-process ####
.fallback_label <- if (.do_reclassify) .reclassify_name else ""
taxonomy_post <- list(
  #### asv_all_tax_prob ####
  # Inject closed-scope priors only for sequences won by a closed classifier.
  asv_all_tax_prob = tar_fst_tbl(
    asv_all_tax_prob,
    {
      tax <- all_tax_prob
      ids <- seqbatch$seq_idx
      if ("classifier" %in% names(tax)) {
        fallback_ids <- unique(
          tax$seq_idx[tax$classifier == !!.fallback_label]
        )
        primary_ids <- setdiff(ids, fallback_ids)
        priors <- dplyr::bind_rows(
          if (!optimotu.pipeline::classifier_is_open("primary")) {
            optimotu.pipeline::scope_prior_taxonomy(
              primary_ids,
              id_col = "seq_idx",
              scope_ranks = !!optimotu.pipeline::classifier_scope_ranks(
                "primary"
              ),
              scope_taxa = !!optimotu.pipeline::classifier_scope_taxa(
                "primary"
              ),
              ranks = !!optimotu.pipeline::tax_ranks()
            )
          },
          if (
            isTRUE(!!.do_reclassify) &&
              !optimotu.pipeline::classifier_is_open("reclassify")
          ) {
            optimotu.pipeline::scope_prior_taxonomy(
              fallback_ids,
              id_col = "seq_idx",
              scope_ranks = !!optimotu.pipeline::classifier_scope_ranks(
                "reclassify"
              ),
              scope_taxa = !!optimotu.pipeline::classifier_scope_taxa(
                "reclassify"
              ),
              ranks = !!optimotu.pipeline::tax_ranks()
            )
          }
        )
      } else if (!optimotu.pipeline::classifier_is_open("primary")) {
        priors <- optimotu.pipeline::scope_prior_taxonomy(
          ids,
          id_col = "seq_idx",
          scope_ranks = !!optimotu.pipeline::classifier_scope_ranks("primary"),
          scope_taxa = !!optimotu.pipeline::classifier_scope_taxa("primary"),
          ranks = !!optimotu.pipeline::tax_ranks()
        )
      } else {
        priors <- optimotu.pipeline::scope_prior_taxonomy(
          integer(),
          id_col = "seq_idx",
          scope_ranks = character(),
          scope_taxa = character(),
          ranks = !!optimotu.pipeline::tax_ranks()
        )
      }
      dplyr::full_join(
        priors,
        tax,
        by = c("seq_idx", "rank", "taxon", "parent_taxonomy", "prob")
      ) |>
        dplyr::inner_join(asv_names, by = "seq_idx") |>
        dplyr::select(seq_id, everything() & !seq_idx) |>
        dplyr::arrange(seq_id, dplyr::desc(rank), dplyr::desc(prob))
    },
    pattern = map(seqbatch, all_tax_prob),
    resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
  ),

  #### asv_tax_prob ####
  asv_tax_prob = tar_fst_tbl(
    asv_tax_prob,
    {
      x <- asv_all_tax_prob |>
        dplyr::summarize(
          prob = max(prob, na.rm = TRUE),
          .by = c(taxon, rank, seq_id)
        ) |>
        dplyr::arrange(dplyr::desc(prob)) |>
        dplyr::summarize(
          taxon = dplyr::first(taxon),
          prob = dplyr::first(prob),
          .by = c(seq_id, rank)
        )
      if ("classifier" %in% names(asv_all_tax_prob)) {
        cls <- asv_all_tax_prob |>
          dplyr::arrange(dplyr::desc(prob)) |>
          dplyr::summarize(
            classifier = dplyr::first(classifier),
            .by = c(seq_id, rank)
          )
        x <- dplyr::left_join(x, cls, by = c("seq_id", "rank"))
      }
      x
    },
    pattern = map(asv_all_tax_prob),
    resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
  ),

  #### asv_unknown_prob ####
  asv_unknown_prob = tar_fst_tbl(
    asv_unknown_prob,
    asv_all_tax_prob |>
      dplyr::summarize(
        novel_prob = sum(prob[taxon == "unk"]),
        known_prob = max(prob[taxon != "unk"], 0),
        known_taxon = if (any(!is.na(taxon)) && known_prob > 0) {
          taxon[taxon != "unk" & prob == known_prob][1]
        } else {
          NA_character_
        },
        .by = c(seq_id, rank)
      ) |>
      tidyr::complete(
        seq_id,
        rank,
        fill = list(novel_prob = 0, known_prob = 0)
      ) |>
      dplyr::filter(
        !rank %in% !!optimotu.pipeline::known_ranks()
      ) |>
      dplyr::arrange(seq_id, desc(rank)) |>
      dplyr::mutate(
        novel_prob = cumsum(novel_prob),
        .by = seq_id
      ),
    pattern = map(asv_all_tax_prob),
    resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
  )
)

taxonomy_plan <- c(
  taxonomy_setup,
  primary_taxonomy,
  fallback_taxonomy,
  taxonomy_post
)

optimotu_plan <- c(optimotu_plan, taxonomy_plan)
