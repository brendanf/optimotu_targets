library(tarchetypes)

# Supplemental ASV pipeline
# Do phase 2 analysis for supplemental ASVs, and combine them with the native
# ASVs to produce a combined set for phase 3.

if (optimotu.pipeline::do_supp_asv()) {
  #### supp_meta ####
  # Table of supplemental sets.
  supp_meta <- tibble::tibble(
    .set_name = optimotu.pipeline::supp_asv_set_names(),
    .set_id = make.unique(make.names(.set_name)),
    .sequences = optimotu.pipeline::supp_asv_sequences(.set_name),
    .sample_table = optimotu.pipeline::supp_asv_sample_table(.set_name),
    .has_sample_table = !is.na(.sample_table),
    .taxonomy_mode = optimotu.pipeline::supp_asv_taxonomy_mode(.set_name),
    .taxonomy_file = optimotu.pipeline::supp_asv_taxonomy_file(.set_name)
  )

  #### supp_meta_sample ####
  # Table of supplemental sets with sample tables.
  supp_meta_sample <- dplyr::filter(supp_meta, .has_sample_table) |>
    dplyr::mutate(
      .supp_asv_names = rlang::syms(sprintf("supp_asv_names_%s", .set_id))
    )

  #### supp_meta_tax_file ####
  # Table of supplemental sets with taxonomy files.
  # Includes target symbol to map `supp_asv_names_{.set_id}` from
  # `supplemental_base_map`.
  supp_meta_tax_file <- dplyr::filter(supp_meta, .taxonomy_mode == "file") |>
    dplyr::mutate(
      .supp_asv_names = rlang::syms(sprintf("supp_asv_names_%s", .set_id))
    )

  #### supp_meta_tax_header ####
  # Table of supplemental sets where taxonomy is encoded in FASTA
  # headers.
  # Includes target symbols to map `supp_asv_names_{.set_id}` and
  # `supp_sequences_file_{.set_id}` from `supplemental_base_map`.
  supp_meta_tax_header <- dplyr::filter(
    supp_meta,
    .taxonomy_mode == "header"
  ) |>
    dplyr::mutate(
      .supp_asv_names = rlang::syms(sprintf("supp_asv_names_%s", .set_id)),
      .supp_sequences_file = rlang::syms(
        sprintf("supp_sequences_file_%s", .set_id)
      ),
    )

  #### supp_meta_tax_classify ####
  # Table of supplemental sets where taxonomy is not provided.
  # Includes target symbols to map `supp_asv_seq_{.set_id}`,
  # `supp_asv_seq_index_{.set_id}`, `supp_asv_seqbatch_{.set_id}`, and
  # `supp_asv_aligned_seq_{.set_id}` from `supplemental_base_map`.
  supp_meta_tax_classify <- dplyr::filter(
    supp_meta,
    .taxonomy_mode == "none"
  ) |>
    dplyr::mutate(
      .supp_asv_names = rlang::syms(
        sprintf("supp_asv_names_%s", .set_id)
      ),
      .supp_asv_seq = rlang::syms(
        sprintf("supp_asv_seq_%s", .set_id)
      ),
      .supp_asv_seq_index = rlang::syms(
        sprintf("supp_asv_seq_index_%s", .set_id)
      ),
      .supp_asv_seqbatch = rlang::syms(
        sprintf("supp_asv_seqbatch_%s", .set_id)
      ),
      .supp_asv_aligned_seq = rlang::syms(
        sprintf("supp_asv_aligned_seq_%s", .set_id)
      )
    )

  #### supp_base_plan ####
  # Static map over supplemental sets defined in `supp_meta`.
  # Defines per-set sequence metadata and standardized sequence files used in
  # downstream taxonomy and clustering.
  supp_base_plan <- tar_map(
    values = supp_meta,
    names = .set_id,
    list(
      ##### supp_sequences_file_{.set_id} #####
      # `character`: source sequence file path for one supplemental set.
      tar_file(supp_sequences_file, .sequences, deployment = "main"),

      ##### supp_set_info_{.set_id} #####
      # `tibble`:
      #   set_name: original supplemental set name from options.
      #   taxonomy_mode: one of "none", "file", "header".
      #   has_sample_table: whether this set provides read abundances.
      tar_fst_tbl(
        supp_set_info,
        tibble::tibble(
          set_name = .set_name,
          taxonomy_mode = .taxonomy_mode,
          has_sample_table = .has_sample_table
        ),
        deployment = "main"
      ),
      ##### supp_asv_names_{.set_id} #####
      # `tibble`:
      #   set_name: supplemental set name.
      #   seq_id_in: normalized source sequence ID for this set.
      #   seq_id: namespaced sequence ID used in combined outputs.
      tar_fst_tbl(
        supp_asv_names,
        {
          seq_ids <- names(Biostrings::fasta.seqlengths(supp_sequences_file))
          seq_ids <- sub("[[:space:];|].*$", "", seq_ids)
          if (any(!nzchar(seq_ids))) {
            stop(
              "Supplemental ASV set '",
              .set_name,
              "' contains FASTA headers with empty sequence IDs ",
              "after header normalization."
            )
          }
          seq_ids <- unique(seq_ids)
          tibble::tibble(
            set_name = .set_name,
            seq_id_in = seq_ids,
            # Always namespace supplemental IDs by set identifier to avoid
            # collisions.
            seq_id = sprintf("%s__%s", .set_id, seq_id_in)
          )
        },
        deployment = "main"
      ),
      ##### supp_asv_seq_{.set_id} #####
      # `character`: path to normalized, namespaced supplemental FASTA.gz.
      # Header normalization strips metadata and deduplicates identical IDs.
      tar_file(
        supp_asv_seq,
        {
          seqs <- Biostrings::readDNAStringSet(supp_sequences_file)
          seq_tbl <- tibble::tibble(
            seq_id_in = sub("[[:space:];|].*$", "", names(seqs)),
            seq = as.character(seqs)
          )
          if (any(!nzchar(seq_tbl$seq_id_in))) {
            stop(
              "Supplemental ASV set '",
              .set_name,
              "' contains FASTA headers with empty sequence IDs ",
              "after header normalization."
            )
          }
          duplicated_ids <- unique(seq_tbl$seq_id_in[duplicated(
            seq_tbl$seq_id_in
          )])
          if (length(duplicated_ids) > 0L) {
            conflicting_ids <- vapply(
              duplicated_ids,
              function(id) {
                length(unique(seq_tbl$seq[seq_tbl$seq_id_in == id])) > 1L
              },
              logical(1)
            )
            if (any(conflicting_ids)) {
              conflict <- duplicated_ids[which(conflicting_ids)[1]]
              stop(
                "Supplemental ASV set '",
                .set_name,
                "' has duplicate normalized sequence ID '",
                conflict,
                "' with different sequences."
              )
            }
            seq_tbl <- seq_tbl[!duplicated(seq_tbl$seq_id_in), , drop = FALSE]
          }
          seqs <- Biostrings::DNAStringSet(seq_tbl$seq)
          names(seqs) <- sprintf("%s__%s", .set_id, seq_tbl$seq_id_in)
          optimotu.pipeline::write_sequence(
            seqs,
            file.path(
              !!optimotu.pipeline::asv_path(),
              sprintf("supplemental_asv_%s.fasta.gz", .set_id)
            ),
            compress = TRUE
          )
        },
        deployment = "main"
      ),
      ##### supp_asv_seq_index_{.set_id} #####
      # `character`: index file for supp_asv_seq.
      tar_target(
        supp_asv_seq_index,
        fastqindexr::create_index(supp_asv_seq),
        deployment = "main"
      ),
      ##### supp_asv_seqbatch_{.set_id} #####
      # `tibble`: batching map for supplemental ASV sequences.
      tar_group_size(
        supp_asv_seqbatch,
        tibble::tibble(
          seq_idx = seq_len(optimotu.pipeline::sequence_size(supp_asv_seq))
        ),
        size = optimotu.pipeline::max_batchsize(),
        deployment = "main",
        format = "fst_tbl"
      ),
      if (isTRUE(optimotu.pipeline::do_model_align())) {
        list(
          ##### supp_asv_aligned_seq_{.set_id} #####
          # `character`: path to per-set aligned supplemental FASTA.gz.
          # Uses CM or HMM alignment according to `amplicon_model_type`.
          if (identical(optimotu.pipeline::amplicon_model_type(), "CM")) {
            tar_file(
              supp_asv_aligned_seq,
              inferrnal::cmalign(
                amplicon_model_file,
                fastqindexr::extract_sequences(
                  index = supp_asv_seq_index,
                  seq_idx = supp_asv_seqbatch$seq_idx,
                  file = supp_asv_seq,
                  return = "seq"
                ),
                global = TRUE,
                notrunc = TRUE,
                dnaout = TRUE,
                cpu = optimotu.pipeline::local_cpus()
              ) |>
                optimotu.pipeline::consensus_columns() |>
                optimotu.pipeline::write_sequence(
                  fname = file.path(
                    !!optimotu.pipeline::aligned_path(),
                    sprintf(
                      "supplemental_asv_aligned_%s_%d.fasta.gz",
                      .set_id,
                      supp_asv_seqbatch$tar_group[1]
                    )
                  ),
                  compress = TRUE
                ),
              pattern = map(supp_asv_seqbatch),
              resources = tar_resources(
                crew = tar_resources_crew(controller = "wide")
              )
            )
          } else if (
            identical(optimotu.pipeline::amplicon_model_type(), "HMM")
          ) {
            tar_file(
              supp_asv_aligned_seq,
              optimotu.pipeline::hmmalign(
                seqs = supp_asv_seq_index,
                hmm = amplicon_model_file,
                files = supp_asv_seq,
                seq_idx = supp_asv_seqbatch$seq_idx,
                outfile = file.path(
                  !!optimotu.pipeline::aligned_path(),
                  sprintf(
                    "supplemental_asv_aligned_%s_%d.fasta.gz",
                    .set_id,
                    supp_asv_seqbatch$tar_group[1]
                  )
                ),
                outformat = "A2M",
                hmmalign = !!optimotu.pipeline::find_hmmalign()
              ),
              pattern = map(supp_asv_seqbatch),
              resources = tar_resources(
                crew = tar_resources_crew(controller = "wide")
              )
            )
          } else {
            stop(
              "Invalid model type for supplemental alignment: ",
              optimotu.pipeline::amplicon_model_type()
            )
          }
        )
      } else {
        NULL
      }
    )
  )

  #### supp_sample_plan ####
  # Map over sets with sample tables defined in `supp_meta_sample` to parse
  # long-format supplemental reads.
  supp_sample_plan <- if (nrow(supp_meta_sample) > 0L) {
    tar_map(
      values = supp_meta_sample,
      names = .set_id,

      ##### supp_sample_table_file_{.set_id} #####
      # `character`: source sample-table file path for one set.
      tar_file(
        supp_sample_table_file,
        .sample_table,
        deployment = "main"
      ),

      ##### supp_sample_table_{.set_id} #####
      # `tibble`:
      #   sample, seq_id_in, nread, seqrun, set_name.
      tar_fst_tbl(
        supp_sample_table,
        optimotu.pipeline::read_long_sequence_table(
          path = supp_sample_table_file,
          default_seqrun = .set_name
        ) |>
          dplyr::rename(seq_id_in = seq_id) |>
          dplyr::left_join(.supp_asv_names, by = "seq_id_in") |>
          dplyr::select(seqrun, sample, seq_id, nread),
        deployment = "main"
      )
    )
  } else {
    list()
  }

  # Taxonomy targets which are required for downstream analyses are
  # supp_best_hit_taxon, supp_tax_prob, and supp_unknown_prob.
  # For sets where taxonomy is provided, these are parsed from the taxonomy
  # file. For sets where taxonomy is not provided, these are calculated using
  # the same taxonomic classifier as the main ASVs.

  #### supp_tax_file_plan ####
  # Map over sets with explicit taxonomy files defined in `supp_meta_tax_file`
  # to parse taxonomy.
  supp_tax_file_plan <- if (nrow(supp_meta_tax_file) > 0L) {
    tar_map(
      values = supp_meta_tax_file,
      names = .set_id,

      ##### supp_taxonomy_file_{.set_id} #####
      # `character`: source taxonomy file path for one set.
      tar_file(
        supp_taxonomy_file,
        .taxonomy_file,
        deployment = "main"
      ),

      ##### supp_best_hit_taxon_{.set_id} (file mode) #####
      # tibble:
      #   `seq_id` character: normalized source sequence ID.
      #   `{rank...}` character: taxonomic name at `rank`.
      # supp_best_hit_taxon and supp_tax_prob contain the same information
      # for "file mode" sets. Because supp_best_hit_taxon is in wide format and
      # has no probability column, it is more convenient to parse directly from
      # the taxonomy file.
      tar_fst_tbl(
        supp_best_hit_taxon,
        optimotu.pipeline::parse_reference_taxonomy(
          supp_taxonomy_file,
          ranks = !!optimotu.pipeline::tax_ranks()
        ) |>
          purrr::reduce2(
            !!optimotu.pipeline::known_ranks(),
            !!optimotu.pipeline::known_taxa(),
            .init = _,
            .f = \(d, rank, taxon) {
              if (!rank %in% names(d)) {
                d[[rank]] <- taxon
              }
              d
            }
          ) |>
          dplyr::rename(seq_id_in = seq_id) |>
          dplyr::mutate(set_name = .set_name) |>
          dplyr::left_join(.supp_asv_names, by = c("set_name", "seq_id_in")) |>
          dplyr::select(
            seq_id,
            !!!optimotu.pipeline::known_rank_vars()
          ),
        deployment = "main"
      ),

      ##### supp_unknown_prob_{.set_id} #####
      # tibble:
      #   `seq_id` character: normalized source sequence ID.
      #   `rank` factor: taxonomic rank.
      #  `novel_prob` numeric : cumulative probability that the ASV belongs to
      #    any novel taxon at `rank`. Always `NA` for reference classification.
      #  `known_prob` numeric : maximum probability that the ASV belongs to any
      #    one known taxon at `rank`. For reference-based classification, `1.0`
      #    if the ASV is assigned to a known taxon, `0.0` if not.
      #  `known_taxon` character : if `known_prob` is nonzero, the name of
      #    a known taxon which the ASV belongs to with probability
      #    `known_prob`. When `known_prob` < 0.5, it is possible for there to be
      #    more than one such taxon, but only one is given. `NA` if
      #    `known_prob` is 0.
      tar_fst_tbl(
        supp_unknown_prob,
        supp_best_hit_taxon |>
          dplyr::pivot_longer(
            cols = any_of(!!optimotu.pipeline::tax_ranks()),
            names_to = "rank",
            values_to = "known_taxon",
            names_transform = \(x) {
              optimotu.pipeline::rank2factor(
                x,
                !!optimotu.pipeline::tax_ranks()
              )
            }
          ) |>
          dplyr::transmute(
            seq_id,
            rank,
            novel_prob = NA_real_,
            known_prob = ifelse(is.na(known_taxon), 0.0, 1.0),
            known_taxon
          ),
        deployment = "main"
      ),

      ##### supp_tax_prob_{.set_id} (file mode) #####
      # `tibble`:
      #   seq_id (character): normalized source sequence ID.
      #   rank (factor): taxonomic rank.
      #   taxon (character): taxonomic name at `rank`.
      #   prob (numeric): probability that the ASV in `seq_id_in` belongs to
      #     `taxon` at `rank`.
      tar_fst_tbl(
        supp_tax_prob,
        dplyr::filter(supp_unknown_prob, !is.na(known_taxon)) |>
          dplyr::select(
            seq_id,
            rank,
            taxon = known_taxon,
            prob = known_prob
          ),
        deployment = "main"
      )
    )
  } else {
    list()
  }

  #### supp_tax_header_plan ####
  # Map over sets where taxonomy is encoded in FASTA headers defined in
  # `supp_meta_tax_header` to parse taxonomy.
  supp_tax_header_plan <- if (nrow(supp_meta_tax_header) > 0L) {
    tar_map(
      values = supp_meta_tax_header,
      names = .set_id,

      ##### supp_best_hit_taxon_{.set_id} (header mode) #####
      # `tibble`: parsed authoritative taxonomy keyed by seq_id.
      #
      # supp_best_hit_taxon and supp_all_tax_prob contain the same information
      # for "header mode" sets. Because supp_best_hit_taxon is in wide format
      # and has no probability column, it is more convenient to parse directly
      # from the FASTA headers.
      tar_fst_tbl(
        supp_best_hit_taxon,
        optimotu.pipeline::parse_reference_taxonomy(
          .supp_sequences_file,
          ranks = !!optimotu.pipeline::tax_ranks()
        ) |>
          purrr::reduce2(
            !!optimotu.pipeline::known_ranks(),
            !!optimotu.pipeline::known_taxa(),
            .init = _,
            .f = \(d, rank, taxon) {
              if (!rank %in% names(d)) {
                d[[rank]] <- taxon
              }
              d
            }
          ) |>
          dplyr::rename(seq_id_in = seq_id) |>
          dplyr::mutate(set_name = .set_name) |>
          dplyr::left_join(.supp_asv_names, by = c("set_name", "seq_id_in")) |>
          dplyr::select(
            seq_id,
            !!!optimotu.pipeline::tax_rank_vars()
          ),
        deployment = "main"
      ),

      ##### supp_unknown_prob_{.set_id} #####
      # tibble:
      #   `seq_id` character: normalized source sequence ID.
      #   `rank` factor: taxonomic rank.
      #   `novel_prob` numeric: cumulative probability that the ASV belongs to
      #     any novel taxon at `rank`. Always `NA` for reference-based
      #     classification.
      #   `known_prob` numeric: maximum probability that the ASV belongs to any
      #     one known taxon at `rank`. For reference-based classification, `1.0`
      #     if the ASV is assigned to a known taxon, `0.0` if not.
      #   `known_taxon` character: if `known_prob` is nonzero, the name of a
      #    known taxon which the ASV belongs to with probability `known_prob`.
      #    When `known_prob` < 0.5, it is possible for there to be more than
      #    one such taxon, but only one is given. `NA` if `known_prob` is 0.
      tar_fst_tbl(
        supp_unknown_prob,
        supp_best_hit_taxon |>
          tidyr::pivot_longer(
            cols = any_of(!!optimotu.pipeline::tax_ranks()),
            names_to = "rank",
            values_to = "known_taxon",
            names_transform = \(x) {
              optimotu.pipeline::rank2factor(
                x,
                !!optimotu.pipeline::tax_ranks()
              )
            }
          ) |>
          dplyr::transmute(
            seq_id,
            rank,
            novel_prob = NA_real_,
            known_prob = ifelse(is.na(known_taxon), 0.0, 1.0),
            known_taxon
          ),
        deployment = "main"
      ),

      ##### supp_all_tax_prob_{.set_id} #####
      # `tibble`: all taxonomic assignments for this set.
      tar_fst_tbl(
        supp_all_tax_prob,
        dplyr::filter(supp_unknown_prob, !is.na(known_taxon)) |>
          dplyr::select(
            seq_id,
            rank,
            taxon = known_taxon,
            prob = known_prob
          ),
        deployment = "main"
      )
    )
  } else {
    list()
  }

  #### supp_tax_classify_plan ####
  # Map over sets with taxonomy_mode "none" defined in `supp_meta_tax_classify`
  # to classify sequences.
  supp_tax_classify_plan <- if (nrow(supp_meta_tax_classify) > 0L) {
    tar_map(
      values = supp_meta_tax_classify,
      names = .set_id,
      ##### supp_all_tax_prob_{.set_id} #####
      # `tibble`:
      #   seq_id (character): unique ASV ID, as in `supp_asv_names_{.set_id}`.
      #   rank (factor): taxonomic rank.
      #   parent_taxonomy (character): comma-separated taxonomy of parent to
      #     `taxon`.
      #   taxon (character): taxonomic name at `rank`.
      #   prob (numeric): probability that the ASV in `seq_id` belongs to
      #     `taxon` at `rank`.
      # This target is calculated according to the same taxonomic classifier
      # as `asv_all_tax_prob` for the native ASVs.

      if (isTRUE(optimotu.pipeline::do_sintax())) {
        ###### supp_all_tax_prob_{.set_id} (sintax) ######
        tar_fst_tbl(
          supp_all_tax_prob,
          optimotu.pipeline::sintax(
            query = .supp_asv_seq_index,
            seq_idx = .supp_asv_seqbatch$seq_idx,
            files = .supp_asv_seq,
            ref = sintax_ref_file,
            ncpu = optimotu.pipeline::local_cpus(),
            id_is_int = FALSE,
            vsearch = !!optimotu.pipeline::find_vsearch()
          ),
          resources = tar_resources(
            crew = tar_resources_crew(controller = "wide")
          )
        )
      } else if (isTRUE(optimotu.pipeline::do_bayesant())) {
        if (optimotu.pipeline::bayesant_aligned()) {
          ###### supp_all_tax_prob_{.set_id} (bayesant aln) ######
          tar_fst_tbl(
            supp_all_tax_prob,
            optimotu.pipeline::bayesant(
              query = .supp_asv_aligned_seq,
              model = !!optimotu.pipeline::read_bayesant_model(),
              ncpu = optimotu.pipeline::local_cpus(),
              id_is_int = FALSE
            ),
            pattern = map(.supp_asv_aligned_seq),
            resources = tar_resources(
              crew = tar_resources_crew(controller = "wide")
            )
          )
        } else {
          ###### supp_all_tax_prob_{.set_id} (bayesant unaln) ######
          tar_fst_tbl(
            supp_all_tax_prob,
            optimotu.pipeline::bayesant(
              query = .supp_asv_seq_index,
              file = .supp_asv_seq,
              seq_idx = .supp_asv_seqbatch$seq_idx,
              model = !!optimotu.pipeline::read_bayesant_model(),
              ncpu = optimotu.pipeline::local_cpus(),
              id_is_int = FALSE,
              hash = .supp_asv_seqbatch_hash
            ),
            pattern = map(.supp_asv_seqbatch, .supp_asv_seqbatch_hash),
            resources = tar_resources(
              crew = tar_resources_crew(controller = "wide")
            )
          )
        }
      } else if (isTRUE(optimotu.pipeline::do_protax())) {
        if (isTRUE(optimotu.pipeline::protax_aligned())) {
          ###### supp_all_tax_prob_{.set_id} (protax aligned) ######
          tar_fst_tbl(
            supp_all_tax_prob,
            optimotu.pipeline::run_protax_animal(
              .supp_asv_aligned_seq,
              modeldir = protax_dir,
              id_is_int = FALSE,
              min_p = 0.02,
              info = TRUE,
              options = c("-m", "300"),
              executable = !!optimotu.pipeline::find_executable("classify_info")
            ) |>
              dplyr::transmute(
                seq_id,
                rank = optimotu.pipeline::rank2factor(
                  (!!optimotu.pipeline::unknown_ranks())[rank],
                  !!optimotu.pipeline::tax_ranks()
                ),
                parent_taxonomy = paste(
                  paste(!!optimotu.pipeline::known_taxa(), collapse = ","),
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
              ),
            pattern = map(.supp_asv_aligned_seq),
            resources = tar_resources(
              crew = tar_resources_crew(controller = "wide")
            )
          )
        } else {
          ###### supp_all_tax_prob_{.set_id} (protax unaln) ######
          tar_fst_tbl(
            supp_all_tax_prob,
            {
              outdir <- withr::local_tempdir()
              out_files <- optimotu.pipeline::run_protax(
                seqs = .supp_asv_seq,
                outdir = outdir,
                modeldir = protax_model,
                script = file.path(script_dir, "runprotax")
              )
              optimotu.pipeline::parse_protax_nameprob(
                grep("query\\d.nameprob", out_files, value = TRUE),
                id_is_int = FALSE
              ) |>
                dplyr::select(seq_id, rank, parent_taxonomy, taxon, prob)
            },
            resources = tar_resources(
              crew = tar_resources_crew(controller = "wide")
            )
          )
        }
      } else if (isTRUE(optimotu.pipeline::do_epa())) {
        ###### supp_all_tax_prob_{.set_id} (epa) ######
        # unlike the main ASVs, this is not split into targets for each step.
        tar_fst_tbl(
          supp_all_tax_prob,
          {
            outdir <- withr::local_tempdir()
            jplace <- optimotu.pipeline::epa_ng(
              ref_msa = !!optimotu.pipeline::epa_ref(),
              tree = !!optimotu.pipeline::epa_tree(),
              query = .supp_asv_aligned_seq,
              outdir = outdir,
              model = !!optimotu.pipeline::epa_params(),
              strip_inserts = TRUE,
              epa_ng = !!optimotu.pipeline::find_epa_ng()
            )
            optimotu.pipeline::gappa_assign(
              jplace = jplace,
              taxonomy = !!optimotu.pipeline::epa_taxonomy(),
              outgroup = !!optimotu.pipeline::epa_outgroup(),
              ncpu = !!optimotu.pipeline::local_cpus(),
              id_is_int = FALSE,
              gappa = !!optimotu.pipeline::find_gappa()
            )
          },
          pattern = map(.supp_asv_aligned_seq),
          resources = tar_resources(
            crew = tar_resources_crew(controller = "wide")
          )
        )
      } else {
        stop("No supported taxonomy classifier is active.")
      },
      ##### supp_best_hit_taxon_{.set_id} #####
      # `tibble`: best-hit taxonomy derived from authoritative references for
      # this set.
      #   seq_id (character): unique ASV ID, as in `supp_asv_names_{.set_id}`.
      #   {rank...} (character): taxonomic name at `rank`.
      #
      # This target is calculated according to the same search method as
      # `asv_best_hit_taxon` for the native ASVs.
      supp_best_hit_taxon = if (
        identical(optimotu.pipeline::cluster_dist_config()$method, "hamming")
      ) {
        ###### supp_best_hit_taxon_{.set_id} (hamming) ######
        tar_fst_tbl(
          supp_best_hit_taxon,
          optimotu::seq_search(
            query = .supp_asv_aligned_seq,
            ref = outgroup_aligned,
            threshold = 0.5,
            dist_config = !!(optimotu.pipeline::cluster_dist_config()$call),
            parallel_config = optimotu::parallel_concurrent(
              optimotu.pipeline::local_cpus()
            )
          ) |>
            dplyr::slice_min(dist, by = seq_id, with_ties = FALSE) |>
            dplyr::left_join(outgroup_taxonomy, by = "ref_id") |>
            dplyr::select(
              seq_id,
              !!!optimotu.pipeline::known_rank_vars()
            ),
          pattern = cross(.supp_asv_aligned_seq, outgroup_aligned),
          resources = tar_resources(
            crew = tar_resources_crew(controller = "wide")
          )
        )
      } else {
        ###### supp_best_hit_taxon_{.set_id} (other) ######
        tar_fst_tbl(
          supp_best_hit_taxon,
          optimotu::seq_search(
            query = fastqindexr::extract_sequences(
              index = .supp_asv_seq_index,
              seq_idx = .supp_asv_seqbatch$seq_idx,
              file = .supp_asv_seq,
              return = "seq"
            ),
            ref = unaligned_ref_seqs,
            threshold = 0.5,
            dist_config = !!(optimotu.pipeline::cluster_dist_config()$call),
            parallel_config = optimotu::parallel_concurrent(
              optimotu.pipeline::local_cpus()
            )
          ) |>
            dplyr::slice_min(dist, by = seq_id, with_ties = FALSE) |>
            dplyr::left_join(outgroup_taxonomy, by = "ref_id") |>
            dplyr::select(
              seq_id,
              !!!optimotu.pipeline::known_rank_vars()
            ),
          pattern = map(.supp_asv_seqbatch),
          resources = tar_resources(
            crew = tar_resources_crew(controller = "wide")
          )
        )
      },

      ##### supp_tax_prob_{.set_id} #####
      # tibble:
      #   `set_name` character: supplemental set name.
      #   `seq_id_in` character: normalized source sequence ID.
      #   `rank` factor: taxonomic rank.
      #   `taxon` character: taxonomic name at `rank`.
      #   `prob` numeric: probability that the ASV in `seq_id_in` belongs to
      #     `taxon` at `rank`.
      tar_fst_tbl(
        supp_tax_prob,
        tidyr::crossing(
          seq_id = .supp_asv_names$seq_id[.supp_asv_seqbatch$seq_idx],
          tibble::tibble(
            rank = optimotu.pipeline::rank2factor(
              !!optimotu.pipeline::known_ranks(),
              !!optimotu.pipeline::tax_ranks()
            ),
            taxon = !!optimotu.pipeline::known_taxa(),
            prob = 1.0
          ) |>
            dplyr::arrange(dplyr::desc(rank)) |>
            dplyr::mutate(
              parent_taxonomy = purrr::accumulate(
                .x = dplyr::lag(taxon, default = NA_character_),
                .f = \(x, y) if (is.na(x)) y else paste(x, y, sep = ",")
              )
            )
        ) |>
          dplyr::full_join(
            supp_all_tax_prob,
            by = c("seq_id", "rank", "taxon", "parent_taxonomy", "prob")
          ) |>
          dplyr::arrange(dplyr::desc(prob)) |>
          dplyr::summarize(
            taxon = dplyr::first(taxon),
            prob = dplyr::first(prob),
            .by = c(seq_id, rank)
          ),
        pattern = map(.supp_asv_seqbatch, supp_all_tax_prob),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      ),

      ##### supp_unknown_prob_{.set_id} #####
      # tibble:
      #   `seq_id` character: normalized source sequence ID.
      #   `rank` factor: taxonomic rank.
      #   `novel_prob` numeric: cumulative probability that the ASV belongs to
      #     any novel taxon at `rank`.
      #   `known_prob` numeric: maximum probability that the ASV belongs to any
      #     one known taxon at `rank`.
      #   `known_taxon` character: if `known_prob` is nonzero, the name of a
      #     known taxon which the ASV belongs to with probability `known_prob`.
      #     When `known_prob` < 0.5, it is possible for there to be more than
      #     one such taxon, but only one is given. `NA` if `known_prob` is 0.
      tar_fst_tbl(
        supp_unknown_prob,
        supp_all_tax_prob |>
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
          dplyr::filter(!rank %in% !!optimotu.pipeline::known_ranks()) |>
          dplyr::arrange(seq_id, desc(rank)) |>
          dplyr::mutate(
            novel_prob = cumsum(novel_prob),
            .by = seq_id
          ),
        pattern = map(supp_all_tax_prob),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      )
    )
  } else {
    list()
  }

  supp_tax_plan <- purrr::reduce(
    list(
      supp_tax_file_plan,
      supp_tax_header_plan,
      supp_tax_classify_plan
    ),
    optimotu.pipeline::tar_merge
  )

  #### bind_supp_{} helpers ####
  # Symbolic aggregations of per-set targets consumed by downstream combined
  # targets.
  bind_supp_set_info <- optimotu.pipeline::tar_map_bind_rows(
    supp_base_plan,
    "supp_set_info"
  )
  bind_supp_asv_names <- optimotu.pipeline::tar_map_bind_rows(
    supp_base_plan,
    "supp_asv_names"
  )
  bind_supp_asv_seq <- optimotu.pipeline::tar_map_c(
    supp_base_plan,
    "supp_asv_seq"
  )

  bind_supp_asv_aligned_seq <- if (
    isTRUE(!!optimotu.pipeline::do_model_align())
  ) {
    optimotu.pipeline::tar_map_c(
      supp_base_plan,
      "supp_asv_aligned_seq"
    )
  } else {
    quote(character())
  }

  bind_supp_sample_table <- if (length(supp_sample_plan) > 0) {
    optimotu.pipeline::tar_map_bind_rows(
      supp_sample_plan,
      "supp_sample_table"
    )
  } else {
    quote(tibble::tibble(
      sample = character(),
      seq_id_in = character(),
      nread = integer(),
      seqrun = character(),
      set_name = character()
    ))
  }

  bind_supp_tax_prob <- optimotu.pipeline::tar_map_bind_rows(
    supp_tax_plan,
    "supp_tax_prob"
  )

  bind_supp_best_hit_taxon <- optimotu.pipeline::tar_map_bind_rows(
    supp_tax_plan,
    "supp_best_hit_taxon"
  )

  bind_supp_unknown_prob <- optimotu.pipeline::tar_map_bind_rows(
    supp_tax_plan,
    "supp_unknown_prob"
  )

  supp_plan <- c(
    supp_base_plan,
    supp_sample_plan,
    supp_tax_plan,
    list(
      ###### supp_set_info ######
      # `tibble`: union of per-set metadata from supp_base_plan.
      tar_fst_tbl(
        supp_set_info,
        !!bind_supp_set_info,
        deployment = "main"
      ),
      ###### supp_asv_names ######
      # `tibble`: union of per-set normalized and namespaced ASV IDs.
      tar_fst_tbl(
        supp_asv_names,
        !!bind_supp_asv_names,
        deployment = "main"
      ),
      ###### combo_asv_table ######
      # `tibble`: supplemental long-format read table in namespaced seq_id
      # space.
      tar_fst_tbl(
        combo_asv_table,
        dplyr::bind_rows(
          asv_table,
          !!bind_supp_sample_table,
        ),
        deployment = "main"
      ),
      ###### combo_tax_prob ######
      # `tibble`: union of classified supplemental taxonomy probabilities.
      tar_fst_tbl(
        combo_tax_prob,
        dplyr::bind_rows(
          asv_tax_prob,
          !!bind_supp_tax_prob,
        ),
        deployment = "main"
      ),

      ###### combo_best_hit_taxon ######
      # `tibble`: union of best-hit taxonomy tables from supp_tax_plan.
      tar_fst_tbl(
        combo_best_hit_taxon,
        dplyr::bind_rows(
          asv_best_hit_taxon,
          !!bind_supp_best_hit_taxon,
        ),
        deployment = "main"
      ),

      ###### combo_unknown_prob ######
      # `tibble`: union of unknown probability tables from supp_tax_plan.
      tar_fst_tbl(
        combo_unknown_prob,
        dplyr::bind_rows(
          asv_unknown_prob,
          !!bind_supp_unknown_prob,
        ),
        deployment = "main"
      ),

      ###### combo_seq ######
      # `character` vector : paths to all ASV sequences.
      tar_file(
        combo_seq,
        c(asv_seq, !!bind_supp_asv_seq),
        deployment = "main"
      ),

      ###### combo_seq_index ######
      # `character`: index file for combined_asv_seq.
      tar_target(
        combo_seq_index,
        fastqindexr::create_index(combo_seq),
        deployment = "main"
      )
    ),
    if (isTRUE(optimotu.pipeline::do_model_align())) {
      list(
        ###### combo_aligned_seq ######
        # `character`: paths to all ASV sequences aligned to the model.
        tar_file(
          combo_aligned_seq,
          c(asv_aligned_seq, !!bind_supp_asv_aligned_seq),
          deployment = "main"
        ),
        ###### combo_aligned_seq_index ######
        # `character`: index file for combo_aligned_seq.
        tar_target(
          combo_aligned_seq_index,
          fastqindexr::create_index(combo_aligned_seq),
          deployment = "main"
        )
      )
    } else {
      list()
    }
  )

  optimotu_plan <- c(optimotu_plan, supp_plan)

  #### "final_asv" pointers ####
  final_asv_unaln_seq <- quote(combo_seq)
  final_asv_unaln_seq_index <- quote(combo_seq_index)
  if (isTRUE(optimotu.pipeline::do_model_align())) {
    final_asv_seq <- quote(combo_aligned_seq)
    final_asv_seq_index <- quote(combo_aligned_seq_index)
  } else {
    final_asv_seq <- quote(combo_seq)
    final_asv_seq_index <- quote(combo_seq_index)
  }
  final_asv_tax_prob <- quote(combo_tax_prob)
  final_asv_table <- quote(combo_asv_table)
  final_asv_best_hit_taxon <- quote(combo_best_hit_taxon)
  final_asv_unknown_prob <- quote(combo_unknown_prob)
} else {
  final_asv_unaln_seq <- quote(asv_seq)
  final_asv_unaln_seq_index <- quote(asv_seq_index)
  if (isTRUE(optimotu.pipeline::do_model_align())) {
    final_asv_seq <- quote(aligned_seq)
    final_asv_seq_index <- quote(aligned_seq_index)
  } else {
    final_asv_seq <- quote(asv_seq)
    final_asv_seq_index <- quote(asv_seq_index)
  }
  final_asv_tax_prob <- quote(asv_tax_prob)
  final_asv_table <- quote(asv_table)
  final_asv_best_hit_taxon <- quote(asv_best_hit_taxon)
  final_asv_unknown_prob <- quote(asv_unknown_prob)
}
