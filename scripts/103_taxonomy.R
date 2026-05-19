# Taxonomy assignment

# This script defines the targets for taxonomy assignment.
#
# The first part of the plan dispatches to the appropriate taxonomy assignment
# algorithm, defining `all_tax_prob` as the main output. For some classifiers
# this is the only target, but for others there are also targets for
# input or intermediate files/results. All `all_tax_prob` targets are
# dynamically mapped over sequence batches as defined in `seqbatch`.
#
# The second part of the plan post-processes and combines the results into
# the final tibbles `asv_all_tax_prob`, `asv_tax_prob`, and `asv_unknown_prob`.
# This is also per sequence batch, in order to avoid potential memory issues
# from trying to load all `all_tax_prob` results at once; for large datasets
# there may be more than 1M candidate ASVs, and classifiers that return
# alternative assignments may return many rows per ASV.
#
# `asv_tax_prob` and `asv_unknown_prob` are more compact, and can be safely
# loaded in their entirety.

library(tarchetypes)

if (optimotu.pipeline::do_epa()) {
  epa_path <- "data/intermediate/epa"
  if (!dir.exists(epa_path)) {
    dir.create(epa_path, recursive = TRUE)
  }
}

taxonomy_plan <- c(
  if (optimotu.pipeline::do_protax()) {
    list(
      #### protax_dir ####
      # character : directory name
      #
      # the main Protax directory (often a symlink). Here to be sure that it is
      # present and has not changed
      protax_dir = tar_file(
        protax_dir,
        !!optimotu.pipeline::protax_location(),
        deployment = "main"
      )
    )
  },

  if (optimotu.pipeline::protax_aligned()) {
    list(
      #### aligned protax ####
      ##### all_tax_prob #####
      # tibble:
      #  `seq_idx` integer : index of sequence in seq_all_trim
      #  `rank` ordered factor : rank of taxonomic assignment (default:
      #    kingdom ... species)
      #  `parent_taxonomy` character : comma-separated taxonomy of parent to
      #    this taxon.
      #  `taxon` character : name of the taxon. Should never be `NA`. When no
      #    assignment was made then the row is dropped. An assignment to a
      #    novel taxon is indicated by `"unk"`.
      #  `prob` numeric : probability that the asv given by `seq_idx` belongs to
      #    `taxon` at `rank`. Should never be `NA`.
      #  `best_id` character : ID of best matching reference sequence. `NA` if
      #    `taxon` is `"unk"`.
      #  `best_dist` numeric : distance to best matching reference sequence.
      #    `NA` if `taxon` is `"unk"`.
      #  `second_id` character : ID of second best matching reference sequence.
      #    `NA` if `taxon` is `"unk"`.
      #  `second_dist` numeric : distance to second best matching reference
      #     sequence. `NA` if `taxon` is `"unk"`.
      #
      # In contrast to the unaligned case, each ASV may or may not have at least
      # one row at each rank; if no assignment at all was made at that rank,
      # then it will be missing.
      # When alternative assignments are each above the probability threshold
      # then all are included on different rows.
      all_tax_prob = tar_target(
        all_tax_prob,
        optimotu.pipeline::run_protax_animal(
          aln_seqs = seq_model_align,
          modeldir = protax_dir,
          id_is_int = TRUE,
          min_p = 0.02,
          info = TRUE,
          options = c("-m", "300"),
          ncpu = optimotu.pipeline::local_cpus()
        ) |>
          dplyr::transmute(
            seq_idx,
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
        pattern = map(seq_model_align), # per seqbatch
        resources = tar_resources(
          crew = tar_resources_crew(controller = "wide")
        )
      )
    )
  } else if (optimotu.pipeline::protax_unaligned()) {
    #### unaligned protax ####
    list(
      ##### protax_script #####
      # character: path and file name (executable)
      #
      # main protax script.  Slightly modified to accept various directories as
      # command line arguments
      protax_script = tar_file(
        protax_script,
        file.path(script_dir, "runprotax"),
        deployment = "main"
      ),
      ##### protax #####
      # character of length 24 : path and filename for all protax output files
      protax = tar_file(
        protax,
        withr::with_tempfile(
          "tempout",
          fileext = ".fasta",
          {
            protax_dir # dependency
            protax_script # dependency
            optimotu.pipeline::run_protax(
              seqs = optimotu.pipeline::fastx_gz_extract(
                infile = !!seq_all_trim,
                index = seq_index,
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
        pattern = map(seqbatch, seqbatch_hash), # per seqbatch
        iteration = "list",
        resources = tar_resources(
          crew = tar_resources_crew(controller = "wide")
        )
      ),

      ##### all_tax_prob #####
      # tibble:
      #  `seq_idx` integer : index of sequence in seq_all_trim
      #  `rank` ordered factor : rank of taxonomic assignment (phylum ...
      #    species)
      #  `parent_taxonomy` character : comma-separated taxonomy of parent to
      #    this taxon
      #  `taxon` character : name of the taxon
      #  `prob` numeric : probability that the asv in `seq_id` belongs to
      #    `taxon`
      #
      # Each ASV should have at least one row at each rank; if no assignment was
      # made at that rank, then `taxon` will be `NA`, `parent_taxon` may be
      # `NA`, and `prob` will be 0. When alternative assignments are each above
      # the probability threshold then all are included on different rows.
      all_tax_prob = tar_fst_tbl(
        all_tax_prob,
        grep("query\\d.nameprob", protax, value = TRUE) |>
          optimotu.pipeline::parse_protax_nameprob(id_is_int = TRUE),
        pattern = map(protax),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      )
    )
  } else if (optimotu.pipeline::do_sintax()) {
    list(
      #### sintax_ref_file ####
      # character: file name
      sintax_ref_file = tar_file(
        sintax_ref_file,
        !!optimotu.pipeline::sintax_ref(),
        deployment = "main"
      ),

      #### all_tax_prob ####
      # tibble:
      #  `seq_idx` integer : index of sequence in seq_all_trim
      #  `rank` ordered factor : rank of taxonomic assignment (phylum ...
      #    species)
      #  `parent_taxonomy` character : comma-separated taxonomy of parent to
      #    this taxon
      #  `taxon` character : name of the taxon
      #  `prob` numeric : probability that the asv in `seq_id` belongs to
      #    `taxon`
      all_tax_prob = tar_fst_tbl(
        all_tax_prob,
        optimotu.pipeline::sintax(
          query = optimotu.pipeline::fastx_gz_extract(
            infile = !!seq_all_trim,
            index = seq_index,
            i = seqbatch$seq_idx,
            outfile = withr::local_tempfile(fileext = ".fasta"),
            hash = seqbatch_hash
          ),
          ref = sintax_ref_file,
          ncpu = local_cpus(),
          id_is_int = TRUE
        ),
        pattern = map(seqbatch, seqbatch_hash),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "wide")
        )
      )
    )
  } else if (optimotu.pipeline::do_bayesant()) {
    list(
      if (!is.null(optimotu.pipeline::bayesant_ref())) {
        list(
          #### bayesant_ref_file ####
          # character: file name
          bayesant_ref_file = tar_file(
            bayesant_ref_file,
            !!optimotu.pipeline::bayesant_ref(),
            deployment = "main"
          )
        )
      },
      # filename is not given, store it as an R object.
      bayesant_model = if (is.null(optimotu.pipeline::bayesant_model())) {
        #### bayesant_model ####
        # object of class BayesANT
        tar_target(
          bayesant_model,
          BayesANT::read.BayesANT.data(
            fasta.file = bayesant_ref_file,
            rank = !!length(optimotu.pipeline::unknown_ranks()),
            rank_names = optimotu.pipeline::unknown_ranks()
          ) |>
            BayesANT::BayesANT(
              typeseq = !!(if (optimotu.pipeline::bayesant_aligned()) {
                "aligned"
              } else {
                " not aligned"
              })
            ),
          resources = tar_resources(
            crew = tar_resources_crew(controller = "wide") # for memory
          )
        )
      } else {
        # filename is given, store it as a file and read it in the target
        #### bayesant_model ####
        # `character`: file name
        if (file.exists(optimotu.pipeline::bayesant_model())) {
          # if the file exists, use it
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
                rank = !!length(optimotu.pipeline::unknown_ranks()),
                rank_names = optimotu.pipeline::unknown_ranks()
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
              crew = tar_resources_crew(controller = "wide") # for memory
            )
          )
        }
      },
      #### all_tax_prob ####
      # tibble:
      #  `seq_idx` integer : index of sequence in seq_all_trim
      #  `rank` ordered factor : rank of taxonomic assignment (phylum ...
      #    species)
      #  `parent_taxonomy` character : comma-separated taxonomy of parent to
      #    this taxon
      #  `taxon` character : name of the taxon
      #  `prob` numeric : probability that the asv in `seq_idx` belongs to
      #    `taxon`
      all_tax_prob = if (optimotu.pipeline::bayesant_aligned()) {
        tar_fst_tbl(
          all_tax_prob,
          optimotu.pipeline::bayesant(
            query = seq_model_align,
            model = !!optimotu.pipeline::read_bayesant_model(),
            ncpu = local_cpus(),
            id_is_int = TRUE
          ),
          pattern = map(seq_model_align),
          resources = tar_resources(
            crew = tar_resources_crew(controller = "wide")
          )
        )
      } else {
        tar_fst_tbl(
          all_tax_prob,
          optimotu.pipeline::bayesant(
            query = seq_index,
            file = seq_all_trim,
            seq_idx = seqbatch$seq_idx,
            model = !!optimotu.pipeline::read_bayesant_model(),
            ncpu = local_cpus(),
            id_is_int = TRUE,
            hash = seqbatch_hash
          ),
          pattern = map(seqbatch, seqbatch_hash),
          resources = tar_resources(
            crew = tar_resources_crew(controller = "wide")
          )
        )
      }
    )
  } else if (optimotu.pipeline::do_epa()) {
    #### epa-ng ####
    list(
      ##### epa_ref_file #####
      # character: file name
      tar_file(
        epa_ref_file,
        !!optimotu.pipeline::epa_ref(),
        deployment = "main"
      ),
      ##### epa_taxonomy_file #####
      # character: file name
      tar_file(
        epa_taxonomy_file,
        !!optimotu.pipeline::epa_taxonomy(),
        deployment = "main"
      ),
      ##### epa_tree_file #####
      # character: file name
      tar_file(
        epa_tree_file,
        !!optimotu.pipeline::epa_tree(),
        deployment = "main"
      ),
      if (file.exists(optimotu.pipeline::epa_params())) {
        ##### epa_params #####
        # character: file name
        tar_file(
          epa_params,
          !!optimotu.pipeline::epa_params(),
          deployment = "main"
        )
      } else {
        ##### epa_params #####
        # character: file name
        tar_target(
          epa_params,
          optimotu.pipeline::epa_params(),
          deployment = "main"
        )
      },
      if (file.exists(optimotu.pipeline::epa_outgroup())) {
        ##### epa_outgroup #####
        # character: file name
        tar_file(
          epa_outgroup,
          !!optimotu.pipeline::epa_outgroup(),
          deployment = "main"
        )
      } else {
        ##### epa_outgroup #####
        # character: outgroup(s)
        tar_target(
          epa_outgroup,
          optimotu.pipeline::epa_outgroup(),
          deployment = "main"
        )
      },
      ##### epa_ng #####
      # character: file name
      tar_file(
        epa_ng,
        optimotu.pipeline::epa_ng(
          ref_msa = epa_ref_file,
          tree = epa_tree_file,
          query = seq_model_align,
          outdir = file.path(epa_path, tar_name()),
          model = epa_params,
          strip_inserts = TRUE
        ),
        pattern = map(seq_model_align),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "wide")
        )
      ),

      ##### all_tax_prob #####
      # tibble:
      #  `seq_idx` integer : index of sequence in seq_all_trim
      #  `rank` ordered factor : rank of taxonomic assignment (phylum ...
      #    species)
      #  `parent_taxonomy` character : comma-separated taxonomy of parent to
      #    this taxon
      #  `taxon` character : name of the taxon
      #  `prob` numeric : probability that the asv in `seq_idx` belongs to
      #    `taxon`
      tar_fst_tbl(
        all_tax_prob,
        optimotu.pipeline::gappa_assign(
          jplace = epa_ng,
          taxonomy = epa_taxonomy_file,
          outgroup = epa_outgroup,
          ncpu = optimotu.pipeline::local_cpus(),
          id_is_int = TRUE
        ),
        pattern = map(epa_ng),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "wide")
        )
      )
    )
  } else {
    stop("No taxonomy assignment method selected")
  },
  list(
    #### asv_all_tax_prob ####
    # tibble:
    #  `seq_id` character : unique asv id
    #  `rank` ordered factor : rank of taxonomic assignment (phylum ... species)
    #  `parent_taxonomy` character : comma-separated taxonomy of parent to this
    #     taxon
    #  `taxon` character : name of the taxon
    #  `prob` numeric : probability that the asv in `seq_id` belongs to `taxon`
    #  `...` : additional columns from `all_tax_prob`, which vary by classifier
    #
    # This tibble is the foundation for the subsequent taxonomic tables.
    # It includes the "known" ranks with probability 1.0, and the remaining
    # ranks with the probability assigned by the classifier.
    asv_all_tax_prob = tar_fst_tbl(
      asv_all_tax_prob,
      tidyr::crossing(
        seq_idx = seqbatch$seq_idx,
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
          all_tax_prob,
          by = c("seq_idx", "rank", "taxon", "parent_taxonomy", "prob")
        ) |>
        dplyr::inner_join(asv_names, by = "seq_idx") |>
        dplyr::select(seq_id, everything() & !seq_idx) |>
        dplyr::arrange(seq_id, dplyr::desc(rank), dplyr::desc(prob)),
      pattern = map(seqbatch, all_tax_prob),
      resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
    ),

    #### asv_tax_prob ####
    # tibble:
    #  `seq_id` character : unique ASV id
    #  `rank` character : taxonomic rank (e.g., kingdom...species)
    #  `taxon` character : name of taxon assigned at rank
    #  `prob` numeric : probability that taxon assignment is correct
    #
    # Long format: most probable assignment and associated probability per rank.
    # `taxon` and `prob` should never be `NA`. When no assignment was made then
    # the row is dropped.
    asv_tax_prob = tar_fst_tbl(
      asv_tax_prob,
      asv_all_tax_prob |>
        # ensure known taxa are included
        dplyr::full_join(
          tidyr::crossing(
            seq_id = unique(asv_all_tax_prob$seq_id),
            tibble::tibble(
              rank = !!optimotu.pipeline::known_ranks(),
              taxon = !!optimotu.pipeline::known_taxa(),
              prob = 1.0
            )
          ),
          by = c("seq_id", "rank", "taxon", "prob")
        ) |>
        # collapse to most probable assignment per rank per seq_id
        dplyr::summarize(
          prob = max(prob, na.rm = TRUE),
          .by = c(taxon, rank, seq_id)
        ) |>
        dplyr::arrange(dplyr::desc(prob)) |>
        dplyr::summarize(
          taxon = dplyr::first(taxon),
          prob = dplyr::first(prob),
          .by = c(seq_id, rank)
        ),
      pattern = map(asv_all_tax_prob),
      resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
    ),

    #### asv_unknown_prob ####
    # tibble:
    #  `seq_id` character : unique ASV ID
    #  `rank` ordered factor : taxonomic rank (e.g., kingdom...species)
    #  `novel_prob` numeric : cumulative probability that the ASV belongs to any
    #    novel taxon at `rank`. May be `NA` for classifiers which cannot assign
    #    novelty.
    #  `known_prob` numeric : maximum probability that the ASV belongs to any
    #    one known taxon at `rank`. Should never be `NA`.
    #  `known_taxon` character : if `known_prob` is nonzero, the name of a known
    #    taxon which the ASV belongs to with probability `known_prob`. When
    #    `known_prob` < 0.5, it is possible for there to be more than one such
    #    taxon, but only one is given. `NA` if `known_prob` is 0.
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
        dplyr::filter(!rank %in% !!optimotu.pipeline::known_ranks()) |>
        dplyr::arrange(seq_id, desc(rank)) |>
        dplyr::mutate(
          novel_prob = cumsum(novel_prob),
          .by = seq_id
        ),
      pattern = map(asv_all_tax_prob),
      resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
    ),

    #### asv_tax_prob_reads ####
    # tibble:
    #  `seq_id` character : unique asv ID
    #  `rank` character : taxonomic rank (e.g., kingdom...species)
    #  `taxon` character : name of taxon assigned at rank
    #  `prob` numeric : probability that taxon assignment is correct
    #  `nread` integer : number of reads for the ASV
    asv_tax_prob_reads = tar_fst_tbl(
      asv_tax_prob_reads,
      dplyr::full_join(
        tidyr::pivot_longer(
          asv_tax,
          c(!!!optimotu.pipeline::tax_rank_vars()),
          names_to = "rank",
          values_to = "taxon"
        ),
        tidyr::pivot_longer(
          asv_tax_prob,
          c(!!!optimotu.pipeline::tax_rank_vars()),
          names_to = "rank",
          values_to = "prob"
        ),
        by = c("seq_id", "rank")
      ) |>
        dplyr::inner_join(asv_reads, by = "seq_id"),
      deployment = "main"
    )
  )
)

optimotu_plan <- c(optimotu_plan, taxonomy_plan)
