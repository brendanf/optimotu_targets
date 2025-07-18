library(tarchetypes)

if (optimotu.pipeline::do_epa()) {
  epa_path <- "data/intermediate/epa"
  if (!dir.exists(epa_path)) {
    dir.create(epa_path, recursive = TRUE)
  }
}

taxonomy_plan <- list(

  if (optimotu.pipeline::do_protax()) {
    #### protax_dir ####
    # character : directory name
    #
    # the main Protax directory (often a symlink). Here to be sure that it is
    # present and has not changed
    tar_file(
      protax_dir,
      !!optimotu.pipeline::protax_location(),
      deployment = "main"
    )
  },

  if (optimotu.pipeline::protax_aligned()) {
    #### aligned protax ####
    ##### all_tax_prob #####
    # tibble:
    #  `seq_idx` integer : index of sequence in seq_all_trim
    #  `rank` ordered factor : rank of taxonomic assignment (phylum ... species)
    #  `parent_taxonomy` character : comma-separated taxonomy of parent to this taxon
    #  `taxon` character : name of the taxon
    #  `prob` numeric : probability that the asv in `seq_id` belongs to `taxon`
    #
    # In contrast to the unaligned case, each ASV may or may not have at least
    # one row at each rank; if no assignment at all was made at that rank,
    # then it will be missing.  If `taxon` is `NA`, this indicates an actual
    # prediction of "unknown taxon at this rank", and has an associated
    # `parent_taxon` and `prob`.
    # When alternative assignments are each above the probability threshold (10%)
    # then all are included on different rows.
    tar_target(
      all_tax_prob,
      withr::with_tempfile(
        "td",
        optimotu.pipeline::fastx_split(
          asv_model_align,
          n = optimotu.pipeline::local_cpus(),
          outroot = optimotu.pipeline::ensure_directory(tempfile(tmpdir = td))
        ) |>
          optimotu.pipeline::run_protax_animal(
            modeldir = protax_dir,
            id_is_int = TRUE,
            min_p = 0.02,
            info = TRUE,
            options = c("-m", "300")
          ) |>
          dplyr::transmute(
            seq_idx,
            rank = optimotu.pipeline::int2rankfactor(rank),
            parent_taxonomy = paste(
              paste(optimotu.pipeline::known_taxa(), collapse = ","),
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
      ),
      pattern = map(asv_model_align), # per seqbatch
      resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
    )
  } else if (optimotu.pipeline::protax_unaligned()) {
    #### unaligned protax ####
    list(
      ##### protax_script #####
      # character: path and file name (executable)
      #
      # main protax script.  Slightly modified to accept various directories as
      # command line arguments
      tar_file(
        protax_script,
        file.path(script_dir, "runprotax"),
        deployment = "main"
      ),
      ##### protax #####
      # character of length 24 : path and filename for all protax output files
      tar_file(
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
              outdir = file.path(!!optimotu.pipeline::protax_path(), tar_name()),
              modeldir = protax_model,
              script = protax_script
            )
          }
        ),
        pattern = map(seqbatch, seqbatch_hash), # per seqbatch
        iteration = "list",
        resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
      ),

      ##### all_tax_prob #####
      # tibble:
      #  `seq_idx` integer : index of sequence in seq_all_trim
      #  `rank` ordered factor : rank of taxonomic assignment (phylum ... species)
      #  `parent_taxonomy` character : comma-separated taxonomy of parent to this taxon
      #  `taxon` character : name of the taxon
      #  `prob` numeric : probability that the asv in `seq_id` belongs to `taxon`
      #
      # Each ASV should have at least one row at each rank; if no assignment was
      # made at that rank, then `taxon` will be `NA`, `parent_taxon` may be `NA`,
      # and `prob` will be 0.
      # When alternative assignments are each above the probability threshold (10%)
      # then all are included on different rows.
      tar_fst_tbl(
        all_tax_prob,
        grep("query\\d.nameprob", protax, value = TRUE) |>
          optimotu.pipeline::parse_protax_nameprob(id_is_int = TRUE),
        pattern = map(protax),
        resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
      )
    )
  } else if (optimotu.pipeline::do_sintax()) {
    list(
      #### sintax_ref_file ####
      # character: file name
      tar_file(
        sintax_ref_file,
        !!optimotu.pipeline::sintax_ref(),
        deployment = "main"
      ),

      #### all_tax_prob ####
      # tibble:
      #  `seq_idx` integer : index of sequence in seq_all_trim
      #  `rank` ordered factor : rank of taxonomic assignment (phylum ... species)
      #  `parent_taxonomy` character : comma-separated taxonomy of parent to this taxon
      #  `taxon` character : name of the taxon
      #  `prob` numeric : probability that the asv in `seq_id` belongs to `taxon`
      #
      # Each ASV should have at least one row at each rank; if no assignment was
      # made at that rank, then `taxon` will be `NA`, `parent_taxon` may be `NA`,
      # and `prob` will be 0.
      # When alternative assignments are each above the probability threshold (10%)
      # then all are included on different rows.
      tar_fst_tbl(
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
        resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
      )
    )
  } else if (optimotu.pipeline::do_bayesant()) {
    list(
      if (train_bayesant) {
        list(
          #### bayesant_ref_file ####
          # character: file name
          tar_file(
            bayesant_ref_file,
            !!optimotu.pipeline::bayesant_ref(),
            deployment = "main"
          ),
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
                typeseq = !!(
                  if (optimotu.pipeline::bayesant_aligned()) "aligned" else " not aligned"
                )
              ),
            resources = tar_resources(
              crew = tar_resources_crew(controller = "wide") # for memory
            )
          )
        )
      } else {
        #### bayesant_model ####
        # `character`: file name
        tar_file(
          bayesant_model,
          !!optimotu.pipeline::bayesant_model(),
          deployment = "main"
        )
      },
      #### all_tax_prob ####
      # tibble:
      #  `seq_idx` integer : index of sequence in seq_all_trim
      #  `rank` ordered factor : rank of taxonomic assignment (phylum ... species)
      #  `parent_taxonomy` character : comma-separated taxonomy of parent to this taxon
      #  `taxon` character : name of the taxon
      #  `prob` numeric : probability that the asv in `seq_idx` belongs to `taxon`
      tar_fst_tbl(
        all_tax_prob,
        optimotu.pipeline::bayesant(
          query = optimotu.pipeline::fastx_gz_extract(
            infile = !!seq_all_trim,
            index = seq_index,
            i = seqbatch$seq_idx,
            outfile = withr::local_tempfile(fileext = ".fasta"),
            hash = seqbatch_hash
          ),
          model = bayesant_model,
          ncpu = local_cpus(),
          id_is_int = TRUE
        ),
        pattern = map(seqbatch, seqbatch_hash),
        resources = tar_resources(crew = tar_resources_crew(controller = "wide") )
      )
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
          query = asv_model_align,
          outdir = file.path(epa_path, tar_name()),
          model = epa_params,
          strip_inserts = TRUE
        ),
        pattern = map(asv_model_align),
        resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
      ),

      ##### all_tax_prob #####
      # tibble:
      #  `seq_idx` integer : index of sequence in seq_all_trim
      #  `rank` ordered factor : rank of taxonomic assignment (phylum ... species)
      #  `parent_taxonomy` character : comma-separated taxonomy of parent to this taxon
      #  `taxon` character : name of the taxon
      #  `prob` numeric : probability that the asv in `seq_idx` belongs to `taxon`
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
        resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
      )

    )
  } else {
    stop("No taxonomy assignment method selected")
  },

  #### asv_all_tax_prob ####
  # tibble:
  #  `seq_id` character : unique asv id
  #  `rank` ordered factor : rank of taxonomic assignment (phylum ... species)
  #  `parent_taxonomy` character : comma-separated taxonomy of parent to this taxon
  #  `taxon` character : name of the taxon
  #  `prob` numeric : probability that the asv in `seq_id` belongs to `taxon`
  tar_fst_tbl(
    asv_all_tax_prob,
    all_tax_prob |>
      dplyr::inner_join(asv_names, by = "seq_idx") |>
      dplyr::select(seq_id, everything() & !seq_idx),
    pattern = map(all_tax_prob),
    resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
  ),

  #### asv_tax ####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  {ROOT_RANK} character : taxonomic assignment at ROOT_RANK (e.g. kingdom)
  #  ... character: taxonomic assignments at intermediate ranks
  #  {TIP_RANK} character : taxonomic assignment at TIP_RANK (e.g. species)
  #
  # The most probable assignment for each ASV at each rank.  NA if there was no
  # assignment above Protax's reporting threshold
  tar_fst_tbl(
    asv_tax,
    asv_all_tax_prob |>
      dplyr::summarize(taxon = dplyr::first(taxon), .by = c(rank, seq_id)) |>
      tidyr::pivot_wider(names_from = rank, values_from = taxon, names_expand = TRUE) |>
      purrr::reduce2(
        !!optimotu.pipeline::known_ranks(),
        !!optimotu.pipeline::known_taxa(),
        .init = _,
        .f = \(d, rank, taxon) {d[[rank]] <- taxon; d}
      ) |>
      dplyr::select("seq_id", !!!optimotu.pipeline::tax_rank_vars()),
    pattern = map(asv_all_tax_prob),
    resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
  ),

  #### asv_tax_prob ####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  {ROOT_RANK} numeric : probability for taxonomic assignment at ROOT_RANK
  #    (e.g., kingdom)
  #  ... numeric : probability for taxonomic assignments at intermediate ranks
  #  {TIP_RANK} numeric : probability for taxonomic assignment at TIP_RANK (e.g.,
  #    species)
  #
  # Associated probaility for the most probable assignment for each ASV at each
  # rank.  0 if there was no assignment above Protax's reporting threshold
  tar_fst_tbl(
    asv_tax_prob,
    asv_all_tax_prob |>
      dplyr::summarize(prob = dplyr::first(prob), .by = c(rank, seq_id)) |>
      tidyr::pivot_wider(
        names_from = rank,
        values_from = prob,
        names_expand = TRUE
      ) |>
      purrr::reduce2(
        optimotu.pipeline::known_ranks(),
        optimotu.pipeline::known_taxa(),
        .init = _,
        .f = \(d, rank, taxon) {d[[rank]] = 1.0; d}
      ) |>
      dplyr::select("seq_id", !!!optimotu.pipeline::tax_ranks()),
    pattern = map(asv_all_tax_prob),
    resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
  ),

  #### asv_unknown_prob ####
  tar_fst_tbl(
    asv_unknown_prob,
    asv_all_tax_prob |>
      dplyr::summarize(
        novel_prob = sum(prob[taxon == "unk"]),
        known_prob = max(prob[taxon != "unk"], 0),
        .by = c(seq_id, rank)
      ) |>
      tidyr::complete(seq_id, rank, fill = list(novel_prob = 0, known_prob = 0)) |>
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
  tar_fst_tbl(
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

optimotu_plan <- c(optimotu_plan, taxonomy_plan)
