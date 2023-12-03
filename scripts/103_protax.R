library(tarchetypes)

protax_plan <- list(

  #### protax_dir ####
  # character : directory name
  #
  # the main Protax directory (often a symlink). Here to be sure that it is
  # present and has not changed
  tar_file_fast(
    protax_dir,
    pipeline_options$protax_dir,
    deployment = "main"
  ),

  #### protax_script ####
  # character: path and file name (executable)
  #
  # main protax script.  Slightly modified to accept various directories as
  # command line arguments
  tar_file_fast(
    protax_script,
    file_path(protax_dir, "classify")
  ),

  #### protax ####
  # character of length 24 : path and filename for all protax output files
  tar_target(
    protax,
    fastx_split(hmm_align, n = local_cpus()) |>
      run_protax_animal(modeldir = protax_dir, strip_inserts = TRUE),
    pattern = map(hmm_align) # per seqbatch
  ),

  #### asv_all_tax_prob ####
  # tibble:
  #  `seq_idx` character : unique asv id
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
    asv_all_tax_prob,
    # lapply(protax, grep, pattern = "query\\d.nameprob", value = TRUE) |>
    #   lapply(parse_protax_nameprob) |>
    #   purrr::map2_dfr(
    #     dplyr::group_split(seqbatch_key, tar_group, .keep = FALSE),
    #     dplyr::left_join,
    #     by = "seq_id",
    #   ) |>
    #   dplyr::select(-seq_id) %>%
    #   dplyr::left_join(
    #     .[,"i"] |>
    #       unique() |>
    #       dplyr::arrange(i) |>
    #       name_seqs(prefix = "ASV", id_col = "seq_id"),
    #     by = "i"
    #   ) |>
    #   dplyr::select(seq_id, everything() & !i)
    parse_protaxAnimal_output(protax),
    pattern = map(protax)
  ),

  #### asv_tax ####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  `kingdom` character : taxonomic kingdom assignment
  #  `phylum` character : taxonomic phylum assignment
  #  `class` character : taxonomic class assignment
  #  `order` character : taxonomic order assignment
  #  `family` character : taxonomic family assignment
  #  `genus` character : taxonomic genus assignment
  #  `species` character : taxonomic species assignment
  #
  # The most probable assignment for each ASV at each rank.  NA if there was no
  # assignment above Protax's reporting threshold
  tar_fst_tbl(
    asv_tax,
    asv_all_tax_prob %>%
      dplyr::group_by(rank, seq_id) %>%
      dplyr::summarize(taxon = dplyr::first(taxon), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = rank, values_from = taxon) %>%
      dplyr::mutate(kingdom = "Fungi") %>%
      dplyr::select("seq_id", "kingdom", "phylum", "class", "order", "family", "genus", "species"),
    deployment = "main"
  ),

  #### asv_tax_prob ####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  `kingdom` numeric : probability for taxonomic kingdom assignment
  #  `phylum` numeric : probability for taxonomic phylum assignment
  #  `class` numeric : probability for taxonomic class assignment
  #  `order` numeric : probability for taxonomic order assignment
  #  `family` numeric : probability for taxonomic family assignment
  #  `genus` numeric : probability for taxonomic genus assignment
  #  `species` numeric : probability for taxonomic species assignment
  #
  # Associated probaility for the most probable assignment for each ASV at each
  # rank.  0 if there was no assignment above Protax's reporting threshold
  tar_fst_tbl(
    asv_tax_prob,
    asv_all_tax_prob %>%
      dplyr::group_by(rank, seq_id) %>%
      dplyr::summarize(prob = dplyr::first(prob), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = rank, values_from = prob) %>%
      dplyr::mutate(kingdom = 1) %>%
      dplyr::select("seq_id", "kingdom", "phylum", "class", "order", "family", "genus", "species"),
    deployment = "main"
  ),

  #### reference_tree ####
  # file name
  # reference tree for phylogenetic placement
  tar_file_fast(
    reference_tree,
    "data/reftree/split_seqs_genus_train_iqtree.treefile",
    deployment = "main"
  ),

  #### reference_tree_taxa ####
  # file name
  # taxonomy for the tips of the reference tree
  tar_file_fast(
    reference_tree_taxa,
    "data/reftree/split_seqs_genus_taxonomy.txt",
    deployment = "main"
  ),

  #### reference_tree_outgroup ####
  # file name
  # names of outgroup tips in the reference tree
  tar_file_fast(
    reference_tree_outgroup,
    "data/reftree/split_seqs_genus_outgroup.txt",
    deployment = "main"
  ),

  #### reference_tree_alignment ####
  # file name
  # alignment for the reference tree
  tar_file_fast(
    reference_tree_alignment,
    "data/reftree/split_seqs_genus_train.fasta",
    deployment = "main"
  ),

  #### reference_tree_model ####
  # file name
  # iqtree log file for the reference tree
  tar_file_fast(
    reference_tree_log,
    "data/reftree/split_seqs_genus_train_iqtree.log",
    deployment = "main"
  ),

  #### epa_gappa ####
  # `tibble`:
  #  `seq_idx` character : unique asv id
  #  `rank` ordered factor : rank of taxonomic assignment (phylum ... species)
  #  `parent_taxonomy` character : comma-separated taxonomy of parent to this taxon
  #  `taxon` character : name of the taxon
  #  `prob` numeric : probability that the asv in `seq_id` belongs to `taxon`
  tar_fst_tbl(
    epa_gappa,
    epa_ng(
      ref_msa = reference_tree_alignment,
      tree = reference_tree,
      query = hmm_align,
      model = parse_iqtree_model(reference_tree_log),
      ncpu = 1,
      strip_inserts = TRUE
    ) |>
      gappa_assign(
        taxonomy = reference_tree_taxa,
        outgroup = reference_tree_outgroup,
        ranks = c("class", "order", "family", "subfamily", "tribe", "genus"),
        ncpu = local_cpus(),
        allow_file_overwriting = TRUE,
        id_is_int = TRUE,
      ),
    pattern = map(hmm_align)
  ),

  #### asv_tax_seq ####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  `kingdom` character : taxonomic kingdom assignment
  #  `phylum` character : taxonomic phylum assignment
  #  `class` character : taxonomic class assignment
  #  `order` character : taxonomic order assignment
  #  `family` character : taxonomic family assignment
  #  `genus` character : taxonomic genus assignment
  #  `species` character : taxonomic species assignment
  #  `seq` character : sequence
  #
  # combine taxonomy and sequence
  # don't include subsequence duplicates
  tar_fst_tbl(
    asv_tax_seq,
    dplyr::left_join(asv_tax, asv_seq, by = "seq_id"),
    deployment = "main"
  ),

  #### asv_unknown_prob ####
  tar_fst_tbl(
    asv_unknown_prob,
    asv_all_tax_prob %>%
      dplyr::filter(!is.na(taxon)) %>%
      dplyr::group_by(seq_id, rank) %>%
      dplyr::summarise(prob_unk = 1-sum(prob))
  ),

  #### asv_tax_prob_reads ####
  # tibble:
  #  `seq_id` character : unique asv ID
  #  `rank` character : taxonomic rank (kingdom...species)
  #  `taxon` character : name of taxon assigned at rank
  #  `prob` numeric : probability that taxon assignment is correct
  #  `nread` integer : number of reads for the ASV
  tar_fst_tbl(
    asv_tax_prob_reads,
    dplyr::full_join(
      tidyr::pivot_longer(asv_tax, kingdom:species, names_to = "rank", values_to = "taxon"),
      tidyr::pivot_longer(asv_tax_prob, kingdom:species, names_to = "rank", values_to = "prob"),
      by = c("seq_id", "rank")
    ) %>%
      dplyr::inner_join(asv_reads, by = "seq_id"),
    deployment = "main"
  )
)
