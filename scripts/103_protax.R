library(tarchetypes)

protax_plan <- list(
  tar_file(
    protax_dir,
    "protaxFungi",
    deployment = "main"
  ),
  tar_file(
    protax_script,
    "scripts/runprotax"
  ),
  
  tar_file(
    protax,
    {
      protax_dir
      protax_script
      run_protax(
        seqs = primer_trim,
        outdir = file.path(protax_path, tar_name()),
        modeldir = protax_model
      )
    },
    pattern = map(primer_trim),
    iteration = "list"
  ),
  
  tar_fst_tbl(
    protax_spikelist,
    protax[basename(protax) == "spikeout"] %>%
      lapply(
        readr::read_tsv,
        col_names = c("seq_id", "size", "spike", "match"),
        col_types = "cicd"
      ) %>%
      dplyr::bind_rows(),
    deployment = "main"
  ),
  
  tar_fst_tbl(
    asv_all_tax_prob,
    lapply(protax, grep, pattern = "query\\d.nameprob", value = TRUE) |>
      lapply(parse_protax_nameprob) |>
      purrr::map2_dfr(seqbatch_key, dplyr::left_join, by = "seq_id") |>
      dplyr::select(-seq_id, seq_id = i) |>
      name_seqs("ASV"),
    pattern = map(protax, seqbatch_key)
  ),
  
  tar_fst_tbl(
    asv_tax,
    asv_all_tax_prob %>%
      dplyr::anti_join(protax_spikelist, by = "seq_id") %>%
      dplyr::group_by(rank, seq_id) %>%
      dplyr::summarize(taxon = dplyr::first(taxon), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = rank, values_from = taxon) %>%
      dplyr::mutate(kingdom = "Fungi") %>%
      dplyr::select("seq_id", "kingdom", "phylum", "class", "order", "family", "genus", "species"),
    deployment = "main"
  ),
  
  #### asv_tax_prob ####
  # read probabilities of taxonomic assignments
  tar_fst_tbl(
    asv_tax_prob,
    asv_all_tax_prob %>%
      dplyr::anti_join(protax_spikelist, by = "seq_id") %>%
      dplyr::group_by(rank, seq_id) %>%
      dplyr::summarize(prob = dplyr::first(prob), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = rank, values_from = prob) %>%
      dplyr::mutate(kingdom = 1) %>%
      dplyr::select("seq_id", "kingdom", "phylum", "class", "order", "family", "genus", "species"),
    deployment = "main"
  ),
  
  #### asv_tax_seq ####
  # combine taxonomy and sequence
  # don't include subsequence duplicates
  tar_fst_tbl(
    asv_tax_seq,
    dplyr::left_join(asv_tax, asv_seq, by = "seq_id"),
    deployment = "main"
  ),
  
  #### asv_tax_prob_reads ####
  tar_fst_tbl(
    asv_tax_prob_reads,
    dplyr::full_join(
      tidyr::pivot_longer(asv_tax, kingdom:species, names_to = "rank", values_to = "taxon"),
      tidyr::pivot_longer(asv_tax_prob, kingdom:species, names_to = "rank", values_to = "prob"),
      by = c("seq_id", "rank")
    ) %>%
      dplyr::inner_join(asv_reads),
    deployment = "main"
  )
)
