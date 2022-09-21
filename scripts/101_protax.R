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
  
  tar_fst_tbl(
    asv_seq,
    tibble::tibble(
      ASV = sprintf(
        sprintf("ASV%%0%dd", ceiling(log10(ncol(asvtable)))),
        seq_len(ncol(asvtable))
      ),
      seq = colnames(asvtable)
    ),
    deployment = "main"
  ),
  
  tar_group_count(
    grouped_asv_seq,
    asv_seq,
    count = 12,
    deployment = "main"
  ),
  
  tar_file(
    protax,
    {
      protax_dir
      protax_script
      run_protax(grouped_asv_seq, file.path(protax_path, tar_name()))
    },
    pattern = map(grouped_asv_seq)
  ),
  
  tar_fst_tbl(
    protax_spikelist,
    protax[basename(protax) == "spikeout"] %>%
      lapply(
        readr::read_tsv,
        col_names = c("ASV", "size", "spike", "match"),
        col_types = "cicd"
      ) %>%
      dplyr::bind_rows(),
    deployment = "main"
  ),
  
  tar_fst_tbl(
    asv_all_tax_prob,
    protax[grepl("query\\d.nameprob", basename(protax))] %>%
      set_names(., basename(.)) %>%
      lapply(readLines) %>%
      tibble::enframe() %>%
      tidyr::extract(
        name,
        into = "rank",
        regex = "query(\\d+)\\.nameprob",
        convert = TRUE
      ) %>%
      tidyr::unchop(value) %>%
      dplyr::mutate(
        rank = rank2factor(TAXRANKS[rank]),
        value = gsub("([^\t]+)\t([0-9.]+)", "\\1:\\2", value) %>%
          gsub("(:[0-9.]+)\t", "\\1;", .)
      ) %>%
      tidyr::separate(value, into = c("ASV", "nameprob"), sep = "\t") %>%
      tidyr::separate_rows(nameprob, sep = ";") %>%
      tidyr::separate(nameprob, into = c("name", "prob"), sep = ":", convert = TRUE) %>%
      tidyr::extract(name, into = c("parent_taxonomy", "taxon"), regex = "(.+),([^,]+)$") %>%
      dplyr::mutate(
        taxon = dplyr::na_if(taxon, "unk"),
        prob = ifelse(is.na(taxon), 0, prob)
      ) %>%
      dplyr::arrange(ASV, rank, dplyr::desc(prob)),
    deployment = "main"
  ),
  
  tar_fst_tbl(
    asv_tax,
    asv_all_tax_prob %>%
      dplyr::anti_join(protax_spikelist, by = "ASV") %>%
      dplyr::group_by(rank, ASV) %>%
      dplyr::summarize(taxon = dplyr::first(taxon), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = rank, values_from = taxon) %>%
      dplyr::mutate(kingdom = "Fungi") %>%
      dplyr::select("ASV", "kingdom", "phylum", "class", "order", "family", "genus", "species"),
    deployment = "main"
  ),
  
  #### asv_tax_prob ####
  # read probabilities of taxonomic assignments
  tar_fst_tbl(
    asv_tax_prob,
    asv_all_tax_prob %>%
      dplyr::anti_join(protax_spikelist, by = "ASV") %>%
      dplyr::group_by(rank, ASV) %>%
      dplyr::summarize(prob = dplyr::first(prob), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = rank, values_from = prob) %>%
      dplyr::mutate(kingdom = 1) %>%
      dplyr::select("ASV", "kingdom", "phylum", "class", "order", "family", "genus", "species"),
    deployment = "main"
  ),
  
  #### asv_tax_seq ####
  # combine taxonomy and sequence
  # don't include subsequence duplicates
  tar_fst_tbl(
    asv_tax_seq,
    dplyr::left_join(asv_tax, asv_seq, by = "ASV"),
    deployment = "main"
  ),
  
  tar_fst_tbl(
    asv_table,
    asvtable %>%
      dplyr::na_if(0L) %>%
      tibble::as_tibble(rownames = "sample") %>%
      tidyr::pivot_longer(-1, names_to = "seq", values_to = "nread", values_drop_na = TRUE) %>%
      dplyr::left_join(asv_seq) %>%
      dplyr::transmute(
        ASV = ASV,
        sample = sub("CCDB-\\d{3}_", "", basename(sample)),
        nread = nread
      ) %>%
      dplyr::arrange(ASV, sample),
    deployment = "main"
  ),
  
  #### asv_reads ####
  tar_fst_tbl(
    asv_reads,
    asv_table %>%
      dplyr::group_by(ASV) %>%
      dplyr::summarize(nread = sum(nread)) %>%
      dplyr::semi_join(asv_tax, by = "ASV"),
    deployment = "main"
  ),
  
  #### asv_tax_prob_reads ####
  tar_fst_tbl(
    asv_tax_prob_reads,
    dplyr::full_join(
      tidyr::pivot_longer(asv_tax, kingdom:species, names_to = "rank", values_to = "taxon"),
      tidyr::pivot_longer(asv_tax_prob, kingdom:species, names_to = "rank", values_to = "prob"),
      by = c("ASV", "rank")
    ) %>%
      dplyr::inner_join(asv_reads),
    deployment = "main"
  )
)
