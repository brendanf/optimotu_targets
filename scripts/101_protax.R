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
  
  tar_file(
    protax_spikelist,
    protax[basename(protax) == "spikeout"] %>%
      lapply(readLines) %>%
      unlist() %>%
      write_and_return_file(file.path(protax_path, "spikelist")),
    deployment = "main"
  ),
  
  tar_file(
    protax_asv_protax_levels,
    vapply(
      2:7,
      function(i) {
        protax[basename(protax) == sprintf("query%d.nameprob")] %>%
          lapply(readLines) %>%
          unlist() %>%
          write_and_return_file(
            filepath(protax_path, "ASVprotaxLevel%d.txt", i)
          )
      },
      ""
    ),
    deployment = "main"
  ),
  
  tar_file(
    protax_tables,
    {
      output_root <- file.path(protax_path, "asvprotax")
      result = system2(
        "perl",
        c(
          "protaxFungi/scripts/make_taxtable.pl",
          protax_spikelist,
          Biobase::lcPrefix(protax_asv_protax_levels),
          output_root
        )
      )
      stopifnot(result == 0)
      paste0(output_root, c("_names.txt", "_prob.txt"))
    },
    deployment = "main"
  ),
  
  tar_fst_tbl(
    asv_tax,
    readr::read_tsv(protax_tables[1], col_type = "c") %>%
      dplyr::filter(kingdom != "Spike"),
    deployment = "main"
  ),
  
  #### asv_tax_prob ####
  # read probabilities of taxonomic assignments
  tar_fst_tbl(
    asv_tax_prob,
    readr::read_tsv(protax_tables[2], col_types = "cddddddd"),
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
        sample = sub( "CCDB-\\d{5}_", "", basename(sample)),
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
    purrr::reduce(
      list(
        tidyr::pivot_longer(asv_tax, kingdom:species, names_to = "rank", values_to = "taxon"),
        tidyr::pivot_longer(asv_tax_prob, kingdom:species, names_to = "rank", values_to = "prob"),
        asv_reads
      ),
      dplyr::left_join
    ),
    deployment = "main"
  )
  
  
)
