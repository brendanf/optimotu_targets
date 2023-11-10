target_taxa_plan <- if (length(target_taxa) > 0) {
  tar_map(
    values = dplyr::mutate(
      reliability_meta,
      taxon_table_fungi = rlang::syms(paste0("taxon_table_fungi_", .conf_level)),
      otu_taxonomy = rlang::syms(paste0("otu_taxonomy_", .conf_level))
    ),
    names = .conf_level,
    tar_fst_tbl(
      target_otus,
      find_target_taxa(
        target_taxa,
        asv_all_tax_prob,
        taxon_table_fungi,
        otu_taxonomy
      )
    ),
    tar_file_fast(
      write_target_otus,
      write_and_return_file(
        target_otus,
        sprintf("output/target_taxon_otus_%s.rds", .conf_level),
        type = "rds"
      )
    ),
    tar_file_fast(
      target_otu_seqs,
      vapply(
        target_taxa,
        \(taxon) dplyr::filter(target_otus_reliable, protax_taxon == taxon) |>
          dplyr::left_join(asv_seq, by = c("asv_seq_id" = "seq_id")) |>
          glue::glue_data(">{seq_id}_{asv_seq_id};p={protax_prob};species={taxon}\n{seq}") |>
          write_and_return_file(sprintf("output/%s_%s.fasta", taxon, .conf_level)),
        ""
      )
    ),
    tar_file_fast(
      target_otu_matches,
      vapply(
        target_otu_seqs,
        \(seqfile) {
          vsearch_usearch_global_dbmatched(
            query = seqfile,
            ref = file.path(protax_model, "sintaxits2train.fa"),
            output = sub(".fasta", "_matches.fasta", seqfile, fixed = TRUE),
            threshold = 0.9,
            global = FALSE
          )
        },
        ""
      )
    )
  )
} else {
  list()
}
