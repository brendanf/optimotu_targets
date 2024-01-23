target_taxa_plan <- if (length(target_taxa) > 0) {
  tar_map(
    values = post_cluster_meta,
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
    )
  )
} else {
  list()
}
