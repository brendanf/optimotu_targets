if (length(target_taxa) > 0) {
  target_taxa_plan <- tar_map(
    values = post_cluster_meta,
    names = .conf_level,
    tar_fst_tbl(
      target_otus,
      optimotu.pipeline::find_target_taxa(
        target_taxa,
        asv_all_tax_prob,
        taxon_table_ingroup,
        otu_taxonomy
      ),
      deployment = "main"
    ),
    tar_file(
      write_target_otus,
      optimotu.pipeline::write_and_return_file(
        target_otus,
        sprintf("%s/target_taxon_otus_%s.rds", !!optimotu.pipeline::output_path(), .conf_level),
        type = "rds"
      ),
      deployment = "main"
    )
  )
  optimotu_plan <- c(optimotu_plan, target_taxa_plan)
}
