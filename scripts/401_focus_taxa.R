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
    if (length(optimotu.pipeline::output_table_formats()) > 0L) {
      tar_file(
        write_target_otus,
        optimotu.pipeline::write_tabular_outputs(
          target_otus,
          file.path(
            !!optimotu.pipeline::output_path(),
            sprintf("target_taxon_otus_%s", .conf_level)
          ),
          formats = !!optimotu.pipeline::output_table_formats()
        ),
        deployment = "main"
      )
    }
  )
  optimotu_plan <- c(optimotu_plan, target_taxa_plan)
}
