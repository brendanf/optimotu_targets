if (length(target_taxa) > 0) {
  all_tax_prob_meta <-
    if (optimotu.pipeline::do_supp_asv()) {
      tibble::tibble(
        .all_tax_prob = c(
          rlang::syms("asv_all_tax_prob"),
          optimotu.pipeline::tar_map_symbols(
            supp_tax_file_plan,
            "supp_tax_prob"
          ),
          optimotu.pipeline::tar_map_symbols(
            supp_tax_header_plan,
            "supp_tax_prob"
          ),
          optimotu.pipeline::tar_map_symbols(
            supp_tax_classify_plan,
            "supp_all_tax_prob"
          ),
        ),
        .set_id = c(
          "native",
          supp_meta_tax_file$.set_id,
          supp_meta_tax_header$.set_id,
          supp_meta_tax_classify$.set_id
        )
      )
    } else {
      tibble::tibble(
        .all_tax_prob = rlang::syms("asv_all_tax_prob"),
        .set_id = "native"
      )
    }

  target_taxa_plan <-
    tar_map(
      values = all_tax_prob_meta,
      names = .set_id,
      #### focus_hit ####
      # `tibble` with a single column:
      #  - seq_id: character vector of OTU identifiers
      #
      #  Finds all OTUs in each batch which contain an ASV with *any*
      #  probability of being one of the target taxa.
      focus_otus = tar_fst_tbl(
        focus_otus,
        .all_tax_prob |>
          dplyr::filter(taxon %in% !!target_taxa) |>
          dplyr::distinct(seq_id) |>
          dplyr::mutate(
            seq_id = dplyr::recode_values(
              seq_id,
              from = asv_otu_map$ASV,
              to = asv_otu_map$OTU,
              default = NA_character_
            )
          ) |>
          dplyr::drop_na(seq_id) |>
          dplyr::distinct(seq_id),
        pattern = map(.all_tax_prob),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      ),

      target_otus = tar_fst_tbl(
        target_otus,
        dplyr::inner_join(
          dplyr::select(.all_tax_prob, !any_of("tar_group")),
          focus_asvs,
          by = "seq_id"
        ),
        pattern = map(.all_tax_prob),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      )
    )

  target_taxa_plan <- c(
    target_taxa_plan,
    list(
      focus_asvs = tar_fst_tbl(
        focus_asvs,
        dplyr::inner_join(
          (!!optimotu.pipeline::tar_map_bind_rows(
            target_taxa_plan,
            "focus_otus"
          )) |>
            dplyr::select(!any_of("tar_group")) |>
            unique(),
          asv_otu_map,
          by = c("seq_id" = "OTU")
        ) |>
          dplyr::select(OTU = seq_id, seq_id = ASV),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      ),

      target_otus = tar_fst_tbl(
        target_otus,
        !!optimotu.pipeline::tar_map_bind_rows(
          target_taxa_plan,
          "target_otus"
        ),
        deployment = "main"
      )
    ),

    if (length(optimotu.pipeline::output_table_formats()) > 0L) {
      list(
        tar_file(
          write_target_otus,
          optimotu.pipeline::write_tabular_outputs(
            target_otus,
            file.path(
              !!optimotu.pipeline::output_path(),
              !!if (optimotu.pipeline::do_rarefy()) {
                quote(sprintf(
                  "target_taxon_otus_%s_%s",
                  .conf_level,
                  .rarefy_text
                ))
              } else {
                quote(sprintf("target_taxon_otus_%s", .conf_level))
              }
            ),
            formats = !!optimotu.pipeline::output_table_formats()
          ),
          deployment = "main"
        )
      )
    }
  )

  target_taxa_plan <- tar_map(
    values = post_cluster_meta,
    names = .conf_level,
    target_taxa_plan
  )

  optimotu_plan <- c(optimotu_plan, target_taxa_plan)
}
