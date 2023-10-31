guild_plan <- list(
  #### funguild_db ####
  tar_fst_tbl(
    funguild_db,
    FUNGuildR::get_funguild_db()
  ),

  #### map over confidence levels ####
  tar_map(
    # also map over some previously mapped targets
    values = tibble::tibble(
      .conf_level = c("plausible", "reliable"),
      otu_abund_table_sparse = paste0("otu_abund_table_sparse_", .conf_level) %>%
        rlang::syms(),
      otu_taxonomy = paste0("otu_taxonomy_", .conf_level) %>%
        rlang::syms(),
      asv_otu_map = paste0("asv_otu_map_", .conf_level) %>%
        rlang::syms()
    ),
    names = .conf_level,

    ##### otu_guild_{.conf_level} #####
    tar_fst_tbl(
      otu_guild,
      otu_taxonomy |>
        dplyr::mutate(
          dplyr::across(
            genus:species,
            sub,
            pattern = "([A-Z].+)_[0-9]+",
            replacement = "\\1"
          )
        ) |>
        tidyr::unite("Taxonomy", kingdom:species, sep = ",") |>
        FUNGuildR::funguild_assign(db = funguild_db) |>
        dplyr::select(OTU, guild)
    ),
    ##### write_otu_guild_{.conf_level} #####
    tar_file(
      write_otu_guild,
      write_and_return_file(
        otu_guild,
        sprintf("output/otu_guilds_%s.tsv", .conf_level),
        type = "tsv"
      )
    )
  )
)
