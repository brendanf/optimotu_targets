if (optimotu.pipeline::do_guilds()) {
  guild_dbs <- optimotu.pipeline::guild_databases()
  download_dbs <- dplyr::filter(guild_dbs, source == "download")
  file_dbs <- dplyr::filter(guild_dbs, source != "download")

  guild_plan <- list()
  if (nrow(download_dbs) > 0L) {
    guild_plan <- c(
      guild_plan,
      tar_map(
        values = dplyr::transmute(download_dbs, .guild = name),
        names = .guild,
        #### guild_db_{.guild} ####
        guild_db = tar_fst_tbl(
          guild_db,
          optimotu.pipeline::load_guild_database("download"),
          deployment = "main"
        )
      )
    )
  }
  if (nrow(file_dbs) > 0L) {
    guild_plan <- c(
      guild_plan,
      tar_map(
        values = dplyr::transmute(
          file_dbs,
          .guild = name,
          .source = source,
          .path = path
        ),
        names = .guild,
        #### guild_db_file_{.guild} ####
        guild_db_file = tar_file(
          guild_db_file,
          .path,
          deployment = "main"
        ),
        #### guild_db_{.guild} ####
        guild_db = tar_fst_tbl(
          guild_db,
          optimotu.pipeline::load_guild_database(
            source = .source,
            path = guild_db_file
          ),
          deployment = "main"
        )
      )
    )
  }

  guild_plan <- c(
    guild_plan,
    #### map over confidence levels ####
    tar_map(
      # also map over some previously mapped targets
      values = tibble::tibble(
        .conf_level = c("plausible", "reliable"),
        otu_abund_table_sparse = paste0(
          "otu_abund_table_sparse_",
          .conf_level
        ) |>
          rlang::syms(),
        otu_taxonomy = paste0("otu_taxonomy_", .conf_level) |>
          rlang::syms(),
        asv_otu_map = paste0("asv_otu_map_", .conf_level) |>
          rlang::syms()
      ),
      names = .conf_level,

      tar_map(
        values = tibble::tibble(
          .guild_db = rlang::syms(paste0("guild_db_", guild_dbs$name)),
          .guild = guild_dbs$name
        ),
        names = .guild,

        ###### otu_guild_{.guild}_{.conf_level} ######
        tar_fst_tbl(
          otu_guild,
          optimotu.pipeline::prepare_guild_taxonomy(
            otu_taxonomy,
            ranks = !!optimotu.pipeline::tax_ranks()
          ) |>
            FUNGuildR::funguild_assign(db = .guild_db) |>
            dplyr::select(seq_id, guild),
          deployment = "main"
        ),
        ###### write_otu_guild_{.guild}_{.conf_level} ######
        if (length(optimotu.pipeline::output_table_formats()) > 0L) {
          tar_file(
            write_otu_guild,
            optimotu.pipeline::write_tabular_outputs(
              otu_guild,
              file.path(
                !!optimotu.pipeline::output_path(),
                !!(if (optimotu.pipeline::do_rarefy()) {
                  quote(sprintf(
                    "otu_guilds_%s_%s_%s",
                    .guild,
                    .conf_level,
                    .rarefy_text
                  ))
                } else {
                  quote(sprintf("otu_guilds_%s_%s", .guild, .conf_level))
                })
              ),
              formats = !!optimotu.pipeline::output_table_formats()
            ),
            deployment = "main"
          )
        }
      )
    )
  )

  optimotu_plan <- c(optimotu_plan, guild_plan)
}
