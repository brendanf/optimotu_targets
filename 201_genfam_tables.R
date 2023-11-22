genfam_plan <- list(
  #### map over confidence levels ####
  tar_map(
    # also map over some previously mapped targets
    values = tibble::tibble(
      .conf_level = c("plausible", "reliable"),
      otu_taxonomy = paste0("otu_taxonomy_", .conf_level) %>%
        rlang::syms(),
      otu_guild = paste0("otu_guild_", .conf_level) %>%
        rlang::syms(),
      otu_table_sparse = paste0("otu_table_sparse_", .conf_level) %>%
        rlang::syms()
    ),
    names = .conf_level,
    tar_fst_tbl(
      nonlichen_otu_table,
      dplyr::filter(otu_table_sparse, nread > 0) |>
        dplyr::anti_join(
          dplyr::filter(
            otu_guild,
            guild %in% c("Lichenized", "Lichen Parasite", "Lichen Parasite-Lichenized"),
          ),
          by = "seq_id"
        )
    ),
    tar_fst_tbl(
      genus_otu_table,
      nonlichen_otu_table |>
        dplyr::left_join(
          dplyr::select(otu_taxonomy, seq_id, genus),
          by = "seq_id"
        ) |>
        dplyr::summarize(nread = sum(nread), .by = c(genus, sample))
    ),
    tar_fst_tbl(
      family_otu_table,
      nonlichen_otu_table |>
        dplyr::left_join(
          dplyr::select(otu_taxonomy, seq_id, family),
          by = "seq_id"
        ) |>
        dplyr::summarize(nread = sum(nread), .by = c(family, sample))
    ),
    tar_file(
      write_nonlichen_otu_table,
      nonlichen_otu_table %>%
        dplyr::mutate(sample = factor(sample, levels = sample_table$sample)) %>%
        tidyr::pivot_wider(
          names_from = seq_id,
          values_from = nread,
          values_fill = list(nread = 0L)
        ) %>%
        tidyr::complete(sample) %>%
        dplyr::mutate(dplyr::across(where(is.integer), tidyr::replace_na, 0L)) %>%
        tibble::column_to_rownames("sample") %>%
        t() %>% {
          c(
            write_and_return_file(., sprintf("output/nonlichen_otu_table_%s.rds", .conf_level)),
            write_and_return_file(tibble::as_tibble(., rownames = "OTU"),
                                  sprintf("output/nonlichen_otu_table_%s.tsv", .conf_level),
                                  "tsv")
          )
        }
    ),
    tar_file(
      write_genus_otu_table,
      genus_otu_table %>%
        dplyr::mutate(sample = factor(sample, levels = sample_table$sample)) %>%
        tidyr::pivot_wider(
          names_from = genus,
          values_from = nread,
          values_fill = list(nread = 0L)
        ) %>%
        tidyr::complete(sample) %>%
        dplyr::mutate(dplyr::across(where(is.integer), tidyr::replace_na, 0L)) %>%
        tibble::column_to_rownames("sample") %>%
        t() %>% {
          c(
            write_and_return_file(., sprintf("output/genus_otu_table_%s.rds", .conf_level)),
            write_and_return_file(tibble::as_tibble(., rownames = "genus"),
                                  sprintf("output/genus_otu_table_%s.tsv", .conf_level),
                                  "tsv")
          )
        }
    ),
    tar_file(
      write_family_otu_table,
      family_otu_table %>%
        dplyr::mutate(sample = factor(sample, levels = sample_table$sample)) %>%
        tidyr::pivot_wider(
          names_from = family,
          values_from = nread,
          values_fill = list(nread = 0L)
        ) %>%
        tidyr::complete(sample) %>%
        dplyr::mutate(dplyr::across(where(is.integer), tidyr::replace_na, 0L)) %>%
        tibble::column_to_rownames("sample") %>%
        t() %>% {
          c(
            write_and_return_file(., sprintf("output/family_otu_table_%s.rds", .conf_level)),
            write_and_return_file(tibble::as_tibble(., rownames = "family"),
                                  sprintf("output/family_otu_table_%s.tsv", .conf_level),
                                  "tsv")
          )
        }
    )
  )
)