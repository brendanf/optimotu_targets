occurrence_plan <- list(
  #### climate_zone_file ####
  # `character` representing a filename
  # TSV giving climate zone for each site
  tar_file(
    climate_zone_file,
    "metadata/cz.tsv"
  ),
  #### climate_zones ####
  # `tibble` with columns:
  #  `site` character: 3-character site abbreviation
  #  `meta.CZ` character: one of "Polar-Continental", "Temperate", or
  #    "Tropical-Subtropical"
  tar_fst_tbl(
    climate_zones,
    readr::read_tsv(climate_zone_file, col_types = "cc")
  ),
  
  #### map over confidence levels ####
  tar_map(
    # also map over some previously mapped targets
    values = tibble::tibble(
      .conf_level = c("plausible", "reliable"),
      otu_table_sparse = paste0("otu_table_sparse_", .conf_level) %>%
        rlang::syms(),
      otu_taxonomy = paste0("otu_taxonomy_", .conf_level) %>%
        rlang::syms(),
      asv_otu_map = paste0("asv_otu_map_", .conf_level) %>%
        rlang::syms()
    ),
    names = .conf_level,
    ##### otu_table_sparse_site_{.conf_level} #####
    # `tibble` with columns:
    #  
    tar_fst_tbl(
      otu_table_sparse_site,
      dplyr::mutate(otu_table_sparse, site = substr(sample, 1, 3)) |>
        dplyr::group_by(sample) |>
        dplyr::mutate(fread = nread/sum(nread)) |>
        dplyr::group_by(site) |>
        dplyr::mutate(nsample = dplyr::n_distinct(sample)) |>
        dplyr::group_by(seq_id, site) |>
        dplyr::summarize(
          nread = sum(nread),
          fread = sum(fread/nsample),
          nocc_sample = dplyr::n(),
          .groups = "drop"
        )
    ),
    tar_fst_tbl(
      otu_table_sparse_cz,
      dplyr::inner_join(
        otu_table_sparse_site,
        climate_zones,
        by = c("site" = "meta.site")
      ) %>%
        dplyr::group_by(meta.CZ) %>%
        dplyr::mutate(nsite = dplyr::n_distinct(site)) %>%
        dplyr::group_by(OTU, meta.CZ) %>%
        dplyr::summarize(
          nread = sum(nread),
          fread = sum(fread/nsite),
          nocc_sample = sum(nocc_sample),
          nocc_site = dplyr::n(),
          .groups = "drop"
        ) %>%
        dplyr::group_by(meta.CZ) %>%
        dplyr::mutate(
          focc_sample = nocc_sample / sum(nocc_sample),
          focc_site = nocc_site / sum(nocc_site)
        ) %>%
        dplyr::ungroup()
    ),
    #### otu_unknown_prob_{.conf_level} ####
    tar_fst_tbl(
      otu_unknown_prob,
      dplyr::select(otu_taxonomy, -ref_seq_id, -nsample, -nread) |>
        tidyr::pivot_longer(
          kingdom:species,
          names_to = "rank",
          values_to = "taxon",
          names_transform = list(rank = rank2factor)
        ) |>
        dplyr::filter(rank != "kingdom") |>
        dplyr::left_join(asv_otu_map, by = c("seq_id" = "OTU"), multiple = "all") |>
        dplyr::left_join(asv_unknown_prob, by = c("ASV" = "seq_id", "rank")) |>
        tidyr::replace_na(list(prob_unk = 1)) |>
        dplyr::group_by(seq_id, rank, taxon) |>
        dplyr::summarize(
          prob_unknown = min(prob_unk),
          nasv = dplyr::n()
        )
      ),
    
    #### otu_unknown_by_cz ####
    tar_fst_tbl(
      otu_unknown_by_cz,
      otu_unknown_prob |>
        dplyr::inner_join(
          otu_table_sparse_cz,
          by = "seq_id"
        ) |>
        dplyr::group_by(rank, taxon, meta.CZ) |>
        dplyr::summarize(
          notu = dplyr::n(),
          prob_unknown = min(prob_unknown),
          dplyr::across(
            c(nasv, nread, fread, nocc_sample, nocc_site,
              focc_sample, focc_site),
            sum
          )
        )
    )
  ),
  #### otu_unknown_plot ####
  tar_target(
    otu_unknown_plot,
    otu_unknown_by_cz_reliable %>%
      dplyr::mutate(
        type = dplyr::case_when(
          !startsWith(taxon, "pseudo") ~ "Known taxa",
          prob_unknown >= 0.9 ~ "Novel taxa",
          TRUE ~ "Uncertain"
        ) |>
          factor(levels = c("Known taxa", "Uncertain", "Novel taxa"))
      ) |>
      dplyr::group_by(rank, meta.CZ, type) |>
      dplyr::summarize(
        dplyr::across(tidyselect::matches("[nf](otu|read|asv|occ).*"), sum)
      ) |>
      dplyr::mutate(
        meta.CZ = dplyr::recode(
          meta.CZ,
          "Polar-Continental" = "Pol-Cont",
          "Temperate" = "Temp",
          "Tropical-Subtropical" = "Trop-Sub"
        ),
        rank = forcats::fct_rev(rank) |>
          forcats::fct_relabel(stringr::str_to_title)
      ) |>
      ggplot(aes(x = meta.CZ, weight = focc_sample, group = type,
                 fill = type)) +
      geom_bar() +
      facet_wrap(~rank, nrow = 1) +
      scale_fill_manual(
        values = list(
          "Known taxa" = "forestgreen",
          "Uncertain" = "tomato",
          "Novel taxa" = "steelblue4"
        ),
        name = NULL
      ) +
      theme(
        axis.text.x = element_text(angle = -45, hjust = 0),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank()
      ) +
      xlab(NULL) +
      scale_y_continuous(
        name = "Fraction of occurences",
        expand = expansion(),
        labels = scales::label_percent()
      ),
    packages = "ggplot2"
  ),
  tar_map(
    values = list(format = c("pdf", "svg", "png")),
    tar_file(
      otu_unknown_plot_file,
      write_and_return_file(
        otu_unknown_plot,
        sprintf("output/known_unknown_fig.%s", format),
        width = 6, height = 3
      )
    )
  ),
)
