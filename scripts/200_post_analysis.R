occurrence_plan <- list(
  #### funguild_db ####
  tar_fst_tbl(
    funguild_db,
    FUNGuildR::get_funguild_db()
  ),
  
  #### lifestyle_db_file ####
  tar_file_fast(
    lifestyle_db_file,
    "data/lifestyle/Fung_LifeStyle_Data.RDS",
    deployment = "main"
  ),
  
  #### lifestyle_db ####
  tar_fst_tbl(
    lifestyle_db,
    taxonomy_new |>
      # take genera
      dplyr::filter(rank == 6) |>
      #split classification string into ranks
      tidyr::separate(
        classification,
        c("kingdom", "phylum", "class", "order", "family", "genus"),
        sep = ","
      ) |>
      # genera have mycobank number appaneded to name; remove it
      dplyr::mutate(genus = sub("_[0-9]+", "", genus)) |>
      # in some cases there are multiple entries for each genus.
      # take the one with the highest prior (i.e. highest number of species in MB)
      dplyr::filter(seq_along(prior) == which.max(prior), .by = genus) |>
      
      dplyr::inner_join(
        readRDS(lifestyle_db_file) |>
          dplyr::mutate(genus = sub(" .*", "", taxon)),
        by = "genus",
        multiple = "all"
      ) |>
      (
        \(x) dplyr::bind_rows(
          dplyr::transmute(
            x,
            taxon,
            taxonomicLevel = ifelse(grepl(" ", taxon, fixed = TRUE), 20L, 13L),
            trophicMode = NA_character_,
            guild = chartr(" ", ",", guild),
            citationSource,
            searchkey = paste0("@", sub("[_ ]", "@", taxon), "@")
          ),
          dplyr::summarize(
            x,
            guild = paste(
              setdiff(
                unique(unlist(strsplit(guild, "[ ,]"))),
                c("NA", NA_character_)
              ),
              collapse = ","
            ),
            .by = genus
          ) |>
            dplyr::transmute(
              taxon = genus,
              taxonomicLevel = 13L,
              trophicMode = NA_character_,
              guild,
              citationSource = "combined from species-level annotations",
              searchkey = paste0("@", taxon, "@")
            ) |>
            dplyr::anti_join(x, by = "taxon"),
          dplyr::filter(x, !startsWith(family, "dummy")) |>
            dplyr::summarize(
              guild = paste(
                setdiff(
                  unique(unlist(strsplit(guild, "[ ,]"))),
                  c("NA", NA_character_)
                ),
                collapse = ","
              ),
              .by = family
            ) |>
            dplyr::transmute(
              taxon = family,
              taxonomicLevel = 9L,
              trophicMode = NA_character_,
              guild,
              citationSource = "combined from genus-level annotations",
              searchkey = paste0("@", taxon, "@")
            )
        )
      )(),
    deployment = "main"
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
    
    tar_map(
      values = tibble::tibble(
        .guild_db = rlang::syms(c("funguild_db", "lifestyle_db")),
        .guild = c("funguild", "carlos")
      ),
      names = .guild,
    
      ###### otu_guild_{.guild_db}_{.conf_level} ######
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
          FUNGuildR::funguild_assign(db = .guild_db) |>
          dplyr::select(seq_id, guild),
        deployment = "main"
      ),
      ###### write_otu_guild_{.guild_db}_{.conf_level} ######
      tar_file_fast(
        write_otu_guild,
        write_and_return_file(
          otu_guild,
          sprintf("output/otu_guilds_%s_%s.tsv", .guild, .conf_level),
          type = "tsv"
        ),
        deployment = "main"
      )
    ),    
    ##### otu_unknown_prob_{.conf_level} #####
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
        ),
      deployment = "main"
    )        
  )
)
