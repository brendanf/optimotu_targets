rank_meta <- tibble::tibble(
  .rank = tail(TAXRANKS, -1), # phylum:species
  .rank_sym = rlang::syms(.rank),
  .parent_rank = head(TAXRANKS, -1), # kingdom:genus
  .parent_rank_sym = rlang::syms(.parent_rank),
  .super_ranks = purrr::accumulate(.parent_rank, c),
  .parent_taxa = rlang::syms(paste0("taxon_table_", .parent_rank))
)

#### rank_plan ####
# this ends up inside the reliablility_plan
rank_plan <- tar_map(
  values = rank_meta,
  names = .rank,
  
  #### known_taxon_table_{.rank}_{.conf_level} ####
  # taxonomy as known before we start clustering at this rank
  tar_fst_tbl(
    known_taxon_table,
    filter_asv_tax_prob_reads %>%
      dplyr::filter(
        rank == .rank,
        prob >= .prob_threshold
      ) %>%
      dplyr::select(ASV, .rank_sym := taxon) %>%
      dplyr::left_join(.parent_taxa, ., by = "ASV")
  ),
  
  #### preclosed_taxon_table_{.rank}_{.conf_level} ####
  # find only the groups which need to be closed-ref clustered
  # i.e. they have some known and some unknown
  tar_fst_tbl(
    preclosed_taxon_table,
    known_taxon_table %>%
      dplyr::group_by(.parent_rank_sym) %>%
      dplyr::filter(any(is.na(.rank_sym)) & !all(is.na(.rank_sym))) %>%
      tar_group(),
    iteration = "group"
  ),
  
  #### thresholds_{.rank}_{.conf_level} ####
  tar_target(
    thresholds,
    calc_taxon_thresholds(
      rank = .parent_rank,
      conf_level = .conf_level,
      taxon_table = known_taxon_table,
      fmeasure_optima = fmeasure_optima
    )
  ),
  
  #### min_threshold_{.rank}_{.conf_level} ####
  tar_target(
    min_threshold,
    min(thresholds)
  ),
  
  #### cluster_closed_ref_{.rank}_{.conf_level} ####
  tar_fst_tbl(
    clusters_closed_ref,
    {
      unknowns <- is.na(preclosed_taxon_table[[.rank]])
      taxon <- preclosed_taxon_table[[.parent_rank]][1]
      if (any(unknowns) && !all(unknowns)) {
        vsearch_usearch_global_closed_ref(
          query = select_sequence(asv_seq, preclosed_taxon_table$ASV[unknowns]),
          ref = select_sequence(asv_seq, preclosed_taxon_table$ASV[!unknowns]),
          threshold = thresholds[taxon]/100
        )
      } else {
        tibble::tibble(ASV = character(), cluster = character())
      }
    },
    pattern = map(preclosed_taxon_table)
  ),
  
  #### closedref_taxon_table_{.rank}_{.conf_level} ####
  # incorporates information from the closed-ref clustering into the
  # taxon table
  tar_fst_tbl(
    closedref_taxon_table,
    dplyr::left_join(
      known_taxon_table,
      clusters_closed_ref,
      by = "ASV"
    ) %>%
      dplyr::left_join(
        dplyr::select(known_taxon_table, cluster = ASV, cluster_taxon = .rank_sym),
        by = "cluster"
      ) %>%
      dplyr::mutate(
        .rank_sym := dplyr::coalesce(.rank_sym, cluster_taxon)
      ) %>%
      dplyr::select(-cluster, -cluster_taxon)
  ),
  
  
  #### predenovo_taxon_table_{.rank}_{.conf_level} ####
  tar_fst_tbl(
    predenovo_taxon_table,
    closedref_taxon_table %>%
      dplyr::filter(is.na(.rank_sym)) %>%
      dplyr::group_by(.parent_rank_sym) %>%
      dplyr::filter(dplyr::n() > 1) %>%
      tar_group(),
    iteration = "group"
  ),
  
  #### clusters_denovo_{.rank}_{.conf_level} ####
  tar_target(
    clusters_denovo,
    dplyr::left_join(predenovo_taxon_table, asv_seq, by = "ASV") %$%
      blastclust_usearch(
        seq = seq,
        seq_id = ASV,
        threshold = tryCatch(
          thresholds[[unique(.parent_rank_sym)]],
          error = function(e) thresholds[["_NA_"]]
        ),
        usearch = "bin/usearch"
      ),
    pattern = map(predenovo_taxon_table)
  ),
  
  #### taxon_table_{.rank}_{.conf_level} ####
  tar_fst_tbl(
    taxon_table,
    tibble::tibble(
      ASV = c(
        trimws(clusters_denovo),
        dplyr::filter(closedref_taxon_table, is.na(.rank_sym)) %>%
          dplyr::group_by(.parent_rank_sym) %>%
          dplyr::filter(dplyr::n() == 1) %>%
          dplyr::pull(ASV)
      ) %>%
        magrittr::extract(order(nchar(.), decreasing = TRUE)),
      cluster = sprintf(
        "pseudo%s_%s",
        .rank,
        formatC(
          seq_along(ASV),
          width = ceiling(log10(length(ASV))) + 1,
          flag = "0"
        )
      )
    ) %>%
      tidyr::separate_rows(ASV, sep = " ") %>%
      dplyr::right_join(closedref_taxon_table, by = "ASV") %>%
      dplyr::mutate(.rank_sym := dplyr::coalesce(.rank_sym, cluster)) %>%
      dplyr::select(-cluster)
  )
)

#### reliability_plan ####

reliability_meta <- c(
  plausible = 0.5,
  reliable = 0.9
) %>%
  tibble::enframe(name = ".conf_level", value = ".prob_threshold")

reliability_plan <- tar_map(
  values = reliability_meta,
  names = .conf_level,
  
#  #### PROTAX_unassigned_phylum ####
#  tar_target(
#    PROTAX_unknown_phylum,
#    asv_tax_prob_reads %>%
#      dplyr::filter(
#        rank == "phylum",
#        prob < threshold_meta$prob_threshold
#      ),
#    pattern = map(threshold_meta)
#  ),
  
  #### taxon_table_kingdom_{.conf_level} ####
  # values for other ranks are calculated recursively
  # this should be everything, because PROTAX-fungi assigns all sequences
  # 100% probability of being fungi
  tar_fst_tbl(
    taxon_table_kingdom,
    filter_asv_tax_prob_reads %>%
      dplyr::filter(rank == "kingdom") %>%
      dplyr::mutate(
        taxon = ifelse(prob < .prob_threshold, NA_character_, taxon)
      ) %>%
      dplyr::select(ASV, kingdom = taxon)
  ),
  
  #### taxon_table_fungi_{.conf_level} ####
  tar_fst_tbl(
    taxon_table_fungi,
    taxon_table_species %>%
      dplyr::mutate(
        known_nonfungus = ASV %in% sh_known_nonfungi$ASV,
        known_fungus = ASV %in% sh_known_fungi$ASV,
        unknown_kingdom = ASV %in% sh_unknown_kingdom$ASV
      ) %>%
      dplyr::group_by(phylum) %>%
      dplyr::filter(
        !startsWith(phylum, "pseudophylum") |
          sum(known_fungus) > sum(known_nonfungus) + sum(unknown_kingdom)
      )
  ),
  
  rank_plan,
  
  #### min_threshold_{.conf_level} ####
  tar_combine(
    min_threshold,
    rank_plan$min_threshold,
    command = min(!!!.x),
    use_names = FALSE
  ),
  
  #### chosen_taxonomy_{.conf_level} ####
  tar_fst_tbl(
    chosen_taxonomy,
    taxon_table_species %>%
      dplyr::group_by(phylum) %>%
      dplyr::filter(!any(ASV %in% sh_known_nonfungi$ASV)) %>%
      dplyr::arrange(as.numeric(substr(ASV, start = 4, stop = 100)))
  ),
  #### write_taxonomy_{.conf_level} ####
  tar_file(
    write_taxonomy,
    tibble::column_to_rownames(chosen_taxonomy, "ASV") %>%
      write_and_return_file(sprintf("output/asv2tax_%s.rds", .conf_level), type = "rds")
  ),
  #### duplicate_species_{.conf_level} ####
  tar_file(
    duplicate_species,
    dplyr::group_by(chosen_taxonomy, species) %>%
      dplyr::filter(dplyr::n_distinct(phylum, class, order, family, genus) > 1) %>%
      dplyr::left_join(asv_seq, by = "ASV") %>%
      dplyr::mutate(
        classification = paste(phylum, class, order, family, genus, sep = ";") %>%
          ifelse(
            length(.) > 0L,
            sub(Biobase::lcPrefix(.), "", .),
            .
          ),
        name = sprintf("%s (%s) %s", species, classification, ASV)
      ) %>%
      dplyr::arrange(name) %>%
      dplyr::ungroup() %>%
      dplyr::select(name, seq) %>%
      tibble::deframe() %>%
      Biostrings::DNAStringSet() %>%
      write_and_return_file(sprintf("output/duplicates_%s.fasta", .conf_level))
  )
)

clust_plan <- list(
  
  #### threshold_meta ####
  # tar_fst_tbl(
  #    threshold_meta,
  #    dplyr::select(fmeasure_optima, rank, superrank, supertaxon, threshold, conf_level) %>%
  #       dplyr::left_join(reliability, by = "conf_level") %>%
  #       dplyr::arrange(conf_level, threshold) %>%
  #       dplyr::left_join(parent_rank, by = "rank")
  # ),
  
  #### sh_known_nonfungi ####
  tar_fst_tbl(
    sh_known_nonfungi,
    dplyr::filter(unite_matches_out_97, kingdom != "Fungi" | genus == "Ciliophora") %>%
      dplyr::rename(ASV = seq_accno)
  ),
  
  #### sh_known_fungi ####
  tar_fst_tbl(
    sh_known_fungi,
    dplyr::filter(unite_matches_out_97, kingdom == "Fungi" & genus != "Ciliophora") %>%
      dplyr::rename(ASV = seq_accno)
  ),
  
  #### sh_unknown_kingdom ####
  tar_target(
    sh_unknown_kingdom,
    unite_matches_out_97 %>%
      dplyr::filter(
        is.na(kingdom) |
          (kingdom == "Eukaryota_kingdom_incertae_sedis" &
             phylum == "unidentified")
      ) %>%
      dplyr::rename(ASV = seq_accno)
  ),
  
  #### filter_asv_tax_prob_reads ####
  tar_fst_tbl(
    filter_asv_tax_prob_reads,
    asv_tax_prob_reads %>%
      # remove ASVs which were excluded by Unite as too short or chimeric
      dplyr::anti_join(unite_excluded, by = c("ASV" = "seq_accno"))
    # remove ASVs which match a non-fungal SH
    #dplyr::anti_join(sh_known_nonfungi, by = "ASV") %>%
    # remove ASVs which do not match an SH with a kingdom, and which
    # PROTAX could not assign to (fungal) phylum
    # dplyr::anti_join(
    #    dplyr::inner_join(sh_unknown_kingdom, PROTAX_unknown_phylum,
    #                      by = "ASV"),
    #    by = "ASV"
    # )
  ),
  
  reliability_plan
)
