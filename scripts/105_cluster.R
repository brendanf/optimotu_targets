rank_meta <- tibble::tibble(
  .rank = tail(TAXRANKS, -1), # phylum:species
  .rank_sym = rlang::syms(.rank),
  .parent_rank = head(TAXRANKS, -1), # kingdom:genus
  .parent_rank_sym = rlang::syms(.parent_rank),
  .super_ranks = purrr::accumulate(.parent_rank, c),
  .parent_taxa = rlang::syms(paste0("taxon_table_", .parent_rank)),
  .parent_pseudotaxa = rlang::syms(paste0("pseudotaxon_table_", .parent_rank))
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
    asv_tax_prob_reads %>%
      dplyr::filter(
        rank == .rank,
        prob >= .prob_threshold
      ) %>%
      dplyr::select(seq_id, .rank_sym := taxon) %>%
      dplyr::left_join(.parent_taxa, ., by = "seq_id")
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
      conf_level = "plausible",
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
          query = select_sequence(asv_seq, preclosed_taxon_table$seq_id[unknowns]),
          ref = select_sequence(asv_seq, preclosed_taxon_table$seq_id[!unknowns]),
          threshold = thresholds[taxon]/100
        )
      } else {
        tibble::tibble(seq_id = character(), cluster = character())
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
      by = "seq_id"
    ) %>%
      dplyr::left_join(
        dplyr::select(known_taxon_table, cluster = seq_id, cluster_taxon = .rank_sym),
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
  
  tar_target(
    denovo_thresholds,
    calc_subtaxon_thresholds(
      rank = .parent_rank,
      conf_level = "plausible",
      taxon_table = predenovo_taxon_table,
      fmeasure_optima = fmeasure_optima,
    ),
  ),
  
  #### clusters_denovo_{.rank}_{.conf_level} ####
  tar_target(
    clusters_denovo,
    dplyr::left_join(predenovo_taxon_table, asv_seq, by = "seq_id") %$%
      optimotu::usearch_single_linkage(
        seq = seq,
        seq_id = seq_id,
        thresholds = tryCatch(
          denovo_thresholds[[unique(.parent_rank_sym)]],
          error = function(e) denovo_thresholds[["_NA_"]]
        ),
        usearch = "bin/usearch"
      ) %>%
      t() %>%
      dplyr::as_tibble() %>%
      dplyr::bind_cols(dplyr::select(predenovo_taxon_table, -.rank_sym, -tar_group), .),
    pattern = map(predenovo_taxon_table)
  ),
  
  #### taxon_table_{.rank}_{.conf_level} ####
  tar_fst_tbl(
    taxon_table,
    dplyr::filter(closedref_taxon_table, !is.na(.rank_sym))
  ),
  
  #### pseudotaxon_table_{.rank}_{.conf_level} ####
  tar_fst_tbl(
    pseudotaxon_table,
    dplyr::bind_rows(
      clusters_denovo,
      .parent_pseudotaxa
    ) %>%
      dplyr::mutate(
        .rank_sym := paste(.parent_rank_sym, .rank_sym) %>%
          forcats::fct_relabel(
            ~names(name_seqs(., paste0("pseudo", .rank, "_")))
          ) %>%
          as.character()
      ) 
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
  
  #### taxon_table_kingdom_{.conf_level} ####
  # values for other ranks are calculated recursively
  # this should be everything, because PROTAX-fungi assigns all sequences
  # 100% probability of being fungi
  tar_fst_tbl(
    taxon_table_kingdom,
    asv_tax_prob_reads %>%
      dplyr::filter(rank == "kingdom") %>%
      dplyr::mutate(
        taxon = ifelse(prob < .prob_threshold, NA_character_, taxon)
      ) %>%
      dplyr::select(seq_id, kingdom = taxon)
  ),
  
  #### pseudotaxon_table_kingdom_{.conf_level} ####
  # this is required because pseudotaxon_table_phylum will try to access its
  # parent, but thre are no pseudotaxa at the kingdom level, so it is empty.
  # they are combined with dplyr::bind_row(), which will be fine if we give it
  # NULL instead of a 0-row tibble with the correct columns.
  tar_target(
    pseudotaxon_table_kingdom,
    NULL
  ),
  
  #### taxon_table_fungi_{.conf_level} ####
  tar_fst_tbl(
    taxon_table_fungi,
    dplyr::bind_rows(
      taxon_table_species,
      pseudotaxon_table_species
    ) %>%
      dplyr::mutate(
        known_nonfungus = seq_id %in% asv_known_nonfungi$seq_id,
        known_fungus = seq_id %in% asv_known_fungi$seq_id,
        unknown_kingdom = seq_id %in% asv_unknown_kingdom$seq_id
      ) %>%
      dplyr::group_by(phylum) %>%
      dplyr::filter(
        !startsWith(phylum, "pseudophylum") |
          sum(known_fungus) > sum(known_nonfungus) + sum(unknown_kingdom)
      ) %>%
      dplyr::select(!where(is.logical)) %>%
      dplyr::arrange(seq_id)
  ),
  
  rank_plan,
  
  #### min_threshold_{.conf_level} ####
  tar_combine(
    min_threshold,
    rank_plan$min_threshold,
    command = min(!!!.x),
    use_names = FALSE
  ),
  
  #### write_taxonomy_{.conf_level} ####
  tar_file(
    write_taxonomy,
    tibble::column_to_rownames(taxon_table_fungi, "seq_id") %>%
      write_and_return_file(sprintf("output/asv2tax_%s.rds", .conf_level), type = "rds")
  ),
  #### duplicate_species_{.conf_level} ####
  tar_file(
    duplicate_species,
    dplyr::group_by(taxon_table_fungi, species) %>%
      dplyr::filter(dplyr::n_distinct(phylum, class, order, family, genus) > 1) %>%
      dplyr::left_join(asv_seq, by = "seq_id") %>%
      dplyr::mutate(
        classification = paste(phylum, class, order, family, genus, sep = ";") %>%
          ifelse(
            length(.) > 0L,
            sub(Biobase::lcPrefix(.), "", .),
            .
          ),
        name = sprintf("%s (%s) %s", species, classification, seq_id)
      ) %>%
      dplyr::arrange(name) %>%
      dplyr::ungroup() %>%
      dplyr::select(name, seq) %>%
      tibble::deframe() %>%
      Biostrings::DNAStringSet() %>%
      write_and_return_file(sprintf("output/duplicates_%s.fasta", .conf_level))
  ),
  
  #### otu_taxonomy_{.conf_level} ####
  tar_fst_tbl(
    otu_taxonomy,
    asv_table %>%
      dplyr::group_by(seq_id) %>%
      dplyr::mutate(asv_nsample = dplyr::n(), asv_nread = sum(nread)) %>%
      dplyr::inner_join(taxon_table_fungi, by = "seq_id") %>%
      dplyr::group_by(dplyr::across(kingdom:species)) %>%
      dplyr::arrange(dplyr::desc(asv_nsample), dplyr::desc(asv_nread)) %>%
      dplyr::summarize(
        nsample = dplyr::n_distinct(sample),
        nread = sum(nread),
        ref_seq_id = dplyr::first(seq_id)
      ) %>%
      dplyr::arrange(dplyr::desc(nsample), dplyr::desc(nread)) %>%
      name_seqs("OTU", "seq_id") %>%
      dplyr::select(seq_id, ref_seq_id, nsample, nread, everything())
  ),
  #### write_taxonomy_{.conf_level} ####
  tar_file(
    write_otu_taxonomy,
    tibble::column_to_rownames(otu_taxonomy, "seq_id") %>%
      write_and_return_file(sprintf("output/otu_taxonomy_%s.rds", .conf_level), type = "rds")
  ),
  
  #### otu_table_{.conf_level} ####
  tar_fst_tbl(
    otu_table_sparse,
    asv_table %>%
      dplyr::inner_join(taxon_table_fungi, by = "seq_id") %>%
      dplyr::inner_join(
        dplyr::select(otu_taxonomy, OTU = seq_id, kingdom:species),
        by = TAXRANKS
      ) %>%
      dplyr::group_by(OTU, sample) %>%
      dplyr::summarise(nread = sum(nread), .groups = "drop") %>%
      dplyr::rename(seq_id = OTU)
  ),
  
  #### otu_table_dense_{.conf_level} ####
  tar_file(
    otu_table_dense,
    otu_table_sparse %>%
      dplyr::mutate(sample = factor(sample, levels = sample_table$sample)) %>%
      tidyr::pivot_wider(names_from = seq_id, values_from = nread, values_fill = list(nread = 0L)) %>%
      tidyr::complete(sample) %>%
      dplyr::mutate(dplyr::across(where(is.integer), tidyr::replace_na, 0L)) %>%
      tibble::column_to_rownames("sample") %>%
      t() %>% {
        c(
          write_and_return_file(., sprintf("output/otu_table_%s.rds", .conf_level)),
          write_and_return_file(tibble::as_tibble(., rownames = "OTU"),
                                sprintf("output/otu_table_%s.tsv", .conf_level),
                                "tsv")
        )
      }
  ),
  
  #### otu_refseq_{.conf_level} ####
  tar_file(
    otu_refseq,
    otu_taxonomy %>%
      dplyr::ungroup() %>%
      dplyr::left_join(asv_seq, by = c("ref_seq_id" = "seq_id")) %>%
      dplyr::select(seq_id, seq) %>%
      tibble::deframe() %>%
      Biostrings::DNAStringSet() %>%
      write_and_return_file(
        sprintf("output/otu_%s.fasta.gz", .conf_level),
        compress = TRUE
      )
  ),
  
  #### read_counts_{.conf_level} ####
  tar_fst_tbl(
    read_counts,
    dada2_meta %>%
      dplyr::mutate(fastq_file = file.path(raw_path, fastq_R1)) %>%
      dplyr::left_join(raw_read_counts, by = "fastq_file") %>%
      dplyr::left_join(trim_read_counts, by = "trim_R1") %>%
      dplyr::left_join(filt_read_counts, by = "filt_R1") %>%
      dplyr::mutate(filt_key = sub("_R[12]_filt\\.fastq\\.gz", "", filt_R1)) %>%
      dplyr::left_join(denoise_read_counts, by = "filt_key") %>%
      dplyr::left_join(nochim1_read_counts, by = "filt_key") %>%
      dplyr::left_join(
        nochim2_read_counts %>%
          dplyr::summarize(dplyr::across(everything(), sum), .by = filt_key),
        by = "filt_key"
      ) %>%
      dplyr::left_join(
        nospike_read_counts %>%
          dplyr::summarize(dplyr::across(everything(), sum), .by = filt_key),
        by = "filt_key"
      ) %>%
      dplyr::left_join(
        dplyr::group_by(otu_table_sparse, sample) %>%
          dplyr::summarize(fungi_nread = sum(nread)),
        by = "sample"
      ) %>%
      tidyr::replace_na(list(fungi_nread = 0L)) %>%
      dplyr::select(sample, raw_nread, trim_nread, filt_nread, denoise_nread,
                    nochim1_nread, nochim2_nread, nospike_nread, fungi_nread)
  ),
  #### read_counts_file_{.conf_level} ####
  tar_file(
    read_counts_file,
    c(
      write_and_return_file(
        read_counts,
        sprintf("output/read_counts_%s.rds", .conf_level),
        "rds"
      ),
      write_and_return_file(
        read_counts,
        sprintf("output/read_counts_%s.tsv", .conf_level),
        "tsv"
      )
    )
  )
)

clust_plan <- list(
  
  #### asv_known_nonfungi ####
  tar_fst_tbl(
    asv_known_nonfungi,
    dplyr::filter(
      asv_unite_kingdom,
      !is.na(kingdom),
      !kingdom %in% c("Fungi", "unspecified", "Eukaryota_kgd_Incertae_sedis")
    )
  ),
  
  #### asv_known_fungi ####
  tar_fst_tbl(
    asv_known_fungi,
    dplyr::filter(asv_unite_kingdom, kingdom == "Fungi")
  ),
  
  #### asv_unknown_kingdom ####
  tar_target(
    asv_unknown_kingdom,
    dplyr::filter(
      asv_unite_kingdom,
      is.na(kingdom) |
        kingdom %in% c("unspecified", "Eukaryota_kgd_Incertae_sedis")
    )
  ),
  
  reliability_plan
)
