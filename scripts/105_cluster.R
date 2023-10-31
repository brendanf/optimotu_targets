# values for several variables to be used inside the tar_map
rank_meta <- tibble::tibble(
  .rank = tail(TAXRANKS, -1), # phylum:species
  .rank_sym = rlang::syms(.rank),
  .parent_rank = head(TAXRANKS, -1), # kingdom:genus
  .parent_rank_sym = rlang::syms(.parent_rank),
  .super_ranks = purrr::accumulate(.parent_rank, c),
  .parent_taxa = rlang::syms(paste0("taxon_table_", .parent_rank)), # for recursion
  .parent_pseudotaxa = rlang::syms(paste0("pseudotaxon_table_", .parent_rank)) # for recursion
)

#### rank_plan ####
# this ends up inside the reliablility_plan
rank_plan <- tar_map(
  values = rank_meta,
  names = .rank,

  ##### known_taxon_table_{.rank}_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV ID
  #  `kingdom` character : taxon assigned at the kingdom rank
  #  ... character : taxon assigned at additional ranks (down to .rank)
  #
  # taxonomy as known before we start clustering at this rank
  tar_fst_tbl(
    known_taxon_table,
    asv_tax_prob_reads %>%
      dplyr::filter(
        rank == .rank,
        prob >= .prob_threshold
      ) %>%
      dplyr::select(seq_id, .rank_sym := taxon) %>%
      dplyr::left_join(.parent_taxa, ., by = "seq_id") # results from previous rank
  ),

  ##### preclosed_taxon_table_{.rank}_{.conf_level} #####
  # grouped tibble:
  #  `seq_id` character : unique ASV ID
  #  `kingdom` character : taxon assigned at the kingdom rank
  #  ... character : taxon assigned at additional ranks (down to .rank)
  #  `tar_group` integer : grouping variable for dispatch to multiple dynamic
  #    jobs (matches values in .parent_rank)
  #
  # find only the groups which need to be closed-ref clustered
  # i.e. they have some known and some unknown
  tar_fst_tbl(
    preclosed_taxon_table,
    {
      out <- known_taxon_table %>%
        dplyr::group_by(.parent_rank_sym) %>%
        dplyr::filter(any(is.na(.rank_sym)) & !all(is.na(.rank_sym))) %>%
        tar_group()
      # we can't dynamically map over an empty data frame
      # so give a single row.
      if (nrow(out) == 0) {
        known_taxon_table[1,] %>%
          dplyr::group_by(.parent_rank_sym) %>%
          tar_group()
      } else {
        out
      }
    },
    iteration = "group"
  ),

  ##### thresholds_{.rank}_{.conf_level} #####
  # named numeric : optimal clustering threshold at .rank for different taxa
  #   at the rank of .parent_rank (taxa are given by names)
  tar_target(
    thresholds,
    calc_taxon_thresholds(
      rank = .parent_rank,
      conf_level = "plausible",
      taxon_table = known_taxon_table,
      fmeasure_optima = fmeasure_optima
    )
  ),

  ##### clusters_closed_ref_{.rank}_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV id of an "unknown" sequence
  #  `cluster` character : unique ASV id of a "known" sequence
  #
  # find matches within the chosen threshold between "unknown" and "known"
  # sequences
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
    pattern = map(preclosed_taxon_table) # per taxon at rank .parent_rank
  ),

  ##### closedref_taxon_table_{.rank}_{.conf_level} #####
  # grouped tibble:
  #  `seq_id` character : unique ASV ID
  #  `kingdom` character : taxon assigned at the kingdom rank
  #  ... character : taxon assigned at additional ranks (down to .rank)
  #
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


  ##### predenovo_taxon_table_{.rank}_{.conf_level} #####
  # grouped tibble:
  #  `seq_id` character : unique ASV ID
  #  `kingdom` character : taxon assigned at the kingdom rank
  #  ... character : taxon assigned at additional ranks (down to .rank)
  #  `tar_group` integer : grouping variable for dispatch to multiple dynamic
  #    jobs (matches values in .parent_rank)
  #
  # find only the ASVs which need to be de novo clustered
  # i.e. they are still unknown
  tar_fst_tbl(
    predenovo_taxon_table,
    {
      out <- closedref_taxon_table %>%
        dplyr::filter(is.na(.rank_sym)) %>%
        dplyr::group_by(.parent_rank_sym) %>%
        dplyr::filter(dplyr::n() > 1) %>%
        tar_group()
      # we can't dynamically map over an empty data frame
      # so give a single row.
      if (nrow(out) == 0) {
        closedref_taxon_table[1,] %>%
          dplyr::group_by(.parent_rank_sym) %>%
          tar_group()
      } else {
        out
      }
    },
    iteration = "group"
  ),

  ##### denovo_thresholds_{.rank}_{.conf_level} #####
  # named list of named numeric
  #  list names: taxa at .parent_rank
  #  numeric values: optimal clustering thresholds
  #  numeric names: ranks from .rank to species
  #
  # once we are doing de novo clustering at one rank, we will have to do it at
  # all subranks as well, and there is no new information that can change our
  # thresholds.  So it is easiest to just calculate the distance matrix once and
  # do all the clustering now.
  tar_target(
    denovo_thresholds,
    calc_subtaxon_thresholds(
      rank = .parent_rank,
      conf_level = "plausible",
      taxon_table = predenovo_taxon_table,
      fmeasure_optima = fmeasure_optima,
    ),
  ),

  ##### clusters_denovo_{.rank}_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  `kingdom` character : taxonomic assignment at kingdom level
  #  ... character : additional taxonomic assignments down to .parent_rank
  #  ... integer : unique cluster index, for ranks from .rank to species
  tar_target(
    clusters_denovo,
    if (nrow(predenovo_taxon_table) > 1) {
      dplyr::left_join(predenovo_taxon_table, asv_seq, by = "seq_id") %$%
        optimotu::seq_cluster_usearch(
          seq = seq,
          seq_id = seq_id,
          threshold_config = optimotu::threshold_set(
            tryCatch(
              denovo_thresholds[[unique(.parent_rank_sym)]],
              error = function(e) denovo_thresholds[["_NA_"]]
            )
          ),
          clust_config = optimotu::clust_tree(),
          parallel_config = optimotu::parallel_concurrent(2),
          usearch = "bin/usearch",
          usearch_ncpu = local_cpus()
        ) %>%
        t() %>%
        dplyr::as_tibble() %>%
        dplyr::bind_cols(dplyr::select(predenovo_taxon_table, -.rank_sym, -tar_group), .)
    } else {
      c(
        c("seq_id", superranks(.rank)) %>%
          magrittr::set_names(., .) %>%
          purrr::map(~character(0)),
        c(.rank, subranks(.rank)) %>%
          magrittr::set_names(., .) %>%
          purrr::map(~integer(0))
      ) %>%
        tibble::as_tibble()
    },
    pattern = map(predenovo_taxon_table) # per taxon at .parent_rank
  ),

  ##### taxon_table_{.rank}_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV ID
  #  `kingdom` character : taxon assigned at the kingdom rank
  #  ... character : taxon assigned at additional ranks (down to .rank)
  #
  # taxonomy for all named taxa down to .rank
  tar_fst_tbl(
    taxon_table,
    dplyr::filter(closedref_taxon_table, !is.na(.rank_sym))
  ),

  ##### pseudotaxon_table_{.rank}_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  `kingdom` character : taxonomic assignment at kingdom level
  #  ... character : additional taxonomic assignments down to .rank
  #  ... integer : unique cluster index, for ranks below .rank
  #
  # assign "pseudotaxon" names to de novo clusters at this .rank
  tar_fst_tbl(
    pseudotaxon_table,
    dplyr::bind_rows(
      clusters_denovo,
      .parent_pseudotaxa # results from previous rank
    ) %>%
      dplyr::arrange(seq_id) %>% # pseudotaxon numbers are ordered by ASV numbers
      dplyr::mutate(
        .rank_sym := paste(.parent_rank_sym, .rank_sym) %>%
          forcats::fct_inorder() %>%
          forcats::fct_relabel(
            ~names(name_seqs(., paste0("pseudo", .rank, "_")))
          ) %>%
          as.character()
      )
  )
)

#### reliability_plan ####
# repeats entire clustering process for different assignment probability
# thresholds

reliability_meta <- c(
  plausible = 0.5,
  reliable = 0.9
) %>%
  tibble::enframe(name = ".conf_level", value = ".prob_threshold")

reliability_plan <- tar_map(
  values = reliability_meta,
  names = .conf_level,

  ##### taxon_table_kingdom_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV ID
  #  `kingdom` character : taxon assigned at the kingdom rank
  #
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

  ##### pseudotaxon_table_kingdom_{.conf_level} #####
  # NULL
  #
  # this is required because pseudotaxon_table_phylum will try to access its
  # parent, but there are no pseudotaxa at the kingdom level, so it is empty.
  # they are combined with dplyr::bind_row(), which will be fine if we give it
  # NULL instead of a 0-row tibble with the correct columns.
  tar_target(
    pseudotaxon_table_kingdom,
    NULL
  ),

  rank_plan,

  ##### taxon_table_fungi_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  `kingdom` character : taxonomic kingdom assignment
  #  `phylum` character : taxonomic phylum assignment
  #  `class` character : taxonomic class assignment
  #  `order` character : taxonomic order assignment
  #  `family` character : taxonomic family assignment
  #  `genus` character : taxonomic genus assignment
  #  `species` character : taxonomic species assignment
  #
  # combine the taxon table and pseudotaxon table, and select only phyla
  # which a) are named fungal phyla; or b) contain more known fungi (based on
  # unite) than known non-fungi or ASVs with no known kingdom.
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

  ##### write_taxonomy_{.conf_level} #####
  # character : path and file name (.rds)
  #
  # write the ASV taxonomy to a file in the output directory
  tar_file(
    write_taxonomy,
    tibble::column_to_rownames(taxon_table_fungi, "seq_id") %>%
      write_and_return_file(sprintf("output/asv2tax_%s.rds", .conf_level), type = "rds")
  ),

  ##### asv_otu_map_{.conf_level} #####
  tar_fst_tbl(
    asv_otu_map,
    dplyr::semi_join(
      taxon_table_fungi,
      asv_table,
      by = "seq_id"
    ) |>
      dplyr::left_join(
        dplyr::select(otu_taxonomy, OTU = seq_id, species),
        by = "species"
      ) |>
      dplyr::select(ASV = seq_id, OTU)
  ),
  ##### duplicate_species_{.conf_level} #####
  # character : path and file name
  #
  # for testing purposes, write any species which exist in multiple places in
  # the taxonomy.  This file should be empty if everything has gone correctly.
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

  ##### otu_taxonomy_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique OTU ID
  #  `ref_seq_id` character : unique ASV ID of the most prevalent/abundant ASV
  #    in this OTU
  #  `nsample` integer : the number of samples in which this OTU occurs
  #  `nread` integer : the number of reads of this OTU across all samples
  #  `kingdom` character : taxonomic kingdom assignment
  #  `phylum` character : taxonomic phylum assignment
  #  `class` character : taxonomic class assignment
  #  `order` character : taxonomic order assignment
  #  `family` character : taxonomic family assignment
  #  `genus` character : taxonomic genus assignment
  #  `species` character : taxonomic species assignment
  tar_fst_tbl(
    otu_taxonomy,
    asv_table %>%
      dplyr::group_by(seq_id) %>%
      dplyr::mutate(asv_nsample = dplyr::n(), asv_nread = sum(nread)) %>%
      dplyr::inner_join(taxon_table_fungi, by = "seq_id") %>%
      dplyr::group_by(dplyr::across(kingdom:species)) %>%
      dplyr::arrange(dplyr::desc(asv_nsample), dplyr::desc(asv_nread)) %>%
      dplyr::summarize(
        nsample = as.integer(dplyr::n_distinct(sample)),
        nread = sum(nread),
        ref_seq_id = dplyr::first(seq_id)
      ) %>%
      dplyr::arrange(dplyr::desc(nsample), dplyr::desc(nread)) %>%
      name_seqs("OTU", "seq_id") %>%
      dplyr::select(seq_id, ref_seq_id, nsample, nread, everything())
  ),

  ##### write_taxonomy_{.conf_level} #####
  # character : path and file name
  #
  # write the otu taxonomy to a file in the output directory
  tar_file(
    write_otu_taxonomy,
    tibble::column_to_rownames(otu_taxonomy, "seq_id") %>%
      write_and_return_file(sprintf("output/otu_taxonomy_%s.rds", .conf_level), type = "rds")
  ),

  ##### otu_table_sparse_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique OTU id
  #  `sample` character : sample name
  #  `nread` integer : number of reads
  #
  # OTU sample/abundance matrix, in sparse format (0's are not included)
  tar_fst_tbl(
    otu_table_sparse,
    asv_table %>%
      dplyr::inner_join(taxon_table_fungi, by = "seq_id") |>
      dplyr::inner_join(asv_otu_map, by = c("seq_id" = "ASV")) |>
      dplyr::group_by(OTU, sample) |>
      dplyr::summarise(nread = sum(nread), .groups = "drop") |>
      dplyr::select(seq_id = OTU, sample, nread)
  ),

  ##### otu_abund_table_sparse #####
  tar_fst_tbl(
    otu_abund_table_sparse,
    otu_table_sparse |>
      dplyr::left_join(read_counts, by = "sample") |>
      dplyr::left_join(sample_table, by = "sample") |>
      dplyr::group_by(sample) |>
      dplyr::transmute(
        seq_id,
        nread,
        fread = nread/sum(nread),
        w = nread/(nochim2_nread - nospike_nread + 1) * spike_weight,
        .keep = "none"
      ) |>
      dplyr::ungroup()
  ),

  ##### write_otu_table_sparse_{.conf_level} #####
  # character : path and file name (.tsv)
  #
  # write the otu table as a sparse tsv
  tar_file(
    write_otu_table_sparse,
    write_and_return_file(
      dplyr::rename(otu_abund_table_sparse, OTU = seq_id),
      sprintf("output/otu_table_sparse_%s.tsv", .conf_level),
      type = "tsv"
    )
  ),

  ##### otu_table_dense_{.conf_level} #####
  # character (length 2) : path and file name (.rds and .tsv)
  #
  # output the otu table in "dense" format, as required by most community
  # ecology analysis software
  tar_file(
    otu_table_dense,
    otu_table_sparse %>%
      dplyr::mutate(sample = factor(sample, levels = sample_table$sample)) %>%
      tidyr::pivot_wider(names_from = seq_id, values_from = nread, values_fill = list(nread = 0L)) %>%
      tidyr::complete(sample) %>%
      dplyr::mutate(dplyr::across(where(is.integer), \(x) tidyr::replace_na(x, 0L))) %>%
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

  ##### otu_refseq_{.conf_level} #####
  # character : path and file name (.fasta.gz)
  #
  # reference sequence for each OTU
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

  ##### read_counts_{.conf_level} #####
  # tibble:
  #  `sample` character : sample name
  #  `raw_nread` integer : number of read pairs in input files
  #  `trim_nread` integer : number of read pairs remaining after adapter trimming
  #  `filt_nread` integer : number of read pairs remaining after quality filtering
  #  `denoise_nread` numeric? : number of merged reads remaining after denoising
  #  `nochim1_nread` numeric? : number of merged reads remaining after de novo
  #    chimera removal
  #  `nochim2_nread` numeric? : number of merged reads remaining after reference
  #    based chimera removal
  #  `nospike_nread` numeric? : number of merged reads remaining after spike
  #    removal
  #  `full_length` numeric? : number of merged reads remaining after CM scan for
  #    full-length amplicons
  #  `fungi_nread` numeric? : number of merged reads remaining after non-fungi
  #    removal
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
        full_length_read_counts %>%
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
  ##### read_counts_file_{.conf_level} #####
  # character : path and file name (.rds and .tsv)
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

#### clust_plan ####
# the outer plan which contains the reliability_plan

clust_plan <- list(

  ##### asv_known_nonfungi #####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  `kingdom` character : kingdom identification
  #
  # ASVs whose best UNITE match is to a species of known, non-fungal kingdom
  tar_fst_tbl(
    asv_known_nonfungi,
    dplyr::filter(
      asv_unite_kingdom,
      !is.na(kingdom),
      !kingdom %in% c("Fungi", "unspecified", "Eukaryota_kgd_Incertae_sedis")
    )
  ),

  #### asv_known_fungi ####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  `kingdom` character : kingdom identification
  #
  # ASVs whose best UNITE match is to a fungus
  tar_fst_tbl(
    asv_known_fungi,
    dplyr::filter(asv_unite_kingdom, kingdom == "Fungi")
  ),

  #### asv_unknown_kingdom ####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  `kingdom` character : kingdom identification
  #
  # ASVs whose best UNITE match is to a species whose kingdom is unknown
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
