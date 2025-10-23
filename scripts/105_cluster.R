# values for several variables to be used inside the tar_map
rank_meta <- tibble::tibble(
  .rank = tail(optimotu.pipeline::tax_ranks(), -1), # e.g. phylum:species
  .rank_sym = rlang::syms(.rank),
  .parent_rank = head(optimotu.pipeline::tax_ranks(), -1), # e.g. kingdom:genus
  .parent_rank_sym = rlang::syms(.parent_rank),
  .super_ranks = purrr::accumulate(.parent_rank, c),
  .parent_taxa = rlang::syms(paste0("taxon_table_", .parent_rank)), # for recursion
  .parent_pseudotaxa = rlang::syms(paste0("pseudotaxon_table_", .parent_rank)) # for recursion
)

taxon_table_TIP_RANK <- rlang::sym(
  sprintf("taxon_table_%s", optimotu.pipeline::tip_rank())
)
pseudotaxon_table_TIP_RANK <- rlang::sym(
  sprintf("pseudotaxon_table_%s", optimotu.pipeline::tip_rank())
)

seq_to_cluster_file <- quote(asv_taxsort_seq)
seq_to_cluster_file_index <- quote(asv_taxsort_seq_index)

if (optimotu.pipeline::do_model_align()) {
  seq_to_cluster_file <- quote(aligned_taxsort_seq)
  seq_to_cluster_file_index <- quote(aligned_taxsort_seq_index)
}

#### rank_plan ####
# this ends up inside the reliablility_plan
rank_plan <- tar_map(
  values = rank_meta,
  names = .rank,

  ##### known_taxon_table_{.rank}_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV ID
  #  {ROOT_RANK} character : taxon assigned at ROOT_RANK (e.g. kingdom)
  #  ... character : taxon assigned at additional ranks (down to .rank)
  #
  # taxonomy as known before we start clustering at this rank
  tar_fst_tbl(
    known_taxon_table,
    asv_tax_prob_reads |>
      dplyr::filter(
        rank == .rank,
        prob >= .prob_threshold,
        taxon != "unk"
      ) |>
      dplyr::select(seq_id, .rank_sym := taxon) |>
      dplyr::left_join(.parent_taxa, y = _, by = "seq_id"), # results from previous rank
    deployment = "main"
  ),

  ##### preclosed_taxon_table_large_{.rank}_{.conf_level} #####
  # grouped tibble:
  #  `seq_id` character : unique ASV ID
  #  {ROOT_RANK} character : taxon assigned at ROOT_RANK (e.g. kingdom)
  #  ... character : taxon assigned at additional ranks (down to .rank)
  #  `seq_idx` integer: index of the sequence in asv_taxsort_seq
  #  `tar_group` integer : grouping variable for dispatch to multiple dynamic
  #    jobs (matches values in .parent_rank)
  #
  # find only the groups which need to be closed-ref clustered
  # i.e. they have some known and some unknown
  # This version returns "large" taxa, i.e. those with enough sequences that
  # parallelization is worthwhile.
  tar_fst_tbl(
    preclosed_taxon_table_large,
    optimotu.pipeline::large_preclosed_taxon_table(
      known_taxon_table = known_taxon_table,
      asv_taxsort = asv_taxsort,
      rank = .rank,
      parent_rank = .parent_rank,
      tax_ranks = !!optimotu.pipeline::tax_ranks()
    ),
    iteration = "group",
    resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
  ),

  ##### preclosed_taxon_table_small_{.rank}_{.conf_level} #####
  # grouped tibble:
  #  `seq_id` character : unique ASV ID
  #  {ROOT_RANK} character : taxon assigned at ROOT_RANK (e.g. kingdom)
  #  ... character : taxon assigned at additional ranks (down to .rank)
  #  `seq_idx` integer: index of the sequence in asv_taxsort_seq
  #  `tar_group` integer : grouping variable for dispatch to multiple dynamic
  #    jobs (matches values in .parent_rank)
  #
  # find only the groups which need to be closed-ref clustered
  # i.e. they have some known and some unknown
  # This version returns "small" taxa, i.e. those with too few sequences for
  # parallelization to be worthwhile, in batches to keep the total number of
  # targets smaller.
  tar_fst_tbl(
    preclosed_taxon_table_small,
    optimotu.pipeline::small_preclosed_taxon_table(
      known_taxon_table = known_taxon_table,
      asv_taxsort = asv_taxsort,
      rank = .rank,
      parent_rank = .parent_rank,
      tax_ranks = !!optimotu.pipeline::tax_ranks()
    ),
    iteration = "group",
    resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
  ),

  ##### thresholds_{.rank}_{.conf_level} #####
  # named numeric : optimal clustering threshold at .rank for different taxa
  #   at the rank of .parent_rank (taxa are given by names)
  tar_target(
    thresholds,
    optimotu::calc_taxon_thresholds(
      rank = .parent_rank,
      taxon_table = known_taxon_table,
      optima = cluster_optima,
      ranks = !!optimotu.pipeline::tax_ranks(),
      measure = !!optimotu.pipeline::cluster_measure()
    ),
    deployment = "main"
  ),

  ##### clusters_closed_ref_large_{.rank}_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV id of an "unknown" sequence
  #  `ref_id` character : unique ASV id of a "known" sequence
  #  `dist` double : distance between the two sequences
  #
  # Find matches within the chosen threshold between "unknown" and "known"
  # sequences.  This version for "large" clusters which should use parallel
  # distance calculations.
  tar_fst_tbl(
    clusters_closed_ref_large,
    optimotu.pipeline::do_closed_ref_cluster(
      preclosed_taxon_table = preclosed_taxon_table_large,
      rank = .rank,
      parent_rank = .parent_rank,
      seq_file = !!seq_to_cluster_file,
      seq_file_index = !!seq_to_cluster_file_index,
      thresholds = thresholds,
      dist_config = !!(
        if (optimotu.pipeline::cluster_dist_config()$method == "usearch") {
          substitute(
            update(dc, usearch_ncpu = optimotu.pipeline::local_cpus()),
            list(dc = optimotu.pipeline::cluster_dist_config())
          ))
        } else {
          optimotu.pipeline::cluster_dist_config()
        }
      ),
      parallel_config = !!(
        if (optimotu.pipeline::cluster_dist_config()$method == "usearch") {
          quote(optimotu::parallel_concurrent(2))
        } else {
          quote(optimotu::parallel_concurrent(optimotu.pipeline::local_cpus()))
        }
      )
    ),
    pattern = map(preclosed_taxon_table_large), # per taxon at rank .parent_rank
    resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
  ),

  ##### clusters_closed_ref_small_{.rank}_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV id of an "unknown" sequence
  #  `ref_id` character : unique ASV id of a "known" sequence
  #  `dist` double : distance between the two sequences
  #
  # Find matches within the chosen threshold between "unknown" and "known"
  # sequences.  This version for "small" clusters which should use serial
  # distance calculations.
  tar_fst_tbl(
    clusters_closed_ref_small,
    optimotu.pipeline::do_closed_ref_cluster(
      preclosed_taxon_table = preclosed_taxon_table_small,
      rank = .rank,
      parent_rank = .parent_rank,
      seq_file = !!seq_to_cluster_file,
      seq_file_index = !!seq_to_cluster_file_index,
      thresholds = thresholds,
      dist_config = !!(
        if (optimotu.pipeline::cluster_dist_config()$method == "usearch") {
          substitute(
            update(dc, usearch_ncpu = 1),
            list(dc = optimotu.pipeline::cluster_dist_config())
          ))
        } else {
          optimotu.pipeline::cluster_dist_config()
        }
      ),
      parallel_config = optimotu::parallel_concurrent(1)
    ),
    pattern = map(preclosed_taxon_table_small), # per batch of taxon at rank .parent_rank
    resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
  ),

  ##### closedref_taxon_table_{.rank}_{.conf_level} #####
  # grouped tibble:
  #  `seq_id` character : unique ASV ID
  #  {ROOT_RANK} character : taxonomic assignment at ROOT_RANK (e.g. kingdom)
  #  ... character : taxon assigned at additional ranks (down to .rank)
  #
  # incorporates information from the closed-ref clustering into the
  # taxon table
  tar_fst_tbl(
    closedref_taxon_table,
    dplyr::left_join(
      known_taxon_table,
      dplyr::bind_rows(
        clusters_closed_ref_small,
        clusters_closed_ref_large
      ),
      by = "seq_id"
    ) |>
      dplyr::left_join(
        dplyr::select(known_taxon_table, ref_id = seq_id, cluster_taxon = .rank_sym),
        by = "ref_id"
      ) |>
      dplyr::mutate(
        .rank_sym := dplyr::coalesce(.rank_sym, cluster_taxon)
      ) |>
      dplyr::select(-ref_id, -cluster_taxon),
    deployment = "main"
  ),


  ##### predenovo_taxon_table_small_{.rank}_{.conf_level} #####
  # grouped tibble:
  #  `seq_id` character : unique ASV ID
  #  `seq_idx` integer : taxonomically sorted index of the sequence
  #  {ROOT_RANK} character : taxon assigned at ROOT_RANK (e.g. kingdom)
  #  ... character : taxon assigned at additional ranks (down to .rank)
  #  `tar_group` integer : grouping variable for dispatch to multiple dynamic
  #    jobs (matches values in .parent_rank)
  #
  # find only the ASVs which need to be de novo clustered
  # i.e. they are still unknown.
  # This version returns "small" taxa, i.e. those with too few sequences for
  # parallelization to be worthwhile, in batches to keep the total number of
  # targets smaller.
  tar_fst_tbl(
    predenovo_taxon_table_small,
    optimotu.pipeline::small_predenovo_taxon_table(
      closedref_taxon_table = closedref_taxon_table,
      asv_taxsort = asv_taxsort,
      rank = .rank,
      parent_rank = .parent_rank,
      tax_ranks <- !!optimotu.pipeline::tax_ranks()
    ),
    iteration = "group",
    deployment = "main"
  ),


  ##### predenovo_taxon_table_large_{.rank}_{.conf_level} #####
  # grouped tibble:
  #  `seq_id` character : unique ASV ID
  #  `seq_idx` integer : taxonomically sorted index of the sequence
  #  {ROOT_RANK} character : taxon assigned at ROOT_RANK (e.g. kingdom)
  #  ... character : taxon assigned at additional ranks (down to .rank)
  #  `tar_group` integer : grouping variable for dispatch to multiple dynamic
  #    jobs (matches values in .parent_rank)
  #
  # find only the ASVs which need to be de novo clustered
  # i.e. they are still unknown.
  # This version returns "large" taxa, i.e. those with many sequences so that
  # parallelization is worthwhile.
  tar_fst_tbl(
    predenovo_taxon_table_large,
    optimotu.pipeline::large_predenovo_taxon_table(
      closedref_taxon_table = closedref_taxon_table,
      asv_taxsort = asv_taxsort,
      rank = .rank,
      parent_rank = .parent_rank,
      tax_ranks <- !!optimotu.pipeline::tax_ranks()
    ),
    iteration = "group",
    deployment = "main"
  ),

  ##### denovo_thresholds_{.rank}_{.conf_level} #####
  # named list of named numeric
  #  list names: taxa at .parent_rank
  #  numeric values: optimal clustering thresholds
  #  numeric names: ranks from .rank to {TIP_RANK} (usually species)
  #
  # once we are doing de novo clustering at one rank, we will have to do it at
  # all subranks as well, and there is no new information that can change our
  # thresholds.  So it is easiest to just calculate the distance matrix once and
  # do all the clustering now.
  tar_target(
    denovo_thresholds,
    optimotu::calc_subtaxon_thresholds(
      rank = .parent_rank,
      taxon_table = dplyr::bind_rows(
        predenovo_taxon_table_small,
        predenovo_taxon_table_large
      ),
      optima = cluster_optima,
      ranks = !!optimotu.pipeline::tax_ranks(),
      measure = !!optimotu.pipeline::cluster_measure()
    ),
    deployment = "main"
  ),

  ##### clusters_denovo_small_{.rank}_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  {ROOT_RANK} character : taxon assigned at ROOT_RANK (e.g. kingdom)
  #  ... character : additional taxonomic assignments down to .parent_rank
  #  ... integer : unique cluster index, for ranks from .rank to TIP_RANK (usually species)
  tar_target(
    clusters_denovo_small,
    optimotu.pipeline::do_denovo_cluster(
      predenovo_taxon_table = predenovo_taxon_table_small,
      seq_file = !!seq_to_cluster_file,
      seq_file_index = !!seq_to_cluster_file_index,
      rank = .rank,
      parent_rank = .parent_rank,
      tax_ranks = !!optimotu.pipeline::tax_ranks(),
      denovo_thresholds = denovo_thresholds,
      dist_config = !!(
        if (optimotu.pipeline::cluster_dist_config()$method == "usearch") {
          quote(
            update(dc, usearch_ncpu = 1),
            list(dc = optimotu.pipeline::cluster_dist_config())
          ))
        } else {
          optimotu.pipeline::cluster_dist_config()
        }
      ),
      parallel_config = optimotu::parallel_concurrent(1)
    ),
    pattern = map(predenovo_taxon_table_small), # per taxon at .parent_rank
    resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
  ),

  ##### clusters_denovo_large_{.rank}_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  {ROOT_RANK} character : taxon assigned at ROOT_RANK (e.g. kingdom)
  #  ... character : additional taxonomic assignments down to .parent_rank
  #  ... integer : unique cluster index, for ranks from .rank to TIP_RANK (usually species)
  tar_target(
    clusters_denovo_large,
    optimotu.pipeline::do_denovo_cluster(
      predenovo_taxon_table = predenovo_taxon_table_large,
      seq_file = !!seq_to_cluster_file,
      seq_file_index = !!seq_to_cluster_file_index,
      rank = .rank,
      parent_rank = .parent_rank,
      tax_ranks = !!optimotu.pipeline::tax_ranks(),
      denovo_thresholds = denovo_thresholds,
      dist_config = !!(
        if (optimotu.pipeline::cluster_dist_config()$method == "usearch") {
          substitute(
            update(dc, usearch_ncpu = optimotu.pipeline::local_cpus()),
            list(dc = optimotu.pipeline::cluster_dist_config())
          ))
        } else {
          optimotu.pipeline::cluster_dist_config()
        }
      ),
      parallel_config = !!(
        if (optimotu.pipeline::cluster_dist_config()$method == "usearch") {
          quote(optimotu::parallel_concurrent(2))
        } else {
          quote(optimotu::parallel_concurrent(optimotu.pipeline::local_cpus()))
        }
      )
    ),
    pattern = map(predenovo_taxon_table_large), # per taxon at .parent_rank
    resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
  ),

  ##### taxon_table_{.rank}_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV ID
  #  {ROOT_RANK} character : taxon assigned at ROOT_RANK (e.g. kingdom)
  #  ... character : taxon assigned at additional ranks (down to .rank)
  #
  # taxonomy for all named taxa down to .rank
  tar_fst_tbl(
    taxon_table,
    dplyr::filter(closedref_taxon_table, !is.na(.rank_sym)),
    deployment = "main"
  ),

  ##### pseudotaxon_table_{.rank}_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  {ROOT_RANK} character : taxon assigned at ROOT_RANK (e.g. kingdom)
  #  ... character : additional taxonomic assignments down to .rank
  #  ... integer : unique cluster index, for ranks below .rank
  #
  # assign "pseudotaxon" names to de novo clusters at this .rank
  tar_fst_tbl(
    pseudotaxon_table,
    dplyr::bind_rows(
      clusters_denovo_small,
      clusters_denovo_large,
      .parent_pseudotaxa # results from previous rank
    ) |>
      dplyr::arrange(seq_id) |> # pseudotaxon numbers are ordered by ASV numbers
      dplyr::mutate(
        .rank_sym := paste(.parent_rank_sym, .rank_sym) |>
          forcats::fct_inorder() |>
          forcats::fct_relabel(
            ~names(optimotu.pipeline::name_seqs(., paste0("pseudo", .rank, "_")))
          ) |>
          as.character()
      ),
    deployment = "main"
  )
)

#### reliability_plan ####
# repeats entire clustering process for different assignment probability
# thresholds

reliability_meta <- c(
  plausible = 0.5,
  reliable = 0.9
) |>
  tibble::enframe(name = ".conf_level", value = ".prob_threshold")

reliability_plan <- tar_map(
  values = reliability_meta,
  names = .conf_level,

  ##### taxon_table_{ROOT_RANK}_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV ID
  #  {ROOT_RANK} character : taxon assigned at ROOT_RANK (e.g. kingdom)
  #
  # values for other ranks are calculated recursively
  # this should be everything, because PROTAX-fungi assigns all sequences
  # 100% probability of being fungi
  tar_target_raw(
    sprintf("taxon_table_%s", optimotu.pipeline::root_rank()),
    substitute(
      asv_tax_prob_reads |>
        dplyr::filter(rank == ROOT_RANK) |>
        dplyr::mutate(
          taxon = ifelse(prob < .prob_threshold, NA_character_, taxon)
        ) |>
        dplyr::select(seq_id, ROOT_RANK_VAR := taxon),
      list(
        ROOT_RANK = optimotu.pipeline::root_rank(),
        ROOT_RANK_VAR = optimotu.pipeline::root_rank_var()
      )
    ),
    format = "fst_tbl",
    deployment = "main"
  ),

  ##### pseudotaxon_table_{ROOT_RANK}_{.conf_level} #####
  # NULL
  #
  # this is required because pseudotaxon_table_{SECOND_RANK} will try to access its
  # parent, but there are no pseudotaxa at the root level, so it is empty.
  # they are combined with dplyr::bind_row(), which will be fine if we give it
  # NULL instead of a 0-row tibble with the correct columns.
  tar_target_raw(
    sprintf("pseudotaxon_table_%s", optimotu.pipeline::root_rank()),
    NULL,
    deployment = "main"
  ),

  rank_plan,

  ##### taxon_table_ingroup_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  {ROOT_RANK} character : taxonomic assignment at ROOT_RANK (e.g. kingdom)
  #  ... character : taxonomic assignments at all intermediate levels
  #  {TIP_RANK} character : taxonomic assignment at TIP_RANK (e.g. species)
  #
  # combine the taxon table and pseudotaxon table, and select only phyla
  # which a) are named ingroup phyla; or b) contain more known ingroup (based on
  # reference) than known outgroup or ASVs with no known ROOT_RANK.
  tar_fst_tbl(
    taxon_table_ingroup,
    dplyr::bind_rows(
      !!(taxon_table_TIP_RANK),
      !!(pseudotaxon_table_TIP_RANK)
    ) |>
      dplyr::mutate(
        known_outgroup = seq_id %in% asv_known_outgroup$seq_id,
        known_ingroup = seq_id %in% asv_known_ingroup$seq_id,
        unknown_outin = seq_id %in% asv_unknown_outin$seq_id
      ) |>
      dplyr::group_by(dplyr::across(!!optimotu.pipeline::second_rank_var())) |>
      dplyr::filter(
        !startsWith(!!optimotu.pipeline::second_rank_var(), "pseudo") |
          sum(known_ingroup) > sum(known_outgroup) + sum(unknown_outin)
      ) |>
      dplyr::select(!where(is.logical)) |>
      dplyr::arrange(seq_id),
    deployment = "main"
  ),

  ##### asv_otu_map_{.conf_level} #####
  tar_fst_tbl(
    asv_otu_map,
    dplyr::semi_join(
      taxon_table_ingroup,
      asv_table,
      by = "seq_id"
    ) |>
      dplyr::left_join(
        dplyr::select(otu_taxonomy, OTU = seq_id, !!optimotu.pipeline::tip_rank_var()),
        by = !!optimotu.pipeline::tip_rank()
      ) |>
      dplyr::select(ASV = seq_id, OTU),
    deployment = "main"
  ),

  ##### otu_taxonomy_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique OTU ID
  #  `ref_seq_id` character : unique ASV ID of the most prevalent/abundant ASV
  #    in this OTU
  #  `nsample` integer : the number of samples in which this OTU occurs
  #  `nread` integer : the number of reads of this OTU across all samples
  #  {ROOT_RANK} character : taxonomic assignment at ROOT_RANK (e.g. kingdom)
  #  ... character : taxonomic assignments at all intermediate levels
  #  {TIP_RANK} character : taxonomic assignment at TIP_RANK (e.g. species)
  tar_fst_tbl(
    otu_taxonomy,
    asv_table |>
      dplyr::mutate(asv_nsample = dplyr::n(), asv_nread = sum(nread), .by = seq_id) |>
      dplyr::inner_join(taxon_table_ingroup, by = "seq_id") |>
      dplyr::arrange(dplyr::desc(asv_nsample), dplyr::desc(asv_nread)) |>
      dplyr::summarize(
        nsample = as.integer(dplyr::n_distinct(sample)),
        nread = sum(nread),
        ref_seq_id = dplyr::first(seq_id),
        .by = c(!!!optimotu.pipeline::tax_rank_vars())
      ) |>
      dplyr::arrange(dplyr::desc(nsample), dplyr::desc(nread)) |>
      optimotu.pipeline::name_seqs("OTU", "seq_id") |>
      dplyr::select(seq_id, ref_seq_id, nsample, nread, everything()),
    deployment = "main"
  ),

  ##### otu_table_sparse_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique OTU id
  #  `sample` character : sample name
  #  `seqrun` character : name of sequencing run
  #  `nread` integer : number of reads
  #
  # OTU sample/abundance matrix, in sparse format (0's are not included)
  tar_fst_tbl(
    otu_table_sparse,
    asv_table |>
      dplyr::inner_join(taxon_table_ingroup, by = "seq_id") |>
      dplyr::inner_join(asv_otu_map, by = c("seq_id" = "ASV")) |>
      dplyr::group_by(OTU, sample, seqrun) |>
      dplyr::summarise(nread = sum(nread), .groups = "drop") |>
      dplyr::select(seq_id = OTU, sample, seqrun, nread),
    deployment = "main"
  )
)

#### clust_plan ####
# the outer plan which contains the reliability_plan

clust_plan <- c(
  list(
    ##### asv_known_outgroup #####
    # tibble:
    #  `seq_id` character : unique ASV id
    #  {KNOWN_RANKS} character : taxonomic assignment at KNOWN_RANKS (e.g. kingdom)
    #
    # ASVs whose best match is to a species of known outgroup
    asv_known_outgroup = tar_fst_tbl(
      asv_known_outgroup,
      {
        out <- asv_best_hit_taxon
        outgroup_cols <- character(length(!!optimotu.pipeline::known_ranks()))
        for (i in seq_along(!!optimotu.pipeline::known_ranks())) {
          rank_i = (!!optimotu.pipeline::known_ranks())[i]
          taxon_i <- (!!optimotu.pipeline::known_taxa())[i]
          outgroup_cols[i] <- paste0(rank_i, "_outgroup")
          out[[outgroup_cols[i]]] <-
            !is.na(out[[rank_i]]) &
            !out[[rank_i]] %in% c(taxon_i, "unspecified",
                                  "Eukaryota_kgd_Incertae_sedis", "None")
        }
        dplyr::filter(out, dplyr::if_any(all_of(outgroup_cols))) |>
          dplyr::select("seq_id", !!!optimotu.pipeline::known_ranks())
      },
      deployment = "main"
    ),

    ##### asv_known_ingroup #####
    # tibble:
    #  `seq_id` character : unique ASV id
    #  {KNOWN_RANKS} character : taxonomic assignment at ROOT_RANK (e.g. kingdom)
    #
    # ASVs whose best match is to an ingroup
    asv_known_ingroup = tar_fst_tbl(
      asv_known_ingroup,
      dplyr::filter(
        asv_best_hit_taxon,
        !!optimotu.pipeline::ingroup_rank_var() == !!optimotu.pipeline::ingroup_taxon()) |>
        dplyr::select(seq_id, !!!optimotu.pipeline::known_rank_vars()),
      deployment = "main"
    ),

    ##### asv_unknown_outin #####
    # tibble:
    #  `seq_id` character : unique ASV id
    #  {ROOT_RANK} character : taxonomic assignment at ROOT_RANK (e.g. kingdom)
    #
    # ASVs whose best match is to a species whose identity at ROOT_RANK is unknown
    asv_unknown_outin = tar_target(
      asv_unknown_outin,
      asv_best_hit_taxon |>
        dplyr::anti_join(asv_known_outgroup, by = "seq_id") |>
        dplyr::anti_join(asv_known_ingroup, by = "seq_id") |>
        dplyr::select(seq_id, !!!optimotu.pipeline::known_rank_vars()),
      deployment = "main"
    )
  ),

  reliability_plan
)

optimotu_plan <- c(optimotu_plan, clust_plan)

post_cluster_meta <-
  purrr::map_dfc(reliability_plan, tar_select_names, everything()) |>
  dplyr::mutate_all(rlang::syms) |>
  dplyr::bind_cols(reliability_meta)
