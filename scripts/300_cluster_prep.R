# Final preparation of sequences for clustering

# We have relatively fast random access lookup of sequences by ASV id,
# but in various stages it is still probably faster to sort sequences by
# taxonomy, so that similar sequences are grouped together.

precluster_plan <- list(
  #### final_asv_taxsort ####
  # tibble:
  #  `seq_id` character : uniques ASV id, as in final_asv_seq
  #  `seq_idx` integer : index of the sequence after taxonomic sorting
  final_asv_taxsort = tar_fst_tbl(
    final_asv_taxsort,
    dplyr::left_join(
      Biostrings::fasta.seqlengths(!!final_asv_seq) |>
        names() |>
        tibble::tibble(seq_id = _),
      tidyr::pivot_wider(
        (!!final_asv_tax_prob) |>
          dplyr::mutate(
            rank = optimotu.pipeline::rank2factor(
              rank,
              !!optimotu.pipeline::tax_ranks()
            )
          ),
        id_cols = "seq_id",
        names_from = "rank",
        names_expand = TRUE,
        values_from = "taxon"
      )
    ) |>
      dplyr::arrange(!!!optimotu.pipeline::tax_rank_vars()) |>
      tibble::rowid_to_column("seq_idx") |>
      dplyr::select(seq_idx, seq_id),
    deployment = "main"
  ),
  #### final_asv_taxsort_seq ####
  # `character` filename
  # Sequences for each ASV, sorted by assigned taxonomy.
  # The file is a FASTA file with gzip compression.
  # Sequences are not renamed after sorting; they have the same names as in
  # asv_seq and combo_seq.
  # Because we will be rewriting all of the sequences, there is no need to use
  # `fastqindexr`, we will have to load it all into memory anyway.
  final_asv_taxsort_seq = tar_file(
    final_asv_taxsort_seq,
    Biostrings::readBStringSet(!!final_asv_seq)[final_asv_taxsort$seq_id] |>
      optimotu.pipeline::write_sequence(
        fname = file.path(
          !!optimotu.pipeline::asv_path(),
          !!(if (optimotu.pipeline::do_rarefy()) {
            quote(sprintf("asv_taxsort_%s.fasta.gz", .rarefy_text))
          } else {
            "asv_taxsort.fasta.gz"
          })
        ),
        compress = TRUE
      ),
    # use wide controller for memory
    resources = tar_resources(
      crew = tar_resources_crew(controller = "wide")
    )
  ),
  #### final_asv_taxsort_seq_index ####
  # fastqindexr_index object
  # Index for fast access to sequences in final_asv_taxsort_seq using
  # `fastqindexr`
  final_asv_taxsort_seq_index = tar_target(
    final_asv_taxsort_seq_index,
    fastqindexr::create_index(final_asv_taxsort_seq),
    deployment = "main"
  )
)

optimotu_plan <- c(optimotu_plan, precluster_plan)
