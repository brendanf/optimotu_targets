# If we need to rarefy, wrap the whole plan in a map

if (optimotu.pipeline::do_rarefy()) {

  # some targets do not need to be re-run for each rarefaction
  # these are constant values and much of phase 2.

  # TODO: trim and filter also do not need to be rarefied

  outside_rarefy_names <- intersect(
    c(
      # from 100_dada.R
      "readwise_meta", "readwise_meta_fwd", "readwise_meta_rev",
      "raw_R1", "raw_R1_fwd", "raw_R1_rev",
      "raw_R2", "raw_R2_fwd", "raw_R2_rev",
      "raw_read_counts", "raw_read_counts_fwd", "raw_read_counts_rev",
      "trim", "trim_fwd", "trim_rev",
      "trim_read_counts", "trim_read_counts_fwd", "trim_read_counts_rev",
      "filter_pairs", "filter_pairs_fwd", "filter_pairs_rev",
      "filt_read_counts", "filt_read_counts_fwd", "filt_read_counts_rev",
      "sample_table",
      "seq_all",
      "errfun",
      # from 101_asv_filtering.R
      "seq_trim", "seq_index", "seqbatch", "seqbatch_hash",
      "unaligned_ref_seqs", "ref_chimeras", "spikes", "pos_controls",
      "amplicon_model_file", "amplicon_model_length", "amplicon_model_match",
      "asv_model_align", "asv_cm_align", "numts", "asv_full_length",
      "unaligned_ref_index", "outgroup_seqbatch", "outgroup_aligned",
      "outgroup_taxonomy", "best_hit", "best_hit_taxon", "best_hit_udb",
      # from 102_protax_refseqs.R
      "taxonomy_default_file", "taxonomy_default",
      "taxonomy_ascii7_default_file", "taxonomy_ascii7_default",
      "new_refseq_file", "new_refseq", "new_refseq_metadata_file",
      "new_refseq_metadata", "taxonomy_new", "write_protax_taxonomy_new",
      "write_its2_new", "write_sintaxits2_new", "write_its2udb_new",
      "write_sintaxits2udb_new", "write_amptksynmockudb2",
      "write_protax_taxonomy.ascii7_new", "write_protax_tax",
      "write_protax_ref.tax", "write_protax_rseqs", "custom_protax",
      "protax_model",
      # from 103_taxonomy.R
      "protax_dir", "all_tax_prob", "protax_script", "protax", "sintax_ref_file",
      "bayesant_ref_file", "bayesant_model", "epa_ref_file", "epa_taxonomy_file",
      "epa_tree_file", "epa_params", "epa_outgroup", "epa_ng",
      # from 105_cluster.R
      # from 106_focus_taxa.R
      # from 200_output.R
      # from 201_guilds.R
      "funguild_db", "lifestyle_db_file", "lifestyle_db",
      # from 202_krona.R
      "krona_script", "krona_shortcut_icon", "krona_hiddenimage",
      "krona_loadingimage", "krona_logo"
    ),
    names(optimotu_plan)
  )
  outside_rarefy <- optimotu_plan[outside_rarefy_names]
  optimotu_plan[outside_rarefy_names] <- NULL
  optimotu_plan <- c(
    outside_rarefy,
    tar_map(
      values = optimotu.pipeline::rarefy_meta(),
      names = .rarefy_text,
      optimotu_plan
    )
  )

  # seq_all must be modified (actually just re-initialized) to include all
  # sequences from all rarefactions.
  optimotu_plan$seq_all <-
    tar_file(
      seq_all,
      {
        old_seqs <-
          if (file.exists(seq_all_file)) {
            Biostrings::readDNAStringSet(seq_all_file)
          } else {
            Biostrings::DNAStringSet()
          }
        uniqs <- unique(!!optimotu.pipeline::tar_map_c(
          optimotu_plan[purrr::keep(names(optimotu_plan), startsWith, "seq_merged")]
        ))
        seqs <- c(
          old_seqs,
          Biostrings::DNAStringSet(uniqs[is.na(BiocGenerics::match(uniqs, old_seqs))])
        )
        names(seqs) <- seq_along(seqs)
        optimotu.pipeline::write_and_return_file(
          seqs,
          file = seq_all_file,
          compress = "gzip",
          compression_level = 9
        )
      },
      resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
    )
}
