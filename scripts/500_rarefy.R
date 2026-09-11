# If we need to rarefy, wrap the whole plan in a map

if (optimotu.pipeline::do_rarefy()) {
  # Rarefaction is applied at sample-wise denoising. Targets that do not
  # depend on which reads survived rarefaction are pulled out so they run
  # once. That includes readwise trim, quality filter, and UNOISE pair
  # merging, plus constant files and most of phase 2.

  outside_rarefy_names <- intersect(
    c(
      # from 100_readwise.R / 102_runwise.R
      "readwise_meta",
      "readwise_meta_fwd",
      "readwise_meta_rev",
      "raw_R1",
      "raw_R1_fwd",
      "raw_R1_rev",
      "raw_R2",
      "raw_R2_fwd",
      "raw_R2_rev",
      "raw_read_counts",
      "raw_read_counts_fwd",
      "raw_read_counts_rev",
      "trim",
      "trim_fwd",
      "trim_rev",
      "trim_read_counts",
      "trim_read_counts_fwd",
      "trim_read_counts_rev",
      "filter_pairs",
      "filter_pairs_fwd",
      "filter_pairs_rev",
      "filt_read_counts",
      "filt_read_counts_fwd",
      "filt_read_counts_rev",
      "premerge_seqs",
      "premerge_seqs_fwd",
      "premerge_seqs_rev",
      "merge_read_counts",
      "merge_read_counts_fwd",
      "merge_read_counts_rev",
      "errfun",
      # from 103_phase1_combine.R
      "sample_table",
      "sample_table_key",
      "seq_all",
      # from 200_asv_filtering.R
      "seq_trim",
      "seq_index",
      "seqbatch",
      "seqbatch_hash",
      "unaligned_ref_seqs",
      "ref_chimeras",
      "spikes",
      "pos_controls",
      "amplicon_model_file",
      "amplicon_model_length",
      "amplicon_model_match",
      "seq_model_align",
      "seq_cm_align",
      "numts",
      "seq_full_length",
      "unaligned_ref_index",
      "outgroup_seqbatch",
      "outgroup_aligned",
      "outgroup_taxonomy",
      "best_hit",
      "best_hit_taxon",
      "best_hit_udb",
      # from 201_refseqs.R
      "taxonomy_default_file",
      "taxonomy_default",
      "taxonomy_ascii7_default_file",
      "taxonomy_ascii7_default",
      "new_refseq_file",
      "new_refseq",
      "new_refseq_metadata_file",
      "new_refseq_metadata",
      "taxonomy_new",
      "write_protax_taxonomy_new",
      "write_its2_new",
      "write_sintaxits2_new",
      "write_its2udb_new",
      "write_sintaxits2udb_new",
      "write_amptksynmockudb",
      "write_protax_taxonomy.ascii7_new",
      "write_protax_tax",
      "write_protax_ref.tax",
      "write_protax_rseqs",
      "custom_protax",
      "protax_model",
      # from 202_taxonomy.R
      "protax_dir",
      "all_tax_prob",
      "protax_script",
      "protax",
      "sintax_ref_file",
      "bayesant_ref_file",
      "bayesant_model",
      "epa_ref_file",
      "epa_taxonomy_file",
      "epa_tree_file",
      "epa_params",
      "epa_outgroup",
      "epa_ng",
      # from 203_supplemental_asv.R (not combo_* : those join native ASVs)
      "supp_sequences_file",
      "supp_set_info",
      "supp_asv_names",
      "supp_asv_seq",
      "supp_asv_seq_index",
      "supp_asv_seqbatch",
      "supp_asv_aligned_seq",
      "supp_sample_table_file",
      "supp_sample_table",
      "supp_taxonomy_file",
      "supp_best_hit_taxon",
      "supp_unknown_prob",
      "supp_tax_prob",
      "supp_all_tax_prob",
      # from 301_optimize_thresholds.R (file/reference training inputs)
      "threshold_train_file",
      "threshold_refseq_file",
      "threshold_refseq_index",
      # from 402_guilds.R
      optimotu.pipeline::guild_db_target_names(),
      # from 403_krona.R
      "krona_script",
      "krona_shortcut_icon",
      "krona_hiddenimage",
      "krona_loadingimage",
      "krona_logo"
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
        uniqs <- unique(
          !!optimotu.pipeline::tar_map_c(
            optimotu_plan[purrr::keep(
              names(optimotu_plan),
              startsWith,
              "seq_merged"
            )]
          )
        )
        seqs <- c(
          old_seqs,
          Biostrings::DNAStringSet(uniqs[is.na(BiocGenerics::match(
            uniqs,
            old_seqs
          ))])
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
