#### phase1_plan ####
phase1_plan <- c(
  list(
    ##### sample_table #####
    # `tibble`:
    #   `sample_key` character: unique identifier for the sample
    #   `sample` character: sample name
    #   `orient` character: orientation of the reads
    #   `seqrun` character: seqrun name
    #   `fastq_R1` character: path to the forward read file
    #   `fastq_R2` character: path to the reverse read file
    #   `trim_R1` character: path to the trimmed forward read file
    #   `trim_R2` character: path to the trimmed reverse read file
    #   `filt_R1` character: path to the filtered forward read file
    #   `filt_R2` character: path to the filtered reverse read file
    #   As well as optional columns to pass sample-specific trimming and filtering parameters
    #
    # The sample table will already have been calculated and cached in `optimotu.pipeline`
    # before it is called here, because it is required for static branching, etc.
    # However, the cached version is not available on workers, and it is better to
    # have it saved as a formal target for dependency tracking.  `sample_table_hash()` is used
    # to ensure that the target is re-run if the sample table has changed.
    sample_table = tar_target(
      sample_table,
      optimotu.pipeline::sample_table(!!optimotu.pipeline::sample_table_hash()),
      deployment = "main" # sample table not available on workers
    ),

    ##### sample_table_key #####
    # `tibble`:
    #   `sample_key` character: unique identifier for the sample
    #   `sample` character: sample name
    #   `seqrun` character: seqrun name
    sample_table_key = tar_target(
      sample_table_key,
      dplyr::select(sample_table, sample_key, sample, seqrun) |>
        dplyr::distinct(),
      deployment = "main"
    )
  ),

  seqrun_plan,

  list(
    ##### seq_all #####
    # `character` vector
    #
    # all unique ASV sequences, across all seqruns
    seq_all = tar_file(
      seq_all,
      {
        old_seqs <-
          if (file.exists(seq_all_file)) {
            Biostrings::readDNAStringSet(seq_all_file)
          } else {
            Biostrings::DNAStringSet()
          }
        uniqs <- unique(!!optimotu.pipeline::tar_map_c(seqrun_plan$seq_merged))
        seqs <- c(
          old_seqs,
          Biostrings::DNAStringSet(
            uniqs[is.na(BiocGenerics::match(uniqs, old_seqs))]
          )
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
    ),

    ##### denovo_chimeras #####
    # `integer` vector - index of ASV sequences which were found to be chimeric
    #
    # calculate consensus chimera calls across all seqruns
    denovo_chimeras = tar_target(
      denovo_chimeras,
      optimotu.pipeline::combine_bimera_denovo_tables(
        !!optimotu.pipeline::tar_map_bind_rows(seqrun_plan$bimera_table)
      ),
      resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
    ),

    ##### seqtable_merged #####
    # `tibble`:
    #   `sample (character) - sample name as given in sample_table$sample_key
    #   `seq_idx` (integer) - index of a sequence in seq_all
    #   `nread` (integer) number of reads
    seqtable_merged = tar_fst_tbl(
      seqtable_merged,
      !!optimotu.pipeline::tar_map_bind_rows(
        seqrun_plan,
        seqtable_final_name
      ),
      resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
    )
  ),

  if (isTRUE(optimotu.pipeline::do_lulu())) {
    list(
      ##### lulu_asv_map_{.rarefaction?}_{.replicate?} #####
      # tibble:
      #  `seq_idx` integer: index of sequence in seq_all
      #  `lulu_idx` integer: index of the denoised "parent" sequence in seq_all
      lulu_asv_map = tar_fst_tbl(
        lulu_asv_map,
        optimotu.pipeline::lulu_map(
          !!optimotu.pipeline::tar_map_bind_rows(
            seqrun_plan,
            "seqtable_raw"
          ),
          match_table = (!!optimotu.pipeline::tar_map_bind_rows(
            seqrun_plan,
            "lulu_match_table"
          )) |>
            dplyr::filter(
              dist <= !!optimotu.pipeline::lulu_max_dist(),
              n_gap <=
                !!(if (optimotu.pipeline::lulu_max_gap_total() >= 1) {
                  optimotu.pipeline::lulu_max_gap_total()
                } else {
                  substitute(
                    m * align_length,
                    list(m = optimotu.pipeline::lulu_max_gap_total())
                  )
                }),
              max_gap <=
                !!(if (optimotu.pipeline::lulu_max_gap_length() >= 1) {
                  optimotu.pipeline::lulu_max_gap_length()
                } else {
                  substitute(
                    m * align_length,
                    list(m = optimotu.pipeline::lulu_max_gap_length())
                  )
                })
            ),
          max_dist = !!optimotu.pipeline::lulu_max_dist(),
          min_abundance_ratio = !!optimotu.pipeline::lulu_min_abundance_ratio(),
          min_cooccurrence_ratio = !!optimotu.pipeline::lulu_min_cooccurrence_ratio(),
          use_mean_abundance_ratio = !!optimotu.pipeline::lulu_use_mean_abundance_ratio(),
          id_is_sorted = FALSE
        ),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      )
    )
  },
  list(
    ##### nochim1_read_counts #####
    # tibble:
    #  `sample_key` character: as `sample_table$sample_key`
    #  `nochim1_nread` integer: number of sequences in the sample after first
    #    chimera filtering
    nochim1_read_counts = tar_fst_tbl(
      nochim1_read_counts,
      dplyr::filter(seqtable_merged, !seq_idx %in% denovo_chimeras) |>
        dplyr::summarize(nochim1_nread = sum(nread), .by = sample) |>
        dplyr::rename(sample_key = sample),
      resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
    )
  )
)

optimotu_plan <- c(optimotu_plan, phase1_plan)
