## Filter unwanted sequences (chimeras, spikes) from DADA2 output
## and produce a final ASV table
## Brendan Furneaux

asv_plan <- list(
  #### seqtable_dedup ####
  # Merge no-mismatch pairs
  tar_target(
    seqtable_dedup,
    collapseNoMismatch_vsearch(seqtable_nochim)
  ),
  
  #### seqbatch ####
  # Divide the unique ASV sequences from DADA2 into batches for further
  # processing.
  # In order to save time when the pipeline is re-run after new sequences are
  # added, cache and re-use the previous batch assignments.  Then sequences
  # which already existed 
  tar_fst_tbl(
    seqbatch,
    {
      batches_file <- "data/seqbatches.fst"
      new_batches <- tibble::tibble(seq = colnames(seqtable_dedup))
      if (file.exists(batches_file)) {
        batches <- fst::read_fst(batches_file)
        # check that we have not lost sequences; if we have we should scrap the
        # results and start over.
        if (all(batches$seq %in% new_batches$seq)) {
          maxbatch <- max(batches$tar_group)
          mean_batchsize <- nrow(batches) / maxbatch
          new_batches <- dplyr::anti_join(new_batches, batches, by = "seq")
          if (nrow(new_batches) > 0L) {
            # new_batches$seq_id <- seqhash(new_batches$seq)
            if (maxbatch < n_seqrun) {
              # if the cached number of batches is less than the number of seqruns
              # then assign the remaining sequences to the remaining seqruns
              new_batches$tar_group <- rep_len(seq.int(n_seqrun - maxbatch) + maxbatch, nrow(new_batches))
            } else {
              # otherwise add new batches of the same average size as the old ones
              n_batches <- min(round(nrow(new_batches) / mean_batchsize), n_seqrun)
              new_batches$tar_group <- rep_len(seq.int(n_batches) + maxbatch)
            }
            new_batches <- dplyr::group_by(new_batches, tar_group) |>
              dplyr::mutate(seq_id = as.character(seq.int(dplyr::n()))) |>
              dplyr::ungroup()
            new_batches <- rbind(batches, new_batches)
            fst::write_fst(new_batches, batches_file)
          } else {
            new_batches <- batches
          }
        }
      }
      if (!"tar_group" %in% names(new_batches)) {
        new_batches$tar_group <- rep_len(seq.int(n_seqrun), nrow(new_batches))
        new_batches <- dplyr::group_by(new_batches, tar_group) |>
          dplyr::mutate(seq_id = as.character(seq.int(dplyr::n()))) |>
          dplyr::ungroup()
        # new_batches$seq_id <- seqhash(new_batches$seq_id)
        fst::write_fst(new_batches, batches_file)
      }
      new_batches
    },
    iteration = "group"
  ),
  
  #### seqbatch_key ####
  # conversion between `i`, the column number in the seqtable (which may change
  # whenever new data is added) and `tar_group` the batch number plus `seq_id`
  # the within-batch index.
  tar_fst_tbl(
    seqbatch_key,
    tibble::tibble(
      seq = colnames(seqtable_dedup),
      i = seq_along(seq)
    ) |>
      dplyr::right_join(seqbatch, by = "seq") |>
      dplyr::select(-seq),
    iteration = "group"
  ),
  
  #### seqtable_batch ####
  # Slice of the seqtable containing only the sequences named in the current
  # batch.
  # The column names (full sequences) can be dropped to keep the size down
  tar_target(
    seqtable_batch,
    magrittr::set_colnames(seqtable_dedup, NULL)[,seqbatch_key$i, drop = FALSE],
    iteration = "list",
    pattern = map(seqbatch_key)
  ),
  
  #### ref_chimeras ####
  # Find reference-based chimeras in the current seqbatch.
  tar_fst_tbl(
    ref_chimeras,
    vsearch_uchime_ref(
      query = seqbatch,
      ref = "data/sh_matching_data/sanger_refs_sh.fasta",
      ncpu = local_cpus()
    ),
    pattern = map(seqbatch)
  ),
  
  #### nochim2_read_counts ####
  tar_target(
    nochim2_read_counts,
    tibble::enframe(
      rowSums(seqtable_batch[,-as.integer(ref_chimeras$seq_id), drop = FALSE]),
      name = "filt_key",
      value = "nochim2_nread"
    ),
    pattern = map(seqtable_batch, ref_chimeras)
  ),
  
  #### spikes ####
  # find spike sequences in the current seqbatch
  tar_target(
    spikes,
    seqbatch |>
      dplyr::anti_join(ref_chimeras, by = "seq_id") |>
      vsearch_usearch_global(
          "protaxFungi/addedmodel/amptk_synmock.udb",
          global = FALSE,
          threshold = 0.9
        ),
    pattern = map(seqbatch, ref_chimeras)
  ),
  
  #### nospike_read_counts ####
  tar_fst_tbl(
    nospike_read_counts,
    tibble::enframe(
      rowSums(
        seqtable_batch[,-as.integer(c(ref_chimeras$seq_id, spikes$seq_id)),
                       drop = FALSE]
      ),
      name = "filt_key",
      value = "nospike_nread"
    ),
    pattern = map(seqtable_batch, ref_chimeras, spikes)
  ),
  
  #### primer_trim ####
  tar_target(
    primer_trim,
    seqbatch |>
      dplyr::anti_join(ref_chimeras, by = "seq_id") |>
      dplyr::anti_join(spikes, by = "seq_id") |>
      trim_primer(
        primer = "GCATCGATGAAGAACGCAGC...GCATATCAATAAGCGGAGGA",
        max_err = 0.2,
        min_overlap = 10
      ),
    pattern = map(seqbatch, ref_chimeras, spikes),
    iteration = "list"
  ),
  
  #### unite_udb ####
  tar_file(
    unite_udb,
    build_filtered_udb(
      infile = "data/sh_matching_data/sanger_refs_sh.fasta",
      outfile = "sequences/filtered_sanger_refs_sh.udb",
      blacklist = c(
        "SH1154235.09FU", # chimeric; partial matches to two different fungi but labeled as a fern
        "SH1240531.09FU" # chimera of two fungi, labeled as a plant
      ),
      usearch = Sys.which("vsearch")
    )
  ),
  
  #### unite_match ####
  tar_fst_tbl(
    unite_match,
    vsearch_usearch_global(
      query = primer_trim,
      ref = unite_udb,
      threshold = 0.8,
      global = FALSE
    ),
    pattern = map(primer_trim),
    iteration = "list"
  ),
  
  #### asv_unite_kingdom ####
  tar_fst_tbl(
    asv_unite_kingdom,
    dplyr::mutate(seqbatch_key, seq_id = as.character(seq_id)) |>
      dplyr::group_split(tar_group, .keep = FALSE) |>
      purrr::map2(
        primer_trim,
        dplyr::semi_join,
        by = "seq_id"
      ) |>
      purrr::map2_dfr(
        unite_match,
        dplyr::full_join,
        by = "seq_id"
      ) |>
      dplyr::arrange(i) |>
      name_seqs("ASV", "seq_id") |>
      tidyr::separate(cluster, c("ref_id", "sh_id"), sep = "_") |>
      dplyr::left_join(
        readr::read_tsv(
          "data/sh_matching_data/shs_out.txt",
          col_names = c("sh_id", "taxonomy"),
          col_types = "cc-------"
        ),
        by = "sh_id"
      ) |>
      dplyr::transmute(
        seq_id = seq_id,
        kingdom = sub(";.*", "", taxonomy) |> substr(4, 100)
      )
  ),
  
  #### asv_table ####
  tar_fst_tbl(
    asv_table,
    dplyr::mutate(seqbatch_key, seq_id = as.character(seq_id)) |>
      dplyr::group_split(tar_group, .keep = FALSE) |>
      purrr::map2_dfr(
        primer_trim,
        dplyr::semi_join,
        by = "seq_id"
      ) |>
      dplyr::pull(i) |>
      sort() |>
      `[`(x = seqtable_dedup, ,j=_, drop = FALSE) |>
      name_seqs(prefix = "ASV") |>
      apply(2, dplyr::na_if, 0L) |>
      tibble::as_tibble(rownames = "filt_key") |>
      dplyr::left_join(sample_table[,c("seqrun", "sample", "filt_key")], by = "filt_key") |>
      dplyr::select(sample, seqrun, everything() & !filt_key) |>
      tidyr::pivot_longer(-(1:2), names_to = "seq_id", values_to = "nread", values_drop_na = TRUE) |>
      dplyr::arrange(seq_id, seqrun, sample),
    deployment = "main"
  ),
  
  #### asv_reads ####
  tar_fst_tbl(
    asv_reads,
    asv_table %>%
      dplyr::group_by(seq_id) %>%
      dplyr::summarize(nread = sum(nread)) %>%
      dplyr::semi_join(asv_tax, by = "seq_id"),
    deployment = "main"
  ),
  
  #### write_asvtable ####
  tar_file(
    write_asvtable,
    file.path(asv_path, "asv_tab.rds") %T>%
      saveRDS(asv_table, .),
    deployment = "main"
  ),
  
  #### asv_seq ####
  tar_fst_tbl(
    asv_seq,
    dplyr::mutate(seqbatch_key, seq_id = as.character(seq_id)) |>
      dplyr::group_split(tar_group, .keep = FALSE) |>
      purrr::map2_dfr(primer_trim, dplyr::inner_join, by = "seq_id") |>
      dplyr::arrange(i) |>
      name_seqs(prefix="ASV", id_col = "seq_id") |>
      dplyr::select(-i),
    deployment = "main"
  )
)
