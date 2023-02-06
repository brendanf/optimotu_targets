## Filter unwanted sequences (chimeras, spikes) from DADA2 output
## and produce a final ASV table
## Brendan Furneaux

asv_plan <- list(
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
      new_batches <- tibble::tibble(seq = colnames(seqtable_nochim))
      if (file.exists(batches_file)) {
        batches <- fst::read_fst(batches_file)
        # check that we have not lost sequences; if we have we should scrap the
        # results and start over.
        if (all(batches$seq %in% new_batches$seq)) {
          maxbatch <- max(batches$tar_group)
          mean_batchsize <- nrow(batches) / maxbatch
          new_batches <- dplyr::anti_join(new_batches, batches, by = "seq_id")
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
              dplyr::mutate(seq_id = seq.int(dplyr::n())) |>
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
          dplyr::mutate(seq_id = seq.int(dplyr::n())) |>
          dplyr::ungroup()
        # new_batches$seq_id <- seqhash(new_batches$seq_id)
        fst::write_fst(new_batches, batches_file)
      }
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
      seq = colnames(seqtable_nochim),
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
    unname(seqtable_nochim)[,seqbatch_key$i, drop = FALSE],
    iteration = "list"
  ),
  
  #### ref_chimeras ####
  # Find reference-based chimeras in the current seqbatch.
  tar_fst_tbl(
    ref_chimeras,
    vsearch_uchime_ref(
      query = named_seqbatch,
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
    vsearch_usearch_global(
        dplyr::anti_join(named_seqbatch, ref_chimeras, by = "seq_id"),
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
      trim_primers(
        primer = "GCATCGATGAAGAACGCAGC...GCATATCAATAAGCGGAGGA",
        max_err = 0.2,
        min_overlap = 10
      )
  ),
  
  #### seqtable_dedup ####
  # Merge no-mismatch pairs
  tar_target(
    seqtable_dedup,
    collapseNoMismatch_vsearch(seqtable_nospike)
  ),
  
  #### write_asvtable ####
  tar_file(
    write_asvtable,
    file.path(asv_path, "asv_tab.rds") %T>%
      saveRDS(asvtable, .),
    deployment = "main"
  )
)