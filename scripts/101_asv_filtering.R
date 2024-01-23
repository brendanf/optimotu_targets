## Filter unwanted sequences (chimeras, spikes) from DADA2 output
## and produce a final ASV table
## Brendan Furneaux

asv_plan <- list(
  #### seqtable_dedup ####
  # dada2 sequence table; integer matrix of read counts with column names as
  # sequences and row names as "samples" (i.e. sample_table$filt_key)
  #
  # Merge no-mismatch pairs
  tar_target(
    seqtable_dedup,
    collapseNoMismatch_vsearch(seqtable_nochim) |>
      sort_seq_table()
  ),

  #### seqbatch ####
  # grouped tibble:
  #  `seq` character: sequence
  #  `tar_group` integer: which batch this sequence is assigned to
  #  `seq_id` character: within-batch index of this sequence
  #
  # Divide the unique ASV sequences from DADA2 into batches for further
  # processing.
  # In order to save time when the pipeline is re-run after new sequences are
  # added, cache and re-use the previous batch assignments.  Then sequences
  # which already existed do not need to be re-run
  #
  # seqbatches.fst: data.frame with same columns as seqbatch
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
  # grouped tibble:
  #  `i` integer: column number in seqtable_dedup (which may change when new
  #    data is added)
  #  `tar_group` integer: the batch number
  #  `seq_id` character: the within-batch index.
  tar_fst_tbl(
    seqbatch_key,
    tibble::tibble(
      seq = colnames(seqtable_dedup),
      i = seq_along(seq)
    ) |>
      dplyr::right_join(seqbatch, by = "seq") |>
      dplyr::arrange(as.numeric(seq_id)) |>
      dplyr::select(-seq),
    iteration = "group"
  ),

  #### seqtable_batch ####
  # modified dada2 sequence table; integer matrix of read counts with no column
  # names and row names as "samples" (i.e. sample_table$filt_key)
  #
  # Slice of the seqtable containing only the sequences named in the current
  # batch.
  # The column names (full sequences) are dropped to keep the size down
  tar_target(
    seqtable_batch,
    magrittr::set_colnames(seqtable_dedup, NULL)[,seqbatch_key$i, drop = FALSE],
    iteration = "list",
    pattern = map(seqbatch_key) # per seqbatch
  ),

  #### ref_chimeras ####
  # tibble:
  #  `seq_id` character: within-batch index
  #  `seq` character: sequence
  #
  # Find reference-based chimeras in the current seqbatch.
  tar_fst_tbl(
    ref_chimeras,
    vsearch_uchime_ref(
      query = seqbatch,
      ref = "data/sh_matching_data/sanger_refs_sh.fasta",
      ncpu = local_cpus()
    ),
    pattern = map(seqbatch) # per seqbatch
  ),

  #### nochim2_read_counts ####
  # tibble:
  #  `filt_key` character: as `sample_table$filt_key`
  #  `nochim2_nread` integer: number of sequences in the sample after second
  #    chimera filtering
  tar_target(
    nochim2_read_counts,
    tibble::enframe(
      rowSums(seqtable_batch[,-as.integer(ref_chimeras$seq_id), drop = FALSE]),
      name = "filt_key",
      value = "nochim2_nread"
    ),
    pattern = map(seqtable_batch, ref_chimeras) # per seqbatch
  ),

  #### spikes ####
  # tibble:
  #  `seq_id` character: within-batch index
  #  `cluster` character: name of matching spike sequence
  #
  # find spike sequences in the current seqbatch
  tar_fst_tbl(
    spikes,
    seqbatch |>
      dplyr::anti_join(ref_chimeras, by = "seq_id") |>
      vsearch_usearch_global(
          "protaxFungi/addedmodel/amptk_synmock.udb",
          global = FALSE,
          threshold = 0.9
        ),
    pattern = map(seqbatch, ref_chimeras) # per seqbatch
  ),

  #### nospike_read_counts ####
  # tibble:
  #  `filt_key` character: as `sample_table$filt_key`
  #  `nospike_nread` integer: number of sequences in the sample after spike
  #    removal.
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
    pattern = map(seqtable_batch, ref_chimeras, spikes) # per seqrun
  ),

  #### primer_trim ####
  # tibble:
  #  `seq_id` character: within-batch index
  #  `seq` character: trimmed sequence
  #
  # If primers were already removed in the initial step, then this is not needed now
  if (trim_options$action %in% c("retain", "lowercase", "none")) {
    tar_fst_tbl(
      primer_trim,
      seqbatch |>
        dplyr::anti_join(ref_chimeras, by = "seq_id") |>
        dplyr::anti_join(spikes, by = "seq_id") |>
        trim_primer(
          primer = trim_primer_merged,
          cutadapt_options(max_err = 0.2, min_overlap = 10)
        ),
      pattern = map(seqbatch, ref_chimeras, spikes), # per seqbatch
      iteration = "list"
    )
  } else {
    tar_fst_tbl(
      primer_trim,
      seqbatch |>
        dplyr::anti_join(ref_chimeras, by = "seq_id") |>
        dplyr::anti_join(spikes, by = "seq_id"),
      pattern = map(seqbatch, ref_chimeras, spikes), # per seqbatch
      iteration = "list"
    )
  },

  #### amplicon_cm_file ####
  tar_file_fast(
    amplicon_cm_file,
    "data/ITS3_ITS4.cm"
  ),

  #### amplicon_cm_match ####
  tar_fst_tbl(
    amplicon_cm_match,
    {
      sfile <- tempfile(fileext = ".dat")
      inferrnal::cmalign(
        amplicon_cm_file,
        tibble::deframe(primer_trim),
        global = TRUE,
        notrunc = TRUE,
        cpu = local_cpus(),
        sfile = sfile
      )
      read_sfile(sfile)
    },
    pattern = map(primer_trim),
    iteration = "list"
  ),

  #### asv_full_length ####
  tar_fst_tbl(
    asv_full_length,
    dplyr::filter(
      amplicon_cm_match,
      bit_sc > 50,
      cm_from < 5,
      cm_to > 310
    ) |>
      dplyr::semi_join(
        primer_trim,
        y = _,
        by = "seq_id"
      ),
    pattern = map(primer_trim, amplicon_cm_match),
    iteration = "list"
  ),

  #### full_length_read_counts ####
  tar_fst_tbl(
    full_length_read_counts,
    tibble::enframe(
      rowSums(
        seqtable_batch[,as.integer(asv_full_length$seq_id),
                       drop = FALSE]
      ),
      name = "filt_key",
      value = "full_length_nread"
    ),
    pattern = map(seqtable_batch, asv_full_length) # per seqrun
  ),

  #### unite_udb ####
  # character: path and file name for udb of Unite sanger reference sequences
  #
  # build a udb index for fast vsearch
  tar_file_fast(
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
  # tibble:
  #  `seq_id` character: within batch index
  #  `cluster` character: name of best Unite match
  tar_fst_tbl(
    unite_match,
    vsearch_usearch_global(
      query = asv_full_length,
      ref = unite_udb,
      threshold = 0.8,
      global = FALSE
    ),
    pattern = map(asv_full_length), # per seqbatch
    iteration = "list"
  ),

  #### asv_unite_kingdom ####
  # tibble:
  #  `seq_id` character: within batch index
  #  `kingdom` character: kingdom of best Unite match
  #
  # combine seqbatches and look up the kingdom for the best Unite matches
  tar_fst_tbl(
    asv_unite_kingdom,
    dplyr::mutate(seqbatch_key, seq_id = as.character(seq_id)) |>
      dplyr::group_split(tar_group, .keep = FALSE) |>
      purrr::map2(
        asv_full_length,
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

  #### all_spikes ####
  # tibble:
  #  `i` integer: column index in seqtable_dedup
  #  `spike_id` character: name of best hit spike sequence
  #
  # map spike lists back to indices in the seqtable and combine them
  tar_fst_tbl(
    all_spikes,
    dplyr::left_join(spikes, seqbatch_key, by = "seq_id") |>
      dplyr::select(i, spike_id = cluster),
    pattern = map(spikes, seqbatch_key)
  ),

  #### spike_table ####
  # tibble:
  #  `sample` character: sample name (as in sample_table$sample)
  #  `seqrun` character: sequencing run (as in sample_table$seqrun)
  #  `seq_id` character: unique spike ASV id, in format "Spike[0-9]+". numbers
  #    are 0-padded
  #  `spike_id` character: name of best hit spike sequence
  #  `nread` integer: number of reads
  #
  # global table of ASVs which are predicted to be spikes
  tar_fst_tbl(
    spike_table,
    seqtable_dedup[,all_spikes$i] |>
      `colnames<-`(all_spikes$i) |>
      apply(2, dplyr::na_if, 0L) |>
      tibble::as_tibble(rownames = "filt_key") |>
      tidyr::pivot_longer(
        -1,
        names_to = "i",
        names_transform = as.integer,
        values_to = "nread",
        values_drop_na = TRUE
      ) |>
      dplyr::left_join(sample_table, by = "filt_key") |>
      dplyr::summarize(nread = sum(nread), .by = c(sample, seqrun, i)) |>
      dplyr::left_join(
        dplyr::arrange(all_spikes, i)
        |> name_seqs("Spike", "seq_id"),
        by = "i"
      ) |>
      dplyr::select(sample, seqrun, seq_id, spike_id, nread)
  ),

  tar_fst_tbl(
    spike_seqs,
    dplyr::arrange(all_spikes, i) |>
      name_seqs("Spike", "seq_id") |>
      dplyr::mutate(seq = colnames(seqtable_dedup)[i]) |>
      dplyr::left_join(
        dplyr::summarize(spike_table, nread = sum(nread), nsample = dplyr::n_distinct(sample), nseqrun = dplyr::n_distinct(seqrun), .by = seq_id),
        by = "seq_id") |>
      dplyr::arrange(seq_id)
  ),

  #### asv_table ####
  # tibble:
  #  `sample` character: sample name (as in sample_table$sample)
  #  `seqrun` character: sequencing run (as in sample_table$seqrun)
  #  `seq_id` character: unique ASV id, in format "ASV[0-9]+". numbers are
  #    0-padded
  #  `nread` integer: number of reads
  #
  # combine batches to form a sparse global ASV table
  tar_fst_tbl(
    asv_table,
    dplyr::mutate(seqbatch_key, seq_id = as.character(seq_id)) |>
      dplyr::group_split(tar_group, .keep = FALSE) |>
      purrr::map2_dfr(
        asv_full_length,
        dplyr::semi_join,
        by = "seq_id"
      ) |>
      dplyr::pull(i) |>
      sort() |>
      `[`(x = seqtable_dedup, ,j=_, drop = FALSE) |>
      name_seqs(prefix = "ASV") |>
      apply(2, dplyr::na_if, 0L) |>
      tibble::as_tibble(rownames = "filt_key") |>
      dplyr::left_join(
        unique(sample_table[,c("seqrun", "sample", "filt_key")]),
        by = "filt_key"
      ) |>
      dplyr::select(sample, seqrun, everything() & !filt_key) |>
      tidyr::pivot_longer(-(1:2), names_to = "seq_id", values_to = "nread", values_drop_na = TRUE) |>
      dplyr::arrange(seq_id, seqrun, sample),
    deployment = "main"
  ),

  #### asv_reads ####
  # tibble:
  #  `seq_id` character: unique ASV id
  #  `nread` integer: total reads across all samples
  #
  # calculate total read counts for all ASVs (at least those present in asv_tax)
  tar_fst_tbl(
    asv_reads,
    asv_table %>%
      dplyr::group_by(seq_id) %>%
      dplyr::summarize(nread = sum(nread)) %>%
      dplyr::semi_join(asv_tax, by = "seq_id"),
    deployment = "main"
  ),

  #### asv_seq ####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  `seq` character: sequence
  tar_fst_tbl(
    asv_seq,
    dplyr::mutate(seqbatch_key, seq_id = as.character(seq_id)) |>
      dplyr::group_split(tar_group, .keep = FALSE) |>
      purrr::map2_dfr(asv_full_length, dplyr::inner_join, by = "seq_id") |>
      dplyr::arrange(i) |>
      name_seqs(prefix="ASV", id_col = "seq_id") |>
      dplyr::select(-i),
    deployment = "main"
  ),

  #### seqbatch_result_map ####
  tar_fst_tbl(
    seqbatch_result_map,
    seqbatch_key |>
      dplyr::left_join(
        dplyr::transmute(ref_chimeras, seq_id, nochim2 = FALSE),
        by = "seq_id"
      ) |>
      dplyr::left_join(
        dplyr::transmute(spikes, seq_id, nonspike = FALSE),
        by = "seq_id"
      ) |>
      dplyr::left_join(
        dplyr::transmute(asv_full_length, seq_id, model_match = TRUE),
        by = "seq_id"
      ) |>
      dplyr::mutate(
        result = as.raw(
          0x0F * dplyr::coalesce(nochim2, TRUE) +
            0x10 * dplyr::coalesce(nonspike, nochim2, TRUE) +
            0x20 * dplyr::coalesce(model_match, FALSE)
        )
      ) |>
      dplyr::select(i, result),
    pattern = map(seqbatch_key, ref_chimeras, spikes, asv_full_length),
    deployment = "main"
  ),

  #### asv_map ####
  tar_fst_tbl(
    asv_map,
    dplyr::left_join(
      seqbatch_result_map,
      dplyr::filter(seqbatch_result_map, result == 0x1f) |>
        dplyr::arrange(i) |>
        name_seqs(prefix = "ASV"),
      by = c("i", "result")
    ) |>
      dplyr::left_join(
        attr(seqtable_dedup, "map"),
        y = _,
        by = c("seq_id_out" = "i")
      ) |>
      dplyr::select(seq_id = seq_id_in, result, ASV)
  )
)
