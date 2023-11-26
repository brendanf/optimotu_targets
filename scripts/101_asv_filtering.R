## Filter unwanted sequences (chimeras, spikes) from DADA2 output
## and produce a final ASV table
## Brendan Furneaux

asv_plan <- list(
  #### seq_all_index ####
  # character filename
  # index file for fast access to sequences in seq_all
  tar_file_fast(
    seq_all_index,
    fastx_gz_index(seq_all)
  ),

  #### seqtable_dedup ####
  # dada2 sequence table; integer matrix of read counts with column names as
  # sequences and row names as "samples" (i.e. sample_table$filt_key)
  #

  # Discard ASVs with frame shifts or stop codons
  tar_target(
    seqtable_numtFilt,
    numts_filter(seqtable_nochim)
  ),

  # # Keep ASVs with invertebrate genetic code (5) [wrapping ORFfinder]
  # tar_target(
  #   seqtable_ORF,
  #   ORFfinder_run(seqtable_numtFilt)
  # ),


  # Merge no-mismatch ASVs
  tar_target(
    seqtable_dedup,
    collapseNoMismatch_vsearch(seqtable_numtFilt) |>
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
      new_batchkey <- tibble::tibble(seq_idx = seq_len(sequence_size(seq_all)))
      if (file.exists(batches_file)) {
        old_batchkey <- fst::read_fst(batches_file)
        nbatch_old <- max(batches$tar_group)
        mean_old_batchsize <- nrow(old_batchkey) / nbatch_old
        new_batchkey <- dplyr::anti_join(new_batchkey, old_batchkey, by = "seq_idx")
        if (nrow(new_batchkey) > 0L) {
          # new_batches$seq_id <- seqhash(new_batches$seq)
          min_nbatch_new <- ceiling(nrow(new_batchkey) / max_batchsize)
          if (min_nbatch_new + nbatch_old < n_seqrun * jobs_per_seqrun) {
            # if the cached number of batches is less than the target number of
            # jobs, and we can fit all the new sequences in the target number of
            # jobs, then do that.
            # This is the normal case when adding new seqruns to the analysis
            # (if the seqruns are small enough that the batchsize is less than
            #  the maximum)
            nbatch_new <- n_seqrun * jobs_per_seqrun - nbatch_old
          } else {
            # otherwise add new batches of the same average size as the old ones
            nbatch_new <- ceiling(round(nrow(new_batchkey) / mean_old_batchsize))
          }
          new_batchkey$tar_group <- rep(
            seq_len(nbatch_new) + nbatch_old,
            each = ceiling(nrow(new_batchkey) / nbatch_new),
            length.out = nrow(new_batchkey)
          )
          new_batchkey <- dplyr::mutate(
            new_batchkey,
            seq_id = as.character(seq_len(dplyr::n())),
            .by = tar_group
          )
          new_batchkey <- rbind(old_batchkey, new_batches)
          fst::write_fst(new_batches, batches_file)
        } else {
          new_batchkey <- old_batchkey
        }
      }
      if (!"tar_group" %in% names(new_batchkey)) {
        nbatch_new <- ceiling(max(
          nrow(new_batchkey) / max_batchsize, # maximum size batches
          n_seqrun * jobs_per_seqrun # one batch per job
        ))
        new_batchkey$tar_group <- rep(
          seq_len(nbatch_new),
          each = ceiling(nrow(new_batchkey) / nbatch_new),
          length.out = nrow(new_batchkey)
        )
        new_batchkey <-
          dplyr::mutate(
            new_batchkey,
            seq_id = as.character(seq_len(dplyr::n())),
            .by = tar_group
          )
        fst::write_fst(new_batchkey, batches_file)
      }
      new_batchkey
    },
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

  #### amplicon_hmm_file ####
  tar_file_fast(
    amplicon_hmm_file,
    "protaxAnimal/refs.hmm"
  ),

  #### hmm_align ####
  tar_file_fast(
    hmm_align,
    fastx_gz_extract(
      infile = seq_all,
      index = seq_all_index,
      i = seqbatch$seq_idx,
      outfile = withr::local_tempfile(fileext=".fasta")
    ) |>
      fastx_split(n = local_cpus(), outroot = withr::local_tempfile()) |>
      hmmalign(
        hmm = amplicon_hmm_file,
        outfile = sprintf("sequences/05_aligned/batch%05i.fasta.gz", seqbatch$tar_group[1])
      ),
    pattern = map(seqbatch)
  ),

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

  tar_file_fast(
    write_spike_seqs,
    dplyr::transmute(
      spike_seqs,
      seq_id = glue::glue("{seq_id};{spike_id};nsample={nsample};nseqrun={nseqrun};nread={nread}"),
      seq
    ) |>
      write_sequence("output/spike_asvs.fasta")
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
      dplyr::left_join(sample_table[,c("seqrun", "sample", "filt_key")], by = "filt_key") |>
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

  #### write_asvtable ####
  # character: path + file name
  #
  # write the sparse ASV table to the output directory
  tar_file_fast(
    write_asvtable,
    file.path(asv_path, "asv_tab.rds") %T>%
      saveRDS(asv_table, .),
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
