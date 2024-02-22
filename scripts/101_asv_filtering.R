## Filter unwanted sequences (chimeras, spikes) from DADA2 output
## and produce a final ASV table
## Brendan Furneaux

if (trim_options$action == "trim") {
  seq_all_trim <- quote(seq_all)
} else {
  seq_all_trim <- quote(seq_trim)
}

seq_dedup_file <- "sequences/04_denoised/seq_all_dedup.fasta.gz"

asv_plan <- list(
  if (trim_options$action %in% c("retain", "lowercase", "none")) {
    #### seq_trim ####
    # character filename
    # all sequences in seq_all, after trimming primers
    # this may include duplicates
    tar_file_fast(
      seq_trim,
      cutadapt_filter_trim(
        seq_all,
        primer = trim_primer_merged,
        options = cutadapt_options(
          max_err = trim_options$max_err,
          min_overlap = trim_options$min_overlap,
          action = "trim",
          discard_untrimmed = FALSE
        ),
        ncpu = local_cpus(),
        trim = "sequences/04_denoised/seq_all_trim.fasta.gz"
      )
    )
  },

  #### duplicate_seqs ####
  # tibble::tibble:
  #  `query` integer: index of a (lower-priority) sequence which matches a
  #    higher-priority sequence
  #  `hit` integer: index of the higher-priority sequence which is matched
  tar_fst_tbl(
    duplicate_seqs,
    nomismatch_hits_vsearch(seqtable_merged, !!seq_all_trim)
  ),

  #### seqtable_dedup ####
  # `tibble`:
  #   `sample (character) - sample name as given in sample_table$sample_key
  #   `seq_idx` (integer) - index of a sequence in seq_all
  #   `nread` (integer) number of reads
  #
  # Merge no-mismatch pairs
  tar_fst_tbl(
    seqtable_dedup,
    deduplicate_seqtable(seqtable_merged, duplicate_seqs),
  ),

  #### seqs_dedup ####
  # `character` file name (fasta.gz)
  # deduplicated sequences
  tar_file_fast(
    seq_dedup,
    deduplicate_seqs(
      seqs = !!seq_all_trim,
      hits = duplicate_seqs,
      outfile = seq_dedup_file
    ),
    deployment = "main"
  ),

  #### denovo_chimeras_dedup ####
  # `integer` : seq_index in seqs_dedup for denovo chimeras
  tar_target(
    denovo_chimeras_dedup,
    unique(deduplicate_seq_idx(denovo_chimeras, duplicate_seqs)),
    deployment = "main"
  ),

  #### seq_index ####
  # character filename
  # index file for fast access to sequences in seq_dedup
  tar_file_fast(
    seq_index,
    fastx_gz_index(seq_dedup),
    deployment = "main"
  ),

  #### seqbatch ####
  # grouped tibble:
  #  `seq_idx` integer: index of this sequence in seq_dedup
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
      new_batchkey <- tibble::tibble(seq_idx = seq_len(sequence_size(seq_dedup)))
      if (
        file.exists(batches_file) &&
        nrow(old_batchkey <- fst::read_fst(batches_file)) <= nrow(new_batchkey)
      ) {
        old_batchkey <- fst::read_fst(batches_file)
        nbatch_old <- max(old_batchkey$tar_group)
        mean_old_batchsize <- nrow(old_batchkey) / nbatch_old
        new_batchkey <- dplyr::anti_join(new_batchkey, old_batchkey, by = "seq_idx")
        if (nrow(new_batchkey) > 0L) {
          # new_batches$seq_id <- seqhash(new_batches$seq)
          min_nbatch_new <- ceiling(nrow(new_batchkey) / max_batchsize)
          if (min_nbatch_new + nbatch_old < n_workers) {
            # if the cached number of batches is less than the target number of
            # workers, and we can fit all the new sequences in the target number
            # of workers, then do that. This is the normal case when adding new
            # seqruns to the analysis (if the seqruns are small enough that the
            # batchsize is less than the maximum)
            nbatch_new <- n_workers - nbatch_old
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
          new_batchkey <- rbind(old_batchkey, new_batchkey)
          fst::write_fst(new_batchkey, batches_file)
        } else {
          new_batchkey <- old_batchkey
        }
      }
      if (!"tar_group" %in% names(new_batchkey)) {
        nbatch_new <- ceiling(max(
          nrow(new_batchkey) / max_batchsize, # maximum size batches
          n_workers # one batch per worker
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

  #### seqbatch_hash ####
  tar_target(
    seqbatch_hash,
    fastx_gz_hash(
      infile = seq_dedup,
      index = seq_index,
      start = min(seqbatch$seq_idx),
      n = nrow(seqbatch)
    ),
    pattern = map(seqbatch)
  ),

  #### unaligned_ref_seqs ####
  # character filename
  # sequences to use as reference for uchime
  tar_file_fast(
    unaligned_ref_seqs,
    "data/sh_matching_data/sanger_refs_sh.fasta"
  ),

  #### ref_chimeras ####
  # `integer` vector: indices in `seq_all` which are reference-based chimeras
  #
  # Find reference-based chimeras in the current seqbatch.
  tar_target(
    ref_chimeras,
    vsearch_uchime_ref(
      query = fastx_gz_extract(
        infile = seq_dedup_file, # actual file not a dependency
        index = seq_index,
        i = seqbatch$seq_idx,
        outfile = withr::local_tempfile(fileext=".fasta.gz"),
        hash = seqbatch_hash
      ),
      ref = unaligned_ref_seqs,
      ncpu = local_cpus(),
      id_only = TRUE,
      id_is_int = TRUE
    ),
    pattern = map(seqbatch, seqbatch_hash) # per seqbatch
  ),

  #### nochim2_read_counts ####
  # tibble:
  #  `sample_key` character: as `sample_table$sample_key`
  #  `nochim2_nread` integer: number of sequences in the sample after second
  #    chimera filtering
  tar_target(
    nochim2_read_counts,
    dplyr::filter(
      seqtable_dedup,
      !seq_idx %in% denovo_chimeras_dedup,
      !seq_idx %in% ref_chimeras
    ) |>
      dplyr::summarize(nochim2_nread = sum(nread), .by = sample) |>
      dplyr::rename(sample_key = sample)
  ),

  #### spikes ####
  # tibble:
  #  `seq_idx` integer: index of sample in seq_dedup
  #  `cluster` character: name of matching spike sequence
  #
  # find spike sequences in the current seqbatch
  tar_fst_tbl(
    spikes,
    vsearch_usearch_global(
      fastx_gz_extract(
        infile = seq_dedup_file, # actual file not a dependency
        index = seq_index,
        i = seqbatch$seq_idx,
        outfile = withr::local_tempfile(fileext=".fasta.gz"),
        hash = seqbatch_hash
      ),
      "protaxFungi/addedmodel/amptk_synmock.udb",
      global = FALSE,
      threshold = 0.9,
      id_is_int = TRUE
    ),
    pattern = map(seqbatch, ref_chimeras) # per seqbatch
  ),

  #### nospike_read_counts ####
  # tibble:
  #  `sample_key` character: as `sample_table$sample_key`
  #  `nospike_nread` integer: number of sequences in the sample after spike
  #    removal.
  tar_fst_tbl(
    nospike_read_counts,
    dplyr::filter(
      seqtable_dedup,
      !seq_idx %in% denovo_chimeras_dedup,
      !seq_idx %in% ref_chimeras
    ) |>
      dplyr::anti_join(spikes, by = "seq_idx") |>
      dplyr::summarize(nospike_nread = sum(nread), .by = sample) |>
      dplyr::rename(sample_key = sample)
  ),

  #### amplicon_cm_file ####
  # `character`: file name of CM file
  tar_file_fast(
    amplicon_cm_file,
    "data/ITS3_ITS4.cm"
  ),

  #### amplicon_cm_match ####
  # `tibble`:
  #  `idx` integer: index of the sequence _in this query_
  #  `seq_id` character: index of the sequence in seq_dedup
  #  `match_len` integer: length of the CM match
  #  `cm_from` integer: first matching position in the CM
  #  `cm_to` integer: last matching postion in the CM
  #  `trunc` character: is this a truncated match (whole CM is not present)
  #  `bit_sc` numeric: bit score of the match; higher is better.
  #  `avg_pp` numeric: average per-nucleotide posterior probability
  #  `time_band_calc` numeric: time in seconds to do banding calculations
  #  `time_alignment` numeric: time in seconds to do alignment
  #  `time_total` numeric: total time in seconds
  #  `mem_mb` numeric: maximum memory used
  tar_fst_tbl(
    amplicon_cm_match,
    {
      sfile <- tempfile(fileext = ".dat")
      inferrnal::cmalign(
        amplicon_cm_file,
        fastx_gz_extract(
          infile = seq_dedup_file, # actual file not a dependency
          index = seq_index,
          i = seqbatch$seq_idx,
          outfile = withr::local_tempfile(fileext=".fasta"),
          hash = seqbatch_hash
        ),
        global = TRUE,
        notrunc = TRUE,
        cpu = local_cpus(),
        sfile = sfile
      )
      read_sfile(sfile)
    },
    pattern = map(seqbatch, seqbatch_hash)
  ),

  #### asv_full_length ####
  # `integer`: index of sequences in seqs_dedup which are full-length CM matches
  tar_target(
    asv_full_length,
    dplyr::filter(
      amplicon_cm_match,
      bit_sc > 50,
      cm_from < 5,
      cm_to > 310
    )$seq_id |>
      as.integer(),
    pattern = map(amplicon_cm_match)
  ),

  #### full_length_read_counts ####
  tar_fst_tbl(
    full_length_read_counts,
    dplyr::filter(
      seqtable_dedup,
      !seq_idx %in% denovo_chimeras_dedup,
      !seq_idx %in% ref_chimeras,
      seq_idx %in% asv_full_length
    ) |>
      dplyr::anti_join(spikes, by = "seq_idx") |>
      dplyr::summarize(full_length_nread = sum(nread), .by = sample) |>
      dplyr::rename(sample_key = sample)
  ),

  #### best_hit_udb ####
  # character: path and file name for udb of reference sequences
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

  #### best_hit_kingdom ####
  # tibble:
  #  `seq_idx` integer: index in seqtable_dedup
  #  `ref_id` character: reference sequence id of best hit
  #  `sh_id` character: species hypothesis of best hit
  #  `kingdom` character: kingdom of best hit
  tar_fst_tbl(
    best_hit_kingdom,
    vsearch_usearch_global(
      query = fastx_gz_extract(
        infile = seq_dedup_file, # actual file not a dependency
        index = seq_index,
        i = seqbatch$seq_idx,
        outfile = withr::local_tempfile(fileext=".fasta"),
        hash = seqbatch_hash
      ),
      ref = unite_udb,
      threshold = 0.8,
      global = FALSE,
      id_is_int = TRUE
    ) |>
      dplyr::arrange(seq_idx) |>
      tidyr::separate(cluster, c("ref_id", "sh_id"), sep = "_") |>
      dplyr::left_join(
        readr::read_tsv(
          "data/sh_matching_data/shs_out.txt",
          col_names = c("sh_id", "taxonomy"),
          col_types = "cc-------"
        ),
        by = "sh_id"
      ) |>
      dplyr::mutate(
        kingdom = sub(";.*", "", taxonomy) |> substr(4, 100),
        .keep = "unused"
      ),
    pattern = map(seqbatch, seqbatch_hash) # per seqbatch
  ),

  #### spike_table ####
  # tibble:
  #  `sample` character: sample name (as in sample_table$sample)
  #  `seqrun` character: sequencing run (as in sample_table$seqrun)
  #  `seq_id` character: unique spike ASV id, in format "Spike[0-9]+". numbers
  #    are 0-padded
  #  'seq_idx` integer: index of sequence in seqs_dedup
  #  `spike_id` character: name of best hit spike sequence
  #  `nread` integer: number of reads
  #
  # global table of ASVs which are predicted to be spikes
  # (this does include chimeras)
  tar_fst_tbl(
    spike_table,
    dplyr::arrange(spikes, seq_idx) |>
      name_seqs("Spike", "seq_id") |>
      dplyr::left_join(seqtable_dedup, by = "seq_idx") |>
      dplyr::rename(spike_id = cluster, sample_key = sample) |>
      dplyr::left_join(sample_table_key, by = "sample_key") |>
      dplyr::select(sample, seqrun, seq_id, seq_idx, spike_id, nread)
  ),

  #### asv_names ####
  # tibble:
  #  `seq_idx` integer : index of a sequence in seqs_dedup
  #  `seq_id` character : unique ASV id, in format "ASV[0-9]+". numbers are
  #    0-padded
  tar_fst_tbl(
    asv_names,
    tibble::tibble(
      seq_idx = unname(asv_full_length) |>
        setdiff(denovo_chimeras_dedup) |>
        setdiff(ref_chimeras) |>
        setdiff(spikes$seq_idx) |>
        sort()
    ) |>
      name_seqs("ASV", "seq_id")
  ),

  #### asv_table ####
  # tibble:
  #  `sample` character: sample name (as in sample_table$sample)
  #  `seqrun` character: sequencing run (as in sample_table$seqrun)
  #  `seq_id` character: unique ASV id, in format "ASV[0-9]+". numbers are
  #    0-padded
  #  `seq_idx` character: index of ASV sequence in seqs_dedup
  #  `nread` integer: number of reads
  tar_fst_tbl(
    asv_table,
    seqtable_dedup |>
      dplyr::filter(
        !seq_idx %in% denovo_chimeras_dedup,
        !seq_idx %in% ref_chimeras,
        !seq_idx %in% spikes$seq_idx,
        seq_idx %in% asv_full_length
      ) |>
      dplyr::left_join(asv_names, by = "seq_idx") |>
      dplyr::rename(sample_key = sample) |>
      dplyr::left_join(sample_table_key, by = "sample_key") |>
      dplyr::select(sample, seqrun, seq_id, seq_idx, nread),
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
  # `character` filename
  # sequence for each ASV
  tar_file_fast(
    asv_seq,
    fastx_gz_extract(
      seq_dedup,
      seq_index,
      asv_names$seq_idx,
      "sequences/04_denoised/asv.fasta.gz"
    ) |>
      name_seqs(prefix = "ASV"),
    deployment = "main"
  ),

  #### asv_seq_index ####
  # `character` filename
  # sequence for each ASV
  tar_file_fast(
    asv_seq_index,
    fastx_gz_index(asv_seq),
    deployment = "main"
  ),

  #### asv_best_hit_kingdom ####
  # `tibble`:
  #  `seq_id` character: unique ASV identifier
  #  `ref_id` character: reference sequence id of best hit
  #  `sh_id` character: species hypothesis of best hit
  #  `kingdom` character: kingdom of best hit
  tar_fst_tbl(
    asv_best_hit_kingdom,
    dplyr::left_join(
      asv_names,
      best_hit_kingdom,
      by = "seq_idx"
    ) |>
      dplyr::select(-seq_idx),
    deployment = "main"
  ),

  #### asv_taxsort ####
  # `tibble`:
  #  `seq_idx` integer: inded of sequence in asv_taxsort_seq
  #  `seq_idx_in` integer: index of sequence in asv_seq
  tar_fst_tbl(
    asv_taxsort,
    tibble::rowid_to_column(asv_tax, "seq_idx_in") |>
      dplyr::arrange(dplyr::across(all_of(c(TAXRANKS, "seq_id")))) |>
      tibble::rowid_to_column("seq_idx") |>
      dplyr::select(seq_idx, seq_idx_in),
    deployment = "main"
  ),

  #### asv_taxsort_seq ####
  # `character` filename
  # sequence for each ASV
  tar_file_fast(
    asv_taxsort_seq,
    fastx_gz_extract(
      asv_seq,
      asv_seq_index,
      i = asv_taxsort$seq_idx_in,
      "sequences/04_denoised/asv_taxsort.fasta.gz"
    ),
    deployment = "main"
  ),

  #### asv_taxsort_seq_index ####
  # `character` filename
  # sequence for each ASV
  tar_file_fast(
    asv_taxsort_seq_index,
    fastx_gz_index(asv_taxsort_seq),
    deployment = "main"
  ),

  #### seqbatch_result_map ####
  tar_fst_tbl(
    seqbatch_result_map,
    seqbatch |>
      dplyr::transmute(
        seq_idx,
        result = as.raw(
          0x10 * (!seq_idx %in% denovo_chimeras_dedup) +
            0x20 * (!seq_idx %in% ref_chimeras) +
            0x40 * (!seq_idx %in% spikes$seq_idx) +
            0x80 * (seq_idx %in% asv_full_length)
          )
      ),
    pattern = map(seqbatch, ref_chimeras, spikes, asv_full_length),
    deployment = "main"
  ),

  #### asv_map ####
  tar_fst_tbl(
    asv_map,
    tibble::tibble(
      seq_idx_in = seq_len(sequence_size(!!seq_all_trim)),
      seq_idx = deduplicate_seq_idx(seq_idx_in, duplicate_seqs, merge = FALSE)
    ) |>
      dplyr::left_join(seqbatch_result_map, by = "seq_idx") |>
      dplyr::left_join(asv_names, by = "seq_idx") |>
      dplyr::select(seq_idx = seq_idx_in, result, seq_id)
  )
)

optimotu_plan <- c(optimotu_plan, asv_plan)
