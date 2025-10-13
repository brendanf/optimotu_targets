## Filter unwanted sequences (chimeras, spikes) from DADA2 output
## and produce a final ASV table
## Brendan Furneaux

seq_trim_file <- file.path(optimotu.pipeline::asv_path(), "seq_all_trim.fasta.gz")
if (optimotu.pipeline::trim_options()$action == "trim") {
  seq_all_trim <- quote(seq_all)
  seq_all_trim_file <- seq_all_file
} else {
  seq_all_trim <- quote(seq_trim)
  seq_all_trim_file <- seq_trim_file
}
seq_index_file <- paste0(seq_all_trim_file, ".fqi")

asv_plan <- c(
  list(
    seq_trim = if (optimotu.pipeline::trim_options()$action %in% c("retain", "lowercase", "none")) {
      #### seq_trim ####
      # character filename
      # all sequences in seq_all, after trimming primers
      # this may include duplicates
      tar_file(
        seq_trim,
        optimotu.pipeline::cutadapt_filter_trim(
          seq_all,
          primer = !!optimotu.pipeline::trim_primer_merged(),
          options = optimotu.pipeline::cutadapt_options(
            max_err = !!(optimotu.pipeline::trim_options()$max_err),
            min_overlap = !!(optimotu.pipeline::trim_options()$min_overlap),
            action = "trim",
            discard_untrimmed = FALSE
          ),
          ncpu = optimotu.pipeline::local_cpus(),
          trim = !!seq_trim_file
        ),
        resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
      )
    },

    #### seq_index ####
    # character filename
    # index file for fast access to sequences in seq_all_trim
    seq_index = tar_file(
      seq_index,
      optimotu.pipeline::fastx_gz_index(!!seq_all_trim),
      deployment = "main"
    ),

    #### seqbatch ####
    # grouped tibble:
    #  `seq_idx` integer: index of this sequence in seq_all_trim
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
    seqbatch = tar_fst_tbl(
      seqbatch,
      {
        batches_file <- optimotu.pipeline::ensure_directory("data/seqbatches.fst")
        new_batchkey <- tibble::tibble(seq_idx = seq_len(optimotu.pipeline::sequence_size(!!seq_all_trim)))
        if (
          file.exists(batches_file) &&
          nrow(old_batchkey <- fst::read_fst(batches_file)) <= nrow(new_batchkey)
        ) {
          old_batchkey <- fst::read_fst(batches_file)
          nbatch_old <- max(old_batchkey$tar_group)
          mean_old_batchsize <- nrow(old_batchkey) / nbatch_old
          new_batchkey <- dplyr::anti_join(new_batchkey, old_batchkey, by = "seq_idx")
          if (nrow(new_batchkey) > 0L) {
            # new_batches$seq_id <- optimotu.pipeline::seqhash(new_batches$seq)
            min_nbatch_new <- ceiling(nrow(new_batchkey) / !!optimotu.pipeline::max_batchsize())
            if (min_nbatch_new + nbatch_old < !!optimotu.pipeline::n_workers()) {
              # if the cached number of batches is less than the target number of
              # workers, and we can fit all the new sequences in the target number
              # of workers, then do that. This is the normal case when adding new
              # seqruns to the analysis (if the seqruns are small enough that the
              # batchsize is less than the maximum)
              nbatch_new <- !!optimotu.pipeline::n_workers() - nbatch_old
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
            nrow(new_batchkey) / !!optimotu.pipeline::max_batchsize(), # maximum size batches
            !!optimotu.pipeline::n_workers() # one batch per worker
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
      iteration = "group",
      resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
    ),

    #### seqbatch_hash ####
    seqbatch_hash = tar_target(
      seqbatch_hash,
      optimotu.pipeline::fastx_gz_hash(
        infile = !!seq_all_trim,
        index = seq_index,
        start = min(seqbatch$seq_idx),
        n = nrow(seqbatch)
      ),
      pattern = map(seqbatch),
      resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
    ),

    #### unaligned_ref_seqs ####
    # character filename
    # sequences to use as reference for uchime
    unaligned_ref_seqs = tar_file(
      unaligned_ref_seqs,
      !!optimotu.pipeline::outgroup_reference(),
      resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
    ),

    #### ref_chimeras ####
    # `integer` vector: indices in `seq_all` which are reference-based chimeras
    #
    # Find reference-based chimeras in the current seqbatch.
    ref_chimeras = tar_target(
      ref_chimeras,
      withr::with_tempfile(
        "outfile",
        fileext = ".fasta.gz",
        optimotu.pipeline::vsearch_uchime_ref(
          query = optimotu.pipeline::fastx_gz_extract(
            infile = seq_all_trim_file, # actual file not a dependency
            index = seq_index_file, # actual file not a dependency
            i = seqbatch$seq_idx,
            outfile = outfile,
            hash = seqbatch_hash # this is where the dependency is tracked
          ),
          ref = unaligned_ref_seqs,
          ncpu = optimotu.pipeline::local_cpus(),
          id_only = TRUE,
          id_is_int = TRUE
        )
      ),
      pattern = map(seqbatch, seqbatch_hash), # per seqbatch
      resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
    ),

    #### nochim2_read_counts ####
    # tibble:
    #  `sample_key` character: as `sample_table$sample_key`
    #  `nochim2_nread` integer: number of sequences in the sample after second
    #    chimera filtering
    nochim2_read_counts = tar_target(
      nochim2_read_counts,
      dplyr::filter(
        seqtable_merged,
        !seq_idx %in% denovo_chimeras,
        !seq_idx %in% ref_chimeras
      ) |>
        dplyr::summarize(nochim2_nread = sum(nread), .by = sample) |>
        dplyr::rename(sample_key = sample),
      resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
    )
  ),
  if (optimotu.pipeline::do_spike()) {
    list(
      #### spikes ####
      # tibble:
      #  `seq_idx` integer: index of sequence in seq_all_trim
      #  `cluster` character: name of matching spike sequence
      #
      # find spike sequences in the current seqbatch
      spikes = tar_fst_tbl(
        spikes,
        withr::with_tempfile(
          "outfile",
          fileext = ".fasta.gz",
          optimotu.pipeline::vsearch_usearch_global(
            optimotu.pipeline::fastx_gz_extract(
              infile = seq_all_trim_file, # actual file not a dependency
              index = seq_index,
              i = seqbatch$seq_idx,
              outfile = outfile,
              hash = seqbatch_hash
            ),
            !!optimotu.pipeline::spike_file(),
            global = FALSE,
            threshold = 0.9,
            id_is_int = TRUE
          )
        ),
        pattern = map(seqbatch, ref_chimeras), # per seqbatch
        resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
      ),

      #### nospike_read_counts ####
      # tibble:
      #  `sample_key` character: as `sample_table$sample_key`
      #  `nospike_nread` integer: number of sequences in the sample after spike
      #    removal.
      nospike_read_counts = tar_fst_tbl(
        nospike_read_counts,
        dplyr::filter(
          seqtable_merged,
          !seq_idx %in% denovo_chimeras,
          !seq_idx %in% ref_chimeras
        ) |>
          dplyr::anti_join(spikes, by = "seq_idx") |>
          dplyr::summarize(nospike_nread = sum(nread), .by = sample) |>
          dplyr::rename(sample_key = sample),
        resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
      ),

      #### spike_read_counts ####
      # tibble:
      #  `sample_key` character: as `sample_table$sample_key`
      #  `spike_nread` integer: number of sequences in the sample which were spikes
      spike_read_counts = tar_fst_tbl(
        spike_read_counts,
        dplyr::filter(
          seqtable_merged,
          seq_idx %in% spikes$seq_idx,
          !seq_idx %in% denovo_chimeras,
          !seq_idx %in% ref_chimeras
        ) |>
          dplyr::summarize(spike_nread = sum(nread), .by = sample) |>
          dplyr::rename(sample_key = sample),
        resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
      )
    )
  },

  if (optimotu.pipeline::do_pos_control()) {
    list(
      #### positive controls ####
      # tibble:
      #  `seq_idx` integer: index of sequence in seq_all_trim
      #  `cluster` character: name of matching spike sequence
      #
      # find positive control sequences in the current seqbatch
      pos_controls = tar_fst_tbl(
        pos_controls,
        withr::with_tempfile(
          "outfile",
          fileext = ".fasta.gz",
          optimotu.pipeline::vsearch_usearch_global(
            optimotu.pipeline::fastx_gz_extract(
              infile = seq_all_trim_file, # actual file not a dependency
              index = seq_index,
              i = seqbatch$seq_idx,
              outfile = outfile,
              hash = seqbatch_hash
            ),
            !!optimotu.pipeline::pos_control_file(),
            global = FALSE,
            threshold = 0.9,
            id_is_int = TRUE
          )
        ),
        pattern = map(seqbatch, ref_chimeras), # per seqbatch
        resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
      ),

      #### nocontrol_read_counts ####
      # tibble:
      #  `sample_key` character: as `sample_table$sample_key`
      #  `nocontrol_nread` integer: number of sequences in the sample after
      #    positive control removal.
      nocontrol_read_counts = tar_fst_tbl(
        nocontrol_read_counts,
        dplyr::filter(
          seqtable_merged,
          !seq_idx %in% denovo_chimeras,
          !seq_idx %in% ref_chimeras,
          !! (
            if (optimotu.pipeline::do_spike()) {
              quote(!seq_idx %in% spikes$seq_idx)
            } else {
              TRUE
            }
          )
        ) |>
          dplyr::anti_join(pos_controls, by = "seq_idx") |>
          dplyr::summarize(nocontrol_nread = sum(nread), .by = sample) |>
          dplyr::rename(sample_key = sample),
        resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
      ),

      #### control_read_counts ####
      # tibble:
      #  `sample_key` character: as `sample_table$sample_key`
      #  `control_nread` integer: number of sequences in the sample which were
      #    positive controls
      control_read_counts = tar_fst_tbl(
        control_read_counts,
        dplyr::filter(
          seqtable_merged,
          seq_idx %in% pos_controls$seq_idx,
          !seq_idx %in% denovo_chimeras,
          !seq_idx %in% ref_chimeras,
          !! (
            if (optimotu.pipeline::do_spike()) {
              quote(!seq_idx %in% spikes$seq_idx)
            } else {
              TRUE
            }
          )
        ) |>
          dplyr::summarize(control_nread = sum(nread), .by = sample) |>
          dplyr::rename(sample_key = sample),
        resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
      )
    )
  },

  #### Amplicon models ####
  if (!identical(optimotu.pipeline::amplicon_model_type(), "none")) {
    c(
      list(
        ##### amplicon_model_file #####
        # `character`: file name of CM or HMM file
        amplicon_model_file = tar_file(
          amplicon_model_file,
          optimotu.pipeline::amplicon_model_file(),
          deployment = "main"
        )
      ),
      if (identical(optimotu.pipeline::amplicon_model_type(), "CM")) {
        c(
          ##### CM #####
          list(
            ###### amplicon_model_length ######
            amplicon_model_length = tar_target(
              amplicon_model_length,
              readLines(amplicon_model_file, n = 100) |>
                purrr::keep(startsWith, "CLEN") |>
                readr::parse_number() |>
                as.integer(),
              deployment = "main"
            )
          ),

          if (optimotu.pipeline::do_model_filter_only()) {
            list(
              ###### do_model_filter_only ######
              ####### amplicon_model_match #######
              # `tibble`:
              #  `seq_idx` character: index of the sequence in seq_all_trim
              #  `model_from` integer: first matching position in the CM
              #  `model_to` integer: last matching postion in the CM
              #  `bit_score` numeric: bit score of the match; higher is better.
              amplicon_model_match = tar_fst_tbl(
                amplicon_model_match,
                withr::with_tempfile(
                  "sfile",
                  fileext = ".dat",
                  {
                    withr::with_tempfile(
                      "outfile",
                      fileext = ".fasta",
                      inferrnal::cmalign(
                        amplicon_model_file,
                        optimotu.pipeline::fastx_gz_extract(
                          infile = seq_all_trim_file, # actual file not a dependency
                          index = seq_index,
                          i = seqbatch$seq_idx,
                          outfile = outfile,
                          hash = seqbatch_hash
                        ),
                        global = TRUE,
                        notrunc = TRUE,
                        cpu = optimotu.pipeline::local_cpus(),
                        sfile = sfile
                      )
                    )
                    optimotu.pipeline::read_sfile(sfile) |>
                      dplyr::transmute(
                        seq_idx = as.integer(seq_id),
                        model_from = cm_from,
                        model_to = cm_to,
                        bit_score = bit_sc
                      )
                  }
                ),
                pattern = map(seqbatch, seqbatch_hash),
                resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
              )
            )
          } else if (optimotu.pipeline::do_model_align_only()) {
            list(
              ###### do_model_align_only ######
              ####### asv_model_align #######
              asv_model_align = tar_file(
                asv_model_align,
                withr::with_tempfile(
                  "tempout",
                  fileext = ".fasta",
                  inferrnal::cmalign(
                    amplicon_model_file,
                    optimotu.pipeline::fastx_gz_extract(
                      infile = seq_all_trim_file, # actual file not a dependency
                      index = seq_index,
                      i = seqbatch$seq_idx,
                      outfile = tempout,
                      hash = seqbatch_hash
                    ),
                    global = TRUE,
                    notrunc = TRUE,
                    dnaout = TRUE,
                    cpu = optimotu.pipeline::local_cpus()
                  ) |>
                    optimotu.pipeline::consensus_columns() |>
                    optimotu.pipeline::write_sequence(
                      fname = file.path(
                        !!optimotu.pipeline::aligned_path(),
                        sprintf("batch%05i.fasta.gz", seqbatch$tar_group[1])
                      ),
                      compress = TRUE
                    )
                ),
                pattern = map(seqbatch, seqbatch_hash),
                resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
              )
            )
          } else if (do_model_both) {
            ###### do_model_both ######
            list(
              ####### asv_cm_align #######
              # `character`: two file names, for the alignment and the alignment
              # stats
              asv_cm_align = tar_file(
                asv_cm_align,
                withr::with_tempfile(
                  "tempout",
                  fileext = ".fasta",
                  {
                    sfile <- file.path(
                      !!optimotu.pipeline::aligned_path(),
                      sprintf("batch%05i.sfile", seqbatch$tar_group[1])
                    )
                    optimotu.pipeline::ensure_directory(sfile)
                    inferrnal::cmalign(
                      amplicon_model_file,
                      optimotu.pipeline::fastx_gz_extract(
                        infile = seq_all_trim_file, # actual file not a dependency
                        index = seq_index,
                        i = seqbatch$seq_idx,
                        outfile = tempout,
                        hash = seqbatch_hash
                      ),
                      global = TRUE,
                      notrunc = TRUE,
                      dnaout = TRUE,
                      cpu = optimotu.pipeline::local_cpus(),
                      sfile = sfile
                    ) |>
                      optimotu.pipeline::consensus_columns() |>
                      optimotu.pipeline::write_sequence(
                        fname = file.path(
                          !!optimotu.pipeline::aligned_path(),
                          sprintf("batch%05i.fasta.gz", seqbatch$tar_group[1])
                        ),
                        compress = TRUE
                      ) |>
                      c(sfile)
                  }
                ),
                pattern = map(seqbatch, seqbatch_hash),
                resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
              ),

              ####### asv_model_align #######
              asv_model_align = tar_file(
                asv_model_align,
                asv_cm_align[1],
                pattern = map(asv_cm_align),
                deployment = "main"
              ),

              ####### amplicon_model_match #######
              # `tibble`:
              #  `seq_idx` character: index of the sequence in seq_all_trim
              #  `model_from` integer: first matching position in the CM
              #  `model_to` integer: last matching postion in the CM
              #  `bit_score` numeric: bit score of the match; higher is better.
              amplicon_model_match = tar_fst_tbl(
                amplicon_model_match,
                optimotu.pipeline::read_sfile(asv_cm_align[2]) |>
                  dplyr::transmute(
                    seq_idx = as.integer(seq_id),
                    model_from = cm_from,
                    model_to = cm_to,
                    bit_score = bit_sc
                  ),
                pattern = map(asv_cm_align),
                deployment = "main"
              )
            )
          }
        )
      } else if (identical(optimotu.pipeline::amplicon_model_type(), "HMM")) {
        ##### HMM #####
        c(
          list(
            ###### amplicon_model_length ######
            amplicon_model_length = tar_target(
              amplicon_model_length,
              readLines(amplicon_model_file, n = 100) |>
                purrr::keep(startsWith, "LENG") |>
                readr::parse_number() |>
                as.integer(),
              deployment = "main"
            )
          ),

          if (optimotu.pipeline::do_model_filter()) {
            list(
              ###### amplicon_model_match ######
              # `tibble`:
              #  `seq_idx` (integer) : index of the sequence in seq_all_trim
              #  `hmm_from` (integer) : start position of the match in the HMM
              #  `hmm_to` (integer) : end position of the match in the HMM
              #  `bit_score` (numeric) : score for the match (higher is better)
              amplicon_model_match = tar_fst_tbl(
                amplicon_model_match,
                withr::with_tempfile(
                  "tempout",
                  fileext = ".fasta",
                  optimotu.pipeline::nhmmer(
                    seqs = optimotu.pipeline::fastx_gz_extract(
                      infile = seq_all_trim_file,
                      index = seq_index,
                      i = seqbatch$seq_idx,
                      outfile = tempout,
                      hash = seqbatch_hash
                    ),
                    hmm = amplicon_model_file
                  ) |>
                    dplyr::transmute(
                      seq_idx = as.integer(seq_name),
                      model_from = hmm_from,
                      model_to = hmm_to,
                      bit_score
                    )
                ),
                pattern = map(seqbatch, seqbatch_hash),
                resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
              )
            )
          },

          if (optimotu.pipeline::do_model_align()) {
            list(
              ###### asv_model_align ######
              asv_model_align = tar_file(
                asv_model_align,
                withr::with_tempfile(
                  "tempout",
                  fileext = ".fasta",
                  optimotu.pipeline::fastx_gz_extract(
                    infile = seq_all_trim_file,
                    index = seq_index,
                    i = seqbatch$seq_idx,
                    outfile = tempout,
                    hash = seqbatch_hash
                  ) |>
                    optimotu.pipeline::fastx_split(
                      n = optimotu.pipeline::local_cpus(),
                      outroot = optimotu.pipeline::ensure_directory(
                        tempfile(tmpdir = withr::local_tempdir())
                      ),
                      compress = TRUE
                    ) |>
                    optimotu.pipeline::hmmalign(
                      hmm = amplicon_model_file,
                      outfile = file.path(
                        !!optimotu.pipeline::aligned_path(),
                        sprintf("batch%05i.fasta.gz", seqbatch$tar_group[1])
                      )
                    )
                ),
                pattern = map(seqbatch, seqbatch_hash),
                resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
              )
            )
          },

          if (optimotu.pipeline::do_numt_filter()) {
            list(
              ###### numts ######
              tar_fst_tbl(
                numts,
                optimotu.pipeline::detect_numts(asv_model_align, id_is_int = TRUE),
                pattern = map(asv_model_align),
                deployment = "main"
              )
            )
          }
        )
      } else {
        stop("invalid value for amplicon_model_type: ",
             optimotu.pipeline::amplicon_model_type())
      },

      if (optimotu.pipeline::do_model_filter()) {
        list(
          ##### asv_full_length #####
          # `integer`: index of sequences in seqs_dedup which are full-length model matches
          asv_full_length = tar_target(
            asv_full_length,
            dplyr::filter(
              amplicon_model_match,
              bit_score >= !!optimotu.pipeline::min_model_score(),
              model_from <= !!optimotu.pipeline::max_model_start(),
              model_to >= !!optimotu.pipeline::min_model_end()
            )$seq_idx,
            pattern = map(amplicon_model_match),
            resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
          ),

          ##### full_length_read_counts #####
          full_length_read_counts = tar_fst_tbl(
            full_length_read_counts,
            dplyr::filter(
              seqtable_merged,
              !seq_idx %in% denovo_chimeras,
              !seq_idx %in% ref_chimeras,
              seq_idx %in% asv_full_length,
              !! (if (optimotu.pipeline::do_spike()) {
                quote(!seq_idx %in% spikes$seq_idx)
              } else {
                TRUE
              }),
              !! (if (optimotu.pipeline::do_pos_control()) {
                quote(!seq_idx %in% pos_controls$seq_idx)
              } else {
                TRUE
              })
            ) |>
              dplyr::summarize(full_length_nread = sum(nread), .by = sample) |>
              dplyr::rename(sample_key = sample),
            resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
          )
        )
      }
    )
  },

  if (optimotu.pipeline::do_model_align()) {
    #### aligned ####
    list(
      ##### unaligned_ref_index #####
      unaligned_ref_index =
        if (endsWith(optimotu.pipeline::outgroup_reference(), ".gz")) {
          tar_file(
            unaligned_ref_index,
            optimotu.pipeline::fastx_gz_index(unaligned_ref_seqs),
            resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
          )
        } else {
          tar_fst(
            unaligned_ref_index,
            Biostrings::fasta.index(unaligned_ref_seqs),
            resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
          )
        },

      ##### outgroup_seqbatch #####
      # tibble:
      #  `batch` integer : index of batch
      #  `batch_id` character : name of batch
      #  `from` integer : index in outgroup reference file to start this batch
      #  `to` integer : index in outgroup reference file to end this batch
      outgroup_seqbatch = tar_fst_tbl(
        outgroup_seqbatch,
        {
          n_seq <- optimotu.pipeline::sequence_size(unaligned_ref_seqs)
          n_batch <- as.integer(ceiling(n_seq / !!optimotu.pipeline::max_batchsize()))
          batchsize <- as.integer(ceiling(n_seq / n_batch))
          tibble::tibble(
            batch = seq_len(n_batch),
            batch_id = optimotu.pipeline::make_seq_names(n_batch, "ref"),
            from = (batch - 1L) * batchsize + 1L,
            to = dplyr::lead(from, 1L, default = n_seq + 1L) - 1L
          )
        },
        resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
      ),

      ##### outgroup_aligned #####
      # character: path and file name for fasta file representing a batch of
      #   aligned reference sequences
      outgroup_aligned = tar_file(
        outgroup_aligned,
        withr::with_tempfile(
          "tempout",
          fileext = ".fasta",
          {
            !!if (endsWith(optimotu.pipeline::outgroup_reference(), ".gz")) {
              quote(
                optimotu.pipeline::fastx_gz_extract(
                  infile = unaligned_ref_seqs,
                  index = unaligned_ref_index,
                  i = seq(outgroup_seqbatch$from, outgroup_seqbatch$to),
                  outfile = tempout
                )
              )
            } else {
              quote(
                Biostrings::readDNAStringSet(
                  unaligned_ref_index[with(outgroup_seqbatch, from:to),]
                ) |>
                  optimotu.pipeline::write_sequence(tempout, width=19999L)
              )
            }
            withr::with_tempfile(
              "outroot",
              optimotu.pipeline::fastx_split(
                tempout,
                n = optimotu.pipeline::local_cpus(),
                outroot = outroot,
                compress = FALSE
              ) |>
                optimotu.pipeline::hmmalign(
                  hmm = amplicon_model_file,
                  outfile = file.path(
                    !!optimotu.pipeline::aligned_path(),
                    sprintf("%s.fasta.gz", outgroup_seqbatch$batch_id)
                  )
                )
            )
          }
        ),
        pattern = map(outgroup_seqbatch),
        resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
      ),
      ##### outgroup_taxonomy #####
      # tibble:
      #  `ref_id` character: reference sequence id
      #  {TAX_RANKS} character: taxonomy of reference sequence
      outgroup_taxonomy = tar_fst_tbl(
        outgroup_taxonomy,
        names(Biostrings::fasta.seqlengths(unaligned_ref_seqs)) |>
          tibble::tibble(name = _) |>
          tidyr::separate_wider_delim(name, delim = "|", names_sep = "_") |>
          dplyr::select(ref_id = 1, taxonomy = last_col()) |>
          tidyr::separate(taxonomy, !!optimotu.pipeline::tax_ranks(), sep = ",", extra = "drop"),
        resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
      ),
      ##### best_hit #####
      # tibble:
      #  `seq_idx` integer: index in seq_all_trim
      #  `ref_id` character: reference sequence id of best hit
      #  `dist` numeric: distance to best hit
      best_hit = tar_fst_tbl(
        best_hit,
        optimotu::seq_search(
          query = asv_model_align,
          ref = outgroup_aligned,
          threshold = 0.5,
          dist_config = optimotu::dist_hamming(min_overlap = 300, ignore_gaps = FALSE),
          parallel_config = optimotu::parallel_concurrent(optimotu.pipeline::local_cpus())
        ) |>
          dplyr::mutate(
            seq_idx = as.integer(seq_id),
            .keep = "unused",
            .before = 1
          ),
        pattern = cross(asv_model_align, outgroup_aligned),
        resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
      ),
      ##### best_hit_taxon #####
      # tibble:
      #  `seq_idx` integer: index in seq_all_trim
      #  `ref_id` character: reference sequence id of best hit
      #  `dist` numeric: distance to best hit
      #  `sh_id` character: species hypothesis of best hit
      #  {INGROUP_RANK} character: taxon of best hit at {INGROUP_RANK} (e.g., kingdom)
      best_hit_taxon = tar_fst_tbl(
        best_hit_taxon,
        dplyr::slice_min(best_hit, dist, by = seq_idx, with_ties = FALSE) |>
          dplyr::left_join(outgroup_taxonomy, by = "ref_id"),
        deployment = "main"
      )
    )
  } else {
    #### unaligned ####
    list(
      ##### best_hit_udb #####
      # character: path and file name for udb of reference sequences
      #
      # build a udb index for fast vsearch
      best_hit_udb = tar_file(
        best_hit_udb,
        optimotu.pipeline::build_filtered_udb(
          infile = unaligned_ref_seqs,
          outfile = !!file.path(optimotu.pipeline::seq_path(), "outgroup_reference.udb"),
          blacklist = c(
            "SH1154235.09FU", # chimeric; partial matches to two different fungi but labeled as a fern
            "SH1240531.09FU" # chimera of two fungi, labeled as a plant
          ),
          usearch = Sys.which("vsearch")
        ),
        resources = tar_resources(crew = tar_resources_crew(controller = "wide")) # memory
      ),

      ##### best_hit_taxon #####
      # tibble:
      #  `seq_idx` integer: index in seq_all_trim
      #  `ref_id` character: reference sequence id of best hit
      #  `dist` numeric: distance to best hit
      #  `sh_id` character: species hypothesis of best hit
      #  {KNOWN_RANKS} character: taxon of best hit at {KNOWN_RANKS} (e.g., kingdom)
      best_hit_taxon = if (is.null(optimotu.pipeline::outgroup_taxonomy())) {
        tar_fst_tbl(
          best_hit_taxon,
          withr::with_tempfile(
            "tempout",
            fileext = ".fasta",
            optimotu.pipeline::vsearch_usearch_global(
              query = optimotu.pipeline::fastx_gz_extract(
                infile = seq_all_trim_file, # actual file not a dependency
                index = seq_index,
                i = seqbatch$seq_idx,
                outfile = tempout,
                hash = seqbatch_hash
              ),
              ref = best_hit_udb,
              threshold = 0.8,
              global = FALSE,
              id_is_int = TRUE
            ) |>
              dplyr::arrange(seq_idx) |>
              tidyr::separate(cluster, c("ref_id", "sh_id", "taxonomy"), sep = "[|]") |>
              tidyr::separate(taxonomy, !!optimotu.pipeline::tax_ranks(), sep = ",", fill = "right")
          ),
          pattern = map(seqbatch, seqbatch_hash), # per seqbatch
          resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
        )
      } else {
        tar_fst_tbl(
          best_hit_taxon,
          withr::with_tempfile(
            "tempout",
            fileext = ".fasta",
            optimotu.pipeline::vsearch_usearch_global(
              query = optimotu.pipeline::fastx_gz_extract(
                infile = seq_all_trim_file, # actual file not a dependency
                index = seq_index,
                i = seqbatch$seq_idx,
                outfile = tempout,
                hash = seqbatch_hash
              ),
              ref = best_hit_udb,
              threshold = 0.8,
              global = FALSE,
              id_is_int = TRUE
            ) |>
              dplyr::arrange(seq_idx) |>
              tidyr::separate(cluster, c("ref_id", "sh_id"), sep = "_") |>
              dplyr::left_join(
                readr::read_tsv(
                  !!optimotu.pipeline::outgroup_taxonomy(),
                  col_names = c("sh_id", "taxonomy"),
                  col_types = "cc-------"
                ),
                by = "sh_id"
              ) |>
              dplyr::mutate(
                !!optimotu.pipeline::ingroup_rank_var() := sub(";.*", "", taxonomy)
                |> substr(4, 100),
                .keep = "unused"
              )
          ),
          pattern = map(seqbatch, seqbatch_hash), # per seqbatch
          resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
        )
      }
    )
  },

  if (optimotu.pipeline::do_spike()) {
    list(
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
      spike_table = tar_fst_tbl(
        spike_table,
        dplyr::arrange(spikes, seq_idx) |>
          optimotu.pipeline::name_seqs("Spike", "seq_id") |>
          dplyr::left_join(seqtable_merged, by = "seq_idx") |>
          dplyr::rename(spike_id = cluster, sample_key = sample) |>
          dplyr::left_join(sample_table_key, by = "sample_key") |>
          dplyr::select(sample, seqrun, seq_id, seq_idx, spike_id, nread),
        resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
      )
    )
  },

  if (optimotu.pipeline::do_pos_control()) {
    list(
      #### control_table ####
      # tibble:
      #  `sample` character: sample name (as in sample_table$sample)
      #  `seqrun` character: sequencing run (as in sample_table$seqrun)
      #  `seq_id` character: unique control ASV id, in format "Control[0-9]+". numbers
      #    are 0-padded
      #  'seq_idx` integer: index of sequence in seqs_dedup
      #  `control_id` character: name of best hit positive control sequence
      #  `nread` integer: number of reads
      #
      # global table of ASVs which are predicted to be control sequences
      # (this does include chimeras)
      control_table = tar_fst_tbl(
        control_table,
        dplyr::arrange(pos_controls, seq_idx) |>
          optimotu.pipeline::name_seqs("Control", "seq_id") |>
          dplyr::left_join(seqtable_merged, by = "seq_idx") |>
          dplyr::rename(control_id = cluster, sample_key = sample) |>
          dplyr::left_join(sample_table_key, by = "sample_key") |>
          dplyr::select(sample, seqrun, seq_id, seq_idx, control_id, nread),
        resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
      )
    )
  },

  list(
    #### asv_names ####
    # tibble:
    #  `seq_idx` integer : index of a sequence in seqs_dedup
    #  `seq_id` character : unique ASV id, in format "ASV[0-9]+". numbers are
    #    0-padded
    asv_names = tar_fst_tbl(
      asv_names,
      tibble::tibble(
        seq_idx =
          (!!(
            if (optimotu.pipeline::do_model_filter())
              quote(unname(asv_full_length))
            else
              quote(seq_len(optimotu.pipeline::sequence_size(!!seq_all_trim)))
          )) |>
          setdiff(denovo_chimeras) |>
          setdiff(ref_chimeras) |>
          setdiff(!!(
            if (optimotu.pipeline::do_spike())
              quote(spikes$seq_idx)
            else
              integer()
          )) |>
          setdiff(!!(
            if (optimotu.pipeline::do_pos_control())
              quote(pos_controls$seq_idx)
            else
              integer()
          )) |>
          setdiff(!!(
            if (optimotu.pipeline::do_numt_filter())
              quote(numts$seq_idx)
            else
              integer()
          )) |>
          intersect(seqtable_merged$seq_idx) |>
          sort()
      ) |>
        optimotu.pipeline::name_seqs("ASV", "seq_id"),
      deployment = "main"
    ),

    #### asv_table ####
    # tibble:
    #  `sample` character: sample name (as in sample_table$sample)
    #  `seqrun` character: sequencing run (as in sample_table$seqrun)
    #  `seq_id` character: unique ASV id, in format "ASV[0-9]+". numbers are
    #    0-padded
    #  `seq_idx` character: index of ASV sequence in seqs_dedup
    #  `nread` integer: number of reads
    asv_table = tar_fst_tbl(
      asv_table,
      seqtable_merged |>
        dplyr::filter(
          !seq_idx %in% denovo_chimeras,
          !seq_idx %in% ref_chimeras,
          !!(if (optimotu.pipeline::do_spike())
            quote(!seq_idx %in% spikes$seq_idx) else TRUE),
          !!(if (optimotu.pipeline::do_model_filter())
            quote(seq_idx %in% asv_full_length) else TRUE),
          !!(if (optimotu.pipeline::do_numt_filter())
            quote(!seq_idx %in% numts$seq_idx) else TRUE)
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
    asv_reads = tar_fst_tbl(
      asv_reads,
      asv_table |>
        dplyr::group_by(seq_id) |>
        dplyr::summarize(nread = sum(nread)) |>
        dplyr::semi_join(asv_tax, by = "seq_id"),
      deployment = "main"
    ),

    #### asv_seq ####
    # `character` filename
    # sequence for each ASV
    asv_seq = tar_file(
      asv_seq,
      optimotu.pipeline::write_sequence(
        Biostrings::readDNAStringSet(!!seq_all_trim)[asv_names$seq_idx] |>
          optimotu.pipeline::name_seqs(prefix = "ASV"),
        file.path(
          !!optimotu.pipeline::asv_path(),
          !!(if(optimotu.pipeline::do_rarefy()) {
            quote(sprintf("asv_%s.fasta.gz", .rarefy_text))
          } else {
            "asv.fasta.gz"
          })
        ),
        compress = TRUE
      ),
      deployment = "main"
    ),

    #### asv_seq_index ####
    # `character` filename
    # sequence for each ASV
    asv_seq_index = tar_file(
      asv_seq_index,
      optimotu.pipeline::fastx_gz_index(asv_seq),
      deployment = "main"
    ),

    #### asv_best_hit_taxon ####
    # `tibble`:
    #  `seq_id` character: unique ASV identifier
    #  `ref_id` character: reference sequence id of best hit
    #  `sh_id` character: species hypothesis of best hit
    #  `dist` numeric: distance to best hit
    #  {{INGROUP_RANK}} character: taxon at rank INGROUP_RANK (e.g. kingdom) of
    #    best hit
    asv_best_hit_taxon = tar_fst_tbl(
      asv_best_hit_taxon,
      dplyr::left_join(
        asv_names,
        best_hit_taxon,
        by = "seq_idx"
      ) |>
        dplyr::select(-seq_idx),
      deployment = "main"
    ),

    #### asv_taxsort ####
    # `tibble`:
    #  `seq_idx` integer: index of sequence in asv_taxsort_seq
    #  `seq_idx_in` integer: index of sequence in asv_seq
    asv_taxsort = tar_fst_tbl(
      asv_taxsort,
      tibble::rowid_to_column(asv_tax, "seq_idx_in") |>
        dplyr::arrange(
          dplyr::across(all_of(!!c(optimotu.pipeline::tax_ranks(), "seq_id")))
        ) |>
        tibble::rowid_to_column("seq_idx") |>
        dplyr::select(seq_idx, seq_idx_in),
      deployment = "main"
    ),

    #### asv_taxsort_seq ####
    # `character` filename
    # sequence for each ASV
    asv_taxsort_seq = tar_file(
      asv_taxsort_seq,
      optimotu.pipeline::write_sequence(
        Biostrings::readDNAStringSet(asv_seq)[asv_taxsort$seq_idx_in],
        file.path(
          !!optimotu.pipeline::asv_path(),
          !!(if (optimotu.pipeline::do_rarefy()) {
            quote(sprintf("asv_taxsort_%s.fasta.gz", .rarefy_text))
          } else {
            "asv_taxsort.fasta.gz"
          })
        ),
        compress = TRUE
      ),
      deployment = "main"
    ),

    #### asv_taxsort_seq_index ####
    # `character` filename
    # sequence for each ASV
    asv_taxsort_seq_index = tar_file(
      asv_taxsort_seq_index,
      optimotu.pipeline::fastx_gz_index(asv_taxsort_seq),
      deployment = "main"
    )
  ),

  if (optimotu.pipeline::do_model_align()) {
    list(
      #### aligned_taxsort_seq ####
      # `character` filename
      # aligned sequence for eah ASV, sorted by protax taxonomy
      aligned_taxsort_seq = tar_file(
        aligned_taxsort_seq,
        optimotu.pipeline::write_sequence(
          # Use BString instead of DNAString because it will preserve case
          Biostrings::readBStringSet(asv_model_align)[asv_taxsort$seq_idx_in],
          file.path(
            !!optimotu.pipeline::aligned_path(),
            !!(if (optimotu.pipeline::do_rarefy()) {
              quote(sprintf("aligned_taxsort_%s.fasta.gz", .rarefy_text))
            } else {
              "aligned_taxsort.fasta.gz"
            })
          ),
          compress = TRUE
        )
      ),

      #### aligned_taxsort_seq_index ####
      # `character` filename
      # sequence for each ASV
      aligned_taxsort_seq_index = tar_file(
        aligned_taxsort_seq_index,
        optimotu.pipeline::fastx_gz_index(aligned_taxsort_seq),
        deployment = "main"
      )
    )
  },

  list(
    #### seqbatch_result_map ####
    seqbatch_result_map = tar_fst_tbl(
      seqbatch_result_map,
      seqbatch |>
        dplyr::transmute(
          seq_idx,
          result = as.raw(
            0x10 * (!seq_idx %in% denovo_chimeras) +
              0x20 * (!seq_idx %in% ref_chimeras) +
              !!( if (optimotu.pipeline::do_spike())
                quote(0x40 * (!seq_idx %in% spikes$seq_idx))
                else 0
              ) +
              !!(
                if (optimotu.pipeline::do_model_filter())
                  quote(0x80 * (seq_idx %in% asv_full_length))
                else
                  0
              )
          )
        ),
      pattern = !!{
        p <- quote(map(seqbatch, ref_chimeras))
        if (optimotu.pipeline::do_spike())
          p[[length(p) + 1]] <- quote(spikes)
        if (optimotu.pipeline::do_pos_control())
          p[[length(p) + 1]] <- quote(pos_controls)
        if (optimotu.pipeline::do_model_filter())
          p[[length(p) + 1]] <- quote(asv_full_length)
        p
      },
      deployment = "main"
    ),

    #### asv_map ####
    asv_map = tar_fst_tbl(
      asv_map,
      seqbatch_result_map |>
        dplyr::left_join(asv_names, by = "seq_idx") |>
        dplyr::select(seq_idx, result, seq_id),
      deployment = "main"
    )
  )
)

optimotu_plan <- c(optimotu_plan, asv_plan)
