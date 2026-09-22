# If not using LULU, then run uncross on the raw seqtable.
# otherwise, run LULU then uncross the LULU output
seqtable_pre_uncross <- quote(seqtable_raw)
if (optimotu.pipeline::do_lulu()) {
  seqtable_pre_uncross <- quote(seqtable_lulu)
}

# the "final" seqtable depends on whether we are doing tag-jump filtering,
# LULU, both, or neither
seqtable_final <- quote(seqtable_raw)
seqtable_final_name <- "seqtable_raw"
if (isTRUE(optimotu.pipeline::do_tag_jump())) {
  seqtable_final <- quote(seqtable_uncross)
  seqtable_final_name <- "seqtable_uncross"
} else if (optimotu.pipeline::do_lulu()) {
  seqtable_final <- quote(seqtable_lulu)
  seqtable_final_name <- "seqtable_lulu"
}

#### orientation plan ####
# the orientation plan consists of targets which need to be run within each
# sequencing run separately for different read orientations.
# There are minor variants for sequencing runs which are entirely either
# forward or reverse oriented, vs. those which contain both orientations.
#
# Library class is `orient_group` (fwd / rev / mixfwd / mixrev). File-path
# `orient` stays fwd/rev. Disjoint groups annotate read_map at batch scope;
# mixfwd/mixrev merge per aligned batch, then annotate.

# Single-orientation slice: used for half-runs and for disjoint libraries
# inside both-runs. Annotate read_map here when LULU and/or UNCROSS run.
orientation_plan_single <- c(
  readwise_plan,
  samplewise_plan
)
if (
  optimotu.pipeline::do_lulu() ||
    isTRUE(optimotu.pipeline::do_tag_jump())
) {
  if (optimotu.pipeline::do_unoise()) {
    orientation_plan_single[["read_map"]] <-
      tar_target(
        read_map,
        !!optimotu.pipeline::with_read_map_annotate(quote(
          optimotu.pipeline::unoise_read_map(
            sample = samplewise_meta$sample_key,
            fq_raw = samplewise_meta$fastq_R1,
            fq_trim = samplewise_meta$trim_R1,
            fq_merged = predenoise_merged,
            uc = unoise,
            denoise_map = denoise_map,
            vsearch = !!optimotu.pipeline::find_vsearch()
          )
        )),
        pattern = map(samplewise_meta, predenoise_merged, unoise, denoise_map),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "wide")
        )
      )
  } else {
    orientation_plan_single[["read_map"]] <-
      tar_target(
        read_map,
        !!optimotu.pipeline::with_read_map_annotate(quote(
          optimotu.pipeline::dada2_read_map(
            sample = samplewise_meta$sample_key,
            fq_raw = samplewise_meta$fastq_R1,
            fq_trim = samplewise_meta$trim_R1,
            fq_filt = samplewise_meta$filt_R1,
            dadaF = denoise_R1,
            derepF = derep_R1,
            dadaR = denoise_R2,
            derepR = derep_R2,
            merged = merged,
            denoise_map = denoise_map
          )
        )),
        pattern = map(
          samplewise_meta,
          denoise_R1,
          derep_R1,
          denoise_R2,
          derep_R2,
          merged,
          denoise_map
        ),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "wide")
        )
      )
  }
}

# Both-run disjoint fwd/rev libraries (orient_group == file orient).
orientation_plan_both_disjoint <- tar_map(
  values = tibble::tibble(
    .orient_group = c("fwd", "rev"),
    .orient = c("fwd", "rev"),
    .keep_empty = c(FALSE, FALSE)
  ),
  names = .orient_group,
  orientation_plan_single
)

# Mixed libraries: unannotated batch read_maps, merged below per batch.
# .orient is the file/primer orientation; names use orient_group.
readwise_plan_both_mixed <- readwise_plan
readwise_plan_both_mixed[["trim"]] <- tar_file(
      trim,
      optimotu.pipeline::trim_raw_pairs(
        pairs_meta = readwise_meta,
        seqrun = .seqrun,
        orient = .orient,
        trim_options = !!optimotu.pipeline::trim_options(),
        primer_R1 = !!optimotu.pipeline::trim_primer_R1(),
        primer_R2 = !!optimotu.pipeline::trim_primer_R2(),
        cutadapt = !!optimotu.pipeline::find_cutadapt(),
        log_tag = .orient_group,
        raw_R1 = raw_R1,
        raw_R2 = raw_R2
      ),
      pattern = map(readwise_meta, raw_R1, raw_R2),
      resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
    )
  
  samplewise_plan_both_mixed <- samplewise_plan
  samplewise_plan_both_mixed[["samplewise_meta"]] <- tar_target(
    samplewise_meta,
    dplyr::select(readwise_meta, orient, readwise_key) |>
      dplyr::left_join(sample_table, by = c("orient", "readwise_key")) |>
      dplyr::filter(
        !!(if (optimotu.pipeline::do_rarefy()) {
          quote(rarefy_text == .rarefy_text)
        } else {
          TRUE
        })
      ) |>
      dplyr::select(
        seqrun,
        sample,
        readwise_key,
        sample_key,
        fastq_R1,
        trim_R1,
        filt_R1,
        filt_R2,
        to_denoise_R1,
        to_denoise_R2,
        any_of(c(
          "rarefy_text",
          "numerator",
          "denominator",
          "number",
          "tar_seed"
        ))
      ),
    pattern = map(readwise_meta)
  )

orientation_plan_both_mixed <- tar_map(
  values = tibble::tibble(
    .orient_group = c("mixfwd", "mixrev"),
    .orient = c("fwd", "rev"),
    .keep_empty = c(TRUE, TRUE)
  ),
  names = .orient_group,
  readwise_plan_both_mixed,
  samplewise_plan_both_mixed
)

has_mixed_orient_keys <- any(
  optimotu.pipeline::sample_table()$orient_group %in% c("mixfwd", "mixrev")
)

#### seq_all_trim ####
# This is the file which is used for all subsequent processing.
# It is `seq_all` if primers were trimmed prior to denoising, otherwise it is
# `seq_trim`.
# `seq_all_trim` is the quoted name of the target which generates the file;
# `seq_all_trim_file` is the actual file name.
# This is *not* itself a target, it is just used for convenience.
# It should always be pre-evaluated with !!

seq_all_file <- file.path(
  optimotu.pipeline::asv_path(),
  "seq_all.fasta.gz"
)

seq_trim_file <- file.path(
  optimotu.pipeline::asv_path(),
  "seq_all_trim.fasta.gz"
)

if (optimotu.pipeline::trim_options()$action == "trim") {
  seq_all_trim <- quote(seq_all)
  seq_all_trim_file <- seq_all_file
} else {
  seq_all_trim <- quote(seq_trim)
  seq_all_trim_file <- seq_trim_file
}
seq_index_file <- paste0(seq_all_trim_file, ".index.qs2")

#### seqrun_plan ####

# the seqrun plan consists of steps that are run once per sequencing run.

seqrun_targets <- list(
  ##### denoise_read_counts_{.seqrun}_{.rarefaction?}_{.replicate?} #####
  # tibble:
  #  `sample_key` character: as `sample_table$sample_key`
  #  `denoise_nread` integer: number of sequences in the sample after denoising
  denoise_read_counts = tar_fst_tbl(
    denoise_read_counts,
    dplyr::summarize(
      seqtable_raw,
      denoise_nread = sum(nread),
      .by = sample
    ) |>
      dplyr::rename(sample_key = sample),
    resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
  ),

  ##### bimera_table_{.seqrun}_{.rarefaction?}_{.replicate?} #####
  # tibble:
  #  `nflag` integer: number of samples in which the sequence was considered
  #    chimeric
  #  `nsam` integer: number of samples in which the sequence occurred
  #  `seq` character: sequence
  #
  # find denovo chimeric sequences in each sample independently
  tar_fst_tbl(
    bimera_table,
    optimotu.pipeline::bimera_denovo_table(
      !!(seqtable_final),
      seq_all,
      allowOneOff = TRUE,
      multithread = optimotu.pipeline::local_cpus()
    ),
    resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
  )
)

if (optimotu.pipeline::do_dada2()) {
  seqrun_targets <- c(
    list(
      ##### errfun_{.seqrun} #####
      errfun = tar_target(
        errfun,
        optimotu.pipeline::choose_dada_error_function(raw_R2),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      )
    ),
    seqrun_targets
  )
}

if (isTRUE(optimotu.pipeline::do_lulu())) {
  seqrun_targets <- c(
    seqrun_targets,
    list(
      ##### seqrun_sentinel_{.seqrun}_{.rarefaction?}_{.replicate?} #####
      # character: a hash value
      #
      # This sentinel exists to ensure that lulu_table is calculated with an
      # updated seq_all_trim_file and seq_index_file, without introducing those
      # files as dependencies for lulu_table, because by design changes to
      # those files should not break targets calculated on earlier sequencing
      # runs.
      # Mentioning `seq_index` orders this target after the index (and thus
      # after seq_all_trim), but the value depends only on `seqtable_raw`, so
      # a rebuilt index does not invalidate LULU for unchanged runs.
      seqrun_sentinel = tar_target(
        seqrun_sentinel,
        {
          seq_index
          targets:::hash_object(seqtable_raw)
        },
        pattern = map(seqtable_raw),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      ),

      ##### lulu_match_table{.seqrun}_{.rarefaction?}_{.replicate?} #####
      # tibble:
      #  `seq_idx1` integer: index of first sequence in seq_all
      #  `seq_idx2` integer: index of second sequence in seq_all
      #  `nread1` integer: number of reads for the first sequence
      #  `nread2` integer: number of reads for the second sequence
      #  `dist` numeric: pairwise distance between the two sequences in [0,1]
      #
      # pairwise distances between ASVs in each sample (gap-filtered)
      lulu_match_table = if (
        optimotu.pipeline::lulu_dist_config()$method == "hamming"
      ) {
        if (!optimotu.pipeline::do_model_align()) {
          stop(
            "lulu_dist_config method is 'hamming' but do_model_align is ",
            "FALSE. This is not a valid configuration because the Hamming ",
            "distance requires sequences to be aligned."
          )
        }
        tar_fst_tbl(
          lulu_match_table,
          dplyr::reframe(
            seqtable_raw,
            optimotu.pipeline::lulu_distmx(
              seqall_file = seq_model_align, # this one does trigger dependency
              seqtable = dplyr::pick(seq_idx, nread),
              threshold = !!optimotu.pipeline::lulu_max_dist(),
              dist_config = !!(optimotu.pipeline::lulu_dist_config()$call),
              sentinel = seqrun_sentinel
            ),
            .by = sample
          ) |>
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
            ) |>
            dplyr::select(seq_idx1, seq_idx2, nread1, nread2, dist),
          pattern = map(seqtable_raw, seqrun_sentinel),
          resources = tar_resources(
            crew = tar_resources_crew(controller = "wide")
          )
        )
      } else {
        tar_fst_tbl(
          lulu_match_table,
          dplyr::reframe(
            seqtable_raw,
            optimotu.pipeline::lulu_distmx(
              seqall_file = seq_all_trim_file, # does not trigger dependency
              seqall_index = !!seq_index_file, # does not trigger dependency
              seqtable = dplyr::pick(seq_idx, nread),
              threshold = !!optimotu.pipeline::lulu_max_dist(),
              dist_config = !!(optimotu.pipeline::lulu_dist_config()$call),
              sentinel = seqrun_sentinel
            ),
            .by = sample
          ) |>
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
            ) |>
            dplyr::select(seq_idx1, seq_idx2, nread1, nread2, dist),
          pattern = map(seqtable_raw, seqrun_sentinel),
          resources = tar_resources(
            crew = tar_resources_crew(controller = "wide")
          )
        )
      },

      ##### lulu_seq_map_batch_{.seqrun} #####
      # Sparse non-identity LULU parents for seqs that occur in only this
      # seqtable_raw / lulu_match_table group (one batch).
      lulu_seq_map_batch = tar_fst_tbl(
        lulu_seq_map_batch,
        optimotu.pipeline::lulu_map_scoped(
          stats = lulu_seq_stats,
          match_table = lulu_match_table,
          scope = "batch",
          max_dist = !!optimotu.pipeline::lulu_max_dist(),
          min_abundance_ratio = !!optimotu.pipeline::lulu_min_abundance_ratio(),
          min_cooccurrence_ratio = !!optimotu.pipeline::lulu_min_cooccurrence_ratio(),
          use_mean_abundance_ratio = !!optimotu.pipeline::lulu_use_mean_abundance_ratio()
        ),
        pattern = map(lulu_match_table),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "wide")
        ),
        retrieval = "none"
      ),

      ##### lulu_seq_map_seqrun_{.seqrun} #####
      # Sparse non-identity LULU parents for seqs that occur in multiple
      # batches of this sequencing run but not in other seqruns.
      lulu_seq_map_seqrun = tar_fst_tbl(
        lulu_seq_map_seqrun,
        optimotu.pipeline::lulu_map_scoped(
          stats = lulu_seq_stats,
          match_table = lulu_match_table,
          scope = "seqrun",
          max_dist = !!optimotu.pipeline::lulu_max_dist(),
          min_abundance_ratio = !!optimotu.pipeline::lulu_min_abundance_ratio(),
          min_cooccurrence_ratio = !!optimotu.pipeline::lulu_min_cooccurrence_ratio(),
          use_mean_abundance_ratio = !!optimotu.pipeline::lulu_use_mean_abundance_ratio()
        ),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "wide")
        ),
        retrieval = "none"
      ),

      ##### seqtable_lulu_{.seqrun}_{.rarefaction?}_{.replicate?} #####
      # `tibble` with columns:
      #   `sample` (character) sample name as given in sample_table$sample_key
      #   `seq_idx` (integer) index of a sequence in seq_all
      #   `nread` (integer) number of reads
      seqtable_lulu = tar_fst_tbl(
        seqtable_lulu,
        optimotu.pipeline::lulu_table(lulu_seq_map, seqtable_raw),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      )
    )
  )
}

if (isTRUE(optimotu.pipeline::do_tag_jump())) {
  seqrun_targets <- c(
    seqrun_targets,
    list(
      ##### uncross_{.seqrun}_{.rarefaction?}_{.replicate?} #####
      # `tibble`:
      #   `sample (character) - sample name as given in sample_table$filt_key
      #   `nread` (integer) - number of reads in that sample (for this ASV)
      #   `total` (integer) - number of reads of that ASV across all samples (for
      #     this ASV)
      #   `uncross` (numeric) - UNCROSS score
      #   `is_tag_jump` (logical) - whether this occurrence is considered a likely
      #     tag jump
      #   `seq_idx` (integer) - index of a sequence in seq_all; same id as the
      #     input table (`seqtable_lulu` when LULU is on, else `seqtable_raw`).
      #     Row order also matches that input.
      #
      # remove tag-jumps (UNCROSS2)
      uncross = tar_fst_tbl(
        uncross,
        optimotu.pipeline::remove_tag_jumps(
          !!(seqtable_pre_uncross), # raw_ASV_table
          !!optimotu.pipeline::tag_jump_f(), # f-value (expected cross-talk rate)
          !!optimotu.pipeline::tag_jump_p(), # p-value (power to rise the exponent)
          "seq_idx" # name of column which uniquely identifies the sequence
        ),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        ),
        cue = tar_cue()
      ),

      ##### uncross_summary_{.seqrun}_{.rarefaction?}_{.replicate?} #####
      # `tibble`:
      #   `sample (character) - sample name as given in sample_table$filt_key
      #   `Total_reads` (integer) - total reads in the sample (all ASVs)
      #   `Number_of_TagJump_Events` (integer) - number of ASVs in the sample which
      #     are considered tag jumps.
      #   `TagJump_reads` (integer) - number of reads in the sample which belong to
      #     tag-jump ASVs.
      uncross_summary = tar_fst_tbl(
        uncross_summary,
        optimotu.pipeline::summarize_uncross(uncross),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      ),

      ##### uncross_read_counts_{.seqrun}_{.rarefaction?}_{.replicate?} #####
      # tibble:
      #  `sample_key` character: as `sample_table$sample_key`
      #  `uncross_nread` integer: number of sequences in the sample after denoising
      uncross_read_counts = tar_fst_tbl(
        uncross_read_counts,
        uncross_summary |>
          dplyr::transmute(
            sample_key = sample,
            uncross_nread = Total_reads - TagJump_reads
          ),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      ),
      ##### seqtable_uncross_{.seqrun}_{.rarefaction?}_{.replicate?} #####
      # `tibble`:
      #   `sample (character) - sample name as given in sample_table$sample_key
      #   `seq_idx` (integer) - index of a sequence in seq_all
      #   `nread` (integer) number of reads
      seqtable_uncross = tar_fst_tbl(
        seqtable_uncross,
        (!!seqtable_pre_uncross)[!uncross$is_tag_jump, ],
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      )
    )
  )
}

# There are a few seqrun targets which differ between forward, reverse, and "both" runs
seqrun_forward_targets <- c(
  seqrun_targets,
  list(
    ##### seq_merged_{.seqrun}_{.rarefaction?}_{.replicate?} #####
    seq_merged = if (optimotu.pipeline::do_unoise()) {
      tar_target(
        seq_merged,
        unique(as.character(unlist(lapply(unoise, \(x) x$clusters$seq)))),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      )
    } else {
      tar_target(
        seq_merged,
        unique(as.character(unlist(lapply(merged, \(x) x$sequence)))),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      )
    }
  )
)

seqrun_reverse_targets <- c(
  seqrun_targets,
  list(
    ##### seq_merged_{.seqrun}_{.rarefaction?}_{.replicate?} #####
    seq_merged = if (optimotu.pipeline::do_unoise()) {
      tar_target(
        seq_merged,
        unique(as.character(unlist(lapply(unoise, \(x) {
          dada2::rc(x$clusters$seq)
        })))),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      )
    } else {
      tar_target(
        seq_merged,
        unique(as.character(unlist(lapply(merged, \(x) {
          dada2::rc(x$sequence)
        })))),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      )
    }
  )
)

seqrun_both_seq_merged <- if (optimotu.pipeline::do_unoise()) {
  if (has_mixed_orient_keys) {
    quote(unique(as.character(unlist(c(
      lapply(unoise_fwd, \(x) x$clusters$seq),
      lapply(unoise_rev, \(x) dada2::rc(x$clusters$seq)),
      lapply(unoise_mixfwd, \(x) x$clusters$seq),
      lapply(unoise_mixrev, \(x) dada2::rc(x$clusters$seq))
    )))))
  } else {
    quote(unique(as.character(unlist(c(
      lapply(unoise_fwd, \(x) x$clusters$seq),
      lapply(unoise_rev, \(x) dada2::rc(x$clusters$seq))
    )))))
  }
} else if (has_mixed_orient_keys) {
  quote(unique(as.character(unlist(c(
    lapply(merged_fwd, \(x) x$sequence),
    lapply(merged_rev, \(x) dada2::rc(x$sequence)),
    lapply(merged_mixfwd, \(x) x$sequence),
    lapply(merged_mixrev, \(x) dada2::rc(x$sequence))
  )))))
} else {
  quote(unique(as.character(unlist(c(
    lapply(merged_fwd, \(x) x$sequence),
    lapply(merged_rev, \(x) dada2::rc(x$sequence))
  )))))
}

seqrun_both_seqtable_raw <- if (has_mixed_orient_keys) {
  quote(
    dplyr::bind_rows(
      seqtable_raw_fwd,
      seqtable_raw_rev,
      seqtable_raw_mixfwd,
      seqtable_raw_mixrev
    ) |>
      dplyr::summarize(nread = sum(nread), .by = c(sample, seq_idx)) |>
      dplyr::mutate(
        tar_group = (as.integer(factor(sample)) - 1L) %/% 96L + 1L
      )
  )
} else {
  quote(
    dplyr::bind_rows(seqtable_raw_fwd, seqtable_raw_rev) |>
      dplyr::summarize(nread = sum(nread), .by = c(sample, seq_idx)) |>
      dplyr::mutate(
        tar_group = (as.integer(factor(sample)) - 1L) %/% 96L + 1L
      )
  )
}

seqrun_both_targets <- c(
  seqrun_targets,
  list(
    ##### seq_merged_{.seqrun}_{.rarefaction?}_{.replicate?} #####
    seq_merged = tar_target(
      seq_merged,
      !!seqrun_both_seq_merged,
      resources = tar_resources(
        crew = tar_resources_crew(controller = "thin")
      )
    ),

    ##### seqtable_raw_{.seqrun}_{.rarefaction?}_{.replicate?} #####
    # `tibble`:
    #   `sample (character) - sample name as given in sample_table$sample_key
    #   `seq_idx` (integer) - index of a sequence in seq_all
    #   `nread` (integer) number of reads
    #
    # Combines disjoint and (when present) mixed orientation seqtables.
    seqtable_raw = tar_fst_tbl(
      seqtable_raw,
      !!seqrun_both_seqtable_raw,
      iteration = "group",
      resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
    )
  ),
  if (has_mixed_orient_keys) {
    list(
      ##### read_map_mixed_{.seqrun}_{.rarefaction?}_{.replicate?} #####
      # Batch merge of mixfwd/mixrev fate maps (aligned by sample_key groups).
      # No seqrun-level combined read_map is stored.
      read_map_mixed = tar_target(
        read_map_mixed,
        !!optimotu.pipeline::with_read_map_annotate(quote(
          optimotu.pipeline::merge_read_maps(read_map_mixfwd, read_map_mixrev)
        )),
        pattern = map(read_map_mixfwd, read_map_mixrev),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "wide")
        )
      )
    )
  } else {
    list()
  }
)

# with "both" orientation we also need to consider all orientation slices of raw_R2
if (optimotu.pipeline::do_dada2()) {
  seqrun_both_errfun <- if (has_mixed_orient_keys) {
    quote(
      optimotu.pipeline::choose_dada_error_function(
        unique(c(raw_R2_fwd, raw_R2_rev, raw_R2_mixfwd, raw_R2_mixrev))
      )
    )
  } else {
    quote(
      optimotu.pipeline::choose_dada_error_function(
        unique(c(raw_R2_fwd, raw_R2_rev))
      )
    )
  }
  seqrun_both_targets$errfun <- tar_target(
    errfun,
    !!seqrun_both_errfun,
    resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
  )
}


# one row for each sequencing run
seqrun_orient_meta <- dplyr::summarize(
  optimotu.pipeline::sample_table(),
  .orient = dplyr::case_when(
    all(orient == "fwd") ~ "fwd",
    all(orient == "rev") ~ "rev",
    TRUE ~ "both"
  ),
  .by = seqrun
)

# one row for each sequencing run which has _only_ forward orientation
seqrun_forward_meta <- dplyr::filter(seqrun_orient_meta, .orient == "fwd") |>
  dplyr::select(.seqrun = seqrun, .orient) |>
  dplyr::mutate(.orient_group = .orient, .keep_empty = FALSE)
seqrun_forward_plan <- tar_map(
  values = seqrun_forward_meta,
  names = .seqrun,
  orientation_plan_single,
  seqrun_forward_targets
)


#  `seqrun` character; name of sequencing run (directory in sequences/01_raw)
#  `sample` character; name of sample, based on parsing file name
#  `readwise_key`character; common prefix of trim_R1, trim_R2, filt_R1 and
#      filt_R2
#  `sample_key` character; as `sample_key` but also with rarefy_text
#  `fastq_R1` character; file name with path for raw R1 file
#  `trim_R1` character; file name with path for trimmed R1 file
#  `filt_R1` character; file name with path for filtered R1 file
#  `filt_R2` character; file name with path for filtered R2 file
#  `to_denoise_R1` character; file name with path for R1 to denoise
#  `to_denoise_R2` character; file name with path for R2 to denoise
#  `rarefy_text` (optional) character; specification of the rarefaction
#  `numerator` (optional) integer; numerator for fractional rarefacation
#  `denominator` (optional) integer; denominator for fractional rarefacation
#  `number` (optional) integer; count for count.based rarefaction
#  `tar_seed` integer; random seed to use for dereplication
samplewise_dummy <- tibble::tibble(
  seqrun = character(),
  sample = character(),
  readwise_key = character(),
  sample_key = character(),
  fastq_R1 = character(),
  trim_R1 = character(),
  filt_R1 = character(),
  filt_R2 = character(),
  to_denoise_R1 = character(),
  to_denoise_R2 = character(),
  merged = character(),
  to_denoise_merged = character()
)

merge_read_counts_dummy <- tibble::tibble(
  merged = character(),
  merge_nread = integer()
)
filt_read_counts_dummy <- tibble::tibble(
  filt_R1 = character(),
  filt_nread = integer()
)

# add "dummy" values if there were no forward-only seqruns
if (nrow(seqrun_forward_meta) == 0) {
  seqrun_forward_plan$samplewise_meta <- list(tar_fst_tbl(
    samplewise_meta_dummy_fwd,
    samplewise_dummy,
    deployment = "main"
  ))
  seqrun_forward_plan$raw_read_counts <- list(tar_fst_tbl(
    raw_read_counts_dummy_fwd,
    tibble::tibble(
      fastq_file = character(),
      raw_nread = integer()
    ),
    deployment = "main"
  ))
  seqrun_forward_plan$trim_read_counts <- list(tar_fst_tbl(
    trim_read_counts_dummy_fwd,
    tibble::tibble(
      trim_R1 = character(),
      trim_nread = integer()
    ),
    deployment = "main"
  ))
  if (optimotu.pipeline::do_dada2()) {
    seqrun_forward_plan$filt_read_counts <- list(tar_fst_tbl(
      filt_read_counts_dummy_fwd,
      filt_read_counts_dummy,
      deployment = "main"
    ))
  }
  if (optimotu.pipeline::do_unoise()) {
    seqrun_forward_plan$merge_read_counts <- list(tar_fst_tbl(
      merge_read_counts_dummy_fwd,
      merge_read_counts_dummy,
      deployment = "main"
    ))
  }
}

# one row for each sequencing run which has _only_ reverse orientation
seqrun_reverse_meta <- dplyr::filter(seqrun_orient_meta, .orient == "rev") |>
  dplyr::select(.seqrun = seqrun, .orient) |>
  dplyr::mutate(.orient_group = .orient, .keep_empty = FALSE)
seqrun_reverse_plan <- tar_map(
  values = seqrun_reverse_meta,
  names = .seqrun,
  orientation_plan_single,
  seqrun_reverse_targets
)

# add "dummy" values if there were no reverse-only seqruns
if (nrow(seqrun_reverse_meta) == 0) {
  seqrun_reverse_plan$samplewise_meta <- list(tar_fst_tbl(
    samplewise_meta_dummy_rev,
    samplewise_dummy,
    deployment = "main"
  ))
  seqrun_reverse_plan$raw_read_counts <- list(tar_fst_tbl(
    raw_read_counts_dummy_rev,
    tibble::tibble(
      fastq_file = character(),
      raw_nread = integer()
    ),
    deployment = "main"
  ))
  seqrun_reverse_plan$trim_read_counts <- list(tar_fst_tbl(
    trim_read_counts_dummy_rev,
    tibble::tibble(
      trim_R1 = character(),
      trim_nread = integer()
    ),
    deployment = "main"
  ))
  if (optimotu.pipeline::do_dada2()) {
    seqrun_reverse_plan$filt_read_counts <- list(tar_fst_tbl(
      filt_read_counts_dummy_rev,
      filt_read_counts_dummy,
      deployment = "main"
    ))
  }
  if (optimotu.pipeline::do_unoise()) {
    seqrun_reverse_plan$merge_read_counts <- list(tar_fst_tbl(
      merge_read_counts_dummy_rev,
      merge_read_counts_dummy,
      deployment = "main"
    ))
  }
}

# one row for each sequencing run which has both orientations
seqrun_both_meta <- dplyr::filter(seqrun_orient_meta, .orient == "both") |>
  dplyr::select(.seqrun = seqrun)
seqrun_both_plan <- if (has_mixed_orient_keys) {
  tar_map(
    values = seqrun_both_meta,
    names = .seqrun,
    orientation_plan_both_disjoint,
    orientation_plan_both_mixed,
    seqrun_both_targets
  )
} else {
  tar_map(
    values = seqrun_both_meta,
    names = .seqrun,
    orientation_plan_both_disjoint,
    seqrun_both_targets
  )
}

# add "dummy" values if there were no both-orientation seqruns
if (nrow(seqrun_both_meta) == 0) {
  seqrun_both_plan$samplewise_meta_fwd <- list(tar_fst_tbl(
    samplewise_meta_fwd_dummy_both,
    samplewise_dummy,
    deployment = "main"
  ))
  seqrun_both_plan$samplewise_meta_rev <- list(tar_fst_tbl(
    samplewise_meta_rev_dummy_both,
    samplewise_dummy,
    deployment = "main"
  ))
  seqrun_both_plan$raw_read_counts_fwd <- list(tar_fst_tbl(
    raw_read_counts_fwd_dummy_both,
    tibble::tibble(
      fastq_file = character(),
      raw_nread = integer()
    ),
    deployment = "main"
  ))
  seqrun_both_plan$raw_read_counts_rev <- list(tar_fst_tbl(
    raw_read_counts_rev_dummy_both,
    tibble::tibble(
      fastq_file = character(),
      raw_nread = integer()
    ),
    deployment = "main"
  ))
  seqrun_both_plan$trim_read_counts_fwd <- list(tar_fst_tbl(
    trim_read_counts_fwd_dummy_both,
    tibble::tibble(
      trim_R1 = character(),
      trim_nread = integer()
    ),
    deployment = "main"
  ))
  seqrun_both_plan$trim_read_counts_rev <- list(tar_fst_tbl(
    trim_read_counts_rev_dummy_both,
    tibble::tibble(
      trim_R1 = character(),
      trim_nread = integer()
    ),
    deployment = "main"
  ))
  if (optimotu.pipeline::do_dada2()) {
    seqrun_both_plan$filt_read_counts_fwd <- list(tar_fst_tbl(
      filt_read_counts_fwd_dummy_both,
      filt_read_counts_dummy,
      deployment = "main"
    ))
    seqrun_both_plan$filt_read_counts_rev <- list(tar_fst_tbl(
      filt_read_counts_rev_dummy_both,
      filt_read_counts_dummy,
      deployment = "main"
    ))
  }
  if (optimotu.pipeline::do_unoise()) {
    seqrun_both_plan$merge_read_counts_fwd <- list(tar_fst_tbl(
      merge_read_counts_fwd_dummy_both,
      merge_read_counts_dummy,
      deployment = "main"
    ))
    seqrun_both_plan$merge_read_counts_rev <- list(tar_fst_tbl(
      merge_read_counts_rev_dummy_both,
      merge_read_counts_dummy,
      deployment = "main"
    ))
  }
}

# Mix-orientation dummies when the study has no mixfwd/mixrev libraries (or no
# both-runs). Output read-count assembly always binds these names.
if (!has_mixed_orient_keys || nrow(seqrun_both_meta) == 0) {
  seqrun_both_plan$samplewise_meta_mixfwd <- list(tar_fst_tbl(
    samplewise_meta_mixfwd_dummy_both,
    samplewise_dummy,
    deployment = "main"
  ))
  seqrun_both_plan$samplewise_meta_mixrev <- list(tar_fst_tbl(
    samplewise_meta_mixrev_dummy_both,
    samplewise_dummy,
    deployment = "main"
  ))
  seqrun_both_plan$raw_read_counts_mixfwd <- list(tar_fst_tbl(
    raw_read_counts_mixfwd_dummy_both,
    tibble::tibble(
      fastq_file = character(),
      raw_nread = integer()
    ),
    deployment = "main"
  ))
  seqrun_both_plan$raw_read_counts_mixrev <- list(tar_fst_tbl(
    raw_read_counts_mixrev_dummy_both,
    tibble::tibble(
      fastq_file = character(),
      raw_nread = integer()
    ),
    deployment = "main"
  ))
  seqrun_both_plan$trim_read_counts_mixfwd <- list(tar_fst_tbl(
    trim_read_counts_mixfwd_dummy_both,
    tibble::tibble(
      trim_R1 = character(),
      trim_nread = integer()
    ),
    deployment = "main"
  ))
  seqrun_both_plan$trim_read_counts_mixrev <- list(tar_fst_tbl(
    trim_read_counts_mixrev_dummy_both,
    tibble::tibble(
      trim_R1 = character(),
      trim_nread = integer()
    ),
    deployment = "main"
  ))
  if (optimotu.pipeline::do_dada2()) {
    seqrun_both_plan$filt_read_counts_mixfwd <- list(tar_fst_tbl(
      filt_read_counts_mixfwd_dummy_both,
      filt_read_counts_dummy,
      deployment = "main"
    ))
    seqrun_both_plan$filt_read_counts_mixrev <- list(tar_fst_tbl(
      filt_read_counts_mixrev_dummy_both,
      filt_read_counts_dummy,
      deployment = "main"
    ))
  }
  if (optimotu.pipeline::do_unoise()) {
    seqrun_both_plan$merge_read_counts_mixfwd <- list(tar_fst_tbl(
      merge_read_counts_mixfwd_dummy_both,
      merge_read_counts_dummy,
      deployment = "main"
    ))
    seqrun_both_plan$merge_read_counts_mixrev <- list(tar_fst_tbl(
      merge_read_counts_mixrev_dummy_both,
      merge_read_counts_dummy,
      deployment = "main"
    ))
  }
}

seqrun_plan <- optimotu.pipeline::tar_merge(
  seqrun_forward_plan,
  seqrun_reverse_plan
) |>
  optimotu.pipeline::tar_merge(seqrun_both_plan)
