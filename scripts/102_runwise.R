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

# for single orientation (fwd or rev) we can add LULU/uncross when the
# read_map is created. For multi-orientation (both) we need to do it later,
# when read_map_fwd and read_map_rev are merged.
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

# for multiple orientations, we duplicate the readwise and samplewise plans
# for the two orientations
orientation_plan_multi <- tar_map(
  values = list(.orient = c("fwd", "rev")),
  names = .orient,
  readwise_plan,
  samplewise_plan
)

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
      # This sentinal exists to ensure that lulu_table is calculated with an
      # updated seq_all_trim_file and seq_index_file, without introducing those
      # files as dependencies for lulu_table, because by design changes to
      # those files should not break targets calculated on earlier sequencing
      # runs.
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
      #  `seq_id1` integer: index of first sequence in seq_all
      #  `seq_id2` integer: index of second sequence in seq_all
      #  `dist` numeric: pairwise distance between the two sequences in [0,1]
      #  `len` integer: length of the alignment between the two sequences
      #  `n_gap` integer: number of gaps in the alignment
      #  `max_gap` integer: length of the longest gap in the alignment
      #
      # pairwise distances between ASVs in each sample
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
          ),
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
              seqall_index = seq_index_file, # does not trigger dependency
              seqtable = dplyr::pick(seq_idx, nread),
              threshold = !!optimotu.pipeline::lulu_max_dist(),
              dist_config = !!(optimotu.pipeline::lulu_dist_config()$call),
              sentinel = seqrun_sentinel
            ),
            .by = sample
          ),
          pattern = map(seqtable_raw, seqrun_sentinel),
          resources = tar_resources(
            crew = tar_resources_crew(controller = "wide")
          )
        )
      },

      ##### seqtable_lulu_{.seqrun}_{.rarefaction?}_{.replicate?} #####
      # `tibble` with columns:
      #   `sample` (character) sample name as given in sample_table$sample_key
      #   `seq_idx` (integer) index of a sequence in seq_all
      #   `nread` (integer) number of reads
      seqtable_lulu = tar_fst_tbl(
        seqtable_lulu,
        optimotu.pipeline::lulu_table(lulu_asv_map, seqtable_raw),
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

seqrun_both_targets <- c(
  seqrun_targets,
  list(
    ##### seq_merged_{.seqrun}_{.rarefaction?}_{.replicate?} #####
    seq_merged = if (optimotu.pipeline::do_unoise()) {
      tar_target(
        seq_merged,
        unique(as.character(unlist(c(
          lapply(unoise_fwd, \(x) x$clusters$seq),
          lapply(unoise_rev, \(x) dada2::rc(x$clusters$seq))
        )))),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      )
    } else {
      tar_target(
        seq_merged,
        unique(as.character(unlist(c(
          lapply(merged_fwd, \(x) x$sequence),
          lapply(merged_rev, \(x) dada2::rc(x$sequence))
        )))),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      )
    },

    ##### seqtable_raw_{.seqrun}_{.rarefaction?}_{.replicate?} #####
    # `tibble`:
    #   `sample (character) - sample name as given in sample_table$sample_key
    #   `seq_idx` (integer) - index of a sequence in seq_all
    #   `nread` (integer) number of reads
    #
    # This combines seqtable_raw_fwd_{.seqrun} and seqtable_raw_rev_{.seqrun}
    seqtable_raw = tar_fst_tbl(
      seqtable_raw,
      dplyr::bind_rows(seqtable_raw_fwd, seqtable_raw_rev) |>
        dplyr::summarize(nread = sum(nread), .by = c(sample, seq_idx)) |>
        dplyr::mutate(
          tar_group = (as.integer(factor(sample)) - 1L) %/% 96L + 1L
        ),
      iteration = "group",
      resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
    ),

    ##### read_map_{.seqrun}_{.rarefaction?}_{.replicate?} #####
    # `tibble`:
    #   `sample (character) - sample name as given in sample_table$sample_key
    #   `raw_idx` (integer) - index of read in the un-rarified fastq file
    #   `seq_idx` (integer) - index of the current community-table ASV in
    #     seq_all (LULU parent when LULU ran)
    #   `prelulu_idx` (integer) - denoise-time ASV in seq_all; present only
    #     when LULU ran. A daughter is prelulu_idx != seq_idx.
    #   `flags` (raw) - bits for presence after each processing stage:
    #    0x01: trim
    #    0x02: filter
    #    0x04: denoise & merge
    #    0x08: survived tag-jump removal (if performed)
    #    0x10-0x80: reserved for asv_map$result (not set here)
    #
    # This combines read_map_fwd_{.seqrun} and read_map_rev_{.seqrun}
    #
    # If LULU and/or tag-jump removal is performed, it remaps seq_idx to the
    # LULU parent and/or adds the uncross information.
    read_map = tar_fst_tbl(
      read_map,
      !!optimotu.pipeline::with_read_map_annotate(quote(
        optimotu.pipeline::merge_read_maps(read_map_fwd, read_map_rev)
      )),
      resources = tar_resources(
        crew = tar_resources_crew(controller = "wide")
      )
    )
  )
)

# with "both" orientation we also need to consider both versions of raw_R2
if (optimotu.pipeline::do_dada2()) {
  seqrun_both_targets$errfun <- tar_target(
    errfun,
    optimotu.pipeline::choose_dada_error_function(
      unique(c(raw_R2_fwd, raw_R2_rev))
    ),
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
  dplyr::select(.seqrun = seqrun, .orient)
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
  dplyr::select(.seqrun = seqrun, .orient)
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
seqrun_both_plan <- tar_map(
  values = seqrun_both_meta,
  names = .seqrun,
  orientation_plan_multi,
  seqrun_both_targets
)

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

seqrun_plan <- optimotu.pipeline::tar_merge(
  seqrun_forward_plan,
  seqrun_reverse_plan
) |>
  optimotu.pipeline::tar_merge(seqrun_both_plan)
