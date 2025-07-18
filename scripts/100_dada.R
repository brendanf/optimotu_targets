# DADA2 quality filtering and denoising for Illumina paired-end metabarcoding data
# Brendan Furneaux
# Based on DADA2 analysis for GSSP from Jenni Hultman
# edits by Sten Anslan - account for reverse complementary oriented sequences and add UNCROSS2 tag-jumps filtering per run

library(targets)
library(tarchetypes)

#### inner_dada_plan ####
# the "inner" dada plan consists of targets which need to be done separately
# for each orientation, if not all reads are in the same orientation.
# These targets also need to be done separately for each sequencing run.

inner_dada_plan <- list(
  ##### dada2_meta_{.orient?}_{.seqrun} #####
  # grouped tibble:
  #  `seqrun` character; name of sequencing run (directory in sequences/01_raw)
  #  `sample` character; name of sample, based on parsing file name
  #  `fastq_R1` character; file name with path for raw R1 file
  #  `fastq_R2` character; file name with path for raw R2 file
  #  `trim_R1` character; file name with path for trimmed R1 file
  #  `trim_R2` character; file name with path for trimmed R2 file
  #  `filt_R1` character; file name with path for filtered R1 file
  #  `filt_R2` character; file name with path for filtered R2 file
  #  `sample_key`character; common prefix of trim_R1, trim_R2, filt_R1 and
  #      filt_R2; used as sample name by dada2 functions and read counts
  dada2_meta = tar_group_size(
    dada2_meta,
    sample_table |>
      dplyr::filter(orient == .orient, seqrun == .seqrun) |>
      dplyr::select(seqrun, sample, fastq_R1, fastq_R2, trim_R1, trim_R2,
                    filt_R1, filt_R2, sample_key,
                    any_of(optimotu.pipeline::cutadapt_paired_option_names),
                    any_of("maxEE")),
    size = 96
  ),

  ##### raw_R2_{.orient?}_{.seqrun} #####
  # character: path and file name
  # raw reads, for dependency tracking
  raw_R2 = tar_file(
    raw_R2,
    file.path(!!optimotu.pipeline::raw_path(), dada2_meta$fastq_R2),
    pattern = map(dada2_meta),
    resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
  ),

  ##### raw_read_counts_{.orient?}_{.seqrun} #####
  # tibble:
  #  `fastq_file` character: file name of raw R2 file
  #  `raw_nread` integer: number of sequences in the file
  raw_read_counts = tar_fst_tbl(
    raw_read_counts,
    tibble::tibble(
      fastq_file = raw_R2,
      raw_nread = optimotu.pipeline::sequence_size(fastq_file)
    ),
    pattern = map(raw_R2),
    resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
  ),


  ##### trim_{.orient?}_{.seqrun} #####
  # character: file names with path of trimmed read files (fastq.gz)
  #
  # remove adapters and barcodes
  # also do some preliminary quality filtering
  trim = tar_file(
    trim,
    optimotu.pipeline::trim_raw_pairs(
      pairs_meta = dada2_meta,
      seqrun = .seqrun,
      orient = .orient,
      trim_options = !!optimotu.pipeline::trim_options(),
      primer_R1 = !!optimotu.pipeline::trim_primer_R1(),
      primer_R2 = !!optimotu.pipeline::trim_primer_R2(),
      raw_path = !!optimotu.pipeline::raw_path()
    ),
    pattern = map(dada2_meta),
    resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
  ),

  ##### trim_read_counts_{.orient?}_{.seqrun} #####
  # tibble:
  #  `trim_R1` character: file name with path of trimmed R1 file
  #  `trim_nread` integer: number of sequences in the file
  #
  # count of reads per sample after adapter trimming
  trim_read_counts = tar_fst_tbl(
    trim_read_counts,
    tibble::tibble(
      trim_R1 = purrr::keep(trim, endsWith, "_R1_trim.fastq.gz"),
      trim_nread = optimotu.pipeline::sequence_size(trim_R1)
    ),
    pattern = map(trim),
    resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
  ),

  # DADA2 quality filtering on read-pairs

  ##### filter_pairs_{.orient?}_{.seqrun} #####
  # character: file names with path of filtered read files (fastq.gz)
  #
  # additional quality filtering on read-pairs
  filter_pairs = tar_file(
    filter_pairs,
    dada2_meta |>
      dplyr::mutate(
        trim_R1 = purrr::keep(trim, endsWith, "_R1_trim.fastq.gz"),
        trim_R2 = purrr::keep(trim, endsWith, "_R2_trim.fastq.gz")
      ) |>
      dplyr::group_by(
        dplyr::pick(any_of(c("seqrun", optimotu.pipeline::filter_option_names)))
      ) |>
      dplyr::group_map(
        ~ optimotu.pipeline::filterAndTrim(
          fwd = .x$trim_R1,
          filt = .x$filt_R1,
          rev = .x$trim_R2,
          filt.rev = .x$filt_R2,
          maxEE = update(!!optimotu.pipeline::dada2_maxEE(), .y), # max expected errors (fwd, rev)
          rm.phix = TRUE, #remove matches to phiX genome
          compress = TRUE, # write compressed files
          multithread = optimotu.pipeline::local_cpus(),
          verbose = TRUE
        )
      ) |>
      unlist(),
    pattern = map(dada2_meta, trim),
    resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
  ),

  ##### filt_read_counts_{.orient?}_{.seqrun} #####
  # tibble:
  #  `filt_R1` character: file name with path of filtered R1 file
  #  `filt_nread` integer: number of sequences in the file
  #
  # count of reads after filtering
  filt_read_counts = tar_fst_tbl(
    filt_read_counts,
    tibble::tibble(
      filt_R1 = purrr::keep(filter_pairs, endsWith, "_R1_filt.fastq.gz"),
      filt_nread = optimotu.pipeline::sequence_size(filt_R1)
    ),
    pattern = map(filter_pairs),
    resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
  ),

  ##### map over R1 and R2 #####
  # inside the tar_map, every occurrence of `read` is replaced by "R1" or "R2"
  # the read name is also appended to all target names
  # so all of this gets done separately for forward and reverse reads.
  # pattern=map() means we are also keeping the different sequencing runs
  # separate
  tar_map(
    values = list(read = c("R1", "R2")),

    ###### filtered_{read}_{.orient?}_{.seqrun} ######
    # character: path and file name of filtered reads; fastq.gz
    #
    # select only the files corresponding to the read we are working on
    filtered = tar_file(
      filtered,
      purrr::keep(filter_pairs, endsWith, paste0(read, "_filt.fastq.gz")),
      pattern = map(filter_pairs),
      deployment = "main"
    ),

    ###### derep_{read}_{.orient?}_{.seqrun} ######
    # list of dada2 `derep` objects
    #
    # dereplicate
    derep = tar_target(
      derep,
      optimotu.pipeline::derepFastq(
        filtered,
        verbose = TRUE,
        names = optimotu.pipeline::file_to_sample_key(filtered)
      ),
      pattern = map(filtered),
      resources = tar_resources(crew = tar_resources_crew(controller = "wide")) # memory
    ),

    ###### err_{read}_{.orient?}_{.seqrun} ######
    # list: see dada2::LearnErrors
    #
    # fit error profile
    err = tar_target(
      err,
      optimotu.pipeline::learnErrors(
        purrr::discard(filtered, grepl, pattern = "BLANK|NEG"),
        errorEstimationFunction = errfun,
        multithread = optimotu.pipeline::local_cpus(),
        verbose = TRUE
      ),
      resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
    ),

    ###### denoise_{read}_{.orient?}_{.seqrun} ######
    # list of dada2 `dada` objects
    denoise = tar_target(
      denoise,
      optimotu.pipeline::dada(
        derep,
        err = err,
        errorEstimationFunction = errfun,
        multithread = optimotu.pipeline::local_cpus(),
        verbose = TRUE
      ),
      pattern = map(derep),
      resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
    )
  ),

  ##### merged_{.orient?}_{.seqrun} #####
  # list of data.frame; see dada2::mergePairs
  #
  # Merge paired reads and make a sequence table for each sequencing run
  merged = tar_target(
    merged,
    optimotu.pipeline::mergePairs(
      denoise_R1,
      derep_R1,
      denoise_R2,
      derep_R2,
      minOverlap = 10,
      maxMismatch = 1,
      verbose = TRUE
    ),
    pattern = map(denoise_R1, derep_R1, denoise_R2, derep_R2),
    resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
  ),

  ##### seqtable_raw_{.orient?}_{.seqrun} #####
  # `tibble` with columns:
  #   `sample` (character) sample name as given in sample_table$sample_key
  #   `seq_idx` (integer) index of a sequence in seq_all
  #   `nread` (integer) number of reads
  seqtable_raw = tar_fst_tbl(
    seqtable_raw,
    optimotu.pipeline::make_mapped_sequence_table(
      merged,
      seq_all,
      rc = .orient == "rev"
    ),
    pattern = map(merged),
    resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
  ),

  ##### dada_map_{.orient?}_{.seqrun} #####
  # map the raw reads to the nochim ASVs
  # indexes are per-sample
  #
  # a tibble:
  #  sample: character identifying the sample, as in sample_table
  #  raw_idx: integer index of read in the un-rarified fastq file
  #  seq_idx: integer index of ASV in seq_all
  #  flags: raw, bits give presence/absence of the read after different stages:
  #    0x01: trim
  #    0x02: filter
  #    0x04: denoise & merge
  #    0x08: tag-jump removal (if performed)
  dada_map = tar_target(
    dada_map,
    mapply(
      FUN = optimotu.pipeline::seq_map,
      sample = dada2_meta$sample_key,
      fq_raw = file.path(!!optimotu.pipeline::raw_path(), dada2_meta$fastq_R1),
      fq_trim = dada2_meta$trim_R1,
      fq_filt = dada2_meta$filt_R1,
      dadaF = denoise_R1,
      derepF = derep_R1,
      dadaR = denoise_R2,
      derepR = derep_R2,
      merged = merged,
      MoreArgs = list(
        seq_all = seq_all,
        rc = .orient == "rev"
      ),
      SIMPLIFY = FALSE
    ) |>
      purrr::list_rbind(
        ptype = tibble::tibble(
          sample = character(),
          raw_idx = integer(),
          seq_idx = integer(),
          flags = raw()
        )
      ),
    pattern = map(dada2_meta, denoise_R1, derep_R1, denoise_R2, derep_R2, merged),
    resources = tar_resources(crew = tar_resources_crew(controller = "wide")) # for memory
  )
)

# for single orientation (fwd or rev) we can add the uncross information when the
# dada_map is created. For multi-orientation (both) we need to do it later, when
# dada_map_fwd and dada_map_rev are merged.
inner_dada_plan_single <- inner_dada_plan
if (optimotu.pipeline::do_tag_jump()) {
  inner_dada_plan_single[["dada_map"]] <-
    tar_target(
      dada_map,
      mapply(
        FUN = optimotu.pipeline::seq_map,
        sample = dada2_meta$sample_key,
        fq_raw = file.path(!!optimotu.pipeline::raw_path(), dada2_meta$fastq_R1),
        fq_trim = dada2_meta$trim_R1,
        fq_filt = dada2_meta$filt_R1,
        dadaF = denoise_R1,
        derepF = derep_R1,
        dadaR = denoise_R2,
        derepR = derep_R2,
        merged = merged,
        MoreArgs = list(
          seq_all = seq_all,
          rc = .orient == "rev"
        ),
        SIMPLIFY = FALSE
      ) |>
        purrr::list_rbind(
          ptype = tibble::tibble(
            sample = character(),
            raw_idx = integer(),
            seq_idx = integer(),
            flags = raw()
          )
        ) |>
        optimotu.pipeline::add_uncross_to_seq_map(seqtable_raw, uncross),
      pattern = map(dada2_meta, denoise_R1, derep_R1, denoise_R2, derep_R2, merged),
      resources = tar_resources(crew = tar_resources_crew(controller = "wide")) # for memory
    )
}

# duplicate the inner plan for each orientation
inner_dada_plan_multi <- tar_map(
  values = list(.orient = c("fwd", "rev")),
  names = .orient,
  inner_dada_plan
)

#### seqrun_plan ####

# the seqrun plan consists of steps that are run once per sequencing run.

seqrun_targets <- list(
  ##### errfun_{.seqrun} #####
  # function to use in dada or learnErrors
  #
  # only scans R2 because it is more likely to contain the lowest bin.
  errfun = tar_target(
    errfun,
    optimotu.pipeline::choose_dada_error_function(raw_R2),
    resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
  ),

  ##### denoise_read_counts_{.seqrun} #####
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

  ##### bimera_table_{.seqrun} #####
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
      !!(
        if (isTRUE(optimotu.pipeline::do_tag_jump())) quote(seqtable_uncross)
        else quote(seqtable_raw)
      ),
      seq_all,
      allowOneOff = TRUE,
      multithread = optimotu.pipeline::local_cpus()
    ),
    resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
  )
)

if (isTRUE(optimotu.pipeline::do_tag_jump())) {
  seqrun_targets <- c(
    seqrun_targets,
    list(
      ##### uncross_{.seqrun} #####
      # `tibble`:
      #   `sample (character) - sample name as given in sample_table$filt_key
      #   `nread` (integer) - number of reads in that sample (for this ASV)
      #   `total` (integer) - number of reads of that ASV across all samples (for
      #     this ASV)
      #   `uncross` (numeric) - UNCROSS score
      #   `is_tag_jump` (logical) - whether this occurrence is considered a likely
      #     tag jump
      #
      # `seq_idx` is not explicitly included; however the order of rows matches
      # `seqtable_raw`.
      #
      # remove tag-jumps (UNCROSS2)
      uncross = tar_fst_tbl(
        uncross,
        optimotu.pipeline::remove_tag_jumps(
          seqtable_raw, # raw_ASV_table
          !!optimotu.pipeline::tag_jump_f(), # f-value (expected cross-talk rate)
          !!optimotu.pipeline::tag_jump_p(), # p-value (power to rise the exponent)
          "seq_idx" # name of column which uniquely identifies the sequence
        ),
        resources = tar_resources(crew = tar_resources_crew(controller = "thin")),
        cue = tar_cue()
      ),

      ##### uncross_summary_{.seqrun} #####
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
        resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
      ),

      ##### uncross_read_counts_{.seqrun} #####
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
        resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
      ),
      ##### seqtable_uncross_{.seqrun} #####
      # `tibble`:
      #   `sample (character) - sample name as given in sample_table$sample_key
      #   `seq_idx` (integer) - index of a sequence in seq_all
      #   `nread` (integer) number of reads
      seqtable_uncross = tar_fst_tbl(
        seqtable_uncross,
        seqtable_raw[!uncross$is_tag_jump,],
        resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
      )
    )
  )
}

# There are a few seqrun targets which differ between forward, reverse, and "both" runs
seqrun_forward_targets <- c(
  seqrun_targets,
  list(
    ##### seq_merged_{.seqrun} #####
    # `character` vector
    #
    # all unique merged ASV sequences for each seqrun
    seq_merged = tar_target(
      seq_merged,
      unique(as.character(unlist(lapply(merged, \(x) x$sequence)))),
      resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
    )
  )
)

seqrun_reverse_targets <- c(
  seqrun_targets,
  list(
    ##### seq_merged_{.seqrun} #####
    # `character` vector
    #
    # all unique merged ASV sequences for each seqrun
    seq_merged = tar_target(
      seq_merged,
      unique(as.character(unlist(lapply(merged, \(x) dada2::rc(x$sequence))))),
      resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
    )
  )
)

seqrun_both_targets <- c(
  seqrun_targets,
  list(
    ##### seq_merged_{.seqrun} #####
    # `character` vector
    #
    # all unique merged ASV sequences for each seqrun
    seq_merged = tar_target(
      seq_merged,
      unique(as.character(unlist(c(
        lapply(merged_fwd, \(x) x$sequence),
        lapply(merged_rev, \(x) dada2::rc(x$sequence))
      )))),
      resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
    ),

    ##### seqtable_raw_{.seqrun} #####
    # `tibble`:
    #   `sample (character) - sample name as given in sample_table$sample_key
    #   `seq_idx` (integer) - index of a sequence in seq_all
    #   `nread` (integer) number of reads
    #
    # This combines seqtable_raw_fwd_{.seqrun} and seqtable_raw_rev_{.seqrun}
    seqtable_raw = tar_target(
      seqtable_raw,
      dplyr::bind_rows(seqtable_raw_fwd, seqtable_raw_rev) |>
        dplyr::summarize(nread = sum(nread), .by = c(sample, seq_idx)),
      resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
    ),

    ##### dada_map_{.seqrun} #####
    # `tibble`:
    #   `sample (character) - sample name as given in sample_table$sample_key
    #   `raw_idx` (integer) - index of read in the un-rarified fastq file
    #   `seq_idx` (integer) - index of ASV in seq_all
    #   `flags` (raw) - bits give presence/absence of the read after different stages:
    #    0x01: trim
    #    0x02: filter
    #    0x04: denoise & merge
    #    0x08: tag-jump removal (if performed)
    #
    # This combines dada_map_fwd_{.seqrun} and dada_map_rev_{.seqrun}
    #
    # If tag-jump removal is performed, it also adds the uncross information.
    dada_map = if (isTRUE(optimotu.pipeline::do_tag_jump())) {
      tar_fst_tbl(
        dada_map,
        optimotu.pipeline::merge_seq_maps(dada_map_fwd, dada_map_rev) |>
          optimotu.pipeline::add_uncross_to_seq_map(seqtable_raw, uncross),
        resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
      )
    } else {
      tar_fst_tbl(
        dada_map,
        optimotu.pipeline::merge_seq_maps(dada_map_fwd, dada_map_rev),
        resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
      )
    }
  )
)

# with "both" orientation we also need to consider both versions of raw_R2
seqrun_both_targets$errfun = tar_target(
  errfun,
  optimotu.pipeline::choose_dada_error_function(unique(c(raw_R2_fwd, raw_R2_rev))),
  resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
)


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
  inner_dada_plan_single,
  seqrun_forward_targets
)

# add "dummy" values if there were no forward-only seqruns
if (nrow(seqrun_forward_meta) == 0) {
  seqrun_forward_plan$dada2_meta <- list(tar_fst_tbl(
    dada2_meta_dummy_fwd,
    tibble::tibble(
      sample_key = character(),
      fastq_R1 = character(),
      fastq_R2 = character(),
      filt_R1 = character(),
      filt_R2 = character(),
      trim_R1 = character(),
      trim_R2 = character(),
    ),
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
  seqrun_forward_plan$filt_read_counts <- list(tar_fst_tbl(
    filt_read_counts_dummy_fwd,
    tibble::tibble(
      filt_R1 = character(),
      filt_nread = integer()
    ),
    deployment = "main"
  ))
}

# one row for each sequencing run which has _only_ reverse orientation
seqrun_reverse_meta <- dplyr::filter(seqrun_orient_meta, .orient == "rev") |>
  dplyr::select(.seqrun = seqrun, .orient)
seqrun_reverse_plan <- tar_map(
  values = seqrun_reverse_meta,
  names = .seqrun,
  inner_dada_plan_single,
  seqrun_reverse_targets
)

# add "dummy" values if there were no reverse-only seqruns
if (nrow(seqrun_reverse_meta) == 0) {
  seqrun_reverse_plan$dada2_meta <- list(tar_fst_tbl(
    dada2_meta_dummy_rev,
    tibble::tibble(
      sample_key = character(),
      fastq_R1 = character(),
      fastq_R2 = character(),
      filt_R1 = character(),
      filt_R2 = character(),
      trim_R1 = character(),
      trim_R2 = character(),
    ),
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
  seqrun_reverse_plan$filt_read_counts <- list(tar_fst_tbl(
    filt_read_counts_dummy_rev,
    tibble::tibble(
      filt_R1 = character(),
      filt_nread = integer()
    ),
    deployment = "main"
  ))
}

# one row for each sequencing run which has both orientations
seqrun_both_meta <- dplyr::filter(seqrun_orient_meta, .orient == "both") |>
  dplyr::select(.seqrun = seqrun)
seqrun_both_plan <- tar_map(
  values = seqrun_both_meta,
  names = .seqrun,
  inner_dada_plan_multi,
  seqrun_both_targets
)

# add "dummy" values if there were no both-orientation seqruns
if (nrow(seqrun_both_meta) == 0) {
  seqrun_both_plan$dada2_meta_fwd <- list(tar_fst_tbl(
    dada2_meta_fwd_dummy_both,
    tibble::tibble(
      sample_key = character(),
      fastq_R1 = character(),
      fastq_R2 = character(),
      filt_R1 = character(),
      filt_R2 = character(),
      trim_R1 = character(),
      trim_R2 = character()
    ),
    deployment = "main"
  ))
  seqrun_both_plan$dada2_meta_rev <- list(tar_fst_tbl(
    dada2_meta_rev_dummy_both,
    tibble::tibble(
      sample_key = character(),
      fastq_R1 = character(),
      fastq_R2 = character(),
      filt_R1 = character(),
      filt_R2 = character(),
      trim_R1 = character(),
      trim_R2 = character()
    ),
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
  seqrun_both_plan$filt_read_counts_fwd <- list(tar_fst_tbl(
    filt_read_counts_fwd_dummy_both,
    tibble::tibble(
      filt_R1 = character(),
      filt_nread = integer()
    ),
    deployment = "main"
  ))
  seqrun_both_plan$filt_read_counts_rev <- list(tar_fst_tbl(
    filt_read_counts_rev_dummy_both,
    tibble::tibble(
      filt_R1 = character(),
      filt_nread = integer()
    ),
    deployment = "main"
  ))
}

seqrun_plan <- optimotu.pipeline::tar_merge(
  seqrun_forward_plan,
  seqrun_reverse_plan
) |>
  optimotu.pipeline::tar_merge(seqrun_both_plan)

seq_all_file <- file.path(optimotu.pipeline::asv_path(), "all_asv.fasta.gz")

#### dada_plan ####
dada_plan <- c(
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
          Biostrings::DNAStringSet(uniqs[is.na(BiocGenerics::match(uniqs, old_seqs))])
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
        if (isTRUE(optimotu.pipeline::do_tag_jump())) "seqtable_uncross" else "seqtable_raw"
      ),
      resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
    ),

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

optimotu_plan <- c(optimotu_plan, dada_plan)
