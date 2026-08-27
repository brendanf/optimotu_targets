# Phase 1 analysis;
# Steps which are independant between reads, samples, and sequencing runs.
# Brendan Furneaux

# Based on DADA2 analysis for GSSP from Jenni Hultman
# edits by Sten Anslan - account for reverse complementary oriented sequences and add UNCROSS2 tag-jumps filtering per run

library(targets)
library(tarchetypes)

#### readwise_plan ####
# The "readwise" plan consists of targets which, in principle, could be run on
# each read pair individually. In practice they are still run on batches of
# files for each target, but the results for different read pairs within the
# same file are independent of one another. This means that these targets can be
# run prior to "raw read" rarefaction.
# These targets do need to be done separately for each orientation, if not all
# reads are in the same orientation. For convenience they are also done
# separately within each sequencing run.

readwise_plan <- c(
  list(
    ##### readwise_meta_{.orient?}_{.seqrun} #####
    # grouped tibble:
    #  `seqrun` character; name of sequencing run (directory in sequences/01_raw)
    #  `sample` character; name of sample, based on parsing file name
    #  `fastq_R1` character; file name with path for raw R1 file
    #  `fastq_R2` character; file name with path for raw R2 file
    #  `trim_R1` character; file name with path for trimmed R1 file
    #  `trim_R2` character; file name with path for trimmed R2 file
    #  `filt_R1` character; file name with path for filtered R1 file
    #  `filt_R2` character; file name with path for filtered R2 file
    #  `readwise_key`character; common prefix of trim_R1, trim_R2, filt_R1 and
    #      filt_R2
    readwise_meta = tar_group_size(
      readwise_meta,
      sample_table |>
        dplyr::filter(
          orient == .orient,
          seqrun == .seqrun
        ) |>
        dplyr::select(
          seqrun,
          sample,
          fastq_R1,
          fastq_R2,
          trim_R1,
          trim_R2,
          filt_R1,
          filt_R2,
          tidyselect::all_of("merged"),
          readwise_key,
          orient,
          any_of(optimotu.pipeline::cutadapt_paired_option_names),
          any_of(optimotu.pipeline::filter_option_names),
          any_of(optimotu.pipeline::merged_filter_option_names)
        ) |>
        dplyr::distinct() |>
        dplyr::arrange(readwise_key),
      size = 48,
      resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
    ),
    ##### raw_R1_{.orient?}_{.seqrun} #####
    # character: path and file name
    # raw reads, for dependency tracking
    raw_R1 = tar_file(
      raw_R1,
      unlist(strsplit(readwise_meta$fastq_R1, ",")),
      pattern = map(readwise_meta),
      resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
    ),
    ##### raw_R2_{.orient?}_{.seqrun} #####
    # character: path and file name
    # raw reads, for dependency tracking
    raw_R2 = tar_file(
      raw_R2,
      unlist(strsplit(readwise_meta$fastq_R2, ",")),
      pattern = map(readwise_meta),
      resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
    ),

    ##### raw_read_counts_{.orient?}_{.seqrun}_{.rarefaction?}_{.replicate?} #####
    # tibble:
    #  `fastq_file` character: file name of raw R2 file
    #  `raw_nread` integer: number of sequences in the file
    raw_read_counts = tar_fst_tbl(
      raw_read_counts,
      tibble::tibble(
        fastq_file = raw_R1,
        raw_nread = optimotu.pipeline::sequence_size(fastq_file)
      ),
      pattern = map(raw_R1),
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
        pairs_meta = readwise_meta,
        seqrun = .seqrun,
        orient = .orient,
        trim_options = !!optimotu.pipeline::trim_options(),
        primer_R1 = !!optimotu.pipeline::trim_primer_R1(),
        primer_R2 = !!optimotu.pipeline::trim_primer_R2(),
        cutadapt = !!optimotu.pipeline::find_cutadapt(),
        raw_R1 = raw_R1,
        raw_R2 = raw_R2
      ),
      pattern = map(readwise_meta, raw_R1, raw_R2),
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
    )
  ),

  if (optimotu.pipeline::do_dada2()) {
    list(
      ##### filter_pairs_{.orient?}_{.seqrun} #####
      # character: file names with path of filtered read files (fastq.gz)
      #
      # additional quality filtering on read-pairs
      filter_pairs = tar_file(
        filter_pairs,
        readwise_meta |>
          dplyr::mutate(
            trim_R1 = purrr::keep(trim, endsWith, "_R1_trim.fastq.gz"),
            trim_R2 = purrr::keep(trim, endsWith, "_R2_trim.fastq.gz")
          ) |>
          dplyr::group_by(
            dplyr::pick(any_of(c(
              "seqrun",
              optimotu.pipeline::filter_option_names
            )))
          ) |>
          dplyr::group_map(
            ~ optimotu.pipeline::filterAndTrim(
              fwd = .x$trim_R1,
              filt = .x$filt_R1,
              rev = .x$trim_R2,
              filt.rev = .x$filt_R2,
              maxEE = update(!!optimotu.pipeline::dada2_maxEE(), .y),
              rm.phix = TRUE,
              compress = TRUE,
              multithread = optimotu.pipeline::local_cpus(),
              verbose = TRUE
            )
          ) |>
          unlist(),
        pattern = map(readwise_meta, trim),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "wide")
        )
      ),

      ##### filt_read_counts_{.orient?}_{.seqrun} #####
      # tibble:
      #  `filt_R1` character: file name with path of filtered R1 file
      #  `filt_nread` integer: number of sequences in the file
      filt_read_counts = tar_fst_tbl(
        filt_read_counts,
        tibble::tibble(
          filt_R1 = purrr::keep(filter_pairs, endsWith, "_R1_filt.fastq.gz"),
          filt_nread = optimotu.pipeline::sequence_size(filt_R1)
        ),
        pattern = map(filter_pairs),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      )
    )
  },

  if (optimotu.pipeline::do_unoise()) {
    list(
      ##### premerge_seqs_{.orient?}_{.seqrun} #####
      # character: merged and quality-filtered FASTQ files
      #
      # Merge-then-filter for UNOISE (and future merge-first denoisers).
      premerge_seqs = tar_file(
        premerge_seqs,
        {
          meta <- readwise_meta |>
            dplyr::mutate(
              trim_R1 = purrr::keep(trim, endsWith, "_R1_trim.fastq.gz"),
              trim_R2 = purrr::keep(trim, endsWith, "_R2_trim.fastq.gz")
            )
          meta |>
            dplyr::group_by(
              dplyr::pick(
                any_of(c(
                  "seqrun",
                  optimotu.pipeline::merged_filter_option_names
                ))
              )
            ) |>
            dplyr::group_map(
              ~ optimotu.pipeline::vsearch_fastq_merge_pairs(
                seq_R1 = .x$trim_R1,
                seq_R2 = .x$trim_R2,
                seq_out = .x$merged,
                min_overlap = !!optimotu.pipeline::merge_min_overlap(),
                max_mismatch = !!optimotu.pipeline::merge_max_mismatch(),
                filter_options = stats::update(
                  !!optimotu.pipeline::merged_filter_options(),
                  .y
                ),
                threads = 1L,
                shards = optimotu.pipeline::local_cpus(),
                vsearch = !!optimotu.pipeline::find_vsearch()
              )
            ) |>
            unlist()
        },
        pattern = map(readwise_meta, trim),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "wide")
        )
      ),

      ##### merge_read_counts_{.orient?}_{.seqrun} #####
      # tibble:
      #  `merged` character: file name with path of merged FASTQ
      #  `merge_nread` integer: number of merged sequences in the file
      merge_read_counts = tar_fst_tbl(
        merge_read_counts,
        tibble::tibble(
          merged = premerge_seqs,
          merge_nread = optimotu.pipeline::sequence_size(merged)
        ),
        pattern = map(premerge_seqs),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      )
    )
  }
)
