# DADA2 quality filtering and denoising for Sonja's spruce log metabarcoding data
# Brendan Furneaux
# Based on DADA2 analysis for GSSP from Jenni Hultman
  # edits by Sten Anslan - account for reverse complementary oriented sequences and add UNCROSS2 tag-jumps filtering per run

############################################################
## TODO: ##
# - skip RC if no rc seqs! Currently = ERROR, if no rc seqs.
# - remove empty files after cutadapt in 02_trim/
############################################################

library(magrittr)
library(targets)
library(tarchetypes)

#### DADA2 analysis based on Jenni Hultman's pipeline for GSSP

dada_plan_mixed <- list(
  #### dada2_meta ####
  # grouped tibble:
  #  `seqrun` character; name of sequencing run (directory in sequences/01_raw)
  #  `sample` character; name of sample, based on parsing file name
  #  `fastq_R1` character; file name with path for raw R1 file
  #  `fastq_R2` character; file name with path for raw R2 file
  #  `trim_R1` character; file name with path for trimmed R1 file
  #  `trim_R2` character; file name with path for trimmed R2 file
  #  `filt_R1` character; file name with path for filtered R1 file
  #  `filt_R2` character; file name with path for filtered R2 file
  #  `filt_key`character; common prefix of filt_R1 and filt_R2; used as sample
  #      name by dada2 functions and read counts
  #
  # grouping structure here will lead to separate dada2 error models
  # `sample_table` is defined in scripts/010_load_samples.R
  tar_target(
    dada2_meta,
    sample_table |>
      dplyr::select(seqrun, sample, fastq_R1, fastq_R2, trim_R1, trim_R2,
                    filt_R1, filt_R2, filt_key) |>
      dplyr::group_by(seqrun) |>
      tar_group(),
    iteration = "group",
    deployment = "main"
  ),

  #### raw_read_counts ####
  # tibble:
  #  `fastq_file` character: file name of raw R1 file
  #  `raw_nread` integer: number of sequences in the file
  tar_fst_tbl(
    raw_read_counts,
    tibble::tibble(
      fastq_file = file.path(raw_path, dada2_meta$fastq_R1),
      raw_nread = sequence_size(fastq_file)
    ),
    pattern = map(dada2_meta) # per seqrun
  ),

  #### trim, round 1 - clip fwd primer in R1 and rev primer in R2 ####
  # character: file names with path of trimmed read files (fastq.gz)
  #
  # remove adapters and barcodes
  # also do some preliminary quality filtering
  tar_file_fast(
    trim,
    purrr::pmap(
      dplyr::transmute(
        dada2_meta,
        file_R1 = file.path(raw_path, fastq_R1),
        file_R2 = file.path(raw_path, fastq_R2),
        trim_R1 = trim_R1,
        trim_R2 = trim_R2
      ),
      cutadapt_paired_filter_trim,
      primer_R1 = trim_primer_R1,
      primer_R2 = trim_primer_R2,
      options = trim_options,
      ncpu = local_cpus(),
    ) %>%
      unlist(),
    pattern = map(dada2_meta), # per seqrun
    iteration = "list"
  ),

  #### trim, round 2 - clip rev primer in R1 and fwd primer in R2 ####
     # needed when seqs are in revcomp orientation
  tar_target(
    dada2_meta_rc,
    sample_table_rc |>
      dplyr::select(seqrun, sample, fastq_R1, fastq_R2, trim_R1, trim_R2,
                    filt_R1, filt_R2, filt_key) |>
      dplyr::group_by(seqrun) |>
      tar_group(),
    iteration = "group",
    deployment = "main"
  ),

  tar_file_fast(
    trim_rc,
    purrr::pmap(
      dplyr::transmute(
        dada2_meta_rc,
        file_R1 = file.path(raw_path, fastq_R1),
        file_R2 = file.path(raw_path, fastq_R2),
        trim_R1 = trim_R1,
        trim_R2 = trim_R2
      ),
      cutadapt_paired_filter_trim,
      primer_R1 = trim_primer_R2, # reverse primer as 5’ primer in R1; for trimming primers in reverse complementary seqs
      primer_R2 = trim_primer_R1,
      options = trim_options,
      ncpu = local_cpus(),
    ) %>%
      unlist(),
    pattern = map(dada2_meta_rc), # per seqrun
    iteration = "list",
    cue = tar_cue(format = FALSE)
  ),

  #### trim_read_counts ####
  # tibble:
  #  `trim_R1` character: file name with path of trimmed R1 file
  #  `trim_nread` integer: number of sequences in the file
  #
  # count of reads per sample after adapter trimming
  tar_fst_tbl(
    trim_read_counts,
    tibble::tibble(
      trim_R1 = purrr::keep(c(trim, trim_rc), endsWith, "_R1_trim.fastq.gz"),
      trim_nread = sequence_size(trim_R1)
    ),
    pattern = map(trim, trim_rc)
  ),
  tar_target(
    dada2_meta_up,
    rbind(dada2_meta, dada2_meta_rc),
    iteration = "group",
    deployment = "main"
  ),

  #### filter_pairs ####
  # character: file names with path of filtered read files (fastq.gz)
  #
  # DADA2 quality filtering on read-pairs
  tar_file_fast(
    filter_pairs,
    {
      file.create(c(dada2_meta_up$filt_R1, dada2_meta_up$filt_R2))
      dada2::filterAndTrim(
        fwd = purrr::keep(c(trim, trim_rc), endsWith, "_R1_trim.fastq.gz"),
        filt = dada2_meta_up$filt_R1,
        rev = purrr::keep(c(trim, trim_rc), endsWith, "_R2_trim.fastq.gz"),
        filt.rev = dada2_meta_up$filt_R2,
        #maxN = 0, # max 0 ambiguous bases (done in cutadapt)
        maxEE = dada2_maxEE, # max expected errors (fwd, rev)
        #truncQ = 2, # truncate at first base with quality <= 2 (done in cutadapt)
        rm.phix = TRUE, #remove matches to phiX genome
        #minLen = 100, # remove reads < 100bp (done by cutadapt)
        compress = TRUE, # write compressed files
        multithread = local_cpus(),
        verbose = TRUE
      )
      # return file names for samples where at least some reads passed
      c(dada2_meta_up$filt_R1, dada2_meta_up$filt_R2) %>%
        purrr::keep(file.exists)
    },
    pattern = map(trim, trim_rc, dada2_meta_up), # per seqrun
    iteration = "list"
  ),

  #### filt_read_counts ####
  # tibble:
  #  `filt_R1` character: file name with path of filtered R1 file
  #  `filt_nread` integer: number of sequences in the file
  #
  # count of reads after filtering
  tar_fst_tbl(
    filt_read_counts,
    tibble::tibble(
      filt_R1 = purrr::keep(filter_pairs, endsWith, "_R1_filt.fastq.gz"),
      filt_nread = sequence_size(filt_R1)
    ),
    pattern = map(filter_pairs) # per seqrun
  ),

  # inside the tar_map, every occurrence of `read` is replaced by "R1" or "R2"
  # the read name is also appended to all target names
  # so all of this gets done separately for forward and reverse reads.
  # pattern=map() means we are also keeping the different sequencing runs
  # separate
  tar_map(
    values = list(read = c("R1", "R2")),

    #### filtered_{read} ####
    # character: path and file name of filtered reads; fastq.gz
    #
    # select only the files corresponding to the read we are working on
    tar_file_fast(
      filtered,
      purrr::keep(filter_pairs, endsWith, paste0(read, "_filt.fastq.gz")),
      pattern = map(filter_pairs), # per seqrun × read
      iteration = "list",
      deployment = "main"
    ),
    tar_file_fast(
      filtered_rc,
      purrr::keep(filter_pairs, endsWith, paste0(read, "_filt_rc.fastq.gz")),
      pattern = map(filter_pairs), # per seqrun × read
      iteration = "list",
      deployment = "main"
    ),

    #### derep_{read} ####
    # list of dada2 `derep` objects
    #
    # dereplicate
    tar_target(
      derep,
      dada2::derepFastq(filtered, verbose = TRUE) %>%
        set_names(sub("_R[12]_filt\\.fastq\\.gz", "", filtered)),
      pattern = map(filtered), # per seqrun × read
      iteration = "list"
    ),
    # dereplicate files that have revcomp seqs
    tar_target(
      derep_rc,
      dada2::derepFastq(filtered_rc, verbose = TRUE) %>%
        set_names(sub("_R[12]_filt_rc\\.fastq\\.gz", "", filtered_rc)),
      pattern = map(filtered_rc), # per seqrun × read
      iteration = "list"
    ),

    #### err_{read} ####
    # list: see dada2::LearnErrors
    #
    # fit error profile
    tar_target(
      err,
      dada2::learnErrors(
        purrr::discard(filtered, grepl, pattern = "BLANK|NEG"),
        multithread = local_cpus(),
        verbose = TRUE
      ),
      pattern = map(filtered), # per seqrun × read
      iteration = "list"
    ),
    tar_target(
      err_rc,
      dada2::learnErrors(
        purrr::discard(filtered_rc, grepl, pattern = "BLANK|NEG"),
        multithread = local_cpus(),
        verbose = TRUE
      ),
      pattern = map(filtered_rc), # per seqrun × read
      iteration = "list"
    ),

    #### denoise_{read} ####
    # list of dada2 `dada` objects
    tar_target(
      denoise,
      dada2::dada(derep, err = err, multithread = local_cpus(), verbose = TRUE),
      pattern = map(derep, err), # per seqrun × read
      iteration = "list"
    ),
    tar_target(
      denoise_rc,
      dada2::dada(derep_rc, err = err_rc, multithread = local_cpus(), verbose = TRUE),
      pattern = map(derep_rc, err_rc), # per seqrun × read
      iteration = "list"
    )
  ),

  #### merged ####
  # list of data.frame; see dada2::mergePairs
  #
  # Merge paired reads and make a sequence table for each sequencing run, and
  # make dada2 sequence table for each sequencing run; integer matrix of read counts with column names as
  # sequences and row names as "samples" (i.e. sample_table$filt_key).

  tar_target(
    merged,
    dada2::mergePairs(
      denoise_R1,
      derep_R1,
      denoise_R2,
      derep_R2,
      minOverlap = 10,
      maxMismatch = 1,
      verbose = TRUE
    ),
    pattern = map(denoise_R1, derep_R1, denoise_R2, derep_R2),
    iteration  = "list"
  ),

  tar_target(
    merged_rc,
    dada2::mergePairs(
      denoise_rc_R1,
      derep_rc_R1,
      denoise_rc_R2,
      derep_rc_R2,
      minOverlap = 10,
      maxMismatch = 1,
      verbose = TRUE
    ),
    pattern = map(denoise_rc_R1, derep_rc_R1, denoise_rc_R2, derep_rc_R2),
    iteration  = "list"
  ),

  #### seq_merged ####
  # `character` vector
  #
  # all unique merged ASV sequences for each seqrun
  tar_target(
    seq_merged,
    unique(
      unlist(
        c(
          lapply(merged, \(x) x$sequence),
          lapply(merged_rc, \(x) dada2::rc(x$sequence))
          )
      )
    ),
    pattern = map(merged, merged_rc) # per seqrun
  ),

  #### seq_all ####
  # `character` vector
  #
  # all unique ASV sequences, across all seqruns
  tar_target(
    seq_all,
    unique(seq_merged)
  ),

  #### seqtable_raw_fwd ####
  # `tibble` with columns:
  #   `sample` (character) sample name as given in sample_table$filt_key
  #   `seq_idx` (integer) index of a sequence in seq_all
  #   `nread` (integer) number of reads
  tar_fst_tbl(
    seqtable_raw_fwd,
    make_mapped_sequence_table(merged, seq_all, rc = FALSE),
    pattern = map(merged), # per seqrun
    iteration = "list"
  ),

  #### seqtable_raw_rc ####
  # `tibble` with columns:
  #   `sample` (character) sample name as given in sample_table$filt_key
  #   `seq_idx` (integer) index of a sequence in seq_all
  #   `nread` (integer) number of reads
  tar_fst_tbl(
    seqtable_raw_rc,
    make_mapped_sequence_table(merged_rc, seq_all, rc = TRUE),
    pattern = map(merged_rc), # per seqrun
    iteration = "list"
  ),

  #### seqtable_raw ####
  # `tibble`:
  #   `sample (character) - sample name as given in sample_table$filt_key
  #   `seq_idx` (integer) - index of a sequence in seq_all
  #   `nread` (integer) number of reads
  tar_fst_tbl(
    seqtable_raw,
    rbind(seqtable_raw_fwd, seqtable_raw_rc) |>
      dplyr::summarize(nread = sum(nread), .by = c(sample, seq_idx)),
    pattern = map(seqtable_raw_fwd, seqtable_raw_rc) # per seqrun
  ),

  #### denoise_read_counts ####
  # tibble:
  #  `filt_key` character: as `sample_table$filt_key`
  #  `denoise_nread` integer: number of sequences in the sample after denoising
  tar_fst_tbl(
    denoise_read_counts,
    dplyr::summarize(
      seqtable_raw,
      denoise_nread = sum(nread),
      .by = sample
    ) |>
      dplyr::rename(filt_key = sample),
    pattern = map(seqtable_raw) # per seqrun
  ),

  #### uncross ####
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
  tar_fst_tbl(
    uncross,
    remove_tag_jumps(
      seqtable_raw, # raw_ASV_table
      pipeline_options$f, # f-value (expected cross-talk rate)
      pipeline_options$p, # p-value (power to rise the exponent)
      "seq_idx" # name of column which uniquely identifies the sequence
    ),
    pattern = map(seqtable_raw) # per seqrun
  ),

  #### uncross_summary ####
  # `tibble`:
  #   `sample (character) - sample name as given in sample_table$filt_key
  #   `Total_reads` (integer) - total reads in the sample (all ASVs)
  #   `Number_of_TagJump_Events` (integer) - number of ASVs in the sample which
  #     are considered tag jumps.
  #   `TagJump_reads` (integer) - number of reads in the sample which belong to
  #     tag-jump ASVs.
  tar_fst_tbl(
    uncross_summary,
    summarize_uncross(uncross),
    pattern = map(uncross),
    iteration = "list"
  ),

  #### uncross_read_counts ####
  # `tibble`:
  #   `filt_key` character: as `sample_table$filt_key`
  #   `uncross_nread` integer: number of sequences in the sample after uncrossing
  tar_fst_tbl(
    uncross_read_counts,
    uncross_summary |>
      dplyr::transmute(
        filt_key = sample,
        uncross_nread = Total_reads - TagJump_reads
      ),
    pattern = map(uncross_summary) # per seqrun
  ),

  #### bimera_table ####
  # tibble:
  #  `nflag` integer: number of samples in which the sequence was considered
  #    chimeric
  #  `nsam` integer: number of samples in which the sequence occurred
  #  `seq` character: sequence
  #
  # find denovo chimeric sequences in each sample independently
  tar_fst_tbl(
    bimera_table,
    bimera_denovo_table(
      seqtable_raw[!uncross$is_tag_jump,],
      seq_all,
      allowOneOff=TRUE,
      multithread=local_cpus()
    ),
    pattern = map(seqtable_raw, uncross) # per seqrun
  ),

  #### seq_nochim ####
  # `character` vector - ASV sequences which were found not to be chimeric
  #
  # calculate consensus chimera calls across all seqruns
  tar_target(
    seq_nochim,
    combine_bimera_denovo_tables(bimera_table)
  ),

  #### seqtable_nochim ####
  # long sequence table; `data.frame`integer matrix of read counts with column names as
  # sequences and row names as "samples" (i.e. sample_table$filt_key)
  #
  # combine sequence tables and remove consensus bimeras by combining
  # results from each seqrun.
  tar_target(
    seqtable_nochim,
    dplyr::filter(seqtable_raw, !uncross$is_tag_jump) |>
      dplyr::filter(seq_idx %in% (seq_all %in% seq_nochim))
  ),

  #### dada_map ####
  # map the raw reads to seq_all
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

  tar_target(
    dada_map,
    list(
      fwd = mapply(
        FUN = seq_map,
        sample = dada2_meta$sample,
        fq_raw = file.path(raw_path, dada2_meta$fastq_R1),
        fq_trim = dada2_meta$trim_R1,
        fq_filt = dada2_meta$filt_R1,
        dadaF = denoise_R1,
        derepF = derep_R1,
        dadaR = denoise_R2,
        derepR = derep_R2,
        merged = merged,
        MoreArgs = list(seq_all = seq_all),
        SIMPLIFY = FALSE
      ),
      rev = mapply(
        FUN = seq_map,
        sample = dada2_meta_rc$sample,
        fq_raw = file.path(raw_path, dada2_meta_rc$fastq_R1),
        fq_trim = dada2_meta_rc$trim_R1,
        fq_filt = dada2_meta_rc$filt_R1,
        dadaF = denoise_rc_R1,
        derepF = derep_rc_R1,
        dadaR = denoise_rc_R2,
        derepR = derep_rc_R2,
        merged = merged_rc,
        MoreArgs = list(seq_all = seq_all, rc = TRUE),
        SIMPLIFY = FALSE
      )
    ) |>
      purrr::pmap_dfr(
        .f = \(fwd, rev) dplyr::full_join(
          fwd,
          rev,
          by = c("sample", "raw_idx"),
          suffix = c("_fwd", "_rev")
        ) |>
          dplyr::transmute(
            seq_idx = dplyr::coalesce(seq_idx_fwd, seq_idx_rev),
            flags = flags_fwd | flags_rev
          )
      ),
    pattern = map(dada2_meta, denoise_R1, derep_R1, denoise_R2, derep_R2, merged,
                  dada2_meta_rc, denoise_rc_R1, derep_rc_R1, denoise_rc_R2,
                  derep_rc_R2, merged_rc)
  ),
  #### nochim1_read_counts ####
  # tibble:
  #  `filt_key` character: as `sample_table$filt_key`
  #  `nochim1_nread` integer: number of sequences in the sample after first
  #    chimera filtering
  tar_fst_tbl(
    nochim1_read_counts,
    dplyr::summarize(seqtable_nochim, nochim1_nread = sum(nread), .by = sample) |>
      dplyr::rename(filt_key = sample)
  )
)
