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

orient_meta <- tibble::tibble(
  .orient = unique(sample_table$orient)
)

#### inner_dada_plan ####
# the "inner" dada plan consists of targets which need to be done separately
# for each orientation, if not all reads are in the same orientation.

inner_dada_plan <- list(
  ##### dada2_meta_{.orient}_{.seqrun} #####
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
  # `sample_table` is defined in scripts/010_load_samples.R
  tar_fst_tbl(
    dada2_meta,
    sample_table |>
      dplyr::filter(orient == .orient, seqrun == .seqrun) |>
      dplyr::select(seqrun, sample, fastq_R1, fastq_R2, trim_R1, trim_R2,
                    filt_R1, filt_R2, filt_key, any_of(cutadapt_option_names))
  ),

  ##### raw_R1_{.orient}_{.seqrun} #####
  # character: path and file name
  # raw reads, for dependency tracking
  tar_file_fast(
    raw_R1,
    file.path(raw_path, dada2_meta$fastq_R1)
  ),

  ##### raw_read_counts_{.orient}_{.seqrun} #####
  # tibble:
  #  `fastq_file` character: file name of raw R1 file
  #  `raw_nread` integer: number of sequences in the file
  tar_fst_tbl(
    raw_read_counts,
    tibble::tibble(
      fastq_file = raw_R1,
      raw_nread = sequence_size(fastq_file)
    )
  ),


  ##### trim_{.orient}_{.seqrun} #####
  # character: file names with path of trimmed read files (fastq.gz)
  #
  # remove adapters and barcodes
  # also do some preliminary quality filtering
  tar_file_fast(
    trim,
    dplyr::group_by(
      dada2_meta,
      dplyr::pick(any_of(c("seqrun", cutadapt_option_names)))
    ) |>
    dplyr::group_map(
      ~ purrr::pmap(
        dplyr::transmute(
          .x,
          file_R1 = file.path(raw_path, fastq_R1),
          file_R2 = file.path(raw_path, fastq_R2),
          trim_R1 = trim_R1,
          trim_R2 = trim_R2
        ),
        cutadapt_paired_filter_trim,
        primer_R1 = ifelse(.orient == "fwd", trim_primer_R1, trim_primer_R2),
        primer_R2 = ifelse(.orient == "fwd", trim_primer_R2, trim_primer_R1),
        options = update(trim_options, .x),
        ncpu = local_cpus(),
      ) |>
        unlist(),
      .keep = TRUE
    ) |>
      unlist()
  ),

  ##### trim_read_counts_{.orient}_{.seqrun} #####
  # tibble:
  #  `trim_R1` character: file name with path of trimmed R1 file
  #  `trim_nread` integer: number of sequences in the file
  #
  # count of reads per sample after adapter trimming
  tar_fst_tbl(
    trim_read_counts,
    tibble::tibble(
      trim_R1 = purrr::keep(trim, endsWith, "_R1_trim.fastq.gz"),
      trim_nread = sequence_size(trim_R1)
    )
  ),

  # DADA2 quality filtering on read-pairs

  ##### filter_pairs_{.orient}_{.seqrun} #####
  # character: file names with path of filtered read files (fastq.gz)
  #
  # additional quality filtering on read-pairs
  tar_file_fast(
    filter_pairs,
    if (nrow(dada2_meta) > 0L) {
      file.create(c(dada2_meta$filt_R1, dada2_meta$filt_R2))
      dada2::filterAndTrim(
        fwd = purrr::keep(trim, endsWith, "_R1_trim.fastq.gz"),
        filt = dada2_meta$filt_R1,
        rev = purrr::keep(trim, endsWith, "_R2_trim.fastq.gz"),
        filt.rev = dada2_meta$filt_R2,
        maxEE = dada2_maxEE, # max expected errors (fwd, rev)
        rm.phix = TRUE, #remove matches to phiX genome
        compress = TRUE, # write compressed files
        multithread = local_cpus(),
        verbose = TRUE
      )
      # return file names for samples where at least some reads passed
      c(dada2_meta$filt_R1, dada2_meta$filt_R2) %>%
        purrr::keep(file.exists)
    } else {
      character()
    }
  ),

  ##### filt_read_counts_{.orient}_{.seqrun} #####
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
    )
  ),

  ##### map over R1 and R2 #####
  # inside the tar_map, every occurrence of `read` is replaced by "R1" or "R2"
  # the read name is also appended to all target names
  # so all of this gets done separately for forward and reverse reads.
  # pattern=map() means we are also keeping the different sequencing runs
  # separate
  tar_map(
    values = list(read = c("R1", "R2")),

    ###### filtered_{read}_{.orient}_{.seqrun} ######
    # character: path and file name of filtered reads; fastq.gz
    #
    # select only the files corresponding to the read we are working on
    tar_file_fast(
      filtered,
      purrr::keep(filter_pairs, endsWith, paste0(read, "_filt.fastq.gz")),
      deployment = "main"
    ),

    ###### derep_{read}_{.orient}_{.seqrun} ######
    # list of dada2 `derep` objects
    #
    # dereplicate
    tar_target(
      derep,
      dada2::derepFastq(filtered, verbose = TRUE) %>%
        set_names(sub("_(fwd|rev)_R[12]_filt\\.fastq\\.gz", "", filtered)),
      iteration = "list"
    ),

    ###### err_{read}_{.orient}_{.seqrun} ######
    # list: see dada2::LearnErrors
    #
    # fit error profile
    tar_target(
      err,
      if (sum(!grepl("BLANK|NEG", filtered)) > 0L) {
        dada2::learnErrors(
          purrr::discard(filtered, grepl, pattern = "BLANK|NEG"),
          multithread = local_cpus(),
          verbose = TRUE
        )
      } else {
        NULL
      }
    ),

    ###### denoise_{read}_{.orient}_{.seqrun} ######
    # list of dada2 `dada` objects
    tar_target(
      denoise,
      if (is.null(err)) {
        NULL
      } else {
        dada2::dada(derep, err = err, multithread = local_cpus(), verbose = TRUE)
      }
    )
  ),

  ##### merged_{.orient}_{.seqrun} #####
  # list of data.frame; see dada2::mergePairs
  #
  # Merge paired reads and make a sequence table for each sequencing run
  tar_target(
    merged,
    if (is.null(denoise_R1) | is.null(denoise_R2)) {
      list()
    } else {
      dada2::mergePairs(
        denoise_R1,
        derep_R1,
        denoise_R2,
        derep_R2,
        minOverlap = 10,
        maxMismatch = 1,
        verbose=TRUE
      )
    }
  ),

  ##### dada_map_{.orient}_{.seqrun} #####
  # map the raw reads to the nochim ASVs
  # indexes are per-sample
  #
  # a tibble:
  #  sample: character identifying the sample, as in sample_table
  #  read_in_sample: integer index of read in the un-rarified fastq file
  #  flags: raw, bits give presence/absence of the read after different stages:
  #    0x01: trim
  #    0x02: filter
  #    0x04: denoise
  #    0x08: chimera check
  #  nochim_id: integer index of ASV in columns of seqtable_nochim

  tar_target(
    dada_map,
    mapply(
      FUN = nochim_map,
      sample = dada2_meta$sample,
      fq_raw = file.path(raw_path, dada2_meta$fastq_R1),
      fq_trim = dada2_meta$trim_R1,
      fq_filt = dada2_meta$filt_R1,
      dadaF = denoise_R1,
      derepF = derep_R1,
      dadaR = denoise_R2,
      derepR = derep_R2,
      merged = merged,
      MoreArgs = list(
        seqtable_nochim = seqtable_nochim,
        orient = .orient
      ),
      SIMPLIFY = FALSE
    ) |>
      purrr::list_rbind()
  )
)

if (nrow(orient_meta) > 1L) {
  inner_dada_plan <- tar_map(
    values = orient_meta,
    inner_dada_plan
  )
} else {
  .orient <- orient_meta$.orient
}

seqrun_meta <- tibble::tibble(
  .seqrun = unique(sample_table$seqrun)
)

seqrun_plan <- tar_map(
  values = seqrun_meta,
  inner_dada_plan,

  #### seqtable_raw_{.seqrun} ####
  # dada2 sequence table; integer matrix of read counts with column names as
  # sequences and row names as "samples" (i.e. sample_table$filt_key)
  #
  # Make sequence table for each sequencing run
  # these may contain some sequences which are no-mismatch pairs, i.e. only
  # differ by length
  if (nrow(orient_meta) == 2) {
    tar_target(
      seqtable_raw,
      if (length(merged_fwd) > 0 && length(merged_rev) > 0) {
      dada2::mergeSequenceTables(
        dada2::makeSequenceTable(merged_fwd),
        dada2::makeSequenceTable(merged_rev) |>
          (\(x) magrittr::set_colnames(x, dada2::rc(colnames(x))))(),
        repeats = "sum"
      )
      } else if (length(merged_fwd) > 0) {
        dada2::makeSequenceTable(merged_fwd)
      } else if (length(merged_rev) > 0) {
        dada2::makeSequenceTable(merged_rev) |>
          (\(x) magrittr::set_colnames(x, dada2::rc(colnames(x))))()
      } else {
        dada2::makeSequenceTable(list())
      }
    )
  } else {
    tar_target(
      seqtable_raw,
      dada2::makeSequenceTable(merged)
    )
  },

  #### denoise_read_counts_{.seqrun} ####
  # tibble:
  #  `filt_key` character: as `sample_table$filt_key`
  #  `denoise_nread` integer: number of sequences in the sample after denoising
  tar_fst_tbl(
    denoise_read_counts,
    tibble::enframe(
      rowSums(seqtable_raw),
      name = "filt_key",
      value = "denoise_nread"
    )
  ),

  #### bimera_table_{.seqrun} ####
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
      seqtable_raw,
      allowOneOff=TRUE,
      multithread=local_cpus()
    )
  )
)

dada_plan <- list(
  seqrun_plan,

  #### seqtable_nochim ####
  # dada2 sequence table; integer matrix of read counts with column names as
  # sequences and row names as "samples" (i.e. sample_table$filt_key)
  #
  # combine sequence tables and remove consensus bimeras by combining
  # results from each seqrun.
  tar_target(
    seqtable_nochim,
    remove_bimera_denovo_tables(
        seqtabs = !!tar_map_list(seqrun_plan$seqtable_raw),
        bimdf = !!tar_map_bind_rows(seqrun_plan$bimera_table)
    )
  ),

  #### nochim1_read_counts ####
  # tibble:
  #  `filt_key` character: as `sample_table$filt_key`
  #  `nochim1_nread` integer: number of sequences in the sample after first
  #    chimera filtering
  tar_fst_tbl(
    nochim1_read_counts,
    tibble::enframe(
      rowSums(seqtable_nochim),
      name = "filt_key",
      value = "nochim1_nread"
    )
  )
)
