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

dada_plan <- list(
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
    dplyr::group_by(sample_table, seqrun) %>%
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
  tar_file(
    trim,
    purrr::pmap(
      dplyr::transmute(
        dada2_meta,
        file_R1 = file.path(raw_path, fastq_R1),
        file_R2 = file.path(raw_path, fastq_R2),
        trim_R1 = trim_R1,
        trim_R2 = trim_R2
      ),
      cutadapt_paired_filter_trim_rc,
      max_err = pipeline_options$max_err, 
      min_overlap = pipeline_options$min_overlap, 
      truncQ_R1 = c(pipeline_options$truncQ_R1_f,pipeline_options$truncQ_R1_r), 
      truncQ_R2 = c(pipeline_options$truncQ_R2_f,pipeline_options$truncQ_R2_r), 
      max_n = pipeline_options$max_n, 
      min_length = pipeline_options$min_length,
      primer_R1 = pipeline_options$forward_primer, 
      primer_R2 = pipeline_options$reverse_primer, 
      cut_R2 = pipeline_options$cut_R2,
      action = pipeline_options$action, 
      discard_untrimmed = TRUE, #discard sequences that do not contain primers
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
    dplyr::group_by(sample_table_rc, seqrun) %>%
      tar_group(),
    iteration = "group",
    deployment = "main"
  ),

  tar_file(
    trim_rc,
    purrr::pmap(
      dplyr::transmute(
        dada2_meta_rc,
        file_R1 = file.path(raw_path, fastq_R1),
        file_R2 = file.path(raw_path, fastq_R2),
        trim_R1 = trim_R1,
        trim_R2 = trim_R2
      ),
      cutadapt_paired_filter_trim_rc,
      max_err = pipeline_options$max_err, 
      min_overlap = pipeline_options$min_overlap, 
      truncQ_R1 = c(pipeline_options$truncQ_R1_f,pipeline_options$truncQ_R1_r), 
      truncQ_R2 = c(pipeline_options$truncQ_R2_f,pipeline_options$truncQ_R2_r), 
      max_n = pipeline_options$max_n, 
      min_length = pipeline_options$min_length,
      primer_R1 = pipeline_options$reverse_primer, # reverse primer as 5’ primer in R1; for trimming primers in reverse complementary seqs 
      primer_R2 = pipeline_options$forward_primer, # forward primer as 5’ primer in R2; for trimming primers in reverse complementary seqs
      cut_R2 = pipeline_options$cut_R2,
      action = pipeline_options$action, 
      discard_untrimmed = TRUE, #discard sequences that do not contain primers
      ncpu = local_cpus(),
    ) %>%
      unlist(), 
    pattern = map(dada2_meta_rc), # per seqrun
    iteration = "list"
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

  #### filter_pairs ####
  # character: file names with path of filtered read files (fastq.gz)
  #
  # DADA2 quality filtering on read-pairs
  tar_target(
    dada2_meta_up,
    dplyr::group_by(rbind(sample_table, sample_table_rc), seqrun) %>%
      tar_group(),
    iteration = "group",
    deployment = "main"
  ),

  tar_file(
    filter_pairs,
    {
      file.create(c(dada2_meta_up$filt_R1, dada2_meta_up$filt_R2))
      dada2::filterAndTrim(
        fwd = purrr::keep(c(trim, trim_rc), endsWith, "_R1_trim.fastq.gz"),
        filt = dada2_meta_up$filt_R1,
        rev = purrr::keep(c(trim, trim_rc), endsWith, "_R2_trim.fastq.gz"),
        filt.rev = dada2_meta_up$filt_R2,
        #maxN = 0, # max 0 ambiguous bases (done in cutadapt)
        maxEE = c(pipeline_options$Fmaxee, pipeline_options$Rmaxee), # max expected errors (fwd, rev)
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
    tar_file(
      filtered,
      purrr::keep(filter_pairs, endsWith, paste0(read, "_filt.fastq.gz")),
      pattern = map(filter_pairs), # per seqrun × read
      iteration = "list",
      deployment = "main"
    ),
    tar_file(
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
    seqtable_raw,
    dada2::mergeSequenceTables(
      dada2::makeSequenceTable(
        dada2::mergePairs(
        denoise_R1, 
        derep_R1, 
        denoise_R2, 
        derep_R2, 
        minOverlap = 10, maxMismatch = 1, verbose = TRUE)
      ),
      dada2::makeSequenceTable(
        dada2::mergePairs(
        denoise_rc_R1, 
        derep_rc_R1, 
        denoise_rc_R2, 
        derep_rc_R2, 
        minOverlap = 10, maxMismatch = 1, verbose = TRUE)
      ),
      repeats = "error",
      tryRC = TRUE
    ),
    pattern = map(denoise_R1, derep_R1, denoise_R2, derep_R2, denoise_rc_R1, derep_rc_R1, denoise_rc_R2, derep_rc_R2), # per seqrun
    iteration = "list"
  ),

  #### remove tag-jumps (UNCROSS2) ####
  tar_target(
    seqtable_tagFilt,
    remove_tag_jumps(seqtable_raw, pipeline_options$f, pipeline_options$p),  #raw_ASV_table, f-value (expected cross-talk rate), p-value (power to rise the exponent)
    pattern = map(seqtable_raw) # per seqrun
  ),
 
  #### denoise_read_counts ####
  # tibble:
  #  `filt_key` character: as `sample_table$filt_key`
  #  `denoise_nread` integer: number of sequences in the sample after denoising
  tar_fst_tbl(
    denoise_read_counts,
    tibble::enframe(
      rowSums(seqtable_raw),
      name = "filt_key",
      value = "denoise_nread"
    ),
    pattern = map(seqtable_raw) # per seqrun
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
      seqtable_raw,
      allowOneOff=TRUE,
      multithread=local_cpus()
    ),
    pattern = map(seqtable_raw) # per seqrun
  ),
  
  #### seqtable_nochim ####
  # dada2 sequence table; integer matrix of read counts with column names as
  # sequences and row names as "samples" (i.e. sample_table$filt_key)
  #
  # combine sequence tables and remove consensus bimeras by combining
  # results from each seqrun.
  tar_target(
    seqtable_nochim,
    remove_bimera_denovo_tables(seqtable_tagFilt, bimera_table)
  )
  
  # #### nochim1_read_counts ####
  # # tibble:
  # #  `filt_key` character: as `sample_table$filt_key`
  # #  `nochim1_nread` integer: number of sequences in the sample after first
  # #    chimera filtering
  # tar_fst_tbl(
  #   nochim1_read_counts,
  #   tibble::enframe(
  #     rowSums(seqtable_nochim),
  #     name = "filt_key",
  #     value = "nochim1_nread"
  #   )
  # )
)
