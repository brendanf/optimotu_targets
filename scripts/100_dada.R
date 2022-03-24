# DADA2 quality filtering and denoising for Sonja's spruce log metabarcoding data
# Brendan Furneaux
# Based on DADA2 analysis for GSSP from Jenni Hultman

library(dada2)
library(tidyverse)
library(magrittr)
library(targets)
library(tarchetypes)

#### DADA2 analysis based on Jenni Hultman's pipeline for GSSP

dada_plan <- list(
  #### dada2_meta ####
  # grouping structure here will lead to separate dada2 error models
  tar_target(
    dada2_meta,
    dplyr::group_by(sample_table, "seqrun") %>%
      tar_group(),
    iteration = "group"
  ),
  #### fastq_r1 ####
  # declare a file target for dependency tracking
  tar_file(
    fastq_R1,
    file.path(raw_path, dada2_meta$fastq_R1),
    pattern = map(dada2_meta),
    iteration = "list"
  ),
  #### fastq_r2 ####
  # declare a file target for dependency tracking
  tar_file(
    fastq_R2,
    file.path(raw_path, dada2_meta$fastq_R2),
    pattern = map(dada2_meta),
    iteration = "list"
  ),

  #### all_filtered ####
  tar_file(
    all_filtered,
    {
      filterAndTrim(
        fwd = fastq_R1,
        filt = dada2_meta$filt_R1,
        rev = fastq_R2,
        filt.rev = dada2_meta$filt_R2,
        maxN = 0, # max 0 ambiguous bases
        maxEE = c(3, 5), # max expected errors (fwd, rev)
        truncQ = 2, # truncate at first base with quality <= 2
        rm.phix = TRUE, #remove matches to phiX genome
        minLen = 50, # remove reads < 50bp
        compress = TRUE, # write compressed files
        multithread = local_cpus(),
        verbose = TRUE
      )
      # return file names so targets know what was created
      c(dada2_meta$filt_R1, dada2_meta$filt_R2)
    },
    pattern = map(dada2_meta, fastq_R1, fastq_R2),
    iteration = "list"
  ),

  
  # inside the tar_map, every occurrence of read is replaced by "R1" or "R2"
  # the read name is also appended to all target names
  # so all of this gets done separately for forward and reverse reads.
  # pattern=map() means we are also keeping the different sequencing runs
  # separate
  tar_map(
    values = list(read = c("R1", "R2")),
    #### filtered_{read} ####
    tar_file(
      filtered,
      purrr::keep(all_filtered, endsWith, paste0(read, "_filt.fastq.gz")),
      pattern = map(all_filtered),
      iteration = "list"
    ),
    
    #### derep_{read} ####
    tar_target(
      derep,
      derepFastq(filtered, verbose = TRUE) %>%
        set_names(dada2_meta$sample),
      pattern = map(dada2_meta, filtered),
      iteration = "list"
    ),
    
    #### err_{read} ####
    tar_target(
      err,
      learnErrors(filtered, multithread = local_cpus(), verbose = TRUE),
      pattern = map(filtered),
      iteration = "list"
    ),
    
    #### denoise_{read} ####
    tar_target(
      denoise,
      dada(derep, err = err, multithread = local_cpus(), verbose = TRUE),
      pattern = map(derep, err),
      iteration = "list"
    )
  ),
  #### merged ####
  # Merge paired reads and make a sequence table for each sequencing run
  tar_target(
    merged,
    mergePairs(denoise_R1, derep_R1, denoise_R2, derep_R2,
               minOverlap = 10, maxMismatch = 1, verbose=TRUE),
    pattern = map(denoise_R1, derep_R1, denoise_R2, derep_R2),
    iteration = "list"
  ),
  
  #### seqtable_dup ####
  # Make sequence table for each sequencing run
  # these may contain some sequences which are no-mismatch pairs, i.e. only
  # differ by length
  tar_target(
    seqtable_dup,
    makeSequenceTable(merged),
    pattern = map(merged),
    iteration = "list"
  ),
  
  #### seqtable ####
  # Merge no-mismatch pairs
  tar_target(
    seqtable,
    collapseNoMismatch(seqtable_dup, minOverlap = 50, verbose = TRUE),
    pattern = map(seqtable_dup),
    iteration = "list"
  ),
  
  #### asvtable_dup ####
  # combine all the remaining sequence tables into a master table
  # again, it's possible we have no-mismatch pairs between the different
  # sequencing runs.
  tar_target(
    asvtable_dup,
    mergeSequenceTables(tables = seqtable),
  ),
  
  #### asvtable ####
  # Merge no-mismatch pairs
  tar_target(
    asvtable,
    collapseNoMismatch(asvtable_dup)
  ),
  
  #### write_asvtable ####
  tar_file(
    write_asvtable,
    file.path(asv_path, "asv_tab.rds") %T>%
      saveRDS(asvtable, .)
  )
)
