# DADA2 quality filtering and denoising for Sonja's spruce log metabarcoding data
# Brendan Furneaux
# Based on DADA2 analysis for GSSP from Jenni Hultman

library(dada2)
library(tidyverse)
library(magrittr)
library(targets)
library(tarchetypes)

#### DADA2 analysis based on Jenni Hultman's pipeline for GSSP

list(
  #### dada2_meta ####
  # grouping structure hear will lead to separate dada2 error models
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
    file_path(raw_path, dada2_meta$fastq_R1),
    pattern = map(dada2_meta),
    iteration = "list"
  ),
  #### fastq_r2 ####
  # declare a file target for dependency tracking
  tar_file(
    fastq_R2,
    file_path(raw_path, dada2_meta$fastq_R2),
    pattern = map(dada2_meta),
    iteration = "list"
  ),

  #### filtered ####
  tar_file(
    filtered,
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
        multithread = ncpus,
        verbose = TRUE
      )
      # return file names so targets know what was created
      c(dada2_meta$filt_R1, dada2_meta$filt_R2)
    },
    pattern = map(dada2_meta, fastq_R1, fastq_R2),
    iteration = "list"
  ),

#### Fit error model ####
errR1 <- sample_table %$%
  learnErrors(file.path(filt_path, filt_R1), multithread = TRUE)
errR2 <- sample_table %$%
  learnErrors(file.path(filt_path, filt_R2), multithread = TRUE)

#### Dereplicate reads ####
derepR1 <- sample_table %$%
  derepFastq(file.path(filt_path, filt_R1), verbose = TRUE)
names(derepR1) <- sample_table$sample

derepR2 <- sample_table %$%
  derepFastq(file.path(filt_path, filt_R2), verbose = TRUE)
names(derepR2) <- sample_table$sample

#### Denoise ####
dadaR1 <- dada(derepR1, err = errR1, multithread = TRUE, verbose = TRUE)
dadaR2 <- dada(derepR2, err = errR2, multithread = TRUE, verbose = TRUE)

#### Merge paired reads ####
merged <- mergePairs(dadaR1, derepR1, dadaR2, derepR2,
                                minOverlap = 10, maxMismatch = 1, verbose=TRUE)

#### Merge no-mismatch pairs ####



#### Generate ASV table ####
asvtab <- makeSequenceTable(mergers)
saveRDS(file.path(asv_path, "asv_tab.rds"))
save.image(file.path(asv_path, "denoising.Rdata"))
