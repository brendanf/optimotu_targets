# DADA2 quality filtering and denoising for Sonja's spruce log metabarcoding data
# Brendan Furneaux
# Based on DADA2 analysis for GSSP from Jenni Hultman

# ensure we have all necessary packages
# if already installed (e.g. on CSC) this will copy them once per user to a
# cache in ~/.cache/R
renv::hydrate()

library(dada2)
library(tidyverse)
library(magrittr)

#how many cpu cores do we have?
ncpus <- Sys.getenv("SLURM_CPUS_PER_TASK")
if (nchar(ncpus) == 0) {
  ncpus <- parallel::detectCores()
} else {
  ncpus <- as.integer(ncpus)
}

cat("Running DADA2 with", ncpus, "cores.\n" )

#define paths
path <- "."
seq_path <- file.path(path, "sequences")
raw_path <- file.path(seq_path, "01_raw")
filt_path <- file.path(seq_path, "02_filter")

# create paths if missing
if (!dir.exists(filt_path)) dir.create(filt_path, recursive = TRUE)

# find files
sample_table <- tibble(
  fastq_R1 = list.files(raw_path, ".*R1(_001)?.fastq.gz"),
  fastq_R2 = list.files(raw_path, ".*R2(_001)?.fastq.gz")
) %>%
# parse filenames
  tidyr::extract(
    fastq_R1,
    into = "sample",
    regex = paste0(
      "(?:MI_M06648_\\d+\\.001\\.N\\d+--S\\d+\\.)?",
      "(19[KFLS][EAU]\\d{3}|CCDB-\\d+NEG(?:EXT|PCR[12]))",
      "_(?:S\\d+_L001_)?R1(?:_001)?.fastq.gz"
    ),
    remove = FALSE
  ) %>%
# generate filenames for filtered reads
  dplyr::mutate(
    filt_R1 = file.path(filt_path, paste0(sample, "_R1_filt.fastq.gz")),
    filt_R2 = file.path(filt_path, paste0(sample, "_R2_filt.fastq.gz"))
  ) %>%
  dplyr::filter(!is.na(sample))

cat("Found", nrow(sample_table), "samples.\n")

out2 <- sample_table %$%
  filterAndTrim(
    fwd = file.path(raw_path, fastq_R1),
    filt = file.path(filt_path, filt_R1),
    rev = file.path(raw_path, fastq_R2),
    filt.rev = file.path(filt_path, filt_R2),
    maxN = 0, # max 0 ambiguous bases
    maxEE = c(3, 5), # max expected errors (fwd, rev)
    truncQ = 2, # truncate at first base with quality <= 2
    rm.phix = TRUE, #remove matches to phix genome
    minLen = 50, # remove reads < 50bp
    compress = TRUE, # write compressed files
    multithread = ncpus,
    verbose = TRUE
  )
