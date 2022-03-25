library(magrittr)

#### find all the sequencing files and load the sequencing metadata

#### Get the samples organized ####

#define paths
path <- "."
meta_path <- file.path(path, "metadata")
seq_path <- file.path(path, "sequences")
raw_path <- file.path(seq_path, "01_raw")
trim_path <- file.path(seq_path, "02_trim")
filt_path <- file.path(seq_path, "03_filter")
asv_path <- file.path(seq_path, "04_denoised")

# create paths if missing
if (!dir.exists(filt_path)) dir.create(filt_path, recursive = TRUE)
if (!dir.exists(trim_path)) dir.create(trim_path, recursive = TRUE)
if (!dir.exists(asv_path)) dir.create(asv_path, recursive = TRUE)

# Get sequencing metadata
metadata_files <- list.files(path = meta_path, pattern = "*.xlsm", full.names = TRUE)
metadata <- purrr::map_dfr(
  metadata_files,
  readxl::read_xlsx,
  skip = 2,
  col_types = c("text", "text", rep("skip", 14))
) %>%
  tidyr::separate(col = "Sample Locator", into = c("seqrun", "well"), sep = " ") %>%
  dplyr::filter(!is.na(`BOLD Sample IDs`))

# find files
sample_table <- tibble::tibble(
  fastq_R1 = list.files(raw_path, ".*R1(_001)?.fastq.gz"),
  fastq_R2 = list.files(raw_path, ".*R2(_001)?.fastq.gz")
) %>%
  # parse filenames
  tidyr::extract(
    fastq_R1,
    into = "sample",
    regex = paste0(
      "(?:MI_M06648_\\d+\\.001\\.N\\d+--S\\d+\\.)?",
      "(19[KFLS][EAU]\\d{3}|CCDB-\\d+NEG(?:EXT|PCR[12])|BLANK\\d)",
      "_(?:S\\d+_L001_)?R1(?:_001)?.fastq.gz"
    ),
    remove = FALSE
  ) %>%
  # generate filenames for trimmed and filtered reads
  dplyr::mutate(
    trim_R1 = file.path(trim_path, paste0(sample, "_R1_trim.fastq.gz")),
    trim_R2 = file.path(trim_path, paste0(sample, "_R2_trim.fastq.gz")),
    filt_R1 = file.path(filt_path, paste0(sample, "_R1_filt.fastq.gz")),
    filt_R2 = file.path(filt_path, paste0(sample, "_R2_filt.fastq.gz"))
  ) %>%
  # dplyr::filter(!is.na(sample)) %>%
  dplyr::left_join(metadata, by = c("sample" = "BOLD Sample IDs")) %>%
  dplyr::mutate(seqrun = dplyr::coalesce(seqrun, substring(sample, 1, 10)))

cat("Found", nrow(sample_table), "samples.\n")