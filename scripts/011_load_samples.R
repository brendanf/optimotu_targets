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
protax_path <- file.path(seq_path, "05_protax")

# create paths if missing
if (!dir.exists(filt_path)) dir.create(filt_path, recursive = TRUE)
if (!dir.exists(trim_path)) dir.create(trim_path, recursive = TRUE)
if (!dir.exists(asv_path)) dir.create(asv_path, recursive = TRUE)
if (!dir.exists(protax_path)) dir.create(protax_path, recursive = TRUE)

if (!is.null(pipeline_options$custom_sample_table)) {
  #todo: support sample table in other formats (csv, excel, ...)
  sample_table <- readr::read_tsv(
    pipeline_options$custom_sample_table,
    col_types = readr::cols(
      seqrun = readr::col_character(),
      sample = readr::col_character(),
      fastq_R1 = readr::col_character(),
      fastq_R2 = readr::col_character(),
      .default = readr::col_guess()
    )
  )
  checkmate::check_data_frame(
    sample_table,
    col.names = "named"
  )
  checkmate::check_names(
    names(sample_table),
    must.include = c("sample", "seqrun", "fastq_R1", "fastq_R2")
  )
  checkmate::check_character(sample_table$sample, any.missing = FALSE, unique = TRUE)
  checkmate::check_character(sample_table$seqrun, any.missing = FALSE)
  checkmate::check_file_exists(file.path(raw_path, sample_table$fastq_R1), access = "r")
  checkmate::check_file_exists(file.path(raw_path, sample_table$fastq_R2), access = "r")
} else {
  # find files
  sample_table <- tibble::tibble(
    fastq_R1 = sort(list.files(raw_path, ".*R1(_001)?.fastq.gz", recursive = TRUE)),
    fastq_R2 = sort(list.files(raw_path, ".*R2(_001)?.fastq.gz", recursive = TRUE))
  ) %>%
    # parse filenames
    tidyr::extract(
      fastq_R1,
      into = c("seqrun", "sample"),
      regex = "([^/]+)/(?:.*/)?(.+?)_(?:S\\d+_L001_)?R1(?:_001)?.fastq.gz",
      remove = FALSE
  ) %>%
  dplyr::mutate(
    sample = ifelse(
      startsWith(sample, "BLANK"),
      paste(seqrun, sample, sep = "_"),
      sample
    )
  )
}
sample_table <- sample_table %>%
  # generate filenames for trimmed and filtered reads
  dplyr::mutate(
    trim_R1 = file.path(trim_path,
                        paste(seqrun, sample, "R1_trim.fastq.gz", sep = "_")),
    trim_R2 = file.path(trim_path,
                        paste(seqrun, sample, "R2_trim.fastq.gz", sep = "_")),
    filt_key = file.path(filt_path, paste(seqrun, sample, sep = "_")),
    filt_R1 = paste(filt_key, "R1_filt.fastq.gz", sep = "_"),
    filt_R2 = paste(filt_key, "R2_filt.fastq.gz", sep = "_")
  )

# spike_strength is used along with the nonspike/spike ratio to convert from
# read number to "weight"
if (!("spike_weight") %in% names(sample_table))
  sample_table$spike_weight <- 1

assertthat::assert_that(
  !any(is.na(sample_table$seqrun)),
  !any(is.na(sample_table$sample)),
  is.numeric(sample_table$spike_weight),
  !any(duplicated(sample_table$fastq_R1)),
  !any(duplicated(sample_table$fastq_R2)),
  !any(duplicated(sample_table$trim_R1)),
  !any(duplicated(sample_table$trim_R2)),
  !any(duplicated(sample_table$filt_R1)),
  !any(duplicated(sample_table$filt_R2))
)

n_seqrun <- dplyr::n_distinct(sample_table$seqrun)

cat("Found", nrow(sample_table), "samples in", n_seqrun, "runs.\n",
    "sample_table hash is:", digest::digest(sample_table), "\n",
    "sample_table targets hash is:", targets:::digest_obj64(sample_table), "\n"
)
