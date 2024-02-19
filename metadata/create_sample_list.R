#define paths
path <- "."
meta_path <- file.path(path, "metadata")
seq_path <- file.path(path, "sequences")
raw_path <- file.path(seq_path, "01_raw")
trim_path <- file.path(seq_path, "02_trim")
filt_path <- file.path(seq_path, "03_filter")

cyclone_sample_regex <- "^C[0-9A-Z]{5}$"
cyclone_neg_regex <- "^((GSSP|CCDB)-\\d{5})?(NEGEXT|NEGPCR[12]|control[123])$"
soil_sample_regex <- "^LIFEP-GSSP-([12])-(S[A-Z0-9]{5})$"
soil_neg_regex <- "^LIFEP-GSSP-([12])-((Neg|PCR)\\d{1,2})$"

sample_table <- tibble::tibble(
    fastq_R1 = sort(list.files(raw_path, ".*R1(_001)?.fastq.gz", recursive = TRUE)),
    fastq_R2 = sort(list.files(raw_path, ".*R2(_001)?.fastq.gz", recursive = TRUE))
  ) |>
  # parse filenames
  tidyr::extract(
    fastq_R1,
    into = c("seqrun", "sample", "lane"),
    regex = "([^/]+)/(?:.*/)?(.+?)_(?:S\\d+)_(L\\d+)_R1(?:_001)?.fastq.gz",
    remove = FALSE
  ) |>
  dplyr::mutate(
    seqrun = if (dplyr::n_distinct(lane) > 1L) paste(seqrun, lane, sep = "_") else seqrun,
    .by = seqrun,
    .keep = "unused"
  ) |>
  dplyr::filter(
    grepl(cyclone_sample_regex, sample) |
    grepl(cyclone_neg_regex, sample) |
    grepl(soil_sample_regex, sample) |
    grepl(soil_neg_regex, sample)
  ) |>
  dplyr::mutate(
    neg_control = grepl(cyclone_neg_regex, sample) | grepl(soil_neg_regex, sample),
    sample = dplyr::case_when(
      grepl(cyclone_sample_regex, sample) ~ sample,
      grepl(cyclone_neg_regex, sample) ~ paste(seqrun, sub(cyclone_neg_regex, "\\3", sample), sep = "_"),
      grepl(soil_sample_regex, sample) ~ sub(soil_sample_regex, "\\2_Rep\\1", sample),
      grepl(soil_neg_regex, sample) ~ paste(seqrun, sub(soil_neg_regex, "\\2_Rep\\1", sample), sep = "_")
    ),
    cut_R2 = ifelse(startsWith(seqrun, "BIONAME"), "0", "16")
  ) |>
  dplyr::select(seqrun, sample, neg_control, cut_R2, fastq_R1, fastq_R2)

readr::write_tsv(sample_table, file.path(meta_path, "sample_table.tsv"))


