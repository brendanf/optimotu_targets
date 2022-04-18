# test the performance of different clustering thresholds for replicating
# "true" taxonomic groups.

threshold_file <- file.path(meta_path, "GSSP_thresholds.tsv")
fmeasure_optima <- readr::read_tsv(threshold_file, col_types = "ccccdd")
