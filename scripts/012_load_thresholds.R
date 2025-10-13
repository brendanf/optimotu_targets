# test the performance of different clustering thresholds for replicating
# "true" taxonomic groups.

fmeasure_optima <- readr::read_tsv(optimotu.pipeline::cluster_thresholds(),
                                   col_types = "ccccdd")
