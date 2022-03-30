library(tarchetypes)

protax_plan <- list(
  tar_fst_tbl(
    asv_seq,
    tibble::tibble(
      ASV = sprintf(
        sprintf("ASV%%0%dd", ceiling(log10(ncol(asvtable)))),
        seq_len(ncol(asvtable))
      ),
      seq = colnames(asvtable)
    ),
    deployment = "main"
  ),
  
  tar_group_count(
    grouped_asv_seq,
    asv_seq,
    count = 12,
    deployment = "main"
  ),
  
  tar_file(
    protax,
    run_protax(grouped_asv_seq, file.path(protax_path, tar_name())),
    pattern = map(grouped_asv_seq)
  )
)
