readr::read_tsv(list.files("metadata/UMIMaps", "LIFEPLAN_.*_UMIMap\\.txt", full.names = TRUE)[1:2], id = "umi_file") |>
  dplyr::transmute(
    sample = `Sample ID`,
    seqrun = substr(basename(umi_file), 1, 14) |> chartr(old = "_", new = "-"),
    fastq_R1 = sprintf("%s/%s.R1.fastq.gz", seqrun, sample),
    fastq_R2 = sprintf("%s/%s.R2.fastq.gz", seqrun, sample),
    orient = ifelse(startsWith(Forward, "BF3"), "fwd", "rev")
  ) |>
  readr::write_tsv("metadata/samples.tsv")
