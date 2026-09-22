readr::read_tsv(list.files("metadata/UMIMaps", "LIFEPLAN_.*_UMIMap\\.txt", full.names = TRUE), id = "umi_file") |>
  dplyr::transmute(
    sample = `Sample ID` |>
      sub("BLANK_(Plate[01][0-9])_(H0[789])", "AMPtk_CTRL_\\2_\\1", x = _),
    seqrun = substr(basename(umi_file), 1, 14) |> chartr(old = "_", new = "-"),
    sample = ifelse(
      seqrun == "LIFEPLAN-00001", 
      sub("BLANK_(Plate[01][0-9])_(H1[012])", "BLANKDNA_CTRL_\\2_\\1", sample),
      sample
    ),
    fastq_R1 = sprintf("%s/%s.R1.fastq.gz", seqrun, sample),
    fastq_R2 = sprintf("%s/%s.R2.fastq.gz", seqrun, sample),
    orient = ifelse(startsWith(Forward, "BF3"), "fwd", "rev")
  ) |>
  readr::write_tsv("metadata/samples.tsv")
