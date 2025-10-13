#### find all the sequencing files and load the sequencing metadata

local({
  # the sample table may be several Mb, so avoid putting it in the global
  # namespace. It is cached in `optimotu.pipeline::sample_table()`, so it does
  # not need to be recomputed every time, which can be especially slow on
  # networked file systems.
  sample_table <- optimotu.pipeline::sample_table()
  cat("Found", dplyr::n_distinct(sample_table$sample, sample_table$seqrun),
      "samples in", optimotu.pipeline::n_seqrun(), "runs.\n",
      "sample_table targets hash is:", targets:::hash_object(sample_table), "\n"
  )
  for (n in colnames(sample_table)) {
    cat(sprintf("sample_table$%s hash: %s\n", n, targets:::hash_object(sample_table[[n]])))
  }
})
