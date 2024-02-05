output_targets <- rlang::syms(
  tarchetypes::tar_select_names(optimotu_plan, starts_with("write_"))
  )

optimotu_plan <- c(
  optimotu_plan,
  list(
    tar_file_fast(
      zip_output,
      {
        outfile <- sprintf(
          "%s/%s_%s.zip",
          "output",
          project_name,
          strftime(Sys.Date(), "%Y%m%d")
        )
        result <- zip(outfile, c(!!!output_targets), zip = "zip")
        stopifnot(result == 0)
        outfile
      }
    )
  )
)
