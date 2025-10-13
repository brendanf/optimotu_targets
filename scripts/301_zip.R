output_targets <- rlang::syms(
  tarchetypes::tar_select_names(optimotu_plan, starts_with("write_"))
  )

optimotu_plan <- c(
  optimotu_plan,
  list(
    zip_output = tar_file(
      zip_output,
      {
        outfile <- file.path(
          !!optimotu.pipeline::output_path(),
          sprintf(
            "%s_%s.zip",
            !!optimotu.pipeline::project_name(),
            strftime(Sys.Date(), "%Y%m%d")
          )
        )
        i <- 1
        while (file.exists(outfile)) {
          outfile <- file.path(
            !!optimotu.pipeline::output_path(),
            sprintf("%s_%s_%i.zip",
              !!optimotu.pipeline::project_name(),
              strftime(Sys.Date(), "%Y%m%d"),
              i
            )
          )
          i <- i + 1
        }
        result <- zip(outfile, c(!!!output_targets), zip = "zip", flags = "-j9X")
        stopifnot(result == 0)
        outfile
      },
      deployment = "main"
    )
  )
)
