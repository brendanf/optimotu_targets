output_targets <- rlang::syms(
  tarchetypes::tar_select_names(optimotu_plan, starts_with("write_"))
  )

optimotu_plan <- c(
  optimotu_plan,
  list(
    #### git_diff ####
    # `character` file name
    # file giving changes that have been made to the pipeline code since the
    # last commit, or a
    git_diff = tar_file(
      git_diff,
      {
        outfile <-
          optimotu.pipeline::ensure_directory("output/pipeline_changes.diff")
        result <- processx::run(
          "git",
          args = "diff",
          error_on_status = FALSE,
          stdout = outfile,
          wd = !!dirname(targets::tar_config_get("script"))
        )
        if (result$status != 0L) {
          writeLines(
            c(
              "'git diff' failed, no tracking information available",
              "If you want this information, make sure git is installed and ",
              "that you have checked out the project from git (rather than just",
              "downloading it as a .zip, for instance."
            ),
            outfile
          )
        }
        outfile
      },
      deployment = "main",
      cue = tar_cue("always")
    ),
    #### git_commit ####
    # `character` file name
    # file giving the full 40-digit hash of the current commit
    git_commit = tar_file(
      git_commit,
      {
        outfile <-
          optimotu.pipeline::ensure_directory("output/pipeline_version.txt")
        result <- processx::run(
          "git",
          args = c("rev-parse", "HEAD"),
          error_on_status = FALSE,
          stdout = outfile,
          wd = !!dirname(targets::tar_config_get("script"))
        )
        if (result$status != 0L) {
          writeLines(
            c(
              "'git rev-parse' failed, no commit information available",
              "If you want this information, make sure git is installed and ",
              "that you have checked out the project from git (rather than just",
              "downloading it as a .zip, for instance."
            ),
            outfile
          )
        }
        outfile
      },
      deployment = "main",
      cue = tar_cue("always")
    ),

    session_info = tar_file(
      session_info,
      {
        outfile <-
          file.path(!!optimotu.pipeline::output_path(), "sessionInfo.txt")
        # load all the package specified in the conda environment file
        file.path(
          !!dirname(targets::tar_config_get("script")),
          "conda/OptimOTU_v6.yaml"
        ) |>
          readLines() |>
          grep("- (bioconducto)?r-", x = _, value = TRUE) |>
          sub("^ *- (bioconducto)?r-", "", x = _) |>
          sub("=.+$", "", x = _) |>
          purrr::walk(
            \(x) if (requireNamespace(x, quietly = TRUE)) loadNamespace(x)
          )
        sink(outfile)
        print(sessionInfo())
        sink()
        outfile
      },
      deployment = "main",
      cue = tar_cue("always")
    ),

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
        result <- zip(outfile, c(git_commit, git_diff, session_info,
          "pipeline_options.yaml",
          !!optimotu.pipeline::custom_sample_table(), !!!output_targets),
          zip = "zip", flags = "-j9X")
        stopifnot(result == 0)
        outfile
      },
      deployment = "main"
    )
  )
)
