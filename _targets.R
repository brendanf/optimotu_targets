# OptimOTU pipeline
#   by Brendan Furneaux
# Based on:
#   DADA2 analysis for GSSP from Jenni Hultman
#   PROTAX taxonomy assignment for GSSP by Panu Somervuo
#   Clustering for GSSP by Brendan Furneaux

library(targets)
library(tarchetypes)
library(qs2)
library(fst)

tar_option_set(
  format = "qs",
  memory = "transient",
  garbage_collection = TRUE,
  priority = 0.5,
  workspace_on_error = TRUE
)

min_pipeline_version <- "0.5.2.9011"

if (packageVersion("optimotu.pipeline") < min_pipeline_version) {
  stop("optimotu.pipeline version ", packageVersion("optimotu.pipeline"),
  " is too old.  Please update to version ", min_pipeline_version, " or later.")
}

min_optimotu_version <- "0.9.3.9003"

if (packageVersion("optimotu") < min_pipeline_version) {
  stop("optimotu version ", packageVersion("optimotu"),
       " is too old.  Please update to version ", min_optimotu_version, " or later.")
}

optimotu_plan <- list()

# if the pipeline code is inside a container, then the script and bin directories
# are not in the working directory
script_dir <- file.path(dirname(targets::tar_config_get("script")), "scripts")
bin_dir <- file.path(dirname(targets::tar_config_get("script")), "bin")
Sys.setenv(OPTIMOTU_BIN_DIR = bin_dir)

# Numbered R scripts define the targets plan.
# They are numbered in the order they are used.

# certain runners (all except run_node) run the 0** scripts, because they need
# information from the configuration file.  Only run them now if they have not
# been run before.
if (!isTRUE(optimotu.pipeline::did_pipeline_options())) {
  for (f in list.files(script_dir, "^0[[:digit:]]{2}_.+[.]R$", full.names = TRUE)) {
    source(f)
  }
}

# the 1** and higher scripts are the ones which actually define the plan.
for (f in list.files(script_dir, "^[1-9][[:digit:]]{2}_.+[.]R$", full.names = TRUE)) {
  source(f)
}

cat("Detected", optimotu.pipeline::local_cpus(), "cores for main process.\n" )

# Make sure log directory exists
if (!dir.exists("logs")) dir.create("logs")

# End this file with a list of target objects.
invisible(optimotu_plan)
