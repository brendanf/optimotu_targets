# Bioinformatics pipeline for GSSP
#   by Brendan Furneaux, 2021-2023

library(targets)
library(tarchetypes)
library(magrittr)
library(qs)
library(fst)

tar_option_set(
  format = "qs",
  memory = "transient",
  garbage_collection = TRUE,
  priority = 0.5,
  workspace_on_error = TRUE
)

# Numbered R scripts define the targets plan.
# They are numbered in the order they are used.
for (f in list.files("scripts", "^\\d{3}_.+.R$", full.names = TRUE)) {
  source(f)
}

cat("Detected", local_cpus(), "cores for main process.\n" )

# End this file with a list of target objects.
invisible(
  list(
    dada_plan,
    asv_plan,
    protax_plan,
    threshold_plan,
    clust_plan
  )
)
