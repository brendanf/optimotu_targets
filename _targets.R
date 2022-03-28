# Bioinformatics pipeline for Sonja Saine's spruce log metabarcoding project
#   by Brendan Furneaux
# Based on:
#   DADA2 analysis for GSSP from Jenni Hultman
#   PROTAX taxonomy assignment for GSSP by Panu Somervuo
#   Clustering for GSSP by Brendan Furneaux

# ensure we have all necessary packages
# if already installed (e.g. on CSC) this will copy them once per user to a
# cache in ~/.cache/R
renv::hydrate()

library(targets)
library(tarchetypes)
library(magrittr)
library(qs)
library(fst)

# Numbered R scripts define the targets plan.
# They are numbered in the order they are used.
for (f in list.files("scripts", "\\d{3}_.+.R", full.names = TRUE)) {
  source(f)
}

cat("Running pipeline with", local_cpus(), "cores.\n" )

tar_option_set(format = "qs", memory = "transient", garbage_collection = TRUE)

# End this file with a list of target objects.
list(
  dada_plan
)
