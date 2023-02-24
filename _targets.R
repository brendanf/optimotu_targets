# Bioinformatics pipeline for Sonja Saine's spruce log metabarcoding project
#   by Brendan Furneaux
# Based on:
#   DADA2 analysis for GSSP from Jenni Hultman
#   PROTAX taxonomy assignment for GSSP by Panu Somervuo
#   Clustering for GSSP by Brendan Furneaux

library(targets)
library(tarchetypes)
library(magrittr)
library(qs)
library(fst)

# Numbered R scripts define the targets plan.
# They are numbered in the order they are used.
for (f in list.files("scripts", "^\\d{3}_.+.R$", full.names = TRUE)) {
  source(f)
}

cat("Detected", local_cpus(), "cores for main process.\n" )

tar_option_set(
  format = "qs",
  memory = "transient",
  garbage_collection = TRUE,
  priority = 0.5
)

# End this file with a list of target objects.
list(
  dada_plan,
  asv_plan,
  refseq_plan,
  protax_plan,
#  SH_plan,
  clust_plan,
  target_taxa_plan
)
