# parse Unite SH assignments
# copied from GSSP 
library(magrittr)
library(targets)
library(tarchetypes)

unite_thresholds <- c(
  "97" = 97,
  "975" = 97.5,
  "98" = 98,
  "985" = 98.5,
  "99" = 99
) %>%
  tibble::enframe()

jobnumber <- 1
sh_infile <- sprintf("indata/source_%d", jobnumber)
sh_outfile <- sprintf("outdata/source_%d.zip", jobnumber)
sh_datafile <- "sh_matching_data_0_5.zip"
sh_dataurl <- "https://files.plutof.ut.ee/public/orig/E4/7C/E47CE4EEBC48A22618FBDF07218E4EC8DC32170CE75D55914765FBF3CA455CB2.zip"

SH_plan <- list(
  #### asvs_to_unite ####
  # write a fasta file of non-spike ASV sequences for SH matching
  tar_file(
    asvs_to_unite,
    write_sequence(asv_seq, sh_infile),
    deployment = "main",
    priority = 1
  ),
  #### sh_matching_script ####
  tar_file(
    sh_matching_script,
    "sh_matching_pub/sh_matching_analysis/run_pipeline.sh"
  ),
  #### sh_matching_analysis ####
  tar_file(
    sh_matching_analysis,
    list.files(path="sh_matching_pub/sh_matching_analysis/scripts", full.names = TRUE)
  ),
  #### sh_matching_data ####
  tar_file(
    sh_matching_data,
    {
      download.file(
        url = sh_matching_dataurl,
        destfile = sh_matching_datafile,
        quiet = TRUE
      )
      out <- unzip(
        zipfile = sh_matching_datafile,
        exdir = "data/sh_matching_data",
        junkpaths = TRUE,
        overwrite = TRUE
      )
      unlink(sh_matching_datafile)
      out
    }
  ),
  #### sh_matching ####
  # run the SH matching pipeline locally
  tar_file(
    sh_matching,
    {
      # mention these here just for dependency tracking
      asvs_to_unite
      sh_matching_analysis
      sh_matching_data
      # remove the output file
      unlink(sh_outfile, force = TRUE)
      # run the pipeline
      system2(sh_matching_script, c(jobnumber, "its2"))
      # return the output file (for dependency tracking)
      sh_outfile
    },
    priority = 1
  ),
  
  #### unite_excluded ####
  # look up the sequences which Unite excluded, and why.
  # typically these will be based on chimeras or length
  tar_fst_tbl(
    unite_excluded,
    dplyr::left_join(
      readr::read_tsv(
        unz(sh_matching, !!sprintf("excluded_%d.txt", jobnumber)),
        col_types = "ccc",
        col_names = c("seq_id_tmp", "step", "reason")
      ),
      readr::read_tsv(
        unz(sh_matching, !!sprintf("source_%s_names", jobnumber)),
        col_types = "cci",
        col_names = c("seq_accno", "seq_id_tmp", "n")
      ),
      by = "seq_id_tmp"
    )
  ),
  tar_map(
    values = unite_thresholds,
    names = name,
    #### unite_matches_out ####
    # parse the Unite matches_out* files
    # matches_out is for sequences with at least 75% similarity to a sequence
    # already in Unite.
    # matches_out_1 is for sequences which do NOT have a match at 75% similarity
    # to a sequence already in Unite.
    tar_fst_tbl(
      unite_matches_out,
      purrr::map_dfr(
        c("matches", "matches_1"),
        # for some reason readr::read_tsv had parsing errors
        ~ readLines(
          unz(sh_matching, sprintf("matches/%s_out_%s.csv", .x, name))
        ) %>%
          tibble::tibble(data = .) %>%
          tidyr::separate(
            col = data,
            into = unlist(strsplit(.$data[1], "\t")),
            sep = "\t",
            extra = "drop"
          ) %>%
          dplyr::slice(-1)
      ) %>%
        # the Unite column names include the threshold (as dissimilarity),
        # but we are taking care of that by keeping the data separate.
        dplyr::rename_with(.fn = sub, pattern = " \\(\\d\\.\\d\\)",
                           replacement = "") %>%
        tidyr::separate(
          `SH/compound taxonomy`,
          into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
          sep = ";",
          fill = "right",
          remove = FALSE
        ) %>%
        dplyr::mutate(
          dplyr::across(kingdom:species, sub, pattern = "^[kpcofgs]__", replacement = ""),
          dplyr::across(kingdom:species, dplyr::na_if, "unidentified"),
          species = ifelse(endsWith(species, "_sp"), NA_character_, species)
        )
    )
  )
)
