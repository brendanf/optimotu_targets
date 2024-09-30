parse_project_name <- function(pipeline_options) {
  if (!("project_name" %in% names(pipeline_options))
      || length(pipeline_options$project_name) == 0) {
    warning(
      "Missing project name in 'pipeline_options.yaml'.\n",
      "Using default name 'metabarcoding_project'."
    )
    return("metabarcoding_project")
  } else if (length(pipeline_options$project_name) > 1) {
    stop("Project can only have one name (file: pipeline_options.yaml)")
  } else if (pipeline_options$project_name == "metabarcoding_project") {
    message(
      "Option 'project_name' is the default value 'metabarcoding_project'.\n",
      "You can change it by editing the file 'pipeline_options.yaml'"
    )
  } else if (!grepl("^[[:alnum:]_-]+$", pipeline_options$project_name)) {
    stop("Project name should consist of alphanumeric characters, '_', and '-'. (file:pipeline_options.yaml)")
  }
  pipeline_options$project_name
}

parse_primer_R1 <- function(pipeline_options) {
  checkmate::assert_string(
    pipeline_options$forward_primer,
    null.ok = TRUE,
    min.chars = 10,
    pattern = "[ACGTSWRYMKBDHVIN]+",
    ignore.case = TRUE
  )
  if (is.null(pipeline_options$forward_primer)) {
    message("Forward primer string missing (file: pipeline_options.yaml)\n",
            "Using default: GCATCGATGAAGAACGCAGC")
    "GCATCGATGAAGAACGCAGC"
  } else {
    pipeline_options$forward_primer
  }
}

parse_primer_R2 <- function(pipeline_options) {
  checkmate::assert_string(
    pipeline_options$reverse_primer,
    null.ok = TRUE,
    min.chars = 10,
    pattern = "[ACGTSWRYMKBDHVIN]+",
    ignore.case = TRUE
  )
  if (is.null(pipeline_options$reverse_primer)) {
    message("Reverse primer string missing (file: pipeline_options.yaml)\n",
            "Using default: TCCTCCGCTTATTGATATGC")
    "TCCTCCGCTTATTGATATGC"
  } else {
    pipeline_options$reverse_primer
  }
}

parse_amplicon_model_options <- function(amplicon_model_options) {
  checkmate::assert_list(amplicon_model_options)
  amplicon_model_options <- unnest_yaml_list(amplicon_model_options)
  checkmate::assert_names(
    names(amplicon_model_options),
    must.include = "model_type"
  )
  checkmate::assert_string(amplicon_model_options$model_type)
  checkmate::assert_subset(
    amplicon_model_options$model_type,
    c("CM", "HMM", "none")
  )
  amplicon_model_type <<- amplicon_model_options$model_type
  if (!identical(amplicon_model_type, "none")) {
    checkmate::assert_names(
      names(amplicon_model_options),
      must.include = "model_file"
    )
    if ("model_trim" %in% names(amplicon_model_options)) {
      checkmate::assert_flag(amplicon_model_options$model_trim)
      do_model_trim <<- amplicon_model_options$model_trim
    }
    if ("model_refine" %in% names(amplicon_model_options)) {
      checkmate::assert_flag(amplicon_model_options$model_refine)
      do_model_refine <<- amplicon_model_options$model_refine
    }

    checkmate::assert_string(amplicon_model_options$model_file)
    model_file <<- amplicon_model_options$model_file
    if (isFALSE(do_model_refine) && isFALSE(do_model_trim)) {
      checkmate::assert_file_exists(model_file)
    } else {
      checkmate::assert_path_for_output(model_file, overwrite = TRUE)
    }

    if ("model_seed" %in% names(amplicon_model_options)) {
      model_seed <<- amplicon_model_options$model_seed
      checkmate::assert_file_exists(model_seed, "r")
    }

    ##### amplicon model filtering settings #####
    if (!is.null(amplicon_model_options$model_filter)) {
      do_model_filter <- TRUE
      checkmate::assert_list(amplicon_model_options$model_filter, min.len = 1)
      checkmate::assert_names(
        names(amplicon_model_options$model_filter),
        subset.of = c("max_model_start", "min_model_end", "min_model_score")
      )
      model_filter <<- unnest_yaml_list(amplicon_model_options$model_filter)
      if ("max_model_start" %in% names(model_filter)) {
        checkmate::assert_number(model_filter$max_model_start)
      } else {
        model_filter$max_model_start <<- Inf
      }

      if ("min_model_end" %in% names(model_filter)) {
        checkmate::assert_number(model_filter$min_model_end)
      } else {
        model_filter$min_model_end <<- -Inf
      }

      if ("min_model_score" %in% names(model_filter)) {
        checkmate::assert_number(model_filter$min_model_score)
      } else {
        model_filter$min_model_score <<- -Inf
      }
    }

    #### amplicon alignment settings ###
    if (!is.null(amplicon_model_options$model_align)) {
      checkmate::assert_flag(amplicon_model_options$model_align)
      do_model_align <<- amplicon_model_options$model_align
    }

    #### NuMt detection settings ####
    if ("numt_filter" %in% names(amplicon_model_options)) {
      checkmate::assert_logical(amplicon_model_options$numt_filter)
      do_numt_filter <<- amplicon_model_options$numt_filter
      if (do_numt_filter && amplicon_model_type != "HMM")
        stop("NuMt filter is only valid when HMM alignment is used")
    }
  }
}

parse_taxonomy_ranks <- function(rank_options) {
  checkmate::assert(
    checkmate::check_list(
      rank_options,
      types = c("character", "list"),
      min.len = 1
    ),
    checkmate::check_character(
      rank_options,
      unique = TRUE,
      min.len = 1
    )
  )
  KNOWN_RANKS <<- purrr::keep(
    rank_options,
    \(x) dplyr::cumall(checkmate::test_list(x))
  ) |>
    unlist()
  UNKNOWN_RANKS <<- purrr::discard(
    rank_options,
    \(x) dplyr::cumall(checkmate::test_list(x))
  ) |>
    unlist()
  if (length(UNKNOWN_RANKS) == 0 || !is.null(names(UNKNOWN_RANKS))) {
    stop(
      "Option 'taxonomy':'ranks' should start from the most inclusive rank (e.g. kingdom)\n",
      "  and continue to the least inclusive rank (e.g. species).  Optionally the first\n",
      "  rank(s) may be defined (e.g. '- kingdom: Fungi') but subsequent ranks must be \n",
      "  undefined (e.g. '- class')."
    )
  }
  ROOT_TAXON <<- unname(KNOWN_RANKS[1])
  KNOWN_TAXA <<- unname(KNOWN_RANKS)
  KNOWN_RANKS <<- names(KNOWN_RANKS)
  TAX_RANKS <<- c(KNOWN_RANKS, UNKNOWN_RANKS)
  found_ranks <<- TRUE
}

parse_protax_options <- function(protax_options) {
  checkmate::assert_list(protax_options)
  protax_options <- unnest_yaml_list(protax_options)

  ##### protax version #####
  if ("aligned" %in% names(protax_options)) {
    checkmate::assert_flag(protax_options$aligned)
    protax_aligned <<- protax_options$aligned
    protax_unaligned <<- !protax_aligned
  } else {
    message("Using unaligned protax by default.")
  }

  ##### protax location #####
  if ("location" %in% names(protax_options)) {
    checkmate::assert_directory_exists(protax_options$location)
    protax_root <- protax_options$location
  } else {
    message("Using default protax directory: ", protax_root)
  }

  if ("ranks" %in% names(protax_options)) {
    parse_taxonomy_ranks(protax_options$ranks)
  }
}

parse_sintax_options <- function(sintax_options) {
  checkmate::assert_list(sintax_options)
  sintax_options <- unnest_yaml_list(sintax_options)
  do_sintax <<- TRUE
  do_protax_aligned <<- FALSE
  do_protax_unaligned <<- FALSE
  checkmate::assert_file_exists(sintax_options$reftax, "r")
  sintax_ref <<- sintax_options$reftax
}

parse_taxonomy_options <- function(taxonomy_options) {
  checkmate::assert_list(taxonomy_options)
  taxonomy_options <- unnest_yaml_list(taxonomy_options)
  if ("protax" %in% names(taxonomy_options)) {
    if ("sintax" %in% names(taxonomy_options)) {
      stop("Only one of options 'taxonomy:protax' and 'taxonomy:sintax'",
           " may be given. (File: pipeline_options.yaml)")
    }
    #### taxonomy protax ####
    parse_protax_options(taxonomy_options$protax)
  } else if ("sintax" %in% names(taxonomy_options)) {
    #### taxonomy sintax ####
    parse_sintax_options(taxonomy_options$sintax)
  }

  ##### taxonomy ranks #####
  if ("ranks" %in% names(taxonomy_options)) {
    parse_taxonomy_ranks(taxonomy_options$ranks)
  }
}
