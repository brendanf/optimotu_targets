############################################################
## TODO: ##
# - other sanity checks for the loaded settings
############################################################

# don't import the whole package, but let's use the null default operator
`%||%` <- rlang::`%||%`

#### load options ####
if (file.exists("pipeline_options.yaml")) {
  pipeline_options <- yaml::read_yaml("pipeline_options.yaml")
} else {
  warning(
    "Options file 'pipeline_options.yaml' is missing!\n",
    "Using defaults for all parameters."
  )
  pipeline_options <- list()
}

#### project_name ####

project_name <- parse_project_name(pipeline_options)

#### custom_sample_table ####
checkmate::assert(
  checkmate::check_null(pipeline_options$custom_sample_table),
  checkmate::check_false(pipeline_options$custom_sample_table),
  checkmate::check_file_exists(pipeline_options$custom_sample_table)
)
pipeline_options$custom_sample_table <-
  pipeline_options$custom_sample_table %||% FALSE

#### parallelism ####
checkmate::assert_count(pipeline_options$max_batchsize, na.ok = TRUE, null.ok = TRUE)
max_batchsize <- NULL
if (checkmate::test_count(pipeline_options$max_batchsize, positive = TRUE))
  max_batchsize <- pipeline_options$max_batchsize

workers_per_seqrun <- 1L
checkmate::assert_count(pipeline_options$workers_per_seqrun, positive = TRUE, null.ok = TRUE)
checkmate::assert_count(pipeline_options$jobs_per_seqrun, positive = TRUE, null.ok = TRUE)
if (!is.null(pipeline_options$workers_per_seqrun)) {
  if (!is.null(pipeline_options$jobs_per_seqrun)) {
    warning(
      "both 'workers_per_seqrun' and 'jobs_per_seqrun' (deprecated) were given",
      " in 'pipeline_options.yaml'. Using 'workers_per_seqrun' (=",
      pipeline_options$workers_per_seqrun, ")")
  }
  workers_per_seqrun <- pipeline_options$workers_per_seqrun
} else if (!is.null(pipeline_options$jobs_per_seqrun)) {
  message(
    "Option 'jobs_per_seqrun' is deprecated in 'pipeline_options.yaml'.",
    "Please use 'workers_per_seqrun' instead."
  )
  workers_per_seqrun <- pipeline_options$jobs_per_seqrun
}

min_workers <- 1L
checkmate::assert_int(pipeline_options$min_workers, lower = 1L, null.ok = TRUE)
if (!is.null(pipeline_options$min_workers))
  min_workers <- pipeline_options$min_workers

max_workers <- Inf
checkmate::assert_int(pipeline_options$max_workers, lower = min_workers, null.ok = TRUE)
if (!is.null(pipeline_options$max_workers))
  max_workers <- pipeline_options$max_workers

#### file_extension ####
checkmate::assert(
  checkmate::check_null(pipeline_options$file_extension),
  checkmate::check_string(pipeline_options$file_extension)
)

if (is.null(pipeline_options$file_extension)) {
  pipeline_options$file_extension <- "f(ast)?q([.]gz)?"
} else if (!is.null(pipeline_options$custom_sample_table) &&
           !isFALSE(pipeline_options$custom_sample_table)) {
  warning("Both 'custom_sample_table' and 'file_extension' options given ",
          "in 'pipeline_options.yaml'.\nIgnoring 'file_extension'.")
}


#### added_reference ####
if (!is.null(pipeline_options$added_reference)) {
  checkmate::assert_list(pipeline_options$added_reference)
  pipeline_options$added_reference <-
    unnest_yaml_list(pipeline_options$added_reference)
  checkmate::assert(
    checkmate::check_null(pipeline_options$added_reference$fasta),
    checkmate::check_file_exists(pipeline_options$added_reference$fasta)
  )

  checkmate::assert(
    checkmate::check_null(pipeline_options$added_reference$table),
    checkmate::check_file_exists(pipeline_options$added_reference$table)
  )

  if (xor(is.null(pipeline_options$added_reference$fasta),
          is.null(pipeline_options$added_reference$table))) {
    stop(
      "If one of 'added_reference_fasta' and 'added_reference_table' is given ",
      "in 'pipeline_options.yaml', then both must be given."
    )
  }
}

#### primers ####
primer_R1 <- parse_primer_R1(pipeline_options)
primer_R2 <- parse_primer_R2(pipeline_options)

# these are the primer sequences to send to cutadapt
trim_primer_R1 <- paste0(primer_R1, "...", dada2::rc(primer_R2), ";optional")
trim_primer_R2 <- paste0(primer_R2, "...", dada2::rc(primer_R1), ";optional")
trim_primer_merged <- paste0(primer_R1, "...", dada2::rc(primer_R2))

#### primer trim settings ####
checkmate::assert_list(pipeline_options$trimming, null.ok = TRUE)
if (is.null(pipeline_options$trimming)) {
  message("No 'trimming' options given in 'pipeline_options.yaml'\n",
          "Using defaults.")
  trim_options <- cutadapt_paired_options()
} else {
  trim_options <- do.call(
    cutadapt_paired_options,
    unnest_yaml_list(pipeline_options$trimming)
  )
}

#### filtering settings ####
dada2_maxEE <- c(2, 2)
checkmate::assert_list(pipeline_options$filtering, null.ok = TRUE)
if (is.null(pipeline_options$filtering)) {
  message("No 'filtering' options given in 'pipeline_options.yaml'\n",
          "Using defaults.")
} else {
  pipeline_options$filtering <- unnest_yaml_list(pipeline_options$filtering)
  checkmate::assert_names(
    names(pipeline_options$filtering),
    subset.of = c("maxEE_R1", "maxEE_R2")
  )
  checkmate::assert_number(
    pipeline_options$filtering$maxEE_R1,
    lower = 0,
    finite = TRUE,
    null.ok = TRUE
  )
  if (!is.null(pipeline_options$filtering$maxEE_R1))
    dada2_maxEE[1] <- pipeline_options$filtering$maxEE_R1
  checkmate::assert_number(
    pipeline_options$filtering$maxEE_R2,
    lower = 0,
    finite = TRUE,
    null.ok = TRUE
  )
  if (!is.null(pipeline_options$filtering$maxEE_R2))
    dada2_maxEE[2] <- pipeline_options$filtering$maxEE_R2
}

#### tag_jump settings ####
checkmate::assert(
  checkmate::check_list(pipeline_options$tag_jump, null.ok = TRUE),
  checkmate::check_false(pipeline_options$tag_jump)
)
if (is.null(pipeline_options$tag_jump) || isFALSE(pipeline_options$tag_jump)) {
  do_uncross <- FALSE
} else {
  do_uncross <- TRUE
  tagjump_options <- list(
    f = 0.05,
    p = 1.0
  )
  pipeline_options$tag_jump <- unnest_yaml_list(pipeline_options$tag_jump)
  checkmate::assert_names(
    names(pipeline_options$tag_jump),
    subset.of = c("f", "p")
  )
  checkmate::assert_number(
    pipeline_options$tag_jump$f,
    lower = 0,
    upper = 1,
    finite = TRUE,
    null.ok = TRUE
  )
  if (!is.null(pipeline_options$tag_jump$f))
    tagjump_options$f <- pipeline_options$tag_jump$f
  checkmate::assert_number(
    pipeline_options$tag_jump$p,
    lower = 0,
    finite = TRUE,
    null.ok = TRUE
  )
  if (!is.null(pipeline_options$tag_jump$p))
    tagjump_options$p <- pipeline_options$tag_jump$p
}

#### amplicon model settings ####
amplicon_model_type <- "none"
model_file <- character(0)
do_model_trim <- FALSE
do_model_refine <- FALSE
do_model_filter <- FALSE
model_filter <- list()
do_model_align <- FALSE
do_numt_filter <- FALSE
model_seed <- character(0)

if (!is.null(pipeline_options$amplicon_model)) {
  parse_amplicon_model_options(pipeline_options$amplicon_model)
}

do_model_align_only <- do_model_align && !do_model_filter
do_model_filter_only <- do_model_filter && !do_model_align
do_model_both <- do_model_align && do_model_filter

#### taxonomy settings ####
protax_aligned <- FALSE
protax_unaligned <- FALSE
protax_root <- "protaxFungi"

do_sintax <- FALSE
sintax_ref <- character()

KNOWN_RANKS <- TAX_RANKS[1]
UNKNOWN_RANKS <- TAX_RANKS[-1]
ROOT_TAXON <- "Fungi"
KNOWN_TAXA <- ROOT_TAXON
found_ranks <- FALSE

if (!is.null(pipeline_options$taxonomy)) {
  parse_taxonomy_options(pipeline_options$taxonomy)
} else if ("protax" %in% names(pipeline_options)) {
  parse_protax_options(pipeline_options$protax)
}
if (isFALSE(found_ranks)) {
  message("Using default ranks: ", paste(TAX_RANKS, collapse = ", "))
}

# these values (and TAX_RANKS) are treated as global variables in the sense that
# they are freely used inside functions where they are not passed as arguments.
ROOT_RANK <- TAX_RANKS[1]
ROOT_RANK_VAR <- rlang::sym(ROOT_RANK)
SECOND_RANK <- TAX_RANKS[2]
SECOND_RANK_VAR <- rlang::sym(SECOND_RANK)
INGROUP_RANK <- TAX_RANKS[length(KNOWN_RANKS)]
INGROUP_RANK_VAR <- rlang::sym(INGROUP_RANK)
INGROUP_TAXON <- KNOWN_TAXA[length(KNOWN_TAXA)]
RANK_OFFSET <- length(KNOWN_RANKS)
TIP_RANK <- TAX_RANKS[length(TAX_RANKS)]
TIP_RANK_VAR <- rlang::sym(TIP_RANK)

#### outgroup reference settings ####
outgroup_reference_file <- "data/sh_matching_data/sanger_refs_sh.fasta"
outgroup_taxonomy_file <- "data/sh_matching_data/shs_out.txt"
if (!is.null(pipeline_options$outgroup_reference)) {
  if (!is.null(pipeline_options$outgroup_reference$sequences)) {
    checkmate::assert_file_exists(pipeline_options$outgroup_reference$sequences)
    outgroup_reference_file <- pipeline_options$outgroup_reference$sequences
    outgroup_taxonomy_file <- NULL
  }
  if (!is.null(pipeline_options$outgroup_reference$taxonomy)) {
    checkmate::assert_file_exists(pipeline_options$outgroup_reference$taxonomy)
    outgroup_taxonomy_file <- pipeline_options$outgroup_reference$taxonomy
  }
}

#### cluster ####
optimize_thresholds <- FALSE
threshold_file <- "metadata/thresholds.tsv"
dist_config <- optimotu::dist_usearch(find_executable("usearch"))

if (!is.null(pipeline_options$cluster)) {

  ##### cluster dist_config #####
  checkmate::assert_list(pipeline_options$cluster$dist_config, null.ok = TRUE)
  if (is.null(pipeline_options$cluster$dist_config)) {
    dist_config <- optimotu::dist_usearch()
  } else {
    dist_config <- do.call(
      optimotu::dist_config,
      unnest_yaml_list(pipeline_options$cluster$dist_config)
    )
  }

  #### cluster optimize ####
  if (!is.null(pipeline_options$cluster$optimize)) {
    checkmate::assert(
      checkmate::check_flag(pipeline_options$cluster$optimize),
      checkmate::check_list(pipeline_options$cluster$optimize)
    )
    if (isTRUE(pipeline_options$cluster$optimize)) {
      stop(
        "Invalid value for `cluster:optimize` (file: 'pipeline_options.yaml')\n",
        "Option must be missing, FALSE, or a list of additional options.\n",
        "(Minimally cluster:optimize:test_data)"
      )
    } else if (is.list(pipeline_options$cluster$optimize)) {
      optimize_thresholds <- TRUE
      ##### cluster optimize train_data #####
      checkmate::assert_string(pipeline_options$cluster$optimize$train_data)
      checkmate::assert(
        checkmate::check_choice(
          pipeline_options$cluster$optimize$train_data,
          choices = c("reference", "self")
        ),
        checkmate::check_file_exists(
          pipeline_options$cluster$optimize$train_data,
          access = "r"
        )
      )
      optimize_train_data <- pipeline_options$cluster$optimize$train_data

      ##### cluster optimize min_conf #####
      checkmate::assert_number(
        pipeline_options$cluster$optimize$min_conf,
        lower = 0.5,
        upper = 1.0,
        null.ok = TRUE
      )
      if (!is.null(pipeline_options$cluster$optimize$min_conf) &&
          optimize_train_data != "self") {
        warning(
          "Option cluster:optimize:min_conf is ignored when option\n",
          "cluster:optimize:train_data is not 'self' (file: 'pipeline_options.yaml)"
        )
      }
      optimize_min_conf <- pipeline_options$cluster$optimize$min_conf %||% 0.5

      ##### cluster optimize dist_max #####
      checkmate::assert_number(
        pipeline_options$cluster$optimize$dist_max,
        lower = 0.0,
        upper = 1.0,
        null.ok = TRUE
      )
      optimize_max_dist <- pipeline_options$cluster$optimize$max_dist %||% 0.3

      ##### cluster optimize dist_step #####
      checkmate::assert_number(
        pipeline_options$cluster$optimize$dist_step,
        lower = .Machine$double.eps,
        upper = optimize_max_dist,
        null.ok = TRUE
      )
      optimize_dist_step <- pipeline_options$cluster$optimize$dist_step %||% 1e-3
    }

  }
  threshold_file <- pipeline_options$cluster$thresholds
  if (isFALSE(optimize_thresholds)) {
    checkmate::assert_file_exists(threshold_file)
  } else {
    checkmate::assert_path_for_output(threshold_file)
  }
}

#### guilds settings ####
do_guilds <- TRUE
if (!is.null(pipeline_options$guilds)) {
  checkmate::assert_flag(pipeline_options$guilds)
  do_guilds <- pipeline_options$guilds
}

#### OTU table output settings ####
do_dense_otu_table <- FALSE
if (!is.null(pipeline_options$dense_table)) {
  checkmate::assert_flag(pipeline_options$dense_table)
  do_dense_otu_table <- pipeline_options$dense_table
}
