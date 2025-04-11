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
optimotu.pipeline::parse_project_name(pipeline_options)

#### custom_sample_table ####
checkmate::assert(
  checkmate::check_null(pipeline_options$custom_sample_table),
  checkmate::check_false(pipeline_options$custom_sample_table),
  checkmate::check_file_exists(pipeline_options$custom_sample_table)
)
pipeline_options$custom_sample_table <-
  pipeline_options$custom_sample_table %||% FALSE

#### parallelism ####
checkmate::assert_count(pipeline_options$local_threads, positive = TRUE, null.ok = TRUE)
if (!is.null(pipeline_options$local_threads)) {
  options(optimotu_num_threads = pipeline_options$local_threads)
}

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
    optimotu.pipeline::unnest_yaml_list(pipeline_options$added_reference)
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
optimotu.pipeline::parse_forward_primer(pipeline_options)
optimotu.pipeline::parse_reverse_primer(pipeline_options)

#### primer trim settings ####
checkmate::assert_list(pipeline_options$trimming, null.ok = TRUE)
if (is.null(pipeline_options$trimming)) {
  message("No 'trimming' options given in 'pipeline_options.yaml'\n",
          "Using defaults.")
  trim_options <- optimotu.pipeline::cutadapt_paired_options()
} else {
  trim_options <- do.call(
    optimotu.pipeline::cutadapt_paired_options,
    optimotu.pipeline::unnest_yaml_list(pipeline_options$trimming)
  )
}

#### filtering settings ####
dada2_maxEE <- optimotu.pipeline::dada2_filter_options(2, 2)
checkmate::assert_list(pipeline_options$filtering, null.ok = TRUE)
if (is.null(pipeline_options$filtering)) {
  message("No 'filtering' options given in 'pipeline_options.yaml'\n",
          "Using defaults.")
} else {
  pipeline_options$filtering <- optimotu.pipeline::unnest_yaml_list(pipeline_options$filtering)
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
    update(dada2_maxEE, list(maxEE_R1 = pipeline_options$filtering$maxEE_R1))

  checkmate::assert_number(
    pipeline_options$filtering$maxEE_R2,
    lower = 0,
    finite = TRUE,
    null.ok = TRUE
  )
  if (!is.null(pipeline_options$filtering$maxEE_R2))
    update(dada2_maxEE, list(maxEE_R2 = pipeline_options$filtering$maxEE_R2))
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
  pipeline_options$tag_jump <- optimotu.pipeline::unnest_yaml_list(pipeline_options$tag_jump)
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
if (!is.null(pipeline_options$amplicon_model)) {
  optimotu.pipeline::parse_amplicon_model_options(pipeline_options$amplicon_model)
}
if (!optimotu.pipeline::do_model_filter()) {
  full_length_read_counts <- tibble::tibble(sample_key = character())
}

#### control sequence settings ####
do_spike <- FALSE
spike_file <- NULL
do_pos_control <- FALSE
pos_control_file <- NULL
spike_read_counts <- nospike_read_counts <- control_read_counts <-
  nocontrol_read_counts <- tibble::tibble(sample_key = character())
if (!is.null(pipeline_options$control)) {
  checkmate::assert_list(pipeline_options$control)
  pipeline_options$control <- optimotu.pipeline::unnest_yaml_list(pipeline_options$control)
  checkmate::assert_names(
    names(pipeline_options$control),
    subset.of = c("spike", "positive")
  )

  if ("spike" %in% names(pipeline_options$control)) {
    checkmate::assert(
      checkmate::check_file_exists(pipeline_options$control$spike),
      checkmate::check_flag(pipeline_options$control$spike, null.ok = TRUE)
    )
    if (isTRUE(pipeline_options$control$spike)) {
      stop("Option 'control':'spike' should be a file path, evaluate to FALSE,",
           "or be left blank")
    }
    if (is.character(pipeline_options$control$spike)) {
      spike_file <- pipeline_options$control$spike
      do_spike <- TRUE
      remove(spike_read_counts, nospike_read_counts)
    }
  }
  if ("positive" %in% names(pipeline_options$control)) {
    checkmate::assert(
      checkmate::check_file_exists(pipeline_options$control$positive),
      checkmate::check_flag(pipeline_options$control$positive, null.ok = TRUE)
    )
    if (isTRUE(pipeline_options$control$positive)) {
      stop("Option 'control':'positive' should be a file path, evaluate to FALSE,",
           "or be left blank")
    }
    if (is.character(pipeline_options$control$positive)) {
      pos_control_file <- pipeline_options$control$positive
      do_pos_control <- TRUE
      remove(control_read_counts, nocontrol_read_counts)
    }
  }
}

#### protax settings ####
if (!is.null(pipeline_options$protax)) {
  optimotu.pipeline::parse_protax_options(pipeline_options$protax)
}

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


#### clustering settings ####
threshold_file <- "metadata/GSSP_thresholds.tsv"
if (!is.null(pipeline_options$cluster_thresholds)) {
  checkmate::assert_file_exists(pipeline_options$cluster_thresholds)
  threshold_file <- pipeline_options$cluster_thresholds
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
