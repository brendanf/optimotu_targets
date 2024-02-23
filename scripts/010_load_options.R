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
if (!("project_name" %in% names(pipeline_options))
    || length(pipeline_options$project_name) == 0) {
  warning(
    "Missing project name in 'pipeline_options.yaml'.\n",
    "Using default name 'metabarcoding_project'."
  )
  pipeline_options$project_name <- "metabarcoding_project"
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
project_name <-pipeline_options$project_name

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
checkmate::assert_string(
  pipeline_options$forward_primer,
  null.ok = TRUE,
  min.chars = 10,
  pattern = "[ACGTSWRYMKBDHVIN]+",
  ignore.case = TRUE
)
if (is.null(pipeline_options$forward_primer)) {
  primer_R1 <- "GCATCGATGAAGAACGCAGC"
  message("Forward primer string missing (file: pipeline_options.yaml)\n",
          "Using default: GCATCGATGAAGAACGCAGC")
} else {
  primer_R1 <- pipeline_options$forward_primer
}

checkmate::assert_string(
  pipeline_options$reverse_primer,
  null.ok = TRUE,
  min.chars = 10,
  pattern = "[ACGTSWRYMKBDHVIN]+",
  ignore.case = TRUE
)
if (is.null(pipeline_options$reverse_primer)) {
  primer_R2 <- "TCCTCCGCTTATTGATATGC"
  message("Reverse primer string missing (file: pipeline_options.yaml)\n",
          "Using default: TCCTCCGCTTATTGATATGC")
} else {
  primer_R2 <- pipeline_options$reverse_primer
}

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
