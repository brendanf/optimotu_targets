############################################################
## TODO: ##
# - add primer sequence validations check (ATGCRYSKMBDHVIN only)
# - other sanity checks for the loaded settings
############################################################

# don't import the whole package, but let's use the null default operator
`%||%` <- rlang::`%||%`

if (file.exists("pipeline_options.yaml")) {
  pipeline_options <- yaml::read_yaml("pipeline_options.yaml")
} else {
  warning(
    "Options file 'pipeline_options.yaml' is missing!\n",
    "Using defaults for all parameters."
  )
  pipeline_options <- list()
}

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
}

checkmate::assert(
  checkmate::check_null(pipeline_options$added_reference_fasta),
  checkmate::check_file_exists(pipeline_options$added_reference_fasta)
)

checkmate::assert(
  checkmate::check_null(pipeline_options$added_reference_table),
  checkmate::check_file_exists(pipeline_options$added_reference_table)
)

if (xor(is.null(pipeline_options$added_reference_fasta),
        is.null(pipeline_options$added_reference_table))) {
  stop(
    "If one of 'added_reference_fasta' and 'added_reference_table' is given ",
    "in 'pipeline_options.yaml', then both must be given."
    )
}

pipeline_options$custom_sample_table <-
  pipeline_options$custom_sample_table %||% FALSE
checkmate::assert(
  checkmate::check_flag(pipeline_options$custom_sample_table)
)
if (isTRUE(pipeline_options$custom_sample_table)) {
  warning(
    "Pipeline option 'custom_sample_table' is not yet supported.\n",
    "Defaulting to implicit sample table from file names."
  )
}

### cutadapt settings
if (length(pipeline_options$forward_primer) == 0) {
  stop("ERROR: forward primer string missing (file: pipeline_options.yaml)")
} else if (length(pipeline_options$forward_primer) > 1) {
  stop("ERROR: specify only one forward primer (file: pipeline_options.yaml).")
}
if (length(pipeline_options$reverse_primer) == 0) {
  stop("ERROR: reverse primer string missing (file: pipeline_options.yaml)")
} else if (length(pipeline_options$reverse_primer) > 1) {
  stop("ERROR: specify only one reverse primer (file: pipeline_options.yaml).")
}