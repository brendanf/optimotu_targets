############################################################
## TODO: ##
# - other sanity checks for the loaded settings
############################################################

# don't import the whole package, but let's use the null default operator
`%||%` <- rlang::`%||%`

optimotu.pipeline::parse_pipeline_options()

# default global value to replace non-existent targets
if (!optimotu.pipeline::do_model_filter()) {
  full_length_read_counts <- tibble::tibble(sample_key = character())
}

if (!optimotu.pipeline::do_spike()) {
  spike_read_counts <- nospike_read_counts <-
    tibble::tibble(sample_key = character())
}

if (!optimotu.pipeline::do_pos_control()) {
  control_read_counts <- nocontrol_read_counts <-
    tibble::tibble(sample_key = character())
}
