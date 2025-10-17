# define functions and metadata for the plan
for (f in list.files("scripts", "^0[[:digit:]]{2}_.+[.]R$", full.names = TRUE)) {
  source(f)
}

# controller which requests workers for targets which need a lot of memory,
# can take advantage of internal parallelism, or both
controller_wide <- crew.cluster::crew_controller_slurm(
  name = "wide",
  seconds_launch = 7200,
  workers = optimotu.pipeline::n_workers(),
  tasks_max = 1000,
  seconds_idle = 120,
  garbage_collection = TRUE,
  options_cluster = crew.cluster::crew_options_slurm(
    verbose = TRUE,
    script_lines = readLines("slurm/puhti_crew.tmpl"),
    log_output = "crew_wide-%A.out",
    log_error = NULL,
    memory_gigabytes_per_cpu = 16,
    cpus_per_task = 10,
    time_minutes = 12*60,
    partition = "small"
  ),
  host = Sys.info()["nodename"]
)

# controller for workers which use a single core with moderate memory
# requirements
controller_thin <- crew.cluster::crew_controller_slurm(
  name = "thin",
  seconds_launch = 7200,
  workers = optimotu.pipeline::n_workers(),
  tasks_max = 1000,
  seconds_idle = 120,
  garbage_collection = TRUE,
#  launch_max = 3,
  options_cluster = crew.cluster::crew_options_slurm(
    verbose = TRUE,
    script_lines = readLines("slurm/puhti_crew.tmpl"),
    log_output = "crew_thin-%A.out",
    log_error = NULL,
    memory_gigabytes_per_cpu = 4,
    cpus_per_task = 1,
    time_minutes = 12*60,
    partition = "small"
  ),
  host = Sys.info()["nodename"]
)

targets::tar_option_set(
  # by default, workers run the targets, retrieve their prerequisites, and
  # store their results
  deployment = "worker",
  storage = "worker",
  retrieval = "worker",
  error = "abridge",
  controller = crew::crew_controller_group(controller_wide, controller_thin)
)

cat("Running pipeline with a pool of at most", optimotu.pipeline::n_workers(), "crew workers.\n")

target = strsplit(Sys.getenv("OPTIMOTU_TARGET"), "[, ;]")[[1]]
if (length(target) > 0) {
  targets::tar_make(any_of(target), callr_function=NULL, reporter = "timestamp")
} else {
  targets::tar_make(callr_function=NULL, reporter = "timestamp")
}
