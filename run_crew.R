n_seqrun_dir <- length(list.dirs("sequences/01_raw", recursive = FALSE))

targets::tar_option_set(
  # by default, workers run the targets, retrieve their prerequisites, and
  # store their results
  deployment = "worker",
  storage = "worker",
  retrieval = "worker",
  error = "abridge",
  controller = crew.cluster::crew_controller_slurm(
    name = "OptimOTU_crew",
    seconds_launch = 7200,
    workers = n_seqrun_dir*4,
    tasks_max = 1000,
    garbage_collection = TRUE,
    launch_max = 3,
    verbose = TRUE,
    script_lines = readLines("slurm/puhti_crew.tmpl"),
    slurm_log_output = "crew-%A.out",
    slurm_log_error = NULL,
    slurm_memory_gigabytes_per_cpu = 16,
    slurm_cpus_per_task = 8,
    slurm_time_minutes = 12*60,
    slurm_partition = "small",
    host = Sys.info()["nodename"]
  )
)

target = strsplit(Sys.getenv("OPTIMOTU_TARGET"), "[, ;]")[[1]]
if (length(target) > 0) {
  targets::tar_make(any_of(target), callr_function=NULL, reporter = "timestamp")
} else {
  targets::tar_make(callr_function=NULL, reporter = "timestamp")
}
