targets::tar_option_set(
  # by default, workers run the targets, retrieve their prerequisites, and
  # store their results
  deployment = "worker",
  storage = "worker",
  retrieval = "worker",
  error = "abridge",
  resources = targets::tar_resources(
    clustermq = targets::tar_resources_clustermq(
      template = list(
        job_name = "priority_effects_worker",
        cores = 8,
        time = "12:00:00",
        temp_space=10,
        memory=16384
      )
    )
  )
)
options(
  clustermq.scheduler = "slurm",
  clustermq.template = file.path(getwd(), "slurm", "puhti_clustermq.tmpl")
)

n_seqrun_dir <- length(list.dirs("sequences/01_raw", recursive = FALSE))
targets::tar_make_clustermq(callr_function=NULL, workers = n_seqrun_dir, reporter = "timestamp")
