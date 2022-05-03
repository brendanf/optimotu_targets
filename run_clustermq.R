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
        job_name = "deadwood_restoration_worker",
        cores = 8,
        time = "48:00:00",
        temp_space=0,
        memory=16384
      )
    )
  )
)
options(
  clustermq.scheduler = "slurm",
  clustermq.template = file.path(getwd(), "slurm", "puhti_clustermq.tmpl")
)

targets::tar_make_clustermq(callr_function=NULL, workers = 12, reporter = "timestamp")
