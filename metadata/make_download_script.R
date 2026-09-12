links <- readLines("metadata/download_links.txt")

script <- tibble::tibble(link = links) |>
  tidyr::extract(
    link,
    into = c("lane", "filename"),
    regex = "https://www.dropbox.com/.+/LPLAN([0-9]{2})_NSEQ-[0-9]{5}/([^?]+)\\?.+",
    remove = FALSE
  ) |>
  dplyr::mutate(lane = as.integer(lane)) |>
  with(sprintf(
    'mkdir -p LIFEPLAN-%05d && wget "%s" -O LIFEPLAN-%05d/%s',
    lane, link, lane, filename))

header <- c(
  "#!/usr/bin/env bash",
  "#SBATCH --cpus-per-task=1",
  "#SBATCH --mem=1G",
  "#SBATCH --time=24:00:00",
  "#SBATCH --account=project_2005718",
  "#SBATCH --mail-type=ALL",
  "#SBATCH --partition=small",
  "# Download LIFEPLAN new sequencing data from Dropbox",
  "set -e"
)

writeLines(c(header, script), "metadata/download_script.sh")
