#!/usr/bin/env bash

#filter and assign taxonomy to demultiplexed Illumina reads

#SBATCH --job-name deadwood_restoration
#SBATCH --account project_2003104
#SBATCH --partition small
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --nodes 1
#SBATCH --mem 64G
#SBATCH --time 24:00:00
#SBATCH --mail-type ALL

SINGULARITY_BIND="bin/usearch11.0.667_i86linux32:/sh_matching/programs/usearch"
SINGULARITY_BIND="${SINGULARITY_BIND},sh_matching_pub/sh_matching_analysis/scripts:/sh_matching/scripts"
SINGULARITY_BIND="${SINGULARITY_BIND},sh_matching_pub/sh_matching_analysis/readme.txt:/sh_matching/readme.txt"
SINGULARITY_BIND="${SINGULARITY_BIND},data/sh_matching_data:/sh_matching/data"
export SINGULARITY_BIND
export PATH="$(pwd)/conda/deadwood_restoration/bin:$PATH"
R --no-save -e 'targets::tar_make(callr_function=NULL, reporter="timestamp")'
