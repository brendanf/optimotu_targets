#!/usr/bin/env bash

#filter and assign taxonomy to demultiplexed Illumina reads

#SBATCH --job-name its2_taxonomy_first
#SBATCH --account project_2003104
#SBATCH --partition small
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --nodes 1
#SBATCH --mem 16G
#SBATCH --time 48:00:00
#SBATCH --mail-type ALL

SINGULARITY_BIND="bin/usearch11.0.667_i86linux32:/sh_matching/programs/usearch"
SINGULARITY_BIND="${SINGULARITY_BIND},sh_matching_pub/sh_matching_analysis/scripts:/sh_matching/scripts"
SINGULARITY_BIND="${SINGULARITY_BIND},sh_matching_pub/sh_matching_analysis/readme.txt:/sh_matching/readme.txt" 
SINGULARITY_BIND="${SINGULARITY_BIND},data/sh_matching_data:/sh_matching/data"
export SINGULARITY_BIND
export PATH="/projappl/project_2003156/its2_taxonomy_first/bin:$PATH"
R --vanilla -f run_clustermq.R
