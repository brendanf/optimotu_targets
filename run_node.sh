#!/usr/bin/env bash

#filter and assign taxonomy to demultiplexed Illumina reads

#SBATCH --job-name spruce-log-fungi
#SBATCH --account project_2003104
#SBATCH --partition small
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --nodes 1
#SBATCH --mem-per-cpu 8G
#SBATCH --time 1:00:00
#SBATCH --mail-type ALL

module load r-env-singularity/4.0.5
srun singularity_wrapper exec Rscript scripts/dada.R