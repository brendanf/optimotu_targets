#!/usr/bin/env bash

#filter and assign taxonomy to demultiplexed Illumina reads

#SBATCH --job-name deadwood_restoration
#SBATCH --account project_2003104
#SBATCH --partition small
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --nodes 1
#SBATCH --mem 16G
#SBATCH --time 24:00:00
#SBATCH --mail-type ALL

export PATH=$(pwd)/conda/deadwood_restoration/bin:$PATH
srun R --no-save -f run_clustermq.R
