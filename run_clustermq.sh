#!/usr/bin/env bash

#filter and assign taxonomy to demultiplexed Illumina reads

#SBATCH --job-name predcom_test
#SBATCH --account project_2003156
#SBATCH --partition small
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 16G
#SBATCH --time 12:00:00
#SBATCH --mail-type ALL

export PATH="/projappl/project_2003156/GSSP/bin:$PATH"
R --vanilla -f run_clustermq.R
