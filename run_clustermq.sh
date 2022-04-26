#!/usr/bin/env bash

#filter and assign taxonomy to demultiplexed Illumina reads

#SBATCH --job-name deadwood_restoration
#SBATCH --account project_2003104
#SBATCH --partition small
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --nodes 1
#SBATCH --mem 16G
#SBATCH --time 48:00:00
#SBATCH --mail-type ALL

srun singularity exec -B /scratch deadwood_restoration_bioinformatics.sif R --no-save -f run_clustermq.R
