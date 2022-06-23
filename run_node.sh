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
##SBATCH --gres=nvme:100

export OMP_STACKSIZE=8096
export OMP_THREAD_LIMIT=$SLURM_CPUS_ON_NODE
mkdir -p userdir
SINGULARITY_BIND="sh_matching_pub/sh_matching_analysis/scripts:/sh_matching/scripts"
SINGULARITY_BIND="${SINGULARITY_BIND},sh_matching_pub/sh_matching_analysis/readme.txt:/sh_matching/readme.txt"
SINGULARITY_BIND="${SINGULARITY_BIND},data/sh_matching_data:/sh_matching/data"
SINGULARITY_BIND="${SINGULARITY_BIND},bin:/sh_matching/programs"
if [ -d "$LOCAL_SCRATCH" ] ; then
  export TMPDIR=$LOCAL_SCRATCH
  SINGULARITY_BIND="${SINGULARITY_BIND},${LOCAL_SCRATCH}:$(pwd)/userdir"
fi
export SINGULARITY_BIND
echo "bind paths: $SINGULARITY_BIND"
export PATH="$(pwd)/conda/deadwood_restoration/bin:$PATH"
R --vanilla -e 'targets::tar_make(callr_function=NULL, reporter="timestamp")'
