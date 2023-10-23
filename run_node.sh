#!/usr/bin/env bash

#filter and assign taxonomy to demultiplexed Illumina reads

#SBATCH --job-name LP_ITS
#SBATCH --account project_2005718
#SBATCH --partition small
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --nodes 1
#SBATCH --mem 64G
#SBATCH --time 24:00:00
#SBATCH --mail-type ALL
##SBATCH --gres=nvme:100

export OMP_STACKSIZE=8096
[ -v $SLURM_CPUS_ON_NODE ] && export OMP_THREAD_LIMIT=$SLURM_CPUS_ON_NODE
if [ -d "$LOCAL_SCRATCH" ] ; then
  export TMPDIR=$LOCAL_SCRATCH
  SINGULARITY_BIND="${SINGULARITY_BIND},${LOCAL_SCRATCH}:$(pwd)/userdir"
  export SINGULARITY_BIND
  echo "bind paths: $SINGULARITY_BIND"
fi
export PATH="/projappl/project_2003156/GSSP/bin:$PATH"
if [[ $1 == "test" ]] ; then
if [[ $2 == "" ]] ; then
echo "testing outdated targets..."
R --vanilla --quiet -e 'targets::tar_outdated(callr_function=NULL)'
else
echo "testing outdated targets leading to $2"
R --vanilla --quiet -e "targets::tar_outdated($2, callr_function=NULL)"
fi
elif [[ $1 == "" ]] ; then
echo "building plan"
R --vanilla --quiet -e 'targets::tar_make(callr_function=NULL, reporter="timestamp")'
else
echo "building target '$1'"
R --vanilla --quiet -e "targets::tar_make($1, callr_function=NULL, reporter='timestamp')"
fi
