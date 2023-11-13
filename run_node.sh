#!/usr/bin/env bash

#filter and assign taxonomy to demultiplexed Illumina reads

#SBATCH --job-name OptimOTU
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
if [ -v SLURM_CPUS_ON_NODE ] ; then
  export OMP_THREAD_LIMIT=$SLURM_CPUS_ON_NODE
fi
if [ -d "$LOCAL_SCRATCH" ] ; then
  export TMPDIR=$LOCAL_SCRATCH
  export SINGULARITY_BIND="${LOCAL_SCRATCH}:$(pwd)/userdir"
  echo "bind paths: $SINGULARITY_BIND"
fi
export PATH="/projappl/project_2003156/OptimOTU/bin:$PATH"
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
