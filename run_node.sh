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
#SBATCH --gres=nvme:10

export OMP_STACKSIZE=8096
if [ -v SLURM_CPUS_ON_NODE ] ; then
  export OMP_THREAD_LIMIT=$SLURM_CPUS_ON_NODE
fi
if [ -d "$LOCAL_SCRATCH" ] ; then
  export TMPDIR=$LOCAL_SCRATCH
fi
if [ -L protaxFungi ] ; then
  export SINGULARITY_BIND="$(realpath protaxFungi),$SINGULARITY_BIND"
fi

echo "bind paths: $SINGULARITY_BIND"

CONTAINER="/projappl/project_2005718/bin/OptimOTU_v7.sif"

if [[ $1 == "test" ]] ; then
 if [[ $2 == "" ]] ; then
  echo "Testing outdated targets..."
  $CONTAINER  R --vanilla --quiet --no-echo -e 'targets::tar_outdated(callr_function=NULL)'
 elif [[ $2 == starts_with\(*\) ]] ; then
  echo "Testing outdated targets matching $2"
  $CONTAINER R --vanilla --quiet --no-echo -e "targets::tar_outdated($2, callr_function=NULL)"
 else
  echo "Testing outdated targets leading to $2"
  $CONTAINER R --vanilla --quiet --no-echo -e "targets::tar_outdated(any_of(strsplit('$2', '[ ,;]')[[1]]), callr_function=NULL)"
 fi
elif [[ $1 == "" ]] ; then
 echo "Building plan on local machine"
 $CONTAINER R --vanilla --quiet --no-echo -e 'targets::tar_make(callr_function=NULL, reporter="timestamp")'
elif [[ $1 == starts_with\(*\) ]] ; then
 echo "Building targets matching $1 and their dependencies on local machine"
 $CONTAINER R --vanilla --quiet --no-echo -e "targets::tar_make($1, callr_function=NULL, reporter='timestamp')"
else
 echo "Building target(s) '$1' on local machine"
 $CONTAINER R --vanilla --quiet --no-echo -e "targets::tar_make(any_of(strsplit('$1', '[ ,;]')[[1]]), callr_function=NULL, reporter='timestamp')"
fi
