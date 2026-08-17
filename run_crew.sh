#!/usr/bin/env bash

#filter and assign taxonomy to demultiplexed Illumina reads

#SBATCH --job-name OptimOTU
#SBATCH --account project_2005718
#SBATCH --partition small
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 16G
#SBATCH --time 72:00:00
#SBATCH --mail-type ALL

# map additional bind points so that sbatch can be run from inside the container
export APPTAINER_BIND=/usr/bin/sbatch,/usr/lib64/slurm,/usr/lib/slurm,/run/munge,/usr/lib64/libmunge.so.2:/usr/lib/libmunge.so.2,/usr/lib64/libargos-toml.so:/usr/lib/libargos-toml.so,/etc/slurm,/etc/passwd,/etc/group,/var/spool/slurmd

if [ -L protaxFungi ] ; then
  export APPTAINER_BIND="$(realpath protaxFungi),$APPTAINER_BIND"
fi

echo "bind paths: $APPTAINER_BIND"

shopt -s expand_aliases
alias container_R="apptainer exec /projappl/project_2005718/bin/OptimOTU_v7.sif R"
alias

if [[ $1 == "test" ]] ; then
if [[ $2 == "" ]] ; then
echo "testing outdated targets..."
echo "NOTE: crew is not used for testing, you could have used 'run_node.sh'"
container_R --vanilla --quiet --no-echo -e 'targets::tar_outdated(callr_function=NULL)'
else
echo "testing outdated targets leading to $2"
echo "NOTE: crew is not used for testing, you could have used 'run_node.sh'"
container_R --vanilla --quiet --no-echo -e "targets::tar_outdated($2, callr_function=NULL)"
fi
elif [[ $1 == "" ]] ; then
echo "Building plan using crew"
container_R --vanilla --quiet --no-echo -f run_crew.R
else
echo "Building target '$1' using crew"
OPTIMOTU_TARGET="$1" container_R --vanilla --quiet --no-echo -f run_crew.R
fi
