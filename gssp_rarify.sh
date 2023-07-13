#!/usr/bin/env bash

#SBATCH --job-name gssp_rarify
#SBATCH --account project_2003156
#SBATCH --partition small
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --array=0
#SBATCH --mem 32G
#SBATCH --time 8:00:00
#SBATCH --gres=nvme:300
#SBATCH --mail-type ALL

#13 rarifactions at powers of 2, i.e. 1/2, 1/4, ..., 1/8192
export N_RARIFY=13
export SAMPLE_NUMER=$(bc <<<"2^(${SLURM_ARRAY_TASK_ID}%${N_RARIFY})")
export SAMPLE_DENOM=$(bc <<<"2^${N_RARIFY}")
export SAMPLE_REP=$(bc <<<"${SLURM_ARRAY_TASK_ID}/${N_RARIFY}")

export OMP_STACKSIZE=8096
export OMP_THREAD_LIMIT=$SLURM_CPUS_PER_TASK

OLD_DIR=$(pwd)
export GSSP_ROOT=$LOCAL_SCRATCH/GSSP
mkdir -p "$GSSP_ROOT"
tar -xzf GSSP-optimotu.tar.gz -C "$GSSP_ROOT"
cp -r ../protaxFungi "$LOCAL_SCRATCH/"

wget -nv -O - https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz |
 gunzip -c - >"$GSSP_ROOT/bin/usearch"
chmod +x "$GSSP_ROOT/bin/usearch"

wget -nv https://files.plutof.ut.ee/public/orig/9C/FD/9CFD7C58956E5331F1497853359E874DEB639B17B04DB264C8828D04FA964A8F.zip
unzip -j 9CFD7C58956E5331F1497853359E874DEB639B17B04DB264C8828D04FA964A8F.zip data/shs_out.txt data/sanger_refs_sh.fasta -d "$GSSP_ROOT/data/sh_matching_data"
rm 9CFD7C58956E5331F1497853359E874DEB639B17B04DB264C8828D04FA964A8F.zip

export ALLAS_ROOT="s3allas:2003156_GSSP"
export SEQ_ROOT="$GSSP_ROOT/sequences/01_raw"
source /appl/opt/csc-cli-utils/allas-cli-utils/allas_conf -f -k --mode s3cmd $OS_PROJECT_NAME
rclone lsf -R "$ALLAS_ROOT" |
grep '\.fastq\.gz$' |
xargs -l -P$SLURM_CPUS_PER_TASK ./seq_sample.sh

cd "$GSSP_ROOT"
export PATH="/projappl/project_2003156/its2_taxonomy_first/bin:$PATH"
R --vanilla --quiet -e 'targets::tar_make(callr_function=NULL, reporter="timestamp")'

RARE_DIR="$OLD_DIR/rarified"
mkdir -p "$RARE_DIR"
tar -czf "$RARE_DIR/rarify_output_rep${SAMPLE_REP}_${SAMPLE_NUMER}_per_${SAMPLE_DENOM}.tar.gz" output
