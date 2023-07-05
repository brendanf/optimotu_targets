#!/usr/bin/env bash

#SBATCH --job-name gssp_rarify
#SBATCH --account project_2003156
#SBATCH --partition small
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --array=11
#SBATCH --mem 32G
#SBATCH --time 8:00:00
#SBATCH --mail-type ALL

N_RARIFY=13
SAMPLE_NUM=$(bc <<<"2^(${SLURM_ARRAY_TASK_ID}%${N_RARIFY})")
SAMPLE_DENOM=$(bc <<<"2^${N_RARIFY}")
SAMPLE_REP=$(bc <<<"${SLURM_ARRAY_TASK_ID}/${N_RARIFY}")

export OMP_STACKSIZE=8096
export OMP_THREAD_LIMIT=$SLURM_CPUS_PER_TASK

OLD_DIR=$(pwd)
GSSP_ROOT=$LOCAL_SCRATCH/GSSP
mkdir -p "$GSSP_ROOT"
tar -xzf GSSP_optimotu.tar.gz -C "$GSSP_ROOT"
cp -r ../protaxFungi "$LOCAL_SCRATCH/"

wget -O https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz |
 gunzip -c - >"$GSSP_ROOT/bin/usearch"
chmod +x "$GSSP_ROOT/bin/usearch"

wget https://files.plutof.ut.ee/public/orig/9C/FD/9CFD7C58956E5331F1497853359E874DEB639B17B04DB264C8828D04FA964A8F.zip
unzip -j 9CFD7C58956E5331F1497853359E874DEB639B17B04DB264C8828D04FA964A8F.zip data/shs_out.txt data/sanger_refs_sh.fasta -d "$GSSP_ROOT/data/sh_matching_data"
rm 9CFD7C58956E5331F1497853359E874DEB639B17B04DB264C8828D04FA964A8F.zip

seqdir="$GSSP_ROOT/sequences/01_raw"
for f in $(ls sequences/01_raw/**/*.fastq.gz)
do
 mkdir -p $(dirname "f")
 zcat "$f" |
 paste - - - - |
 awk '{pool[int(NR%denom)]=$0};
 NR%denom==0{
  n=denom;
  for (i=0;i<n;i++) idx_pool[i]=i;
  srand(seed + int((NR-1)/denom));
  for (i=0;i<num;i++) {
   j=int(rand()*n);
   print pool[idx_pool[j]];
   n--;
   idx_pool[j] = idx_pool[n];
  }
 };
 END{
  m = NR%denom;
  if (m != 0) {
   for (i=0;i<m;i++) idx_pool[i]=i+1;
   srand(seed + int((NR-1)/denom));
   n=m;
   for (i=0;i<int(m*num/denom);i++) {
    j=int(rand()*n);
    print pool[idx_pool[j]];
    n--;
    idx_pool[j] = idx_pool[n];
   }
  }
 }' seed=$SAMPLE_REP num=$SAMPLE_NUM denom=$SAMPLE_DENOM |
 tr "\t" "\n" |
 gzip -c - >"$GSSP_ROOT/$f"
done

cd "$GSSP_ROOT"
export PATH="/projappl/project_2003156/its2_taxonomy_first/bin:$PATH"
R --vanilla --quiet -e 'targets::tar_make(callr_function=NULL, reporter="timestamp")'

RARE_DIR="$OLD_DIR/rarified"
mkdir -p "$RARE_DIR"
tar -czf "$RARE_DIR/rarify_output_rep${SAMPLE_REP}_${SAMPLE_NUM}_per_${SAMPLE_DENOM}.tar.gz" output



