#!/usr/bin/env bash
# Download and reproducibly subsample a fastq.gz file
# takes a single command line argument, which is the path to a fastq.gz file inside Allas
# also requires some environmental variables to be set:
# ALLAS_ROOT : path to the bucket in Allas, e.g. "s3allas:2003156_GSSP"
# SEQ_ROOT : path where the subsampled files will be downloaded, e.g. "/local_scratch/GSSP/sequences/01_raw"
#  (Subdirectories below $ALLAS_ROOT are replicated in $SEQ_ROOT.)
# SAMPLE_NUMER : numerator for the fraction of sequences to take, e.g. "32"
# SAMPLE_DENOM : denominator for the fraction of sequences to take, e.g. "1024"
# SAMPLE_REP : replicate index; the same replicate index will always yield the same subsample

# ensure that the output directory exists
mkdir -p $(dirname "$SEQ_ROOT/$1")

# download the file to stdout
rclone cat "$ALLAS_ROOT/$1" |
# decompress
gzip -dc - |
# combine the four lines of each record into a single tab-separated line
paste - - - - |
# actual subsampling happens in awk
awk '
# input separates fields using tabs, but output uses newlines (so it is a valid fastq again)
BEGIN{
 FS="\t"
 OFS="\n"
}
# on every line: add the line to the reservoir
# (this will rewrite the reservoir after the first "denom" sequences)
# also rename the sequences to a (hexadecimal) numerical index
{
 $1=sprintf("@%05x", NR)
 pool[int(NR%denom)]=$0
}
# when the reservoir is full, draw samples and output them
NR%denom==0{
 n=denom;
 # idx_pool is an array of indices which have not yet been samples
 for (i=0;i<n;i++) idx_pool[i]=i;
 # seed the random number generator; different for each rep and each time the reservoir is filled
 srand(seed + int((NR-1)/denom));
 # draw "numer" samples from the reservoir
 for (i=0;i<numer;i++) {
  j=int(rand()*n);
  print pool[idx_pool[j]];
  n--;
  # update the index pool so that there are no repeat indices in the first n elements
  idx_pool[j] = idx_pool[n];
 }
}
# at the end we may probably have an partially filled reservoir
END{
 # m is the number of elements we have in the reservoir
 m = NR%denom;
 if (m != 0) {
  # generate the index pool as above
  for (i=0;i<m;i++) idx_pool[i]=i+1;
  # random seed as above
  srand(seed + int((NR-1)/denom));
  n=m;
  # sample m*numer/denom sequences and output
  for (i=0;i<int(m*numer/denom);i++) {
   j=int(rand()*n);
   print pool[idx_pool[j]];
   n--;
   idx_pool[j] = idx_pool[n];
  }
 }
}' seed=$SAMPLE_REP numer=$SAMPLE_NUMER denom=$SAMPLE_DENOM |
# zip and output the result
gzip -c - >"$SEQ_ROOT/$1"
echo "$(date '+[%Y-%m-%d %H:%M:%S]') finished subsampling $1" >>/dev/stderr
