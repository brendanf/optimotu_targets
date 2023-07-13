#!/usr/bin/env bash
mkdir -p $(dirname "$SEQ_ROOT/$1")
rclone cat "$ALLAS_ROOT/$1" |
gzip -dc - |
paste - - - - |
awk '
{
 pool[int(NR%denom)]=$0
}
NR%denom==0{
 n=denom;
 for (i=0;i<n;i++) idx_pool[i]=i;
 srand(seed + int((NR-1)/denom));
 for (i=0;i<numer;i++) {
  j=int(rand()*n);
  print pool[idx_pool[j]];
  n--;
  idx_pool[j] = idx_pool[n];
 }
}
END{
 m = NR%denom;
 if (m != 0) {
  for (i=0;i<m;i++) idx_pool[i]=i+1;
  srand(seed + int((NR-1)/denom));
  n=m;
  for (i=0;i<int(m*numer/denom);i++) {
   j=int(rand()*n);
   print pool[idx_pool[j]];
   n--;
   idx_pool[j] = idx_pool[n];
  }
 }
}' seed=$SAMPLE_REP numer=$SAMPLE_NUMER denom=$SAMPLE_DENOM |
tr "\t" "\n" |
gzip -c - >"$SEQ_ROOT/$1"
