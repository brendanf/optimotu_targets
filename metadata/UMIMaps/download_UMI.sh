#!/usr/bin/env bash
for N in 00025 00026 00028 {00032..00037} 00041
do
 f=$(rclone ls s3allas:2005718_malaise_meta_1/LIFEPLAN-$N | grep UMI | cut -c11-100)
 rclone copy s3allas:2005718_malaise_meta_1/LIFEPLAN-$N/$f .
 mv $f LIFEPLAN_${N}_UMIMap.txt
done

