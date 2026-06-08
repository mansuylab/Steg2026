#!/bin/bash

# Get input arguments
input_bed="$1"
output_pbc="$2"

# Calculate PBC metrics
cat ${input_bed} | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | \
sort | uniq -c | \
awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} \
END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n", mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${output_pbc}
