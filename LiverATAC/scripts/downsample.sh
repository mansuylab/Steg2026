#!/bin/bash

# downsample.sh <input> <output> <reads> 

# Arguments
INPUT=$1
OUTPUT=$2
READS=$3

# Variables
BNAME=$(basename $INPUT)
NAME=${BNAME%.*}

# Commands
FRACTION=$(samtools idxstats $INPUT | cut -f3 | awk -v ct=$READS 'BEGIN {total=0} {total += $1} END {print ct/total}')
samtools view -b -s ${FRACTION} -@ 8 $INPUT > $OUTPUT$BNAME
