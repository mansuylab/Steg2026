#!/bin/bash

# tracks.sh <input> <output> <log> <binsize> <threads>

# Arguments
INPUT=$1
OUTPUT=$2
LOG=$3
BS=$4
THREADS=$5

# Variables
BNAME=$(basename $INPUT)
NAME=${BNAME%.*}

bamCoverage -b $INPUT -o $OUTPUT$NAME.bw --binSize $BS --normalizeUsing CPM  -p $THREADS 2>$LOG$NAME.log
