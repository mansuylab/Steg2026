#!/bin/bash

# filter.sh <input> <output> <log>

# Arguments
INPUT=$1
OUTPUT=$2
LOG=$3

# Variables
BNAME=$(basename $INPUT)
NAME=${BNAME%.*}

# Function
macs2 callpeak -f BAMPE -g mm --keep-dup auto --nomodel --nolambda -n $NAME -t $INPUT --outdir $OUTPUT 2> $LOG$NAME.log
