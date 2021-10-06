#!/bin/bash
#
# Author: Rohini Gadde
#   Reference: https://datacarpentry.org/wrangling-genomics/03-trimming/index.html
# Usage: ../bin/trimmomatic.sh $1 $2 (run in data directory of project)
#   This script loops through a directory of .fastq.gz files and trims adapter
#   sequences and low-quality sequences.
# 1 - Number of threads
# 2 - File containing adapter/primer sequences
# TODO: Make more general for future use. Create file structure instead of hard
# coding existing paths.

if [[ $# -ne 2 ]]; then
  echo "Incorrect number of arguments"
  exit 1
fi

for file in fastq/*_1.fastq.gz; do
  read=$(basename "${file}" _1.fastq.gz)
  trimmomatic PE -threads "$1" \
    fastq/"${read}_1.fastq.gz" fastq/"${read}_2.fastq.gz" \
    filtered/"${read}_1P.trim.fastq.gz" filtered/"${read}_1U.trim.fastq.gz" \
    filtered/"${read}_2P.trim.fastq.gz" filtered/"${read}_2U.trim.fastq.gz" \
    ILLUMINACLIP:"$2":2:30:10 SLIDINGWINDOW:4:20 MINLEN:25 \
    2> filtered_logs/"${read}_trim.log"
done