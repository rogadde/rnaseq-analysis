#!/bin/bash
#
# Author: Rohini Gadde
#   Reference: https://datacarpentry.org/wrangling-genomics/03-trimming/index.html
# Usage: ../bin/trimmomatic.sh (run in data directory of project)
#   This script loops through a directory of .fastq.gz files and trims adapter
#   sequences and low-quality sequences.
# TODO: Make more general for future use. Create file structure instead of hard
# coding existing paths.

for file in fastq/*_1.fastq.gz; do
  read=$(basename "${file}" _1.fastq.gz)
  trimmomatic PE -threads 10 -summary filtered/trimsum.txt \
    fastq/"${read}"_1.fastq.gz fastq/"${read}"_2.fastq.gz \
    filtered/"${read}"_1p.trim.fastq.gz filtered/"${read}"_1u.trim.fastq.gz \
    filtered/"${read}"_2p.trim.fastq.gz filtered/"${read}"_2u.trim.fastq.gz \
    ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:35
done