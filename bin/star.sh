#!/bin/bash
#
# Author: Rohini Gadde
# Usage: <path>/star.sh $1 $2 $3
#   This script aligns trimmed (paired) reads to a reference genome. It assumes
# the genome indices have already been generated. Alignments are output into 
# subdirectories labelled by accession IDs.
# 1 - Number of threads
# 2 - Path to input file directory (no trailing slash)
# 3 - Path to genome directory

if [[ $# -ne 3 ]]; then
  echo "Incorrect number of arguments. See usage."
  exit 1
fi

thread="$1"
dir_in="$2"
genome_dir="$3"

# Loop over trimmed reads, aligning each read pair to the reference genome
for file in "${dir_in}"/*_1P.trim.fastq.gz; do
  read1=$(basename "${file}")
  acc=$(basename "${read1}" _1P.trim.fastq.gz)
  read2=$(basename "${file}" _1P.trim.fastq.gz)
  read2+="_2P.trim.fastq.gz"
  
  if [[ -d "STAR/${acc}/" ]]; then
    rm -r "STAR/${acc}/" # Remove directory if it already exists
  fi
  mkdir "STAR/${acc}/"
  
  STAR --runMode alignReads --runThreadN "${thread}" --genomeDir "${genome_dir}" \
  --readFilesCommand zcat --genomeLoad LoadAndKeep --limitBAMsortRAM 150000000000 \
  --readFilesIn "${dir_in}/${read1}" "${dir_in}/${read2}" \
  --outFileNamePrefix "STAR/${acc}/${acc}" \
  --outSAMtype BAM Unsorted SortedByCoordinate \
  --quantMode TranscriptomeSAM
done  