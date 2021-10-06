#!/bin/bash
#
# Author: Rohini Gadde
# Usage: ../bin/star.sh $1 $2 $3 $4 (run in data directory of project)
#   This script aligns trimmed (paired) reads to a reference genome. It assumes
# the genome indices have already been generated.
# 1 - Number of threads
# 2 - Path to accession list
# 3 - Path to genome directory
# 4 - Path to input file directory (no trailing slash)

if [[ $# -ne 4 ]]; then
  echo "Incorrect number of arguments"
  exit 1
fi

thread="$1"
srr_list="$2"
genome_dir="$3"
dir_in="$4"

acc_to_list () {
  acc_list=()
  while read acc; do
    acc_list+=("${acc}")
  done < "${srr_list}"
}

acc_to_list
for ((i=0,j=1; i<${#acc_list[@]}; i+=3,j++)); do
  rep1=${acc_list[$i]}
  rep2=${acc_list[$((i+1))]}
  rep3=${acc_list[$((i+2))]}
  
  STAR --runMode alignReads --runThread "${thread}" --genomeDir "${genome_dir}" \
  --readFilesCommand zcat --readFilesIn \
  ${dir_in}/"${rep1}_1P"*,${dir_in}/"${rep2}_1P"*,${dir_in}/"${rep3}_1P"* \
  ${dir_in}/"${rep1}_2P"*,${dir_in}/"${rep2}_2P"*,${dir_in}/"${rep3}_2P"* \
  --outFileNamePrefix STAR/s${j} --outSAMtype BAM Unsorted SortedByCoordinate \
  --quantMode TranscriptomeSAM
done  