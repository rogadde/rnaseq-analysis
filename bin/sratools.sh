#!/bin/bash
#
# Author: Rohini Gadde
# Usage: ./sratool.sh $1 $2
#   This script iterates over an accession list to download and convert all
#   SRA files whose accession numbers are in said list
# $1 - Accession list of SRA files
# $2 - Number of threads to use for fasterq-dump (dflt is 6)
# TODO: Automate script more for future use (e.g. create file structure)
# TODO: Make sure script (specifically fasterq-dump) works regardless of path 
# (symlinked or not)

if [[ $# -ne 2 ]]; then
  echo "Incorrect number of arguments"
  exit 1
fi

acc_list="$1"
thread="$2"

# Clear log before each run
if [[ -e "accession_log.txt" ]]; then
  rm -f ../data/accession_log.txt
fi

while read accession; do
  if [[ -e "${accession}"*.fastq ]]; then
    continue
  else 
    # .sra files output in directory specified in config: ~/ncbi/public/sra
    # fasterq-dump should automatically find accession in public user-repository
    { prefetch "${accession}" \
      && vdb-validate "${accession}" \
      && fasterq-dump "${accession}" --split-files -e "${thread}" --outdir ~/data/Pan2020/fastq/ -p \
      && gzip ../data/fastq/"${accession}"*.fastq \
      && rm -f ../data/sra/"${accession}.sra"; } \
      &>> ../data/accession_log.txt
    if [[ $? -ne 0 ]]; then
      echo "Error: Failed to download and convert ${accession} successfully" \
        >> ../data/accession_log.txt
      echo "{accession}" >> ../data/failed_acc.txt
    fi 
  fi
done < "${acc_list}"

# Delete symlink once all sra files have been converted 
unlink ../data/sra
