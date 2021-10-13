#!/bin/bash
#
# Author: Rohini Gadde
# Usage: <path>/rsem.sh $1 $2 $3 $4
#  This script generates isoform and gene expression results, output into
# subdirectories labelled by accession IDs.
# 1 - Path to gtf annotation file
# 2 - Path to genome fasta file
# 3 - Name of reference used
# 4 - Number of threads

if [[ $# -ne 4 ]]; then
  echo "Incorrect number of arguments. See usage."
  exit 1
fi

gtf="$1"
fasta="$2"
ref="$3"
thread="$4"

# Prepare RSEM reference
rsem-prepare-reference --gtf "${gtf}" "${fasta}" "${ref}"

# Calculate expression levels for each alignment. Input path is the output of
# star.sh.
for dir in STAR/SRR*; do
  acc=$(basename "${dir}")
  
  if [[ -d "rsem/${acc}/" ]]; then
    rm -r "rsem/${acc}/" # Clear any existing directories (reset)
  fi
  mkdir "rsem/${acc}/"
  
  rsem-calculate-expression --alignments --paired-end --no-bam-output \
  -p "${thread}" STAR/${acc}/${acc}Aligned.toTranscriptome.out.bam "${ref}" \
  rsem/${acc}/${acc}
done