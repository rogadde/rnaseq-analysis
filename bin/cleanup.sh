#!/bin/bash
#
# Author: Rohini Gadde
# Usage: ../bin/cleanup.sh (run in data directory of project)
#   This script deletes unnecessary files generated earlier in this analysis.

# Delete unsorted bam files
for dir in STAR/SRR*/; do
  acc=$(basename "${dir}")
  
  if [[ -e "${dir}${acc}Aligned.out.bam" ]]; then
    rm -f "${dir}${acc}Aligned.out.bam"
  fi
done