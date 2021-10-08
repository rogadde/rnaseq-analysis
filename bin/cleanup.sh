#!/bin/bash
#
# Author: Rohini Gadde
# Usage: ../bin/cleanup.sh (run in data directory of project)
#   This script deletes unnecessary files generated earlier in this analysis.

# Delete unsorted bam files
for dir in STAR/SRR*; do
  if [[ -e "${dir}Aligned.out.bam" ]]; then
    rm -f "${dir}Aligned.out.bam"
  fi
done