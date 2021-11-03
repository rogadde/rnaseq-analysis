#!/bin/bash
#
# Author: Rohini Gadde
# Usage: ../bin/quant_matrix.sh (run in data directory of project)
#   This script loops through count information for each RSEM run and compiles 
# the TPM counts in a single matrix. Intended for single-use (hence hard-coded).

cut -f 1 rsem/SRR10079579/SRR10079579.genes.results > rsem/temp.txt 

for dir in rsem/SRR*; do
  acc=$(basename "${dir}")
 
  # Get TPM values
  cut -f 6 "rsem/${acc}/${acc}.genes.results" > "rsem/temp_${acc}.txt"
done

# First column contains gene_ids while following columns contain TPM values for
# each accession in numerical order. Run IDs must be added in Excel; rows and
# columns must be switched before running prcomp()
paste rsem/temp.txt rsem/temp_SRR* > rsem/temp_matrix.txt

# Remove labels
sed '1d' rsem/temp_matrix.txt > rsem/quant_matrix.txt

# Remove temporary files
rm rsem/temp*
