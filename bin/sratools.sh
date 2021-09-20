#!/bin/bash
# Author: Rohini Gadde
# Usage: ./sratool.sh $1
#	This script iterates over an accession list to download and convert all SRA
# 	files whose accession numbers are in said list
# 1 - Accession list of SRA files

if [ $# -ne 1 ]
then
	echo "Incorrect number of arguments"
	exit 1
fi

while read accession || [ -n "$accession" ]
do
	if [ -e "${accession}.sra" -a -e "${accession}"*.fastq ]
	then
		continue
	else 
		# .sra files output in directory specified in config: ~/ncbi/public/sra
		prefetch "$accession"
		# fasterq-dump should automatically find accession in public user-repository
		fasterq-dump "${accession}" --split-3 --skip-technical --outdir ~/data/Pan_TranscriptomicProfiling/fastq
	fi
done < $1
