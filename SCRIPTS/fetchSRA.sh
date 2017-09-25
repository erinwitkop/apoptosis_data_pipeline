#!/bin/bash
#PBS-l nodes=1
#PBS-l walltime=1000:00:00
#PBS -j oe

# This script uses the SRA toolkit to download C.gigas SRA's into fasta format sequences
#	using the list of SRA Accessions downloaded from NCBI that are present in the uploaded project files.

set -e 
echo "START $(date)"
F=/data3/marine_diseases_lab/erin/Bio_project_SRA/
module load SRA-Toolkit/2.8.2-1-centos_linux64

for f in $F/*.txt
do 
  prefetch --option-file $f 
  while read -r LINE; do
    fastq-dump -O $F --readids $LINE
  done < $f 
done   
	
#for paired end reads

for f in $F/*_paired.text
do 
  prefetch --option-file $f 
  while read -r LINE; do
    fastq-dump -O $F --split-files --readids $LINE 
  done < $f 
done   

echo "STOP $(date)"
