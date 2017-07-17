#!/bin/bash
#PBS-l nodes=8
#PBS-l walltime=1000:00:00
#PBS -j oe

# This script uses the SRA toolkit to download C.gigas paired end SRA's using the 
# list of SRA's that are present in the uploaded project files

set -e 
echo "START $(date)"
cd /data3/marine_diseases_lab/erin/Bio_project_SRA/
module load SRA-Toolkit/2.8.2-1-centos_linux64

for f in *_paired.text
do 
  prefetch --option-file $f 
  while read -r LINE; do
    fastq-dump $LINE --split-files --readids --outdir . 
  done < $f 
done   
	


	





echo "STOP $(date)"
