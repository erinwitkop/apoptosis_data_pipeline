#!/bin/bash
#PBS-l nodes=1
#PBS-l walltime=1000:00:00
#PBS -j oe

# This script uses the SRA toolkit to download C.gigas SRA's using the 
# list of SRA's that are present in the uploaded project files

set -e 
echo "START $(date)"
cd /data3/marine_diseases_lab/erin/Bio_project_SRA/
module load SRA-Toolkit/2.8.2-1-centos_linux64

for f in *.txt
do 
  prefetch --option-file $f 
  while read -r LINE; do
    fastq-dump $LINE --readids --outdir . 
  done < $f 
done   
	


echo "STOP $(date)"
