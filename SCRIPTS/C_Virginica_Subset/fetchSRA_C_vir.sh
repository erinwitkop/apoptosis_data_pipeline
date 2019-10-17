#!/bin/bash
#PBS-l nodes=1
#PBS-l walltime=1000:00:00
#PBS -j oe
#PBS -m ae -M erin_roberts@my.uri.edu
#PBS -o out_fetchSRA3


# This script uses the SRA toolkit to download SRA's using the 
# list of SRA's that are present in the uploaded project files

set -e 
echo "START $(date)"
F=/data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files/C_Vir_subset
module load SRA-Toolkit/2.8.2-1-centos_linux64

#For Single End (GX3_F3L SRAs)
#for f in $F/GX3_F3L.txt
#do 
#  prefetch --option-file $f 
#  while read -r LINE; do
#    fastq-dump -O $F --readids $LINE 
#  done < $f 
#done   

#For paired end (UW_transcriptome)
for f in $F/New_Transcriptomes.txt
do 
  prefetch --option-file $f 
  while read -r LINE; do
    fastq-dump -O $F --split-files --readids $LINE
  done < $f 
done 

echo "STOP $(date)"

