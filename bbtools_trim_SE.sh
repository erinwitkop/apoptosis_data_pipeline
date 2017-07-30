#!/bin/bash
#PBS-l nodes=3
#PBS-l walltime=1000:00:00
#PBS -j oe

# This script processes SRA PE end reads with BBtools to find adaptor sequences
#  with BBmerge, and then uses these for adaptor trimming and quality trimming with bbduk.sh.

set -e 
echo "START $(date)"
module load BBMap/37.36-foss-2016b-Java-1.8.0_131

S=/data3/marine_diseases_lab/erin/Bio_project_SRA

# Commands for Single End Read PreProcessing

array3=($(ls $S/*.fastq))
 
for i in ${array3[@]}; do  # @ symbol tells it to go through each item in the array  
   bbduk.sh in1=${i} k=23 ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa stats=${i}.stat out=${i}.out
   echo "STOP" $(date)
done

  # stats.txt will then list the names of adapter sequences found, and their frequency

# Trimming of adaptors found in the previous command
for i in ${array3[@]}; do 
	bbduk.sh in1=${i} out1=${i}.clean ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
	echo "STOP" $(date)
done

echo "STOP" $(date)
