#!/bin/bash
#PBS-l nodes=3
#PBS-l walltime=1000:00:00
#PBS -j oe

# This script performs just Quality filtering for SE reads using BBtools bbduk.sh, using output from previous step to remove adapters

set -e 
echo "START" $(date)
module load BBMap/37.36-foss-2016b-Java-1.8.0_131

S=/data3/marine_diseases_lab/erin/Bio_project_SRA

#Quality Filtering to get rid of entire low quality reads. maq=10 will trim reads that have average quality of less than 10

#Processing steps for SE reads
array3=($(ls $S/*.fastq.clean.trim))

for i in ${array3[@]}; do 
	bbduk.sh in1=${i} out1=${i}.filter maq=10
	echo "STOP" $(date)
done

echo "STOP" $(date)

