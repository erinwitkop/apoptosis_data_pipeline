#!/bin/bash
#PBS-l nodes=3
#PBS-l walltime=1000:00:00
#PBS -j oe

#Processing steps for SE reads. This script just performs Quality trimming and then stops
# not currently up to date
set -e
echo "START" $(date)
module load BBMap/37.36-foss-2016b-Java-1.8.0_131


S=/data3/marine_diseases_lab/erin/Bio_project_SRA

array3=($(ls $S/*.fastq.clean))

#Quality Trimming, of both the left and the right sides to get rid of reads that are less than quality 20

for i in ${array3[@]}; do
        bbduk.sh in1=${i} out1=${i}.trim qtrim=rl trimq=20
        echo "STOP" $(date)
done

echo "STOP" $(date)