#!/bin/bash
#PBS-l nodes=2
#PBS-l walltime=1000:00:00
#PBS -j oe

# This script processes SRA PE and SE reads with BBTOOLS (bbduk.sh) to identify all the adapter
# sequences in your files (if you don't know what they are). 

set -e 
echo "START $(date)"
module load BBMap/37.36-foss-2016b-Java-1.8.0_131

F=/data3/marine_diseases_lab/erin/Bio_project_SRA/PE_fastq
S=/data3/marine_diseases_lab/erin/Bio_project_SRA

#all files are in the home directory and either have ending 
# _1.fq or _2.fq
#going to make two array variables and then iterate through them as an index


##### USE BBDUK.SH TO GET STATS.TXT FILE TELLING YOU WHICH ADAPTERS AND THEIR FREQUENCY ####

# Commands for Paired End read processing

#All paired end reads in my working directory have the following suffix : either _1.fq or _2.fq
array1=($(ls $F/*_1.fq))
array2=($(ls $F/*_2.fq))

for i in ${array1[@]}; do  # @ symbol tells it to go through each item in the array  
   bbduk.sh in1=${i} in2=$(echo ${i}|sed s/_1/_2/) k=23 ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa stats=${i}.stat out=${i}.out
done

  # stats.txt will then list the names of adapter sequences found, and their frequency
  
  # k=23, “k” specifies the kmer size to use (must be at most the length of the adapters)
  # ref= sets the reference
  # stats= sets output name of the stats files

# Commands for Single End Read Processing

#All single end reads in my directory have the following suffix: *.fastq
array3=($(ls $S/*.fastq))
 
for i in ${array3[@]}; do  # @ symbol tells it to go through each item in the array  
   bbduk.sh in1=${i} k=23 ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa stats=${i}.stat out=${i}.out
   echo "STOP" $(date)
done

  # k=23, “k” specifies the kmer size to use (must be at most the length of the adapters)
  # ref= sets the reference
  # stats= sets output name of the stats files
  
echo "STOP" $(date)
