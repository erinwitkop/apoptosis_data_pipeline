#!/bin/bash
#PBS-l nodes=1
#PBS-l walltime=1000:00:00
#PBS -j oe

# This script processes SRA PE and SE with BBtools (bbduk.sh) to perform adaptor trimming.

set -e 
echo "START $(date)"
module load BBMap/37.36-foss-2016b-Java-1.8.0_131

F=/data3/marine_diseases_lab/erin/Bio_project_SRA/PE_fastq
S=/data3/marine_diseases_lab/erin/Bio_project_SRA

# Trimming of adaptors on paired end reads found in the previous 01_C_gigas_RNA_pipeline.sh

array1=($(ls $F/*_1.fq))
array2=($(ls $F/*_2.fq))

for i in ${array1[@]}; do
   bbduk.sh in1=${i} out1=${i}.clean in2=$(echo ${i}|sed s/_1/_2/) out2=$(echo ${i}|sed s/_1/_2/).clean ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
   echo "STOP PE $(date)"
done
	
	#ktrim = r means it will only trim from right side, which is where the adapter should be. (ktrim=l would trim from left)
	#hdist = hamming distance, hdist =1 allows for 1 mismatch
	#flag -tbo specifies to also trim adaptors based on pair overlap detection using BBMerge, which does not require known adapter sequences
	#flag -tpe specified to trim both reads to the same length (if the adapter kmer was only detected in one of them and not other)
	# mink = allows it to use shorter kmers at the ends of the read
	# see manual for additional options: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/
	# instead of the adapters.fa file which has illumina and Nextera adapter sequences, you can put in literal strings for your adapter, literal=ACTGGT,TTTGGTG” for those two literal strings.

# Commands for Single End Read PreProcessing

# Trimming of adaptors from SE reads 
array3=($(ls $S/*.fastq))

for i in ${array3[@]}; do 
	bbduk.sh in1=${i} out1=${i}.clean ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
	echo "STOP SE $(date)"
done


echo "DONE $(date)"
