#!/bin/bash
#PBS-l nodes=2
#PBS-l walltime=1000:00:00
#PBS -j oe

# This script processes SRA PE end reads with BBtools to find adaptor sequences
#  with BBmerge, and then uses these for adaptor trimming and quality trimming with bbduk.sh.

set -e 
echo "START $(date)"
module load BBMap/37.36-foss-2016b-Java-1.8.0_131

F=/data3/marine_diseases_lab/erin/Bio_project_SRA/PE_fastq

#all files are in the home directory and either have ending 
# _1.fq or _2.fq
#going to make two array variables and then iterate through them as an index
#changing all the file endings to .fq to see if bbduk.sh prefers that file name ending

 array1=($(ls $F/*_1.fq))
 array2=($(ls $F/*_2.fq))

 
for i in ${array1[@]}; do  # @ symbol tells it to go through each item in the array  
   bbduk.sh in1=${i} in2=$(echo ${i}|sed s/_1/_2/) k=23 ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa stats=${i}.stat out=${i}.out
done

  # stats.txt will then list the names of adapter sequences found, and their frequency

# Trimming of adaptors found in the previous command

for i in ${array1[@]}; do 
	bbduk.sh in1=${i} out1=${i}.clean in2=$(echo ${i}|sed s/_1/_2/) out2=$(echo ${i}|sed s/_1/_2/).clean ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
done

	#ktrim = r means it will only trim from right side, which is where the adapter should be. (ktrim=l would trim from left)
	#hdist = hamming distance, hdist =1 allows for 1 mismatch
	#flag -tbo specifies to also trim adaptors based on pir overlap detection using BBMerge 
	#which does not require known adapter sequences)
	#flag -tpe specified to trim both reads to the same length (if the adapter kmer was only detected in one of them and not other)


#quality trimming, of both the left and the right sides to get rid of reads that are less than quality 10
	
 for i in ${array1[@]}; do 
	bbduk.sh in1=${i}.clean out1=${i}.trim in2=$(echo ${i}|sed s/_1/_2/).clean out2=$(echo ${i}|sed s/_1/_2/).trim  qtrim=rl trimq=10
 done


#quality filtering to get rid of entire low quality reads. maq=10 will trim reads that have average quality of less than 10

for i in ${array1[@]}; do 
	bbduk.sh in1=${i}.trim out1=${i}.trim.filter in2=$(echo ${i}|sed s/_1/_2/).trim out2=$(echo ${i}|sed s/_1/_2/).trim.filter  maq=10
 done

#histogram generation
 for i in ${array1[@]}; do 
	bbduk.sh in1=${i}.trim.filter in2=$(echo ${i}|sed s/_1/_2/).trim.filter bhist=bhist.txt qhist=qhist.txt gchist=gchist.txt aqhist=aqhist.txt lhist=lhist.txt gcbins=auto
done

# 06 tophat2 and bowtie2-build in 05_C_gigas_RNA_pipeline.sh

# 07 Salmon to perform quantification

#RSEM