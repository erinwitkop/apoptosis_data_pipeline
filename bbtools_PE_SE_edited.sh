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
S=/data3/marine_diseases_lab/erin/Bio_project_SRA

#Commands for Paired End Read Preprocessing

#all files are in the home directory and either have ending 
# _1.fq or _2.fq
#going to make two array variables and then iterate through them as an index
#changing all the file endings to .fq to see if bbduk.sh prefers that file name ending

 #array1=($(ls $F/*_1.fq))
 #array2=($(ls $F/*_2.fq))

 
#for i in ${array1[@]}; do  # @ symbol tells it to go through each item in the array  
   #bbduk.sh in1=${i} in2=$(echo ${i}|sed s/_1/_2/) k=23 ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa stats=${i}.stat out=${i}.out
#done

  # stats.txt will then list the names of adapter sequences found, and their frequency

# Trimming of adaptors found in the previous command
#for i in ${array1[@]}; do 
	#bbduk.sh in1=${i} out1=${i}.clean in2=$(echo ${i}|sed s/_1/_2/) out2=$(echo ${i}|sed s/_1/_2/).clean ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
#done

	#ktrim = r means it will only trim from right side, which is where the adapter should be. (ktrim=l would trim from left)
	#hdist = hamming distance, hdist =1 allows for 1 mismatch
	#flag -tbo specifies to also trim adaptors based on pir overlap detection using BBMerge 
	#which does not require known adapter sequences)
	#flag -tpe specified to trim both reads to the same length (if the adapter kmer was only detected in one of them and not other)

#Quality trimming, of both the left and the right sides to get rid of reads that are less than quality 20
	
 #for i in ${array1[@]}; do 
#	bbduk.sh in1=${i} out1=${i}.trim in2=$(echo ${i}|sed s/_1/_2/) out2=$(echo ${i}|sed s/_1/_2/).trim qtrim=rl trimq=20
 #	echo "STOP" $(date)
 #done

#Quality filtering to get rid of entire low quality reads. maq=10 will trim reads that have average quality of less than 10

#for i in ${array1[@]}; do 
#	bbduk.sh in1=${i}.trim out1=${i}.trim.filter in2=$(echo ${i}|sed s/_1/_2/).trim out2=$(echo ${i}|sed s/_1/_2/).trim.filter maq=10
#	echo "STOP" $(date)
#done

#Histogram generation, only generating for one of the pair (assuming that similar stats will be present). 
#All histogram output contents are combined into one file
 #for i in ${array1[@]}; do 
#	bbduk.sh in1=${i}.trim.filter in2=$(echo ${i}|sed s/_1/_2/).trim.filter bhist=${i}.b.hist qhist=${i}.q.hist gchist=${i}.gc.hist lhist=${i}.l.hist gcbins=auto
 #	cat *${i}*.hist > ${i}.hist.all
 #	echo "STOP" $(date)
 #done
	#lhist = output a read length histogram
	#qhist = per base average quality
	#bhist = output a per-base composition histogram 
	#gchist = output a gc content histogram

# Commands for Single End Read PreProcessing

array3=($(ls $S/*.fastq))
 
for i in ${array3[@]}; do  # @ symbol tells it to go through each item in the array  
   bbduk.sh in1=${i} k=23 ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa stats=${i}.stat out=${i}.out
done

  # stats.txt will then list the names of adapter sequences found, and their frequency

# Trimming of adaptors found in the previous command
for i in ${array3[@]}; do 
	bbduk.sh in1=${i} out1=${i}.clean ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
	echo "STOP" $(date)
done

	#ktrim = r means it will only trim from right side, which is where the adapter should be. (ktrim=l would trim from left)
	#hdist = hamming distance, hdist =1 allows for 1 mismatch
	#flag -tbo specifies to also trim adaptors based on pir overlap detection using BBMerge 
	#which does not require known adapter sequences)
	#flag -tpe specified to trim both reads to the same length (if the adapter kmer was only detected in one of them and not other)


#Quality Trimming, of both the left and the right sides to get rid of reads that are less than quality 20
	
 for i in ${array3[@]}; do 
	bbduk.sh in1=${i}.clean out1=${i}.clean.trim qtrim=rl trimq=20
 	echo "STOP" $(date)
 done


#Quality Filtering to get rid of entire low quality reads. maq=10 will trim reads that have average quality of less than 10

for i in ${array3[@]}; do 
	bbduk.sh in1=${i}.clean.trim out1=${i}.clean.trim.filter maq=10
	echo "STOP" $(date)
done

#histogram generation
 for i in ${array3[@]}; do 
	bbduk.sh in1=${i}.trim.filter bhist=bhist.txt qhist=qhist.txt gchist=gchist.txt aqhist=aqhist.txt lhist=lhist.txt gcbins=auto out=${i}.clean.trim.filter.out
 	echo "STOP" $(date)
 done

echo "STOP" $(date)


