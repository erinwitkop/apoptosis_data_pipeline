#!/bin/bash
#PBS-l nodes=1
#PBS-l walltime=1000:00:00
#PBS -j oe
#PBS -q default
#PBS -o out_BBtools_preprocessing2
#PBS -e err_BBtools_preprocessing2
#PBS -m ae -M erin_roberts@my.uri.edu


#01_C_Vir_Bac_Subset_Preprocessing. This script runs the preprocessing steps on raw SRA transcriptomes
# from GX3 and F3L Single end transcriptomes and UW transcriptomes with RI probiotic challenge. This script skips
# the step to use BBMerge to obtain the stats file about which adapters are present.

set -e 
echo "START $(date)"
module load BBMap/37.36-foss-2016b-Java-1.8.0_131

F=/data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files/C_Vir_subset

#Commands for Paired End Read Preprocessing, all files are in the home directory and either have ending 
# _1.fq or _2.fq
#going to make two array variables and then iterate through them as an index
#changing PE file endings to fq to make array work more simply

array1=($(ls $F/*_1.fq))
array2=($(ls $F/*_2.fq))

 
#for i in ${array1[@]}; do  # @ symbol tells it to go through each item in the array  
   #bbduk.sh in1=${i} in2=$(echo ${i}|sed s/_1/_2/) k=23 ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa stats=${i}.stat out=${i}.out
#done

  # stats.txt will then list the names of adapter sequences found, and their frequency

#Trimming of adaptors found in the previous command
for i in ${array1[@]}; do 
	bbduk.sh in1=${i} out1=${i}.clean in2=$(echo ${i}|sed s/_1/_2/) out2=$(echo ${i}|sed s/_1/_2/).clean ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
	echo "adapter trimming ${i}" $(date)
done

	#ktrim = r means it will only trim from right side, which is where the adapter should be. (ktrim=l would trim from left)
	#hdist = hamming distance, hdist =1 allows for 1 mismatch
	#flag -tbo specifies to also trim adaptors based on pir overlap detection using BBMerge 
	#which does not require known adapter sequences)
	#flag -tpe specified to trim both reads to the same length (if the adapter kmer was only detected in one of them and not other)

#Quality trimming, of both the left and the right sides to get rid of reads that are less than quality 20
for i in ${array1[@]}; do 
	bbduk.sh in1=${i}.clean out1=${i}.clean.trim in2=$(echo ${i}|sed s/_1/_2/).clean out2=$(echo ${i}|sed s/_1/_2/).clean.trim qtrim=rl trimq=20
 	echo "quality trimming ${i}" $(date)
done

#Quality filtering to get rid of entire low quality reads. maq=10 will trim reads that have average quality of less than 10
for i in ${array1[@]}; do 
	bbduk.sh in1=${i}.clean.trim out1=${i}.clean.trim.filter in2=$(echo ${i}|sed s/_1/_2/).clean.trim out2=$(echo ${i}|sed s/_1/_2/).clean.trim.filter maq=10
	echo "STOP" $(date)
	echo "quality filtering ${i}" $(date)
done

#Histogram generation, only generating for one of the pair (assuming that similar stats will be present). 
#All histogram output contents are combined into one file
for i in ${array1[@]}; do
 	 bbduk.sh in1=${i}.clean.trim.filter in2=$(echo ${i}|sed s/_1/_2/).clean.trim.filter  bhist=${i}.b.hist qhist=${i}.q.hist gchist=${i}.gc.hist lhist=${i}.l.hist gcbins=auto
	 echo "STOP" $(date)
     echo ${i} > ${i}.hist.all
     echo "bhist" >> ${i}.hist.all
     cat ${i}.b.hist >> ${i}.hist.all
     echo "qhist" >> ${i}.hist.all
     cat ${i}.q.hist >> ${i}.hist.all
     echo "gchist" >> ${i}.hist.all
     cat ${i}.gc.hist >> ${i}.hist.all
     echo "lhist" >> ${i}.hist.all
     cat ${i}.l.hist >> ${i}.hist.all 
	 echo "histogram DONE" $(date)
done
		#lhist = output a read length histogram
        #qhist = per base average quality
        #bhist = output a per-base composition histogram
        #gchist = output a gc content histogram

# Commands for Single End Read PreProcessing

array3=($(ls $F/*.fastq))
#for i in ${array3[@]}; do  # @ symbol tells it to go through each item in the array  
#   bbduk.sh in1=${i} k=23 ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa stats=${i}.stat out=${i}.out  
#done

  # stats.txt will then list the names of adapter sequences found, and their frequency

# Trimming of adaptors found in the previous command
for i in ${array3[@]}; do 
	bbduk.sh in1=${i} out1=${i}.clean ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
	echo "adapter trimming {i}" $(date)
done

	#ktrim = r means it will only trim from right side, which is where the adapter should be. (ktrim=l would trim from left)
	#hdist = hamming distance, hdist =1 allows for 1 mismatch
	#flag -tbo specifies to also trim adaptors based on pir overlap detection using BBMerge 
	#which does not require known adapter sequences)
	#flag -tpe specified to trim both reads to the same length (if the adapter kmer was only detected in one of them and not other)

#Quality Trimming, of both the left and the right sides to get rid of reads that are less than quality 20
	
 for i in ${array3[@]}; do 
	bbduk.sh in1=${i}.clean out1=${i}.clean.trim qtrim=rl trimq=20
 	echo "quality trimming ${i}" $(date)
 done


#Quality Filtering to get rid of entire low quality reads. maq=10 will trim reads that have average quality of less than 10

for i in ${array3[@]}; do 
	bbduk.sh in1=${i}.clean.trim out1=${i}.clean.trim.filter maq=10
	echo "quality filtering ${i}" $(date)
done

#histogram generation
 #Histogram generation
for i in ${array3[@]}; do
        bbduk.sh in1=${i}.clean.trim.filter bhist=${i}.b.hist qhist=${i}.q.hist gchist=${i}.gc.hist lhist=${i}.l.hist gcbins=auto
        echo "STOP" $(date)
        echo ${i} > ${i}.hist.all
        echo "bhist" >> ${i}.hist.all
        cat ${i}.b.hist >> ${i}.hist.all
        echo "qhist" >> ${i}.hist.all
        cat ${i}.q.hist >> ${i}.hist.all
        echo "gchist" >> ${i}.hist.all
        cat ${i}.gc.hist >> ${i}.hist.all
        echo "lhist" >> ${i}.hist.all
        cat ${i}.l.hist >> ${i}.hist.all
        echo "STOP" $(date)
		echo "histogram ${i}" $(date)
done

echo "DONE" $(date)


