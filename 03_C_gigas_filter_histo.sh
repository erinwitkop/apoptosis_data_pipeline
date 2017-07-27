#!/bin/bash
#PBS-l nodes=3
#PBS-l walltime=1000:00:00
#PBS -j oe

# This script processes SRA PE reads and SE ends with BBtools to perform quality trimming, quality filtering
# and generates histograms using bbduk.sh.

set -e 
echo "START" $(date)
module load BBMap/37.36-foss-2016b-Java-1.8.0_131

F=/data3/marine_diseases_lab/erin/Bio_project_SRA/PE_fastq
S=/data3/marine_diseases_lab/erin/Bio_project_SRA

 array1=($(ls $F/*_1.fq.clean))
 array2=($(ls $F/*_2.fq.clean))

#Quality trimming, of both the left and the right sides to get rid of reads that are less than quality 20
	
 for i in ${array1[@]}; do 
	bbduk.sh in1=${i} out1=${i}.trim in2=$(echo ${i}|sed s/_1/_2/) out2=$(echo ${i}|sed s/_1/_2/).trim qtrim=rl trimq=20
 	echo "STOP" $(date)
 done

#Quality filtering to get rid of entire low quality reads. maq=10 will trim reads that have average quality of less than 10

for i in ${array1[@]}; do 
	bbduk.sh in1=${i}.trim out1=${i}.trim.filter in2=$(echo ${i}|sed s/_1/_2/).trim out2=$(echo ${i}|sed s/_1/_2/).trim.filter maq=10
	echo "STOP" $(date)
done

#Histogram generation, only generating for one of the pair (assuming that similar stats will be present). 
#All histogram output contents are combined into one file
 for i in ${array1[@]}; do 
	bbduk.sh in1=${i}.trim.filter in2=$(echo ${i}|sed s/_1/_2/).trim.filter bhist=${i}.b.hist qhist=${i}.q.hist gchist=${i}.gc.hist lhist=${i}.l.hist gcbins=auto
 	cat *${i}*.hist > ${i}.hist.all
 	echo "STOP" $(date)
 done
	#lhist = output a read length histogram
	#qhist = per base average quality
	#bhist = output a per-base composition histogram 
	#gchist = output a gc content histogram


#Processing steps for SE reads
array3=($(ls $S/*.fastq))

#Quality Trimming, of both the left and the right sides to get rid of reads that are less than quality 20
	
 for i in ${array3[@]}; do 
	bbduk.sh in1=${i}.clean out1=${i}.trim qtrim=rl trimq=20
 	echo "STOP" $(date)
 done


#Quality Filtering to get rid of entire low quality reads. maq=10 will trim reads that have average quality of less than 10

for i in ${array3[@]}; do 
	bbduk.sh in1=${i}.trim out1=${i}.trim.filter maq=10
	echo "STOP" $(date)
done

#histogram generation
 for i in ${array3[@]}; do 
	bbduk.sh in1=${i}.trim.filter bhist=bhist.txt qhist=qhist.txt gchist=gchist.txt aqhist=aqhist.txt lhist=lhist.txt gcbins=auto
 	cat *${i}*.hist > ${i}.hist.all
 	echo "STOP" $(date)
 done

echo "STOP" $(date)

echo "STOP" $(date)