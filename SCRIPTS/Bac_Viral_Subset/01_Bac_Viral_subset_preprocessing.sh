#!/bin/bash
#PBS-l nodes=2
#PBS-l walltime=1000:00:00
#PBS -j oe
#PBS -q default
#PBS -o out_Bac_Viral_HISAT
#PBS -e err_Bac_Viral_HISAT
#PBS -m ae -M erin_roberts@my.uri.edu

#01_BAc_Viral_Subset_Preprocessing. This script runs the preprocessing steps on raw SRA reads just for the 
# OsHV-1 pathogen challenge and the Gram positive and Gram negative Bac challenges. This script skips
# the step to use BBMerge to obtain the stats file about which adapters are present. All SRAs are single end.

set -e 
echo "START $(date)"
module load BBMap/37.36-foss-2016b-Java-1.8.0_131

F=/data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files/Bac_Viral_subset

# Commands for Single End Read PreProcessing

#Trimming of adaptors from SE reads 
array1=($(ls $F/*.fastq))

for i in ${array1[@]}; do 
	bbduk.sh in1=${i} out1=${i}.clean ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
	echo "STOP SE $(date)"
done

#ktrim = r means it will only trim from right side, which is where the adapter should be. (ktrim=l would trim from left)
	#hdist = hamming distance, hdist =1 allows for 1 mismatch
	#flag -tbo specifies to also trim adaptors based on pair overlap detection using BBMerge, which does not require known adapter sequences
	#flag -tpe specified to trim both reads to the same length (if the adapter kmer was only detected in one of them and not other)
	# mink = allows it to use shorter kmers at the ends of the read
	# see manual for additional options: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/
	# instead of the adapters.fa file which has illumina and Nextera adapter sequences, you can put in literal strings for your adapter, literal=ACTGGT,TTTGGTG” for those two literal strings.

echo "DONE $(date)"

#Quality Trimming, of both the left and the right sides to get rid of reads that are less than quality 20
array3=($(ls $F/*.clean))

 for i in ${array3[@]}; do 
	bbduk.sh in1=${i} out1=${i}.trim qtrim=rl trimq=20
 	echo "STOP qual trim $(date)"
 done

#Quality Filtering to get rid of entire low quality reads. maq=10 will trim reads that have average quality of less than 10

for i in ${array3[@]}; do 
	bbduk.sh in1=${i}.trim out1=${i}.trim.filter maq=10
	echo "STOP qual filter $(date)"
done

#histogram generation
 for i in ${array3[@]}; do 
	bbduk.sh in1=${i}.trim.filter bhist=${i}.b.hist qhist=${i}.q.hist gchist=${i}.gc.hist lhist=${i}.l.hist gcbins=auto
 	echo "STOP" $(date)
	echo ${i} > ${i}.hist.all 
	echo "bhist" >> ${i}.hist.all
	cat ${i}.b.hist >> ${i}.hist.all
	echo "qhist" >> ${i}.hist.all
    cat ${i}.q.hist	>> ${i}.hist.all
	echo "gchist" >> ${i}.hist.all
    cat ${i}.gc.hist >> ${i}.hist.all
	echo "lhist" >> ${i}.hist.all
	cat ${i}.l.hist >> ${i}.hist.all
 	echo "histo $(date)"
 done

echo "DONE $(date)"