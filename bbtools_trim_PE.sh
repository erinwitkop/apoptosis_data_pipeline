#!/bin/bash
#PBS-l nodes=1
#PBS-l walltime=1000:00:00
#PBS -j oe

# This script processes SRA PE end reads with BBtools to find adaptor sequences
#  with BBmerge, and then uses these for adaptor trimming and quality trimming with bbduk.sh.

set -e 
echo "START $(date)"
module load BBMap/37.36-foss-2016b-Java-1.8.0_131

#all files are in the home directory and either have ending 
# _1.fastq or _2.fastq
#going to make two array variables and then iterate through them as an index
#changing all the file endings to .fq to see if bbduk.sh prefers that file name ending

array=(
   ERR498201_1.fq
   SRR2601695_1.fq
   SRR2601706_1.fq 
   SRR334213_1.fq
   SRR334339_1.fq
   SRR623231_1.fq
   SRR1060328_1.fq
   SRR2601696_1.fq
   SRR2601707_1.fq
   SRR334214_1.fq
   SRR334340_1.fq
   SRR826503_1.fq
   SRR2601482_1.fq
   SRR2601698_1.fq 
   SRR2601709_1.fq
   SRR334215_1.fq
   SRR497890_1.fq
   SRR957603_1.fq
   SRR2601571_1.fq 
   SRR2601699_1.fq
   SRR2601711_1.fq 
   SRR334216_1.fq  
   SRR497891_1.fq
   SRR960386_1.fq
   SRR2601666_1.fq
   SRR2601701_1.fq
   SRR2601712_1.fq
   SRR334217_1.fq
   SRR5043896_1.fq
   SRR967068_1.fq
   SRR2601689_1.fq  
   SRR2601702_1.fq
   SRR2601714_1.fq
   SRR334218_1.fq 
   SRR5043897_1.fq
   SRR2601690_1.fq
   SRR2601703_1.fq
   SRR2601716_1.fq
   SRR334219_1.fq
   SRR5043898_1.fq
   SRR2601692_1.fq
   SRR2601704_1.fq
   SRR2601718_1.fq
   SRR334220_1.fq
   SRR5043899_1.fq
   SRR2601694_1.fq
   SRR2601705_1.fq
   SRR334212_1.fq
   SRR334221_1.fq
   SRR526975_1.fq
)

array2=(
   ERR498201_2.fq   
   SRR2601695_2.fq  
   SRR2601706_2.fq  
   SRR334213_2.fq  
   SRR334339_2.fq   
   SRR623231_2.fq
   SRR1060328_2.fq  
   SRR2601696_2.fq  
   SRR2601707_2.fq  
   SRR334214_2.fq  
   SRR334340_2.fq  
   SRR826503_2.fq
   SRR2601482_2.fq 
   SRR2601698_2.fq  
   SRR2601709_2.fq  
   SRR334215_2.fq  
   SRR497890_2.fq   
   SRR957603_2.fq
   SRR2601571_2.fq  
   SRR2601699_2.fq  
   SRR2601711_2.fq  
   SRR334216_2.fq  
   SRR497891_2.fq   
   SRR960386_2.fq
   SRR2601666_2.fq  
   SRR2601701_2.fq  
   SRR2601712_2.fq  
   SRR334217_2.fq  
   SRR5043896_2.fq  
   SRR967068_2.fq
   SRR2601689_2.fq  
   SRR2601702_2.fq  
   SRR2601714_2.fq  
   SRR334218_2.fq  
   SRR5043897_2.fq
   SRR2601690_2.fq  
   SRR2601703_2.fq  
   SRR2601716_2.fq  
   SRR334219_2.fq  
   SRR5043898_2.fq
   SRR2601692_2.fq  
   SRR2601704_2.fq  
   SRR2601718_2.fq  
   SRR334220_2.fq  
   SRR5043899_2.fq
   SRR2601694_2.fq  
   SRR2601705_2.fq  
   SRR334212_2.fq   
   SRR334221_2.fq  
   SRR526975_2.fq
)
 
array3=(
)
for index in ${!array[*]}; do    
   bbduk.sh in1=${array[$index]} in2=${array2[$index]} k=23 ref=adapters.fa stats=stats${array[$index]}.txt out=stdout${array[$index]}.fq
done

# stats.txt will then list the names of adapter sequences found, and their frequency

echo "STOP $(date)"
