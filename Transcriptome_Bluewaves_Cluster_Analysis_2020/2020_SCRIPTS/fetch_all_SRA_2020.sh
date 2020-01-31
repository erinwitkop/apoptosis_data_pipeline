#!/bin/bash
#SBATCH -t 1000:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH	-o fetch_SRA_output_1_31_2020
#SBATCH	-e fetch_SRA_error_1_31_2020
#SBATCH	--mail-user=erin_roberts@my.uri.edu

echo "START $(date)"

# Load modules
module load SRA-Toolkit/2.9.0-centos_linux64

# Create variable for each path to raw data folder for each species

C=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/2020_Raw_Transcriptome_Data/
G=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/

####### Load single end read files #######

# Load C_gig_He_2015_OsHV1_SRA_ID
for f in $G/C_gig_He_OsHV1_SRA/C_gig_He_2015_OsHV1_SRA_ID.txt
do
  prefetch --option-file $f
  while read -r LINE; do
    fastq-dump -O $F --readids $LINE
  done < $f
done

echo "C_gig_He_2015_OsHV1_SRA_ID DONE $(date)"

## Load C_gig_Zhang_Vibrio_SRA_ID
for f in $G/C_gig_Zhang_Vibrio_SRA/C_gig_He_2015_OsHV1_SRA_ID.txt
do
  prefetch --option-file $f
  while read -r LINE; do
    fastq-dump -O $F --readids $LINE
  done < $f
done

echo "C_gig_Zhang_Vibrio_SRA_ID DONE $(date)"

## Load C_vir_ROD_SRA_ID
for f in $C/C_vir_ROD_SRA/C_vir_ROD_SRA_ID.txt
do
  prefetch --option-file $f
  while read -r LINE; do
    fastq-dump -O $F --readids $LINE
  done < $f
done

echo "C_vir_ROD_SRA_ID DONE $(date)"


###### Load paired end read files #####

## Load C_gig_deLorgeril_OsHV1
for f in $G/C_gig_deLorgeril_OsHV1_SRA/C_gig_deLorgeril_OsHV1_SRA_ID.txt
do
  prefetch --option-file $f
  while read -r LINE; do
    fastq-dump -O $F --split-files --readids $LINE
  done < $f
done

echo "C_gig_deLorgeril_OsHV1 DONE $(date)"

## Load C_gig_Rubio_Vibrio_SRA_ID
for f in $G/C_gig_Rubio_Vibrio_SRA/C_gig_Rubio_Vibrio_SRA_ID.txt
do
  prefetch --option-file $f
  while read -r LINE; do
    fastq-dump -O $F --split-files --readids $LINE
  done < $f
done

echo "C_gig_Rubio_Vibrio_SRA DONE $(date)"

## Load C_vir_Probiotic_SRA_ID
for f in $C/C_vir_Probiotic_SRA/C_vir_Probiotic_SRA_ID.txt
do
  prefetch --option-file $f
  while read -r LINE; do
    fastq-dump -O $F --split-files --readids $LINE
  done < $f
done

echo "C_vir_Probiotic_SRA_ID DONE $(date)"


echo "STOP $(date)"
