#!/bin/bash
#SBATCH -t 400:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH -o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/extract_IAP_all_seq_5_8_2020
#SBATCH -e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/extract_IAP_all_seq_error_5_8_2020

echo "START $(date)"

# Set paths needed
IAP=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/OrthoFinder_2020/OrthoFinder_Data_Analysis/Results_Mar25/IAP

# script to fetch sequences from all the orthofinder genomes

array1=($(ls $IAP/IAP_ALL_XP_lookup_all_sep.txt)

for i in ${array1[@]}; do
	sed -n '/${i}/,/^>/p' All_genomes_prot.faa | sed '$d' >> $IAP/IAP_all_orthogroups.fa
done

echo "STOP $(date)"
