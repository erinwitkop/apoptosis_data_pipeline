#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=250GB
#SBATCH --export=NONE
#SBATCH -o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/SCRIPTS/Script_out_error_files/module_export_module_new_9_21_2020_out
#SBATCH -e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/SCRIPTS/Script_out_error_files/module_export_module_new_9_21_2020_error

echo "START $(date)"

# Load R version 3.6.0 (3.6.1 doesn't have any relevant bug fixes)
module load R/3.6.0-intel-2019a

Rscript /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/SCRIPTS/TOMsim_WGCNA_DESeq_cluster.R

echo "DONE $(date)"
