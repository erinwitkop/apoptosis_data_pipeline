#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH -o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/SCRIPTS/Script_out_error_files/Dermo_tol_Pro_RE22_net_TOM_out
#SBATCH -e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/SCRIPTS/Script_out_error_files/Dermo_tol_Pro_RE22_net_TOM_error

echo "START $(date)"

# Load R version 3.6.0 (3.6.1 doesn't have any relevant bug fixes)
module load R/3.6.0-intel-2019a


Rscript /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/SCRIPTS/TOMsim_cluster.R

echo "DONE $(date)"
