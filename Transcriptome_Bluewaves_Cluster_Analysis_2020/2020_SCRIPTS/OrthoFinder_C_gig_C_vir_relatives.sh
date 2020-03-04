#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH	-o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/Orthofinder_output_3_4_2020
#SBATCH	-e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/Orthofinder_error_3_4_2020

# Load Orthofinder module and DIAMOND
 module load OrthoFinder/2.3.3-foss-2018b-Python-2.7.15
 module load DIAMOND/0.9.25-foss-2018b

# Set working directory
F=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/OrthoFinder_2020

# Input for OrthoFinder is the full protein sequences from the species, decompressed
# Place protein.faa files for each species (Crassostrea gigas, Crassostrea virginica, Biomphalaria glabrata,
  # Mizuhopecten yessoensis, Octopus bimaculoides) in this folder

# From manual
# 3 SIMPLE USAGE: Run full OrthoFinder analysis on FASTA format proteomes in <dir>
  # $ orthofinder [options] -f <dir>

orthofinder.py -f $F/.
