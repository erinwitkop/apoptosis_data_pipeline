#!/bin/bash
#SBATCH -t 1000:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH	-o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/HISAT_index_build_out_2_4_2020
#SBATCH	-e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/HISAT_index_build_error_2_4_2020

echo "START" $(date)

module load HISAT2/2.1.0-foss-2016b

# Create variable for each path to trimmed and quality filtered data folder
CV=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/Cvir_Genome_and_Indexes
CG=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/Cgig_Genome_and_Indexes

############## BUILDING HISAT2 INDEXES ###################
#Index will be made with reference genome and no annotation file (allowing for novel transcript discovery)
	# create new directory for the HISAT index called genome, and put the genome inside it
	# copy all reads files into this directory as well to ensure easy access by commands

# C. virginica genome index
#Build HISAT index with cvir_edited (this file has extra spaces in header removed so that genome and annotation don't conflict)
hisat2-build -f $CV/cvir_edited.fa $CV/cvir_edited_index
  # -f indicates that the reference input files are FASTA files

# C. gigas genome index
hisat2-build -f $CG/GCF_000297895.1_oyster_v9_genomic.fna   $CG/GCF_000297895.1_oyster_v9_genomic_index
  # -f indicates that the reference input files are FASTA files
