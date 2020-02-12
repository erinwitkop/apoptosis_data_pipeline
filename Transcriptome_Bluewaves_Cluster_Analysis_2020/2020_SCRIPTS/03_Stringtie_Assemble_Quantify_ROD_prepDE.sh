#!/bin/bash
#SBATCH -t 1000:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH	-o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/Stringtie_ROD_prepDE_out_2_11_2020
#SBATCH	-e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/Stringtie_ROD_prepDE_error_2_11_2020

#This script takes bam files from HISAT (processed by SAMtools) and performs StringTie assembly and quantification and converts
# data into a format that is readable as count tables for DESeq2 usage

echo "START" $(date)

module load StringTie/2.1.1-GCCcore-7.3.0 # new version of Stringtie
module load gffcompare/0.11.5-foss-2018b # new version of gffcompare
module load python/2.7.6

# Create variable for each path to trimmed and quality filtered data folder
CV=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/Cvir_Genome_and_Indexes
CG=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/Cgig_Genome_and_Indexes
CP=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/2020_Raw_Transcriptome_Data/C_vir_Probiotic_SRA
CR=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/2020_Raw_Transcriptome_Data/C_vir_ROD_SRA
CD=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/2020_Raw_Transcriptome_Data/Dermo_2015_transcriptomes/Fastq
GLO=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_deLorgeril_OsHV1_SRA
GHO=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_He_OsHV1_SRA
GRV=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA
GZV=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Zhang_Vibrio_SRA

### C_vir_ROD_SRA_ID
# assemble transcripts for each sample with the GFF3 annotation file
array7=($(ls $CR/*.bam))
for i in ${array7[@]}; do
	stringtie -G $CV/ref_C_virginica-3.0_top_level.gff3 -o ${i}.gtf ${i}
	echo "${i} assembled"
	echo "${i}.gtf" >> $CR/ROD_mergelist.txt # Make stringtie mergelist with names of all gtf files with full path
done

#Run StringTie merge, merge transcripts from all samples in single experiment
stringtie --merge -G $CV/ref_C_virginica-3.0_top_level.gff3 -o $CR/ROD_stringtie_merged.gtf $CR/ROD_mergelist.txt
echo "ROD merged"
#gffcompare to compare how transcripts compare to reference annotation
gffcompare -r $CV/ref_C_virginica-3.0_top_level.gff3 -G -o $CR/ROD_stringtie_merged $CR/ROD_stringtie_merged.gtf
echo "ROD gffcompared"
#Re-estimate transcript abundance after merge step
for i in ${array7[@]}; do
		stringtie -A $(echo ${i}|sed "s/\..*//").abd.tab -e -G $CR/ROD_stringtie_merged.gtf -o $(echo ${i}|sed "s/\..*//").merge.gtf ${i}
		echo "${i} ROD transcript abundance re-estimated"
done

# C_vir_ROD_SRA_ID
cd $CR/
for i in *.merge.gtf; do
	# create text file with sample IDs and respective paths
	echo "$(echo $i |sed "s/\..*//") $CR/$i" >> C_vir_ROD_sample_list.txt
done

python prepDE_Oct.2019.py -v -i $CR/C_vir_ROD_sample_list.txt -g ROD_gene_count_matrix.csv -t ROD_transcript_count_matrix.csv
echo "C_vir_ROD_SRA_ID DONE $(date)"
