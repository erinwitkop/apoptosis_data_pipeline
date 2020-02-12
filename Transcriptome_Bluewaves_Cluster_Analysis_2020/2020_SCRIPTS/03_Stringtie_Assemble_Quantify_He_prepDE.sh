#!/bin/bash
#SBATCH -t 1000:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH	-o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/Stringtie_He_prepDE_out_2_11_2020
#SBATCH	-e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/Stringtie_He_prepDE_error_2_11_2020

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

### C_gig_He_2015_OsHV1_SRA_ID
# assemble transcripts for each sample with the GFF3 annotation file
array5=($(ls $GHO/*.bam))
for i in ${array5[@]}; do
	stringtie -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o ${i}.gtf ${i}
	echo "${i} assembled"
	echo "${i}.gtf" >> $GHO/He_OsHV1_mergelist.txt # Make stringtie mergelist with names of all gtf files with full path
done

#Run StringTie merge, merge transcripts from all samples in single experiment
stringtie --merge -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o $GHO/He_OsHV1_stringtie_merged.gtf $GHO/He_OsHV1_mergelist.txt
echo "He OsHV1 merged"
#gffcompare to compare how transcripts compare to reference annotation
gffcompare -r $CG/GCF_000297895.1_oyster_v9_genomic.gff -G -o $GHO/He_OsHV1_stringtie_merged $GHO/He_OsHV1_stringtie_merged.gtf
echo "He OsHV1 gffcompared"
#Re-estimate transcript abundance after merge step
for i in ${array5[@]}; do
		stringtie -A $(echo ${i}|sed "s/\..*//").abd.tab -e -G $GHO/He_OsHV1_stringtie_merged.gtf -o $(echo ${i}|sed "s/\..*//").merge.gtf ${i}
		echo "${i} He OsHV1 transcript abundance re-estimated"
done

echo "He OsHV1 Stringtie complete $(date)"

# C_gig_He_2015_OsHV1_SRA_ID
cd $GHO/
for i in *.merge.gtf; do
	# create text file with sample IDs and respective paths
	echo "$(echo $i |sed "s/\..*//") $GHO/$i" >> C_gig_He_sample_list.txt # keep .. in sed command since these are SE and don't have underscore
done

python prepDE_Oct.2019.py -v  -i C_gig_He_sample_list.txt -g He_gene_count_matrix.csv -t He_transcript_count_matrix.csv
echo "C_gig_He_2015_OsHV1_SRA_ID DONE $(date)"
