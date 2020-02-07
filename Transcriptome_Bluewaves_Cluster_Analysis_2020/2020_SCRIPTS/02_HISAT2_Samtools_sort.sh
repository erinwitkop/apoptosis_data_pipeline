#!/bin/bash
#SBATCH -t 1000:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH	-o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/HISAT2_deLorgeril_Rubio_Pro_out_2_6_2020
#SBATCH	-e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/HISAT2_deLorgeril_Rubio_Pro_error_2_6_2020

echo "START" $(date)

module load HISAT2/2.1.0-foss-2018b
module load SAMtools/1.9-foss-2018b

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

############## BUILDING HISAT2 INDEXES ###################
#Index will be made with reference genome and no annotation file (allowing for novel transcript discovery)
	# create new directory for the HISAT index called genome, and put the genome inside it
	# copy all reads files into this directory as well to ensure easy access by commands

# C. virginica genome index
#Build HISAT index with cvir_edited (this file has extra spaces in header removed so that genome and annotation don't conflict)
#hisat2-build -f $CV/cvir_edited.fa $CV/cvir_edited_index
  # -f indicates that the reference input files are FASTA files

# C. gigas genome index
#hisat2-build -f $CG/GCF_000297895.1_oyster_v9_genomic.fna   $CG/GCF_000297895.1_oyster_v9_genomic_index
  # -f indicates that the reference input files are FASTA files

############# USE HISAT TO ALIGN RNA-READS TO GENOME ##############

####### PE read experiments ######
# PE experiments: C_gig_deLorgeril_OsHV1, C_gig_Rubio_Vibrio_SRA_ID, C_vir_Probiotic_SRA_ID, C_vir_Dermo

# C_gig_deLorgeril_OsHV1
array1=($(ls $GLO/*_1.fastq.gz.clean.trim.filter.gz))
for i in ${array1[@]}; do
    # outputs a single bam file
  	hisat2 --dta -x $CG/GCF_000297895.1_oyster_v9_genomic_index -1 ${i} -2 $(echo ${i}|sed s/_1/_2/) -S ${i}.sam
    echo "HISAT2 PE ${i} $(date)"
    #SAMTOOLS sort to convert the SAM file into a BAM file to be used with StringTie. Stringtie only take sorted bam
    samtools sort ${i}.sam > ${i}.bam
    #Get bam file statistics for percentage aligned with flagstat
    samtools flagstat ${i}.bam > ${i}.bam.stats #get % mapped
    echo "${i} sorted bam done"
done

# -x before the index
# -U before the file to be aligned
# --dta : Report alignments tailored for transcript assemblers including StringTie.
     # With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites.
     #This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computation and memory usage.

echo "C_gig_deLorgeril_OsHV1 DONE $(date)"

# C_gig_Rubio_Vibrio_SRA_ID
array2=($(ls $GRV/*_1.fastq.gz.clean.trim.filter.gz))
for i in ${array2[@]}; do
  # outputs a single bam file
	hisat2 --dta -x $CG/GCF_000297895.1_oyster_v9_genomic_index -1 ${i} -2 $(echo ${i}|sed s/_1/_2/) -S ${i}.sam
	echo "HISAT2 PE ${i} $(date)"
  #SAMTOOLS sort to convert the SAM file into a BAM file to be used with StringTie. Stringtie only take sorted bam
  samtools sort ${i}.sam > ${i}.bam
  #Get bam file statistics for percentage aligned with flagstat
  samtools flagstat ${i}.bam > ${i}.bam.stats #get % mapped
  echo "${i} sorted bam done"
done

echo "C_gig_Rubio_Vibrio_SRA_ID DONE $(date)"

# C_vir_Probiotic_SRA_ID
array3=($(ls $CP/*_1.fastq.gz.clean.trim.filter.gz))
for i in ${array3[@]}; do
  # outputs a single bam file
	hisat2 --dta -x $CV/cvir_edited_index  -1 ${i} -2 $(echo ${i}|sed s/_1/_2/) -S ${i}.sam
	echo "HISAT2 PE ${i} $(date)"
  #SAMTOOLS sort to convert the SAM file into a BAM file to be used with StringTie. Stringtie only take sorted bam
  samtools sort ${i}.sam > ${i}.bam
  #Get bam file statistics for percentage aligned with flagstat
  samtools flagstat ${i}.bam > ${i}.bam.stats #get % mapped
 echo "${i} sorted bam done"
done

echo "C_vir_Probiotic_SRA_ID DONE $(date)"

# C_vir_Dermo
#array4=($(ls $CD/*R1.fastq.gz.clean.trim.filter.gz))
#for i in ${array4[@]}; do
  # outputs a single bam file
#	hisat2 --dta -x $CV/cvir_edited_index  -1 ${i} -2 $(echo ${i}|sed s/R1/R2/) -S ${i}.sam
#	echo "HISAT2 PE ${i} $(date)"
  #SAMTOOLS sort to convert the SAM file into a BAM file to be used with StringTie. Stringtie only take sorted bam
#  samtools sort ${i}.sam > ${i}.bam
  #Get bam file statistics for percentage aligned with flagstat
#  samtools flagstat ${i}.bam > ${i}.bam.stats #get % mapped
#  echo "${i} sorted bam done"
#done

# echo "C_vir_Dermo DONE $(date)"

####### SE read experiments ######
# SE experiments: C_gig_He_2015_OsHV1_SRA_ID, C_gig_Zhang_Vibrio_SRA_ID, C_vir_ROD_SRA_ID

# C_gig_He_2015_OsHV1_SRA_ID
#array5=($(ls $GHO/*.filter.gz))
#for i in ${array5[@]}; do
#        hisat2 --dta -x $CG/GCF_000297895.1_oyster_v9_genomic_index -U ${i} -S ${i}.sam
        #SAMTOOLS sort to convert the SAM file into a BAM file to be used with StringTie. Stringtie only take sorted bam
#        samtools sort ${i}.sam > ${i}.bam
        #Get bam file statistics for percentage aligned with flagstat
#        samtools flagstat ${i}.bam > ${i}.bam.stats #get % mapped
#        echo "${i} sorted bam done"
#done

#echo "C_gig_He_2015_OsHV1_SRA_ID DONE $(date)"

# C_gig_Zhang_Vibrio_SRA_ID
#array6=($(ls $GZV/*.filter.gz))
#for i in ${array6[@]}; do
#        hisat2 --dta -x $CG/GCF_000297895.1_oyster_v9_genomic_index -U ${i} -S ${i}.sam
        #SAMTOOLS sort to convert the SAM file into a BAM file to be used with StringTie. Stringtie only take sorted bam
#        samtools sort ${i}.sam > ${i}.bam
        #Get bam file statistics for percentage aligned with flagstat
#        samtools flagstat ${i}.bam > ${i}.bam.stats #get % mapped
#        echo "${i} sorted bam done"
#done

#echo "C_gig_Zhang_Vibrio_SRA_ID DONE $(date)"

# C_vir_ROD_SRA_ID
#array7=($(ls $CR/*.filter.gz))
#for i in ${array7[@]}; do
#        hisat2 --dta -x $CV/cvir_edited_index -U ${i} -S ${i}.sam
        #SAMTOOLS sort to convert the SAM file into a BAM file to be used with StringTie. Stringtie only take sorted bam
#        samtools sort ${i}.sam > ${i}.bam
        #Get bam file statistics for percentage aligned with flagstat
#        samtools flagstat ${i}.bam > ${i}.bam.stats #get % mapped
#        echo "${i} sorted bam done"
#done

#echo "C_vir_ROD_SRA_ID DONE $(date)"

#echo "full script DONE $(date)"

#reference: Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie, and Ballgown
#https://sequencing.qcfail.com/articles/mapq-values-are-really-useful-but-their-implementation-is-a-mess/
#http://www.htslib.org/doc/samtools.html
