#!/bin/bash
#SBATCH -t 1000:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH	-o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/Stringtie_ALL_out_2_10_2020
#SBATCH	-e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/Stringtie_ALL_error_2_10_2020

#This script takes bam files from HISAT (processed by SAMtools) and performs StringTie assembly and quantification and converts
# data into a format that is readable as count tables for DESeq2 usage


echo "START" $(date)

module load StringTie/2.1.1-GCCcore-7.3.0 # new version of Stringtie
module load gffcompare/0.11.5-foss-2018b # new version of gffcompare

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

### C_gig_deLorgeril_OsHV1
# assemble transcripts for each sample with the GFF3 annotation file
array1=($(ls $GLO/*.bam))
#for i in ${array1[@]}; do
#	stringtie -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o ${i}.gtf -l $(echo ${i}|sed "s/\..*//") ${i}
#	echo "${i} assembled"
#	echo "${i}.gtf" >> $GLO/deLorgeril_mergelist.txt # Make stringtie mergelist with names of all gtf files with full path
#done
	# command structure: $ stringtie <options> -G <reference.gtf or .gff> -o outputname.gtf -l prefix_for_transcripts input_filename.bam
	# -o specifies the output name
	# -G specifies you are aligning with an option GFF or GTF file as well to perform novel transcript discovery
	# -l Sets <label> as the prefix for the name of the output transcripts. Default: STRG
	# don't use -e here if you want it to assemble any novel transcripts

#Run StringTie merge, merge transcripts from all samples in single experiment
stringtie --merge -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o $GLO/deLorgeril_OsHV1_stringtie_merged.gtf $GLO/deLorgeril_mergelist.txt
echo "deLorgeril OsHV1 merged"
		#-G is a flag saying to use the .gff annotation file

#gffcompare to compare how transcripts compare to reference annotation
gffcompare -r $CG/GCF_000297895.1_oyster_v9_genomic.gff -G -o $GLO/deLorgeril_OsHV1_merged $GLO/deLorgeril_OsHV1_stringtie_merged.gtf
echo "deLorgeril OsHV1 gffcompared"
		# -o specifies prefix to use for output files
		# -r followed by the annotation file to use as a reference
	 	# merged.annotation.gtf tells you how well the predicted transcripts track to the reference annotation file
	 	# merged.stats file shows the sensitivity and precision statistics and total number for different features (genes, exons, transcripts)

#Re-estimate transcript abundance after merge step
for i in ${array1[@]}; do
		stringtie -e -A -G $GLO/deLorgeril_OsHV1_stringtie_merged.gtf -o $(echo ${i}|sed "s/\..*//").merge.gtf ${i}
		echo "${i} deLorgeril OsHV1 transcript abundance re-estimated"
done
		# input here is the original set of alignment files
		#-A here creates a gene table output with genomic locations and compiled information that I will need later to fetch gene sequences
			#FROM MANUAL: "If StringTie is run with the -A <gene_abund.tab> option, it returns a file containing gene abundances. "
		# here -G refers to the merged GTF files
		# -e creates more accurate abundance estimations with input transcripts, needed when converting to DESeq2 tables
		# -e Limits the processing of read alignments to only estimate and output the assembled transcripts matching the
				#reference transcripts given with the -G option (requires -G, recommended for -B/-b). With this option,
				# read bundles with no reference transcripts will be entirely skipped, which may provide a considerable
				# speed boost when the given set of reference transcripts is limited to a set of target genes, for example.

echo "deLorgeril Stringtie complete $(date)"

### C_gig_Rubio_Vibrio_SRA_ID
# assemble transcripts for each sample with the GFF3 annotation file
array2=($(ls $GRV/*.bam))
#for i in ${array2[@]}; do
#	stringtie -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o ${i}.gtf -l $(echo ${i}|sed "s/\..*//") ${i}
#	echo "${i} assembled"
#	echo "${i}.gtf" >> $GRV/Rubio_Vibrio_mergelist.txt # Make stringtie mergelist with names of all gtf files with full path
#done

#Run StringTie merge, merge transcripts from all samples in single experiment
stringtie --merge -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o $GRV/Rubio_Vibrio_stringtie_merged.gtf $GRV/Rubio_Vibrio_mergelist.txt
echo "Rubio Vibrio merged"
#gffcompare to compare how transcripts compare to reference annotation
gffcompare -r $CG/GCF_000297895.1_oyster_v9_genomic.gff -G -o $GRV/Rubio_Vibrio_stringtie_merged $GRV/Rubio_Vibrio_stringtie_merged.gtf
echo "Rubio Vibrio gffcompared"
#Re-estimate transcript abundance after merge step
for i in ${array2[@]}; do
		stringtie -e -A -G $GRV/Rubio_Vibrio_stringtie_merged.gtf -o $(echo ${i}|sed "s/\..*//").merge.gtf ${i}
		echo "${i} Rubio Vibrio transcript abundance re-estimated"
done

echo "Rubio Vibrio Stringtie complete $(date)"

### C_vir_Probiotic_SRA_ID
# assemble transcripts for each sample with the GFF3 annotation file
array3=($(ls $CP/*.bam))
# for i in ${array3[@]}; do
#	stringtie -G $CV/ref_C_virginica-3.0_top_level.gff3  -o ${i}.gtf -l $(echo ${i}|sed "s/\..*//") ${i}
#	echo "${i} assembled"
#	echo "${i}.gtf" >> $CP/Probiotic_mergelist.txt # Make stringtie mergelist with names of all gtf files with full path
#done

#Run StringTie merge, merge transcripts from all samples in single experiment
stringtie --merge -G $CV/ref_C_virginica-3.0_top_level.gff3 -o $CP/Probiotic_stringtie_merged.gtf $CP/Probiotic_mergelist.txt
echo "Probiotic merged"
#gffcompare to compare how transcripts compare to reference annotation
gffcompare -r $CV/ref_C_virginica-3.0_top_level.gff3 -G -o $CP/Probiotic_stringtie_merged $CP/Probiotic_stringtie_merged.gtf
echo "Probiotic gffcompared"
#Re-estimate transcript abundance after merge step
for i in ${array3[@]}; do
		stringtie -e -A -G $CP/Probiotic_stringtie_merged.gtf -o $(echo ${i}|sed "s/\..*//").merge.gtf ${i}
		echo "${i} Probiotic transcript abundance re-estimated"
done

echo "Probiotic Stringtie complete $(date)"

### C_vir_Dermo
# assemble transcripts for each sample with the GFF3 annotation file
array4=($(ls $CD/*.bam))
#for i in ${array4[@]}; do
#	stringtie -G $CV/ref_C_virginica-3.0_top_level.gff3  -o ${i}.gtf -l $(echo ${i}|sed "s/\..*//") ${i}
#	echo "${i} assembled"
#	echo "${i}.gtf" >> $CD/Dermo_mergelist.txt # Make stringtie mergelist with names of all gtf files with full path
#done

#Run StringTie merge, merge transcripts from all samples in single experiment
stringtie --merge -G $CV/ref_C_virginica-3.0_top_level.gff3 -o $CD/Dermo_stringtie_merged.gtf $CD/Dermo_mergelist.txt
echo "Dermo merged"
#gffcompare to compare how transcripts compare to reference annotation
gffcompare -r $CV/ref_C_virginica-3.0_top_level.gff3 -G -o $CD/Dermo_stringtie_merged $CD/Dermo_stringtie_merged.gtf
echo "Dermo gffcompared"
#Re-estimate transcript abundance after merge step
for i in ${array4[@]}; do
		stringtie -e -A -G $CD/Dermo_stringtie_merged.gtf -o $(echo ${i}|sed "s/\..*//").merge.gtf ${i}
		echo "${i} Dermo transcript abundance re-estimated"
done

echo "Dermo Stringtie complete $(date)"

### C_gig_He_2015_OsHV1_SRA_ID
# assemble transcripts for each sample with the GFF3 annotation file
array5=($(ls $GHO/*.bam))
#for i in ${array5[@]}; do
#	stringtie -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o ${i}.gtf -l $(echo ${i}|sed "s/\..*//") ${i}
#	echo "${i} assembled"
#	echo "${i}.gtf" >> $GHO/He_OsHV1_mergelist.txt # Make stringtie mergelist with names of all gtf files with full path
#done

#Run StringTie merge, merge transcripts from all samples in single experiment
stringtie --merge -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o $GHO/He_OsHV1_stringtie_merged.gtf $GHO/He_OsHV1_mergelist.txt
echo "He OsHV1 merged"
#gffcompare to compare how transcripts compare to reference annotation
gffcompare -r $CG/GCF_000297895.1_oyster_v9_genomic.gff -G -o $GHO/He_OsHV1_stringtie_merged $GHO/He_OsHV1_stringtie_merged.gtf
echo "He OsHV1 gffcompared"
#Re-estimate transcript abundance after merge step
for i in ${array5[@]}; do
		stringtie -e -A -G $GHO/He_OsHV1_stringtie_merged.gtf -o $(echo ${i}|sed "s/\..*//").merge.gtf ${i}
		echo "${i} He OsHV1 transcript abundance re-estimated"
done

echo "He OsHV1 Stringtie complete $(date)"

### C_gig_Zhang_Vibrio_SRA_ID
# assemble transcripts for each sample with the GFF3 annotation file
array6=($(ls $GZV/*.bam))
for i in ${array6[@]}; do
	stringtie -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o ${i}.gtf -l $(echo ${i}|sed "s/\..*//") ${i}
	echo "${i} assembled"
	echo "${i}.gtf" >> $GZV/Zhang_Vibrio_mergelist.txt # Make stringtie mergelist with names of all gtf files with full path
done

#Run StringTie merge, merge transcripts from all samples in single experiment
stringtie --merge -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o $GZV/Zhang_Vibrio_stringtie_merged.gtf $GZV/Zhang_Vibrio_mergelist.txt
echo "Zhang Vibrio merged"
#gffcompare to compare how transcripts compare to reference annotation
gffcompare -r $CG/GCF_000297895.1_oyster_v9_genomic.gff -G -o $GZV/Zhang_Vibrio_stringtie_merged $GZV/Zhang_Vibrio_stringtie_merged.gtf
echo "Zhang Vibrio gffcompared"
#Re-estimate transcript abundance after merge step
for i in ${array6[@]}; do
		stringtie -e -A -G $GZV/Zhang_Vibrio_stringtie_merged.gtf -o $(echo ${i}|sed "s/\..*//").merge.gtf ${i}
		echo "${i} Zhang Vibrio transcript abundance re-estimated"
done

echo "Zhang Vibrio Stringtie complete $(date)"

### C_vir_ROD_SRA_ID
# assemble transcripts for each sample with the GFF3 annotation file
array7=($(ls $CR/*.bam))
for i in ${array7[@]}; do
	stringtie -G $CV/ref_C_virginica-3.0_top_level.gff3  -o ${i}.gtf -l $(echo ${i}|sed "s/\..*//") ${i}
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
		stringtie -e -A -G $CR/ROD_stringtie_merged.gtf -o $(echo ${i}|sed "s/\..*//").merge.gtf ${i}
		echo "${i} ROD transcript abundance re-estimated"
done

echo "ROD Stringtie complete $(date)"


echo "DONE ALL $(date)"
