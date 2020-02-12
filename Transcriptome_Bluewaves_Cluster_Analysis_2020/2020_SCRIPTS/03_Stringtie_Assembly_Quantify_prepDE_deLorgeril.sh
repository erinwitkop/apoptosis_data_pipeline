#!/bin/bash
#SBATCH -t 1000:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH	-o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/Stringtie_deLorgeril_prepDE_out_2_11_2020
#SBATCH	-e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/Stringtie_deLorgeril_prepDE_error_2_11_2020

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

### C_gig_deLorgeril_OsHV1
# assemble transcripts for each sample with the GFF3 annotation file
array1=($(ls $GLO/*.bam))
for i in ${array1[@]}; do
	stringtie -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o ${i}.gtf ${i}
	echo "${i} assembled"
	echo "${i}.gtf" >> $GLO/deLorgeril_mergelist.txt # Make stringtie mergelist with names of all gtf files with full path
done
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
		stringtie -A $(echo ${i}|sed "s/\..*//").abd.tab -e -G $GLO/deLorgeril_OsHV1_stringtie_merged.gtf -o $(echo ${i}|sed "s/\..*//").merge.gtf ${i}
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

## C_gig_deLorgeril_OsHV1
cd $GLO/
for i in *.merge.gtf; do
	# create text file with sample IDs and respective paths
	echo "$(echo $i |sed "s/\..*//") $GLO/$i" >> C_gig_deLorgeril_sample_list.txt
done

python prepDE_Oct.2019.py -v -i C_gig_deLorgeril_sample_list.txt -g deLorgeril_gene_count_matrix.csv -t deLorgeril_transcript_count_matrix.csv

# -i is the the parent directory of the sample sub-directories or a .txt file listing sample IDs and the paths to GTF files in tab-delimited format
# -g where to output the gene count matrix [default: gene_count_matrix.csv
# -t where to output the transcript count matrix [default: transcript_count_matrix.csv]
# -s STRING, --string=STRING	if a different prefix is used for geneIDs assigned by StringTie [default: MSTRG], don't need however since running now without -l

echo "C_gig_deLorgeril_OsHV1 DONE $(date)"
