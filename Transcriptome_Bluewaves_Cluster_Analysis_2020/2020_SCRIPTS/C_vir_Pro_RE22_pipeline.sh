#!/bin/bash
#SBATCH -t 200:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH	-o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/Probiotic_RE22_pipeline_output_3_2_2020
#SBATCH	-e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/Probiotic_RE22_pipeline_ALL_error_3_2_2020

echo "START $(date)"

# Create variable for each path to raw data folder for each species
CPV=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/2020_Raw_Transcriptome_Data/C_vir_Pro_RE22_SRA
CV=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/Cvir_Genome_and_Indexes

##### DOWNLOAD PAIRED END TRANSCRIPTS #######

# Load modules
module load SRA-Toolkit/2.9.0-centos_linux64

## Load C_vir_Pro_RE22
for f in $CPV/Modak_Pro_RE22_SRA_ID.txt
do
  prefetch --option-file $f
  while read -r LINE; do
    fastq-dump -O $CPV/ --split-files --readids --gzip $LINE
  done < $f
done

vdb-validate --option-file $CPV/C_gig_deLorgeril_OsHV1_SRA/Modak_Pro_RE22_SRA_ID.txt &>  $CPV/C_vir_Pro_RE22_sra_checksum.txt

echo "C_vir_Pro_RE22 download DONE $(date)"

# remove the prefetch files from home
rm -r /home/erin_roberts/ncbi

##### SRA TRIM FILTER #####
module purge
module load BBMap/37.36-foss-2016b-Java-1.8.0_131

array1=($(ls $CPV/*_1.fastq.gz))
for i in ${array1[@]}; do
	bbduk.sh in1=${i} out1=${i}.clean in2=$(echo ${i}|sed s/_1/_2/) out2=$(echo ${i}|sed s/_1/_2/).clean ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
	echo "adapter trimming ${i} $(date)"
	#Quality trimmming, of both the left and the right sides to get rid of reads that are less than quality 20
	bbduk.sh in1=${i}.clean out1=${i}.clean.trim in2=$(echo ${i}|sed s/_1/_2/).clean out2=$(echo ${i}|sed s/_1/_2/).clean.trim qtrim=rl trimq=20
	echo "quality trimming ${i} $(date)"
	#Quality filtering to get rid of entire low quality reads. maq=10 will trim reads that have average quality of less than 10
	bbduk.sh in1=${i}.clean.trim out1=${i}.clean.trim.filter in2=$(echo ${i}|sed s/_1/_2/).clean.trim out2=$(echo ${i}|sed s/_1/_2/).clean.trim.filter maq=10
	rm ${i}.clean
	rm $(echo ${i}|sed s/_1/_2/).clean
	rm ${i}.clean.trim
	rm $(echo ${i}|sed s/_1/_2/).clean.trim
	echo "quality filtering ${i} $(date)"
done

for i in ${array1[@]}; do
	bbduk.sh in1=${i}.clean.trim.filter in2=$(echo ${i}|sed s/_1/_2/).clean.trim.filter  bhist=${i}.b.hist qhist=${i}.q.hist gchist=${i}.gc.hist lhist=${i}.l.hist gcbins=auto
	echo ${i} > ${i}.hist.all
	echo "bhist" >> ${i}.hist.all
	cat ${i}.b.hist >> ${i}.hist.all
	echo "qhist" >> ${i}.hist.all
	cat ${i}.q.hist >> ${i}.hist.all
	echo "gchist" >> ${i}.hist.all
	cat ${i}.gc.hist >> ${i}.hist.all
	echo "lhist" >> ${i}.hist.all
	cat ${i}.l.hist >> ${i}.hist.all
	echo "histogram DONE $(date)"
	rm ${i}.*.hist
	gzip ${i}.clean.trim.filter
	gzip $(echo ${i}|sed s/_1/_2/).clean.trim.filter
done

echo "C_vir_Pro_RE22 trimming and filtering DONE $(date)"

##### BUILD HISAT INDEX #######

#DONE PREVIOUSLY

##### HISAT SAMTOOLS SORT #####
module purge
module load HISAT2/2.1.0-foss-2018b
module load SAMtools/1.9-foss-2018b

array2=($(ls $CPV/*_1.fastq.gz.clean.trim.filter.gz))
for i in ${array2[@]}; do
  # outputs a single bam file
	hisat2 --dta -x $CV/cvir_edited_index  -1 ${i} -2 $(echo ${i}|sed s/_1/_2/) -S ${i}.sam
	echo "HISAT2 PE ${i} $(date)"
  #SAMTOOLS sort to convert the SAM file into a BAM file to be used with StringTie. Stringtie only take sorted bam
  samtools sort ${i}.sam > ${i}.bam
  #Get bam file statistics for percentage aligned with flagstat
  samtools flagstat ${i}.bam > ${i}.bam.stats #get % mapped
 echo "${i} sorted bam done"
done

echo "C_vir_Pro_RE22 HISAT DONE $(date)"

##### Stringtie Assembly Quantify #####
module purge
module load StringTie/2.1.1-GCCcore-7.3.0 # new version of Stringtie
module load gffcompare/0.11.5-foss-2018b # new version of gffcompare
module load python/2.7.6

# assemble transcripts for each sample with the GFF3 annotation file
array3=($(ls $CPV/*.bam))
for i in ${array3[@]}; do
	stringtie -G $CV/ref_C_virginica-3.0_top_level.gff3 -o ${i}.gtf ${i}
	echo "${i} assembled"
	echo "${i}.gtf" >> $CPV/C_vir_Pro_RE22_mergelist.txt # Make stringtie mergelist with names of all gtf files with full path
done

#Run StringTie merge, merge transcripts from all samples in single experiment
stringtie --merge -G $CV/ref_C_virginica-3.0_top_level.gff3 -o $CPV/C_vir_Pro_RE22_stringtie_merged.gtf $CPV/C_vir_Pro_RE22_mergelist.txt
echo "C_vir_Pro_RE22 merged"

#gffcompare to compare how transcripts compare to reference annotation
gffcompare -r $CV/ref_C_virginica-3.0_top_level.gff3 -G -o $CPV/C_vir_Pro_RE22_stringtie_merged $CPV/C_vir_Pro_RE22_stringtie_merged.gtf
echo "C_vir_Pro_RE22 gffcompared"

#Re-estimate transcript abundance after merge step
for i in ${array3[@]}; do
		stringtie -A $(echo ${i}|sed "s/\..*//").abd.tab -e -G $CPV/C_vir_Pro_RE22_stringtie_merged.gtf -o $(echo ${i}|sed "s/\..*//").merge.gtf ${i}
		echo "${i} C_vir_Pro_RE22 transcript abundance re-estimated"
done

echo "C_vir_Pro_RE22 STRINGTIE DONE $(date)"

#### PREP STRINGTIE DESEQ2 #####
cd $CPV/
for i in *.merge.gtf; do
	# create text file with sample IDs and respective paths
	echo "$(echo $i |sed "s/\..*//") $CPV/$i" >> C_vir_Pro_RE22_sample_list.txt
done

python prepDE_Oct.2019.py -v -i C_vir_Pro_RE22_sample_list.txt -g Probiotic_RE22_gene_count_matrix.csv -t Probiotic_RE22_transcript_count_matrix.csv

echo "C_vir_Pro_RE22 DEseq 2 prep DONE $(date)"
echo "C_vir_Pro_RE22 FULL PIPELINE DONE $(date)"
