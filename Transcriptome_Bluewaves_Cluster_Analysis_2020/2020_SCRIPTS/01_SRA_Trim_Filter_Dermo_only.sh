#!/bin/bash
#SBATCH -t 1000:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH	-o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/BBTools_qual_trim_out_Dermo_2_3_2020
#SBATCH	-e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/BBTools_qual_trim_error_Dermo_2_3_2020
#SBATCH	--mail-user=erin_roberts@my.uri.edu

echo "START $(date)"
module load BBMap/37.36-foss-2016b-Java-1.8.0_131

# DERMO QUALITY TRIMMING TEST SCRIPT

CD=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/2020_Raw_Transcriptome_Data/Dermo_2015_transcriptomes/Fastq

## C_vir_Dermo
array4=($(ls $CD/*R1.fastq.gz))
for i in ${array4[@]}; do
		bbduk.sh in1=${i} out1=${i}.clean in2=$(echo ${i}|sed s/R1/R2/) out2=$(echo ${i}|sed s/R1/R2/).clean ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
		echo "adapter trimming ${i}" $(date)
		#Quality trimmming, of both the left and the right sides to get rid of reads that are less than quality 20
		bbduk.sh in1=${i}.clean out1=${i}.clean.trim in2=$(echo ${i}|sed s/R1/R2/).clean out2=$(echo ${i}|sed s/R1/R2/).clean.trim qtrim=rl trimq=20
		echo "quality trimming ${i}" $(date)
		#Quality filtering to get rid of entire low quality reads. maq=10 will trim reads that have average quality of less than 10
		bbduk.sh in1=${i}.clean.trim out1=${i}.clean.trim.filter in2=$(echo ${i}|sed s/R1/R2/).clean.trim out2=$(echo ${i}|sed s/R1/R2/).clean.trim.filter maq=10
		echo "quality filtering ${i}" $(date)
		rm ${i}.clean
		rm $(echo ${i}|sed s/R1/R2/).clean
		rm ${i}.clean.trim
		rm $(echo ${i}|sed s/R1/R2/).clean.trim
done

for i in ${array4[@]}; do
	  bbduk.sh in1=${i}.clean.trim.filter in2=$(echo ${i}|sed s/R1/R2/).clean.trim.filter  bhist=${i}.b.hist qhist=${i}.q.hist gchist=${i}.gc.hist lhist=${i}.l.hist gcbins=auto
		echo "STOP" $(date)
		echo ${i} > ${i}.hist.all
		echo "bhist" >> ${i}.hist.all
	  cat ${i}.b.hist >> ${i}.hist.all
		echo "qhist" >> ${i}.hist.all
		cat ${i}.q.hist >> ${i}.hist.all
		echo "gchist" >> ${i}.hist.all
		cat ${i}.gc.hist >> ${i}.hist.all
		echo "lhist" >> ${i}.hist.all
		cat ${i}.l.hist >> ${i}.hist.all
		rm ${i}.*.hist
		echo "histogram DONE" $(date)
		gzip ${i}.clean.trim.filter
		gzip $(echo ${i}|sed s/R1/R2/).clean.trim.filter
done
