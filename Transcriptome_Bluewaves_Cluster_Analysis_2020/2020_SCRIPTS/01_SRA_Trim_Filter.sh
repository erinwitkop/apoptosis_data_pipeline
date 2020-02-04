#!/bin/bash
#SBATCH -t 1000:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH	-o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/BBTools_qual_trim_out_Zhang_Vibrio_2_4_2020
#SBATCH	-e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/BBTools_qual_trim_error_Zhang_Vibrio_2_4_2020
#SBATCH	--mail-user=erin_roberts@my.uri.edu

echo "START $(date)"
module load BBMap/37.36-foss-2016b-Java-1.8.0_131

# Create variable for each path to raw data folder
CP=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/2020_Raw_Transcriptome_Data/C_vir_Probiotic_SRA
CR=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/2020_Raw_Transcriptome_Data/C_vir_ROD_SRA
CD=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/2020_Raw_Transcriptome_Data/Dermo_2015_transcriptomes/Fastq
GLO=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_deLorgeril_OsHV1_SRA
GHO=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_He_OsHV1_SRA
GRV=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA
GZV=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Zhang_Vibrio_SRA

##### Paired End Read Trimming + Filtering  ######
# PE experiments: C_gig_deLorgeril_OsHV1, C_gig_Rubio_Vibrio_SRA_ID, C_vir_Probiotic_SRA_ID, C_vir_Dermo (done in separate script)
# Create array variables for each set of PE experiment files

## C_gig_deLorgeril_OsHV1
array1=($(ls $GLO/*_1.fastq.gz))
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

	#ktrim = r means it will only trim from right side, which is where the adapter should be. (ktrim=l would trim from left)
	#hdist = hamming distance, hdist =1 allows for 1 mismatch
	#flag -tbo specifies to also trim adaptors based on pir overlap detection using BBMerge
	#which does not require known adapter sequences)
	#flag -tpe specified to trim both reads to the same length (if the adapter kmer was only detected in one of them and not other)

#Histogram generation, only generating for one of the pair (assuming that similar stats will be present).
#All histogram output contents are combined into one file
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
		#lhist = output a read length histogram
        #qhist = per base average quality
        #bhist = output a per-base composition histogram
        #gchist = output a gc content histogram

## C_gig_Rubio_Vibrio_SRA_ID: Trimming of adaptors
#array2=($(ls $GRV/*_1.fastq.gz))
#for i in ${array2[@]}; do
#	bbduk.sh in1=${i} out1=${i}.clean in2=$(echo ${i}|sed s/_1/_2/) out2=$(echo ${i}|sed s/_1/_2/).clean ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
#	echo "adapter trimming ${i} $(date)"
	#Quality trimmming, of both the left and the right sides to get rid of reads that are less than quality 20
#	bbduk.sh in1=${i}.clean out1=${i}.clean.trim in2=$(echo ${i}|sed s/_1/_2/).clean out2=$(echo ${i}|sed s/_1/_2/).clean.trim qtrim=rl trimq=20
#	echo "quality trimming ${i} $(date)"
	#Quality filtering to get rid of entire low quality reads. maq=10 will trim reads that have average quality of less than 10
#	bbduk.sh in1=${i}.clean.trim out1=${i}.clean.trim.filter in2=$(echo ${i}|sed s/_1/_2/).clean.trim out2=$(echo ${i}|sed s/_1/_2/).clean.trim.filter maq=10
#	rm ${i}.clean
#	rm $(echo ${i}|sed s/_1/_2/).clean
#	rm ${i}.clean.trim
#	rm $(echo ${i}|sed s/_1/_2/).clean.trim
#	echo "quality filtering ${i} $(date)"
#done

#for i in ${array2[@]}; do
#	bbduk.sh in1=${i}.clean.trim.filter in2=$(echo ${i}|sed s/_1/_2/).clean.trim.filter  bhist=${i}.b.hist qhist=${i}.q.hist gchist=${i}.gc.hist lhist=${i}.l.hist gcbins=auto
#	echo ${i} > ${i}.hist.all
#	echo "bhist" >> ${i}.hist.all
#	cat ${i}.b.hist >> ${i}.hist.all
#	echo "qhist" >> ${i}.hist.all
#	cat ${i}.q.hist >> ${i}.hist.all
#	echo "gchist" >> ${i}.hist.all
#	cat ${i}.gc.hist >> ${i}.hist.all
#	echo "lhist" >> ${i}.hist.all
#	cat ${i}.l.hist >> ${i}.hist.all
#	echo "histogram DONE $(date)"
#	rm ${i}.*.hist
#	gzip ${i}.clean.trim.filter
#	gzip $(echo ${i}|sed s/_1/_2/).clean.trim.filter
#done

## C_vir_Probiotic_SRA_ID: Trimming of adaptors
#array3=($(ls $CP/*_1.fastq.gz))
#for i in ${array3[@]}; do
#	bbduk.sh in1=${i} out1=${i}.clean in2=$(echo ${i}|sed s/_1/_2/) out2=$(echo ${i}|sed s/_1/_2/).clean ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
#	echo "adapter trimming ${i} $(date)"
	#Quality trimmming, of both the left and the right sides to get rid of reads that are less than quality 20
#	bbduk.sh in1=${i}.clean out1=${i}.clean.trim in2=$(echo ${i}|sed s/_1/_2/).clean out2=$(echo ${i}|sed s/_1/_2/).clean.trim qtrim=rl trimq=20
#	echo "quality trimming ${i} $(date)"
	#Quality filtering to get rid of entire low quality reads. maq=10 will trim reads that have average quality of less than 10
#	bbduk.sh in1=${i}.clean.trim out1=${i}.clean.trim.filter in2=$(echo ${i}|sed s/_1/_2/).clean.trim out2=$(echo ${i}|sed s/_1/_2/).clean.trim.filter maq=10
#	rm ${i}.clean
#	rm $(echo ${i}|sed s/_1/_2/).clean
#	rm ${i}.clean.trim
#	rm $(echo ${i}|sed s/_1/_2/).clean.trim
#	echo "quality filtering ${i} $(date)"
#done

#for i in ${array3[@]}; do
#	bbduk.sh in1=${i}.clean.trim.filter in2=$(echo ${i}|sed s/_1/_2/).clean.trim.filter  bhist=${i}.b.hist qhist=${i}.q.hist gchist=${i}.gc.hist lhist=${i}.l.hist gcbins=auto
#	echo ${i} > ${i}.hist.all
#	echo "bhist" >> ${i}.hist.all
#	cat ${i}.b.hist >> ${i}.hist.all
#	echo "qhist" >> ${i}.hist.all
#	cat ${i}.q.hist >> ${i}.hist.all
#	echo "gchist" >> ${i}.hist.all
#	cat ${i}.gc.hist >> ${i}.hist.all
#	echo "lhist" >> ${i}.hist.all
#	cat ${i}.l.hist >> ${i}.hist.all
#	echo "histogram DONE $(date)"
#	rm ${i}.*.hist
#	gzip ${i}.clean.trim.filter
#	gzip $(echo ${i}|sed s/_1/_2/).clean.trim.filter
#done

## C_vir_Dermo Ran in script by itself "01_SRA_Trim_Filter_Dermo_only.sh"

####### Single End Read Filtering + Trimming ######
# SE experiments: C_gig_He_2015_OsHV1_SRA_ID, C_gig_Zhang_Vibrio_SRA_ID, C_vir_ROD_SRA_ID

## C_gig_He_2015_OsHV1_SRA_ID
#array5=($(ls $GHO/*.fastq.gz))
#for i in ${array5[@]}; do
#	bbduk.sh in1=${i} out1=${i}.clean ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
#	echo "adapter trimming {i} $(date)"
	#Quality Trimming, of both the left and the right sides to get rid of reads that are less than quality 20
#	bbduk.sh in1=${i}.clean out1=${i}.clean.trim qtrim=rl trimq=20
#	echo "quality trimming ${i} $(date)"
 #Quality Filtering to get rid of entire low quality reads. maq=10 will trim reads that have average quality of less than 10
#	bbduk.sh in1=${i}.clean.trim out1=${i}.clean.trim.filter maq=10
#	rm ${i}.clean
#	rm ${i}.clean.trim
#	echo "quality filtering ${i} $(date)"
#done

#histogram generation
#for i in ${array5[@]}; do
#    bbduk.sh in1=${i}.clean.trim.filter bhist=${i}.b.hist qhist=${i}.q.hist gchist=${i}.gc.hist lhist=${i}.l.hist gcbins=auto
#    echo ${i} > ${i}.hist.all
#    echo "bhist" >> ${i}.hist.all
#    cat ${i}.b.hist >> ${i}.hist.all
#    echo "qhist" >> ${i}.hist.all
#    cat ${i}.q.hist >> ${i}.hist.all
#    echo "gchist" >> ${i}.hist.all
#    cat ${i}.gc.hist >> ${i}.hist.all
#    echo "lhist" >> ${i}.hist.all
#    cat ${i}.l.hist >> ${i}.hist.all
 #   rm ${i}.*.hist
 #	 echo "histogram ${i} $(date)"
#		gzip ${i}.clean.trim.filter
#done

## C_gig_Zhang_Vibrio_SRA_ID
array6=($(ls $GZV/*.fastq.gz))
for i in ${array6[@]}; do
	bbduk.sh in1=${i} out1=${i}.clean ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
	echo "adapter trimming {i} $(date)"
	#Quality Trimming, of both the left and the right sides to get rid of reads that are less than quality 20
	bbduk.sh in1=${i}.clean out1=${i}.clean.trim qtrim=rl trimq=20
  echo "quality trimming ${i} $(date)"
  Quality Filtering to get rid of entire low quality reads. maq=10 will trim reads that have average quality of less than 10
  bbduk.sh in1=${i}.clean.trim out1=${i}.clean.trim.filter maq=10
	rm ${i}.clean
	rm ${i}.clean.trim
  echo "quality filtering ${i} $(date)"
done

#histogram generation
for i in ${array6[@]}; do
    bbduk.sh in1=${i}.clean.trim.filter bhist=${i}.b.hist qhist=${i}.q.hist gchist=${i}.gc.hist lhist=${i}.l.hist gcbins=auto
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
		echo "histogram ${i} $(date)"
		gzip ${i}.clean.trim.filter
done

## C_vir_ROD_SRA_ID
array7=($(ls $CR/*.fastq.gz))
for i in ${array7[@]}; do
	bbduk.sh in1=${i} out1=${i}.clean ref=/opt/software/BBMap/37.36-foss-2016b-Java-1.8.0_131/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
	echo "adapter trimming {i} $(date)"
	#Quality Trimming, of both the left and the right sides to get rid of reads that are less than quality 20
	bbduk.sh in1=${i}.clean out1=${i}.clean.trim qtrim=rl trimq=20
	echo "quality trimming ${i} $(date)"
	#Quality Filtering to get rid of entire low quality reads. maq=10 will trim reads that have average quality of less than 10
	bbduk.sh in1=${i}.clean.trim out1=${i}.clean.trim.filter maq=10
	rm ${i}.clean
	rm ${i}.clean.trim
	echo "quality filtering ${i} $(date)"
done

#histogram generation
for i in ${array7[@]}; do
    bbduk.sh in1=${i}.clean.trim.filter bhist=${i}.b.hist qhist=${i}.q.hist gchist=${i}.gc.hist lhist=${i}.l.hist gcbins=auto
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
	  echo "histogram ${i} $(date)"
		gzip ${i}.clean.trim.filter
done

echo "full trim complete $(date)"
