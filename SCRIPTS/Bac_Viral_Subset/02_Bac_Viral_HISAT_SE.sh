#!/bin/bash
#PBS-l nodes=2
#PBS-l walltime=1000:00:00
#PBS -j oe
#PBS -q default
#PBS -o out_Bac_Viral_HISAT
#PBS -e err_Bac_Viral_HISAT
#PBS -m ae -M erin_roberts@my.uri.edu

#02_Bac_Viral_HISAT_SE.sh, 09_07_17 Script to re-do the HISAT alignment steps with just the SRAs
#	from the OsHV-1 challenge and the Gram Negative and Gram positive challenges. Add in critical step
# 	to only use uniquely mapped reads. 

set -e
echo "START" $(date)

module load HISAT2/2.0.4-foss-2016b   
module load SAMtools/1.3.1-foss-2016b

#HISAT2 code
#Indexing a reference genome and no annotation file, make sure everything stays in the same directory
cd /data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files/Bac_Viral_subset
F=/data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files/Bac_Viral_subset

hisat2-build -f /data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files/Crassostrea_gigas_genome.fa  genome_index

# -f indicates that the reference input files are FASTA files

#Stay in the directory created in the previous step

#Aligning single end reads

array1=($(ls $F/*.filter))

for i in ${array1[@]}; do
        hisat2 --dta -x /data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files/genome_index -U ${i} -S ${i}.sam
        echo "${i}_DONE"
done
	
#This runs the HISAT2 aligner, which aligns a set of unpaired reads to the genome region using the index generated in the 

 	#don't need -f because the reads are fastq
	# put -x before the index
	# --dta : Report alignments tailored for transcript assemblers including StringTie.
	 #With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. 
	 #This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computation and memory usage.

#SAMTOOLS sort to convert the SAM file into a BAM file to be used with StringTie
array2=($(ls $F/*.sam))

for i in ${array2[@]}; do
	samtools sort -o ${i}.bam ${i}
	echo "${i}_convert"
done

#put -o before the out.bam and

#SAMTOOLS filter out for only uniquely mapped reads from bam results
#"samtools view -q 10 input > filtered_output

#FILTER OUT ANYTHING THAT DOES NOT HAVE A TMapq score of over 40, will give you reasonable stringency for 
#finding the best, most uniquely mapped reads

#reference: Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie, and Ballgown
#https://sequencing.qcfail.com/articles/mapq-values-are-really-useful-but-their-implementation-is-a-mess/

