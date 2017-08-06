#!/bin/bash
#PBS-l nodes=3
#PBS-l walltime=1000:00:00
#PBS -j oe

set -e
echo "START" $(date)

# 8_5_17

module load HISAT2/2.0.4-foss-2016b   
module load SAMtools/1.3.1-foss-2016b

#HISAT2 code

#Indexing a reference genome
	#create new directory for the index called genome, and put the genome inside it

cd /data3/marine_diseases_lab/erin/Crassostrea_gigas_reference_genome/genome/

hisat2-build -f /data3/marine_diseases_lab/erin/Crassostrea_gigas_reference_genome/genome/Crassostrea_gigas_genome.fa  genome_index

# -f indicates that the reference input files are FASTA files

#Stay in the directory created in the previous step

#aligning paired end reads
P=/data3/marine_diseases_lab/erin/Bio_project_SRA/PE_fastq/PE_revised

for f in $P/*.fq.clean.trim.filter; do
	hisat2 -x /data3/marine_diseases_lab/erin/Crassostrea_gigas_reference_genome/genome/genome_index -1 ${f}_1.fq.clean.trim.filter -2 ${f}_2.fq.clean.trim.filter -S {f}.sam
done
 	#don't need -f because the reads are fastq

#Aligning single end reads
S=/data3/marine_diseases_lab/erin/Bio_project_SRA/SAR

hisat2 -x /data3/marine_diseases_lab/erin/Crassostrea_gigas_reference_genome/genome_index -U /data3/marine_diseases_lab/erin/Bio_project_SRA/reads_1.fa -S eg1.sam

	#put -x right before the genome index base, and the -U right before unpaired reads
	
This runs the HISAT2 aligner, which aligns a set of unpaired reads to the genome region using the index generated in the 


#SAMTOOLS view to convert the SAM file into a BAM file to be used with StringTie

