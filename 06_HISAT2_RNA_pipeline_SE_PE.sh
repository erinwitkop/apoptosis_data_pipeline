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

#Indexing a reference genome and no annotation file
	#create new directory for the index called genome, and put the genome inside it

cd /data3/marine_diseases_lab/erin/Crassostrea_gigas_reference_genome/genome/

hisat2-build -f /data3/marine_diseases_lab/erin/Crassostrea_gigas_reference_genome/genome/Crassostrea_gigas_genome.fa  genome_index

# -f indicates that the reference input files are FASTA files

#Indexing a reference genome with an annotation file 

#Stay in the directory created in the previous step

#Aligning paired end reads
P=/data3/marine_diseases_lab/erin/Bio_project_SRA/PE_fastq/PE_revised

for f in $P/*.fq.clean.trim.filter; do
	hisat2 --dta -x /data3/marine_diseases_lab/erin/Crassostrea_gigas_reference_genome/genome/genome_index -1 ${f}_1.fq.clean.trim.filter -2 ${f}_2.fq.clean.trim.filter -S {f}.sam
done
 	#don't need -f because the reads are fastq
	# put -x before the index
	# --dta : Report alignments tailored for transcript assemblers including StringTie.
	 #With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. 
	 #This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computation and memory usage.
#Aligning single end reads
S=/data3/marine_diseases_lab/erin/Bio_project_SRA/SAR

hisat2 --dta -x /data3/marine_diseases_lab/erin/Crassostrea_gigas_reference_genome/genome_index -U /data3/marine_diseases_lab/erin/Bio_project_SRA/reads_1.fa -S ${f}.sam

	#put -x right before the genome index base, and the -U right before unpaired reads

	
This runs the HISAT2 aligner, which aligns a set of unpaired reads to the genome region using the index generated in the 


#SAMTOOLS sort to convert the SAM file into a BAM file to be used with StringTie

samtools sort -@ 8 -o ERR188044_chrX.bam ERR188044_chrX.sam $ samtools sort -@ 8 -o ERR188104_chrX.bam ERR188104_chrX.sam $ samtools sort -@ 8 -o ERR188234_chrX.bam ERR188234_chrX.sam $ samtools sort -@ 8 -o ERR188245_chrX.bam ERR188245_chrX.sam $ samtools sort -@ 8 -o ERR188257_chrX.bam ERR188257_chrX.sam $ samtools sort -@ 8 -o ERR188273_chrX.bam ERR188273_chrX.sam $ samtools sort -@ 8 -o ERR188337_chrX.bam ERR188337_chrX.sam $ samtools sort -@ 8 -o ERR188383_chrX.bam ERR188383_chrX.sam $ samtools sort -@ 8 -o ERR188401_chrX.bam ERR188401_chrX.sam $ samtools sort -@ 8 -o ERR188428_chrX.bam ERR188428_chrX.sam $ samtools sort -@ 8 -o ERR188454_chrX.bam ERR188454_chrX.sam $ samtools sort -@ 8 -o ERR204916_chrX.bam ERR204916_chrX.sam


