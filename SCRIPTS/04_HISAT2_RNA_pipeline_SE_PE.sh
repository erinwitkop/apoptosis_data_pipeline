#!/bin/bash
#PBS-l nodes=2
#PBS-l walltime=1000:00:00
#PBS -j oe

set -e
echo "START" $(date)

# 8_5_17

module load HISAT2/2.0.4-foss-2016b   
module load SAMtools/1.3.1-foss-2016b

#HISAT2 code

#Indexing a reference genome and no annotation file
	#create new directory for the HISAT index called genome, and put the genome inside it
	# copy all reads files into this directory as well to ensure easy access by commands
cd /data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files
F=/data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files

hisat2-build -f /data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files/Crassostrea_gigas_genome.fa  genome_index

# -f indicates that the reference input files are FASTA files

#Stay in the directory created in the previous step

#Aligning paired end reads
array1=($(ls $F/*_1.fq))

for i in ${array1[@]}; do
	hisat2 --dta -x /data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files/genome_index  -1 ${i} -2 $(echo ${i}|sed s/_1/_2/) -S ${i}.sam
done
 	#don't need -f because the reads are fastq
	# put -x before the index
	# --dta : Report alignments tailored for transcript assemblers including StringTie.
	 #With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. 
	 #This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computation and memory usage.

#Aligning single end reads
array2=($(ls $F/*.fq))

for i in ${array2[@]}; do
        hisat2 --dta -x /data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files/genome_index -U ${i} -S ${i}.sam
        echo "${i}_DONE"
done
	
This runs the HISAT2 aligner, which aligns a set of unpaired reads to the genome region using the index generated in the 


#SAMTOOLS sort to convert the SAM file into a BAM file to be used with StringTie
#SAMTOOLS filter out low quality mapping results from bam file, and keep SAM header (while converting to bam)
array3=($(ls $F/*.sam))
	for i in ${array3[@]}; do
	 #This keeps the @SQ etc headers Stringtie requires
		samtools view -q 40 -b ${i} > ${i}.mapqfilter
		samtools sort ${i}.mapqfilter > ${i}.finalsorted #Stringtie takes as input only sorted bam files
		echo "${i}_filtered"
	done
	
	# -h to include header in the sam output
	# -q is to Skip alignments with MAPQ smaller than INT [0]
	# -b is to output in bam format, doing this first leaves the headers in the file, and sorts it by map quality
	#FILTER OUT ANYTHING THAT DOES NOT HAVE A TMapq score of over 40, will give you reasonable stringency for 
	#finding the best, most uniquely mapped reads


	
#reference: Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie, and Ballgown
