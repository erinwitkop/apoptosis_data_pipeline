#!/bin/bash
#PBS-l nodes=1
#PBS-l walltime=1000:00:00
#PBS -j oe

#7_31_17

#Script to perform read alignments with tophat2 utilizing both bowtie index and GTF file. Discards unmapped reads
# and counts the percentage of mapped reads
set -e 

echo "START" $(date)

module load Bowtie2/2.2.9-foss-2016b 
module load TopHat/2.1.1-foss-2016b
module load SAMtools/1.3.1-foss-2016b

F=/data3/marine_diseases_lab/erin/Bio_project_SRA/PE_fastq
S=/data3/marine_diseases_lab/erin/Bio_project_SRA
C=/data3/marine_diseases_lab/erin/Crassostrea_gigas_reference_genome
C_gigas_genome=$C/Crassostrea_gigas_genome.fa
GFF_file=$C/Crassostrea_gigas.GCA_000297895.1.36.gff3

# First need to create an index from the reference genome using Bowtie2
# Please note that the values in the first column of the provided GTF/GFF file 
# (column which indicates the chromosome or contig on which the feature is located), 
# must match the name of the reference sequence in the Bowtie index you are using with TopHat. 
#Before using a known annotation file with this option make sure that the 
# check with bowtie-inspect command.
 
bowtie2-build $C/Crassostrea_gigas_genome.fa $C/Crassostrea_gigas_bowtie_index

#bowtie-inspect --names $C/Crassostrea_gigas_bowtie_index
	# compare to the GFF3 file first column annotations
#Bowtie2Index= $C/Crassostrea_gigas_bowtie_index

echo "STOP" $(date) 
