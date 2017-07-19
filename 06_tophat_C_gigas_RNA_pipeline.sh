#7_19_17_notes

#Script to process reads with tophat2

Module load  bowtie2/2.2.9-foss-2016b
Module load TopHat/2.1.1-foss-2016b
Module load SAMtools/1.3.1-foss-2016b

F=/data3/marine_diseases_lab/erin/Bio_project_SRA/PE_fastq
C=/data3/marine_diseases_lab/erin/Crassostrea_gigas_reference_genome/

 array1=($(ls $F/*_1.fq.trim.filter))
 array2=($(ls $F/*_2.fq.trim.filter))

#first need to create an index from the reference genome using Bowtie2

bowtie2-build /data3/marine_diseases_lab/erin/Crassostrea_gigas_reference_genome/Crassostrea_gigas.GCA_000297895.1.dna.nonchromosomal.fa Crassostrea_gigas_bowtie_index

for i in ${array1[@]}; do
	tophat -${i}.trim.filter $(echo ${i}|sed s/_1/_2/).trim.filter C/Crassostrea_gigas_bowtie_index

#outputs two SAM format files, one for each of the strands