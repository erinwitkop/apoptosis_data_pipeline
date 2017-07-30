#!/bin/bash
#PBS-l nodes=3
#PBS-l walltime=1000:00:00
#PBS -j oe

#7_19_17

#Script to process reads with tophat2

Module load  bowtie2/2.2.9-foss-2016b
Module load TopHat/2.1.1-foss-2016b
Module load SAMtools/1.3.1-foss-2016b

F=/data3/marine_diseases_lab/erin/Bio_project_SRA/PE_fastq
S=/data3/marine_diseases_lab/erin/Bio_project_SRA
C=/data3/marine_diseases_lab/erin/Crassostrea_gigas_reference_genome/

 array1=($(ls $F/*_1.fq.clean.trim.filter))
 array2=($(ls $F/*_2.fq.clean.trim.filter))
 
 a

#first need to create an index from the reference genome using Bowtie2
# Before using a known annotation file with this option make sure that the 
# 1st column in the annotation file uses the exact same chromosome/contig names (case sensitive) 
# as shown by the bowtie-inspect command.

bowtie2-build /data3/marine_diseases_lab/erin/Crassostrea_gigas_reference_genome/Crassostrea_gigas.GCA_000297895.1.dna.nonchromosomal.fa Crassostrea_gigas_bowtie_index

bowtie-inspect --names Crassostrea_gigas_bowtie_index 


#Command to just prepare the transcriptome index files and then exit before running any reads
#Using this new usage no additional computational time is spent rebuilding the index in the future

tophat -G C/Crassostrea_gigas.GCA_000297895.1.36.gff3 --transcriptome-index C/C_gigas_GTF_index 


#Discussion on Tophat parameters:
	#mate inner distance: expected mean distance between mate pairs, (mean insert-size-2*read_length)
		#default good enough in Rondon et al. 2016
	#max intron length (Trapnell et al. 2013)
	#min intron length; (Rondon et al. 2016)

# Command  for running tophat with PE and SE reads

# Reminder: When running TopHat with paired reads it is critical that the *_1 files an the *_2 
# files appear in separate comma-delimited lists, and that the order of the files in the two lists is the same. 
# TopHat allows the use of additional unpaired reads to be provided after the paired reads. 
# These unpaired reads can be either given at the end of the paired read files on one side (as reads that can no longer be paired with reads from the other side), or they can be given in separate file(s) which are appended (comma delimited) to the list of paired input files on either side e.g.:

for i in ${array1[@]}; do
	tophat --max-intron-length 25000 --min-intron-length 50 --output-dir F --transcriptome-index C/C_gigas_GTF_index C/Crassostrea_gigas_bowtie_index ${i} $(echo ${i}|sed s/_1/_2/) 
done

#outputs two SAM format files, one for each of the strands. File suffix is .bam

#Discard unmapped reads with SAMtools (Rondon et al. 2016)

samtools view -f 4 file.bam > unmapped.sam



	




#References:
	#Rondon, R., F. Akcha, P. Alonso, D. Menard, J. Rouxel, C. Montagnani, G. Mitta, C. Cosseau, and C. Grunau. 2016. Transcriptional changes in Crassostrea gigas oyster spat following a parental exposure to the herbicide diuron. Aquat. Toxicol. 175:47–55. Available from: http://dx.doi.org/10.1016/j.aquatox.2016.03.007
	# http://wiki.bits.vib.be/index.php/Parameters_of_TopHat
	# Trapnell et al. 2013. Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks.
	# https://ccb.jhu.edu/software/tophat/manual.shtml#toph
	# https://ccb.jhu.edu/software/tophat/manual.shtml


