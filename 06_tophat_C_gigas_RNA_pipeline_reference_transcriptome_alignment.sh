#!/bin/bash
#PBS-l nodes=3
#PBS-l walltime=1000:00:00
#PBS -j oe

#7_31_17

#Script to perform read alignments with tophat2 utilizing both bowtie index and GTF file. Discards unmapped reads
# and counts the percentage of mapped reads

module load  bowtie2/2.2.9-foss-2016b
module load TopHat/2.1.1-foss-2016b
module load SAMtools/1.3.1-foss-2016b

F=/data3/marine_diseases_lab/erin/Bio_project_SRA/PE_fastq
S=/data3/marine_diseases_lab/erin/Bio_project_SRA
C=/data3/marine_diseases_lab/erin/Crassostrea_gigas_reference_genome/
C_gigas_genome=$C/Crassostrea_gigas_genome.fa
GFF_file=$C/Crassostrea_gigas.GCA_000297895.1.36.gff3

 array1=($(ls $F/*_1.fq.clean.trim.filter))
 array2=($(ls $F/*_2.fq.clean.trim.filter))
 array3=($(ls $S/*.fastq.clean.trim.filter))
 

# First need to create an index from the reference genome using Bowtie2
# Please note that the values in the first column of the provided GTF/GFF file 
# (column which indicates the chromosome or contig on which the feature is located), 
# must match the name of the reference sequence in the Bowtie index you are using with TopHat. 
# Before using a known annotation file with this option make sure that the 
# check with bowtie-inspect command.
 
bowtie2-build $C/Crassostrea_gigas_genome.fa Crassostrea_gigas_bowtie_index
bowtie2-inspect --names $C/Crassostrea_gigas_bowtie_index
	# compare to the GFF3 file first column annotations
	#They don't match!
Bowtie2Index= $C/Crassostrea_gigas_bowtie_index
	
# Must rearrange the file entries to match, then can use the .GFF file
	#extract name from genome file index
	grep '^>' Crassostrea_gigas_genome.fa > Crassostrea_gigas_genome_headers
	sed 's/\s.*$//' Crassostrea_gigas_genome_headers > Crassostrea_gigas_genome_ID #to remove everything after first space
	sed 's/>//' Crassostrea_gigas_genome_ID > Crassostrea_gigas_genome_ID_string #to remove all the ">" characters
	
	#Reorder GFF 3 file by second column based on matching Crassostrea_gigas_genome_ID_string
	
	

#Prepare transcriptome index file and then exit before running any reads, saves time later


tophat --GTF $GFF_file --transcriptome-index=$C/C_gigas_transcriptome_index/transcriptome $Bowtie2Index
	# this will create a C_gigas_transcriptome folder in the current directory with files known known.gff, known.fa, known.fa.tlst, 
	# known.fa.ver and the known.* Bowtie index files

#Discussion on Tophat parameters:
	#mate inner distance: expected mean distance between mate pairs, (mean insert-size-2*read_length)
		#default good enough in Rondon et al. 2016
	#max intron length (Trapnell et al. 2013)
	#min intron length; (Rondon et al. 2016)

# Running tophat with PE and SE reads

# Reminder: When running TopHat with paired reads it is critical that the *_1 files an the *_2 
# files appear in separate comma-delimited lists, and that the order of the files in the two lists is the same. 
# TopHat allows the use of additional unpaired reads to be provided after the paired reads. 
# These unpaired reads can be either given at the end of the paired read files on one side (as reads that can no longer be paired with reads from the other side), or they can be given in separate file(s) which are appended (comma delimited) to the list of paired input files on either side e.g.:

#Because -G  was run once above to create own transcriptome index to map against for those reads that don't map to the Bowtie index

for i in ${array1[@]}; do
	tophat --max-intron-length 25000 --min-intron-length 50 -o $F --GTF $GFF_file $Bowtie2Index --transcriptome-index=$C/C_gigas_transcriptome_index/transcriptome ${i} $(echo ${i}|sed s/_1/_2/) 
	echo "STOP" $(date)
done

for i in ${array3[@]}; do
	tophat --max-intron-length 25000 --min-intron-length 50 -o $S $Bowtie2Index --transcriptome-index=$C/C_gigas_GTF_index ${i}
	echo "STOP" $(date)
done

#outputs in SAM format files

#Discard unmapped reads with SAMtools (Rondon et al. 2016)

for for i in ${array1[@]}; do 
	samtools view -f 4 ${i}.bam > ${i}.unmapped.sam
	samtools view -b -F 4 ${i}.bam > ${i}mapped.bam
	echo "STOP" $(date)
done
	# -F 4 excludes unmapped reads

for for i in ${array3[@]}; do 
	samtools view -f 4 ${i}.bam > ${i}.unmapped.sam
	samtools view -b -F 4 ${i}.bam > ${i}mapped.bam
	echo "STOP" $(date)
done
	
	#-F 4 excludes unmapped reads

#Count mapped reads percentage
	# look at alignment summary file in any output folder



#References:
	#Rondon, R., F. Akcha, P. Alonso, D. Menard, J. Rouxel, C. Montagnani, G. Mitta, C. Cosseau, and C. Grunau. 2016. Transcriptional changes in Crassostrea gigas oyster spat following a parental exposure to the herbicide diuron. Aquat. Toxicol. 175:47–55. Available from: http://dx.doi.org/10.1016/j.aquatox.2016.03.007
	# http://wiki.bits.vib.be/index.php/Parameters_of_TopHat
	# Trapnell et al. 2013. Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks.
	# https://ccb.jhu.edu/software/tophat/manual.shtml
	# https://ccb.jhu.edu/software/tophat/tutorial.shtml
	# http://khmer-protocols-reftrans.readthedocs.io/en/latest/refTrans/m-tophat.html


# Additional reasoning:
	# Supply TopHat with a set of gene model annotations and/or known transcripts, as a GTF 2.2 or GFF3 formatted file. 
	# If this option is provided, TopHat will first extract the transcript sequences and use Bowtie to align reads to this virtual transcriptome first.
	# Only the reads that do not fully map to the transcriptome will then be mapped on the genome. 
	# The reads that did map on the transcriptome will be converted to genomic mappings (spliced as needed) and
	# merged with the novel mappings and junctions in the final tophat output.
