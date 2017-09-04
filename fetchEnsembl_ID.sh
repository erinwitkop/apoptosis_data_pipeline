#!/bin/bash
#PBS-l nodes=2
#PBS-l walltime=1000:00:00
#PBS -j oe

### Script Author: Erin M. Roberts. 2017. ###

#Script to match MSTRG values from DESeq2 output to their correct gene attributes from the Crassostrea gigas GFF3 file provided by Ensembl

set -e
echo "START" $(date)

#Set WD to working directory, and f to be the .csv file that was output from DESeq2 after performing differential sequence analysis
cd /data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files
f=resoshv1Tran_05_dfSig

#Retrieve first column (MSTRG gene ID) from the list of significantly differentially expressed genes
cut -d, -f1 $f.csv > $f_MSTRG.txt

#Remove quotes from around MSTRG IDs in the file
sed 's/\"//g' $f_MSTRG.txt > $f_MSTRG_edited.txt

#Access Stringtie merge file on the cluster, upload $f_MSTRG_edited.txt file to cluster
 
#Create stringtie.merge.$f_lines.txt file for read array grep output to go in

#Find line with Matching stringtie geneID found in the StringTie --merge output file
readarray a < /data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files/$f_MSTRG_edited.txt
for i in "${a[@]}"; do
	grep "$i" stringtie_merged.gtf >> stringtie.merge.$f_lines.txt
	echo "STOP" $(date)
done

#Extract first column of data from that line, input into new file
sed -e 's/\s.*$//' stringtie.merge.$f_lines.txt > $f_MSTRG_Accession.txt #sed used here removes everything but the first line

#Delete several extraneous "#" from beginning (if necessary)
tail -n +4 $f_MSTRG_Accession.txt > $f_MSTRG_Accession_edited.txt #customize the +4 to be + whatever number of lines from the front that need to be removed
C12748
#checked to make Accession number order of head of both files match
head $f_MSTRG_Accession_edited.txt
head stringtie.merge.$f_lines.txt

#Check the number of lines in each to make sure 
#cat stringtie.merge.$f_lines.txt | wc -l #671973313
#cat resoshv1_05_dfSig_MSTRG_Accession_edited.txt | wc -l # 671973310 (correct because 3 lines were deleted above)
#cat resoshv1_05_dfSig_matching_gff_lines.txt | wc -l #498076 
#cat resoshv1_05_dfSig_MSTRG_Accession_edited_unique.txt | wc -l #5433

#Find in Crassostrea gigas GFF with the same accession
	# Usage grep -F -f file1 file2
grep -F -f $f_MSTRG_Accession_edited.txt Crassostrea_gigas.gff > $f_matching_gff_lines.txt

#only print lines which are for genes or mRNA
#one combined file
awk '$3=="gene"' $f_matching.gff > $f_matching_genes_transcripts.gff
awk '$3=="mRNA"' $f_matching.gff >> $f_matching_genes_transcripts.gff

#one separated by type
awk '$3=="gene"' $f_matching.gff > $f_matching_genes.gff
awk '$3=="mRNA"' $f_matching.gff > $f_matching_transcripts.gff

#Use python to separate out the columns by tab and turn into csv, Get MSTRG IDs back onto the file after getting the GFF ID values
	#Run attributes_by_type.py in python

	

stringtie.merge.$f_lines.txt


#Load output into R


echo "STOP" $(date)


