#!/bin/bash
#PBS-l nodes=2
#PBS-l walltime=1000:00:00
#PBS -j oe

#Script to match MSTRG values from DESeq2 output to their correct gene attributes from the Crassostrea gigas GFF3 file provided by Ensembl

set -e
echo "START" $(date)

#Set WD to working directory, and f to be the .csv file that was output from DESeq2 after performing differential sequence analysis
WD=/data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files
f=resoshv1_05_dfSig

#Retrieve first column (MSTRG gene ID) from the list of significantly differentially expressed genes
cut -d, -f1 $f.csv > $f_MSTRG.txt

#Remove quotes from around MSTRG IDs in the file
sed 's/\"//g' $f_MSTRG.txt > $f_MSTRG_edited.txt

#Access Stringtie merge file on the cluster, upload $f_MSTRG.txt file to cluster
cd WD
 
#Create $f_MSTRG_edited.txt file

#Find line with Matching stringtie geneID found in the StringTie --merge output file
readarray a < /data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files/$f_MSTRG_edited.txt
for i in "${a[@]}"; do
	grep "$i" stringtie_merged.gtf >> stringtie.merge.$f_lines.txt
	echo "STOP" $(date)
done

#Extract first column of data from that line, input into new file
sed -e 's/\s.*$//' stringtie.merge.$f_lines.txt > $f_MSTRG_Accession.txt #sed used here removes everything but the first line

#Delete several extraneous "#" from beginning (if necessary_
#tail -n +4 $f_MSTRG_Accession.txt > $f_MSTRG_Accession_edited.txt #customize the +4 to be + whatever number of lines

#checked to make Accession number order of head of both files match
head $f_MSTRG_Accession_edited.txt
head stringtie.merge.$f_lines.txt

#Remove duplicate Accession entries
awk '!seen[$0]++' $f_MSTRG_Accession_edited.txt > $f_MSTRG_Accession_edited_unique.txt

#Check the number of lines in each to make sure 
#cat stringtie.merge.$f_lines.txt | wc -l #671973313
#cat resoshv1_05_dfSig_MSTRG_Accession_edited.txt | wc -l # 671973310 (correct because 3 lines were deleted above)
#cat resoshv1_05_dfSig_matching_gff_lines.txt | wc -l #498076 
#cat resoshv1_05_dfSig_MSTRG_Accession_edited_unique.txt | wc -l #5433
#after removing duplicate Accessions 
#cat resoshv1_05_dfSig_matching_gff_lines.txt | wc -l #498076 (matches same as above, indicating duplicates before this step don't need to be deleted)

#Find line in Crassostrea gigas GFF with the same accession
	# Usage grep -F -f file1 file2

grep -F -f $f_MSTRG_Accession_edited_unique.txt Crassostrea_gigas.gff > $f_matching_gff_lines.txt

#Extract attributes from GFF3 file using python script extract after pulling out lines based on matching accessions using python parser
python3 extract_attributes.py resoshv1_05_dfSig_matching_gff_lines.txt --print-records #must use the --print-records flag, and python3

resoshv1_05_dfSig_matching_gff_lines.txt


echo "STOP" $(date)


