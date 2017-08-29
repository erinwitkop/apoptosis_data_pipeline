#!/bin/bash
#PBS-l nodes=2
#PBS-l walltime=1000:00:00
#PBS -j oe

#fetchEnsembl_ID.sh, the un-generic version

set -e
echo "START" $(date)

#refer to 8_17_17_notes to see strategy
	#1. Gather all the gene MSTRG values, 
	#2. compare to the stringtie.merge file to find which GFF3 file annotation number they are, starting with >C
	#3. grep these are their corresponding sequences from the genome.fa file, or just use the identifier to get the Ensembl ID
	
	#4. Any MSTRGs that don't match, their sequences need to be gathered and using BLAST2GO they need to be BLASTed after retreiving the sequences from the 
		#original genome
	#5. Add a column with all the Ensembl IDs into the DESeq2 results output files for the significantly differentially
		#expressed ones so that I can do the functional enrichment


#Retrieve first column (MSTRG gene ID) from the list of significantly differentially expressed genes

#cd ~/Documents/PhD_Research/GOMEZCHIARRI2_FILES_DAILYNOTES/DAILYNOTES_DATA/apoptosis_data_pipeline/

#Set f to be the .csv file that was output from DESeq2 after performing differential sequence analysis
f=resoshv1_05_dfSig

cut -d, -f1 $f.csv > $f_MSTRG.txt

cat resoshv1_05_dfSig_MSTRG.txt

#Remove quotes
#sed 's/\"//g' resoshv1_05_dfSig_MSTRG.txt > resoshv1_05_dfSig_MSTRG_edited.txt

#Copy resoshv1_05_dfSig_MSTRG.txt to bluewaves
#scp resoshv1_05_dfSig_MSTRG_edited.txt bluewaves:/data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files

#Access Stringtie merge file on the cluster
#ssh bluewaves #enter password
#qsub -I
#cd /data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files
 
#Create resoshv1_05_dfSig_MSTRG_edited.txt file

#Find line with Matching stringtie geneID
cd /data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files
readarray a < /data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files/resoshv1_05_dfSig_MSTRG_edited.txt

for i in "${a[@]}"; do
	grep "$i" stringtie_merged.gtf >> stringtie.merge.resoshv1Sig_lines.txt
	echo "STOP" $(date)
done
 
cat stringtie.merge.resoshv1Sig_lines.txt

#Extract first column of data from that line, input into new file
#this sed script will remove everything but the first word
sed -e 's/\s.*$//' stringtie.merge.resoshv1Sig_lines.txt > resoshv1_05_dfSig_MSTRG_Accession.txt

#Deleted several extraneous "#" from beginning
tail -n +4 resoshv1_05_dfSig_MSTRG_Accession.txt > resoshv1_05_dfSig_MSTRG_Accession_edited.txt

#checked to make Accession number order of head of both files match
head resoshv1_05_dfSig_MSTRG_Accession_edited.txt
head stringtie.merge.resoshv1Sig_lines.txt

#Remove duplicate Accession entries
awk '!seen[$0]++' resoshv1_05_dfSig_MSTRG_Accession_edited.txt > resoshv1_05_dfSig_MSTRG_Accession_edited_unique.txt

#Check the number of lines
cat stringtie.merge.resoshv1Sig_lines.txt | wc -l #671973313
cat resoshv1_05_dfSig_MSTRG_Accession_edited.txt | wc -l # 671973310
cat resoshv1_05_dfSig_matching_gff_lines.txt | wc -l #498076 
cat resoshv1_05_dfSig_MSTRG_Accession_edited_unique.txt | wc -l #5433
#after removing duplicate Accessions 
cat resoshv1_05_dfSig_matching_gff_lines.txt | wc -l #498076

#Find line in GFF with the same accession
#grep -F -f file1 file2
grep -F -f resoshv1_05_dfSig_MSTRG_Accession_edited_unique.txt Crassostrea_gigas.gff > resoshv1_05_dfSig_matching_gff_lines.txt

#Extract attributes from GFF3 file using python after pulling out lines based on matching accessions using python parser
#gffutils (download via conda with $conda install --channel bioconda gffutils)

echo "STOP" $(date)


