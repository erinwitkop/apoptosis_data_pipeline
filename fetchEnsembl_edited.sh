#!/bin/bash
#PBS-l nodes=2
#PBS-l walltime=1000:00:00
#PBS -j oe

#Script to match MSTRG values from DESeq2 output to their correct gene attributes from the Crassostrea gigas GFF3 file provided by Ensembl
set -e
echo "START" $(date)

#Set WD to working directory, and f to be the .csv file that was output from DESeq2 after performing differential sequence analysis
local=~/Documents/PhD_Research/GOMEZCHIARRI2_FILES_DAILYNOTES/DAILYNOTES_DATA/apoptosis_data_pipeline
cluster=/data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files
f=resoshv1Tran_05_dfSig

#Retrieve first column (MSTRG gene ID) from the list of significantly differentially expressed genes
cut -d, -f1 $f.csv > $/local/$f_MSTRG.txt

#Remove quotes from around MSTRG IDs in the file
sed 's/\"//g' $local/$f_MSTRG.txt > $local/$f_MSTRG_edited.txt

#Access Stringtie merge file on the cluster, upload $f_MSTRG.txt file to cluster

#Find line with Matching stringtie geneID found in the StringTie --merge output file
rm $cluster/stringtie.merge.$f_lines.txt
readarray a < $cluster/$f_MSTRG_edited.txt
for i in "${a[@]}"; do
	grep "$i" $cluster/stringtie_merged.gtf >> $cluster/stringtie.merge.$f_lines.txt
	echo "STOP" $(date)
done

#output
#Cnumbers=accession  exon/transc/etc phase start/end attributes(include mstring 
#need to keep col 1,2, parts of attributes

#Isolate the fields that you want, keeping the accession, type, and the MSTRG gene and MSTRG transcript ID
rm $cluster/stringtie.merge.matching.select.fields.txt
readarray a < $cluster/stringtie.merge.$f_lines.txt
#for each line
for line in "${a[@]}"; do
	#split into an array
	splitline=($line)
	#write to file the items in array that you want, taking account that each space is a separate item in the list
	echo "${splitline[0]}" "${splitline[2]}" "${splitline[8]}" "${splitline[9]}" "${splitline[10]}" "${splitline[11]}" >> stringtie.merge.matching.select.fields.txt
	echo "STOP" $(date)
done

#Find line in Crassostrea gigas GFF with the same accession for genes and transcripts, and keep the MSTRG ID with it
rm $cluster/MSTRG_merged_with_geneinfo_genes.txt
rm $cluster/MSTRG_merged_with_transcriptinfo_transcripts.txt
readarray a < $cluster/stringtie.merge.matching.select.fields.txt
for line in "${a[@]}"; do
		#split into an array
	splitline=($line)
	#write to file the items in array that you want
	gff=$(grep "${splitline[0]}" $cluster/Crassostrea_gigas.gff)
	#split gff
	splitgff=($gff)
	#if ${splitgff[2]} is 'gene'; then 
		echo "${line}" "${gff}" >>  MSTRG_merged_with_geneinfo_genes.txt 
	elif ${splitgff[2]} is 'transcript'; then 
		echo "${line}" "${gff}" >>  MSTRG_merged_with_transcriptinfo_transcripts.txt
	fi
done

echo "STOP" $(date)
