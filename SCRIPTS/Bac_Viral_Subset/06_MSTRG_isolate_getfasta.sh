#!/bin/bash
#PBS-l nodes=2
#PBS-l walltime=1000:00:00
#PBS -j oe
#PBS -q default
#PBS -o out_
#PBS -e err_
#PBS -m ae -M erin_roberts@my.uri.edu

set -e
echo "START" $(date)

module load bio/BEDTools/2.26.0-foss-2016b
module load fastx/0.0.14
module load genometools/1.5.9-foss-2016b 
module load pyfasta/0.5.2
F=/data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files/Bac_Viral_subset

#Using BEDtools getfasta in order to pull out MSTRG gene sequences after subsetting the stringtie.merge.gtf
#file for the sequences that were significantly differentially expressed

#Find line with Matching stringtie geneID found in the StringTie --merge output file
#Make file for all the MSTRGID_Gene_Sig.txt with all the significantly differentially expressed MSTRG transcripts

#cut -f2 OsHv1_MSTRGID_tran_Sig > OsHv1_MSTRGID_tran_Sig.txt
#cut -f2 Bac_MSTRGID_tran_Sig > Bac_MSTRGID_tran_Sig.txt

# USAGE: bash 06_MSTRG_isolate_getfasta.sh > OsHv1_MSTRGID_tran_SIG_merged.gtf

args=($(cat $F/OsHv1_MSTRGID_tran_non_Sig.txt))
while read -r line
do
    for i in ${args[@]}
    do
        case "$line" in
            *"$i"*) echo "$line";;
        esac
    done
done <$F/"Bac_viral_stringtie_merged.gtf"

#args=($(cat $F/Bac_MSTRGID_tran_Sig.txt))
#while read -r line
#do
#    for i in ${args[@]}
#    do
#        case "$line" in
#            *"$i"*) echo "$line";;
#        esac
#    done
#done <$F/"Bac_viral_stringtie_merged.gtf"


#isolate just the transcripts 
awk '$3=="transcript"' $F/OsHv1_MSTRGID_tran_SIG_merged.gtf > $F/OsHv1_MSTRGID_tran_SIG_merged_noexons.gtf

# awk '$3=="transcript"' $F/Bac_MSTRGID_tran_SIG_merged.gtf > $F/Bac_MSTRGID_tran_SIG_merged_noexons.gtf

#Use BEDTools to extract the sequences of MSTRG genes
#getfasta will extract the sequence in the orientation defined in the strand column when the “-s” option is used.
#tip: the headers in the input FASTA file must exactly match the chromosome column in the BED file.
#test this script on a subset first to see if it works properly...

#Make the "name" column the MSTRG transcript ID so the getfasta will retain it#
cut -f9 OsHv1_MSTRGID_tran_SIG_merged_noexons.gtf > OsHv1_MSTRGID_tran_SIG_merged_noexons_MSTRGID.gtf
#This will be merged in R! There is only one sequence per > from output below, so they are in the correct order

bedtools getfasta -name -s -fi $F/Crassostrea_gigas_genome.fa -bed $F/OsHv1_MSTRGID_tran_SIG_merged_noexons.gtf -fo $F/OsHV1_MSTRG_merged_TRAN_SEQUENCE.fa 

	#-name Use the “name” column in the BED file for the FASTA headers in the output FASTA file
	# -s Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented. 

#Will add the MSTRG ID to every line in the correct order back in my R script

#Remove duplicate sequences with genometools (if necessary) -Or can do it in R!
#gt sequniq -o OsHV1_MSTRG_merged_TRAN_SEQUENCE_unique.fa OsHV1_MSTRG_merged_TRAN_SEQUENCE.fa
gt splitfasta -numfiles 5 OsHV1_MSTRG_merged_TRAN_SEQUENCE.fa
gt splitfasta -numfiles 5 Bac_MSTRGID_tran_SIG_merged_SEQUENCE.fa

#IN ORDER TO GET FULL GENE LOCATIONS, YOU NEED TO HAVE GOTTEN THE FEATURE COUNTS SCRIPT OR GET THE -A GENE LOCATION OUTPUT

echo "DONE" $(date) 