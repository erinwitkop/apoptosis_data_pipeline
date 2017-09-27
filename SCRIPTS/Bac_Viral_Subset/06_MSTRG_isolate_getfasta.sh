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
F=/data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files/Bac_Viral_subset

#Using BEDtools getfasta in order to pull out MSTRG gene sequences after subsetting the stringtie.merge.gtf
#file for the sequences that were significantly differentially expressed

#Find line with Matching stringtie geneID found in the StringTie --merge output file
#Make file for all the MSTRGID_Gene_Sig.txt with all the significantly differentially expressed MSTRG transcripts

#cut -f2 OsHv1_MSTRGID_tran_Sig > OsHv1_MSTRGID_tran_Sig.txt
#cut -f2 Bac_MSTRGID_tran_Sig > Bac_MSTRGID_tran_Sig.txt

# USAGE: bash 06_MSTRG_isolate_getfasta.sh > OsHv1_MSTRGID_tran_SIG_merged.gtf

args=($(cat $F/OsHv1_MSTRGID_tran_Sig.txt))
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

bedtools getfasta -fi $F/Crassostrea_gigas_genome.fa -bed $F/OsHv1_MSTRGID_tran_SIG_merged_noexons.gtf -fo $F/OsHV1_MSTRG_merged_SEQUENCE.fa -s

#IN ORDER TO GET FULL GENE LOCATIONS, YOU NEED TO HAVE GOTTEN THE FEATURE COUNTS SCRIPT OR GET THE -A GENE LOCATION OUTPUT


