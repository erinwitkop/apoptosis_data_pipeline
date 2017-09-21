#!/bin/bash
#PBS-l nodes=1
#PBS-l walltime=1000:00:00
#PBS -j oe
#PBS -q default
#PBS -o out_MSTRG_isolate
#PBS -e err_MSTRG_isolate_err
#PBS -m ae -M erin_roberts@my.uri.edu


#isolate MSTRG IDs in the cluster into a smaller file that can be used by R locally 

cd /data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files/Bac_Viral_subset/

#Find line with Matching stringtie geneID found in the StringTie --merge output file
#Make file for all the MSTRGID_Gene_Sig.txt with all the significantly differentially expressed MSTRG transcripts

readarray a < /data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files/Bac_Viral_subset/MSTRGID_Gene_Sig.txt
for i in "${a[@]}"; do
	grep "$i"  Bac_viral_stringtie_merged.gtf >> Bac_viral_Gene_SIG_stringtie.merge.lines.txt
	echo "STOP" $(date)
done

#Make file for all the MSTRGID_gene_non_Sig.txt
readarray a < /data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files/Bac_Viral_subset/MSTRGID_gene_non_Sig.txt
for i in "${a[@]}"; do
	grep "$i"  Bac_viral_stringtie_merged.gtf >> Bac_viral_Gene_NON_SIG_stringtie.merge.lines.txt
	echo "STOP" $(date)
done

echo "DONE" $(date)