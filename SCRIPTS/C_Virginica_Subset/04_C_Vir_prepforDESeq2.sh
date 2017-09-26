#!/bin/bash
#PBS-l nodes=2
#PBS-l walltime=1000:00:00
#PBS -j oe
#PBS -q default
#PBS -o out_C_vir_prep
#PBS -e err_C_vir_prep
#PBS -m ae -M erin_roberts@my.uri.edu

#04_C_Vir_prep.sh

cd /data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files/C_Vir_subset
F=/data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files/C_Vir_subset

array2=($(ls *.merged.gtf))

for i in ${array2[@]}; do
	echo "$(echo ${i}|sed "s/\..*//") $F/${i}" >> C_vir_sample_list.txt
done

python prepDE.py -i C_vir_sample_list.txt
			
echo "STOP" $(date)