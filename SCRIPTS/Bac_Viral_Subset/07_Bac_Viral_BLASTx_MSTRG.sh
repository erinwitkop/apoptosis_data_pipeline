#!/bin/bash
#PBS-l nodes=1
#PBS-l walltime=1000:00:00
#PBS -j oe
#PBS -q default
#PBS -o out_blastx
#PBS -e err_blastx
#PBS -m ae -M erin_roberts@my.uri.edu

#07_Bac_Viral_BLASTx_MSTRG.sh

#Perform BLASTx search of the sequences against the nr database

set -e 
echo "START $(date)"
cd /data3/marine_diseases_lab/erin/Bio_project_SRA/pipeline_files/Bac_Viral_subset

module load BLAST+/2.5.0-foss-2016b-Python-2.7.12

blastx -db nr -query OsHV1_MSTRG_merged_TRAN_SEQUENCE.fa -out OsHV1_MSTRG_merged_TRAN_SEQUENCE_BLASTx.txt -evalue 0.0001 -remote

echo "DONE" $(date)

#evalue of 10^-3 cited in Pales, E., E. Corre, and B. Allam. 2014. Pallial mucus of the oyster Crassostrea virginica regulates the expression of putative virulence genes of its pathogen Perkinsus marinus q. Int. J. Parasitol. 44:305–317. Available from: http://dx.doi.org/10.1016/j.ijpara.2014.01.006

