#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH --nodes 1
#SBATCH --mail-user=erin_roberts@uri.edu
#SBATCH -o /data3/marine_diseases_lab/shared/Cv_fetch_sequences_output
#SBATCH -e /data3/marine_diseases_lab/shared/Cv_fetch_sequences_error
#SBATCH -D /data3/marine_diseases_lab/

#-D submits the start path
echo "START $(date)"


F=/data3/marine_diseases_lab/shared/

array1=($(ls $F/CV_sig_apop_IDs_batch_entrez_feb_28_2018.txt)

for i in ${array1[@]}; do
	sed -n '/${i}/,/^>/p' GCF_002022765.2_C_virginica-3.0_rna.fasta | sed '$d' >> Cv_fasta_from_GENOME.fasta
done

echo "STOP $(date)"

