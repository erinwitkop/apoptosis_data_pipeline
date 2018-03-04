#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH --nodes 1
#SBATCH --mail-user=erin_roberts@my.uri.edu
#SBATCH -o /data3/marine_diseases_lab/erin/CV_apop_clustalw/clustalw_output
#SBATCH -e /data3/marine_diseases_lab/erin/CV_apop_clustalw/clustalw_error
#SBATCH -D /data3/marine_diseases_lab/erin/CV_apop_clustalw

#-D submits the start path
echo "START $(date)"

module load ClustalW2/2.1-foss-2016b
F=/data3/marine_diseases_lab/erin/CV_apop_clustalw/

clustalw2 -infile=$F/BID_conserved_plus_tophit.fa -type=protein -gapopen=10.00 -gapext=0.10 -pwmatrix=GONNET -pwdnamatrix=IUB -outfile=$F/BID_conserved_tophit_clustalw.out -outorder=input


echo "STOP $(date)"
