#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH --nodes 1
#SBATCH --mail-user=erin_roberts@uri.edu
#SBATCH -o /data3/marine_diseases_lab/erin/CV_apop_clustalw/tcoffe_output
#SBATCH -e /data3/marine_diseases_lab/erin/CV_apop_clustalw/tcoffee_error
#SBATCH -D /data3/marine_diseases_lab/erin/CV_apop_clustalw

#-D submits the start path
echo "START $(date)"

module load tcoffee/11.00.ddc7141-foss-2016b
F=/data3/marine_diseases_lab/erin/CV_apop_clustalw/

t_coffee $F/BID_conserved_plus_tophit.fa -mode accurate -outfile $F/BID_conserved_tophit_tcofee.out -email=erin_roberts@uri.edu

echo "STOP $(date)"
