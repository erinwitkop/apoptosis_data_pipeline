#!/bin/bash
#SBATCH -t 400:00:00
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH --exclusive
#SBATCH -o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/MUSCLE_RAxML_5_8_2020
#SBATCH -e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/MUSCLE_RAxML_error_5_8_2020

echo "START $(date)"

# Set paths needed
IAP=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/OrthoFinder_2020/OrthoFinder_Data_Analysis/Results_Mar25/IAP
M=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/MUSCLE_2020
R=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/RAxML_2020

# the two needed modules conflict with one another so I need to load MUSCLE first, then unload and load RAxML
module load MUSCLE/3.8.31-foss-2018b

# Generate alignments of all protein sequences using MUSCLE
muscle -in $IAP/IAP_all_orthogroups.fa -phyiout $M/IAP_all_orthogroups.phy -maxiters 2 -seqtype protein -sv

# -phyiout specifies output in Phylip interleaved format
# -maxiters 2 is suggested for large alignments (several thousand sequences), Running this option generally gives accuracy
    # comparable to TCoffee and speeds much faster than ClustalW
# -seqtype protein just explicitly specifies that the input type is protein
# -sv option is for usage with large number of sequences where memory may be an issue
# can also addi the diags1 option to help optimize for speed because of the very large alignment. Trying first without


echo "MUSCLE done $(date)"

# Unload RAxML and load RAxML
module purge
module load RAxML/8.2.10-goolf-2016b-mpi-avx

# Perform ML search and rapid bootstrapping with one command
mpiexec raxmlHPC-MPI-AVX -T 20 -s $M/IAP_all_orthogroups.phy -n $R/IAP_all_orthogroups_RAxML -m PROTGAMMAAUTO -x 12345 -p 12345 -f a -N autoMRE

# -s is the sequence file name
# -n is the outputFileName
# -m is the substitution model
# -x rapidBootstrapRandomNumberSeed
# -p parsimonyRandomSeed
# -T 20 uses 20 threads
# -f a is rapid Bootstrap analysis and search for best-scoring ML tree in one program run
# -N specifies the number of alternative runs on distinct starting trees to run.
    # -N: In combination with ­ combined with ­-f a -­x a rapid BS search and thereafter a thorough ML search on the original alignment.
    # -N: if you want to use the bootstrapping criteria rather than a set number, specifiy autoMRE with the -x or -b option
    # autoMRE applies the bookstrap convergence criterion

echo "RAxML done $(date)"
