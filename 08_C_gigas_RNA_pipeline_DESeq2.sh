#!/bin/bash
#PBS-l nodes=3
#PBS-l walltime=1000:00:00
#PBS -j oe

#8_2_17

#Script to take Salmon output of read abundances and pipe those into DESeq2 for normalization and Differential Expression Analysis