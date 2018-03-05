#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH --nodes 1
#SBATCH --mail-user=erin_roberts@my.uri.edu
#SBATCH -o /data3/marine_diseases_lab/erin/HHMER_analysis/HMM_output
#SBATCH -e /data3/marine_diseases_lab/erin/HHMER_analysis/HMM_error
#SBATCH -D /data3/marine_diseases_lab/erin/HHMER_analysis

#-D submits the start path
echo "START $(date)"

module load HMMER/3.1b2-foss-2016b
F=/data3/marine_diseases_lab/erin/HHMER_analysis

#Step 1: build a profile HMM with hmmbuild
#input file as Stockholm or FASTA alignments
#It expects Stockholm by default. To read aligned FASTA files, which HMMER calls “afa” format, 
#specify --informat afa on the command line of any program that reads an input alignment

hmmbuild --informat afa $F/BID.hmm $F/BID_CD_alignment.fasta

#Search sequence database with hmmsearch
#hmmsearch accepts any FASTA file as input. It also accepts EMBL/Uniprot text format. 
#It will automatically determine what format your file is in; you don’t have to say. 

hmmsearch $F/BID.hmm /data3/marine_diseases_lab/shared/GCF_002022765.2_C_virginica-3.0_rna.fasta > BID_CV_search.out

#INFO ON INTERPRETTING RESULT
#The second section is the sequence top hits list. It is a list of ranked top hits 
#(sorted by E-value, most significant hit first), formatted in a BLAST-like style
# Things to look at: Full sequence score sums up ALL domains (The full sequence score is summed over all possible alignments), the "best dom" score
#is only for the best scoring domain. 
	#if both E-values are significant (<< 1), the sequence is likely to be homologous to your query.
	#if the fullsequence E-value is significant but the single best domain E-value is not,
		#the target sequence is probably a multidomain remote homolog; but be wary, 
		#and watch out for the case where it’s just a repetitive sequence.
# Doms columns list first the expected # of domains, second column N list the number of domains
	#the software finally went with for annotating and aligning the target sequence
	
#INFO ON DOMAIN ANNOTATION FOR EACH SEQUENCE
#For each sequence in the top hits list, there will be a section containing a table of
#where HMMER3 thinks all the domains are, followed by the alignment inferred for each domain
# ! or ? symbol indicates whether this domain does or does not satisfy both per-sequence and per-domain inclusion thresholds.
#If the independent E-value is significant (<< 1), that means that even this single domain by itself is such a 
#strong hit that it suffices to identify the sequence as a significant homolog with respect to the size of the entire original database search.
#You can be confident that this is a homologous domain.
#Once there’s one or more high-scoring domains in the sequence already, 
#sufficient to decide that the sequence contains homologs of your query,
# you can look (with some caution) at the conditional E-value to decide 
# the statistical significance of additional weak-scoring domains.

echo "STOP $(date)" 

#File for reference on code http://eddylab.org/software/hmmer3/3.0rc2/Userguide.pdf