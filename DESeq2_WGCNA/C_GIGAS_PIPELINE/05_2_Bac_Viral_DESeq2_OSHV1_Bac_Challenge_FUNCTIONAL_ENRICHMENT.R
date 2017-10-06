#05_2_Bac_Viral_DESeq2_OSHV1_Bac_Challenge_FUNCTIONAL_ENRICHMENT.R

#this script utilized topGO to perform Gene Set Ennrichment Analysis (GSEA)

#Load packages
source("http://bioconductor.org/biocLite.R")
biocLite("topGO")
library(topGO)
library(tidyr)
library(dplyr)
library(plyr)
library(genefilter)

####BLAST2GO the MSTRG SIG genes (AND MSTRG NON SIG GENES), then perform gene set enrichment on these and pull out more specific apoptosis genes ####
##subset out only the genes that have "MSTRG ID"
#lookup these lines in the stringtie.merge file to find the sequence for them (SEE 06_MSTRG_isolate_getfasta.sh), then put them into BLAST2GO
#Sig MSTRG's from OsHV1
resoshv1Tran_05_dfSig_Transcript_MSTRG <- resoshv1Tran_05_dfSig[grep("MSTRG", rownames(resoshv1Tran_05_dfSig)), ] #for Significant genes
OsHv1_MSTRGID_tran_Sig <- as.data.frame(rownames(resoshv1Tran_05_dfSig_Transcript_MSTRG))
head(OsHv1_MSTRGID_tran_Sig)
write.table(OsHv1_MSTRGID_tran_Sig[,1], file="OsHv1_MSTRGID_tran_Sig", sep= "\t")

#NON SI MSTRGs OsHV1
resoshv1Tran_05_df_non_Sig_Transcript_MSTRG <- resoshv1Tran_05_df_non_Sig[grep("MSTRG", rownames(resoshv1Tran_05_df_non_Sig)),]
OsHV1_MSTRGID_tran_non_Sig <- as.data.frame(rownames(resoshv1Tran_05_df_non_Sig_Transcript_MSTRG))
head(OsHV1_MSTRGID_tran_non_Sig)
write.table(OsHV1_MSTRGID_tran_non_Sig[,1], file="OsHv1_MSTRGID_tran_non_Sig", sep= "\t")

#Sig MSTRGs from Bacterial Challenge
resBacTran_05_dfSig_Transcript_MSTRG <- resBacTran_05_dfSig[grep("MSTRG", rownames(resBacTran_05_dfSig)), ]
Bac_MSTRGID_tran_Sig <- as.data.frame(rownames(resoshv1Tran_05_dfSig_Transcript_MSTRG))
head(Bac_MSTRGID_tran_Sig)
write.table(Bac_MSTRGID_tran_Sig[,1], file="Bac_MSTRGID_tran_Sig", sep= "\t")

#NON SIG MSTRGs from Bacterial Challenge
resBacTran_05_df_non_Sig_MSTRG <- resBacTran_05_df_non_Sig[grep("MSTRG", rownames(resBacTran_05_df_non_Sig)),]
Bac_MSTRG_tran_non_Sig <- as.data.frame(rownames(resBacTran_05_df_non_Sig_MSTRG))
head(Bac_MSTRG_tran_non_Sig)
write.table(write.table(Bac_MSTRG_tran_non_Sig[,1], file="Bac_MSTRGID_tran_non_Sig", sep= "\t"))

#USE 06_MSTRG_isolate_getfasta.sh script to get the sequence of genes while using the bam files of each and the sequence
#INSERT SIG SEQUENCES SEQUENCES INTO BLAST2GO

####If only using annotated genes!#####


####Data preparation####
#input is a table that contains the genes, the GO mappings, and the p-values
#OsHV1
head(resoshv1Tran_05_df)
OsHV1_topGOsubset <- 
  resoshv1Tran_05_df[grep("transcript:", rownames(resoshv1Tran_05_df)), ]
head(OsHV1_topGOsubset)
nrow(OsHV1_topGOsubset) #rows = 12112
#make a new column with the rownames, called Ensembl_genomes_transcript
OsHV1_topGOsubset$Ensembl_Genomes_Transcript<-rownames(OsHV1_topGOsubset)  
##needs to be made a data frame so that can be joined later by plyr
OsHV1_topGOsubset_df <- as.data.frame(OsHV1_topGOsubset)  
OsHV1_topGOtranscript_subset <- as.data.frame(rownames(OsHV1_topGOsubset))
head(OsHV1_topGOtranscript_subset)
colnames(OsHV1_topGOtranscript_subset) <- "transcript"
OsHV1transcriptIDdf <- OsHV1_topGOtranscript_subset %>%
  separate(transcript, c("transcript", "EKC"), ":")
#write to string
OsHV1_topGO_transcriptIDstring <- toString(OsHV1transcriptIDdf[,2], sep=',')
write(OsHV1_topGO_transcriptIDstring, "OsHV1_topGO_transcriptIDstring", sep = ",")
#cut off commas with sed 's/,//g' INFILE > outfile

#Load annotation file
OsHV1_topGO_annotations <- read.csv("OsHV1_topGO_transcriptIDstring2.csv", sep=",", header=TRUE) #12092 rows
nrow(OsHV1_topGO_annotations)
head(OsHV1_topGO_annotations)
#add "transcript: before every string" in Ensembl Genomes Transcript
OsHV1_topGO_annotations$Ensembl_Genomes_Transcript = paste('transcript', OsHV1_topGO_annotations$Ensembl_Genomes_Transcript, sep=':')
#concatentate gene expression and uniprot table
#Perform join by plyr
OsHV1_topGO_annotations
OsHV1_topGO_annotations_FULL <- join(OsHV1_topGO_annotations, OsHV1_topGOsubset_df, by="Ensembl_Genomes_Transcript")
nrow(OsHV1_topGO_annotations_FULL) #120921 This is what I wanted 

#Bac
head(resBacTran_05_df)
BAC_topGOsubset <- 
  resBacTran_05_df[grep("transcript:", rownames(resBacTran_05_df)), ]
head(BAC_topGOsubset)
#make a new column with the rownames, called Ensembl_genomes_transcript
BAC_topGOsubset$Ensembl_Genomes_Transcript<-rownames(BAC_topGOsubset)  
##needs to be made a data frame so that can be joined later by plyr
BAC_topGOsubset_df <- as.data.frame(BAC_topGOsubset)  #use this for the plyr join!
BAC_topGOtranscript_subset <- as.data.frame(rownames(BAC_topGOsubset))
head(BAC_topGOtranscript_subset)
colnames(BAC_topGOtranscript_subset) <- "transcript"
BACtranscriptIDdf <- BAC_topGOtranscript_subset %>%
  separate(transcript, c("transcript", "EKC"), ":") #edited to make rownames
nrow(BACtranscriptIDdf) #8075
#write to string so I can do lookup
BAC_topGO_transcriptIDstring <- toString(BACtranscriptIDdf[,2], sep=',')
write(BAC_topGO_transcriptIDstring, "BAC_topGO_transcriptIDstring", sep = ",")

#Load annotation file
BAC_topGO_annotations <- read.csv("BAC_topGO_transcriptIDstring.csv", sep=",", header=TRUE)
nrow(BAC_topGO_annotations) #8092
head(BAC_topGO_annotations)
#add "transcript: before every string" in Ensembl Genomes Transcript
BAC_topGO_annotations$Ensembl_Genomes_Transcript = paste('transcript', BAC_topGO_annotations$Ensembl_Genomes_Transcript, sep=':')
#concatentate gene expression and uniprot table
#Perform join by plyr
BAC_topGO_annotations_FULL <- join(BAC_topGO_annotations, BAC_topGOsubset_df, by="Ensembl_Genomes_Transcript")
nrow(BAC_topGO_annotations_FULL) #8092


#####Prepare topGO data ####
#1. filter out low count data

#TOPGO: Running the Enrichment tests




REVIGO.R to get functional enrichment



# Cytoscape word cloud ?

#References: Gene set enrichment analysis with topGO
#Adrian Alexa, Jorg Rahnenfuhrer, April 24, 2017, http://www.mpi-sb.mpg.de/alexa
