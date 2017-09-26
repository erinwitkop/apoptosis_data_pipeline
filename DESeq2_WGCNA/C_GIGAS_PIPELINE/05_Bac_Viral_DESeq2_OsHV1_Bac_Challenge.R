#05_Bac_viral_DESeq2_OsHV1_Bac_challenge

####INPUT DATA GATHERED FROM CRASSOSTREA GIGAS####

#This script takes as input the output Bac_Viral_transcript_count_matrix.csv data prepared from prepDE.py and performs
#differential transcript expression analysis, and subsets out isoforms of GIMAPs and CgIAPs and graphs their
#relative abundance.

#call the DESeq2 library 
source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade") 
biocLite("DESeq2")
library("DESeq2")
#install.packages("fdrtool")
library(fdrtool)
#source("https://bioconductor.org/biocLite.R")
#biocLite("genefilter")
library(genefilter)
install.packages("dplyr")
library(dplyr)
install.packages("tidyr")
library(tidyr)
install.packages("reshape")
library(reshape2)
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
library(ggplot2)
# Construct Bac_Viral_PHENO_DATA.csv that contains SRA run information, such as which contrast, tissue, etc.

####DEG Analysis with TRANSCRIPT Count Matrix ####
#load transcript count matrix and labels
#Bac_Viral_PHENO_DATA.csv file contains metadata on the count table's samples
###Make sure excel PHENODATA is in the same order or these commands will change data to be wrong!!!!###

TranscriptCountData <- as.matrix(read.csv("Bac_Viral_transcript_count_matrix.csv", row.names="transcript_id"))
head(TranscriptCountData)

####Subset Data for OsHV1 ####
#Extract columns you want from the TranscriptCountData, based on which column the correct SRA data is in for the
#Extract columns 1:30 from the FullTranscriptCountData, these are the SRA's from the OsHV-1 experiment
oshv1TranCountData <- as.matrix(TranscriptCountData[ , c(1:30)])
#head(oshv1TranCountData)
TranColData <- read.csv("OsHV1_PHENO_DATA.csv", header=TRUE, sep=",")
oshv1TranColData <- TranColData[, c("sampleID", "condition", "time.h.")]
print(oshv1TranColData)
rownames(oshv1TranColData) <- oshv1TranColData$sampleID
colnames(oshv1TranCountData) <- oshv1TranColData$sampleID
head(oshv1TranCountData)
head(oshv1TranColData)
#Give the stressorLevel column levels
oshv1TranColData$time.h. <- factor(oshv1TranColData$time.h.)
levels(oshv1TranColData$time.h.) #check to see that it has levels 

# Check all sample IDs in oshv1ColData are also in oshv1CountData and match their orders
all(rownames(oshv1TranColData) %in% colnames(oshv1TranCountData))  #Should return TRUE
# returns TRUE
all(rownames(oshv1TranColData) == colnames(oshv1TranCountData))    # should return TRUE
#returns TRUE

#### Create OsHV-1 DESeq Data Set from Matrix ####
# DESeqDataSet from count matrix and labels 
#design purposefully doesn't account for time, not looking at time interaction, just condition
ddsOshv1Tran <- DESeqDataSetFromMatrix(countData = oshv1TranCountData, 
                                       colData = oshv1TranColData, 
                                       design = ~ condition)
#
ddsOshv1Tran <- ddsOshv1Tran[ rowSums(counts(ddsOshv1Tran)) > 1, ]

#if design was design = ~stressorLevel + condition, #this design will gather the effect of condition, accounting for the sample pairing by time
# review how the data set looks
head(ddsOshv1Tran)

#Relevel each to make sure that control is the first level in the treatment factor for each
ddsOshv1Tran$condition <- relevel( ddsOshv1Tran$condition, "control")

#Check we're looking at the right samples
as.data.frame( colData(ddsOshv1Tran) )

#Running the DEG pipeline
ddsOshv1Tran<- DESeq(ddsOshv1Tran) #for designs with interactions, recommends setting betaPrior=FALSE

#Inspect results table
#extract contrasts between control and treatment values
resoshv1Tran<- results(ddsOshv1Tran, contrast = c("condition", "control", "treatment"))
#to extract log2fold change and p values under 0.1 and 0.05
head(resoshv1Tran)

#Order by Log2FC
head( resoshv1Tran[ order( resoshv1Tran$log2FoldChange ), ] ) #head for strongest downregulation
tail( resoshv1Tran[ order( resoshv1Tran$log2FoldChange ), ] ) #tail for strongest up regulation
summary(resoshv1Tran)
sum(resoshv1Tran$padj < 0.1, na.rm=TRUE) #4410
#Change alpha setting in DESeq Results
resoshv1Tran_05 <- results(ddsOshv1Tran,alpha=0.05)
summary(resoshv1Tran_05)
sum(resoshv1Tran_05$padj < 0.05, na.rm=TRUE) #3298

#metadata on meaning of the columns
mcols(resoshv1Tran_05, use.names = TRUE)
#Get more detailed description
mcols(resoshv1Tran_05_df)$description
#shows treatment vs. control with baseMean as the mean of normalized counts for all samples

####p-value correction for all genes in resoshv1Tran ####
#First Visualize histograms
#histogram of P- values to visualize any "hills" or "U shape"
# hill means variance of the null distribution too high, U shape means variance assumed too low
hist(resoshv1Tran_05$pvalue, breaks = 20, col = "grey") #hill

#remove filtered out genes by independent filtering, they have NA adj. pvals
resoshv1Tran_05_df <- resoshv1Tran[ !is.na(resoshv1Tran_05$padj), ]

#remove genes with NA pvals (outliers)
resoshv1Tran_05_df <- resoshv1Tran_05_df[ !is.na(resoshv1Tran_05_df$pvalue), ]

#remove adjsuted pvalues, since we add the fdrtool results later on (based on the correct p-values)
resoshv1Tran_05_df <- resoshv1Tran_05_df[, -which(names(resoshv1Tran_05_df) == "padj")]

#use z-scores as input to FDRtool to re-estimate the p-value
FDR.resoshv1Tran_05_df <- fdrtool(resoshv1Tran_05_df$stat, statistic= "normal", plot = T)

#add values to the results data frame, also ad new BH- adjusted p-values
resoshv1Tran_05_df[,"padj"] <- p.adjust(FDR.resoshv1Tran_05_df$pval, method = "BH")

#replot corrected p-values 
hist(FDR.resoshv1Tran_05_df$pval, col = "royalblue4",
     main = "Correct null model OsHv1 Transcript Count", xlab = "CORRECTED p-values")

#Check how many genes have BH adjusted p values of less than 0.05 after P-value correction?
sum( resoshv1Tran_05_df$padj < 0.05, na.rm=TRUE ) #1459

#Subset the results table to the differentially expressed genes under FDR 0.1, order the Log2FC table first by strongest down regulation
resoshv1Tran_05_dfSig <- resoshv1Tran_05_df[ which(resoshv1Tran_05_df$padj < 0.05 ), ]
head( resoshv1Tran_05_dfSig[ order( resoshv1Tran_05_dfSig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resoshv1Tran_05_dfSig[ order( resoshv1Tran_05_dfSig$log2FoldChange ), ] ) #tail for strongest up regulation
summary(resoshv1Tran_05_dfSig)
resoshv1Tran_05_df_non_Sig <- resoshv1Tran_05_df[ which(resoshv1Tran_05_df$padj > 0.05 ), ]
summary(resoshv1Tran_05_df_non_Sig)

#Visualize Results with Diagnostic Plots#
#MA plot, useful overview for experiment with two-group comparison. Plots log2FC over mean of normalized counts
#genes with adjusted p value 
plotMA(resoshv1Tran_05_dfSig)
plotMA(resoshv1Tran_05_df_non_Sig)

#Export Results to CSV
write.csv( as.data.frame(resoshv1Tran_05_df), file="OsHV1_resoshv1Tran_05_df.csv")
write.csv( as.data.frame(resoshv1Tran_05_dfSig), file="OsHV1_resoshv1Tran_05_dfSig.csv")
write.csv( as.data.frame(resoshv1Tran_05_df_non_Sig), file = "OsHV1_resoshv1Tran_05_df_non_Sig.csv")


####Subset files for only those that have Transcript IDs####
#Extract gene titles of all the significantly differentially expressed genes
OsHV1_withID_subset_resoshv1Tran_05_dfSig <- 
  resoshv1Tran_05_dfSig[grep("transcript:", rownames(resoshv1Tran_05_dfSig)), ]
head(OsHV1_withID_subset_resoshv1Tran_05_dfSig)
transcriptIDdf <- as.data.frame(rownames(OsHV1_withID_subset_resoshv1Tran_05_dfSig))
head(transcriptIDdf)
transcriptID1thru5 <- rownames(head(OsHV1_withID_subset_resoshv1Tran_05_dfSig[1:5,]))
transcriptIDdf= transform(transcriptIDdf, 
        ID = colsplit(rownames(OsHV1_withID_subset_resoshv1Tran_05_dfSig), 
                      split = "\\:", names = c('transcript:', 'EKC33371')))
  #this line splitsup the word transcript and the transcript name

transcriptIDstring <- toString(transcriptIDdf[,3], sep=',')
transcriptIDstring
write(transcriptIDstring, "transcriptIDstring", sep = ",")
#write this to a file and then perform look up on the UniProt website

#Extract gene titles from non Sig genes
OsHV1_withID_subset_resoshv1Tran_05_df_non_Sig <- 
  resoshv1Tran_05_df_non_Sig[grep("transcript:", rownames(resoshv1Tran_05_df_non_Sig)), ]
head(OsHV1_withID_subset_resoshv1Tran_05_df_non_Sig)
transcriptIDdf_nonSig <- as.data.frame(rownames(OsHV1_withID_subset_resoshv1Tran_05_df_non_Sig))
head(transcriptIDdf_nonSig)
transcriptIDdf_nonSig= transform(transcriptIDdf_nonSig, 
                          ID = colsplit(rownames(OsHV1_withID_subset_resoshv1Tran_05_df_non_Sig), 
                                        split = "\\:", names = c('transcript:', 'EKC37466')))
transcriptIDstring_nonSig <- toString(transcriptIDdf_nonSig[,3], sep=',')
transcriptIDstring_nonSig
write(transcriptIDstring_nonSig, "transcriptIDstring_nonSig", sep = ",")
#write this to a file and then perform look up on the UniProt website, using Ensembl Genomes Transcript as look up

####Upload Sig Differentially Expressed transcripts from the UniProt.ws website ####
oshv1_transcriptIDs_UniProt_SIG <- read.csv("OsHV1_resoshv1Tran_dfSig_transcriptIDstring.csv", header=TRUE)
#make sure to upload as csv version so that the protein names load correctly
head(oshv1_transcriptIDs_UniProt_SIG)
oshv1_transcriptIDs_UniProt_SIG["Challenge"] <- "OsHV1"

####Extract GIMAP/IAN proteins and CgIAPs from Significantly Differentially Expressed Genes####
#Significantly differentially Expressed IAPs, using several names
oshv1_transcriptIDs_UniProt_SIG_ProtNames <- oshv1_transcriptIDs_UniProt_SIG$Protein.names
oshv1_IAPs_SIG <- grepl("IAP", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_IAPs_SIG) #3, 355 

oshv1_apoptosis_SIG <- grepl("apoptosis", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_apoptosis_SIG) #296
oshv1_apoptosis_SIG_info <- oshv1_transcriptIDs_UniProt_SIG[296,]

oshv1_inhibitor_SIG <-  grepl("inhibitor", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_inhibitor_SIG) #24, 152, 296, 351 (just 296!)

oshv1_IAPs_SIG_info <- oshv1_transcriptIDs_UniProt_SIG[c(3, 355, 296),]
oshv1_IAPs_SIG_info


#Significant GIMAP Genes
oshv1_IAN_Sig <- grepl("IAN", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_IAN_Sig) #0 TRUE
oshv1_AIG_Sig <- grepl("AIG", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_AIG_Sig) #0 TRUE
oshv1_IMAP_Sig <- grepl("IMAP", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_IMAP_Sig) #53, 65 TRUE

oshv1_GTP_SIG <- grepl("GTP", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_GTP_SIG) # 53  54  65  85 125
oshv1_transcriptIDs_UniProt_SIG[c(53,  54,  65,  85, 125),] # 53, 65 are GIMAP

oshv1_GIMAP_Sig_info <- oshv1_transcriptIDs_UniProt_SIG[c(53,65),]

#Uploaded NON Sig Differentially Expressed transcripts from the UniProt.ws website ####
oshv1_transcriptIDs_UniProt_non_Sig <- read.csv("OsHV1_resoshv1_Tran_df_non_Sig_transcript_ID_string.csv", header=TRUE)
#make sure to upload as csv version so that the protein names load correctly
head(oshv1_transcriptIDs_UniProt_non_Sig)
oshv1_transcriptIDs_UniProt_non_Sig["Challenge"] <- "OsHV1"

####Extract GIMAP and IAP genes from NON significant genes, for comparison ####
##Expressed IAPs
#use grepl to find text strings 
oshv1_transcriptIDs_UniProt_non_Sig_ProtNames <- oshv1_transcriptIDs_UniProt_non_Sig$Protein.names
oshv1_IAPs_non_Sig <- grepl("IAP", oshv1_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_IAPs_non_Sig) 
  #REAL IAPs: 556  1548  3029  5950  6152  7165  8369  8382  8566 10061 10456
oshv1_apoptosis_non_Sig <- grepl("apoptosis", oshv1_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_apoptosis_non_Sig) 
  #REAL IAPS: 1065  1891  2294  2295  2823 6340  6433 7217 8127 8816  8945  9070  

oshv1_inhibitor_nonSig <-  grepl("inhibitor", oshv1_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_inhibitor_nonSig) 
#303   417   454   616   760  1116  1493  1558  1571  1658  1788  1802  1891  2294
#2295  2342  2394  2815  2823  3400  3443  3453  3589  3659  4190  4572  4586  4643
#4796  5168  5215  5271  5390  5962  5963  6089  6198  6257  6302  6340  6433  6714
#6752  6888  6917  7217  7330  7523  8127  8667  8683  8731  8816  8925  8945  9070
# 9079  9311  9566 10365 10420 10455 10907 10943 10974 11195 11273 11282 11433

  #full non Sig IAP list
oshv1_IAPs_non_sig_info <- oshv1_transcriptIDs_UniProt_non_Sig[c(556,  1548,  3029 , 5950,  6152,  7165,  8369,  8382,  8566, 10061, 10456,
                                   1891,  2294,  2295,  2823, 6340 , 6433, 7217, 8127 ,8816,  8945,  9070  ),]
oshv1_IAPs_non_sig_info


##Expressed GIMAP Genes, non significant
oshv1_GTP_non_Sig <- grepl("GTP", oshv1_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_GTP_non_Sig)
oshv1_GIMAP_non_sig <- grepl("IMAP",oshv1_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE)
grep("TRUE", oshv1_GIMAP_non_sig)
#145   532   597  2086  2512  6100  6212  6319  7044  8445  8892  8981 11164
oshv1_IAN_non_Sig <- grepl("IAN", oshv1_transcriptIDs_UniProt_non_Sig$Protein.names) 
grep("TRUE", oshv1_IAN_non_Sig) #0

oshv1_GIMAP_non_sig_info <- oshv1_transcriptIDs_UniProt_non_Sig[c(145,   532,   597,  2086,  2512,  6100,  6212 , 6319,  7044,  8445,  8892 , 8981, 11164),] 
oshv1_GIMAP_non_sig_info


####Link GIMAP and IAN genes, significant and non significant, with Expression Values####
#Sig IAPS
subset_oshv1_IAPs_SIG_info$Ensembl_Genomes_Transcript #EKC24074 EKC42449 EKC20774
subset_oshv1_IAPs_SIG_info <- oshv1_IAPs_SIG_info[,c(4,8,9,10)]
oshv1_IAP_SIG_info_resoshv1Tran_05_dfSig <- grep("EKC24074", rownames(resoshv1Tran_05_dfSig), ignore.case = TRUE) 
oshv1_IAP_SIG_info_resoshv1Tran_05_dfSig <- resoshv1Tran_05_dfSig[4,] #line 4 is what the first grep gave
oshv1_IAP_SIG_info_resoshv1Tran_05_dfSig2 <- grep("EKC42449", rownames(resoshv1Tran_05_dfSig), ignore.case = TRUE) #line 1112
oshv1_IAP_SIG_info_resoshv1Tran_05_dfSig2 <- resoshv1Tran_05_dfSig[1116,]
oshv1_IAP_SIG_info_resoshv1Tran_05_dfSig3 <- grep("EKC20774", rownames(resoshv1Tran_05_dfSig), ignore.case = TRUE) #926
oshv1_IAP_SIG_info_resoshv1Tran_05_dfSig3 <- resoshv1Tran_05_dfSig[930,]
oshv1_IAP_SIG_combined_EXP <- rbind(oshv1_IAP_SIG_info_resoshv1Tran_05_dfSig, oshv1_IAP_SIG_info_resoshv1Tran_05_dfSig2, oshv1_IAP_SIG_info_resoshv1Tran_05_dfSig3)
oshv1_IAP_SIG_combined_FULL <- cbind(oshv1_IAP_SIG_combined_EXP,subset_oshv1_IAPs_SIG_info)
oshv1_IAP_SIG_combined_FULL["Significance"] <- "Significant"
oshv1_IAP_SIG_combined_FULL["Type"] <- "IAP"

#Sig GIMAPs
subset_oshv1_GIMAP_SIG_info <- oshv1_GIMAP_Sig_info[,c(4,8,9,10)]
#EKC41832, EKC30713
oshv1_GIMAP_SIG_info_resoshv1Tran_05_dfSig <- grep("EKC41832", rownames(resoshv1Tran_05_dfSig), ignore.case = TRUE) 
oshv1_GIMAP_SIG_info_resoshv1Tran_05_dfSig <- resoshv1Tran_05_dfSig[147,]
oshv1_GIMAP_SIG_info_resoshv1Tran_05_dfSig2 <- grep("EKC30713", rownames(resoshv1Tran_05_dfSig), ignore.case = TRUE) #170
oshv1_GIMAP_SIG_info_resoshv1Tran_05_dfSig2 <- resoshv1Tran_05_dfSig[170,]
oshv1_GIMAP_SIG_combined_EXP <- rbind(oshv1_GIMAP_SIG_info_resoshv1Tran_05_dfSig, oshv1_GIMAP_SIG_info_resoshv1Tran_05_dfSig2)
oshv1_GIMAP_SIG_combined_FULL <- cbind(oshv1_GIMAP_SIG_combined_EXP, subset_oshv1_GIMAP_SIG_info)
oshv1_GIMAP_SIG_combined_FULL["Significance"] <- "Significant"
oshv1_GIMAP_SIG_combined_FULL["Type"] <- "GIMAP"

#Non Sig IAP
oshv1_IAPs_non_sig_info$Ensembl_Genomes_Transcript 
#EKC32934 EKC38724 EKC34718 EKC30031 EKC24792 EKC20773 EKC18369 EKC20239 EKC42441
#EKC37539 EKC25493 EKC29263 EKC25955 EKC41180 EKC41181 EKC33184 EKC29824 EKC17690
#EKC34720 EKC34022 EKC26454 EKC26950 EKC42442

subset_oshv1_IAPs_non_sig_info <- oshv1_IAPs_non_sig_info[,c(4,8,9,10)]
oshv1_IAP_non_sig_info_resoshv1Tran_05_df2 <- grep("EKC32934", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) #2037
oshv1_IAP_non_sig_info_resoshv1Tran_05_df2 <- resoshv1Tran_05_df_non_Sig[2037,]
oshv1_IAP_non_sig_info_resoshv1Tran_05_df3 <- grep("EKC38724", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) #5678
oshv1_IAP_non_sig_info_resoshv1Tran_05_df3 <- resoshv1Tran_05_df_non_Sig[5678,]
oshv1_IAP_non_sig_info_resoshv1Tran_05_df4 <- grep("EKC34718", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) #11264
oshv1_IAP_non_sig_info_resoshv1Tran_05_df4 <- resoshv1Tran_05_df_non_Sig[11264,]
oshv1_IAP_non_sig_info_resoshv1Tran_05_df5 <- grep("EKC30031", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) #22345
oshv1_IAP_non_sig_info_resoshv1Tran_05_df5 <- resoshv1Tran_05_df_non_Sig[22345,]
oshv1_IAP_non_sig_info_resoshv1Tran_05_df6 <- grep("EKC24792", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) 
oshv1_IAP_non_sig_info_resoshv1Tran_05_df6 <- resoshv1Tran_05_df_non_Sig[23049,]
oshv1_IAP_non_sig_info_resoshv1Tran_05_df7 <- grep("EKC20773", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) 
oshv1_IAP_non_sig_info_resoshv1Tran_05_df7 <- resoshv1Tran_05_df_non_Sig[26798,]
oshv1_IAP_non_sig_info_resoshv1Tran_05_df8 <- grep("EKC18369", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) 
oshv1_IAP_non_sig_info_resoshv1Tran_05_df8 <- resoshv1Tran_05_df_non_Sig[31101,]
oshv1_IAP_non_sig_info_resoshv1Tran_05_df9 <- grep("EKC20239", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) 
oshv1_IAP_non_sig_info_resoshv1Tran_05_df9 <- resoshv1Tran_05_df_non_Sig[31152,]
oshv1_IAP_non_sig_info_resoshv1Tran_05_df10 <- grep("EKC42441", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) 
oshv1_IAP_non_sig_info_resoshv1Tran_05_df10 <- resoshv1Tran_05_df_non_Sig[31856,]
oshv1_IAP_non_sig_info_resoshv1Tran_05_df11 <- grep("EKC37539", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) #
oshv1_IAP_non_sig_info_resoshv1Tran_05_df11 <- resoshv1Tran_05_df_non_Sig[37159,]
oshv1_IAP_non_sig_info_resoshv1Tran_05_df12 <- grep("EKC25493", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) #
oshv1_IAP_non_sig_info_resoshv1Tran_05_df12 <- resoshv1Tran_05_df_non_Sig[38485,]
oshv1_IAP_non_sig_info_resoshv1Tran_05_df15 <- grep("EKC25955", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) #
oshv1_IAP_non_sig_info_resoshv1Tran_05_df15 <- resoshv1Tran_05_df_non_Sig[6982,]
oshv1_IAP_non_sig_info_resoshv1Tran_05_df16 <- grep("EKC41180", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) #
oshv1_IAP_non_sig_info_resoshv1Tran_05_df16 <- resoshv1Tran_05_df_non_Sig[8465,]
oshv1_IAP_non_sig_info_resoshv1Tran_05_df17 <- grep("EKC41181", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) #
oshv1_IAP_non_sig_info_resoshv1Tran_05_df17 <- resoshv1Tran_05_df_non_Sig[8466,]
oshv1_IAP_non_sig_info_resoshv1Tran_05_df18 <- grep("EKC33184", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) #
oshv1_IAP_non_sig_info_resoshv1Tran_05_df18 <- resoshv1Tran_05_df_non_Sig[10471,]
oshv1_IAP_non_sig_info_resoshv1Tran_05_df19 <- grep("EKC29824", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) #
oshv1_IAP_non_sig_info_resoshv1Tran_05_df19 <- resoshv1Tran_05_df_non_Sig[23795,]
oshv1_IAP_non_sig_info_resoshv1Tran_05_df20 <- grep("EKC17690", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) #
oshv1_IAP_non_sig_info_resoshv1Tran_05_df20 <- resoshv1Tran_05_df_non_Sig[24103,]
oshv1_IAP_non_sig_info_resoshv1Tran_05_df21 <- grep("EKC34720", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) #
oshv1_IAP_non_sig_info_resoshv1Tran_05_df21 <- resoshv1Tran_05_df_non_Sig[26969,]
oshv1_IAP_non_sig_info_resoshv1Tran_05_df22 <- grep("EKC34022", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) #
oshv1_IAP_non_sig_info_resoshv1Tran_05_df22 <- resoshv1Tran_05_df_non_Sig[30277,]
oshv1_IAP_non_sig_info_resoshv1Tran_05_df23 <- grep("EKC26454", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) #34851
oshv1_IAP_non_sig_info_resoshv1Tran_05_df23 <- resoshv1Tran_05_df_non_Sig[32754,]
oshv1_IAP_non_sig_info_resoshv1Tran_05_df24 <- grep("EKC26950", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) #34851
oshv1_IAP_non_sig_info_resoshv1Tran_05_df24 <- resoshv1Tran_05_df_non_Sig[33276,]
oshv1_IAP_non_sig_info_resoshv1Tran_05_df25 <- grep("EKC42442", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) #34851
oshv1_IAP_non_sig_info_resoshv1Tran_05_df25 <- resoshv1Tran_05_df_non_Sig[33643,]


oshv1_IAP_non_sig_combined_EXP <- rbind(oshv1_IAP_non_sig_info_resoshv1Tran_05_df2, oshv1_IAP_non_sig_info_resoshv1Tran_05_df3,
                                        oshv1_IAP_non_sig_info_resoshv1Tran_05_df4, oshv1_IAP_non_sig_info_resoshv1Tran_05_df5,
                                        oshv1_IAP_non_sig_info_resoshv1Tran_05_df6, oshv1_IAP_non_sig_info_resoshv1Tran_05_df7,
                                        oshv1_IAP_non_sig_info_resoshv1Tran_05_df8, oshv1_IAP_non_sig_info_resoshv1Tran_05_df9,
                                        oshv1_IAP_non_sig_info_resoshv1Tran_05_df10, oshv1_IAP_non_sig_info_resoshv1Tran_05_df11,
                                        oshv1_IAP_non_sig_info_resoshv1Tran_05_df12,
                                        oshv1_IAP_non_sig_info_resoshv1Tran_05_df15, oshv1_IAP_non_sig_info_resoshv1Tran_05_df16, 
                                        oshv1_IAP_non_sig_info_resoshv1Tran_05_df17, oshv1_IAP_non_sig_info_resoshv1Tran_05_df18,
                                        oshv1_IAP_non_sig_info_resoshv1Tran_05_df19, oshv1_IAP_non_sig_info_resoshv1Tran_05_df20,
                                        oshv1_IAP_non_sig_info_resoshv1Tran_05_df21, oshv1_IAP_non_sig_info_resoshv1Tran_05_df22,
                                        oshv1_IAP_non_sig_info_resoshv1Tran_05_df23, oshv1_IAP_non_sig_info_resoshv1Tran_05_df24,
                                        oshv1_IAP_non_sig_info_resoshv1Tran_05_df24)
oshv1_IAP_non_sig_combined_FULL <- cbind(oshv1_IAP_non_sig_combined_EXP, subset_oshv1_IAPs_non_sig_info)
oshv1_IAP_non_sig_combined_FULL["Significance"] <- "Non_significant"
oshv1_IAP_non_sig_combined_FULL["Type"] <- "IAP"

#Non Sig GIMAP
oshv1_GIMAP_non_sig_info$Ensembl_Genomes_Transcript
# EKC39736 EKC40465 EKC40820 EKC32489 EKC36405 EKC39748 EKC27363 EKC31739 EKC35292 EKC29604 EKC41613 EKC42724 EKC38639
subset_GIMAP_non_sig_info <- oshv1_GIMAP_non_sig_info[,c(4,8,9,10)]
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df <- grep("EKC39736", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE)
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df <- resoshv1Tran_05_df_non_Sig[509,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df2 <- grep("EKC40465", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) 
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df2 <- resoshv1Tran_05_df_non_Sig[1874,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df3 <- grep("EKC40820", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) 
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df3 <- resoshv1Tran_05_df_non_Sig[2137,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df4 <- grep("EKC32489", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) 
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df4 <- resoshv1Tran_05_df_non_Sig[7644,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df5 <- grep("EKC36405", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) 
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df5 <- resoshv1Tran_05_df_non_Sig[9306,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df6 <- grep("EKC39748", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) 
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df6 <- resoshv1Tran_05_df_non_Sig[22765,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df7 <- grep("EKC27363", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) 
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df7 <- resoshv1Tran_05_df_non_Sig[23223,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df8 <- grep("EKC31739", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) 
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df8 <- resoshv1Tran_05_df_non_Sig[23681,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df9 <- grep("EKC35292", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) 
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df9 <- resoshv1Tran_05_df_non_Sig[26534,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df10 <- grep("EKC29604", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE)
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df10 <- resoshv1Tran_05_df_non_Sig[31420,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df11 <- grep("EKC41613", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) 
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df11 <- resoshv1Tran_05_df_non_Sig[33035,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df12 <- grep("EKC42724", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) 
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df12 <- resoshv1Tran_05_df_non_Sig[33389,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df13 <- grep("EKC38639", rownames(resoshv1Tran_05_df_non_Sig), ignore.case = TRUE) 
oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df13 <- resoshv1Tran_05_df_non_Sig[41342,]

oshv1_GIMAP_non_sig_combined_EXP <- rbind(oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df, oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df2,
                                          oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df3, oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df4,
                                          oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df5, oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df6,
                                          oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df7, oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df8,
                                          oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df9, oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df10,
                                          oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df11,oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df12,
                                          oshv1_GIMAP_non_sig_info_resoshv1Tran_05_df13)

oshv1_GIMAP_non_sig_combined_FULL <- cbind(oshv1_GIMAP_non_sig_combined_EXP, subset_GIMAP_non_sig_info)
oshv1_GIMAP_non_sig_combined_FULL["Significance"] <- "Non significant"
oshv1_GIMAP_non_sig_combined_FULL["Type"] <- "GIMAP"

####Combine all Oshv1 IAP and GIMAP data #####
#to rbind the names of all the columns need to be the same, change first column name
oshv1_GIMAP_IAP_combined_FULL <- rbind(oshv1_GIMAP_non_sig_combined_FULL, oshv1_IAP_non_sig_combined_FULL,
                                       oshv1_IAP_SIG_combined_FULL, oshv1_GIMAP_SIG_combined_FULL)


####Bacterial Challenge Differential Gene Expression Analysis ####
#Gram negative SRA's (including control): (SRR796597, SRR796596, SRR796595, SRR796594, SRR796593, SRR796592, SRR796589)
# Gram positive (including control): (SRR796598,  SRR796589)

#Subset Bacterial Challenge Data
BacTranCountData <- as.matrix(TranscriptCountData[ , c(31:38)])
head(BacTranCountData)
BacTranColData <- read.csv("bac_PHENO_DATA.csv", header=TRUE, sep=",")
rownames(BacTranColData) <- BacTranColData$sampleID
colnames(BacTranCountData) <- BacTranColData$sampleID
head(BacTranCountData)
print(BacTranColData)

# Check all sample IDs in BacColData are also in BacCountData and match their orders
all(rownames(BacTranColData) %in% colnames(BacTranCountData))  #Should return TRUE
# returns TRUE
all(rownames(BacTranColData) == colnames(BacTranCountData))    # should return TRUE
#returns TRUE

#Give the "level" column levels
BacTranColData$level <- factor(BacTranColData$level)
levels(BacTranColData$level) #check to see that it has levels 

#### Create Bacterial Challenge DESeq Data Set from Matrix ####
# DESeqDataSet from count matrix and labels 
ddsBacTran <- DESeqDataSetFromMatrix(countData = BacTranCountData, 
                                     colData = BacTranColData, 
                                     design = ~condition)

#this design will gather the effect of condition between control and each treatment
# review how the data set looks
head(ddsBacTran)

#Relevel each to make sure that control is the first level in the condition factor
ddsBacTran$condition <- relevel( ddsBacTran$condition, "control")

#Check we're looking at the right samples
as.data.frame( colData(ddsBacTran) )

#Running the DEG pipeline
ddsBacTran<- DESeq(ddsBacTran) #for designs with interactions, recommends setting betaPrior=FALSE

#Inspect results table
resultsNames(ddsBacTran)

#extract contrasts between control and treatment values
resBacTran <- results(ddsBacTran, contrast = c("condition", "control", "treatment"))
#to extract log2fold change and p values under 0.1 and 0.05
head(resBacTran)
summary(resBacTran)
sum(resBacTran$padj < 0.1, na.rm=TRUE) #1756
#Change alpha setting in DESeq Results
resBacTran_05 <- results(ddsBacTran,alpha=0.05)
summary(resBacTran_05)
sum(resBacTran_05$padj < 0.05, na.rm=TRUE) #1437

#metadata on meaning of the columns
mcols(resBacTran_05, use.names = TRUE)
#shows treatment vs. control with baseMean as the mean of normalized counts for all samples

####p-value correction#### adapted from "Differential expression analysis of RNA-Seq data using DESeq2" Klaus 2014
#First Visualize histograms
#histogram of P- values to visualize any "hills" or "U shape"
# hill means variance of the null distribution too high, U shape means variance assumed too low
hist(resBacTran_05$pvalue, breaks = 20, col = "grey") #hill

#remove filtered out genes by independent filtering, they have NA adj. pvals
resBacTran_05_df <- resBacTran_05[ !is.na(resBacTran_05$padj), ]

#remove genes with NA pvals (outliers)
resBacTran_05_df <- resBacTran_05_df[ !is.na(resBacTran_05_df$pvalue), ]

#remove adjsuted pvalues, since we add the fdrtool results later on (based on the correct p-values)
resBacTran_05_df <- resBacTran_05_df[, -which(names(resBacTran_05_df) == "padj")]

#use z-scores as input to FDRtool to re-estimate the p-value
FDR.resBacTran_05_df <- fdrtool(resBacTran_05_df$stat, statistic= "normal", plot = T)

#add values to the results data frame, also ad new BH- adjusted p-values
resBacTran_05_df[,"padj"] <- p.adjust(FDR.resBacTran_05_df$pval, method = "BH")

#replot corrected p-values
hist(FDR.resBacTran_05_df$pval, col = "royalblue4",
     main = "Correct null model Bacterial Challenge Transcript Counts", xlab = "CORRECTED p-values")

#Check how many genes have BH adjusted p values of less than 0.05 after P-value correction?
sum( resBacTran_05_df$padj < 0.05, na.rm=TRUE ) #1096 (before p-value correction it was 1437)

#Subset the results table to the differentially expressed genes under FDR 0.01, order the Log2FC table first by strongest down regulation
resBacTran_05_dfSig <- resBacTran_05_df[ which(resBacTran_05_df$padj < 0.05 ), ]
head( resBacTran_05_dfSig[ order( resBacTran_05_dfSig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resBacTran_05_dfSig[ order( resBacTran_05_dfSig$log2FoldChange ), ] ) #tail for strongest up regulation
resBacTran_05_df_non_Sig <- resBacTran_05_df[ which(resBacTran_05_df$padj > 0.05), ]
resBacTran_05_df_non_Sig

####Visualize Results with Diagnostic Plots####
#MA plot, useful overview for experiment with two-group comparison. Plots log2FC over mean of normalized counts
#genes with adjusted p value 
plotMA(resBacTran_05_df)
plotMA(resBacTran_05_dfSig)
plotMA(resBacTran_05_df_non_Sig)

#Export Results to CSV
write.csv( as.data.frame(resBacTran_05_df), file="resBacTran_05_df.csv")
write.csv( as.data.frame(resBacTran_05_dfSig), file="resBacTran_05_dfSig.csv")
write.csv( as.data.frame(resBacTran_05_df_non_Sig), file="resBacTran_05_df_non_Sig.csv")

####Subset files for only those that have Transcript IDs####
#Extract gene titles of all the significantly differentially expressed genes
Bac_withID_subset_resBacTran_05_dfSig <- 
  resBacTran_05_dfSig[grep("transcript:", rownames(resBacTran_05_dfSig)), ]
Bac_withID_subset_resBacTran_dfSig #180 rows
BacTranscriptIDdf <- as.data.frame(rownames(Bac_withID_subset_resBacTran_05_dfSig))
head(BacTranscriptIDdf)
BacTranscriptIDdf= transform(BacTranscriptIDdf, 
                          ID = colsplit(rownames(Bac_withID_subset_resBacTran_05_dfSig), 
                                        split = "\\:", names = c('transcript:', 'EKC36484')))
#this line splitsup the word transcript and the transcript name

BacTranscriptIDstring <- toString(BacTranscriptIDdf[,3], sep=',')
write(BacTranscriptIDstring, "BacTranscriptIDstring_sig", sep = ",")
#write this to a file and then perform look up on the UniProt website

#Extract gene titles from non Sig genes
Bac_withID_subset_resBacTran_05_df_non_Sig <- 
  resBacTran_05_df_non_Sig[grep("transcript:", rownames(resBacTran_05_df_non_Sig)), ]
head(Bac_withID_subset_resBacTran_05_df_non_Sig)
BacTranscriptIDdf_nonSig <- as.data.frame(rownames(Bac_withID_subset_resBacTran_05_df_non_Sig))
head(BacTranscriptIDdf_nonSig)
BacTranscriptIDdf_nonSig= transform(BacTranscriptIDdf_nonSig, 
                                 ID = colsplit(rownames(Bac_withID_subset_resBacTran_05_df_non_Sig), 
                                               split = "\\:", names = c('transcript:', 'EKC36261')))
BacTranscriptIDstring_nonSig <- toString(BacTranscriptIDdf_nonSig[,3], sep=',')
BacTranscriptIDstring_nonSig
write(BacTranscriptIDstring_nonSig, "BacTranscriptIDstring_nonSig", sep = ",")
#write this to a file and then perform look up on the UniProt website

####Upload Sig Differentially Expressed transcripts from the UniProt.ws website ####
Bac_transcriptIDs_UniProt_SIG <- read.csv("BacTran_resBacTran_df_Sig_transcriptIDstring.csv", header=TRUE)
#make sure to upload as csv version so that the protein names load correctly
head(Bac_transcriptIDs_UniProt_SIG)
Bac_transcriptIDs_UniProt_SIG["Challenge"] <- "Bacteria"

####Extract GIMAP/IAN proteins and CgIAPs from Significantly Differentially Expressed Genes####
#Significantly differentially Expressed IAPs
Bac_transcriptIDs_UniProt_SIG_ProtNames <- Bac_transcriptIDs_UniProt_SIG$Protein.names
Bac_IAPs_SIG <- grepl("IAP", Bac_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", Bac_IAPs_SIG) #131...IAP 
Bac_apoptosis_SIG <- grepl("apoptosis", Bac_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", Bac_apoptosis_SIG) #44 
Bac_inhibitor_SIG <- grepl("inhibitor", Bac_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", Bac_inhibitor_SIG) #44, 145 (BAX)

Bac_IAP_SIG_info <- Bac_transcriptIDs_UniProt_SIG[c(44,131),]


#Significant GIMAP Genes
Bac_GTP_SIG <- grepl("GTP", Bac_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", Bac_GTP_SIG) #55 114 134 155 169 (none are GIMAP)
Bac_IMAP_SIG <- grepl("IMAP", Bac_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", Bac_IMAP_SIG) #0 
Bac_IAN_Sig <- grepl("IAN", Bac_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", Bac_IAN_Sig) #0 TRUE
Bac_immune_Sig <- grepl("immune", Bac_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", Bac_immune_Sig) #0 TRUE
Bac_AIG_Sig <- grepl("AIG", Bac_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", Bac_AIG_Sig) #0 TRUE
#NO SIG GIMAPS HERE

#Uploaded NON Sig Differentially Expressed transcripts from the UniProt.ws website ####
Bac_transcriptIDs_UniProt_non_Sig <- read.csv("BacTran_resBacTran_df_non_Sig_transcript_ID_string.csv", header=TRUE)
#make sure to upload as csv version so that the protein names load correctly
Bac_transcriptIDs_UniProt_non_Sig["Challenge"] <- "Bacteria"

####Extract GIMAP and IAP genes from NON significant genes, for comparison ####
#Expressed IAPs non sig
#use grepl to find text strings 
Bac_transcriptIDs_UniProt_non_Sig_ProtNames <- Bac_transcriptIDs_UniProt_non_Sig$Protein.names
Bac_IAPs_non_Sig <- grepl("IAP", Bac_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", Bac_IAPs_non_Sig) 
#36 1079 4110 4935 5791 5798 5927 5935 7244 all are BAC IAP

Bac_apoptosis_non_Sig <- grepl("apoptosis", Bac_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", Bac_apoptosis_non_Sig) #202  748 1598 4290 4443 4456 4972 5269 5622 6085 6180 6261 6697 6925
Bac_apoptosis_non_Sig_info <- Bac_transcriptIDs_UniProt_non_Sig[c(202,  748, 1598, 4290, 4443, 4456, 4972, 5269, 5622, 6085, 6180, 6261, 6697, 6925),]
#Actual Apoptosis Inhibitors: 202, 1598, 4443, 4972, 5622, 6085, 6180, 6261

Bac_IAPs_non_sig_info <- Bac_transcriptIDs_UniProt_non_Sig[c(36, 1079, 4110, 4935, 5791, 5798, 5927, 5935, 7244,
                                                             202, 1598, 4443, 4972, 5622, 6085, 6180, 6261 ),]

#Expressed GIMAP Genes, non significant
Bac_IAN_non_Sig <- grepl("IAN", Bac_transcriptIDs_UniProt_non_Sig$Protein.names) 
grep("TRUE", Bac_IAN_non_Sig) #0
Bac_immune_non_Sig <- grepl("immune", Bac_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", Bac_immune_non_Sig) #0 TRUE
Bac_AIG_non_Sig <- grepl("AIG", Bac_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", Bac_AIG_non_Sig) #0 TRUE
Bac_IMAP_non_sig <- grepl("IMAP",Bac_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE)
grep("TRUE", Bac_IMAP_non_sig)
#GIMAPS: 358  401  901 1028 1743 4286 4854 6141 6206 7716 (all are GIMAP)
Bac_GIMAP_non_sig_info <- Bac_transcriptIDs_UniProt_non_Sig[c(358,  401,  901, 1028 ,1743, 4286, 4854, 6141, 6206, 7716),]


####Link GIMAP and IAN genes, significant and non significant, with Expression Values####
#Sig IAPS
subset_Bac_IAPs_SIG_info <- Bac_IAP_SIG_info[,c(4,8,9,10)]
Bac_IAP_SIG_info_resBacTran_05_dfSig <- grep("EKC41180", rownames(resBacTran_05_dfSig), ignore.case = TRUE) # 848
Bac_IAP_SIG_info_resBacTran_05_dfSig <- resBacTran_05_dfSig[224,] 
Bac_IAP_SIG_info_resBacTran_05_dfSig2 <- grep("EKC38618", rownames(resBacTran_05_dfSig), ignore.case = TRUE) # 848
Bac_IAP_SIG_info_resBacTran_05_dfSig2 <- resBacTran_05_dfSig[811,] 
Bac_IAP_SIG_combined_EXP <- rbind(Bac_IAP_SIG_info_resBacTran_05_dfSig,Bac_IAP_SIG_info_resBacTran_05_dfSig2)
Bac_IAP_SIG_combined_FULL <- cbind(Bac_IAP_SIG_combined_EXP,subset_Bac_IAPs_SIG_info)
Bac_IAP_SIG_combined_FULL["Significance"] <- "Significant" #add column about significance
Bac_IAP_SIG_combined_FULL["Type"] <- "IAP"

#Sig GIMAPs
  #NO SIG GIMAPs

#Non Sig IAP
Bac_IAPs_non_sig_info$Ensembl_Genomes_Transcript 
#EKC24074 EKC38724 EKC30031 EKC20773 EKC18369 EKC20239 EKC42441 EKC42449 EKC25493 EKC18368
#EKC41181 EKC17690 EKC34720 EKC34022 EKC26454 EKC26950 EKC42442
subset_Bac_IAPs_non_sig_info <- Bac_IAPs_non_sig_info[,c(4,8,9,10)]
Bac_IAP_non_sig_info_resBacTran_05_df <- grep("EKC24074", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) #
Bac_IAP_non_sig_info_resBacTran_05_df <- resBacTran_05_df_non_Sig[80,]
Bac_IAP_non_sig_info_resBacTran_05_df2 <- grep("EKC38724", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) #
Bac_IAP_non_sig_info_resBacTran_05_df2 <- resBacTran_05_df_non_Sig[3991,]
Bac_IAP_non_sig_info_resBacTran_05_df3 <- grep("EKC30031", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) #
Bac_IAP_non_sig_info_resBacTran_05_df3 <- resBacTran_05_df_non_Sig[15556,]
Bac_IAP_non_sig_info_resBacTran_05_df4 <- grep("EKC20773", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) #
Bac_IAP_non_sig_info_resBacTran_05_df4 <- resBacTran_05_df_non_Sig[18692,]
Bac_IAP_non_sig_info_resBacTran_05_df5 <- grep("EKC18369", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) #
Bac_IAP_non_sig_info_resBacTran_05_df5 <- resBacTran_05_df_non_Sig[21699,]
Bac_IAP_non_sig_info_resBacTran_05_df6 <- grep("EKC20239", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) #
Bac_IAP_non_sig_info_resBacTran_05_df6 <- resBacTran_05_df_non_Sig[21745,]
Bac_IAP_non_sig_info_resBacTran_05_df7 <- grep("EKC42441", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) #
Bac_IAP_non_sig_info_resBacTran_05_df7 <- resBacTran_05_df_non_Sig[22242,]
Bac_IAP_non_sig_info_resBacTran_05_df8 <- grep("EKC42449", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) #
Bac_IAP_non_sig_info_resBacTran_05_df8 <- resBacTran_05_df_non_Sig[22258,]
Bac_IAP_non_sig_info_resBacTran_05_df9 <- grep("EKC25493", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) #
Bac_IAP_non_sig_info_resBacTran_05_df9 <- resBacTran_05_df_non_Sig[26841,]
Bac_IAP_non_sig_info_resBacTran_05_df10 <- grep("EKC18368", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) #
Bac_IAP_non_sig_info_resBacTran_05_df10 <- resBacTran_05_df_non_Sig[763,]
Bac_IAP_non_sig_info_resBacTran_05_df11 <- grep("EKC41181", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) #
Bac_IAP_non_sig_info_resBacTran_05_df11 <- resBacTran_05_df_non_Sig[5925,]
Bac_IAP_non_sig_info_resBacTran_05_df12 <- grep("EKC17690", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) #
Bac_IAP_non_sig_info_resBacTran_05_df12 <- resBacTran_05_df_non_Sig[16780,]
Bac_IAP_non_sig_info_resBacTran_05_df13 <- grep("EKC34720", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) #
Bac_IAP_non_sig_info_resBacTran_05_df13 <- resBacTran_05_df_non_Sig[18808,]
Bac_IAP_non_sig_info_resBacTran_05_df14 <- grep("EKC34022", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) #
Bac_IAP_non_sig_info_resBacTran_05_df14 <- resBacTran_05_df_non_Sig[21134,]
Bac_IAP_non_sig_info_resBacTran_05_df15 <- grep("EKC26454", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) #
Bac_IAP_non_sig_info_resBacTran_05_df15 <- resBacTran_05_df_non_Sig[22842,]
Bac_IAP_non_sig_info_resBacTran_05_df16 <- grep("EKC26950", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) #
Bac_IAP_non_sig_info_resBacTran_05_df16 <- resBacTran_05_df_non_Sig[23217,]
Bac_IAP_non_sig_info_resBacTran_05_df17 <- grep("EKC42442", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) #
Bac_IAP_non_sig_info_resBacTran_05_df17 <- resBacTran_05_df_non_Sig[23451,]

EKC29263

Bac_IAP_non_sig_combined_EXP <- rbind(Bac_IAP_non_sig_info_resBacTran_05_df, Bac_IAP_non_sig_info_resBacTran_05_df2,
                                      Bac_IAP_non_sig_info_resBacTran_05_df3, Bac_IAP_non_sig_info_resBacTran_05_df4,
                                      Bac_IAP_non_sig_info_resBacTran_05_df5,Bac_IAP_non_sig_info_resBacTran_05_df6,
                                      Bac_IAP_non_sig_info_resBacTran_05_df7,Bac_IAP_non_sig_info_resBacTran_05_df8,
                                      Bac_IAP_non_sig_info_resBacTran_05_df9, Bac_IAP_non_sig_info_resBacTran_05_df10,
                                      Bac_IAP_non_sig_info_resBacTran_05_df11,Bac_IAP_non_sig_info_resBacTran_05_df12,
                                      Bac_IAP_non_sig_info_resBacTran_05_df13,Bac_IAP_non_sig_info_resBacTran_05_df14,
                                      Bac_IAP_non_sig_info_resBacTran_05_df15,Bac_IAP_non_sig_info_resBacTran_05_df16,
                                      Bac_IAP_non_sig_info_resBacTran_05_df17)
Bac_IAP_non_sig_combined_FULL <- cbind(Bac_IAP_non_sig_combined_EXP,subset_Bac_IAPs_non_sig_info)
Bac_IAP_non_sig_combined_FULL["Significance"] <- "Non_significant" # add column about Significance
Bac_IAP_non_sig_combined_FULL["Type"] <- "IAP"

#Non Sig GIMAP
Bac_GIMAP_non_sig_info$Ensembl_Genomes_Transcript
#EKC40465 EKC40820 EKC41832 EKC30713 EKC36405 EKC27363 EKC35292 EKC41613 EKC42724 EKC38639
subset_Bac_GIMAP_non_sig_info <- Bac_GIMAP_non_sig_info[,c(4,8,9,10)]
Bac_GIMAP_non_sig_info_resBacTran_05_df <- grep("EKC40465", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) 
Bac_GIMAP_non_sig_info_resBacTran_05_df <- resBacTran_05_df_non_Sig[1321,]
Bac_GIMAP_non_sig_info_resBacTran_05_df2 <- grep("EKC40820", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) 
Bac_GIMAP_non_sig_info_resBacTran_05_df2 <- resBacTran_05_df_non_Sig[1499,]
Bac_GIMAP_non_sig_info_resBacTran_05_df3 <- grep("EKC41832", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE)
Bac_GIMAP_non_sig_info_resBacTran_05_df3 <- resBacTran_05_df_non_Sig[3277,]
Bac_GIMAP_non_sig_info_resBacTran_05_df4 <- grep("EKC30713", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) 
Bac_GIMAP_non_sig_info_resBacTran_05_df4 <- resBacTran_05_df_non_Sig[3761,]
Bac_GIMAP_non_sig_info_resBacTran_05_df5 <- grep("EKC36405", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) 
Bac_GIMAP_non_sig_info_resBacTran_05_df5 <- resBacTran_05_df_non_Sig[6491,]
Bac_GIMAP_non_sig_info_resBacTran_05_df6 <- grep("EKC27363", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) 
Bac_GIMAP_non_sig_info_resBacTran_05_df6 <- resBacTran_05_df_non_Sig[16169,]
Bac_GIMAP_non_sig_info_resBacTran_05_df7 <- grep("EKC35292", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) 
Bac_GIMAP_non_sig_info_resBacTran_05_df7 <- resBacTran_05_df_non_Sig[18513,]
Bac_GIMAP_non_sig_info_resBacTran_05_df8 <- grep("EKC41613", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) 
Bac_GIMAP_non_sig_info_resBacTran_05_df8 <- resBacTran_05_df_non_Sig[23047,]
Bac_GIMAP_non_sig_info_resBacTran_05_df9 <- grep("EKC42724", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) 
Bac_GIMAP_non_sig_info_resBacTran_05_df9 <- resBacTran_05_df_non_Sig[23284,]
Bac_GIMAP_non_sig_info_resBacTran_05_df10 <- grep("EKC38639", rownames(resBacTran_05_df_non_Sig), ignore.case = TRUE) 
Bac_GIMAP_non_sig_info_resBacTran_05_df10 <- resBacTran_05_df_non_Sig[28806,]

Bac_GIMAP_non_sig_combined_EXP <- rbind(Bac_GIMAP_non_sig_info_resBacTran_05_df, Bac_GIMAP_non_sig_info_resBacTran_05_df2,
                                        Bac_GIMAP_non_sig_info_resBacTran_05_df3, Bac_GIMAP_non_sig_info_resBacTran_05_df4,
                                        Bac_GIMAP_non_sig_info_resBacTran_05_df5, Bac_GIMAP_non_sig_info_resBacTran_05_df6,
                                        Bac_GIMAP_non_sig_info_resBacTran_05_df7,Bac_GIMAP_non_sig_info_resBacTran_05_df8,
                                        Bac_GIMAP_non_sig_info_resBacTran_05_df9, Bac_GIMAP_non_sig_info_resBacTran_05_df10)

Bac_GIMAP_non_sig_combined_FULL <- cbind(Bac_GIMAP_non_sig_combined_EXP,subset_Bac_GIMAP_non_sig_info)
Bac_GIMAP_non_sig_combined_FULL["Significance"] <- "Non_significant" #add column about significance
Bac_GIMAP_non_sig_combined_FULL["Type"] <- "GIMAP"

####Combine all Bac IAP and GIMAP data #####
#to rbind the names of all the columns need to be the same
Bac_GIMAP_IAP_combined_FULL <- rbind(Bac_IAP_SIG_combined_FULL, Bac_IAP_non_sig_combined_FULL,
                                     Bac_GIMAP_non_sig_combined_FULL)

#### COMPILE GIMAP AND IAP DATA FROM BOTH OSHV1 AND BAC TRIAL ####
#Check colnames
colnames(oshv1_GIMAP_IAP_combined_FULL)
colnames(Bac_GIMAP_IAP_combined_FULL)

COMBINED_GIMAP_IAP <- rbind(Bac_GIMAP_IAP_combined_FULL, oshv1_GIMAP_IAP_combined_FULL)

####Plotting the combined Significance values ####
COMBINED_GIMAP_IAP_cols <- COMBINED_GIMAP_IAP[,c(2,7,8,9,10,11,12)]

#Download data to simplify Protein Names for Viewing 
write.csv( as.data.frame(COMBINED_GIMAP_IAP_cols), file="COMBINED_GIMAP_IAP_cols.csv")
#Added transcript header
COMBINED_GIMAP_IAP_cols <- read.csv("COMBINED_GIMAP_IAP_cols.csv", header = TRUE)
COMBINED_GIMAP_IAP_cols_Bac <- COMBINED_GIMAP_IAP_cols[1:29,]
COMBINED_GIMAP_IAP_cols_Oshv1 <- COMBINED_GIMAP_IAP_cols[30:69,]
colnames(COMBINED_GIMAP_IAP_cols_Bac)
#Plot the combined data set altogether for comparison per Challenge
  #use geom_col and not geom_bar...geom bar makes the height of the bar proportional to the number
    #of cases in each group

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Plot GIMAP and IAP genes for the Bacterial Challenge Data
label.Bac.df <- data.frame(Transcript= c("transcript:EKC41180", "transcript:EKC38618"), log2FoldChange=c(21.5, 19.5))
COMBINED_GIMAP_IAP_BAC_PLOT <- ggplot(COMBINED_GIMAP_IAP_cols_Bac) + 
  geom_col(aes(x=Transcript, y=log2FoldChange, fill=as.factor(Type))) + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + 
  scale_y_continuous(name ="log2 Fold Change", breaks = scales::pretty_breaks(n = 20)) + 
  ggtitle("GIMAP and IAP Transcript Differential Expression in C. gigas under Bacterial Challenge") + 
  scale_fill_manual("Gene Family", values=(c("GIMAP"="#56B4E9", "IAP"="#009E73"))) +
  theme(plot.title = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5))

COMBINED_GIMAP_IAP_BAC_PLOT2 <- COMBINED_GIMAP_IAP_BAC_PLOT + geom_text(data=label.Bac.df, aes(x=Transcript, y=log2FoldChange), label = c("*")) + geom_hline(yintercept=0) +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  scale_x_discrete(name="Transcript Names (ENSEMBL Transcript ID)", limits=c("transcript:EKC24074",
 "transcript:EKC18369","transcript:EKC30031","transcript:EKC20773","transcript:EKC20239", "transcript:EKC38724",
 "transcript:EKC38618","transcript:EKC42449","transcript:EKC25493","transcript:EKC42441",
"transcript:EKC34022","transcript:EKC41181","transcript:EKC18368","transcript:EKC41180", "transcript:EKC42442",
"transcript:EKC17690", "transcript:EKC26950","transcript:EKC34720", "transcript:EKC26454","transcript:EKC27363",
"transcript:EKC35292","transcript:EKC36405","transcript:EKC38639","transcript:EKC40465","transcript:EKC40820",
"transcript:EKC41613","transcript:EKC41832", "transcript:EKC30713","transcript:EKC42724"),
labels=c("transcript:EKC24074"= "BIR-containing protein 2 (EKC24074)","transcript:EKC18369"="BIR-containing protein 3 (EKC18369)",
"transcript:EKC30031"="BIR-containing protein 3 (EKC30031)","transcript:EKC20773"="BIR-containing protein 3 Fragment (EKC20773)",
"transcript:EKC20239"="BIR-containing protein 6 (EKC20239)", "transcript:EKC38724"="BIR-containing protein 7 (EKC38724)",
"transcript:EKC38618"="BIR-containing protein 7-A (EKC38618)","transcript:EKC42449"="BIR-containing protein 7-A (EKC42449)",
"transcript:EKC25493"="BIR-containing protein 7-B (EKC25493)","transcript:EKC42441"="BIR-containing protein 7-B (EKC42441)",
 "transcript:EKC34022"="IAP 1 (EKC34022)", "transcript:EKC41181"="IAP 1 (EKC41181)","transcript:EKC18368"="IAP 2 (EKC18368)",
"transcript:EKC41180"="IAP 2 (EKC41180)","transcript:EKC42442"="IAP 2 (EKC42442)","transcript:EKC17690"="Putative IAP (EKC17690)",
"transcript:EKC26950"="Putative IAP (EKC26950)","transcript:EKC34720"="Putative IAP (EKC34720)",
"transcript:EKC26454"="Putative IAP ORF42 (EKC26454)","transcript:EKC27363"="GIMAP 4 (EKC27363)",
"transcript:EKC35292"="GIMAP 4 (EKC35292)","transcript:EKC36405"="GIMAP 4 (EKC36405)","transcript:EKC38639"="GIMAP 4 (EKC38639)",
"transcript:EKC40465"="GIMAP 4 (EKC40465)","transcript:EKC40820"="GIMAP 4 (EKC40820)","transcript:EKC41613"="GIMAP 4 (EKC41613)",
"transcript:EKC41832"="GIMAP 4 (EKC41832)","transcript:EKC30713"="GIMAP 7 (EKC30713)","transcript:EKC42724"="GIMAP 7 (EKC42724)"))

#Plot GIMAP and IAP from OsHV1
#Create vector for the significance values! 
label.oshv1Sig.df <- data.frame(Transcript= c("transcript:EKC240741","transcript:EKC424491", "transcript:EKC20774",
   "transcript:EKC418321", "transcript:EKC307131"), log2FoldChange=c(-3.0, -3.0, -24.0, -4.0, -3.0))

COMBINED_GIMAP_IAP_OsHV1_PLOT <- ggplot(COMBINED_GIMAP_IAP_cols_Oshv1) + 
  geom_col(aes(x=Transcript, y=log2FoldChange, fill=as.factor(Type))) +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + 
  scale_y_continuous(name ="log2 Fold Change", breaks = scales::pretty_breaks(n = 20)) + 
  ggtitle("GIMAP and IAP Transcript Differential Expression in C. gigas under OsHV1 Challenge") + 
  scale_fill_manual("Gene Family", values=(c("GIMAP"="#56B4E9", "IAP"="#009E73"))) +
  theme(plot.title = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5))

COMBINED_GIMAP_IAP_OsHV1_PLOT_2 <- COMBINED_GIMAP_IAP_OsHV1_PLOT + geom_text(data=label.oshv1Sig.df, aes(x=Transcript, y=log2FoldChange), label = c("*")) + geom_hline(yintercept=0) +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  scale_x_discrete(name="Transcript Names (ENSEMBL Transcript ID)", limits=c("transcript:EKC240741","transcript:EKC37539",
  "transcript:EKC183691","transcript:EKC300311","transcript:EKC207731", "transcript:EKC24792","transcript:EKC202391",
 "transcript:EKC307131","transcript:EKC387241","transcript:EKC424491", "transcript:EKC254931","transcript:EKC32934",
"transcript:EKC34718", "transcript:EKC424411","transcript:EKC20774","transcript:EKC29824","transcript:EKC340221",
 "transcript:EKC411811","transcript:EKC411801","transcript:EKC269501","transcript:EKC33184","transcript:EKC176901",
"transcript:EKC269502","transcript:EKC347201","transcript:EKC264541","transcript:EKC25955","transcript:EKC31739",
 "transcript:EKC273631","transcript:EKC32489","transcript:EKC352921", "transcript:EKC364051","transcript:EKC386391",
"transcript:EKC39736", "transcript:EKC39748","transcript:EKC404651","transcript:EKC408201","transcript:EKC416131",
"transcript:EKC418321","transcript:EKC427241","transcript:EKC29604"), 
labels=c("transcript:EKC240741"="BIR-containing protein 2 (EKC240741)","transcript:EKC37539"="BIR-containing protein 2 (EKC37539)",
"transcript:EKC183691"="BIR-containing protein 3 (EKC183691)","transcript:EKC300311"="BIR-containing protein 3 (EKC300311)",
 "transcript:EKC207731"="BIR-containing protein 3 Fragment (EKC207731)","transcript:EKC24792"="BIR-containing protein 5 (EKC24792)",
 "transcript:EKC202391"="BIR-containing protein 6 (EKC202391)","transcript:EKC307131"="BIR-containing protein 7 (EKC307131)",
"transcript:EKC387241"="BIR-containing protein 7 (EKC387241)","transcript:EKC424491"="BIR-containing protein 7-A (EKC424491)",
"transcript:EKC254931"="BIR-containing protein 7-B (EKC254931)","transcript:EKC32934"="BIR-containing protein 7-B (EKC329341)",
"transcript:EKC34718"="BIR-containing protein 7-B (EKC34718)","transcript:EKC424411"="BIR-containing protein 7-B (EKC424411)",
 "transcript:EKC20774"="IAP (EKC20774)","transcript:EKC29824"="IAP 1 (EKC29824)", "transcript:EKC340221"="IAP 1 (EKC340221)",
"transcript:EKC411811"="IAP 1 (EKC41181)","transcript:EKC411801"="IAP 2 (EKC41180)","transcript:EKC269501"="IAP 2 (EKC42442)",
  "transcript:EKC33184"="IAP 3 (EKC33184)","transcript:EKC176901"="Putative IAP (EKC176901)","transcript:EKC269502"="Putative IAP (EKC269502)",
"transcript:EKC347201"="Putative IAP (EKC347201)","transcript:EKC264541"="Putative IAP ORF42 (EKC264541)",
"transcript:EKC25955"="TP53-regulated IAP 1 (EKC25955)","transcript:EKC31739"="GIMAP 1 (EKC31739)",
"transcript:EKC273631"="GIMAP 4 (EKC273631)", "transcript:EKC32489"="GIMAP 4 (EKC32489)","transcript:EKC352921"="GIMAP 4 (EKC352921)",
"transcript:EKC364051"="GIMAP 4 (EKC364051)","transcript:EKC386391"="GIMAP 4 (EKC386391)","transcript:EKC39736"="GIMAP 4 (EKC39736)",
"transcript:EKC39748"="GIMAP 4 (EKC39748)","transcript:EKC404651"="GIMAP 4 (EKC404651)","transcript:EKC408201"="GIMAP 4 (EKC408201)",
"transcript:EKC416131"="GIMAP 4 (EKC416131)", "transcript:EKC418321"="GIMAP 4 (EKC418321)","transcript:EKC427241"="GIMAP 7 (EKC427241)",
"transcript:EKC29604"="GIMAP 8 (EKC29604)"))
  
  
#Plot IAP between OsHV1 and Bac
IAP_Bac <- COMBINED_GIMAP_IAP_cols_Bac %>% filter(Type=="IAP")
IAP_OsHV1 <- COMBINED_GIMAP_IAP_cols_Oshv1 %>% filter(Type=="IAP")
IAP_OsHV1_Bac <- rbind(IAP_Bac,IAP_OsHV1)
label.IAP.df <- data.frame(Transcript= c("transcript:EKC41180","transcript:EKC38618","transcript:EKC240741", "transcript:EKC424491","transcript:EKC20774","transcript:EKC41180","transcript:EKC38618","transcript:EKC240741","transcript:EKC424491","transcript:EKC20774"), log2FoldChange=c(21.5, 19.5, -3.0,-3.0,-24.0,21.5,19.2,-3.0,-3.0,-24.0))

IAP_OsHV1_Bac_PLOT <- ggplot(IAP_OsHV1_Bac) + 
  geom_col(aes(x=Transcript, y=log2FoldChange, fill=as.factor(Challenge))) + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + 
  scale_y_continuous(name ="log2 Fold Change", breaks = scales::pretty_breaks(n = 20)) + 
  ggtitle("IAP Transcript Differential Expression in C. gigas under OsHV1 and Bacterial Challenge") + 
  scale_fill_manual("Challenge", values=c("OsHV1"="#0072B2", "Bacteria"= "#D55E00")) +
  theme(plot.title = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5))

IAP_OsHV1_Bac_PLOT + geom_text(data=label.IAP.df, aes(x=Transcript, y=log2FoldChange), label = c("*")) + geom_hline(yintercept=0) +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  scale_x_discrete(name="Transcript Names (ENSEMBL Transcript ID)", 
  limits=c("transcript:EKC24074",
           "transcript:EKC240741",
           "transcript:EKC37539",
           "transcript:EKC18369",
           "transcript:EKC30031",
           "transcript:EKC20773",
          "transcript:EKC183691",
           "transcript:EKC300311",
           "transcript:EKC207731",
           "transcript:EKC24792",
           "transcript:EKC20239",
           "transcript:EKC202391",
           "transcript:EKC38724",
           "transcript:EKC387241",
           "transcript:EKC38618",
          "transcript:EKC42449",
           "transcript:EKC424491",
           "transcript:EKC25493",
           "transcript:EKC42441",
           "transcript:EKC254931",
           "transcript:EKC32934",
           "transcript:EKC34718",
           "transcript:EKC424411",
           "transcript:EKC17690",
           "transcript:EKC26950",
           "transcript:EKC34720",
           "transcript:EKC176901",
           "transcript:EKC269502",
           "transcript:EKC347201",
          "transcript:EKC26454",
           "transcript:EKC264541",
           "transcript:EKC25955",
           "transcript:EKC20774",
           "transcript:EKC34022",
           "transcript:EKC41181",
           "transcript:EKC29824",
           "transcript:EKC340221",
           "transcript:EKC411811",
           "transcript:EKC18368",
           "transcript:EKC41180",
           "transcript:EKC42442",
           "transcript:EKC411801",
           "transcript:EKC269501",
           "transcript:EKC33184"), 
labels=c("transcript:EKC24074"="BIR-containing protein 2 (EKC24074)",
         "transcript:EKC240741"="BIR-containing protein 2 (EKC24074)",
         "transcript:EKC37539"="BIR-containing protein 2 (EKC37539)",
         "transcript:EKC18369"="BIR-containing protein 3 (EKC18369)",
         "transcript:EKC30031"="BIR-containing protein 3 (EKC30031)",
         "transcript:EKC20773"="BIR-containing protein 3 Fragment (EKC20773)",
         "transcript:EKC183691"="BIR-containing protein 3 (EKC18369)",
         "transcript:EKC300311"="BIR-containing protein 3 (EKC30031)",
         "transcript:EKC207731"="BIR-containing protein 3 Fragment (EKC20773)",
         "transcript:EKC24792"="BIR-containing protein 5 (EKC24792)",
         "transcript:EKC20239"="BIR-containing protein 6 (EKC20239)",
         "transcript:EKC202391"="BIR-containing protein 6 (EKC20239)",
         "transcript:EKC38724"="BIR-containing protein 7 (EKC38724)",
         "transcript:EKC387241"="BIR-containing protein 7 (EKC38724)",
         "transcript:EKC38618"="BIR-containing protein 7-A (EKC38618)",
         "transcript:EKC42449"="BIR-containing protein 7-A (EKC42449)",
         "transcript:EKC424491"="BIR-containing protein 7-A (EKC42449)",
         "transcript:EKC25493"="BIR-containing protein 7-B (EKC25493)",
         "transcript:EKC42441"="BIR-containing protein 7-B (EKC42441)",
         "transcript:EKC254931"="BIR-containing protein 7-B (EKC25493)",
         "transcript:EKC32934"="BIR-containing protein 7-B (EKC32934)",
         "transcript:EKC34718"="BIR-containing protein 7-B (EKC34718)",
         "transcript:EKC424411"="BIR-containing protein 7-B (EKC42441)",
         "transcript:EKC17690"="Putative IAP (EKC17690)",
         "transcript:EKC26950"="Putative IAP (EKC26950)",
         "transcript:EKC34720"="Putative IAP (EKC34720)",
         "transcript:EKC176901"="Putative IAP (EKC17690)",
         "transcript:EKC269502"="Putative IAP (EKC26950)",
         "transcript:EKC347201"="Putative IAP (EKC34720)",
         "transcript:EKC26454"="Putative IAP ORF42 (EKC26454)",
         "transcript:EKC264541"="Putative IAP ORF42 (EKC26454)",
         "transcript:EKC25955"="TP53-regulated IAP 1 (EKC25955)",
         "transcript:EKC20774"="IAP (EKC20774)",
         "transcript:EKC34022"="IAP 1 (EKC34022)",
         "transcript:EKC41181"="IAP 1 (EKC41181)",
         "transcript:EKC29824"="IAP 1 (EKC29824)",
         "transcript:EKC340221"="IAP 1 (EKC34022)",
         "transcript:EKC411811"="IAP 1 (EKC41181)",
         "transcript:EKC18368"="IAP 2 (EKC18368)",
         "transcript:EKC41180"="IAP 2 (EKC41180)",
         "transcript:EKC42442"="IAP 2 (EKC42442)",
         "transcript:EKC411801"="IAP 2 (EKC41180)",
         "transcript:EKC269501"="IAP 2 (EKC42442)",
         "transcript:EKC33184"="IAP 3 (EKC33184)"))

#plot IAPs with the same transcript side by side
#Establish which transcripts are shared and which are different
IAP_OsHV1_Bac_duplicated<- duplicated(IAP_OsHV1_Bac$Ensembl_Genomes_Transcript) #none are duplicated!
grep("TRUE", IAP_OsHV1_Bac_duplicated) #21 23 25 26 27 28 30 32 33 36 37 38 39 40 41 42 43

#plot GIMAP between OsHV1 and Bac
GIMAP_Bac <- COMBINED_GIMAP_IAP_cols_Bac %>% filter(Type== "GIMAP")
GIMAP_OsHV1 <- COMBINED_GIMAP_IAP_cols_Oshv1 %>% filter(Type=="GIMAP")
GIMAP_OsHV1_Bac <- rbind(GIMAP_Bac,GIMAP_OsHV1)
label.GIMAP.df <- data.frame(Transcript= c("transcript:EKC418321","transcript:EKC307131"), log2FoldChange=c(-3.5, -2.5))

GIMAP_OsHV1_Bac_PLOT <- ggplot(GIMAP_OsHV1_Bac) + 
  geom_col(aes(x=Transcript, y=log2FoldChange, fill=as.factor(Challenge))) + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + 
  scale_y_continuous(name ="log2 Fold Change", breaks = scales::pretty_breaks(n = 20)) + 
  ggtitle("GIMAP Transcript Differential Expression in C. gigas under OsHV1 and Bacterial Challenge") + 
  scale_fill_manual("Challenge", values=c("Bacteria"="#0072B2", "OsHV1"= "#D55E00")) +
  theme(plot.title = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5)) 

GIMAP_OsHV1_Bac_PLOT2 <- GIMAP_OsHV1_Bac_PLOT + geom_text(data=label.GIMAP.df, aes(x=Transcript, y=log2FoldChange), label = c("*")) + geom_hline(yintercept=0) +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  scale_x_discrete(name="Transcript Names (ENSEMBL Transcript ID)", limits=c("transcript:EKC31739",
                                                                             "transcript:EKC27363",
                                                                             "transcript:EKC35292",
                                                                             "transcript:EKC36405",
                                                                             "transcript:EKC38639",
                                                                             "transcript:EKC40465",
                                                                            "transcript:EKC40820",
                                                                             "transcript:EKC41613",
                                                                             "transcript:EKC41832",
                                                                             "transcript:EKC273631",
                                                                             "transcript:EKC32489",
                                                                             "transcript:EKC352921",
                                                                             "transcript:EKC364051",
                                                                             "transcript:EKC386391",
                                                                             "transcript:EKC39736",
                                                                             "transcript:EKC39748",
                                                                             "transcript:EKC404651",
                                                                             "transcript:EKC408201",
                                                                             "transcript:EKC416131",
                                                                             "transcript:EKC418321",
                                                                             "transcript:EKC307131",
                                                                             "transcript:EKC427241",
                                                                             "transcript:EKC30713",
                                                                             "transcript:EKC42724",
                                                                             "transcript:EKC29604"), 
  labels=c("transcript:EKC31739"="GIMAP 1 (EKC31739)",
           "transcript:EKC27363"="GIMAP 4 (EKC27363)",
           "transcript:EKC35292"="GIMAP 4 (EKC35292)",
           "transcript:EKC36405"="GIMAP 4 (EKC36405)",
           "transcript:EKC38639"="GIMAP 4 (EKC38639)",
           "transcript:EKC40465"="GIMAP 4 (EKC40465)",
           "transcript:EKC40820"="GIMAP 4 (EKC40820)",
           "transcript:EKC41613"="GIMAP 4 (EKC41613)",
           "transcript:EKC41832"="GIMAP 4 (EKC41832)",
           "transcript:EKC273631"="GIMAP 4 (EKC27363)",
           "transcript:EKC32489"="GIMAP 4 (EKC32489)",
           "transcript:EKC352921"="GIMAP 4 (EKC35292)",
           "transcript:EKC364051"="GIMAP 4 (EKC36405)",
           "transcript:EKC386391"="GIMAP 4 (EKC38639)",
           "transcript:EKC39736"="GIMAP 4 (EKC39736)",
           "transcript:EKC39748"="GIMAP 4 (EKC39748)",
           "transcript:EKC404651"="GIMAP 4 (EKC40465)",
           "transcript:EKC408201"="GIMAP 4 (EKC40820)",
           "transcript:EKC416131"="GIMAP 4 (EKC41613)",
           "transcript:EKC418321"="GIMAP 4 (EKC41832)",
           "transcript:EKC307131"="GIMAP 7 (EKC30713)",
           "transcript:EKC427241"="GIMAP 7 (EKC42724)",
           "transcript:EKC30713"="GIMAP 7 (EKC30713)",
           "transcript:EKC42724"="GIMAP 7 (EKC42724)",
           "transcript:EKC29604"="GIMAP 8 (EKC29604)"))

#plot GIMAPs with same accession 
GIMAP_OsHV1_Bac_duplicated<- duplicated(GIMAP_OsHV1_Bac$Ensembl_Genomes_Transcript) 
grep("TRUE", GIMAP_OsHV1_Bac_duplicated) # 12 13 15 17 19 21 22 23 24 25

####Plot Duplicated Transcripts as a heatmap ####
#Using ggplot2
#melt dataframe I'm going to use
#https://www.r-bloggers.com/how-to-create-a-fast-and-easy-heatmap-with-ggplot2/
#sort the data greatest to least
COMBINED_GIMAP_IAP_cols_ordered <- COMBINED_GIMAP_IAP_cols[order(COMBINED_GIMAP_IAP_cols$log2FoldChange),]
COMBINED_GIMAP_IAP_cols_Bac_ordered <- COMBINED_GIMAP_IAP_cols_Bac[order(COMBINED_GIMAP_IAP_cols_Bac$log2FoldChange),]
COMBINED_GIMAP_IAP_cols_Oshv1_ordered <- COMBINED_GIMAP_IAP_cols_Oshv1[order(COMBINED_GIMAP_IAP_cols_Oshv1$log2FoldChange),] 

COMBINED_GIMAP_IAP_cols_ordered_log2FC <- COMBINED_GIMAP_IAP_cols_ordered[,c(2,4,7,9)]
COMBINED_GIMAP_IAP_cols_Bac_ordered_log2FC <- COMBINED_GIMAP_IAP_cols_Bac_ordered[,c(2,4,7,9)]
COMBINED_GIMAP_IAP_cols_Oshv1_ordered_log2FC <- COMBINED_GIMAP_IAP_cols_Oshv1_ordered[,c(2,4,7,9)]

#Name the rows by Protein.Transcript.Names
row.names(COMBINED_GIMAP_IAP_cols_ordered) <- COMBINED_GIMAP_IAP_cols_ordered$Protein.Transcript.Names
row.names(COMBINED_GIMAP_IAP_cols_Bac_ordered) <- COMBINED_GIMAP_IAP_cols_Bac_ordered$Protein.Transcript.Names
row.names(COMBINED_GIMAP_IAP_cols_Oshv1_ordered) <- COMBINED_GIMAP_IAP_cols_Oshv1_ordered$Protein.Transcript.Names

# Elaboration of heatmap (white - steelblue)
ggplot(COMBINED_GIMAP_IAP_cols_ordered_log2FC, aes(Challenge, Protein.Transcript.Names)) +
  geom_tile(aes(fill = log2FoldChange), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  ylab("Challenge") +
  xlab("Protein Transcripts") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "log2FoldChange")


####Investigate Transcripts from Other Gene Families####
##OSHV1## assuming Intrinsic pathway of apoptosis
#Caspases (caspase 2 is important with IAPS)
#Sig caspases (none) OsHV1

####STOPPED EDITING/UPDATING RESULTS HERE #### may not be relevant to do

oshv1_caspase_Sig <- grepl("caspase", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_caspase_Sig) #0
oshv1_protease_Sig <- grepl("protease", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_protease_Sig) #216 26S protease regulatory subunit 6A 
oshv1_caspase_Sig <- grepl("caspase", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_caspase_Sig) #0
oshv1_cysteine_Sig <- grepl("Cysteine", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_cysteinyl_Sig)
oshv1_aspartic_Sig <- grepl("aspartic", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_aspartic_Sig) #0
oshv1_cas_Sig <- grepl("Cas", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_cas_Sig) #136, 244 (neither!)

#Non sig caspases OsHV1
oshv1_caspase_non_sig <- grepl("caspase", oshv1_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE)
grep("TRUE", oshv1_caspase_non_sig)
  # 149  3227  3849  4735  5361  6070  6549  7465  8384  8737  9603  9856  9878 11369


#BcL 2 
#Sig Bcl2
oshv1_Bcl2_Sig <- grepl("Bcl", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_Bcl2_Sig) #245
oshv1_Bcl2_Sig_info <- oshv1_transcriptIDs_UniProt_SIG[245,]
#Non sig Bcl2
oshv1_Bcl2_non_Sig <- grepl("bcl", oshv1_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE)
grep("TRUE", oshv1_Bcl2_non_Sig) #313 1292 5409 6381
oshv1_Bcl2_non_Sig_info <- oshv1_transcriptIDs_UniProt_non_Sig[c( 313 ,6381),]

#BAG
#Sig BAG
oshv1_BAG_Sig <- grepl("BAG", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_BAG_Sig) #0
oshv1_athanogene_Sig <- grepl("athanogene", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_athanogene_Sig ) #0
#Non sig BAG
oshv1_BAG_non_Sig <- grepl("BAG",oshv1_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_BAG_non_Sig ) #2430  5568 11147
oshv1_BAG_non_Sig_info <- oshv1_transcriptIDs_UniProt_non_Sig[c(2430,  5568, 11147),]

#Cytochrome c
#Sig cytochrome 
oshv1_cytochrome_Sig <- grepl("cytochrome", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_cytochrome_Sig) #57 #this is only p 350
#Non sig cytochrome
oshv1_cytochrome_non_Sig <- grepl("cytochrome",oshv1_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_cytochrome_non_Sig )
#221 (cyt c oxidase) ,6438, 4081, 6438, 7902, 9665, 10310, 10311

oshv1_cytochrome_non_Sig_info <- oshv1_transcriptIDs_UniProt_non_Sig[c(221,6438, 4081, 6438, 7902, 9665, 10310, 10311),]

#Apoptotic process#
#Sig
oshv1_apoptotic_Sig <- grepl("apoptotic", oshv1_transcriptIDs_UniProt_SIG$Gene.ontology..GO., ignore.case = TRUE) 
grep("TRUE", oshv1_apoptotic_Sig) #3  45 233 245
oshv1_apoptotis_process_SIG_info <- oshv1_transcriptIDs_UniProt_SIG[c(3 , 45, 233, 245),]
oshv1_apoptotis_process_SIG_info <- oshv1_apoptotis_process_SIG_info %>%
  filter(Protein.names !="Uncharacterized protein")

#Non Sig
oshv1_apoptotic_non_Sig <- grepl("apoptotic",oshv1_transcriptIDs_UniProt_non_Sig$Gene.ontology..GO., ignore.case = TRUE) 
grep("TRUE", oshv1_apoptotic_non_Sig )
oshv1_apoptotic_process_non_Sig_info <-oshv1_transcriptIDs_UniProt_non_Sig[c(29,92,313,388,534,595,651,760,1342,  1389,  1492,  1707,  1742,  1865 , 2007  ,2605,  3096,  3260,  3344,
                                                                              3565,  4047,  4252,  4291,  4411,  4418 , 4857,  4944,  5082,  5263,  5312,  5335 , 5361,  5409,  5448,  5494,  5593,  5637,  5728,
                                                                              6381,  6549,  6565,  7131,  7465 , 7466,  7621,  7747,  7834,  7914 , 7928,  7984,  8320,  8384,  8414 , 8592,  8738,  8888 , 9344,
                                                                              9556,  9580 , 9601 , 9620 , 9952, 10265, 10273, 10441, 10625 ,10821, 10827, 10941, 11043, 11728),]
oshv1_apoptotic_process_non_Sig_info <- oshv1_apoptotic_process_non_Sig_info %>%
  filter(Protein.names != "Uncharacterized protein")


#Fas
#sig 
oshv1_Fas_non_sig <- grepl("Fas", oshv1_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE)
grep("TRUE", oshv1_Fas_non_sig)
oshv1_transcriptIDs_UniProt_non_Sig[c(1363, 1674, 3166, 3770, 3994, 6651, 95560),] # all are just Fas associated 
#non sig
oshv1_FAS_non_sig <- oshv1_transcriptIDs_UniProt_non_Sig[9556,]

###BAC challenge Additional Apoptotic Transcripts####
#Caspases (caspase 2 is important with IAPS)
#Sig caspases (none) OsHV1
oshv1_caspase_Sig <- grepl("caspase", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_caspase_Sig) #0
oshv1_protease_Sig <- grepl("protease", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_protease_Sig) #216 26S protease regulatory subunit 6A 
oshv1_caspase_Sig <- grepl("caspase", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_caspase_Sig) #0
oshv1_cysteine_Sig <- grepl("Cysteine", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_cysteinyl_Sig)
oshv1_aspartic_Sig <- grepl("aspartic", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_aspartic_Sig) #0
oshv1_cas_Sig <- grepl("Cas", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_cas_Sig) #136, 244 (neither!)

#Non sig caspases OsHV1
oshv1_caspase_non_sig <- grepl("caspase", oshv1_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE)
grep("TRUE", oshv1_caspase_non_sig)
# 149  3227  3849  4735  5361  6070  6549  7465  8384  8737  9603  9856  9878 11369


#BcL 2 
#Sig Bcl2
oshv1_Bcl2_Sig <- grepl("Bcl", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_Bcl2_Sig) #245
oshv1_Bcl2_Sig_info <- oshv1_transcriptIDs_UniProt_SIG[245,]
#Non sig Bcl2
oshv1_Bcl2_non_Sig <- grepl("bcl", oshv1_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE)
grep("TRUE", oshv1_Bcl2_non_Sig) #313 1292 5409 6381
oshv1_Bcl2_non_Sig_info <- oshv1_transcriptIDs_UniProt_non_Sig[c( 313 ,6381),]

#BAG
#Sig BAG
oshv1_BAG_Sig <- grepl("BAG", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_BAG_Sig) #0
oshv1_athanogene_Sig <- grepl("athanogene", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_athanogene_Sig ) #0
#Non sig BAG
oshv1_BAG_non_Sig <- grepl("BAG",oshv1_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_BAG_non_Sig ) #2430  5568 11147
oshv1_BAG_non_Sig_info <- oshv1_transcriptIDs_UniProt_non_Sig[c(2430,  5568, 11147),]

#Cytochrome c
#Sig cytochrome 
oshv1_cytochrome_Sig <- grepl("cytochrome", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_cytochrome_Sig) #57 #this is only p 350
#Non sig cytochrome
oshv1_cytochrome_non_Sig <- grepl("cytochrome",oshv1_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_cytochrome_non_Sig )
#221 (cyt c oxidase) ,6438, 4081, 6438, 7902, 9665, 10310, 10311

oshv1_cytochrome_non_Sig_info <- oshv1_transcriptIDs_UniProt_non_Sig[c(221,6438, 4081, 6438, 7902, 9665, 10310, 10311),]

#Apoptotic process#
#Sig
oshv1_apoptotic_Sig <- grepl("apoptotic", oshv1_transcriptIDs_UniProt_SIG$Gene.ontology..GO., ignore.case = TRUE) 
grep("TRUE", oshv1_apoptotic_Sig) #3  45 233 245
oshv1_apoptotis_process_SIG_info <- oshv1_transcriptIDs_UniProt_SIG[c(3 , 45, 233, 245),]
oshv1_apoptotis_process_SIG_info <- oshv1_apoptotis_process_SIG_info %>%
  filter(Protein.names !="Uncharacterized protein")

#Non Sig
oshv1_apoptotic_non_Sig <- grepl("apoptotic",oshv1_transcriptIDs_UniProt_non_Sig$Gene.ontology..GO., ignore.case = TRUE) 
grep("TRUE", oshv1_apoptotic_non_Sig )
oshv1_apoptotic_process_non_Sig_info <-oshv1_transcriptIDs_UniProt_non_Sig[c(29,92,313,388,534,595,651,760,1342,  1389,  1492,  1707,  1742,  1865 , 2007  ,2605,  3096,  3260,  3344,
                                                                             3565,  4047,  4252,  4291,  4411,  4418 , 4857,  4944,  5082,  5263,  5312,  5335 , 5361,  5409,  5448,  5494,  5593,  5637,  5728,
                                                                             6381,  6549,  6565,  7131,  7465 , 7466,  7621,  7747,  7834,  7914 , 7928,  7984,  8320,  8384,  8414 , 8592,  8738,  8888 , 9344,
                                                                             9556,  9580 , 9601 , 9620 , 9952, 10265, 10273, 10441, 10625 ,10821, 10827, 10941, 11043, 11728),]
oshv1_apoptotic_process_non_Sig_info <- oshv1_apoptotic_process_non_Sig_info %>%
  filter(Protein.names != "Uncharacterized protein")


#Fas
#sig 
oshv1_Fas_non_sig <- grepl("Fas", oshv1_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE)
grep("TRUE", oshv1_Fas_non_sig)
oshv1_transcriptIDs_UniProt_non_Sig[c(1363, 1674, 3166, 3770, 3994, 6651, 95560),] # all are just Fas associated 
#non sig
oshv1_FAS_non_sig <- oshv1_transcriptIDs_UniProt_non_Sig[9556,]


#references: 
#https://github.com/sr320/LabDocs/tree/master/code/DESeq
#StringTie manual : http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
