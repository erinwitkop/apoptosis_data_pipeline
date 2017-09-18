#05_Bac_viral_DESeq2_OsHV1_Bac_challenge

#This script takes as input the output Bac_Viral_transcript_count_matrix.csv data prepared from prepDE.py and performs
#differential Gene expression analysis, and subsets out isoforms of GIMAPs and CgIAPs and graphs their
#relative abundance.

#call the DESeq2 library 
source("https://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade") 
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
library(reshape)
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
head(oshv1TranCountData)
TranColData <- read.csv("Bac_Viral_PHENO_DATA.csv", header=TRUE, sep=",")
oshv1TranColData <- TranColData[c(1:30),]
oshv1TranColData <- oshv1TranColData[, c("sampleID", "condition", "stressorLevel")]
print(oshv1TranColData)
rownames(oshv1TranColData) <- oshv1TranColData$sampleID
colnames(oshv1TranCountData) <- oshv1TranColData$sampleID
head(oshv1TranCountData)
head(oshv1TranColData)
#Give the stressorLevel column levels
oshv1TranColData$stressorLevel <- factor(oshv1TranColData$stressorLevel)
levels(oshv1TranColData$stressorLevel) #check to see that it has levels 

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
sum(resoshv1Tran$padj < 0.1, na.rm=TRUE) #4409
resoshv1Tran_05 <- results(ddsOshv1Tran,alpha=0.05)
summary(resoshv1Tran_05)
sum(resoshv1Tran_05$padj < 0.05, na.rm=TRUE) #3282

#metadata on meaning of the columns
mcols(resoshv1Tran_05, use.names = TRUE)
#Get more detailed description
mcols(resoshv1Tran_05_df)$description
#shows treatment vs. control with baseMean as the mean of normalized counts for all samples

####p-value correction for all genes in resoshv1Tran ####
#First Visualize histograms
#histogram of P- values to visualize any "hills" or "U shape"
# hill means variance of the null distribution too high, U shape means variance assumed too low
hist(resoshv1Tran$pvalue, breaks = 20, col = "grey") #hill

#remove filtered out genes by independent filtering, they have NA adj. pvals
resoshv1Tran_df <- resoshv1Tran[ !is.na(resoshv1Tran$padj), ]

#remove genes with NA pvals (outliers)
resoshv1Tran_df <- resoshv1Tran_df[ !is.na(resoshv1Tran_df$pvalue), ]

#remove adjsuted pvalues, since we add the fdrtool results later on (based on the correct p-values)
resoshv1Tran_df <- resoshv1Tran_df[, -which(names(resoshv1Tran_df) == "padj")]

#use z-scores as input to FDRtool to re-estimate the p-value
FDR.resoshv1Tran_df <- fdrtool(resoshv1Tran_df$stat, statistic= "normal", plot = T)

#add values to the results data frame, also ad new BH- adjusted p-values
resoshv1Tran_df[,"padj"] <- p.adjust(FDR.resoshv1Tran_df$pval, method = "BH")

#replot corrected p-values 
hist(FDR.resoshv1Tran_df$pval, col = "royalblue4",
     main = "Correct null model OsHv1 Transcript Count", xlab = "CORRECTED p-values")

#Check how many genes have BH adjusted p values of less than 0.05 after P-value correction?
sum( resoshv1Tran_df$padj < 0.05, na.rm=TRUE ) #1453

#Subset the results table to the differentially expressed genes under FDR 0.1, order the Log2FC table first by strongest down regulation
resoshv1Tran_dfSig <- resoshv1Tran_df[ which(resoshv1Tran_df$padj < 0.05 ), ]
head( resoshv1Tran_dfSig[ order( resoshv1Tran_dfSig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resoshv1Tran_dfSig[ order( resoshv1Tran_dfSig$log2FoldChange ), ] ) #tail for strongest up regulation
summary(resoshv1Tran_dfSig)
resoshv1Tran_df_non_Sig <- resoshv1Tran_df[ which(resoshv1Tran_df$padj > 0.05 ), ]
summary(resoshv1Tran_df_non_Sig)

#Visualize Results with Diagnostic Plots#
#MA plot, useful overview for experiment with two-group comparison. Plots log2FC over mean of normalized counts
#genes with adjusted p value 
plotMA(resoshv1Tran_dfSig)
plotMA(resoshv1Tran_df_non_Sig)

#Export Results to CSV
write.csv( as.data.frame(resoshv1Tran_df), file="OsHV1_resoshv1Tran_df.csv")
write.csv( as.data.frame(resoshv1Tran_dfSig), file="OsHV1_resoshv1Tran_dfSig.csv")
write.csv( as.data.frame(resoshv1Tran_df_non_Sig), file = "OsHV1_resoshv1Tran_df_non_Sig.csv")


####Subset files for only those that have Transcript IDs####
#Extract gene titles of all the significantly differentially expressed genes
OsHV1_withID_subset_resoshv1Tran_dfSig <- 
  resoshv1Tran_dfSig[grep("transcript:", rownames(resoshv1Tran_dfSig)), ]
head(OsHV1_withID_subset_resoshv1Tran_dfSig)
transcriptIDdf <- as.data.frame(rownames(OsHV1_withID_subset_resoshv1Tran_dfSig))
head(transcriptIDdf)
transcriptID1thru5 <- rownames(head(OsHV1_withID_subset_resoshv1Tran_dfSig[1:5,]))
transcriptIDdf= transform(transcriptIDdf, 
        ID = colsplit(rownames(OsHV1_withID_subset_resoshv1Tran_dfSig), 
                      split = "\\:", names = c('transcript:', 'EKC33371')))
  #this line splitsup the word transcript and the transcript name

transcriptIDstring <- toString(transcriptIDdf[,3], sep=',')
transcriptIDstring
write(transcriptIDstring, "transcriptIDstring", sep = ",")
#write this to a file and then perform look up on the UniProt website

#Extract gene titles from non Sig genes
OsHV1_withID_subset_resoshv1Tran_df_non_Sig <- 
  resoshv1Tran_df_non_Sig[grep("transcript:", rownames(resoshv1Tran_df_non_Sig)), ]
head(OsHV1_withID_subset_resoshv1Tran_df_non_Sig)
transcriptIDdf_nonSig <- as.data.frame(rownames(OsHV1_withID_subset_resoshv1Tran_df_non_Sig))
head(transcriptIDdf_nonSig)
transcriptIDdf_nonSig= transform(transcriptIDdf_nonSig, 
                          ID = colsplit(rownames(OsHV1_withID_subset_resoshv1Tran_df_non_Sig), 
                                        split = "\\:", names = c('transcript:', 'EKC37466')))
transcriptIDstring_nonSig <- toString(transcriptIDdf_nonSig[,3], sep=',')
transcriptIDstring_nonSig
write(transcriptIDstring_nonSig, "transcriptIDstring_nonSig", sep = ",")
#write this to a file and then perform look up on the UniProt website

#Add quotes around this 
#transcriptIDparen <- sapply(strsplit(transcriptIDstring, '[, ]+'), function(x) toString(dQuote(x)))
#str(transcriptIDparen)

####For situations only looking at the gene name of a few genes, can look it up locally ####
#Load C. gigas UniProt ID information for C. gigas
source("https://bioconductor.org/biocLite.R")
biocLite("UniProt.ws")
library(UniProt.ws)
browseVignettes("UniProt.ws")
#Taxonomy ID for C. virginica: 29159
availableUniprotSpecies(pattern = "gigas")
CgigasUp <- UniProt.ws(29159)
CgigasUp
keytypes(CgigasUp)
columns(CgigasUp)
keytypes(CgigasUp)
showMethods("keys")
#structure: res <- select(UniProtName, keys, columns, keytype)
#select(up, keys=c("P31946","P62258"), columns=c("PDB","SEQUENCE"), keytype="ENSEMBL_GENOMES TRANSCRIPT)

#Create object to look up UNIPROTKB entries and retreive other gene info
columns <- c("ENSEMBL_PROTEIN","HGNC","SEQUENCE", "GENES","ENTREZ_GENE", "ID", "GO", "GO-ID", "KEGG")
keytype <- "ENSEMBL_GENOMES TRANSCRIPT"
#keys <- transcriptIDdf[,3] #set key to be  vector of all of your transcriptIDs from previous
keys <- transcriptIDparen[1:50] #test to see how this works
res <- select(CgigasUp, keys, columns, keytype)
res

####Otherwise, Uploaded Sig Differentially Expressed transcripts from the UniProt.ws website ####
oshv1_transcriptIDs_UniProt_SIG <- read.csv("OsHV1_resoshv1Tran_dfSig_transcriptIDstring.csv", header=TRUE)
#make sure to upload as csv version so that the protein names load correctly
head(oshv1_transcriptIDs_UniProt_SIG)

####Extract GIMAP/IAN proteins and CgIAPs from Significantly Differentially Expressed Genes####
#Significantly differentially Expressed IAPs, using several names
oshv1_transcriptIDs_UniProt_SIG_ProtNames <- oshv1_transcriptIDs_UniProt_SIG$Protein.names
oshv1_IAPs_SIG <- grepl("IAP", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_IAPs_SIG) #3, 353

oshv1_apoptosis_SIG <- grepl("apoptosis", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_apoptosis_SIG) #3, 353
oshv1_apoptosis_SIG_info <- oshv1_transcriptIDs_UniProt_SIG[294,]

oshv1_inhibitor_SIG <-  grepl("inhibitor", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_inhibitor_SIG) # 24 149 294 349
oshv1_inhibitor_SIG_info <- oshv1_transcriptIDs_UniProt_SIG[c(24, 149, 294, 349),] #only 294!

oshv1_IAPs_SIG_info <- oshv1_transcriptIDs_UniProt_SIG[c(3,353, 294),]
oshv1_IAPs_SIG_info


#Significant GIMAP Genes
oshv1_IAN_Sig <- grepl("IAN", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_IAN_Sig) #0 TRUE
oshv1_AIG_Sig <- grepl("AIG", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_AIG_Sig) #0 TRUE
oshv1_IMAP_Sig <- grepl("IMAP", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_IMAP_Sig) #0 TRUE

oshv1_GTP_SIG <- grepl("GTP", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_GIMAP_SIG)
oshv1_transcriptIDs_UniProt_SIG[c(53,54,65,85,125),] #53 and 65 are GIMAP! called GTPase IMAP
oshv1_GIMAP_Sig_info <- oshv1_transcriptIDs_UniProt_SIG[c(53,65),]
oshv1_GIMAP_Sig_info


#Uploaded NON Sig Differentially Expressed transcripts from the UniProt.ws website ####
oshv1_transcriptIDs_UniProt_non_Sig <- read.csv("OsHV1_resoshv1_Tran_df_non_Sig_transcript_ID_string.csv", header=TRUE)
#make sure to upload as csv version so that the protein names load correctly
head(oshv1_transcriptIDs_UniProt_non_Sig)

####Extract GIMAP and IAP genes from NON significant genes, for comparison ####
##Expressed IAPs
#use grepl to find text strings 
oshv1_transcriptIDs_UniProt_non_Sig_ProtNames <- oshv1_transcriptIDs_UniProt_non_Sig$Protein.names
oshv1_IAPs_non_Sig <- grepl("IAP", oshv1_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_IAPs_non_Sig) 
  #69   575  1584  3107  6112  6316  7355  8579  8592  8784 10334 10736 11378
oshv1_apoptosis_non_Sig <- grepl("apoptosis", oshv1_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_apoptosis_non_Sig) #318  1092  1934  2350  2351  2897  5229  6381  6508  6607  6626  7407  7834  8333  8738  9043
#[17]  9178  9308  9939 10273 10941
oshv1_apoptosis_non_Sig_info <- oshv1_transcriptIDs_UniProt_non_Sig[c(318,  1092,  1934,  2350, 
      2351,  2897,  5229,  6381,  6508,  6607,  6626,  7407,  7834,  8333,  8738,  9043,
                                                  9178,  9308,  9939, 10273, 10941),]
oshv1_apoptosis_non_Sig_info 
#ACTUAL IAP PROTEINS FROM oshv1_apoptosis_non_Sig_info
#318, 1934,2350,2351, 2897, 6508, 6607, 7407, 8333, 9043, 9178, 9308

oshv1_inhibitor_nonSig <-  grepl("inhibitor", oshv1_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_inhibitor_nonSig) 
#  318   435   473   637   782  1143  1527  1594  1608  1696  1831  1845  1934  2350  2351  2401
# 2453  2889  2897  3493  3539  3549  3688  3760  4309  4698  4715  4772  4925  5309  5358  5415
# 5537  6124  6125  6253  6363  6424  6469  6508  6607  6894  6933  7075  7105  7407  7522  7718
# 8333  8887  8903  8953  9043  9157  9178  9308  9317  9556  9820 10643 10700 10735 11192 11228
# 11261 11488 11568 11577 11729

oshv1_inhibitor_non_Sig_info <- oshv1_transcriptIDs_UniProt_non_Sig[c(
  318,   435,   473,   637,   782 , 1143 , 1527,  1594,  1608,  1696,  1831,  1845,1934,2350,  2351,  2401,
  2453,  2889,  2897,  3493,  3539,  3549,  3688,  3760,  4309,  4698,  4715,  4772,  4925,  5309,  5358,  5415,
  5537,  6124,  6125,  6253,  6363,  6424,  6469,  6508 , 6607,  6894,  6933 , 7075,  7105,  7407,  7522,  7718,
  8333,  8887 , 8903,  8953,  9043,  9157,  9178,  9308,  9317,  9556,  9820, 10643, 10700, 10735, 11192, 11228,
  11261, 11488, 11568, 11577, 11729),]
oshv1_inhibitor_non_Sig_info
#Rows that are actual Apoptosis inhibitor proteins IAP: 318, 1934, 2350,2351,2897,6508,6607,
  #7407, 8333, 9043, 9178, 9308, SAME AS ABOVE


#Compiled list from IAP searches above
oshv1_IAPs_non_sig_info <- oshv1_transcriptIDs_UniProt_non_Sig[c(575,1584,3107,
        6112,6316,7355,8579,8592,8784,10334,10736,11378, 318, 1934, 2350,2351,2897,6508,6607,
        7407, 8333, 9043, 9178, 9308),]
oshv1_IAPs_non_sig_info


##Expressed GIMAP Genes, non significant
oshv1_GTP_non_Sig <- grepl("GTP", oshv1_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_GTP_non_Sig)
oshv1_GIMAP_non_sig <- grepl("IMAP",oshv1_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE)
grep("TRUE", oshv1_GIMAP_non_sig)
#151   551   617  2134  2576  6264  6377  6487  7233  8658  9122  9214 11457
oshv1_GIMAP_non_sig_info <- oshv1_transcriptIDs_UniProt_non_Sig[c(151,551,617,2134,
      2576,  6264,  6377,  6487,  7233,  8658,  9122,  9214, 11457),] 
oshv1_GIMAP_non_sig_info
oshv1_IAN_non_Sig <- grepl("IAN", oshv1_transcriptIDs_UniProt_non_Sig$Protein.names) 
grep("TRUE", oshv1_IAN_non_Sig) #0

####Link GIMAP and IAN genes, significant and non significant, with Expression Values####
#Sig IAPS
subset_oshv1_IAPs_SIG_info <- oshv1_IAPs_SIG_info[,c(4,8,9)]
oshv1_IAP_SIG_info_resoshv1Tran_dfSig <- grep("EKC24074", rownames(resoshv1Tran_dfSig), ignore.case = TRUE) 
oshv1_IAP_SIG_info_resoshv1Tran_dfSig <- resoshv1Tran_dfSig[4,] #line 4 is what the first grep gave
oshv1_IAP_SIG_info_resoshv1Tran_dfSig2 <- grep("EKC42449", rownames(resoshv1Tran_dfSig), ignore.case = TRUE) #line 1112
oshv1_IAP_SIG_info_resoshv1Tran_dfSig2 <- resoshv1Tran_dfSig[1112,]
oshv1_IAP_SIG_info_resoshv1Tran_dfSig3 <- grep("EKC20774", rownames(resoshv1Tran_dfSig), ignore.case = TRUE) #926
oshv1_IAP_SIG_info_resoshv1Tran_dfSig3 <- resoshv1Tran_dfSig[926,]
oshv1_IAP_SIG_combined_EXP <- rbind(oshv1_IAP_SIG_info_resoshv1Tran_dfSig, oshv1_IAP_SIG_info_resoshv1Tran_dfSig2, oshv1_IAP_SIG_info_resoshv1Tran_dfSig3)
oshv1_IAP_SIG_combined_FULL <- cbind(oshv1_IAP_SIG_combined_EXP,subset_oshv1_IAPs_SIG_info)

#Sig GIMAPs
subset_oshv1_GIMAP_SIG_info <- oshv1_GIMAP_Sig_info[,c(4,8,9)]
#EKC41832, EKC30713
oshv1_GIMAP_SIG_info_resoshv1Tran_dfSig <- grep("EKC41832", rownames(resoshv1Tran_dfSig), ignore.case = TRUE) 
oshv1_GIMAP_SIG_info_resoshv1Tran_dfSig <- resoshv1Tran_dfSig[147,]
oshv1_GIMAP_SIG_info_resoshv1Tran_dfSig2 <- grep("EKC30713", rownames(resoshv1Tran_dfSig), ignore.case = TRUE) #170
oshv1_GIMAP_SIG_info_resoshv1Tran_dfSig2 <- resoshv1Tran_dfSig[170,]
oshv1_GIMAP_SIG_combined_EXP <- rbind(oshv1_GIMAP_SIG_info_resoshv1Tran_dfSig, oshv1_GIMAP_SIG_info_resoshv1Tran_dfSig2)
oshv1_GIMAP_SIG_combined_FULL <- cbind(oshv1_GIMAP_SIG_combined_EXP, subset_oshv1_GIMAP_SIG_info)


#Non Sig IAP
oshv1_IAPs_non_sig_info$Ensembl_Genome_Transcript 
#EKC32934 EKC38724 EKC34718 EKC30031 EKC24792 EKC20773 EKC18369 EKC20239 EKC42441 EKC37539
#EKC25493 EKC42443
#EKC18368 EKC25955 EKC41180 EKC41181 EKC33184 EKC29824 EKC17690 EKC34720
#EKC34022 EKC26454 EKC26950 EKC42442
subset_oshv1_IAPs_non_sig_info <- oshv1_IAPs_non_sig_info[,c(4,8,9)]
#deleted the first row (69)
oshv1_IAP_non_sig_info_resoshv1Tran_df2 <- grep("EKC32934", rownames(resoshv1Tran_df), ignore.case = TRUE) #2157
oshv1_IAP_non_sig_info_resoshv1Tran_df2 <- resoshv1Tran_df[2157,]
oshv1_IAP_non_sig_info_resoshv1Tran_df3 <- grep("EKC38724", rownames(resoshv1Tran_df), ignore.case = TRUE) #6044
oshv1_IAP_non_sig_info_resoshv1Tran_df3 <- resoshv1Tran_df[6044,]
oshv1_IAP_non_sig_info_resoshv1Tran_df4 <- grep("EKC34718", rownames(resoshv1Tran_df), ignore.case = TRUE) #11980
oshv1_IAP_non_sig_info_resoshv1Tran_df4 <- resoshv1Tran_df[11980,]
oshv1_IAP_non_sig_info_resoshv1Tran_df5 <- grep("EKC30031", rownames(resoshv1Tran_df), ignore.case = TRUE) #23757
oshv1_IAP_non_sig_info_resoshv1Tran_df5 <- resoshv1Tran_df[23757,]
oshv1_IAP_non_sig_info_resoshv1Tran_df6 <- grep("EKC24792", rownames(resoshv1Tran_df), ignore.case = TRUE) #24508
oshv1_IAP_non_sig_info_resoshv1Tran_df6 <- resoshv1Tran_df[24508,]
oshv1_IAP_non_sig_info_resoshv1Tran_df7 <- grep("EKC20773", rownames(resoshv1Tran_df), ignore.case = TRUE) #28522
oshv1_IAP_non_sig_info_resoshv1Tran_df7 <- resoshv1Tran_df[28522,]
oshv1_IAP_non_sig_info_resoshv1Tran_df8 <- grep("EKC18369", rownames(resoshv1Tran_df), ignore.case = TRUE) #33083
oshv1_IAP_non_sig_info_resoshv1Tran_df8 <- resoshv1Tran_df[33083,]
oshv1_IAP_non_sig_info_resoshv1Tran_df9 <- grep("EKC20239", rownames(resoshv1Tran_df), ignore.case = TRUE) #33137
oshv1_IAP_non_sig_info_resoshv1Tran_df9 <- resoshv1Tran_df[33137,]
oshv1_IAP_non_sig_info_resoshv1Tran_df10 <- grep("EKC42441", rownames(resoshv1Tran_df), ignore.case = TRUE) #33888
oshv1_IAP_non_sig_info_resoshv1Tran_df10 <- resoshv1Tran_df[33888,]
oshv1_IAP_non_sig_info_resoshv1Tran_df11 <- grep("EKC37539", rownames(resoshv1Tran_df), ignore.case = TRUE) #39530
oshv1_IAP_non_sig_info_resoshv1Tran_df11 <- resoshv1Tran_df[39530,]
oshv1_IAP_non_sig_info_resoshv1Tran_df12 <- grep("EKC25493", rownames(resoshv1Tran_df), ignore.case = TRUE) #40931
oshv1_IAP_non_sig_info_resoshv1Tran_df12 <- resoshv1Tran_df[40931,]
oshv1_IAP_non_sig_info_resoshv1Tran_df13 <- grep("EKC42443", rownames(resoshv1Tran_df), ignore.case = TRUE) #43642
oshv1_IAP_non_sig_info_resoshv1Tran_df13 <- resoshv1Tran_df[43642,]
oshv1_IAP_non_sig_info_resoshv1Tran_df14 <- grep("EKC18368", rownames(resoshv1Tran_df), ignore.case = TRUE) #1136
oshv1_IAP_non_sig_info_resoshv1Tran_df14 <- resoshv1Tran_df[1136,]
oshv1_IAP_non_sig_info_resoshv1Tran_df15 <- grep("EKC25955", rownames(resoshv1Tran_df), ignore.case = TRUE) #7418
oshv1_IAP_non_sig_info_resoshv1Tran_df15 <- resoshv1Tran_df[7418,]
oshv1_IAP_non_sig_info_resoshv1Tran_df16 <- grep("EKC41180", rownames(resoshv1Tran_df), ignore.case = TRUE) #8987
oshv1_IAP_non_sig_info_resoshv1Tran_df16 <- resoshv1Tran_df[8987,]
oshv1_IAP_non_sig_info_resoshv1Tran_df17 <- grep("EKC41181", rownames(resoshv1Tran_df), ignore.case = TRUE) #8988
oshv1_IAP_non_sig_info_resoshv1Tran_df17 <- resoshv1Tran_df[8988,]
oshv1_IAP_non_sig_info_resoshv1Tran_df18 <- grep("EKC33184", rownames(resoshv1Tran_df), ignore.case = TRUE) #11127
oshv1_IAP_non_sig_info_resoshv1Tran_df18 <- resoshv1Tran_df[11127,]
oshv1_IAP_non_sig_info_resoshv1Tran_df19 <- grep("EKC29824", rownames(resoshv1Tran_df), ignore.case = TRUE) #25305
oshv1_IAP_non_sig_info_resoshv1Tran_df19 <- resoshv1Tran_df[25305,]
oshv1_IAP_non_sig_info_resoshv1Tran_df20 <- grep("EKC17690", rownames(resoshv1Tran_df), ignore.case = TRUE) #25645
oshv1_IAP_non_sig_info_resoshv1Tran_df20 <- resoshv1Tran_df[25645,]
oshv1_IAP_non_sig_info_resoshv1Tran_df21 <- grep("EKC34720", rownames(resoshv1Tran_df), ignore.case = TRUE) #28700
oshv1_IAP_non_sig_info_resoshv1Tran_df21 <- resoshv1Tran_df[28700,]
oshv1_IAP_non_sig_info_resoshv1Tran_df22 <- grep("EKC34022", rownames(resoshv1Tran_df), ignore.case = TRUE) #32211
oshv1_IAP_non_sig_info_resoshv1Tran_df22 <- resoshv1Tran_df[32211,]
oshv1_IAP_non_sig_info_resoshv1Tran_df23 <- grep("EKC26454", rownames(resoshv1Tran_df), ignore.case = TRUE) #34851
oshv1_IAP_non_sig_info_resoshv1Tran_df23 <- resoshv1Tran_df[34851,]
oshv1_IAP_non_sig_info_resoshv1Tran_df24 <- grep("EKC26950", rownames(resoshv1Tran_df), ignore.case = TRUE) #35402
oshv1_IAP_non_sig_info_resoshv1Tran_df24 <- resoshv1Tran_df[35402,]
oshv1_IAP_non_sig_info_resoshv1Tran_df25 <- grep("EKC42442", rownames(resoshv1Tran_df), ignore.case = TRUE) #35792
oshv1_IAP_non_sig_info_resoshv1Tran_df25 <- resoshv1Tran_df[35792,]



oshv1_IAP_non_sig_combined_EXP <- rbind(oshv1_IAP_non_sig_info_resoshv1Tran_df2, 
                                        oshv1_IAP_non_sig_info_resoshv1Tran_df3, oshv1_IAP_non_sig_info_resoshv1Tran_df4,
                                        oshv1_IAP_non_sig_info_resoshv1Tran_df5, oshv1_IAP_non_sig_info_resoshv1Tran_df6, 
                                        oshv1_IAP_non_sig_info_resoshv1Tran_df7, oshv1_IAP_non_sig_info_resoshv1Tran_df8,
                                        oshv1_IAP_non_sig_info_resoshv1Tran_df9, oshv1_IAP_non_sig_info_resoshv1Tran_df10,
                                        oshv1_IAP_non_sig_info_resoshv1Tran_df11, oshv1_IAP_non_sig_info_resoshv1Tran_df12,
                                        oshv1_IAP_non_sig_info_resoshv1Tran_df13, oshv1_IAP_non_sig_info_resoshv1Tran_df14,
                                        oshv1_IAP_non_sig_info_resoshv1Tran_df15, oshv1_IAP_non_sig_info_resoshv1Tran_df16,
                                        oshv1_IAP_non_sig_info_resoshv1Tran_df17, oshv1_IAP_non_sig_info_resoshv1Tran_df18,
                                        oshv1_IAP_non_sig_info_resoshv1Tran_df19, oshv1_IAP_non_sig_info_resoshv1Tran_df20,
                                        oshv1_IAP_non_sig_info_resoshv1Tran_df21, oshv1_IAP_non_sig_info_resoshv1Tran_df22,
                                        oshv1_IAP_non_sig_info_resoshv1Tran_df23, oshv1_IAP_non_sig_info_resoshv1Tran_df24,
                                        oshv1_IAP_non_sig_info_resoshv1Tran_df25)
oshv1_IAP_non_sig_combined_FULL <- cbind(oshv1_IAP_non_sig_combined_EXP, subset_oshv1_IAPs_non_sig_info)


#Non Sig GIMAP
oshv1_GIMAP_non_sig_info$Ensembl_Genome_Transcript
# EKC39736 EKC40465 EKC40820 EKC32489 EKC36405 EKC39748 EKC27363 EKC31739 EKC35292 EKC29604 EKC41613 EKC42724 EKC38639
subset_GIMAP_non_sig_info <- oshv1_GIMAP_non_sig_info[,c(4,8,9)]
oshv1_GIMAP_non_sig_info_resoshv1Tran_df <- grep("EKC39736", rownames(resoshv1Tran_df), ignore.case = TRUE) #538
oshv1_GIMAP_non_sig_info_resoshv1Tran_df <- resoshv1Tran_df[538,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_df2 <- grep("EKC40465", rownames(resoshv1Tran_df), ignore.case = TRUE) #1984
oshv1_GIMAP_non_sig_info_resoshv1Tran_df2 <- resoshv1Tran_df[1984,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_df3 <- grep("EKC40820", rownames(resoshv1Tran_df), ignore.case = TRUE) #2258
oshv1_GIMAP_non_sig_info_resoshv1Tran_df3 <- resoshv1Tran_df[2258,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_df4 <- grep("EKC32489", rownames(resoshv1Tran_df), ignore.case = TRUE) #8113
oshv1_GIMAP_non_sig_info_resoshv1Tran_df4 <- resoshv1Tran_df[8113,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_df5 <- grep("EKC36405", rownames(resoshv1Tran_df), ignore.case = TRUE) #9883
oshv1_GIMAP_non_sig_info_resoshv1Tran_df5 <- resoshv1Tran_df[9883,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_df6 <- grep("EKC39748", rownames(resoshv1Tran_df), ignore.case = TRUE) #24210
oshv1_GIMAP_non_sig_info_resoshv1Tran_df6 <- resoshv1Tran_df[24210,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_df7 <- grep("EKC27363", rownames(resoshv1Tran_df), ignore.case = TRUE) #24691
oshv1_GIMAP_non_sig_info_resoshv1Tran_df7 <- resoshv1Tran_df[24691,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_df8 <- grep("EKC31739", rownames(resoshv1Tran_df), ignore.case = TRUE) #25182
oshv1_GIMAP_non_sig_info_resoshv1Tran_df8 <- resoshv1Tran_df[25182,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_df9 <- grep("EKC35292", rownames(resoshv1Tran_df), ignore.case = TRUE) #28242
oshv1_GIMAP_non_sig_info_resoshv1Tran_df9 <- resoshv1Tran_df[28242,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_df10 <- grep("EKC29604", rownames(resoshv1Tran_df), ignore.case = TRUE) #33417
oshv1_GIMAP_non_sig_info_resoshv1Tran_df10 <- resoshv1Tran_df[33417,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_df11 <- grep("EKC41613", rownames(resoshv1Tran_df), ignore.case = TRUE) #35144
oshv1_GIMAP_non_sig_info_resoshv1Tran_df11 <- resoshv1Tran_df[35144,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_df12 <- grep("EKC42724", rownames(resoshv1Tran_df), ignore.case = TRUE) #35519
oshv1_GIMAP_non_sig_info_resoshv1Tran_df12 <- resoshv1Tran_df[35519,]
oshv1_GIMAP_non_sig_info_resoshv1Tran_df13 <- grep("EKC38639", rownames(resoshv1Tran_df), ignore.case = TRUE) #43959
oshv1_GIMAP_non_sig_info_resoshv1Tran_df13 <- resoshv1Tran_df[43959,]

oshv1_GIMAP_non_sig_combined_EXP <- rbind(oshv1_GIMAP_non_sig_info_resoshv1Tran_df, oshv1_GIMAP_non_sig_info_resoshv1Tran_df2,
                                          oshv1_GIMAP_non_sig_info_resoshv1Tran_df3, oshv1_GIMAP_non_sig_info_resoshv1Tran_df4,
                                          oshv1_GIMAP_non_sig_info_resoshv1Tran_df5, oshv1_GIMAP_non_sig_info_resoshv1Tran_df6,
                                          oshv1_GIMAP_non_sig_info_resoshv1Tran_df7, oshv1_GIMAP_non_sig_info_resoshv1Tran_df8,
                                          oshv1_GIMAP_non_sig_info_resoshv1Tran_df9, oshv1_GIMAP_non_sig_info_resoshv1Tran_df10,
                                          oshv1_GIMAP_non_sig_info_resoshv1Tran_df11,oshv1_GIMAP_non_sig_info_resoshv1Tran_df12,
                                          oshv1_GIMAP_non_sig_info_resoshv1Tran_df13)

oshv1_GIMAP_non_sig_combined_FULL <- cbind(oshv1_GIMAP_non_sig_combined_EXP, subset_GIMAP_non_sig_info)

####Combine all Oshv1 IAP and GIMAP data #####
#to rbind the names of all the columns need to be the same, change first column name
oshv1_GIMAP_IAP_combined_FULL <- rbind(oshv1_GIMAP_non_sig_combined_FULL, oshv1_IAP_non_sig_combined_FULL,
                                       oshv1_IAP_SIG_combined_FULL, oshv1_GIMAP_SIG_combined_FULL)

####OsHV1 Gene Set Enrichment Analysis topGO ####

#Matching the background set
#Get average expressions 
oshv1Tran_BaseMean <- as.matrix(resoshv1Tran_df[, "baseMean", drop=F])
oshv1Tran_backG <- genefinder(oshv1_6_BaseMean, anSig$ensembl_gene_id, 10, method= "manhattan")


####Bacterial Challenge Differential Gene Expression Analysis ####
#Gram negative SRA's (including control): (SRR796597, SRR796596, SRR796595, SRR796594, SRR796593, SRR796592, SRR796589)
# Gram positive (including control): (SRR796598,  SRR796589)

#Subset Bacterial Challenge Data
BacTranCountData <- as.matrix(FullTranscriptCountData[ , c(124:131)])
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
sum(resBacTran$padj < 0.1, na.rm=TRUE) #2322
#going to subset for 0.05 after the p value correction!
sum(resBacTran_05$padj < 0.05, na.rm=TRUE) #1889

#metadata on meaning of the columns
mcols(resBacTran, use.names = TRUE)
#shows treatment vs. control with baseMean as the mean of normalized counts for all samples

####p-value correction#### adapted from "Differential expression analysis of RNA-Seq data using DESeq2" Klaus 2014
#First Visualize histograms
#histogram of P- values to visualize any "hills" or "U shape"
# hill means variance of the null distribution too high, U shape means variance assumed too low
hist(resBacTran$pvalue, breaks = 20, col = "grey") #hill

#remove filtered out genes by independent filtering, they have NA adj. pvals
resBacTran_df <- resBacTran[ !is.na(resBacTran$padj), ]

#remove genes with NA pvals (outliers)
resBacTran_df <- resBacTran_df[ !is.na(resBacTran_df$pvalue), ]

#remove adjsuted pvalues, since we add the fdrtool results later on (based on the correct p-values)
resBacTran_df <- resBacTran_df[, -which(names(resBacTran_df) == "padj")]

#use z-scores as input to FDRtool to re-estimate the p-value
FDR.resBacTran_df <- fdrtool(resBacTran_df$stat, statistic= "normal", plot = T)

#add values to the results data frame, also ad new BH- adjusted p-values
resBacTran_df[,"padj"] <- p.adjust(FDR.resBacTran_df$pval, method = "BH")

#replot corrected p-values
hist(FDR.resBacTran_df$pval, col = "royalblue4",
     main = "Correct null model Bacterial Challenge Transcript Counts", xlab = "CORRECTED p-values")

#Check how many genes have BH adjusted p values of less than 0.05 after P-value correction?
sum( resBacTran_df$padj < 0.05, na.rm=TRUE ) #1398 (before p-value correction it was )

#Subset the results table to the differentially expressed genes under FDR 0.01, order the Log2FC table first by strongest down regulation
resBacTran_dfSig <- resBacTran_df[ which(resBacTran_df$padj < 0.05 ), ]
head( resBacTran_dfSig[ order( resBacTran_dfSig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resBacTran_dfSig[ order( resBacTran_dfSig$log2FoldChange ), ] ) #tail for strongest up regulation
resBacTran_df_non_Sig <- resBacTran_df[ which(resBacTran_df$padj > 0.05 ), ]
summary(resBacTran_df_non_Sig)

####Visualize Results with Diagnostic Plots####
#MA plot, useful overview for experiment with two-group comparison. Plots log2FC over mean of normalized counts
#genes with adjusted p value 
plotMA(resBacTran_df)
plotMA(resBacTran_dfSig)
plotMA(resBacTran_df_non_Sig)

#Export Results to CSV
write.csv( as.data.frame(resBacTran_df), file="resBacTran_df.csv")
write.csv( as.data.frame(resBacTran_dfSig), file="resBacTran_dfSig.csv")
write.csv( as.data.frame(resBacTran_df_non_Sig), file="resBacTran_df_non_Sig.csv")

####Subset files for only those that have Transcript IDs####
#Extract gene titles of all the significantly differentially expressed genes
Bac_withID_subset_resBacTran_dfSig <- 
  resBacTran_dfSig[grep("transcript:", rownames(resBacTran_dfSig)), ]
Bac_withID_subset_resBacTran_dfSig #193 rows
BacTranscriptIDdf <- as.data.frame(rownames(Bac_withID_subset_resBacTran_dfSig))
head(BacTranscriptIDdf)
BacTranscriptIDdf= transform(BacTranscriptIDdf, 
                          ID = colsplit(rownames(Bac_withID_subset_resBacTran_dfSig), 
                                        split = "\\:", names = c('transcript:', 'EKC36484')))
#this line splitsup the word transcript and the transcript name

BacTranscriptIDstring <- toString(BacTranscriptIDdf[,3], sep=',')
write(BacTranscriptIDstring, "BacTranscriptIDstring_sig", sep = ",")
#write this to a file and then perform look up on the UniProt website

#Extract gene titles from non Sig genes
Bac_withID_subset_resBacTran_df_non_Sig <- 
  resBacTran_df_non_Sig[grep("transcript:", rownames(resBacTran_df_non_Sig)), ]
head(Bac_withID_subset_resBacTran_df_non_Sig)
BacTranscriptIDdf_nonSig <- as.data.frame(rownames(Bac_withID_subset_resBacTran_df_non_Sig))
head(BacTranscriptIDdf_nonSig)
BacTranscriptIDdf_nonSig= transform(BacTranscriptIDdf_nonSig, 
                                 ID = colsplit(rownames(Bac_withID_subset_resBacTran_df_non_Sig), 
                                               split = "\\:", names = c('transcript:', 'EKC36261')))
BacTranscriptIDstring_nonSig <- toString(BacTranscriptIDdf_nonSig[,3], sep=',')
BacTranscriptIDstring_nonSig
write(BacTranscriptIDstring_nonSig, "BacTranscriptIDstring_nonSig", sep = ",")
#write this to a file and then perform look up on the UniProt website

#Add quotes around this 
#transcriptIDparen <- sapply(strsplit(transcriptIDstring, '[, ]+'), function(x) toString(dQuote(x)))
#str(transcriptIDparen)

####For situations only looking at the gene name of a few genes, can look it up locally ####
#Load C. gigas UniProt ID information for C. gigas
source("https://bioconductor.org/biocLite.R")
biocLite("UniProt.ws")
library(UniProt.ws)
browseVignettes("UniProt.ws")
#Taxonomy ID for C. virginica: 29159
availableUniprotSpecies(pattern = "gigas")
CgigasUp <- UniProt.ws(29159)
CgigasUp
keytypes(CgigasUp)
columns(CgigasUp)
keytypes(CgigasUp)
showMethods("keys")
#structure: res <- select(UniProtName, keys, columns, keytype)
#select(up, keys=c("P31946","P62258"), columns=c("PDB","SEQUENCE"), keytype="ENSEMBL_GENOMES TRANSCRIPT)

#Create object to look up UNIPROTKB entries and retreive other gene info
columns <- c("ENSEMBL_PROTEIN","HGNC","SEQUENCE", "GENES","ENTREZ_GENE", "ID", "GO", "GO-ID", "KEGG")
keytype <- "ENSEMBL_GENOMES TRANSCRIPT"
#keys <- transcriptIDdf[,3] #set key to be  vector of all of your transcriptIDs from previous
keys <- transcriptIDparen[1:50] #test to see how this works
res <- select(CgigasUp, keys, columns, keytype)
res

####Otherwise, Uploaded Sig Differentially Expressed transcripts from the UniProt.ws website ####
Bac_transcriptIDs_UniProt_SIG <- read.csv("BacTran_resBacTran_df_Sig_transcriptIDstring.csv", header=TRUE)
#make sure to upload as csv version so that the protein names load correctly
head(Bac_transcriptIDs_UniProt_SIG)

####Extract GIMAP/IAN proteins and CgIAPs from Significantly Differentially Expressed Genes####
#Significantly differentially Expressed IAPs
Bac_transcriptIDs_UniProt_SIG_ProtNames <- Bac_transcriptIDs_UniProt_SIG$Protein.names
Bac_IAPs_SIG <- grepl("IAP", Bac_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", Bac_IAPs_SIG) #32...but 32 doesn't have anything 
Bac_apoptosis_SIG <- grepl("apoptosis", Bac_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", Bac_apoptosis_SIG) #126 
Bac_inhibitor_SIG <- grepl("inhibitor", Bac_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", Bac_inhibitor_SIG) #84, 126 , 84 is the BAX inhibitor

Bac_IAP_SIG_info <- Bac_transcriptIDs_UniProt_SIG[126,]


#Significant GIMAP Genes
Bac_GTP_SIG <- grepl("GTP", Bac_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
grep("TRUE", Bac_GTP_SIG) #67 134 157 167, none are GIMAP
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
Bac_transcriptIDs_UniProt_non_Sig

####Extract GIMAP and IAP genes from NON significant genes, for comparison ####
#Expressed IAPs non sig
#use grepl to find text strings 
Bac_transcriptIDs_UniProt_non_Sig_ProtNames <- Bac_transcriptIDs_UniProt_non_Sig$Protein.names
Bac_IAPs_non_Sig <- grepl("IAP", Bac_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", Bac_IAPs_non_Sig) 
#22  367 2484 2528 2707 4848 4850
Bac_IAP_non_Sig_info <- Bac_transcriptIDs_UniProt_non_Sig[c(22,  367, 2484, 2528, 2707, 4848, 4850),]

Bac_apoptosis_non_Sig <- grepl("apoptosis", Bac_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", Bac_apoptosis_non_Sig) #270 1887 1950 2099 2177 2868 3292 3393 3650 4324 4655
Bac_apoptosis_non_Sig_info <- Bac_transcriptIDs_UniProt_non_Sig[c(270,1887, 1950, 2099, 2177, 2868, 3292, 3393, 3650, 4324, 4655),]
#Actual Apoptosis Inhibitors: 1887, 1950, 2099, 2177, 2868, 3292, 3650 
Bac_inhibitor_non_Sig <- grepl("inhibitor", Bac_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", Bac_inhibitor_non_Sig) #190  372  568  656  816 1029 1134 1246 1600 1813 1887 1950 2099 2177 2307 2519 2868 3154 3174
# 3210 3292 3358 3378 3650 3694 3705 3915 3921 4057 4067 4348 4487 4570 4924 5341
Bac_inhibitor_non_Sig_info <- Bac_transcriptIDs_UniProt_non_Sig[c(190,  372,  568,  656,  816, 1029, 1134, 
1246, 1600, 1813, 1887, 1950, 2099, 2177, 2307, 2519, 2868, 3154, 3174, 3210, 3292, 3358, 3378, 3650, 3694,
3705, 3915, 3921 ,4057, 4067, 4348, 4487, 4570, 4924, 5341),] 
#Actual Apoptosis Inhibitors: 1887, 1950, 2099, 2177, 2868, 3292, 3650, same as above!

Bac_IAPs_non_sig_info <- Bac_transcriptIDs_UniProt_non_Sig[c(22,  367, 2484, 2528, 2707, 4848, 4850,
                                                             1887, 1950, 2099, 2177, 2868, 3292, 3650),]

#Expressed GIMAP Genes, non significant
Bac_GTP_non_Sig <- grepl("GTP", Bac_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", Bac_GTP_non_Sig)
#141  310  318  600  841  924  939  985 1009 1325 1338 1535 1596 1626 1685 1702 1704 1879 2020
#2131 2377 2429 2669 2805 3132 3221 3461 3580 4049 4169 4197 4199 4273 4353 4383 4599 4710 4767
#4816 4969 5009 5023 5201 5474
Bac_GTP_non_Sig_info <- Bac_transcriptIDs_UniProt_non_Sig[c(141,  310,  318,  600,  841,  924 , 939,  985, 1009, 1325 ,
                                                            1338, 1535, 1596, 1626, 1685, 1702, 1704, 1879, 2020,
                                                            2131, 2377, 2429, 2669, 2805, 3132, 3221, 3461, 3580,
                                                            4049, 4169, 4197, 4199, 4273, 4353, 4383, 4599, 4710, 4767,
                                                            4816, 4969, 5009, 5023, 5201, 5474),]
#Actual GIMAPS: 141, 318, 1338, 1685, 1879, 3132, 4199, 4969, 5474
Bac_IAN_non_Sig <- grepl("IAN", Bac_transcriptIDs_UniProt_non_Sig$Protein.names) 
grep("TRUE", Bac_IAN_non_Sig) #0
Bac_immune_non_Sig <- grepl("immune", Bac_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", Bac_immune_non_Sig) #0 TRUE
Bac_AIG_non_Sig <- grepl("AIG", Bac_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", Bac_AIG_non_Sig) #0 TRUE
Bac_IMAP_non_sig <- grepl("IMAP",Bac_transcriptIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE)
grep("TRUE", Bac_IMAP_non_sig)
#GIMAPS: 141  318 1338 1685 1879 3132 4199 4969 5474

Bac_GIMAP_non_sig_info <- Bac_transcriptIDs_UniProt_non_Sig[c(141,  318, 1338, 1685, 1879, 3132, 4199, 4969, 5474),]


####Link GIMAP and IAN genes, significant and non significant, with Expression Values####
#Sig IAPS
subset_Bac_IAPs_SIG_info <- Bac_IAP_SIG_info[,c(4,8,9)]
Bac_IAP_SIG_info_resBacTran_dfSig <- grep("EKC41180", rownames(resBacTran_dfSig), ignore.case = TRUE) # 848
Bac_IAP_SIG_info_resBacTran_dfSig <- resBacTran_dfSig[848,] 
Bac_IAP_SIG_combined_FULL <- cbind(Bac_IAP_SIG_info_resBacTran_dfSig,subset_Bac_IAPs_SIG_info)

#Sig GIMAPs
  #NO SIG GIMAPs

#Non Sig IAP
Bac_IAPs_non_sig_info$Ensembl_Genome_Transcript 
#EKC24074 EKC38724 EKC42441 EKC25493 EKC42449 EKC18369 EKC20239 EKC17690 EKC34022 EKC26454 EKC42442
#EKC18368 EKC41181 EKC26950
subset_Bac_IAPs_non_sig_info <- Bac_IAPs_non_sig_info[,c(4,8,9)]
Bac_IAP_non_sig_info_resBacTran_df <- grep("EKC24074", rownames(resBacTran_df_non_Sig), ignore.case = TRUE) #48
Bac_IAP_non_sig_info_resBacTran_df <- resBacTran_df_non_Sig[48,]
Bac_IAP_non_sig_info_resBacTran_df2 <- grep("EKC38724", rownames(resBacTran_df_non_Sig), ignore.case = TRUE) #2294
Bac_IAP_non_sig_info_resBacTran_df2 <- resBacTran_df_non_Sig[2294,]
Bac_IAP_non_sig_info_resBacTran_df3 <- grep("EKC42441", rownames(resBacTran_df_non_Sig), ignore.case = TRUE) #14607
Bac_IAP_non_sig_info_resBacTran_df3 <- resBacTran_df_non_Sig[14607,]
Bac_IAP_non_sig_info_resBacTran_df4 <- grep("EKC25493", rownames(resBacTran_df_non_Sig), ignore.case = TRUE) #14856
Bac_IAP_non_sig_info_resBacTran_df4 <- resBacTran_df_non_Sig[14856,]
Bac_IAP_non_sig_info_resBacTran_df5 <- grep("EKC42449", rownames(resBacTran_df_non_Sig), ignore.case = TRUE) #15901
Bac_IAP_non_sig_info_resBacTran_df5 <- resBacTran_df_non_Sig[15901,]
Bac_IAP_non_sig_info_resBacTran_df6 <- grep("EKC18369", rownames(resBacTran_df_non_Sig), ignore.case = TRUE) #28537
Bac_IAP_non_sig_info_resBacTran_df6 <- resBacTran_df_non_Sig[28537,]
Bac_IAP_non_sig_info_resBacTran_df7 <- grep("EKC20239", rownames(resBacTran_df_non_Sig), ignore.case = TRUE) #28555
Bac_IAP_non_sig_info_resBacTran_df7 <- resBacTran_df_non_Sig[28555,]
Bac_IAP_non_sig_info_resBacTran_df8 <- grep("EKC17690", rownames(resBacTran_df_non_Sig), ignore.case = TRUE) #11190
Bac_IAP_non_sig_info_resBacTran_df8 <- resBacTran_df_non_Sig[11190,]
Bac_IAP_non_sig_info_resBacTran_df9 <- grep("EKC34022", rownames(resBacTran_df_non_Sig), ignore.case = TRUE) #11641
Bac_IAP_non_sig_info_resBacTran_df9 <- resBacTran_df_non_Sig[11641,]
Bac_IAP_non_sig_info_resBacTran_df10 <- grep("EKC26454", rownames(resBacTran_df_non_Sig), ignore.case = TRUE) #12597
Bac_IAP_non_sig_info_resBacTran_df10 <- resBacTran_df_non_Sig[12597,]
Bac_IAP_non_sig_info_resBacTran_df11 <- grep("EKC42442", rownames(resBacTran_df_non_Sig), ignore.case = TRUE) #12949
Bac_IAP_non_sig_info_resBacTran_df11 <- resBacTran_df_non_Sig[12949,]
Bac_IAP_non_sig_info_resBacTran_df12 <- grep("EKC18368", rownames(resBacTran_df_non_Sig), ignore.case = TRUE) #16861
Bac_IAP_non_sig_info_resBacTran_df12 <- resBacTran_df_non_Sig[16861,]
Bac_IAP_non_sig_info_resBacTran_df13 <- grep("EKC41181", rownames(resBacTran_df_non_Sig), ignore.case = TRUE) #19258
Bac_IAP_non_sig_info_resBacTran_df13 <- resBacTran_df_non_Sig[19258,]
Bac_IAP_non_sig_info_resBacTran_df14 <- grep("EKC26950", rownames(resBacTran_df_non_Sig), ignore.case = TRUE) #21340
Bac_IAP_non_sig_info_resBacTran_df14 <- resBacTran_df_non_Sig[21340,]
Bac_IAP_non_sig_combined_EXP <- rbind(Bac_IAP_non_sig_info_resBacTran_df, Bac_IAP_non_sig_info_resBacTran_df2,
                                      Bac_IAP_non_sig_info_resBacTran_df3, Bac_IAP_non_sig_info_resBacTran_df4,
                                      Bac_IAP_non_sig_info_resBacTran_df5,Bac_IAP_non_sig_info_resBacTran_df6,
                                      Bac_IAP_non_sig_info_resBacTran_df7,Bac_IAP_non_sig_info_resBacTran_df8,
                                      Bac_IAP_non_sig_info_resBacTran_df9, Bac_IAP_non_sig_info_resBacTran_df10,
                                      Bac_IAP_non_sig_info_resBacTran_df11,Bac_IAP_non_sig_info_resBacTran_df12,
                                      Bac_IAP_non_sig_info_resBacTran_df13,Bac_IAP_non_sig_info_resBacTran_df14)
Bac_IAP_non_sig_combined_FULL <- cbind(Bac_IAP_non_sig_combined_EXP,subset_Bac_IAPs_non_sig_info)

#Non Sig GIMAP
Bac_GIMAP_non_sig_info$Ensembl_Genome_Transcript
#EKC40820 EKC41832 EKC41613 EKC35292 EKC36405 EKC30713 EKC40465 EKC42724 EKC38639
subset_Bac_GIMAP_non_sig_info <- Bac_GIMAP_non_sig_info[,c(4,8,9)]
Bac_GIMAP_non_sig_info_resBacTran_df <- grep("EKC40820", rownames(resBacTran_df_non_Sig), ignore.case = TRUE) #825
Bac_GIMAP_non_sig_info_resBacTran_df <- resBacTran_df_non_Sig[825,]
Bac_GIMAP_non_sig_info_resBacTran_df2 <- grep("EKC41832", rownames(resBacTran_df_non_Sig), ignore.case = TRUE) 
Bac_GIMAP_non_sig_info_resBacTran_df2 <- resBacTran_df_non_Sig[1896,]
Bac_GIMAP_non_sig_info_resBacTran_df3 <- grep("EKC41613", rownames(resBacTran_df_non_Sig), ignore.case = TRUE)
Bac_GIMAP_non_sig_info_resBacTran_df3 <- resBacTran_df_non_Sig[7998,]
Bac_GIMAP_non_sig_info_resBacTran_df4 <- grep("EKC35292", rownames(resBacTran_df_non_Sig), ignore.case = TRUE) 
Bac_GIMAP_non_sig_info_resBacTran_df4 <- resBacTran_df_non_Sig[10221,]
Bac_GIMAP_non_sig_info_resBacTran_df5 <- grep("EKC36405", rownames(resBacTran_df_non_Sig), ignore.case = TRUE) 
Bac_GIMAP_non_sig_info_resBacTran_df5 <- resBacTran_df_non_Sig[11154,]
Bac_GIMAP_non_sig_info_resBacTran_df6 <- grep("EKC30713", rownames(resBacTran_df_non_Sig), ignore.case = TRUE) 
Bac_GIMAP_non_sig_info_resBacTran_df6 <- resBacTran_df_non_Sig[18385,]
Bac_GIMAP_non_sig_info_resBacTran_df7 <- grep("EKC40465", rownames(resBacTran_df_non_Sig), ignore.case = TRUE) 
Bac_GIMAP_non_sig_info_resBacTran_df7 <- resBacTran_df_non_Sig[24720,]
Bac_GIMAP_non_sig_info_resBacTran_df8 <- grep("EKC42724", rownames(resBacTran_df_non_Sig), ignore.case = TRUE) 
Bac_GIMAP_non_sig_info_resBacTran_df8 <- resBacTran_df_non_Sig[29415,]
Bac_GIMAP_non_sig_info_resBacTran_df9 <- grep("EKC38639", rownames(resBacTran_df_non_Sig), ignore.case = TRUE) 
Bac_GIMAP_non_sig_info_resBacTran_df9 <- resBacTran_df_non_Sig[32500,]

Bac_GIMAP_non_sig_combined_EXP <- rbind(Bac_GIMAP_non_sig_info_resBacTran_df, Bac_GIMAP_non_sig_info_resBacTran_df2,
                                        Bac_GIMAP_non_sig_info_resBacTran_df3, Bac_GIMAP_non_sig_info_resBacTran_df4,
                                        Bac_GIMAP_non_sig_info_resBacTran_df5, Bac_GIMAP_non_sig_info_resBacTran_df6,
                                        Bac_GIMAP_non_sig_info_resBacTran_df7,Bac_GIMAP_non_sig_info_resBacTran_df8,
                                        Bac_GIMAP_non_sig_info_resBacTran_df9)

Bac_GIMAP_non_sig_combined_FULL <- cbind(Bac_GIMAP_non_sig_combined_EXP,subset_Bac_GIMAP_non_sig_info)

####Combine all Bac IAP and GIMAP data #####
#to rbind the names of all the columns need to be the same, change first column name
#Fix column name error for column 9 in SIG set
colnames(Bac_IAP_SIG_combined_FULL)[9] <- "Ensembl_Genome_Transcript"
Bac_GIMAP_IAP_combined_FULL <- rbind(Bac_IAP_SIG_combined_FULL, Bac_IAP_non_sig_combined_FULL,
                                     Bac_GIMAP_non_sig_combined_FULL)

####Bacterial Challenge Gene Set Enrichment Analysis ####
#Matching the background set
#Get average expressions 
#Bac_BaseMean <- as.matrix(resoshv1_6_df[, "baseMean", drop=F])
#Bac_backG <- genefinder(oshv1_6_BaseMean, anSig$ensembl_gene_id, 10, method= "manhattan")


#### COMPILE GIMAP AND IAP DATA FROM BOTH OSHV1 AND BAC TRIAL ####
Bac_GIMAP_IAP_combined_FULL["Challenge"] <- "Bacteria"
#Add column indicating trial
oshv1_GIMAP_IAP_combined_FULL["Challenge"] <- "OsHV1"
colnames(oshv1_GIMAP_IAP_combined_FULL)
colnames(Bac_GIMAP_IAP_combined_FULL)

COMBINED_GIMAP_IAP <- rbind(Bac_GIMAP_IAP_combined_FULL, oshv1_GIMAP_IAP_combined_FULL)

####Plotting the combined Significance values ####
COMBINED_GIMAP_IAP_cols <- COMBINED_GIMAP_IAP[,c(2,7,10)]

#Download data to simplify Protein Names for Viewing 
write.csv( as.data.frame(COMBINED_GIMAP_IAP_cols), file="COMBINED_GIMAP_IAP_cols.csv")
COMBINED_GIMAP_IAP_cols <- read.csv("COMBINED_GIMAP_IAP_cols_edited.csv", header = TRUE)
COMBINED_GIMAP_IAP_cols_Bac <- COMBINED_GIMAP_IAP_cols[1:24,]
COMBINED_GIMAP_IAP_cols_Oshv1 <- COMBINED_GIMAP_IAP_cols[25:66,] 

#Plot the combined data set altogether for comparison
plot(log2FoldChange~factor(Protein.names), COMBINED_GIMAP_IAP_cols, las=2,
     xlab="GIMAP and IAP Transcripts", main="Differentially Expressed GIMAP and IAP Transcripts")

#Plot GIMAP between OsHV1 and Bac

#Plot IAP between OsHV1 and Bac

#references: 
#https://github.com/sr320/LabDocs/tree/master/code/DESeq
#StringTie manual : http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual