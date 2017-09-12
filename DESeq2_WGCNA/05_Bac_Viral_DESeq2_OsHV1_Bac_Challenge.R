#05_Bac_viral_DESeq2_OsHV1_Bac_challenge

#This script takes as input the output Bac_Viral_gff3HIT_subset_transcript_count_matrix.csv data prepared from prepDE.py and performs
#differential Gene expression analysis, and subsets out isoforms of GIMAPs and CgIAPs and graphs their
#relative abundance.

#call the DESeq2 library 
#source("https://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade") 
#biocLite("DESeq2")
library("DESeq2")
#install.packages("fdrtool")
library(fdrtool)
#source("https://bioconductor.org/biocLite.R")
#biocLite("genefilter")
library(genefilter)
# Construct Bac_Viral_PHENO_DATA.csv that contains SRA run information, such as which contrast, tissue, etc.

####DEG Analysis with TRANSCRIPT Count Matrix ####
#load transcript count matrix and labels
#Bac_Viral_PHENO_DATA.csv file contains metadata on the count table's samples
###Make sure excel PHENODATA is in the same order or these commands will change data to be wrong!!!!###

TranscriptCountData <- as.matrix(read.csv("Bac_Viral_gff3HIT_subset_transcript_count_matrix.csv", row.names="transcript_id"))
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
summary(resoshv1Tran)
#sum(resoshv1Tran$padj < 0.1, na.rm=TRUE) #3501
resoshv1Tran_05 <- results(ddsOshv1Tran,alpha=0.05)
summary(resoshv1Tran_05)
sum(resoshv1Tran_05$padj < 0.05, na.rm=TRUE) #2590

#metadata on meaning of the columns
mcols(resoshv1Tran_05, use.names = TRUE)
#shows treatment vs. control with baseMean as the mean of normalized counts for all samples

####p-value correction#### adapted from "Differential expression analysis of RNA-Seq data using DESeq2" Klaus 2014
#First Visualize histograms
#histogram of P- values to visualize any "hills" or "U shape"
# hill means variance of the null distribution too high, U shape means variance assumed too low
hist(resoshv1Tran_05$pvalue, breaks = 20, col = "grey") #hill

#remove filtered out genes by independent filtering, they have NA adj. pvals
resoshv1Tran_05_df <- resoshv1Tran_05[ !is.na(resoshv1Tran_05$padj), ]

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

#Check how many genes have BH adjusted p values of less than 0.01 after P-value correction?
sum( resoshv1Tran_05_df$padj < 0.05, na.rm=TRUE ) #1715 (before p-value correction it was )... now its 2

#Subset the results table to the differentially expressed genes under FDR 0.01, order the Log2FC table first by strongest down regulation
resoshv1Tran_05_dfSig <- resoshv1Tran_05_df[ which(resoshv1Tran_05_df$padj < 0.05 ), ]
head( resoshv1Tran_05_dfSig[ order( resoshv1Tran_05_dfSig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resoshv1Tran_05_dfSig[ order( resoshv1Tran_05_dfSig$log2FoldChange ), ] ) #tail for strongest up regulation

####Visualize Results with Diagnostic Plots####
#MA plot, useful overview for experiment with two-group comparison. Plots log2FC over mean of normalized counts
#genes with adjusted p value 
plotMA(resoshv1Tran_05_df)
plotMA(resoshv1Tran_05_dfSig)

#Export Results to CSV
write.csv( as.data.frame(resoshv1Tran_05_df), file="resoshv1Tran_05_df.csv")
write.csv( as.data.frame(resoshv1Tran_05_dfSig), file="resoshv1Tran_05_dfSig.csv")

####Use Bash script fetchEnsembl_ID.sh to get the Ensembl IDs ####
#Extract gene titles of all the significantly differentially expressed genes

####OsHV1 Gene Set Enrichment Analysis ####
#Matching the background set
#Get average expressions 
oshv1Tran_BaseMean <- as.matrix(resoshv1Tran_df[, "baseMean", drop=F])
oshv1Tran_backG <- genefinder(oshv1_6_BaseMean, anSig$ensembl_gene_id, 10, method= "manhattan")

####BLAST2GO after subsetting genes####

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
resBacTran_05 <- results(ddsBacTran, contrast = c("condition", "control", "treatment"), alpha=0.05)
summary(resBacTran_05)
sum(resBacTran_05$padj < 0.05, na.rm=TRUE) #1889

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

#Check how many genes have BH adjusted p values of less than 0.01 after P-value correction?
sum( resBacTran_05_df$padj < 0.05, na.rm=TRUE ) #1398 (before p-value correction it was )

#Subset the results table to the differentially expressed genes under FDR 0.01, order the Log2FC table first by strongest down regulation
resBacTran_05_dfSig <- resBacTran_05_df[ which(resBacTran_05_df$padj < 0.05 ), ]
head( resBacTran_05_dfSig[ order( resBacTran_05_dfSig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resBacTran_05_dfSig[ order( resBacTran_05_dfSig$log2FoldChange ), ] ) #tail for strongest up regulation

####Visualize Results with Diagnostic Plots####
#MA plot, useful overview for experiment with two-group comparison. Plots log2FC over mean of normalized counts
#genes with adjusted p value 
plotMA(resBacTran_05_df)
plotMA(resBacTran_05_dfSig)

#Export Results to CSV
write.csv( as.data.frame(resBacTran_05_df), file="resBacTran_05_df.csv")
write.csv( as.data.frame(resBacTran_05_dfSig), file="resBacTran_05_dfSig.csv")

####Use Bash script fetchEnsembl_ID.sh to get the Ensembl IDs ####
#Extract gene titles of all the significantly differentially expressed genes

####Bacterial Challenge Gene Set Enrichment Analysis ####
#Matching the background set
#Get average expressions 
#Bac_BaseMean <- as.matrix(resoshv1_6_df[, "baseMean", drop=F])
#Bac_backG <- genefinder(oshv1_6_BaseMean, anSig$ensembl_gene_id, 10, method= "manhattan")

#### Look at Differential Expression of Transcript Isoforms using Transcript Count Matrix ####

####Use Bash script fetchEnsembl_ID.sh to get the Ensembl IDs ####

#references: 
#https://github.com/sr320/LabDocs/tree/master/code/DESeq
#Code adapted from Stephen Roberts DESeq2 Github script called "SCRIPT_DESeq_CLAM_no replication.R" 
# StringTie manual : http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual

#SCRIPT_DESeq_SeaFan_replication.R also consulted, https://github.com/sr320/LabDocs/blob/master/code/DESeq/SCRIPT_DESeq_SeaFan_replication.R
