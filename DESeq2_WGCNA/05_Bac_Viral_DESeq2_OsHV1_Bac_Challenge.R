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
sum( resoshv1Tran_05_df$padj < 0.05, na.rm=TRUE ) #1462 (before p-value correction it was )... now its 2

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
write.csv( as.data.frame(resoshv1Tran_05_df), file="OsHV1_resoshv1Tran_05_df.csv")
write.csv( as.data.frame(resoshv1Tran_05_dfSig), file="OsHV1_resoshv1Tran_05_dfSig.csv")

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
transcriptIDstring <- toString(transcriptIDdf[,3], sep=',')
transcriptIDstring
write(transcriptIDstring, "transcriptIDstring", sep = ",")
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

####Uploaded transcripts from the UniProt.ws website ####
oshv1_transcriptIDs_UniProt <- read.csv("oshv1_transcriptIDstring_UniProtKB.tab", sep = "\t")
head(oshv1_transcriptIDs_UniProt)

####Extract GIMAP/IAN proteins and CgIAPs ####
#Strategy: Based on which ones BLAST to GIMAP proteins?


####OsHV1 Gene Set Enrichment Analysis topGO ####
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

#references: 
#https://github.com/sr320/LabDocs/tree/master/code/DESeq
#StringTie manual : http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual