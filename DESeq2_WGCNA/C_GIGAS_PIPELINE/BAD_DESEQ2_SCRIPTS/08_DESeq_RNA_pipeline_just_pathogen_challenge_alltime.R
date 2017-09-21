#08_DESeq2_RNA_pipeline.R

#This script takes as input the ouput gene_count_matrix.csv data prepared from prepDE.py

#call the DESeq2 library 
source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade") 
biocLite()
biocLite("DESeq2")
library("DESeq2")
install.packages("fdrtool")
library(fdrtool)
source("https://bioconductor.org/biocLite.R")
biocLite("genefilter")
library(genefilter)
# Construct Full_PHENO_DATA.csv that contains SRA run information, such as which contrast, tissue, etc

####DEG Analysis with Gene Count Matrix ####
#load gene count matrix and labels
#Full_PHENO_DATA file contains metadata on the count table's samples
###Make sure excel Full_PHENO_DATA is in the same order or these commands will change data to be wrong!!!!###

FullCountData <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
head(FullCountData)

####Subset Data for OsHV1 ####
#Extract columns you want from the FullCountData, based on which column the correct SRA data is in for the
#Extract columns 1:30 from the FullCountData, these are the SRA's from the OsHV-1 experiment
oshv1CountData <- as.matrix(FullCountData[ , c(1:30)])
oshv1ColData <- read.csv("OsHV1_PHENO_DATA.csv", header=TRUE, sep=",")
rownames(oshv1ColData) <- oshv1ColData$sampleID
colnames(oshv1CountData) <- oshv1ColData$sampleID
head(oshv1CountData)
head(oshv1ColData)
#Give the time.h. column levels
oshv1ColData$time.h. <- factor(oshv1ColData$time.h.)
levels(oshv1ColData$time.h.) #check to see that it has levels 

# Check all sample IDs in oshv1ColData are also in oshv1CountData and match their orders
all(rownames(oshv1ColData) %in% colnames(oshv1CountData))  #Should return TRUE
  # returns TRUE
all(rownames(oshv1ColData) == colnames(oshv1CountData))    # should return TRUE
  #returns TRUE

#### Create OsHV-1 DESeq Data Set from Matrix ####
# DESeqDataSet from count matrix and labels 
#design purposefully doesn't account for time, not looking at time interaction, just condition
ddsOshv1 <- DESeqDataSetFromMatrix(countData = oshv1CountData, 
  colData = oshv1ColData, 
  design = ~condition)

#if design was design = ~time.h. + condition, #this design will gather the effect of condition, accounting for the sample pairing by time
# review how the data set looks
head(ddsOshv1)

#Relevel each to make sure that control is the first level in the treatment factor for each
ddsOshv1$condition <- relevel( ddsOshv1$condition, "control")

#Check we're looking at the right samples
as.data.frame( colData(ddsOshv1) )

#Running the DEG pipeline
ddsOshv1<- DESeq(ddsOshv1) #for designs with interactions, recommends setting betaPrior=FALSE

#Inspect results table
#extract contrasts between control and treatment values
resoshv1 <- results(ddsOshv1, contrast = c("condition", "control", "treatment"))
#to extract log2fold change and p values under 0.1 and 0.05
head(resoshv1)
summary(resoshv1)
sum(resoshv1$padj < 0.1, na.rm=TRUE) #4601
resoshv1_05 <- results(ddsOshv1,alpha=0.05)
summary(resoshv1_05)
sum(resoshv1_05$padj < 0.05, na.rm=TRUE) #3502

#metadata on meaning of the columns
mcols(resoshvl_05, use.names = TRUE)
  #shows treatment vs. control with baseMean as the mean of normalized counts for all samples

####p-value correction#### adapted from "Differential expression analysis of RNA-Seq data using DESeq2" Klaus 2014
#First Visualize histograms
#histogram of P- values to visualize any "hills" or "U shape"
# hill means variance of the null distribution too high, U shape means variance assumed too low
hist(resoshvl_05$pvalue, breaks = 20, col = "grey") #hill

#remove filtered out genes by independent filtering, they have NA adj. pvals
resoshv1_05_df <- resoshvl_05[ !is.na(resoshvl_05$padj), ]

#remove genes with NA pvals (outliers)
resoshv1_05_df <- resoshv1_05_df[ !is.na(resoshv1_05_df$pvalue), ]

#remove adjsuted pvalues, since we add the fdrtool results later on (based on the correct p-values)
resoshv1_05_df <- resoshv1_05_df[, -which(names(resoshv1_05_df) == "padj")]

#use z-scores as input to FDRtool to re-estimate the p-value
FDR.resoshv1_05_df <- fdrtool(resoshv1_05_df$stat, statistic= "normal", plot = T)

#add values to the results data frame, also ad new BH- adjusted p-values
resoshv1_05_df[,"padj"] <- p.adjust(FDR.resoshv1_05_df$pval, method = "BH")

#replot corrected p-values
hist(FDR.resoshv1_05_df$pval, col = "royalblue4",
     main = "Correct null model full OsHV-1", xlab = "CORRECTED p-values")

#Check how many genes have BH adjusted p values of less than 0.01 after P-value correction?
sum( resoshv1_05_df$padj < 0.05, na.rm=TRUE ) #663 (before p-value correction it was )

#Subset the results table to the differentially expressed genes under FDR 0.01, order the Log2FC table first by strongest down regulation
resoshv1_05_dfSig <- resoshv1_05_df[ which(resoshv1_05_df$padj < 0.05 ), ]
head( resoshv1_05_dfSig[ order( resoshv1_05_dfSig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resoshv1_05_dfSig[ order( resoshv1_05_dfSig$log2FoldChange ), ] ) #tail for strongest up regulation

####Visualize Results with Diagnostic Plots####
#MA plot, useful overview for experiment with two-group comparison. Plots log2FC over mean of normalized counts
#genes with adjusted p value 
plotMA(resoshv1_05_df)
plotMA(resoshv1_05_dfSig)

#Export Results to CSV
write.csv( as.data.frame(resoshv1_05_df), file="resoshv1_05_df.csv")
write.csv( as.data.frame(resoshv1_05_dfSig), file="resoshv1_05_dfSig.csv")

####Use Bash script fetchEnsembl_ID.sh to get the Ensembl IDs ####
#Extract gene titles of all the significantly differentially expressed genes

####OsHV1 Gene Set Enrichment Analysis ####
#Matching the background set
#Get average expressions 
oshv1_6_BaseMean <- as.matrix(resoshv1_6_df[, "baseMean", drop=F])
oshv1_6_backG <- genefinder(oshv1_6_BaseMean, anSig$ensembl_gene_id, 10, method= "manhattan")

####BLAST2GO after subsetting genes####

####Bacterial Challenge Differential Gene Expression Analysis ####
#Gram negative SRA's (including control): (SRR796597, SRR796596, SRR796595, SRR796594, SRR796593, SRR796592, SRR796589)
# Gram positive (including control): (SRR796598,  SRR796589)

#Subset Bacterial Challenge Data
BacCountData <- as.matrix(FullCountData[ , c(124:131)])
head(BacCountData)
BacColData <- read.csv("bac_PHENO_DATA.csv", header=TRUE, sep=",")
rownames(BacColData) <- BacColData$sampleID
colnames(BacCountData) <- BacColData$sampleID
head(BacCountData)
print(BacColData)

# Check all sample IDs in BacColData are also in BacCountData and match their orders
all(rownames(BacColData) %in% colnames(BacCountData))  #Should return TRUE
# returns TRUE
all(rownames(BacColData) == colnames(BacCountData))    # should return TRUE
#returns TRUE

#Give the "level" column levels
BacColData$level <- factor(BacColData$level)
levels(BacColData$level) #check to see that it has levels 

#### Create Bacterial Challenge DESeq Data Set from Matrix ####
# DESeqDataSet from count matrix and labels 
ddsBac <- DESeqDataSetFromMatrix(countData = BacCountData, 
                                   colData = BacColData, 
                                   design = ~condition)

#this design will gather the effect of condition between control and each treatment
# review how the data set looks
head(ddsBac)

#Relevel each to make sure that control is the first level in the condition factor
ddsBac$condition <- relevel( ddsBac$condition, "control")

#Check we're looking at the right samples
as.data.frame( colData(ddsBac) )

#Running the DEG pipeline
ddsBac<- DESeq(ddsBac) #for designs with interactions, recommends setting betaPrior=FALSE


#Inspect results table
resultsNames(ddsBac)

#extract contrasts between control and treatment values
resBac <- results(ddsBac, contrast = c("condition", "control", "treatment"))
#to extract log2fold change and p values under 0.1 and 0.05
head(resBac)
summary(resBac)
sum(resBac$padj < 0.1, na.rm=TRUE) #303
resBac_05 <- results(ddsBac,alpha=0.05)
summary(resBac_05)
sum(resBac_05$padj < 0.05, na.rm=TRUE) #196
resultsNames(res)

#metadata on meaning of the columns
mcols(resBac_05, use.names = TRUE)
#shows treatment vs. control with baseMean as the mean of normalized counts for all samples

####p-value correction#### adapted from "Differential expression analysis of RNA-Seq data using DESeq2" Klaus 2014
#First Visualize histograms
#histogram of P- values to visualize any "hills" or "U shape"
# hill means variance of the null distribution too high, U shape means variance assumed too low
hist(resoshvl_05$pvalue, breaks = 20, col = "grey") #hill

#remove filtered out genes by independent filtering, they have NA adj. pvals
resoshv1_05_df <- resoshvl_05[ !is.na(resoshvl_05$padj), ]

#remove genes with NA pvals (outliers)
resoshv1_05_df <- resoshv1_05_df[ !is.na(resoshv1_05_df$pvalue), ]

#remove adjsuted pvalues, since we add the fdrtool results later on (based on the correct p-values)
resoshv1_05_df <- resoshv1_05_df[, -which(names(resoshv1_05_df) == "padj")]

#use z-scores as input to FDRtool to re-estimate the p-value
FDR.resoshv1_05_df <- fdrtool(resoshv1_05_df$stat, statistic= "normal", plot = T)

#add values to the results data frame, also ad new BH- adjusted p-values
resoshv1_05_df[,"padj"] <- p.adjust(FDR.resoshv1_05_df$pval, method = "BH")

#replot corrected p-values
hist(FDR.resoshv1_05_df$pval, col = "royalblue4",
     main = "Correct null model full OsHV-1", xlab = "CORRECTED p-values")

#Check how many genes have BH adjusted p values of less than 0.01 after P-value correction?
sum( resoshv1_05_df$padj < 0.05, na.rm=TRUE ) #663 (before p-value correction it was )

#Subset the results table to the differentially expressed genes under FDR 0.01, order the Log2FC table first by strongest down regulation
resoshv1_05_dfSig <- resoshv1_05_df[ which(resoshv1_05_df$padj < 0.05 ), ]
head( resoshv1_05_dfSig[ order( resoshv1_05_dfSig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resoshv1_05_dfSig[ order( resoshv1_05_dfSig$log2FoldChange ), ] ) #tail for strongest up regulation

####Visualize Results with Diagnostic Plots####
#MA plot, useful overview for experiment with two-group comparison. Plots log2FC over mean of normalized counts
#genes with adjusted p value 
plotMA(resoshv1_05_df)
plotMA(resoshv1_05_dfSig)

#Export Results to CSV
write.csv( as.data.frame(resoshv1_05_df), file="resoshv1_05_df.csv")
write.csv( as.data.frame(resoshv1_05_dfSig), file="resoshv1_05_dfSig.csv")

####Use Bash script fetchEnsembl_ID.sh to get the Ensembl IDs ####
#Extract gene titles of all the significantly differentially expressed genes

####OsHV1 Gene Set Enrichment Analysis ####
#Matching the background set
#Get average expressions 
oshv1_6_BaseMean <- as.matrix(resoshv1_6_df[, "baseMean", drop=F])
oshv1_6_backG <- genefinder(oshv1_6_BaseMean, anSig$ensembl_gene_id, 10, method= "manhattan")

####BLAST2GO after subsetting genes####
# DESeqDataSet from count matrix and labels
ddsBac <- DESeqDataSetFromMatrix(countData = BacCountData, colData = BacColData , design =~condition)

#Check we're looking at the right samples
as.data.frame( colData(ddsBac) )

#Run DEG Analysis
ddsBac <- DESeq(ddsBac)

#Inspect results table
#to extract just log2fold change and p values

#Check how many genes have BH adjusted p values of less than 0.01?
sum( resBac$padj < 0.1, na.rm=TRUE ) # 0

#Export Results to CSV
write.csv( as.data.frame(resBac), file="resBac.csv")

####Use Bash script fetchEnsembl_ID.sh to get the Ensembl IDs ####

####Perform GO Analysis in R on Samples####



#### Look at Differential Expression of Transcript Isoforms using Transcript Count Matrix ####

#load gene count matrix and labels
#Full_PHENO_DATA file contains metadata on the count table's samples
###Make sure excel Full_PHENO_DATA is in the same order or these commands will change data to be wrong!!!!###

FullCountData <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
head(FullCountData)

####Subset Data for OsHV1 ####
#Extract columns you want from the FullCountData, based on which column the correct SRA data is in for the
#Extract columns 1:30 from the FullCountData, these are the SRA's from the OsHV-1 experiment
oshv1CountData <- as.matrix(FullCountData[ , c(1:30)])
oshv1ColData <- read.csv("OsHV1_PHENO_DATA.csv", header=TRUE, sep=",")
rownames(oshv1ColData) <- oshv1ColData$sampleID
colnames(oshv1CountData) <- oshv1ColData$sampleID
head(oshv1CountData)
head(oshv1ColData)
#Give the time.h. column levels
oshv1ColData$time.h. <- factor(oshv1ColData$time.h.)
levels(oshv1ColData$time.h.) #check to see that it has levels 

# Check all sample IDs in oshv1ColData are also in oshv1CountData and match their orders
all(rownames(oshv1ColData) %in% colnames(oshv1CountData))  #Should return TRUE
# returns TRUE
all(rownames(oshv1ColData) == colnames(oshv1CountData))    # should return TRUE
#returns TRUE

#### Create OsHV-1 DESeq Data Set from Matrix ####
# DESeqDataSet from count matrix and labels 
#design purposefully doesn't account for time, not looking at time interaction, just condition
ddsOshv1 <- DESeqDataSetFromMatrix(countData = oshv1CountData, 
                                   colData = oshv1ColData, 
                                   design = ~condition)

#if design was design = ~time.h. + condition, #this design will gather the effect of condition, accounting for the sample pairing by time
# review how the data set looks
head(ddsOshv1)

####OsHV1 Differential Gene Expression Analysis ####
#Relevel each to make sure that control is the first level in the treatment factor for each
ddsOshv1$condition <- relevel( ddsOshv1$condition, "control")

#Check we're looking at the right samples
as.data.frame( colData(ddsOshv1) )

#Running the DEG pipeline
ddsOshv1<- DESeq(ddsOshv1) #for designs with interactions, recommends setting betaPrior=FALSE

#Inspect results table
#extract contrasts between control and treatment values
resoshv1 <- results(ddsOshv1, contrast = c("condition", "control", "treatment"))
#to extract log2fold change and p values under 0.1 and 0.05
head(resoshv1)
summary(resoshv1)
sum(resoshv1$padj < 0.1, na.rm=TRUE) #4601
resoshv1_05 <- results(ddsOshv1,alpha=0.05)
summary(resoshv1_05)
sum(resoshv1_05$padj < 0.05, na.rm=TRUE) #3502

#metadata on meaning of the columns
mcols(resoshvl_05, use.names = TRUE)
#shows treatment vs. control with baseMean as the mean of normalized counts for all samples

####p-value correction#### adapted from "Differential expression analysis of RNA-Seq data using DESeq2" Klaus 2014
#First Visualize histograms
#histogram of P- values to visualize any "hills" or "U shape"
# hill means variance of the null distribution too high, U shape means variance assumed too low
hist(resoshvl_05$pvalue, breaks = 20, col = "grey") #hill

#remove filtered out genes by independent filtering, they have NA adj. pvals
resoshv1_05_df <- resoshvl_05[ !is.na(resoshvl_05$padj), ]

#remove genes with NA pvals (outliers)
resoshv1_05_df <- resoshv1_05_df[ !is.na(resoshv1_05_df$pvalue), ]

#remove adjsuted pvalues, since we add the fdrtool results later on (based on the correct p-values)
resoshv1_05_df <- resoshv1_05_df[, -which(names(resoshv1_05_df) == "padj")]

#use z-scores as input to FDRtool to re-estimate the p-value
FDR.resoshv1_05_df <- fdrtool(resoshv1_05_df$stat, statistic= "normal", plot = T)

#add values to the results data frame, also ad new BH- adjusted p-values
resoshv1_05_df[,"padj"] <- p.adjust(FDR.resoshv1_05_df$pval, method = "BH")

#replot corrected p-values
hist(FDR.resoshv1_05_df$pval, col = "royalblue4",
     main = "Correct null model full OsHV-1", xlab = "CORRECTED p-values")

#Check how many genes have BH adjusted p values of less than 0.01 after P-value correction?
sum( resoshv1_05_df$padj < 0.05, na.rm=TRUE ) #663 (before p-value correction it was )

#Subset the results table to the differentially expressed genes under FDR 0.01, order the Log2FC table first by strongest down regulation
resoshv1_05_dfSig <- resoshv1_05_df[ which(resoshv1_05_df$padj < 0.05 ), ]
head( resoshv1_05_dfSig[ order( resoshv1_05_dfSig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resoshv1_05_dfSig[ order( resoshv1_05_dfSig$log2FoldChange ), ] ) #tail for strongest up regulation

####Visualize Results with Diagnostic Plots####
#MA plot, useful overview for experiment with two-group comparison. Plots log2FC over mean of normalized counts
#genes with adjusted p value 
plotMA(resoshv1_05_df)
plotMA(resoshv1_05_dfSig)

#Export Results to CSV
write.csv( as.data.frame(resoshv1_05_df), file="resoshv1_05_df.csv")
write.csv( as.data.frame(resoshv1_05_dfSig), file="resoshv1_05_dfSig.csv")

####Get .csv of all differentially expressed genes in contrasts####
#Extract gene titles of all the significantly differentially expressed genes

####BLAST2GO after subsetting genes####

####OsHV1 Gene Set Enrichment Analysis ####
#Matching the background set
#Get average expressions 
oshv1_6_BaseMean <- as.matrix(resoshv1_6_df[, "baseMean", drop=F])
oshv1_6_backG <- genefinder(oshv1_6_BaseMean, anSig$ensembl_gene_id, 10, method= "manhattan")

####Bacterial Challenge Differential Gene Expression Analysis ####

#Gram negative SRA's (including control): (SRR796597, SRR796596, SRR796595, SRR796594, SRR796593, SRR796592, SRR796589)
# Gram positive (including control): (SRR796598,  SRR796589)
#make a new DESeq2 matrix for each

BacCountData <- as.matrix(FullCountData[ , c(124:131)])
head(BacCountData)
BacColData <- read.csv("bac_PHENO_DATA.csv", header=TRUE, sep=",")
rownames(BacColData) <- BacColData$sampleID
colnames(BacCountData) <- BacColData$sampleID
head(BacCountData)
print(BacColData)

# Check all sample IDs in gramNegColData are also in gramNegCountData and match their orders
all(rownames(BacColData) %in% colnames(BacCountData))  #Should return TRUE
# returns TRUE
all(rownames(BacColData) == colnames(BacCountData))    # should return TRUE
#returns TRUE

# DESeqDataSet from count matrix and labels
ddsBac <- DESeqDataSetFromMatrix(countData = BacCountData, colData = BacColData , design =~condition)

#Check we're looking at the right samples
as.data.frame( colData(ddsBac) )

#Run DEG Analysis
ddsBac <- DESeq(ddsBac)

#Inspect results table
#to extract just log2fold change and p values

#Check how many genes have BH adjusted p values of less than 0.01?
sum( resBac$padj < 0.1, na.rm=TRUE ) # 0

#Export Results to CSV
write.csv( as.data.frame(resBac), file="resBac.csv")

####Use Bash script fetchEnsembl_ID.sh to get the Ensembl IDs ####

####Perform GO Analysis in R on Samples####




#references: 
#https://github.com/sr320/LabDocs/tree/master/code/DESeq
#Code adapted from Stephen Roberts DESeq2 Github script called "SCRIPT_DESeq_CLAM_no replication.R" 
# StringTie manual : http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual

SCRIPT_DESeq_SeaFan_replication.R also consulted, https://github.com/sr320/LabDocs/blob/master/code/DESeq/SCRIPT_DESeq_SeaFan_replication.R
