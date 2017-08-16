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
# Construct Full_PHENO_DATA.csv that contains SRA run information, such as which contrast, tissue, etc

#load gene(/transcript) count matrix and labels
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

# Check all sample IDs in oshv1ColData are also in oshv1CountData and match their orders
all(rownames(oshv1ColData) %in% colnames(oshv1CountData))  #Should return TRUE
  # returns TRUE

all(rownames(oshv1ColData) == colnames(oshv1CountData))    # should return TRUE
  #returns TRUE

#### Create OsHV-1 DESeq Data Set from Matrix ####
# DESeqDataSet from count matrix and labels 

ddsOshv1 <- DESeqDataSetFromMatrix(countData = oshv1CountData, 
  colData = oshv1ColData, 
  design = ~ condition)
# review how the data set looks
head(ddsOshv1)
####OsHV1 Differential Gene Expression Analysis ####

#subset columns by time point comparison (not comparing across time)
ddsOshv1_6 <- ddsOshv1[ , ddsOshv1$time.h. == 6 ]
ddsOshv1_12 <- ddsOshv1[ , ddsOshv1$time.h. == 12 ]
ddsOshv1_24 <- ddsOshv1[ , ddsOshv1$time.h. == 24 ]
ddsOshv1_48 <- ddsOshv1[ , ddsOshv1$time.h. == 48 ]
ddsOshv1_120 <- ddsOshv1[ , ddsOshv1$time.h. == 120 ]

#Relevel each to make sure that control is the first level in the treatment factor for each
ddsOshv1_6$condition <- relevel( ddsOshv1_6$condition, "control")
ddsOshv1_12$condition <- relevel( ddsOshv1_12$condition, "control")
ddsOshv1_24$condition <- relevel( ddsOshv1_24$condition, "control")
ddsOshv1_48$condition <- relevel( ddsOshv1_48$condition, "control")
ddsOshv1_120$condition <- relevel( ddsOshv1_120$condition, "control")

#Check we're looking at the right samples
as.data.frame( colData(ddsOshv1_6) )
as.data.frame( colData(ddsOshv1_12) )
as.data.frame( colData(ddsOshv1_24) )
as.data.frame( colData(ddsOshv1_48) )
as.data.frame( colData(ddsOshv1_120) )

#Running the DEG pipeline
ddsOshv1_6_DE <- DESeq(ddsOshv1_6)
ddsOshv1_12_DE <- DESeq(ddsOshv1_12)
ddsOshv1_24_DE <- DESeq(ddsOshv1_24)
ddsOshv1_48_DE <- DESeq(ddsOshv1_48)
ddsOshv1_120_DE <- DESeq(ddsOshv1_120)

#Inspect results table
#to extract just log2fold change and p values
resoshv1_6 <- results(ddsOshv1_6_DE)
resoshv1_6_data_frame <- as.data.frame(resoshv1_6)
head(resoshv1_6)
resoshv1_12 <- results(ddsOshv1_12_DE)
resoshv1_12_data_frame <- as.data.frame(resoshv1_12)
head(resoshv1_12)
resoshv1_24 <- results(ddsOshv1_24_DE)
resoshv1_24_data_frame <- as.data.frame(resoshv1_24)
head(resoshv1_24)
resoshv1_48 <- results(ddsOshv1_48_DE)
resoshv1_48_data_frame <- as.data.frame(resoshv1_48)
head(resoshv1_48)
resoshv1_120 <- results(ddsOshv1_120_DE)
resoshv1_120_data_frame <- as.data.frame(resoshv1_120)
head(resoshv1_120)

####p-value correction#### adapted from "Differential expression analysis of RNA-Seq data using DESeq2" Klaus 2014
#First Visualize histograms
#histogram of P- values to visualize any "hills" or "U shape"
# hill means variance of the null distribution too high, U shape means variance assumed too low
hist(resoshv1_6_data_frame$pvalue, breaks = 20, col = "grey") #hill, 
hist(resoshv1_12_data_frame$pvalue, breaks = 20, col = "grey") #hill
hist(resoshv1_24_data_frame$pvalue, breaks = 20, col = "grey") #U shape
hist(resoshv1_48_data_frame$pvalue, breaks = 20, col = "grey") #U shape
hist(resoshv1_120_data_frame$pvalue, breaks = 20, col = "grey") #U shape

#remove filtered out genes by independent filtering, they have NA adj. pvals
resoshv1_6_df <- resoshv1_6_data_frame[ !is.na(resoshv1_6_data_frame$padj), ]
resoshv1_12_df <- resoshv1_12_data_frame[ !is.na(resoshv1_12_data_frame$padj), ]
resoshv1_24_df <- resoshv1_24_data_frame[ !is.na(resoshv1_24_data_frame$padj), ]
resoshv1_48_df <- resoshv1_48_data_frame[ !is.na(resoshv1_48_data_frame$padj), ]
resoshv1_120_df <- resoshv1_120_data_frame[ !is.na(resoshv1_120_data_frame$padj), ]

#remove genes with NA pvals (outliers)
resoshv1_6_df <- resoshv1_6_df[ !is.na(resoshv1_6_df$pvalue), ]
resoshv1_12_df <- resoshv1_12_df[ !is.na(resoshv1_12_df$pvalue), ]
resoshv1_24_df <- resoshv1_24_df[ !is.na(resoshv1_24_df$pvalue), ]
resoshv1_48_df <- resoshv1_48_df[ !is.na(resoshv1_48_df$pvalue), ]
resoshv1_120_df <- resoshv1_120_df[ !is.na(resoshv1_120_df$pvalue), ]

#remove adjsuted pvalues, since we add the fdrtool results later on (based on the correct p-values)
resoshv1_6_df <- resoshv1_6_df[, -which(names(resoshv1_6_df) == "padj")]
resoshv1_12_df <- resoshv1_12_df[, -which(names(resoshv1_12_df) == "padj")]
resoshv1_24_df <- resoshv1_24_df[, -which(names(resoshv1_24_df) == "padj")]
resoshv1_48_df <- resoshv1_48_df[, -which(names(resoshv1_48_df) == "padj")]
resoshv1_120_df <- resoshv1_120_df[, -which(names(resoshv1_120_df) == "padj")]

#use z-scores as input to FDRtool to re-estimate the p-value
FDR.resoshv1_6_df <- fdrtool(resoshv1_6_df$stat, statistic= "normal", plot = T)
FDR.resoshv1_12_df <- fdrtool(resoshv1_12_df$stat, statistic= "normal", plot = T)
FDR.resoshv1_24_df <- fdrtool(resoshv1_24_df$stat, statistic= "normal", plot = T)
FDR.resoshv1_48_df <- fdrtool(resoshv1_48_df$stat, statistic= "normal", plot = T)
FDR.resoshv1_120_df <- fdrtool(resoshv1_120_df$stat, statistic= "normal", plot = T)

#add values to the results data frame, also ad new BH- adjusted p-values
resoshv1_6_df[,"padj"] <- p.adjust(FDR.resoshv1_6_df$pval, method = "BH")
resoshv1_12_df[,"padj"] <- p.adjust(FDR.resoshv1_12_df$pval, method = "BH")
resoshv1_24_df[,"padj"] <- p.adjust(FDR.resoshv1_24_df$pval, method = "BH")
resoshv1_48_df[,"padj"] <- p.adjust(FDR.resoshv1_48_df$pval, method = "BH")
resoshv1_120_df[,"padj"] <- p.adjust(FDR.resoshv1_120_df$pval, method = "BH")

#replot corrected p-values
hist(FDR.resoshv1_6_df$pval, col = "royalblue4",
     main = "Correct null model 6hr OsHV-1", xlab = "CORRECTED p-values")
hist(FDR.resoshv1_12_df$pval, col = "royalblue4",
     main = "Correct null model 12hr OsHV-1", xlab = "CORRECTED p-values")
hist(FDR.resoshv1_24_df$pval, col = "royalblue4",
     main = "Correct null model 24hr OsHV-1", xlab = "CORRECTED p-values")
hist(FDR.resoshv1_48_df$pval, col = "royalblue4",
     main = "Correct null model 48hr OsHV-1", xlab = "CORRECTED p-values")
hist(FDR.resoshv1_120_df$pval, col = "royalblue4",
     main = "Correct null model 120hr OSHV-1", xlab = "CORRECTED p-values")

#Check how many genes have BH adjusted p values of less than 0.01 after P-value correction?
sum( resoshv1_6_df$padj < 0.1, na.rm=TRUE ) #0 (before p-value correction it was 694)
sum( resoshv1_12_df$padj < 0.1, na.rm=TRUE ) # 833 (before p-value correction it was 747)
sum( resoshv1_24_df$padj < 0.1, na.rm=TRUE ) # 881 (before p-value correction it was 4084)
sum( resoshv1_48_df$padj < 0.1, na.rm=TRUE ) # 366 (before p-value correction it was 1026)
sum( resoshv1_120_df$padj < 0.1, na.rm=TRUE ) # 421 (before p-value correction it was 1137)

#Subset the results table to the differentially expressed genes under FDR 0.01, order the Log2FC table first by strongest down regulation
#6hr challenge
resoshv1_6_dfSig <- resoshv1_6_df[ which(resoshv1_6_df$padj < 0.1 ), ]
head( resoshv1_6_dfSig[ order( resoshv1_6_dfSig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resoshv1_6_dfSig[ order( resoshv1_6_dfSig$log2FoldChange ), ] ) #tail for strongest up regulation

#12 hr challenge
resoshv1_12_dfSig <- resoshv1_12_df[ which(resoshv1_12_df$padj < 0.1 ), ]
head( resoshv1_12_dfSig[ order( resoshv1_12_dfSig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resoshv1_12_dfSig[ order( resoshv1_12_dfSig$log2FoldChange ), ] ) #tail for strongest up regulation

#24hr challenge
resoshv1_24_dfSig <- resoshv1_24_df[ which(resoshv1_24_df$padj < 0.1 ), ]
head( resoshv1_24_dfSig[ order( resoshv1_24_dfSig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resoshv1_24_dfSig[ order( resoshv1_24_dfSig$log2FoldChange ), ] ) #tail for strongest up regulation

#48 hr challenge
resoshv1_48_dfSig <- resoshv1_48_df[ which(resoshv1_48_df$padj < 0.1 ), ]
head( resoshv1_48_dfSig[ order( resoshv1_48_dfSig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resoshv1_48_dfSig[ order( resoshv1_48_dfSig$log2FoldChange ), ] ) #tail for strongest up regulation

#120 hr challenge
resoshv1_120_dfSig <- resoshv1_120_df[ which(resoshv1_120_df$padj < 0.1 ), ]
head( resoshv1_120_dfSig[ order( resoshv1_120_dfSig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resoshv1_120_dfSig[ order( resoshv1_120_dfSig$log2FoldChange ), ] ) #tail for strongest up regulation

####Visualize Results with Diagnostic Plots####
#MA plot, useful overview for experiment with two-group comparison. Plots log2FC over mean of normalized counts
#genes with adjusted p value 
plotMA(resoshv1_6_df)
plotMA(resoshv1_12_df)
plotMA(resoshv1_24_df)
plotMA(resoshv1_48_df)
plotMA(resoshv1_120_df)

#Export Results to CSV
write.csv( as.data.frame(resoshv1_6_df), file="resoshv1_6_df.csv")
write.csv( as.data.frame(resoshv1_6_dfSig), file="resoshv1_6_dfSig.csv")
write.csv( as.data.frame(resoshv1_12_df), file="resoshv1_12_df.csv")
write.csv( as.data.frame(resoshv1_12_dfSig), file="resoshv1_12_dfSig.csv")
write.csv( as.data.frame(resoshv1_24_df), file="resoshv1_24_df.csv")
write.csv( as.data.frame(resoshv1_24_dfSig), file="resoshv1_24_dfSig.csv")
write.csv( as.data.frame(resoshv1_48_df), file="resoshv1_48_df.csv")
write.csv( as.data.frame(resoshv1_48_dfSig), file="resoshv1_48_dfSig.csv")
write.csv( as.data.frame(resoshv1_120_df), file="resoshv1_120_df.csv")
write.csv( as.data.frame(resoshv1_120_dfSig), file="resoshv1_120_dfSig.csv")


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


####Perform GO Analysis in R on Samples####






#references: 
#https://github.com/sr320/LabDocs/tree/master/code/DESeq
#Code adapted from Stephen Roberts DESeq2 Github script called "SCRIPT_DESeq_CLAM_no replication.R" 
# StringTie manual : http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual

SCRIPT_DESeq_SeaFan_replication.R also consulted, https://github.com/sr320/LabDocs/blob/master/code/DESeq/SCRIPT_DESeq_SeaFan_replication.R
