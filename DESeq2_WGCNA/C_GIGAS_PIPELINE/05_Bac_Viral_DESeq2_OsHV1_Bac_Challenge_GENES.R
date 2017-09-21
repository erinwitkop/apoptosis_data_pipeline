#05_Bac_viral_DESeq2_OsHV1_Bac_challenge_GENE_SET

####INPUT DATA GATHERED FROM CRASSOSTREA GIGAS####

#This script takes as input the output Bac_Viral_gene_count_matrix.csv data prepared from prepDE.py and performs
#differential Gene expression analysis, performs gene set enrichment, isolates apoptosis-related processes,
#graphs their relative abundance, and creates heatmaps. 

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
install.packages("dplyr")
library(dplyr)
install.packages("tidyr")
library(tidyr)
#install.packages("reshape")
library(reshape)
#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
#library(biomaRt)
library(ggplot2)
#install.packages("pheatmap")
library(pheatmap)
# Construct Bac_Viral_PHENO_DATA.csv that contains SRA run information, such as which contrast, tissue, etc.

####DEG Analysis with GENE Count Matrix ####
#load GENE count matrix and labels
#Bac_Viral_PHENO_DATA.csv file contains metadata on the count table's samples
###Make sure excel PHENODATA is in the same order or these commands will change data to be wrong!!!!###

GeneCountData <- as.matrix(read.csv("Bac_Viral_gene_count_matrix.csv", row.names="gene_id"))
head(GeneCountData)

####Subset Data for OsHV1 ####
#Extract columns you want from the GeneCountData, based on which column the correct SRA data is in for the
#Extract columns 1:30 from the GeneCountData, these are the SRA's from the OsHV-1 experiment
oshv1GeneCountData <- as.matrix(GeneCountData[ , c(1:30)])
head(oshv1GeneCountData)
GeneColData <- read.csv("Bac_Viral_PHENO_DATA.csv", header=TRUE, sep=",")
oshv1GeneColData <- GeneColData[c(1:30),]
oshv1GeneColData <- oshv1GeneColData[, c("sampleID", "condition", "stressorLevel")]
print(oshv1GeneColData)
rownames(oshv1GeneColData) <- oshv1GeneColData$sampleID
colnames(oshv1GeneCountData) <- oshv1GeneColData$sampleID
head(oshv1GeneCountData)
head(oshv1GeneColData)
#Give the stressorLevel column levels
oshv1GeneColData$stressorLevel <- factor(oshv1GeneColData$stressorLevel)
levels(oshv1GeneColData$stressorLevel) #check to see that it has levels 

# Check all sample IDs in oshv1ColData are also in oshv1CountData and match their orders
all(rownames(oshv1GeneColData) %in% colnames(oshv1GeneCountData))  #Should return TRUE
# returns TRUE
all(rownames(oshv1GeneColData) == colnames(oshv1GeneCountData))    # should return TRUE
#returns TRUE

#### Create OsHV-1 DESeq Data Set from Matrix ####
# DESeqDataSet from count matrix and labels 
#design purposefully doesn't account for time, not looking at time interaction, just condition
ddsOshv1Gene <- DESeqDataSetFromMatrix(countData = oshv1GeneCountData, 
                                       colData = oshv1GeneColData, 
                                       design = ~ condition)
#
ddsOshv1Gene <- ddsOshv1Gene[ rowSums(counts(ddsOshv1Gene)) > 1, ]
#if design was design = ~stressorLevel + condition, #this design will gather the effect of condition, accounting for the sample pairing by time
# review how the data set looks
head(ddsOshv1Gene)

#Relevel each to make sure that control is the first level in the treatment factor for each
ddsOshv1Gene$condition <- relevel( ddsOshv1Gene$condition, "control")

#Check we're looking at the right samples
as.data.frame( colData(ddsOshv1Gene) )

#Running the DEG pipeline
ddsOshv1Gene<- DESeq(ddsOshv1Gene) #for designs with interactions, recommends setting betaPrior=FALSE

#Inspect results table
#extract contrasts between control and treatment values
resoshv1Gene<- results(ddsOshv1Gene, contrast = c("condition", "control", "treatment"))
#to extract log2fold change and p values under 0.1 and 0.05
head(resoshv1Gene)

#Order by Log2FC
head( resoshv1Gene[ order( resoshv1Gene$log2FoldChange ), ] ) #head for strongest downregulation
tail( resoshv1Gene[ order( resoshv1Gene$log2FoldChange ), ] ) #tail for strongest up regulation
summary(resoshv1Gene)
#sum(resoshv1Tran$padj < 0.1, na.rm=TRUE) #
#resoshv1Tran_05 <- results(ddsOshv1Tran,alpha=0.05)
#summary(resoshv1Tran_05)
sum(resoshv1Gene$padj < 0.05, na.rm=TRUE) #3776

#metadata on meaning of the columns
mcols(resoshv1Gene, use.names = TRUE)
#Get more detailed description
mcols(resoshv1Gene)$description
#shows treatment vs. control with baseMean as the mean of normalized counts for all samples

####p-value correction for all genes in resoshv1Gene ####
#First Visualize histograms
#histogram of P- values to visualize any "hills" or "U shape"
# hill means variance of the null distribution too high, U shape means variance assumed too low
hist(resoshv1Gene$pvalue, breaks = 20, col = "grey") #hill

#remove filtered out genes by independent filtering, they have NA adj. pvals
resoshv1Gene_df <- resoshv1Gene[ !is.na(resoshv1Gene$padj), ]

#remove genes with NA pvals (outliers)
resoshv1Gene_df <- resoshv1Gene_df[ !is.na(resoshv1Gene_df$pvalue), ]

#remove adjsuted pvalues, since we add the fdrtool results later on (based on the correct p-values)
resoshv1Gene_df <- resoshv1Gene_df[, -which(names(resoshv1Gene_df) == "padj")]

#use z-scores as input to FDRtool to re-estimate the p-value
FDR.resoshv1Gene_df <- fdrtool(resoshv1Gene_df$stat, statistic= "normal", plot = T)

#add values to the results data frame, also ad new BH- adjusted p-values
resoshv1Gene_df[,"padj"] <- p.adjust(FDR.resoshv1Gene_df$pval, method = "BH")

#replot corrected p-values 
hist(FDR.resoshv1Gene_df$pval, col = "royalblue4",
     main = "Correct null model OsHv1 Gene Count", xlab = "CORRECTED p-values")

#Check how many genes have BH adjusted p values of less than 0.05 after P-value correction?
sum( resoshv1Gene_df$padj < 0.05, na.rm=TRUE ) #435

#Subset the results table to the differentially expressed genes under FDR 0.1, order the Log2FC table first by strongest down regulation
resoshv1Gene_dfSig <- resoshv1Gene_df[ which(resoshv1Gene_df$padj < 0.05 ), ]
head( resoshv1Gene_dfSig[ order( resoshv1Gene_dfSig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resoshv1Gene_dfSig[ order( resoshv1Gene_dfSig$log2FoldChange ), ] ) #tail for strongest up regulation
summary(resoshv1Gene_dfSig)
resoshv1Gene_df_non_Sig <- resoshv1Gene_df[ which(resoshv1Gene_df$padj > 0.05 ), ]
summary(resoshv1Gene_df_non_Sig)

#Visualize Results with Diagnostic Plots#
#MA plot, useful overview for experiment with two-group comparison. Plots log2FC over mean of normalized counts
#genes with adjusted p value 
plotMA(resoshv1Gene_dfSig)
plotMA(resoshv1Gene_df_non_Sig)

#Export Results to CSV
write.csv( as.data.frame(resoshv1Gene_df), file="OsHV1_resoshv1Gene_df.csv")
write.csv( as.data.frame(resoshv1Gene_dfSig), file="OsHV1_resoshv1Gene_dfSig.csv")
write.csv( as.data.frame(resoshv1Gene_df_non_Sig), file = "OsHV1_resoshv1Gene_df_non_Sig.csv")


####Subset files for only those that have Gene IDs####
#Extract gene titles of all the significantly differentially expressed genes
OsHV1_withID_subset_resoshv1Gene_dfSig <- 
  resoshv1Gene_dfSig[grep("gene:", rownames(resoshv1Gene_dfSig)), ]
head(OsHV1_withID_subset_resoshv1Gene_dfSig) #0 genes are significant!

#GeneIDdf <- as.data.frame(rownames(OsHV1_withID_subset_resoshv1Gene_dfSig))
#head(GeneIDdf)
#GeneID1thru5 <- rownames(head(OsHV1_withID_subset_resoshv1Gene_dfSig[1:5,]))
#GeneIDdf= transform(GeneIDdf, 
                          ID = colsplit(rownames(OsHV1_withID_subset_resoshv1Gene_dfSig), 
                                        split = "\\:", names = c()))

#this line splitsup the word transcript and the transcript name
#GeneIDstring <- toString(GeneIDdf[,3], sep=',')
#GeneIDstring
#write(GeneIDstring, "GeneIDstring", sep = ",")
#write this to a file and then perform look up on the UniProt website

#Extract gene titles from non Sig genes
OsHV1_withID_subset_resoshv1Gene_df_non_Sig <- 
  resoshv1Gene_df_non_Sig[grep("gene:", rownames(resoshv1Gene_df_non_Sig)), ]
head(OsHV1_withID_subset_resoshv1Gene_df_non_Sig)
GeneIDdf_nonSig <- as.data.frame(rownames(OsHV1_withID_subset_resoshv1Gene_df_non_Sig))
head(GeneIDdf_nonSig)
GeneIDdf_nonSig= transform(GeneIDdf_nonSig, 
                                 ID = colsplit(rownames(OsHV1_withID_subset_resoshv1Gene_df_non_Sig), 
                                               split = "\\:", names = c('gene:', 'CGI_10024332')))
GeneIDstring_nonSig <- toString(GeneIDdf_nonSig[,3], sep=',')
GeneIDstring_nonSig
write(GeneIDstring_nonSig, "GeneIDstring_nonSig", sep = ",")
#write this to a file and then perform look up on the UniProt website

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
#oshv1_GeneIDs_UniProt_SIG <- read.csv("", header=TRUE)
#make sure to upload as csv version so that the protein names load correctly
#head(oshv1_GeneIDs_UniProt_SIG)

####Extract GIMAP/IAN proteins and CgIAPs from Significantly Differentially Expressed Genes####
#Significantly differentially Expressed IAPs, using several names
#oshv1_GeneIDs_UniProt_SIG_ProtNames <- oshv1_GeneIDs_UniProt_SIG$Protein.names
#oshv1_IAP_gene_SIG <- grepl("IAP", oshv1_GeneIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
#grep("TRUE", oshv1_IAP_gene_SIG) 

#oshv1_apoptosis_gene_SIG <- grepl("apoptosis", oshv1_GeneIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
#grep("TRUE", oshv1_apoptosis_geme_SIG) #3, 353
#oshv1_apoptosis_gene_SIG_info <- oshv1_GeneIDs_UniProt_SIG[294,]

#oshv1_inhibitor_SIG <-  grepl("inhibitor", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
#grep("TRUE", oshv1_inhibitor_SIG) # 24 149 294 349
#oshv1_inhibitor_SIG_info <- oshv1_transcriptIDs_UniProt_SIG[c(24, 149, 294, 349),] #only 294!

#oshv1_IAPs_SIG_info <- oshv1_transcriptIDs_UniProt_SIG[c(3,353, 294),]
#oshv1_IAPs_SIG_info


#Significant GIMAP Genes
#oshv1_IAN_Sig <- grepl("IAN", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
#grep("TRUE", oshv1_IAN_Sig) #0 TRUE
#oshv1_AIG_Sig <- grepl("AIG", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
#grep("TRUE", oshv1_AIG_Sig) #0 TRUE
#oshv1_IMAP_Sig <- grepl("IMAP", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
#grep("TRUE", oshv1_IMAP_Sig) #0 TRUE

#oshv1_GTP_SIG <- grepl("GTP", oshv1_transcriptIDs_UniProt_SIG$Protein.names, ignore.case = TRUE) 
#grep("TRUE", oshv1_GIMAP_SIG)
#oshv1_transcriptIDs_UniProt_SIG[c(53,54,65,85,125),] #53 and 65 are GIMAP! called GTPase IMAP
#oshv1_GIMAP_Sig_info <- oshv1_transcriptIDs_UniProt_SIG[c(53,65),]
#oshv1_GIMAP_Sig_info


#Uploaded NON Sig Differentially Expressed transcripts from the UniProt.ws website ####
oshv1_GeneIDs_UniProt_non_Sig <- read.csv("Bac_Viral_resoshv1Gene_df_non_sig_UniprotId.csv", header=TRUE)
#make sure to upload as csv version so that the protein names load correctly
head(oshv1_GeneIDs_UniProt_non_Sig)

####Extract GIMAP and IAP genes from NON significant genes, for comparison ####
##Expressed IAPs
#use grepl to find text strings 
oshv1_GeneIDs_UniProt_non_Sig_ProtNames <- oshv1_GeneIDs_UniProt_non_Sig$Protein.names
oshv1_IAP_gene_non_Sig <- grepl("IAP", oshv1_GeneIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_IAP_gene_non_Sig) #137
oshv1_apoptosis_gene_non_Sig <- grepl("apoptosis", oshv1_GeneIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_apoptosis_gene_non_Sig) #0
oshv1_inhibitor_gene_non_Sig <- grepl("inhibitor", oshv1_GeneIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_inhibitor_gene_non_Sig) #0 none

oshv1_IAP_gene_non_Sig_info <- oshv1_GeneIDs_UniProt_non_Sig[c(137),]
oshv1_IAP_gene_non_Sig_info 
# #only 1 non significant one

##Expressed GIMAP Genes, non significant
oshv1_GTP_gene_non_Sig <- grepl("GTP", oshv1_GeneIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE) 
grep("TRUE", oshv1_GTP_gene_non_Sig) # 138 161 212 (138 and 212 are GIMAP!)
oshv1_GIMAP_gene_non_sig <- grepl("IMAP",oshv1_GeneIDs_UniProt_non_Sig$Protein.names, ignore.case = TRUE)
grep("TRUE", oshv1_GIMAP_gene_non_sig)
#138, 212
oshv1_IAN_gene_non_Sig <- grepl("IAN", oshv1_GeneIDs_UniProt_non_Sig$Protein.names) 
grep("TRUE", oshv1_IAN_gene_non_Sig) #0

oshv1_GIMAP_gene_non_sig_info <- oshv1_GeneIDs_UniProt_non_Sig[c(138,212),] 
oshv1_GIMAP_gene_non_sig_info

####Link GIMAP and IAN genes, significant and non significant, with Expression Values####
#Sig IAPS
#none

#Sig GIMAPs
#none

#Non Sig IAP
oshv1_IAP_gene_non_Sig_info$Gene.names 
# CGI_10022060
subset_oshv1_IAP_gene_non_sig_info <- oshv1_IAP_gene_non_Sig_info[,c(4,5,7)]
#deleted the first row (69)
oshv1_IAP_gene_non_sig_info_resoshv1Gene_df <- grep("CGI_10022060", rownames(resoshv1Gene_df_non_Sig), ignore.case = TRUE) #11542
oshv1_IAP_gene_non_sig_info_resoshv1Gene_df <- resoshv1Gene_df_non_Sig[11542,]
oshv1_IAP_gene_non_sig_combined_FULL <- cbind(subset_oshv1_IAP_gene_non_sig_info,oshv1_IAP_gene_non_sig_info_resoshv1Gene_df)
oshv1_IAP_gene_non_sig_combined_FULL["Significance"] <- "Non_significant"
oshv1_IAP_gene_non_sig_combined_FULL["Challenge"] <- "OsHV1"
oshv1_IAP_gene_non_sig_combined_FULL["Type"] <- "IAP"

#Non Sig GIMAP
oshv1_GIMAP_gene_non_sig_info$Gene.names
#CGI_10027391 CGI_10026503

subset_GIMAP_gene_non_sig_info <- oshv1_GIMAP_gene_non_sig_info[,c(4,5,7)
                                                           ]
oshv1_GIMAP_gene_non_sig_info_resoshv1Tran_df <- grep("CGI_10027391", rownames(resoshv1Gene_df_non_Sig), ignore.case = TRUE) #11547
oshv1_GIMAP_gene_non_sig_info_resoshv1Tran_df <- resoshv1Gene_df_non_Sig[11547,]
oshv1_GIMAP_gene_non_sig_info_resoshv1Tran_df2 <- grep("CGI_10026503", rownames(resoshv1Gene_df_non_Sig), ignore.case = TRUE) #17545
oshv1_GIMAP_gene_non_sig_info_resoshv1Tran_df2 <- resoshv1Gene_df_non_Sig[17545,]

oshv1_GIMAP_gene_non_sig_combined_EXP <- rbind(oshv1_GIMAP_gene_non_sig_info_resoshv1Tran_df,oshv1_GIMAP_gene_non_sig_info_resoshv1Tran_df2)

oshv1_GIMAP_gene_non_sig_combined_FULL <- cbind(subset_GIMAP_gene_non_sig_info, oshv1_GIMAP_gene_non_sig_combined_EXP)
oshv1_GIMAP_gene_non_sig_combined_FULL["Significance"] <- "Non significant"
oshv1_GIMAP_gene_non_sig_combined_FULL["Challenge"] <- "OsHV1"
oshv1_GIMAP_gene_non_sig_combined_FULL["Type"] <- "GIMAP"

####Combine all Oshv1 IAP and GIMAP data #####
#to rbind the names of all the columns need to be the same, change first column name
oshv1_GIMAP_IAP_gene_combined_FULL <- rbind(oshv1_GIMAP_gene_non_sig_combined_FULL, oshv1_IAP_gene_non_sig_combined_FULL)

####OsHV1 Gene Set Enrichment Analysis topGO ####


####Bacterial Challenge Differential Gene Expression Analysis ####
#Gram negative SRA's (including control): (SRR796597, SRR796596, SRR796595, SRR796594, SRR796593, SRR796592, SRR796589)
# Gram positive (including control): (SRR796598,  SRR796589)

#Subset Bacterial Challenge Data
BacGeneCountData <- as.matrix(GeneCountData[ , c(31:38)])
head(BacGeneCountData)
BacGeneColData <- read.csv("bac_PHENO_DATA.csv", header=TRUE, sep=",")
rownames(BacGeneColData) <- BacGeneColData$sampleID
colnames(BacGeneCountData) <- BacGeneColData$sampleID
head(BacGeneCountData)
print(BacGeneColData)

# Check all sample IDs in BacColData are also in BacCountData and match their orders
all(rownames(BacGeneColData) %in% colnames(BacGeneCountData))  #Should return TRUE
# returns TRUE
all(rownames(BacGeneColData) == colnames(BacGeneCountData))    # should return TRUE
#returns TRUE

#Give the "level" column levels
BacGeneColData$level <- factor(BacGeneColData$level)
levels(BacGeneColData$level) #check to see that it has levels 

#### Create Bacterial Challenge DESeq Data Set from Matrix ####
# DESeqDataSet from count matrix and labels 
ddsBacGene <- DESeqDataSetFromMatrix(countData = BacGeneCountData, 
                                     colData = BacGeneColData, 
                                     design = ~condition)

#this design will gather the effect of condition between control and each treatment
# review how the data set looks
head(ddsBacGene)

#Relevel each to make sure that control is the first level in the condition factor
ddsBacGene$condition <- relevel( ddsBacGene$condition, "control")

#Check we're looking at the right samples
as.data.frame( colData(ddsBacGene) )

#Running the DEG pipeline
ddsBacGene <- DESeq(ddsBacGene) #for designs with interactions, recommends setting betaPrior=FALSE

#Inspect results table
resultsNames(ddsBacGene)

#extract contrasts between control and treatment values
resBacGene <- results(ddsBacGene, contrast = c("condition", "control", "treatment"))
#to extract log2fold change and p values under 0.1 and 0.05
head(resBacGene)
summary(resBacGene)
sum(resBacGene$padj < 0.1, na.rm=TRUE) #2322
#going to subset for 0.05 after the p value correction!
sum(resBacGene$padj < 0.05, na.rm=TRUE) #173

#metadata on meaning of the columns
mcols(resBacGene, use.names = TRUE)
#shows treatment vs. control with baseMean as the mean of normalized counts for all samples

####p-value correction#### adapted from "Differential expression analysis of RNA-Seq data using DESeq2" Klaus 2014
#First Visualize histograms
#histogram of P- values to visualize any "hills" or "U shape"
# hill means variance of the null distribution too high, U shape means variance assumed too low
hist(resBacGene$pvalue, breaks = 20, col = "grey") #hill

#remove filtered out genes by independent filtering, they have NA adj. pvals
resBacGene <- resBacGene[ !is.na(resBacGene$padj), ]

#remove genes with NA pvals (outliers)
resBacGene <- resBacGene[ !is.na(resBacGene$pvalue), ]

#remove adjsuted pvalues, since we add the fdrtool results later on (based on the correct p-values)
resBacGene <- resBacGene[, -which(names(resBacGene) == "padj")]

#use z-scores as input to FDRtool to re-estimate the p-value
FDR.resBacGene <- fdrtool(resBacGene$stat, statistic= "normal", plot = T)

#add values to the results data frame, also ad new BH- adjusted p-values
resBacGene[,"padj"] <- p.adjust(FDR.resBacGene$pval, method = "BH")

#replot corrected p-values
hist(FDR.resBacGene$pval, col = "royalblue4",
     main = "Correct null model Bacterial Challenge Gene Counts", xlab = "CORRECTED p-values")

#Check how many genes have BH adjusted p values of less than 0.05 after P-value correction?
sum( resBacGene$padj < 0.05, na.rm=TRUE ) #125 (before p-value correction it was )

#Subset the results table to the differentially expressed genes under FDR 0.01, order the Log2FC table first by strongest down regulation
resBacGene_Sig <- resBacGene[ which(resBacGene$padj < 0.05 ), ]
head( resBacGene_Sig[ order( resBacGene_Sig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resBacGene_Sig[ order( resBacGene_Sig$log2FoldChange ), ] ) #tail for strongest up regulation
resBacGene_non_Sig <- resBacGene[ which(resBacGene$padj > 0.05 ), ]
summary(resBacGene_non_Sig)

####Visualize Results with Diagnostic Plots####
#MA plot, useful overview for experiment with two-group comparison. Plots log2FC over mean of normalized counts
#genes with adjusted p value 
plotMA(resBacGene)
plotMA(resBacGene_Sig)
plotMA(resBacGene_non_Sig)

#Export Results to CSV
write.csv( as.data.frame(resBacGene), file="resBacGene.csv")
write.csv( as.data.frame(resBacGene_Sig), file="resBacGene_Sig.csv")
write.csv( as.data.frame(resBacGene_non_Sig), file="resBacGene_non_Sig.csv")

####Subset files for only those that have Gene IDs####
#Extract gene titles of all the significantly differentially expressed genes
Bac_withID_subset_resBacGene_Sig <- 
  resBacGene_Sig[grep("gene:", rownames(resBacGene_Sig)), ]
Bac_withID_subset_resBacGene_Sig #0 significant ones! 
BacGeneIDdf <- as.data.frame(rownames(Bac_withID_subset_resBacGene_Sig))
head(BacGeneIDdf)
#BacGeneIDdf= transform(BacGeneIDdf, 
                             ID = colsplit(rownames(Bac_withID_subset_resBacGene_Sig), 
                                           split = "\\:", names = c('')))
#this line splitsup the word transcript and the transcript name

#BacTranscriptIDstring <- toString(BacTranscriptIDdf[,3], sep=',')
#write(BacTranscriptIDstring, "BacTranscriptIDstring_sig", sep = ",")
#write this to a file and then perform look up on the UniProt website

#Extract gene titles from non Sig genes
Bac_withID_subset_resBacGene_non_Sig <- 
  resBacGene_non_Sig[grep("gene:", rownames(resBacGene_non_Sig)), ]
head(Bac_withID_subset_resBacGene_non_Sig)
BacGeneID_nonSig <- as.data.frame(rownames(Bac_withID_subset_resBacGene_non_Sig))
head(BacGeneID_nonSig)
BacGeneID_nonSig= transform(BacGeneID_nonSig, 
                                    ID = colsplit(rownames(Bac_withID_subset_resBacGene_non_Sig), 
                                                  split = "\\:", names = c('gene:', 'CGI_10024332')))
BacGeneID_nonSig <- toString(BacGeneID_nonSig[,3], sep=',')
BacGeneID_nonSig
write(BacGeneID_nonSig, "BacGeneID_nonSig", sep = ",")
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
Bac_IAP_SIG_combined_FULL["Significance"] <- "Significant" #add column about significance
#only 1 sig IAP

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
Bac_IAP_non_sig_combined_FULL["Significance"] <- "Non_significant" # add column about Significance

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
Bac_GIMAP_non_sig_combined_FULL["Significance"] <- "Non_significant" #add column about significance

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
COMBINED_GIMAP_IAP_cols <- COMBINED_GIMAP_IAP[,c(2,6,7,10,11)]

#Download data to simplify Protein Names for Viewing 
write.csv( as.data.frame(COMBINED_GIMAP_IAP_cols), file="COMBINED_GIMAP_IAP_cols.csv")
#Added rows for Gene Type, and edited rows for protein 
COMBINED_GIMAP_IAP_cols <- read.csv("COMBINED_GIMAP_IAP_cols_edited.csv", header = TRUE)
COMBINED_GIMAP_IAP_cols_Bac <- read.csv("COMBINED_GIMAP_IAP_cols_Bac.csv", header = TRUE)
COMBINED_GIMAP_IAP_cols_Oshv1 <- read.csv("COMBINED_GIMAP_IAP_cols_OsHV1.csv", header = TRUE)
colnames(COMBINED_GIMAP_IAP_cols_Bac)
#Plot the combined data set altogether for comparison per Challenge
#use geom_col and not geom_bar...geom bar makes the height of the bar proportional to the number
#of cases in each group

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Pot GIMAP and IAP genes for the Bacterial Challenge Data
label.Bac.df <- data.frame(Transcript= c("transcript:EKC41180"), log2FoldChange=c(-22.0))
COMBINED_GIMAP_IAP_BAC_PLOT <- ggplot(COMBINED_GIMAP_IAP_cols_Bac) + 
  geom_col(aes(x=Transcript, y=log2FoldChange, fill=as.factor(Type))) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  scale_y_continuous(name ="log2 Fold Change", breaks = scales::pretty_breaks(n = 20)) + 
  ggtitle("GIMAP and IAP Transcript Differential Expression in C. gigas under Bacterial Challenge") + 
  scale_fill_manual("Gene Family", values=(c("GIMAP"="#56B4E9", "IAP"="#009E73"))) +
  theme(plot.title = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5))

COMBINED_GIMAP_IAP_BAC_PLOT2 <- COMBINED_GIMAP_IAP_BAC_PLOT + geom_text(data=label.Bac.df, aes(x=Transcript, y=log2FoldChange), label = c("*")) + geom_hline(yintercept=0) +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  scale_x_discrete(name="Transcript Names (ENSEMBL Transcript ID)", limits=c("transcript:EKC24074", "transcript:EKC18369", "transcript:EKC20239", "transcript:EKC38724", "transcript:EKC42449",
                                                                             "transcript:EKC42441", "transcript:EKC25493", "transcript:EKC34022", "transcript:EKC42442", "transcript:EKC18368", "transcript:EKC41181",
                                                                             "transcript:EKC41180",  "transcript:EKC17690",  "transcript:EKC26950", "transcript:EKC26454", "transcript:EKC40820", "transcript:EKC41832",
                                                                             "transcript:EKC41613", "transcript:EKC35292", "transcript:EKC36405", "transcript:EKC40465", "transcript:EKC38639", "transcript:EKC30713", "transcript:EKC42724"),
                   labels=c("transcript:EKC24074"="BIR protein 2 (EKC24074)","transcript:EKC18369"="BIR protein 3 (EKC18369)","transcript:EKC20239"="BIRprotein 6 (EKC20239)",
                            "transcript:EKC38724"="BIR protein 7 (EKC38724)","transcript:EKC42449"="BIR protein 7-A (EKC42449)","transcript:EKC42441"="BIR protein 7-B(EKC42441)",
                            "transcript:EKC25493"="BIR protein 7-B (EKC25493)","transcript:EKC34022"="IAP 1 (EKC34022)", "transcript:EKC42442"="IAP 2 (EKC42442)",
                            "transcript:EKC18368"="IAP 2 (EKC18368)","transcript:EKC41181"="IAP 2 (EKC41181)", "transcript:EKC41180"="IAP 2 (EKC41180)", 
                            "transcript:EKC17690"="Putative IAP (EKC17690)","transcript:EKC26950"="Putative IAP (EKC26950)","transcript:EKC26454"="Putative IAP ORF42 (EKC26454)",
                            "transcript:EKC40820"="GIMAP 4 (EKC40820)", "transcript:EKC41832"="GIMAP 4 (EKC41832)", "transcript:EKC41613"="GIMAP 4 (EKC41613)",
                            "transcript:EKC35292"="GIMAP 4 (EKC35292)","transcript:EKC36405"="GIMAP 4 (EKC36405)","transcript:EKC40465"="GIMAP 4 (EKC40465)",
                            "transcript:EKC38639"="GIMAP 4 (EKC38639)", "transcript:EKC30713"="GIMAP 7 (EKC30713)","transcript:EKC42724"="GIMAP 7 (EKC42724)"))

#Plot GIMAP and IAP from OsHV1
#Create vector for the significance values! 
label.oshv1Sig.df <- data.frame(Transcript= c("transcript:EKC240741","transcript:EKC424491", "transcript:EKC20774",
                                              "transcript:EKC418321", "transcript:EKC307131"), log2FoldChange=c(-3.0, -3.0, -24.0, -4.0, -3.0))

COMBINED_GIMAP_IAP_OsHV1_PLOT <- ggplot(COMBINED_GIMAP_IAP_cols_Oshv1) + 
  geom_col(aes(x=Transcript, y=log2FoldChange, fill=as.factor(Type))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  scale_y_continuous(name ="log2 Fold Change", breaks = scales::pretty_breaks(n = 20)) + 
  ggtitle("GIMAP and IAP Transcript Differential Expression in C. gigas under OsHV1 Challenge") + 
  scale_fill_manual("Gene Family", values=(c("GIMAP"="#56B4E9", "IAP"="#009E73"))) +
  theme(plot.title = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5))

COMBINED_GIMAP_IAP_OsHV1_PLOT_2 <- COMBINED_GIMAP_IAP_OsHV1_PLOT + geom_text(data=label.oshv1Sig.df, aes(x=Transcript, y=log2FoldChange), label = c("*")) + geom_hline(yintercept=0) +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  scale_x_discrete(name="Transcript Names (ENSEMBL Transcript ID)", limits=c("transcript:EKC240741","transcript:EKC37539", "transcript:EKC42443", "transcript:EKC30031","transcript:EKC20773", 
                                                                             "transcript:EKC183691","transcript:EKC24792", "transcript:EKC202391","transcript:EKC387241","transcript:EKC424491", "transcript:EKC32934","transcript:EKC34718", 
                                                                             "transcript:EKC424411", "transcript:EKC254931", "transcript:EKC20774", "transcript:EKC411811", "transcript:EKC29824", "transcript:EKC340221",
                                                                             "transcript:EKC25955",  "transcript:EKC183681", "transcript:EKC411801", "transcript:EKC424421","transcript:EKC33184", "transcript:EKC34720", 
                                                                             "transcript:EKC264541", "transcript:EKC269501", "transcript:EKC176901","transcript:EKC31739", "transcript:EKC418321",  "transcript:EKC39736",
                                                                             "transcript:EKC404651","transcript:EKC408201", "transcript:EKC32489", "transcript:EKC364051", "transcript:EKC39748", "transcript:EKC27363" ,
                                                                             "transcript:EKC386391", "transcript:EKC416131", "transcript:EKC352921","transcript:EKC307131", "transcript:EKC427241", "transcript:EKC29604"), 
                   labels=c( "transcript:EKC240741"="BIR protein 2 (EKC240741)", "transcript:EKC37539"="BIR protein 2 (EKC37539)",
                             "transcript:EKC42443"="BIR protein 3 (Fragment) (EKC42443)","transcript:EKC30031"="BIR protein 3 (EKC30031)",
                             "transcript:EKC20773"="BIR protein 3 (Fragment) (EKC20773)","transcript:EKC183691"="BIR protein 3 (EKC183691)",
                             "transcript:EKC24792"="BIR protein 5 (EKC24792)","transcript:EKC202391"="BIR protein 6 (EKC202391)",
                             "transcript:EKC387241"="BIR protein 7 (EKC387241)","transcript:EKC424491"="BIR protein 7-A (EKC424491)",
                             "transcript:EKC32934"="BIR protein 7-B (EKC32934)","transcript:EKC34718"="BIR protein 7-B (EKC34718)",
                             "transcript:EKC424411"="BIR protein 7-B (EKC424411)", "transcript:EKC254931"="BIR protein 7-B (EKC254931)",
                             "transcript:EKC20774"="IAP(EKC20774)","transcript:EKC411811"="IAP 1 (EKC411811)",
                             "transcript:EKC29824"="IAP 1 (EKC29824)", "transcript:EKC340221"="IAP 1 (EKC340221)",
                             "transcript:EKC25955"="TP53-regulated IAP 1 (EKC25955)","transcript:EKC183681"="IAP 2 (EKC183681)",  
                             "transcript:EKC411801"="IAP 2 (EKC411801)", "transcript:EKC424421"="IAP 2 (EKC424421)",
                             "transcript:EKC33184"="IAP 3 (EKC33184)", "transcript:EKC34720"="Putative IAP (EKC34720)",
                             "transcript:EKC264541"="Putative IAP ORF42 (EKC264541)","transcript:EKC269501"="Putative IAP (EKC269501)",
                             "transcript:EKC176901"="Putative IAP (EKC176901)", "transcript:EKC31739"="GIMAP 1 (EKC31739)",
                             "transcript:EKC418321"="GIMAP 4 (EKC418321)","transcript:EKC39736"="GIMAP 4 (EKC39736)",
                             "transcript:EKC404651"="GIMAP 4 (EKC404651)","transcript:EKC408201"="GIMAP 4 (EKC408201)", 
                             "transcript:EKC32489"="GIMAP 4 (EKC32489)","transcript:EKC364051"="GIMAP 4 (EKC364051)",
                             "transcript:EKC39748"="GIMAP 4 (EKC39748)", "transcript:EKC27363"="GIMAP 4 (EKC27363)",
                             "transcript:EKC386391"="GIMAP 4 (EKC386391)","transcript:EKC416131"="GIMAP 4 (EKC416131)",
                             "transcript:EKC352921"="GIMAP 4 (EKC352921)","transcript:EKC307131"="GIMAP 7 (EKC307131)", "transcript:EKC427241"="GIMAP 7 (EKC427241)",
                             "transcript:EKC29604"="GIMAP 8 (EKC29604)"))


#Plot IAP between OsHV1 and Bac
IAP_Bac <- COMBINED_GIMAP_IAP_cols_Bac %>% filter(Type== "IAP")
IAP_OsHV1 <- COMBINED_GIMAP_IAP_cols_Oshv1 %>% filter(Type=="IAP")
IAP_OsHV1_Bac <- rbind(IAP_Bac,IAP_OsHV1)
label.IAP.df <- data.frame(Transcript= c("transcript:EKC41180", "transcript:EKC240741", "transcript:EKC424491", "transcript:EKC20774"), log2FoldChange=c(-22.0, -3.0,-3.0, -24.5))

IAP_OsHV1_Bac_PLOT <- ggplot(IAP_OsHV1_Bac) + 
  geom_col(aes(x=Transcript, y=log2FoldChange, fill=as.factor(Challenge))) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  scale_y_continuous(name ="log2 Fold Change", breaks = scales::pretty_breaks(n = 20)) + 
  ggtitle("IAP Transcript Differential Expression in C. gigas under OsHV1 and Bacterial Challenge") + 
  scale_fill_manual("Challenge", values=c("OsHV1"="#0072B2", "Bacteria"= "#D55E00")) +
  theme(plot.title = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5))

IAP_OsHV1_Bac_PLOT + geom_text(data=label.IAP.df, aes(x=Transcript, y=log2FoldChange), label = c("*")) + geom_hline(yintercept=0) +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  scale_x_discrete(name="Transcript Names (ENSEMBL Transcript ID)", limits=c(
    "transcript:EKC24074", "transcript:EKC240741","transcript:EKC37539", "transcript:EKC18369",
    "transcript:EKC183691", "transcript:EKC30031","transcript:EKC20773","transcript:EKC42443", "transcript:EKC24792",
    "transcript:EKC20239","transcript:EKC202391","transcript:EKC38724","transcript:EKC387241","transcript:EKC42449",
    "transcript:EKC424491","transcript:EKC25493","transcript:EKC254931","transcript:EKC32934","transcript:EKC34718",
    "transcript:EKC424411","transcript:EKC42441","transcript:EKC20774","transcript:EKC29824","transcript:EKC34022",
    "transcript:EKC340221","transcript:EKC411811","transcript:EKC18368","transcript:EKC183681","transcript:EKC41180",
    "transcript:EKC411801","transcript:EKC41181","transcript:EKC424421","transcript:EKC42442","transcript:EKC33184",
    "transcript:EKC17690","transcript:EKC176901","transcript:EKC26950","transcript:EKC269501","transcript:EKC34720",
    "transcript:EKC26454","transcript:EKC264541","transcript:EKC25955"), labels=c(
      "transcript:EKC24074" = "BIR protein 2 (EKC24074)", "transcript:EKC240741"="BIR protein 2 (EKC240741)","transcript:EKC37539"="BIR protein 2 (EKC37539)", "transcript:EKC18369"="BIR protein 3 (EKC18369)",
      "transcript:EKC183691"="BIR protein 3 (EKC183691)", "transcript:EKC30031"="BIR protein 3 (EKC30031)","transcript:EKC20773"="BIR protein 3 (Fragment) (EKC20773)","transcript:EKC42443"="BIR protein 3 (Fragment) (EKC42443)", "transcript:EKC24792"="BIR protein 5 (EKC24792)",
      "transcript:EKC20239"="BIR protein 6 (EKC20239)","transcript:EKC202391"="BIR protein 6 (EKC202391)","transcript:EKC38724"="BIR protein 7 (EKC38724)","transcript:EKC387241"="BIR protein 7 (EKC387241)","transcript:EKC42449"="BIR protein 7-A (EKC42449)",
      "transcript:EKC424491"="BIR protein 7-A (EKC424491)","transcript:EKC25493"="BIR protein 7-B (EKC25493)","transcript:EKC254931"="BIR protein 7-B (EKC254931)","transcript:EKC32934"="BIR protein 7-B (EKC32934)","transcript:EKC34718"="BIR protein 7-B (EKC34718)",
      "transcript:EKC424411"="BIR protein 7-B (EKC424411)","transcript:EKC42441"="BIR protein 7-B(EKC42441)","transcript:EKC20774"="IAP (EKC20774)","transcript:EKC29824"="IAP 1 (EKC29824)","transcript:EKC34022"="IAP 1 (EKC34022)",
      "transcript:EKC340221"="IAP 1 (EKC340221)","transcript:EKC411811"="IAP 1 (EKC411811)","transcript:EKC18368"="IAP 2 (EKC18368)","transcript:EKC183681"="IAP 2 (EKC183681)","transcript:EKC41180"="IAP 2 (EKC41180)",
      "transcript:EKC411801"="IAP 2 (EKC411801)","transcript:EKC41181"="IAP 2 (EKC41181)","transcript:EKC424421"="IAP 2 (EKC424421)","transcript:EKC42442"="IAP 2 (EKC42442)","transcript:EKC33184"="IAP 3 (EKC33184)",
      "transcript:EKC17690"="Putative IAP (EKC17690)","transcript:EKC176901"="Putative IAP (EKC176901)","transcript:EKC26950"="Putative IAP (EKC26950)","transcript:EKC269501"="Putative IAP (EKC269501)","transcript:EKC34720"="Putative IAP (EKC34720)",
      "transcript:EKC26454"="Putative IAP ORF42 (EKC26454)","transcript:EKC264541"="Putative IAP ORF42 (EKC264541)","transcript:EKC25955"="TP53-regulated IAP 1 (EKC25955)"))

#plot IAPs with the same transcript side by side
#Establish which transcripts are shared and which are different
duplicated(IAP_OsHV1_Bac$Transcript) #none are duplicated!

#plot GIMAP between OsHV1 and Bac
GIMAP_Bac <- COMBINED_GIMAP_IAP_cols_Bac %>% filter(Type== "GIMAP")
GIMAP_OsHV1 <- COMBINED_GIMAP_IAP_cols_Oshv1 %>% filter(Type=="GIMAP")
GIMAP_OsHV1_Bac <- rbind(GIMAP_Bac,GIMAP_OsHV1)
label.GIMAP.df <- data.frame(Transcript= c("transcript:EKC418321","transcript:EKC307131"), log2FoldChange=c(-3.5, -2.5))

GIMAP_OsHV1_Bac_PLOT <- ggplot(GIMAP_OsHV1_Bac) + 
  geom_col(aes(x=Transcript, y=log2FoldChange, fill=as.factor(Challenge))) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  scale_y_continuous(name ="log2 Fold Change", breaks = scales::pretty_breaks(n = 20)) + 
  ggtitle("GIMAP Transcript Differential Expression in C. gigas under OsHV1 and Bacterial Challenge") + 
  scale_fill_manual("Challenge", values=c("Bacteria"="#0072B2", "OsHV1"= "#D55E00")) +
  theme(plot.title = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5)) 

GIMAP_OsHV1_Bac_PLOT2 <- GIMAP_OsHV1_Bac_PLOT + geom_text(data=label.GIMAP.df, aes(x=Transcript, y=log2FoldChange), label = c("*")) + geom_hline(yintercept=0) +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  scale_x_discrete(name="Transcript Names (ENSEMBL Transcript ID)", limits=c(
    "transcript:EKC31739", "transcript:EKC40820", "transcript:EKC41832","transcript:EKC41613", 
    "transcript:EKC35292", "transcript:EKC36405", "transcript:EKC40465", "transcript:EKC38639",
    "transcript:EKC39736", "transcript:EKC404651", "transcript:EKC408201","transcript:EKC32489",
    "transcript:EKC364051","transcript:EKC39748","transcript:EKC27363","transcript:EKC352921",
    "transcript:EKC416131","transcript:EKC386391","transcript:EKC418321","transcript:EKC30713",
    "transcript:EKC42724", "transcript:EKC427241","transcript:EKC307131","transcript:EKC29604"), 
    labels=c("transcript:EKC31739"="GIMAP 1 (EKC31739)", "transcript:EKC40820"="GIMAP 4 (EKC40820)",
             "transcript:EKC41832"="GIMAP 4 (EKC41832)", "transcript:EKC41613"="GIMAP 4 (EKC41613)",
             "transcript:EKC35292"="GIMAP 4 (EKC35292)","transcript:EKC36405"="GIMAP 4 (EKC36405)",
             "transcript:EKC40465"="GIMAP 4 (EKC40465)","transcript:EKC38639"="GIMAP 4 (EKC38639)",
             "transcript:EKC39736"="GIMAP 4 (EKC39736)", "transcript:EKC404651"="GIMAP 4 (EKC404651)",
             "transcript:EKC408201"="GIMAP 4 (EKC408201)", "transcript:EKC32489"="GIMAP 4 (EKC32489)",
             "transcript:EKC364051"="GIMAP 4 (EKC364051)", "transcript:EKC39748"="GIMAP 4 (EKC39748)",
             "transcript:EKC27363"="GIMAP 4 (EKC27363)", "transcript:EKC352921"="GIMAP 4 (EKC352921)",
             "transcript:EKC416131"="GIMAP 4 (EKC416131)", "transcript:EKC386391"="GIMAP 4 (EKC386391)",
             "transcript:EKC418321"="GIMAP 4 (EKC418321)","transcript:EKC30713"="GIMAP 7 (EKC30713)",
             "transcript:EKC42724"="GIMAP 7 (EKC42724)","transcript:EKC427241"="GIMAP 7 (EKC427241)",
             "transcript:EKC307131"="GIMAP 7 (EKC307131)","transcript:EKC29604"="GIMAP 8 (EKC29604)"))

#plot GIMAPs with same accession 
duplicated(GIMAP_OsHV1_Bac$Transcript) 
#none of them have been duplicated from trial to trial

#references: 
#https://github.com/sr320/LabDocs/tree/master/code/DESeq
#StringTie manual : http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
