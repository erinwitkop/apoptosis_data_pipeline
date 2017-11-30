#05_C_Vir_Bac_DESeq2_transcript_analysis.R

####INPUT GATHERED FROM CRASSOSTREA VIRGINICA ####

#This script takes as input the output C_Vir_transcript_count_matrix.csv data prepared from prepDE.py and performs
#differential transcript expression analysis, and subsets out isoforms of GIMAPs and IAPs and graphs their
#relative abundance.

#call the DESeq2 library 
source("http://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade") 
biocLite("DESeq2")
library("DESeq2")
#install.packages("fdrtool")
library(fdrtool)
#source("https://bioconductor.org/biocLite.R")
#install.packages("dplyr")
library(dplyr)
#install.packages("tidyr")
library(tidyr)
#install.packages("reshape2")
library(reshape2)
library(ggplot2)
# Construct C_Vir_PHENO_DATA.csv that contains SRA run information, such as which contrast, tissue, etc.

####Match "rna#" value with the Gene LOC name in the stringtie file
#Read in stringtie merged file
C_vir_stringtie <- read.csv(file= "C_Vir_stringtie_merged_RERUN.gtf", sep="\t", header = FALSE)
#set colnames for the attributes
colnames(C_vir_stringtie) <- c('seqname', 'source', 'feature','start','end','score','strand','frame', 'attribute')
#subset for just transcript lines
C_vir_stringtie_transcripts <- C_vir_stringtie %>% filter(feature =="transcript")
#grep out only the lines with "rna" in the line
C_vir_stringtie_transcripts_rna <-filter(C_vir_stringtie_transcripts, grepl("rna", attribute))
#grep out only the lines with "rna" and a "LOC" value
C_vir_stringtie_transcripts_rna_LOC <- filter(C_vir_stringtie_transcripts_rna, grepl("LOC", attribute))
#make into df
C_vir_stringtie_transcripts_rna_LOC <- as.data.frame(C_vir_stringtie_transcripts_rna_LOC)
#separate attributes column with the ; separator in tidyr
#use the extra  argument to control what happens when every row doesn't split into the same number of pieces
C_vir_stringtie_transcripts_rna_LOC_separated <- C_vir_stringtie_transcripts_rna_LOC %>% separate(attribute, c("gene_id","transcript_id","gene_name","ref_gene_id"), ";", extra="merge")   
#clean up white space in the data
C_vir_stringtie_transcripts_rna_LOC_separated_clean <- data.frame(lapply(C_vir_stringtie_transcripts_rna_LOC_separated, trimws),stringsAsFactors = FALSE)
#remove spaces in terminal
# tr ' ' ';'<C_vir_stringtie_transcripts_rna_LOC_separated_clean.csv >C_vir_stringtie_transcripts_rna_LOC_separated_clean2.csv
write.csv(C_vir_stringtie_transcripts_rna_LOC_separated_clean, file ="C_vir_stringtie_transcripts_rna_LOC_separated_clean.csv")

####DEG Analysis with TRANSCRIPT Count Matrix ####
#load transcript count matrix and labels
#C_VirPHENO_DATA.csv file contains metadata on the count table's samples
###Make sure excel PHENODATA is in the same order or these commands will change data to be wrong!!!!###

#Extract correct rows and columns from the PHENO DATA and transcript data
C_vir_TranscriptCountData <- as.data.frame(read.csv("C_vir_RERUN_transcript_count_matrix.csv", row.names="transcript_id"))
head(C_vir_TranscriptCountData)
#duplicate control columns for CGX to be used with the ROD resistant family
#SRR1293904_2, SRR1298417_2, SRR1298710_2
C_vir_TranscriptCountData$SRR1293904_2 = C_vir_TranscriptCountData$SRR1293904
C_vir_TranscriptCountData$SRR1298417_2 = C_vir_TranscriptCountData$SRR1298417
C_vir_TranscriptCountData$SRR1298710_2 = C_vir_TranscriptCountData$SRR1298710

#load Phenotype data with addded column for extra control
#upload pheno data in the correct order with the controls of A and B listed first
#then change the order of the TranscriptCountData to match this order 
C_vir_TranColData <- read.csv("C_vir_PHENO_DATA_edited_reordered.csv", header=TRUE, sep=",")
print(C_vir_TranColData)
#reorder TranscriptCountData
C_vir_TranscriptCountData <- C_vir_TranscriptCountData[c(1,3,4,11, 5,6,7,12,19,20,21,
                                                         8,9,10,2,13,14,15,16,17,18)]
head(C_vir_TranscriptCountData)
#change rownames to match
rownames(C_vir_TranColData) <- C_vir_TranColData$sampleID
colnames(C_vir_TranscriptCountData) <- C_vir_TranColData$sampleID
head(C_vir_TranColData)
head(C_vir_TranscriptCountData)
# Check all sample IDs in C_vir_TranColData are also in C_vir_TranscriptCountData and match their orders
all(rownames(C_vir_TranColData) %in% colnames(C_vir_TranscriptCountData))  #Should return TRUE
# returns TRUE
all(rownames(C_vir_TranColData) == colnames(C_vir_TranscriptCountData))    # should return TRUE
#returns TRUE

####Subset Data for ROD Transcriptomes ####
#Extract columns you want from the C_vir_TranscriptCountData, based on which column the correct SRA data
#NOT USING DAY 30 BECAUSE IN PROBIOTIC CHALLENGE IT WAS 5, 12, AND 16 DAYS
ROD_C_vir_TranColData <- C_vir_TranColData[c(1,2,4,5,6,8,9,10,11,12,13,14),]
ROD_C_vir_TranscriptCountData <- C_vir_TranscriptCountData[,c(1,2,4,5,6,8,9,10,11,12,13,14),]
ROD_C_vir_TranColData
head(ROD_C_vir_TranscriptCountData)

#Give the condition column levels
ROD_C_vir_TranColData$condition <- factor(ROD_C_vir_TranColData$condition)
levels(ROD_C_vir_TranColData$condition) #check to see that it has levels 

#give the treatment column levels
ROD_C_vir_TranColData$treatment <- factor(ROD_C_vir_TranColData$treatment)
levels(ROD_C_vir_TranColData$treatment)

#### Create ROD DESeq Data Set from Matrix ####
# DESeqDataSet from count matrix and labels, separate into resistant and susceptible 
#add an interaction term to compare treatment between two conditions 
#layout used for interactions: https://support.bioconductor.org/p/58162/

ddsRODTran <- DESeqDataSetFromMatrix(countData = ROD_C_vir_TranscriptCountData, 
                                       colData = ROD_C_vir_TranColData, 
                                       design =  ~ condition + treatment + condition:treatment)

ddsRODTran<- ddsRODTran[ rowSums(counts(ddsRODTran)) > 1, ]

# review how the data set looks
head(ddsRODTran)

#Relevel each to make sure that control is the first level in the treatment factor for each
ddsRODTran$condition <- relevel(ddsRODTran$condition, "A")
ddsRODTran$treatment <- relevel( ddsRODTran$treatment, "control")

#Check we're looking at the right samples
as.data.frame( colData(ddsRODTran) )

#Running the DEG pipeline
ddsRODTran<- DESeq(ddsRODTran, betaPrior = FALSE) #for designs with interactions, recommends setting betaPrior=FALSE

#Inspect results
#extract contrasts between control and treatment values for interaction
resRODTran<- results(ddsRODTran)
head(resRODTran)

# this last line is all that is needed because the interaction term is last
  # in the design formula (said my Michael Love)
#Small p-values for the interaction term indicate that the log fold change
  # due to treatment is significantly different for the two conditions.

#summary is just printing a table for you, you need to tell it what threshold you want
help("summary",package="DESeq2")
alpha <- 0.05 #set alpha to 0.05, this will control FDR
summary(resRODTran) #default FDR is still 0.1
summary(resRODTran, alpha) #no showing all genes with FRD < 0.05

#To get the significant genes
#The independent filtering in results() has an argument 'alpha'
#which is used to optimize a cutoff on mean normalized count
#to maximize the number of genes with padj < alpha
resRODTran_05 <- results(ddsRODTran, alpha= alpha) #set FDR to 0.05 now
resRODTran_05_Sig <- resRODTran[which(resRODTran$padj < alpha),]
summary(resRODTran_05) #this is all the genes
summary(resRODTran_05_Sig) #this is the significant ones!
sum(resRODTran_05$padj < 0.05, na.rm=TRUE) #4121 tells you how many genes have expected FDR ≤ 0.05
sum(resRODTran_05_Sig$padj < 0.05, na.rm=TRUE) #4102, differ by 19 genes only 
sig="significant"
resRODTran_05_Sig$Significance <- sig
resRODTran_05_nonSig <- resRODTran[which(resRODTran$padj > alpha),] #create list of nonsig
nonsig <- "non-significant"
resRODTran_05_nonSig$Significance <- nonsig

#add ID column with the rownames so a merge can happen later
resRODTran_05_Sig["ID"] <- rownames(resRODTran_05_Sig) #add a new column with the rownames for match
resRODTran_05_nonSig["ID"] <- rownames(resRODTran_05_nonSig)
resRODTran_05["ID"] <- rownames(resRODTran_05)

#metadata on meaning of the columns
mcols(resRODTran_05_Sig, use.names = TRUE)
mcols(resRODTran_05_Sig)$description
#shows treatment vs. control with baseMean as the mean of normalized counts for all samples

#Order by Log2FC
head( resRODTran_05[ order( resRODTran_05$log2FoldChange ), ] ) #head for strongest downregulation
tail( resRODTran_05[ order( resRODTran_05$log2FoldChange ), ] ) #tail for strongest up regulation

####p-value correction for both sets of results ####
#First Visualize histograms
#histogram of P- values to visualize any "hills" or "U shape"
# hill means variance of the null distribution too high, U shape means variance assumed too low
hist(resRODTran_05$pvalue, breaks= 20, col = "grey")
hist(resRODTran_05_Sig$pvalue, breaks = 20, col = "grey") #hill

#looks good, don't need to do the correction

#Visualize Results with Diagnostic Plots#
#MA plot, useful overview for experiment with two-group comparison. Plots log2FC over mean of normalized counts
#genes with adjusted p value 
plotMA(resRODTran_05)
plotMA(resRODTran_05_Sig)

#Export Results to CSV
write.csv( as.data.frame(resRODTran_05), file="ROD_resRODTran_05.csv")
write.csv( as.data.frame(resRODTran_05_Sig), file="ROD_resRODTran_05_Sig.csv")
write.csv( as.data.frame(resRODTran_05_nonSig), file = "ROD_resRODTran_05_nonSig.csv")

####Match lines for only those that have a "LOC" in the C_vir_stringtie_transcripts_rna_LOC_separated$transcript_id ####
# %in% returns true for every value in the first argument that matches a value in the second argument.
# the order of arguments is important
#load file with the spaces changed to ; so that the column can be separated
C_vir_stringtie_transcripts_rna_LOC_separated_clean2 <- read.csv(file="C_vir_stringtie_transcripts_rna_LOC_separated_clean2.csv", header=TRUE)

#separate the transcript id column with the semicolon
C_vir_stringtie_transcripts_rna_LOC_separated_clean_transcript_id <- C_vir_stringtie_transcripts_rna_LOC_separated_clean2 %>% separate(transcript_id, c("transcript_id", "ID"), ";", extra="merge")
#check duplicate rownames
ROD_dupliate_names<- duplicated(resRODTran_05_Sig$rownames)
grep("TRUE", ROD_dupliate_names) #0 duplicates
#must comvert to data frame
resRODTran_05_Sig_df <- data.frame(resRODTran_05_Sig)
#Merge columns together based on match in the "ID" column 
resRODTran_05_Sig_df_FULL <- merge(resRODTran_05_Sig_df, C_vir_stringtie_transcripts_rna_LOC_separated_clean_transcript_id[,c("ID", "gene_name")], by="ID")
nrow(resRODTran_05_Sig_df_FULL) #1854
#put LOC info in new column
resRODTran_05_Sig_df_FULL <- resRODTran_05_Sig_df_FULL %>% separate(gene_name, c("gene_name", "gene_ID"), ";")

#Repeat for non significant genes
resRODTran_05_nonSig["ID"] <- rownames(resRODTran_05_nonSig) 
resRODTran_05_df_non_Sig_df <- data.frame(resRODTran_05_df_non_Sig)
resRODTran_05_df_non_Sig_FULL <- merge(resRODTran_05_df_non_Sig_df, C_vir_stringtie_transcripts_rna_LOC_separated_clean_transcript_id[,c("ID", "gene_name")], by="ID")
#put LOC info in new column
resRODTran_05_df_non_Sig_FULL <- resRODTran_05_df_non_Sig_FULL %>% separate(gene_name, c("gene_name", "gene_ID"), ";")

####STOPPED EDITING HERE####

####Lookup LOC values using Batch Entrez ####
#dfSig
write.csv(resRODTran_05_dfSig_FULL$gene_ID, file="resRODTran_05_dfSig_FULL_gene_ID.csv")
#put column with LOC info into text file
#perform batch Entrez lookup to the gene database to retrieve IDs
#df_non_Sig
write.csv(resRODTran_05_df_non_Sig_FULL$gene_ID, file="resRODTran_05_df_non_Sig_FULL_gene_ID.csv")
#put column wit LOC into text file
#perform batch Entrez lookup with the gene database to retrieve IDs
#batch Entrez gets rid of duplicates

#####Merge the LOC values with gene name####
#SIG
#convert file to csv using excel 
resRODTran_05_Sig_ENTREZGENE <- read.csv(file="resRODTran_05_df_non_Sig_ENTREZ_RESULTS_FULL.csv", header = TRUE)
#change symbol to gene_ID
colnames(resRODTran_05_Sig_ENTREZGENE)[6] <- "gene_ID"
#merge file based on the LOC column
resRODTran_05_Sig_ENTREZGENE_DATA <-  merge(resRODTran_05_Sig_ENTREZGENE, resRODTran_05_dfSig_FULL[,c("gene_ID", "log2FoldChange","pvalue","padj")], by="gene_ID")

####Extract GIMAP/IAN proteins and CgIAPs from Significantly Differentially Expressed Genes####
#Significantly differentially Expressed IAPs, using several names
ROD_IAPs_SIG <- grepl("IAP", resRODTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", ROD_IAPs_SIG) #184, 1037, 1458, 1459 
ROD_IAPs_SIG_info <- resRODTran_05_Sig_ENTREZGENE_DATA[c(184, 1037, 1458, 1459),]
ROD_IAPs_SIG_info

#Sig apoptosis
ROD_apoptosis_SIG <- grepl("apoptosis", resRODTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE", ROD_apoptosis_SIG) #0


#Significant GIMAP Genes
ROD_IAN_Sig <- grepl("IAN", resRODTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", ROD_IAN_Sig) #346 347 378 379 384 434...NONE
ROD_AIG_Sig <- grepl("AIG", resRODTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE", ROD_AIG_Sig) #0 TRUE
ROD_IMAP_Sig <- grepl("IMAP", resRODTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", ROD_IMAP_Sig) #404 TRUE!

ROD_GTP_SIG <- grepl("GTP", resRODTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", ROD_GTP_SIG) # 9  130  144  150  151  181  233  234  235  283  284  285  286  404  460  461  495  649  650  726  727
# 728  912  913  914 1423 1424 1889 1890
resRODTran_05_Sig_ENTREZGENE_DATA[c(9,  130,  144,  150,  151,  181,  233,  234,  235,
                                    283,  284,  285,  286,  404,  460,  461,  495,  649,
                                    650,  726,  727,
                                    728,  912,  913,  914, 1423, 1424, 1889, 1890),]
ROD_GIMAP_Sig_info <- resRODTran_05_Sig_ENTREZGENE_DATA[404,]

####NON Sig from Entrez for ROD challenge ####
resRODTran_05_df_non_Sig_ENTREZ_RESULTS_FULL
#convert file to csv using excel 
resRODTran_05_non_Sig_ENTREZGENE <- read.csv(file="/data/resRODTran_05_df_non_Sig_ENTREZ_RESULTS_FULL.csv", header = TRUE)
#change symbol to gene_ID
colnames(resRODTran_05_non_Sig_ENTREZGENE)[6] <- "gene_ID"
#merge file based on the LOC column
resRODTran_05_non_Sig_ENTREZGENE_DATA <-  merge(resRODTran_05_non_Sig_ENTREZGENE, resRODTran_05_df_non_Sig_FULL[,c("gene_ID", "log2FoldChange","pvalue","padj")], by="gene_ID")
#remove duplicates
resRODTran_05_non_Sig_ENTREZGENE_DATA[!duplicated(resRODTran_05_non_Sig_ENTREZGENE_DATA),]

####Extract GIMAP/IAN proteins and CgIAPs from non Sig Expressed Genes####
#Non sig IAPS
ROD_IAPs_non_SIG <- grepl("IAP", resRODTran_05_non_Sig_ENTREZGENE_DATA$description , ignore.case = TRUE) 
grep("TRUE", ROD_IAPs_non_SIG) #  947,   948,   949,   950,   959,   960 ,  978,   986,   987,  1005,  1006,  1007,  1008,  1182,  1448,
1469,  1470 , 1984,  1985 , 2087 , 2526,  2676 , 2817,  2820 , 2891  ,2892 , 3011,  3053  ,3054,  3055 , 3318,
3341,  3342 , 3495 , 3721  , 3767 , 3768 , 3901 , 3902 , 3903,  4295  ,4300,  4496 , 4497 ,13132,
17162,22227, 22228, 23039, 28602 
ROD_IAPs_non_Sig_info <- resRODTran_05_non_Sig_ENTREZGENE_DATA[c(947,   948,   949,   950,   959,   960 ,  978,   986,   987,  1005,  1006,  1007,  1008,  1182,  1448,
                                                                 1469,  1470 , 1984,  1985 , 2087 , 2526,  2676 , 2817,  2820 , 2891  ,2892 , 3011,  3053  ,3054,  3055 , 3318,
                                                                 3341,  3342 , 3495 , 3721  , 3767 , 3768 , 3901 , 3902 , 3903,  4295  ,4300,  4496 , 4497 ,13132,
                                                                 17162,22227, 22228, 23039, 28602),]

ROD_apoptosis_non_SIG <- grepl("apoptosis", resRODTran_05_non_Sig_ENTREZGENE_DATA$description , ignore.case = TRUE) 
grep("TRUE", ROD_apoptosis_non_SIG) #294  1359  1360  1361  1446  1978  2324  2762  3284  3671  4059 13035 13036 13534 13725 14788 15726
[18] 15727 15867 15868 15869 16416 20870 21078 21812 22539 24495 25340 28795 29499
ROD_apoptosis_non_SIG_info <- resRODTran_05_non_Sig_ENTREZGENE_DATA[c(294,  1359,  1360,  1361,  1446,  1978,  2324,2762,  
  3284,  3671,  4059, 13035, 13036, 13534, 13725, 14788, 15726,
  15727, 15867, 15868, 15869, 16416, 20870, 21078, 21812, 22539, 24495, 25340, 28795, 29499),]
ROD_inhibitors_of_apoptosis_info <- resRODTran_05_non_Sig_ENTREZGENE_DATA[c(1359,1360,1361,1446,1978, 2324,3284,4059,25340),]

#inhibitors of apoptosis: 1359,1360,1361,1446,1978, 2324,3284,4059,25340


#Non sig GIMAP Genes
ROD_IAN_non_Sig <- grepl("IAN", resRODTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", ROD_IAN_non_Sig) # 635  3371  3372  5463  5464  5926  5997  5998  6501  6595 10469 13278 13341 13342 13494 13596 14076
[18] 14077 16209 16837 17376 18323 20068 21219 21374 23360 23361 23749 23750 23963 23964 24418 25505 25506
[35] 25507 25508 25509 27470 28944 28971 28972 29801 29811 29898 30333 30334 30335 30336
ROD_IAN_non_Sig_info <- resRODTran_05_non_Sig_ENTREZGENE_DATA[c(635,  3371,  3372,  5463,  5464,  5926,  5997,  5998,  6501,  6595, 10469, 13278, 13341, 13342, 13494, 13596, 14076,
                                                                14077, 16209, 16837, 17376, 18323, 20068, 21219, 21374, 23360, 23361, 23749, 23750, 23963, 23964, 24418, 25505, 25506,
                                                                25507, 25508, 25509, 27470, 28944, 28971, 28972, 29801, 29811 ,29898, 30333, 30334, 30335, 30336),]
#None

ROD_AIG_non_Sig <- grepl("AIG", resRODTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", ROD_AIG_non_Sig) #0 TRUE
ROD_IMAP_non_Sig <- grepl("IMAP", resRODTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", ROD_IMAP_non_Sig) # 86  2737  2823  3080  3082  4398  4399  4400  4401  4402  4630  4767  4858  4861  5014  5021  5461
  #6167  6168  6169  6265  6299  6538  6708  6882  7032  7122  7123  7170  7269  7434  7575  8296  8636
  #13308 13309 13310 13964 22616 22618 22619 22861 22862 2286423
ROD_IMAP_non_Sig_info <- resRODTran_05_non_Sig_ENTREZGENE_DATA[c(2737,  2823,  3080,  3082,  4398,  4399,  4400, 4401,  4402,  4630,  4767,  4858,  4861,  5014,  5021,  5461,
                                                                 6167,  6168,  6169,  6265,  6299,  6538,  6708,  6882,  7032,  7122,  7123,  7170,  7269,  7434,  7575,  8296,  8636,
                                                                 13308, 13309, 13310, 13964, 22616, 22618, 22619, 22861, 22862, 2286423),]

ROD_GTPase_non_SIG <- grepl("GTPase IMAP", resRODTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", ROD_GTPase_non_SIG) # 2386  2737  2823  3080  3082  4398  4399  4400  4401  4402  4630  4767  4858  4861  5014  5021  5461
6167  6168  6169  6265  6299  6538  6708  6882  7032  7122  7123  7170  7269  7434  7575  8296  8636
13308 13309 13310 13964 22616 22618 22619 22861 22862 22864
ROD_GTPase_non_Sig_info <- resRODTran_05_non_Sig_ENTREZGENE_DATA[c(2386,  2737,  2823,  3080,  3082 , 4398 , 4399,  4400 , 4401,  4402,  4630 , 4767,  4858 , 4861,  5014 , 5021,  5461,
                                        6167,  6168,  6169,  6265,  6299,  6538 , 6708 , 6882 , 7032,  7122,  7123,  7170 , 7269,  7434 , 7575,  8296 , 8636,
                                        13308, 13309, 13310, 13964, 22616, 22618, 22619 ,22861, 22862, 22864),] 


#### RI PROBIOTIC CHALLENGED Differential Gene Expression Analysis ####

#Subset RIF Challenge Data
RIFTranCountData <- C_vir_TranCountData[ , c(4,5,6,7,8,9)]
head(RIFTranCountData)
RIFTranColData <- C_vir_TranColData[c(4,5,6,7,8,9),]
rownames(RIFTranColData) <- RIFTranColData$sampleID
colnames(RIFTranCountData) <- RIFTranColData$sampleID
head(RIFTranCountData)
print(RIFTranColData)

# Check all sample IDs in BacColData are also in BacCountData and match their orders
all(rownames(RIFTranColData) %in% colnames(RIFTranCountData))  #Should return TRUE
# returns TRUE
all(rownames(RIFTranColData) == colnames(RIFTranCountData))    # should return TRUE
#returns TRUE

#Give the "treatment" column levels
RIFTranColData$treatment <- factor(RIFTranColData$treatment)
levels(RIFTranColData$treatment) #check to see that it has levels 

#### Create Bacterial Challenge DESeq Data Set from Matrix ####
# DESeqDataSet from count matrix and labels 
ddsRIFTran <- DESeqDataSetFromMatrix(countData = RIFTranCountData, 
                                     colData = RIFTranColData, 
                                     design = ~ treatment)

#this design will gather the effect of condition between control and treatment
# review how the data set looks
head(ddsRIFTran)

#Relevel each to make sure that control is the first level in the treatment factor
ddsRIFTran$treatment <- relevel( ddsRIFTran$treatment, "control")
levels(ddsRIFTran$treatment)

#Check we're looking at the right samples
as.data.frame( colData(ddsRIFTran) )

#Running the DEG pipeline
ddsRIFTran<- DESeq(ddsRIFTran) #for designs with interactions, recommends setting betaPrior=FALSE

#Inspect results table
resultsNames(ddsRIFTran)

#extract contrasts between control and treatment values
resRIFTran <- results(ddsRIFTran)
#to extract log2fold change and p values under 0.1 and 0.05
head(resRIFTran)
summary(resRIFTran)
sum(resRIFTran$padj < 0.1, na.rm=TRUE) #390
#Change alpha setting in DESeq Results
resRIFTran_05 <- results(ddsRIFTran,alpha=0.05)
summary(resRIFTran_05)
sum(resRIFTran_05$padj < 0.05, na.rm=TRUE) #306

#metadata on meaning of the columns
mcols(resRIFTran_05, use.names = TRUE)
#shows treatment vs. control with baseMean as the mean of normalized counts for all samples

####p-value correction#### adapted from "Differential expression analysis of RNA-Seq data using DESeq2" Klaus 2014
#First Visualize histograms
#histogram of P- values to visualize any "hills" or "U shape"
# hill means variance of the null distribution too high, U shape means variance assumed too low
hist(resRIFTran_05$pvalue, breaks = 20, col = "grey") #hill

#remove filtered out genes by independent filtering, they have NA adj. pvals
resRIFTran_05_df <- resRIFTran_05[ !is.na(resRIFTran_05$padj), ]

#remove genes with NA pvals (outliers)
resRIFTran_05_df <- resRIFTran_05_df[ !is.na(resRIFTran_05_df$pvalue), ]

#remove adjsuted pvalues, since we add the fdrtool results later on (based on the correct p-values)
resRIFTran_05_df <- resRIFTran_05_df[, -which(names(resRIFTran_05_df) == "padj")]

#use z-scores as input to FDRtool to re-estimate the p-value
FDR.resRIFTran_05_df <- fdrtool(resRIFTran_05_df$stat, statistic= "normal", plot = T)

#add values to the results data frame, also ad new BH- adjusted p-values
resRIFTran_05_df[,"padj"] <- p.adjust(FDR.resRIFTran_05_df$pval, method = "BH")

#replot corrected p-values
hist(FDR.resRIFTran_05_df$pval, col = "royalblue4",
     main = "Correct null model RIF Challenge Transcript Counts", xlab = "CORRECTED p-values")

#Check how many genes have BH adjusted p values of less than 0.05 after P-value correction?
sum( resRIFTran_05_df$padj < 0.05, na.rm=TRUE ) #2036

#Subset the results table to the differentially expressed genes under FDR 0.01, order the Log2FC table first by strongest down regulation
resRIFTran_05_dfSig <- resRIFTran_05_df[ which(resRIFTran_05_df$padj < 0.05 ), ]
head( resRIFTran_05_dfSig[ order( resRIFTran_05_dfSig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resRIFTran_05_dfSig[ order( resRIFTran_05_dfSig$log2FoldChange ), ] ) #tail for strongest up regulation
resRIFTran_05_df_non_Sig <- resRIFTran_05_df[ which(resRIFTran_05_df$padj > 0.05), ]
resRIFTran_05_df_non_Sig

####Visualize Results with Diagnostic Plots####
#MA plot, useful overview for experiment with two-group comparison. Plots log2FC over mean of normalized counts
#genes with adjusted p value 
plotMA(resRIFTran_05_df)
plotMA(resRIFTran_05_dfSig)
plotMA(resRIFTran_05_df_non_Sig)

#Export Results to CSV
write.csv( as.data.frame(resRIFTran_05_df), file="resRIFTran_05_df.csv")
write.csv( as.data.frame(resRIFTran_05_dfSig), file="resRIFTran_05_dfSig.csv")
write.csv( as.data.frame(resRIFTran_05_df_non_Sig), file="resRIFTran_05_df_non_Sig.csv")

####Match lines for only those that have a "LOC" in the C_vir_stringtie_transcripts_rna_LOC_separated$transcript_id ####
# %in% returns true for every value in the first argument that matches a value in the second argument.
  # the order of arguments is important
#load file with the spaces changed to ; so that the column can be separated
C_vir_stringtie_transcripts_rna_LOC_separated_clean2 <- read.csv(file="C_vir_stringtie_transcripts_rna_LOC_separated_clean2.csv", header=TRUE)

#separate the transcript id column with the semicolon
C_vir_stringtie_transcripts_rna_LOC_separated_clean_transcript_id <- C_vir_stringtie_transcripts_rna_LOC_separated_clean2 %>% separate(transcript_id, c("transcript_id", "ID"), ";", extra="merge")
#add ID column with the rownames so a merge can happen later
resRIFTran_05_dfSig["ID"] <- rownames(resRIFTran_05_dfSig) #add a new column with the rownames for match
#check duplicate rownames
dupliate_names<- duplicated(resRIFTran_05_dfSig$rownames)
grep("TRUE", dupliate_names) #0 duplicates
#must comvert to data frame
resRIFTran_05_dfSig_df <- data.frame(resRIFTran_05_dfSig)
#Merge columns together based on match in the "ID" column 
resRIFTran_05_dfSig_FULL <- merge(resRIFTran_05_dfSig_df, C_vir_stringtie_transcripts_rna_LOC_separated_clean_transcript_id[,c("ID", "gene_name")], by="ID")
nrow(resRIFTran_05_dfSig_FULL) #2031
#put LOC info in new column
resRIFTran_05_dfSig_FULL <- resRIFTran_05_dfSig_FULL %>% separate(gene_name, c("gene_name", "gene_ID"), ";")

#Repeat for non significant genes
resRIFTran_05_df_non_Sig["ID"] <- rownames(resRIFTran_05_df_non_Sig) 
resRIFTran_05_df_non_Sig_df <- data.frame(resRIFTran_05_df_non_Sig)
resRIFTran_05_df_non_Sig_FULL <- merge(resRIFTran_05_df_non_Sig_df, C_vir_stringtie_transcripts_rna_LOC_separated_clean_transcript_id[,c("ID", "gene_name")], by="ID")
#put LOC info in new column
resRIFTran_05_df_non_Sig_FULL <- resRIFTran_05_df_non_Sig_FULL %>% separate(gene_name, c("gene_name", "gene_ID"), ";")

####Lookup LOC values using Batch Entrez ####
#dfSig
write.csv(resRIFTran_05_dfSig_FULL$gene_ID, file="resRIFTran_05_dfSig_FULL_gene_ID.csv")
#put column with LOC info into text file
#perform batch Entrez lookup to the gene database to retrieve IDs
#df_non_Sig
write.csv(resRIFTran_05_df_non_Sig_FULL$gene_ID, file="resRIFTran_05_df_non_Sig_FULL_gene_ID.csv")
#put column wit LOC into text file
#perform batch Entrez lookup with the gene database to retrieve IDs

#####Merge the LOC values with gene name####
#sig
#convert file to csv using excel 
resRIFTran_05_Sig_ENTREZGENE <- read.csv(file="resRIFTran_05_dfSig_ENTREZ_RESULTS.csv", header = TRUE)
#change symbol to gene_ID
colnames(resRIFTran_05_Sig_ENTREZGENE)[6] <- "gene_ID"
#merge file based on the LOC column
resRIFTran_05_Sig_ENTREZGENE_DATA <-  merge(resRIFTran_05_Sig_ENTREZGENE, resRIFTran_05_dfSig_FULL[,c("gene_ID", "log2FoldChange","pvalue","padj")], by="gene_ID")

#Non Sig
#convert file to csv using excel 
resRIFTran_05_non_Sig_ENTREZGENE <- read.csv(file="resRIFTran_05_df_non_Sig_ENTREZ_RESULTS_FULL.csv", header = TRUE)
#change symbol to gene_ID
colnames(resRIFTran_05_non_Sig_ENTREZGENE)[6] <- "gene_ID"
#merge file based on the LOC column
resRIFTran_05_non_Sig_ENTREZGENE_DATA <-  merge(resRIFTran_05_non_Sig_ENTREZGENE, resRIFTran_05_df_non_Sig_FULL[,c("gene_ID", "log2FoldChange","pvalue","padj")], by="gene_ID")

####Extract CgIAPs from Significantly Differentially Expressed Genes####
#Significantly differentially Expressed IAPs
RIF_dfSig_IAP <- grepl("IAP", resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_dfSig_IAP) #62 253 289
RIF_dfSig_other_IAP <- grepl("IAP", resRIFTran_05_Sig_ENTREZGENE_DATA$other_designations, ignore.case = TRUE) 
grep("TRUE",RIF_dfSig_other_IAP) #same
RIF_dfSig_inhibitor <- grepl("inhibitor", resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE",RIF_dfSig_inhibitor) #19  260  466  700  701  968 1048 1818
resRIFTran_05_Sig_ENTREZGENE_DATA[c(19,  260,  466,  700,  701,  968, 1048, 1818),]
  #none are inhibitors of apoptosis
RIF_dfSig_IAP_info <- resRIFTran_05_Sig_ENTREZGENE_DATA[c(62, 253),]

#non Sig IAPs
RIF_non_Sig_IAP <- grepl("IAP", resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_non_Sig_IAP) # 1270  1278  1279  1300  1318  1358  1359  1360  1580  1581  1962  1990
#  2643  2644  2800 3469  3470  3589  3674  3850  3854  3958  3959  3960  3961  4104
#  4105  4154  4155  4156  4520  4541  4542  4746  5026  5027  5028  5078  5079  5256  5257  5815
#  5820  6117  6118  6119  6130  6210  7192  9878  9879 18434 23779 24891 26145 30549 30550 30551
# 30552 31681 31682 31683 34258
RIF_non_Sig_IAP <- resRIFTran_05_non_Sig_ENTREZGENE_DATA[c(  1270,  1278,  1279,  1300,  1318,  1358, 1359,  
                                                             1360,  1580,  1581,  1962 , 1990, 2643,  2644 ,
                                                             2800 ,3469,  3470,  3589,  3674,  3850,  3854, 
                                                             3958,  3959,  3960,  3961,  4104, 4105,  4154, 
                                                             4155,  4156,  4520,  4541,  4542 , 4746,  5026,  
                                                             5027,  5028,  5078,  5079,  5256,  5257,  5815, 
                                                             5820,  6117,  6118,  6119,  6130,  6210,  7192, 
                                                             9878 , 9879, 18434, 23779, 24891, 26145, 30549, 
                                                             30550, 30551, 30552, 31681, 31682, 31683, 34258),]
RIF_non_Sig_inhibitor <- grepl("inhibitor", resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE",RIF_non_Sig_inhibitor) 
229 ,  280,   281 ,  375,   376 ,  666,  1154 , 1248,  1735 , 1822,  1823 , 1959,  1960 , 2101,  2102,  2316,  2446 , 2645,  2899 , 2962,
3228,  3260,  3296,  3361,  4477,  4532,  4533,  4865,  4949,  5362,  5471,  5548,  6103,  6265,  6301,  6444,  6679,  6704,  6728,  6729,
6817,  6818,  6901,  7038,  7069,  7786,  7787,  8258 , 8661,  9026 , 9068,  9104 , 9112,  9152,  9303,  9602,  9844,  9931,  9944,  9984,
10052, 10243, 10610, 10860, 10918, 11132, 11341, 11600, 11949, 12827, 12830, 13045, 14138, 14395, 14457, 14579, 14602, 14980, 15302, 15632,
15702, 15808, 15853, 15854 ,15855, 15883, 16168, 16459, 16610, 16611, 16956, 17505, 17690, 18012, 18013, 18036, 18037, 18038, 18039, 19008,
19764, 19967, 19968, 19969, 19970, 20310, 20476, 20763, 21038, 21160, 22030, 22031, 22032, 22033, 22130, 22349, 22701, 22702, 22703, 22770,
23426, 23562, 24519, 26429, 26593, 27467, 28582, 28583, 29740, 29755, 29812, 29813, 29814, 29815, 29816, 30214, 30523, 30524, 30525, 30526,
30527, 30573, 31266, 32699, 32753, 32823, 32858 ,32931, 33534, 33535, 33592, 33768, 33819, 33876, 34685, 34775, 35249, 35254, 35255, 35256,
36121, 36122, 36123, 36124, 36196, 36282, 36283, 37064, 37065, 37066, 37067, 37068, 37069, 38156, 38188, 38497, 40086, 40165, 40649, 40838,
40839, 41471
resRIFTran_05_non_Sig_ENTREZGENE_DATA[c(229 ,  280,   281 ,  375,   376 ,  666,  1154 , 1248,  1735 , 1822,  1823 , 1959,  1960 , 2101,  2102,  2316,  2446 , 2645,  2899 , 2962,
                                        3228,  3260,  3296,  3361,  4477,  4532,  4533,  4865,  4949,  5362,  5471,  5548,  6103,  6265,  6301,  6444,  6679,  6704,  6728,  6729,
                                        6817,  6818,  6901,  7038,  7069,  7786,  7787,  8258 , 8661,  9026 , 9068,  9104 , 9112,  9152,  9303,  9602,  9844,  9931,  9944,  9984,
                                        10052, 10243, 10610, 10860, 10918, 11132, 11341, 11600, 11949, 12827, 12830, 13045, 14138, 14395, 14457, 14579, 14602, 14980, 15302, 15632,
                                        15702, 15808, 15853, 15854 ,15855, 15883, 16168, 16459, 16610, 16611, 16956, 17505, 17690, 18012, 18013, 18036, 18037, 18038, 18039, 19008,
                                        19764, 19967, 19968, 19969, 19970, 20310, 20476, 20763, 21038, 21160, 22030, 22031, 22032, 22033, 22130, 22349, 22701, 22702, 22703, 22770,
                                        23426, 23562, 24519, 26429, 26593, 27467, 28582, 28583, 29740, 29755, 29812, 29813, 29814, 29815, 29816, 30214, 30523, 30524, 30525, 30526,
                                        30527, 30573, 31266, 32699, 32753, 32823, 32858 ,32931, 33534, 33535, 33592, 33768, 33819, 33876, 34685, 34775, 35249, 35254, 35255, 35256,
                                        36121, 36122, 36123, 36124, 36196, 36282, 36283, 37064, 37065, 37066, 37067, 37068, 37069, 38156, 38188, 38497, 40086, 40165, 40649, 40838,
                                        40839, 41471),]
#putative inhibitor of apoptosis: 1822,`1823,1959, 1960, 2962, 3228, 4477, 
RIF_non_Sig_inhibitor <- resRIFTran_05_non_Sig_ENTREZGENE_DATA[c(1822, 1823,1959, 1960, 2962, 3228, 4477),]

RIF_non_Sig_IAP_full_info <- resRIFTran_05_non_Sig_ENTREZGENE_DATA[c(1822, 1823,1959, 1960, 2962, 3228, 4477,
                                                                     1270,  1278,  1279,  1300,  1318,  1358, 1359,  
                                                                     1360,  1580,  1581,  1962 , 1990, 2643,  2644 ,
                                                                     2800 ,3469,  3470,  3589,  3674,  3850,  3854, 
                                                                     3958,  3959,  3960,  3961,  4104, 4105,  4154, 
                                                                     4155,  4156,  4520,  4541,  4542 , 4746,  5026,  
                                                                     5027,  5028,  5078,  5079,  5256,  5257,  5815, 
                                                                     5820,  6117,  6118,  6119,  6130,  6210,  7192, 
                                                                     9878 , 9879, 18434, 23779, 24891, 26145, 30549, 
                                                                     30550, 30551, 30552, 31681, 31682, 31683, 34258),]



#GIMAP
#Sig GIMAP
RIF_sig_GIMAP <- grepl("GTPase",resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE",RIF_sig_GIMAP)
RIF_sig_GIMAP_info <-  resRIFTran_05_Sig_ENTREZGENE_DATA[508,]
RIF_sig_IMAP <- grepl("IMAP", resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE", RIF_sig_IMAP)

#non-Sig GIMAP
RIF_non_sig_GIMAP <- grepl("GTPase", resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE", RIF_non_sig_GIMAP) 
RIF_non_sig_IMAP <- grepl("IMAP", resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE", RIF_non_sig_IMAP)
RIF_non_sig_IMAP_info <- resRIFTran_05_non_Sig_ENTREZGENE_DATA[c(842,   843 , 5846,  5966 , 5967,  5968 , 5969,  5970,  
                                                                 5971,  5974,  6318,  6531,  6643,  6851,  6868,  7431, 
                                                                 8617,  8657,  8983,  9506,
                                                                 9687,  9688,  9842, 11734 ,18676 ,18677 ,18678, 19559 ,
                                                                 31117, 31119, 31441, 31442, 31443, 31445),]
#all are GIMAP!

###RI probiotic challenge Additional Apoptotic Transcripts####
#Caspases (caspase 2 is important with IAPS)
#Sig caspases (none) Bac
RIF_caspase_Sig <- grepl("caspase",resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE", RIF_caspase_Sig) #842  843 1048 1209 1557
RIF_caspase_Sig_info <- resRIFTran_05_Sig_ENTREZGENE_DATA[c(842,  843, 1209, 1557),]

#non sig caspase
RIF_caspase_non_Sig <- grepl("caspase",resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE", RIF_caspase_non_Sig) #145  5600  5601  5649  5650  5651  6548  6722  7135  8068  8269  8789 10882 14776 16985 16987 17151 17152 17153 17154
18711 18712 19502 20065 20066 20067 22030 22031 22032 22033 22364 23127 23128 24871 24872 24935 24936 25000 25842 27451
27452 32565 34700 34703 38137 38138 38871 38872 39310 39536 39537 40930
RIF_caspase_non_Sig_info <- resRIFTran_05_non_Sig_ENTREZGENE_DATA[c(145,  5600,  5601,  5649,  5650,  5651,  6548,  6722,  7135,  8068,  8269,  8789, 10882, 14776, 16985, 16987, 17151, 17152, 17153, 17154,
                                                                  19502, 22030, 22031, 22032, 22033, 22364, 23127, 23128, 24871, 24872, 24935, 24936, 25000, 25842, 27451,
                                                                    27452, 32565, 34700, 34703, 38137, 38138, 38871, 38872, 39310, 39536, 39537, 40930),]
#BcL 2 
#Sig Bcl2
RIF_Bcl2_Sig <- grepl("Bcl", resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_Bcl2_Sig) #0

#non sig Bcl2
RIF_Bcl2_non_Sig <- grepl("Bcl", resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_Bcl2_non_Sig) #12386 12387 36507
RIF_Bcl2_non_Sig_info <- resRIFTran_05_non_Sig_ENTREZGENE_DATA[c(12386, 12387),]

#BAG
#Sig BAG
RIF_BAG_Sig <- grepl("BAG", resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_BAG_Sig) #0
RIF_athanogene_Sig <- grepl("athanogene", resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_athanogene_Sig ) #0

#non sig BAG
RIF_BAG_non_Sig <- grepl("BAG", resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_BAG_non_Sig) #4455 21607 22282 22283 28207
RIF_athanogene_non_Sig <- grepl("athanogene", resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_athanogene_non_Sig ) #0
RIF_BAG_non_Sig_info <- resRIFTran_05_non_Sig_ENTREZGENE_DATA[c(4455, 21607, 22282, 22283, 28207),]


#IFNLP: , IFN-like protein 
#sig 
RIF_IFN_Sig <- grepl("IFN", resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE",RIF_IFN_Sig ) #0
RIF_interferon_Sig <- grepl("interferon", resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE",  RIF_interferon_Sig ) #0

#non IFNLP
RIF_IFN_non_Sig <- grepl("IFN",resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE",RIF_IFN_non_Sig ) #0
RIF_interferon_non_Sig <- grepl("interferon", resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE",  RIF_interferon_non_Sig ) #4839  6641  6801  6806  6852  7097  7189  7441  9033  9035  9193  9194  9195  9196  9197  9198  9501  9591  9711  9714
# 10390 11025 12579 12580 12581 12583 12584 14668 16762 16763 24539 24540 24541 25813 29989 30387 32001 34442 34443 34444
#  35065 36376 36389 36390 39676 40068
RIF_interferon_non_sig_info <- resRIFTran_05_non_Sig_ENTREZGENE_DATA[c(4839,  6641,  6801,  6806,  6852,  7097,  7189,  7441,  9033,  9035,  9193,  9194,  9195,  9196,  9197,  9198,  9501,  9591,  9711,  9714,
                                                                       10390, 11025, 12579, 12580, 12581, 12583, 12584, 14668, 16762, 16763, 24539, 24540, 24541, 25813, 29989, 30387, 32001, 34442 ,34443, 34444,
                                                                       35065, 36376, 36389, 36390, 39676, 40068),]
#NONE are interferon per se 

#Cytochrome c
#Sig cytochrome 
RIF_cytochrome_Sig <- grepl("cytochrome",resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE",RIF_cytochrome_Sig) #241  459  892  894  895 1059 1917
RIF_cytochrome_Sig_info <-resRIFTran_05_Sig_ENTREZGENE_DATA[c(241,  459,  892,  894,  895, 1059, 1917),]
#none are cytochrome c

#non sig cytochrome c
RIF_cytochrome_non_Sig <- grepl("cytochrome c",resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE",RIF_cytochrome_non_Sig) # 4816 12470 16434 17036 18324 19897 19898 20336 20550 20633 20932 25078 26684 27977 32025 32026 33616 33633 34200 35299
# 38664 40196, 
# real cytochrome proteins: 12470, 33633 , 40196
RIF_cytochrome_non_Sig_info <-resRIFTran_05_non_Sig_ENTREZGENE_DATA[c(12470, 33633 , 40196),]

#FasL, FAIM, FADD
#sig
RIF_Fas_Sig <- grepl("Fas",resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_Fas_Sig ) #0

#non sig Fas
RIF_Fas_non_Sig <- grepl("Fas",resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_Fas_non_Sig ) #12143 12144 13193 15152 16826 16913 19951 21151 31224 37230 38188
RIF_FAIM_non_Sig_info <- resRIFTran_05_non_Sig_ENTREZGENE_DATA[38188,] #only 1 is Fas apoptotic inhibitory molecule
RIF_FADD_non_sig_info <- resRIFTran_05_non_Sig_ENTREZGENE_DATA[16913,]


#smac/DIABLO
#sig
RIF_smac_Sig <- grepl("smac",resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_sAC_Sig ) #0
RIF_DIABLO_Sig <- grepl("DIABLO",resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_DIABLO_Sig ) #0
RIF_second_sig <- grepl("second", resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE", RIF_second_sig)
RIF_activator_Sig <- grepl("activator of caspases", resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE", RIF_activator_Sig) #0
RIF_mitochondrial_Sig <- grepl("mitochondrial-derived", resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE", RIF_mitochondrial_Sig) #0
RIF_IAPbinding_Sig <- grepl("IAP-binding",  resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE", RIF_IAPbinding_Sig) #0

#non Sig smac/DIABLO
RIF_smac_non_Sig <- grepl("smac",resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_smac_non_Sig ) #0
RIF_DIABLO_non_Sig <- grepl("DIABLO",resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_DIABLO_non_Sig ) #31525 32575 (kelch like protein diablo..NO)
RIF_second_non_sig <- grepl("second", resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE", RIF_second_non_sig)
RIF_activator_non_Sig <- grepl("activator of caspases", resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE", RIF_activator_non_Sig) #0
RIF_mitochondrial_non_Sig <- grepl("mitochondrial-derived", resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE", RIF_mitochondrial_non_Sig) #0
RIF_IAPbinding_non_Sig <- grepl("IAP-binding", resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE", RIF_IAPbinding_non_Sig) #0

#NONE

#Extrinsic pathway molecules
#PCD and Death domain containing protein (DD), Death effector domain (DED)
#sig 
RIF_PCD_Sig <- grepl("death",resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_PCD_Sig ) #1413
RIF_DD_Sig <- resRIFTran_05_Sig_ENTREZGENE_DATA[1413,]

# non Sig
RIF_PCD_non_Sig <- grepl("death",resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_PCD_non_Sig ) #1125,  1126,  1735,  3375,  4258,  4259,  4331,  5531,  5532,  5553,  5789,  8193,  8194,  8442,  8919,  8920,  9543,  9544, 10893, 10942,
#11451, 11452, 11873, 12292, 12293, 13006, 13007, 13008, 13342, 15295, 15698, 16778, 16779, 16780, 16781, 16782, 16913, 17236, 18893, 18943,
#21131, 21132, 22050, 22540, 22853, 26895, 27193, 27481, 28479, 28528, 29120, 29121, 29122, 29123, 31941, 33413, 35522, 35593, 35669, 36585,
#36586, 36587, 36588, 37184, 37185, 37186, 37941, 37947, 37948, 38035, 39597, 41433, 41434
RIF_death_non_sig_info <- resRIFTran_05_non_Sig_ENTREZGENE_DATA[c(8193,  8194,  8442,  9543,  9544,
  11451, 11452, 13006, 13007, 13008, 18943,
   22540, 26895, 27193, 27481, 29120, 29121, 29122, 29123, 33413, 35522, 35593, 36585,
  36586, 36587, 36588, 37184, 37185, 37186, 37941, 37947, 37948, 38035, 39597, 41433, 41434),]
RIF_PCD_non_sig_info <- resRIFTran_05_non_Sig_ENTREZGENE_DATA[c(17236,21131, 21132, 22050, 31941),]
RIF_DED_non_sig_info <- resRIFTran_05_non_Sig_ENTREZGENE_DATA[22853,]

#TRAF
#sig 
RIF_TRAF_Sig <- grepl("traf",resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_TRAF_Sig ) #1046, 1912 NONE
RIF_tumornecrosis_Sig <- grepl("necrosis",resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_tumornecrosis_Sig  ) #0

#non Sig TRAF
RIF_TRAF_non_Sig <- grepl("traf",resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_TRAF_non_Sig ) #1594,  4527  4528  5235  5338  6599  9222 11309 15949 16631 16644 17700 18143 18740 19696 19990 20005 20937 21009 21446
21960 22573 23418 23419 23420 24299 25129 26154 26910 26911 26912 26913 26914 26915 27654 27655 27656 27657 27658 27786
[41] 29600 29601 30674 30675 30676 33527 33528 33529 33530 33531 37432 37433 37860 38633 38710 38986 39025 39026 39125 40031
[61] 41400
#NONE

RIF_tumornecrosis_non_Sig <- grepl("necrosis",resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_tumornecrosis_non_Sig  ) #[1]    96  4557  4558  4559  4560  4698  4699  5260  5261  5262  5263  5313  5343  5872  5948  6556  6559  7932  9974 10103
# 10703 10704 10705 12203 15580 16060 16117 17682 20490 20564 20565 20768 20769 20770 20771 20772 24922 24925 25843 27136
# 27137 27494 30721 34502 35554 36170 36171 36172 36173 38034 38036 38396 38397
#NONE

#TNF
#sig none
RIF_TNF_Sig <-  grepl("tumor necrosis factor", resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_TNF_Sig) #0

#non Sig 
RIF_TNF_non_Sig <- grepl("tumor necrosis factor",resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_TNF_non_Sig  ) #96  4557  4558  4559  4560  4698  4699 5313  5343  5872  5948  6556  6559  7932  9974 10103
[21] 10703 10704 10705 12203 15580 16060 16117 17682 20490 20768 20769 20770 20771 20772 24922 24925 25843 27136
[41] 27137 27494 34502 35554 36173 38034 38036 38396 38397

RIF_TNF_non_Sig_info <- resRIFTran_05_non_Sig_ENTREZGENE_DATA[c(96,  4557,  4558,  4559,  4560,  4698,  4699, 5313,  5343,  5872,  5948,  6556,  6559,  7932,  9974, 10103,
                                                                10703, 10704, 10705, 12203, 15580, 16060, 16117, 17682, 20490, 20768, 20769, 20770, 20771, 20772, 24922, 24925, 25843, 27136,
                                                                27137, 27494 ,34502, 35554, 36173 ,38034, 38036, 38396, 38397),]


#TNFR
#sig - necrosis found nothing
RIF_TNFR_Sig <- grepl("tumor necrosis factor receptor", resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE",  RIF_TNFR_Sig) #0

#non TNFR
RIF_TNFR_non_Sig <- grepl("tumor necrosis factor receptor",resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_TNFR_non_Sig) # 25843 27136 27137 38396 38397 
RIF_TNFR_non_Sig_info <-  resRIFTran_05_non_Sig_ENTREZGENE_DATA[c(25843, 27136, 27137, 38396, 38397),]


#TRAIL (should have come up with TNF) TNF-related apoptosis inducing ligand
#sig 
RIF_TRAIL_Sig <- grepl("trail",resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_TRAIL_Sig ) #none 

#non TRAIL
#none

#siglec (sialic acid binding immunoglobulin-type lectin)
#sig 
RIF_siglec_Sig <- grepl("siglec",resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE",RIF_siglec_Sig ) #0
RIF_sialic_Sig <- grepl("sialic",resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE",RIF_sialic_Sig ) #0
RIF_immunoglobulin_lectin_Sig <- grepl("immunoglobulin", resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_immunoglobulin_lectin_Sig ) #0

#non siglec
RIF_siglec_non_Sig <- grepl("siglec",resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE",RIF_siglec_non_Sig ) #0
RIF_sialic_non_Sig <- grepl("sialic",resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE",RIF_sialic_non_Sig ) #0
RIF_immunoglobulin_lectin_non_Sig <- grepl("immunoglobulin", resRIFTran_05_non_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_immunoglobulin_lectin_non_Sig ) #1405  3161  3261  6843 15413 15414 16468 16469 24157 28662 35966 35967 35968 35969 37380 41348
RIF_immunoglobulin_lectin_non_Sig_info <- resRIFTran_05_non_Sig_ENTREZGENE_DATA[c(1405,  3161,  3261,  6843, 15413, 15414,
                                                                                  16468, 16469, 24157, 28662, 35966, 35967, 35968, 35969 ,37380,
                                                                                  41348),] 
#none

#cgBTG 1
#sig 
RIF_BTG_Sig <- grepl("BTG",resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_BTG_Sig ) #0
RIF_Bcell_Sig <- grepl("translocation",resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE) 
grep("TRUE", RIF_Bcell_Sig ) #0

#non cgBTG-1

#p53 (can induce apoptosis, sokolova 2009)
#sig
RIF_p53_Sig <- grepl("p53", resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE", RIF_p53_Sig ) #0

#non p53

#TLR
RIF_TLR_Sig <- grepl("Toll",resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE", RIF_TLR_Sig) #499,596
RIF_TLR_Sig_info <- resRIFTran_05_Sig_ENTREZGENE_DATA[c(499,596),]
                                                      
#non TLR

#cAMP - can induce apoptosis in C. gigas hemocytes (sokolova 2009)
#Sig
RIF_cAMP_Sig <- grepl("cAMP",resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE", RIF_cAMP_Sig) #743  744  745  850 1026 , not cAMP

RIF_cyclic_Sig <- grepl("cyclic",resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE", RIF_cyclic_Sig) 
RIF_adenosine_Sig <- grepl("adenosine",resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE", RIF_adenosine_Sig) #561 562, none

#non sig cAMP

#Nitric oxide may be important in apoptosis regulation! (sokolova 2009)
#sig
RIF_nitric_Sig <- grepl("nitric",resRIFTran_05_Sig_ENTREZGENE_DATA$description, ignore.case = TRUE)
grep("TRUE", RIF_nitric_Sig) #0

#non Sig nitric oxide

#COMPILE AND GRAPH THE APOPTOSIS GENES THAT WERE SIGNIFICANT
#Add type column
RIF_dfSig_IAP_info['Type'] <- "IAP"
RIF_caspase_Sig_info['Type'] <- "Caspase"
RIF_TLR_Sig_info["Type"] <- "TLR"

RIF_Sig_APOPTOSIS_GENES_DATA <- rbind(RIF_dfSig_IAP_info, RIF_caspase_Sig_info, RIF_TLR_Sig_info)


#subet relevant columns
RIF_Sig_APOPTOSIS_GENES_DATA_subset <- RIF_Sig_APOPTOSIS_GENES_DATA[,c("gene_ID","description","pvalue","padj","log2FoldChange","Type")]

RIF_Sig_APOPTOSIS_GENE_DATA_PLOT <- ggplot(RIF_Sig_APOPTOSIS_GENES_DATA_subset, aes(x=gene_ID, y=log2FoldChange, fill=Type)) + 
  theme_set(theme_classic(base_size=15))+
  geom_bar(stat="identity", position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1, size=15, face="bold")) + 
  scale_y_continuous(name ="log2 Fold Change", breaks = scales::pretty_breaks(n = 10)) + 
  ggtitle("Log2 Fold Change of C. virginica Significant\nApoptosis Transcripts under Bacterial Challenge") + 
  theme(plot.title = element_text(size=15, face="bold")) + theme(plot.title = element_text(hjust = 0.5, size=15)) 


RIF_Sig_APOPTOSIS_GENE_DATA_PLOT2 <- RIF_Sig_APOPTOSIS_GENE_DATA_PLOT + geom_hline(yintercept=0) +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  scale_x_discrete(name="Gene Name", limits=c("LOC111100471", "LOC111104430", "LOC111118288", "LOC111118290",
                                              "LOC111124891", "LOC111130912", "LOC111110045", "LOC111112158"), 
                   labels=c("LOC111100471"="IAP protein 3-like", "LOC111104430"="IAP protein 3-like", "LOC111118288"="caspase-3-like", "LOC111118290"="caspase-3-like",
                            "LOC111124891"="caspase-7-like", "LOC111130912"="caspase-3-like", "LOC111110045"="toll-like receptor 6", "LOC111112158"="toll-like receptor 6"))



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
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4743157/
