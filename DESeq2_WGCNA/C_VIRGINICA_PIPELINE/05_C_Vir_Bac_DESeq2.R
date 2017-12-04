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
install.packages("tm")
library(tm)
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
C_vir_TranscriptCountData <- C_vir_TranscriptCountData[c(1,3,4,11,5,6,7,12,19,20,21,
                                                         8,9,10,2,13,16,17,18,14,15)]
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
resRODTran_05_nonSig_df <- data.frame(resRODTran_05_nonSig)
resRODTran_05_df_nonSig_FULL <- merge(resRODTran_05_nonSig_df, C_vir_stringtie_transcripts_rna_LOC_separated_clean_transcript_id[,c("ID", "gene_name")], by="ID")
#put LOC info in new column
resRODTran_05_df_nonSig_FULL <- resRODTran_05_df_nonSig_FULL %>% separate(gene_name, c("gene_name", "gene_ID"), ";")

####Lookup LOC values using Batch Entrez ####
#Sig
write.csv(resRODTran_05_Sig_df_FULL$gene_ID, file="resRODTran_05_Sig_df_FULL_gene_ID.csv")
#put column with LOC info into text file
#perform batch Entrez lookup to the gene database to retrieve IDs

#non_Sig
write.csv(resRODTran_05_df_nonSig_FULL$gene_ID, file="resRODTran_05_df_nonSig_FULL_gene_ID.csv")
#put column wit LOC into text file
#perform batch Entrez lookup with the gene database to retrieve IDs
#batch Entrez gets rid of duplicates

#####Merge the LOC values with gene name####
#SIG
#convert file to csv using excel 
resRODTran_05_Sig_ENTREZGENE <- read.csv(file="resRODTran_05_df_Sig_ENTREZ_RESULTS_FULL.csv", header = TRUE)
#change symbol to gene_ID
colnames(resRODTran_05_Sig_ENTREZGENE)[6] <- "gene_ID"
#merge file based on the LOC column
resRODTran_05_Sig_ENTREZGENE_DATA <-  merge(resRODTran_05_Sig_ENTREZGENE, resRODTran_05_Sig_df_FULL[,c("gene_ID", "log2FoldChange","pvalue","padj")], by="gene_ID")
head(resRODTran_05_Sig_ENTREZGENE_DATA)

#NON SIG
resRODTran_05_NON_Sig_ENTREZGENE <- read.csv(file="resRODTran_05_df_non_Sig_ENTREZ_RESULTS_FULL.csv", header = TRUE)
#change symbol to gene_ID
colnames(resRODTran_05_NON_Sig_ENTREZGENE)[6] <- "gene_ID"
#merge file based on the LOC column
resRODTran_05_NON_Sig_ENTREZGENE_DATA <-  merge(resRODTran_05_NON_Sig_ENTREZGENE, resRODTran_05_df_nonSig_FULL[,c("gene_ID", "log2FoldChange","pvalue","padj")], by="gene_ID")
head(resRODTran_05_Sig_ENTREZGENE_DATA)

#add columns for significance and merge the two files so only one search needs to be done
resRODTran_05_Sig_ENTREZGENE_DATA$Significance <- sig
resRODTran_05_NON_Sig_ENTREZGENE_DATA$Significance <- nonsig

resRODTran_05_FULL_ENTREZGENE_DATA <- rbind(resRODTran_05_Sig_ENTREZGENE_DATA,resRODTran_05_NON_Sig_ENTREZGENE_DATA)

####Add full list of apoptosis related genes and aliases to grep for
#https://stackoverflow.com/questions/17130129/find-matches-of-a-vector-of-strings-in-another-vector-of-strings
#https://stackoverflow.com/questions/37684179/how-to-search-for-multiple-words-in-the-same-string-using-r

# full data frame
#make the description column a character
resRODTran_05_FULL_ENTREZGENE_DATA$description <- as.character(resRODTran_05_FULL_ENTREZGENE_DATA$description)

# Sample vector of keywords or phrases
  #Search terms come from the strings that found hits in the C_virginica genome annotation view in NCBI:
  #all terms with no hits in the annotation on NCBI were INCLUDED IN THIS LIST..just in case, and all aliases searched for

apoptosis_keywords <- as.character(c("IP3R","inositol","RhoA","RhoB","RhoC","rho","RhoE","ras","Rac1","Cdc42",
 "PKC","IkB","kappa","NFKB","NF-kappa-B","GPCR","Apo3","TWEAK","DR3LG","TNFSF12","Apo","TNFSF10","TNFSF","TRAIL",
"Apo2L","c-FLIP","FLICE", "Casper","FLAME","FLAME-1","CASH","CLARP","MRIT","FADD-like","IL-1β","I-FLICE",
 "RIP","RIPK1","RIPK","TNFR-interacting","TNFR","AMP","adenosine","cyclic","kinase","PKC","adenylate","CREB",
"TNFRSF12","Apo3","TRAMP","DDR3","LARD","WSL","TNFRSF","DR4","TNFRSF10A","TRAILR1",
 "TRAILR","APO2","TNFRS10B","TRAIL-R2","TRICK2","ZTNFR9","KILLER","Fas","FASL","TNFSF6","FasR",
 "TNF-alpha","alpha","Fas","FAIM","TNF","TRADD","TNFRSF1A","Fas","casp","casp2","casp8","caspase",
"procasp8","caspase","caspase","procaspase","procasp","siglec","immunoglobulin","sialic",
"Cgsiglec-1","IFNLP","IFN","interferon","interferon-like","TNF","TNFRSF1A","TRAF",
"mitotic","ICAD","DFF-45","Dnase","CAD","DFF-40","CPAN","caspase","casp6","caspase","caspase",
"casp3/7","caspase","casp7","caspase","casp3","caspase","casp14","caspase",
"caspase","casp9","caspase","APAF","APIP","apoptotic","APAF-1","protease","SAC",
 "soluble","AC","BI-1", "bcl-2", "bcl","PDRP","p53","cyto","apoptosis","APAF-1","APAF",
  "procaspase","p55-BLK","Blk","MGC10442","myc","oncogene","aven","NY-CO-13","trp53","p53",
"smac","DIABLO","IAP","IAP-binding","xL","bcl2","bcl-2", "xS","bcl-w","myeloid","mcl-1",
"mcl","leukemia","A1","Bak","BCL2L7","antagonist", "Bcl-XL","BAD","BAG","BH3","Bid",
"tBid","BIM","Bik", "JFY1","puma","ceramide", "MEK","MAPK","p38","phosphatidyl","cyclin",
 "JUN","p39","c-jun","AP-1","proto-oncogene", "AP1","JUN","CD151","BTG","caspase",
"IAP","apoptosis","CAAP","caspase","NR13","bcl-2","PRKC","WT1","PPP1R13B",
 "p53","IFI44","interferon","IFI6","TLR","MyD88","myeloid","CCAR", "GIMAP", "ovarian", "bok"))

apoptosis_multiple_keywords <- as.character(c("protein kinase c","NF-kappa-B inhibitor", "nuclear factor","G protein","Apo3 ligand",
"Apo2 ligand","Receptor interacting","threonine kinase", "serine threonine kinase","Protein kinase A","adenylate cyclase","cyclic responsive",
"cAMP responsive", "death receptor","death receptor 4", "death receptor", "death receptor 5","fatty acid synthase","tumor necrosis factor",
"tumor necrosis factor","tumor necrosis factor","apoptotic inhibitory","death domain","receptor associated","Fas associated death domain","mitotic apparatus",
"nuclear mitotic","SP-H antigen","adenylyl cyclase","bax inhibitor","DNA damage","cytochrome c","apoptosis inducing","peptidase activity","protease activity",
"apoptotic protease","apoptotic pedtidase","B lymphoid", "tyrosine kinase","tumor antigen","tumor protein","low PI","mitochondrial activator",
"activator caspase","cell differentiation","cell death","antagonist killer","bax","bcl-2","ovarian","bok","ovarian killer","associated death promoter",
"binding component","tBid","interacting killer","modulator of apoptosis", "MEK Kinase", "endonuclease G", "akt-kinase","protein kinase B",
"kinase B","N-terminal kinase","cyclin dependent cyclase 5 activator","cyclase 5","proto-oncogene","transcription factor","B cell translocation",
"programmed cell death","WT1 regulator", "apoptosis regulator", "apoptosis-stimulating","cell division","apoptosis regulator","death agonist"))


#remove duplicates
apoptosis_keywords_unique <- unique(unlist(strsplit(apoptosis_keywords, " ")))
apoptosis_multiple_keywords_no_space_unique <- unique(unlist(strsplit(apoptosis_multiple_keywords_no_space, " ")))

#could paste all "keywords" together and separate them with the pipe character (|) which will work like an "or" statement
resRODTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_SINGLE_WORDS <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl(paste(apoptosis_keywords_unique, collapse="|"), ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
#remove all lines containing the word "uncharacterized"
resRODTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_SINGLE_WORDS <-  resRODTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_SINGLE_WORDS[!grepl("uncharacterized", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_SINGLE_WORDS$description),]
write.csv(resRODTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_SINGLE_WORDS, file="resRODTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_SINGLE_WORDS.csv")

#searching for each one individually:
IP3R <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("inositol", resRODTran_05_FULL_ENTREZGENE_DATA$description) & 
                                             grepl("phosphate", resRODTran_05_FULL_ENTREZGENE_DATA$description) &  grepl("receptor", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
IP3Ra <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("IP3R", resRODTran_05_FULL_ENTREZGENE_DATA$description),] #8 
Rho <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("ras", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) & 
                                             grepl("GTP", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) &  grepl("rho", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
Rac1 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("ras", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) & 
                                           grepl("botulinum", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) &  grepl("C3", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
cdc42 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("cdc42", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                              !grepl("activated", ignore.case=TRUE,resRODTran_05_FULL_ENTREZGENE_DATA$description),]

PKCiota <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("protein", resRODTran_05_FULL_ENTREZGENE_DATA$description) & 
                                            grepl("kinase", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl("iota", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
PKCdelta <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("protein", resRODTran_05_FULL_ENTREZGENE_DATA$description) & 
                                                 grepl("kinase", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl("delta", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                                 !grepl("calmodulin", resRODTran_05_FULL_ENTREZGENE_DATA$description) & !grepl("ribosomal", resRODTran_05_FULL_ENTREZGENE_DATA$description) ,]

IkB <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("NF-kappa-B", ignore.case= TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) & 
                                            grepl( "inhibitor", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) & !grepl("ras-like", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
NFkB <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("NF-kappa-B", ignore.case= TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) & 
                                             grepl( "nuclear", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) & !grepl("inhibitor", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]

GPCRsevenpass <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("G-type", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl("seven-pass", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
GPCR <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("G-protein" , resRODTran_05_FULL_ENTREZGENE_DATA$description) & 
                                             grepl("coupled", resRODTran_05_FULL_ENTREZGENE_DATA$description) &  grepl("receptor", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

# 0 hits
Apo3 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("Apo3",ignore.case= TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
# 0 hits
Apo2 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("Apo2",ignore.case= TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
# 0 hits
ApoL <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("Apo" , ignore.case= TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl("ligand", ignore.case= TRUE,resRODTran_05_FULL_ENTREZGENE_DATA$description),]
# 0 hits
TWEAK <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("TWEAK",ignore.case= TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
# 0 hits
TRAIL <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("TRAIL",ignore.case= TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
# 0 hits
TRAILb <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("TNF-related", ignore.case= TRUE,  resRODTran_05_FULL_ENTREZGENE_DATA$description) & 
                                             grepl("apoptosis",ignore.case= TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) &  grepl("ligand",ignore.case= TRUE,  resRODTran_05_FULL_ENTREZGENE_DATA$description),]
#0 hits
TRAILc <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("TNFSF12",ignore.case= TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
# 0 hits
cFLIP <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("c-FLIP",ignore.case= TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
# 0 hits
cFLIP2 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("FLICE",ignore.case= TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
cFLIP3 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("CASPER",ignore.case= TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
cFLIP4 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("FLAME",ignore.case= TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
cFLIP5 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("CASH",ignore.case= TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
cFLIP6 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("CLARP",ignore.case= TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
cFLIP7 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("MRIT",ignore.case= TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
cFLIP8 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("FADD-like",ignore.case= TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) &
   grepl("inhibitory",ignore.case= TRUE,  resRODTran_05_FULL_ENTREZGENE_DATA$description),]

RIP <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("receptor", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl("interacting", resRODTran_05_FULL_ENTREZGENE_DATA$description)
                                          & !grepl("thyroid", resRODTran_05_FULL_ENTREZGENE_DATA$description)
                                          & !grepl("glutamate", resRODTran_05_FULL_ENTREZGENE_DATA$description)
                                          & !grepl("cannabinoid", resRODTran_05_FULL_ENTREZGENE_DATA$description)
                                          & !grepl("AH", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
#non are cAMP specifically
cAMP <-  resRODTran_05_FULL_ENTREZGENE_DATA[grepl("cAMP",ignore.case= TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]

#no hits
cAMP2 <-  resRODTran_05_FULL_ENTREZGENE_DATA[grepl("cyclic",ignore.case= TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) 
                                             & grepl("adenosine", ignore.case = TRUE,resRODTran_05_FULL_ENTREZGENE_DATA$description),]
#none are PKA
PKA <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("protein", resRODTran_05_FULL_ENTREZGENE_DATA$description) & 
                                                grepl("kinase", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl("A", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

AC <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("adenylate", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl("cyclase", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
CREB <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("cAMP-responsive", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]

#0 hits
DR <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("death",ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl("receptor", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
DR1 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("TRAMP", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
DR2 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("WSL",ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
DR3 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("LARD",ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
DR4 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("DDR3",ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]

Fas <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("fatty", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl( "acid", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl( "synthase", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                            !grepl("acyl",resRODTran_05_FULL_ENTREZGENE_DATA$description),]


FasL <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("FasL", ignore.case = TRUE ,resRODTran_05_FULL_ENTREZGENE_DATA$description),]
FASL2 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("TNF", ignore.case = TRUE ,resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                              grepl("superfamily", ignore.case = TRUE ,resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                              grepl("6", ignore.case = TRUE ,resRODTran_05_FULL_ENTREZGENE_DATA$description),]
#none are what I want
Fas2 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("Fas", ignore.case = TRUE ,resRODTran_05_FULL_ENTREZGENE_DATA$description),]

#no DR5
DR5 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("TNFRS10B", ignore.case = TRUE ,resRODTran_05_FULL_ENTREZGENE_DATA$description),]
DR5a <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("TRAIL-R2", ignore.case = TRUE ,resRODTran_05_FULL_ENTREZGENE_DATA$description),]
DR5b <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("TRICK2", ignore.case = TRUE ,resRODTran_05_FULL_ENTREZGENE_DATA$description),]
DR5c <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("KILLER", ignore.case = TRUE ,resRODTran_05_FULL_ENTREZGENE_DATA$description),]
DR5d <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("ZTNFR9", ignore.case = TRUE ,resRODTran_05_FULL_ENTREZGENE_DATA$description),]


TNF <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("tumor", resRODTran_05_FULL_ENTREZGENE_DATA$description) & 
                                            grepl( "necrosis", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl( "factor", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                            !grepl("C1q", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                            !grepl("alpha-induced",resRODTran_05_FULL_ENTREZGENE_DATA$description),]

FAIM <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("fas", resRODTran_05_FULL_ENTREZGENE_DATA$description) & 
                                             grepl("apoptotic", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("inhibitory", resRODTran_05_FULL_ENTREZGENE_DATA$description) ,]

#non are TRAD
TRADD <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("death", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl(  "domain", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl(  "receptor", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

FADD <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("FAS-associated", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("death", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
caspase <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("caspase", ignore.case = TRUE ,resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                                !grepl("activity",resRODTran_05_FULL_ENTREZGENE_DATA$description),]
procaspase <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("procaspase", ignore.case = TRUE ,resRODTran_05_FULL_ENTREZGENE_DATA$description),]
                                                
#no sialic acid
siglec <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("siglec", ignore.case = TRUE ,resRODTran_05_FULL_ENTREZGENE_DATA$description),]
sialic <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("sialic", ignore.case = TRUE ,resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                               grepl("immunoglobulin", ignore.case = TRUE ,resRODTran_05_FULL_ENTREZGENE_DATA$description),]
#none are IFNLP
IFNLP <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("IFNLP", ignore.case = TRUE ,resRODTran_05_FULL_ENTREZGENE_DATA$description),]
IFNLPa <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("interferon", ignore.case = TRUE ,resRODTran_05_FULL_ENTREZGENE_DATA$description),]

TRAF <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("TNF", ignore.case = TRUE ,resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("receptor-associated", ignore.case = TRUE ,resRODTran_05_FULL_ENTREZGENE_DATA$description) ,]

#no NuMA
NuMa <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("mitotic", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "apparatus", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
NuMAa <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("nuclear", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "mitotic", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
NuMAb<- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("SP-H", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

ICAD <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("DNA", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl("fragmentation", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("alpha", resRODTran_05_FULL_ENTREZGENE_DATA$description) ,]
CAD <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("DNA", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl("fragmentation", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl("beta", resRODTran_05_FULL_ENTREZGENE_DATA$description) ,]
CADb <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("caspase", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl("activated", resRODTran_05_FULL_ENTREZGENE_DATA$description) ,]

#no APIP
APIP <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("apoptotic", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "peptidase", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
APIPb <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("apoptotic", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "protease", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
APIPc <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("APAF", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) ,]
APIPc <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("APIP", ignore.case=TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) ,]

#no sAC
sAC <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("adenylyl", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "cyclase", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

#no BI-1
BI <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("bax", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "inhibitor", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

PDRP <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("p53", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "damage", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

cytoc <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("cytochrome", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "c-like", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

AIF <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("apoptosis", resRODTran_05_FULL_ENTREZGENE_DATA$description) & 
                                            grepl( "inducing", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl( "factor", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
#no Blk
Blk <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("blk", ignore.case=TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
Blkb <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("B", ignore.case=TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "lymphoid", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
Blkc <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("p55", ignore.case=TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]

#no myc
myc <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("myc", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
myc <-  resRODTran_05_FULL_ENTREZGENE_DATA[grepl("oncogene", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

aven <-   resRODTran_05_FULL_ENTREZGENE_DATA[grepl("aven", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                               !grepl("scavenger", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
#no p53..
p53 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("tumor", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "antigen", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
p53 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("p53", resRODTran_05_FULL_ENTREZGENE_DATA$description) ,]
p53 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("tumor", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "protein", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

#no smac
smac <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("low", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "PI", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
smac <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("mitochondrial", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "activator", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
smac1 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("IAP", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "direct", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

bcl <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("bcl",ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]

#no mcl
mcl <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("myeloid",ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl("leukemia", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
mcl <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("cell", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "differentiation", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

#no A1
A1 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("bcl-2",ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]

#no Bak
Bak <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("antagonist", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( " killer", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
Bak2 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("bcl2l7",ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]

BAX <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("bax-like",ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]

#no Bok
Bok <-  resRODTran_05_FULL_ENTREZGENE_DATA[grepl("ovarian", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "killer", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
Bok <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("ovarian", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

#no BAD
BAD <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("death", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "promoter", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
BAD2 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("binding", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "component", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
BAD3 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("BAD", ignore.case=TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]

BAG <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("BAG",ignore.case=TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]

#no BH3
BH3 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("BH3",ignore.case=TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
Bid <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("Bid",ignore.case=TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]

#no BIM
BIM <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("BIM",ignore.case=TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]

# no Bik
Bik <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("interacting", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "killer", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

#no Puma
Puma <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("modulator", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "apoptosis", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
Puma <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("JFY", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

ceramide <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("ceramide", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                                 !grepl("galactosylceramide",resRODTran_05_FULL_ENTREZGENE_DATA$description)&
                                                 !grepl("kinase",resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                                 grepl("synthase",resRODTran_05_FULL_ENTREZGENE_DATA$description ) &
                                                 !grepl("phospho", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

MEKK_MAPK <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("mitogen-activated", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl( "kinase", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                             !grepl("kinase-binding", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

#no p38 MAPK
p38_MAPK <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("mitogen-activated", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                                 grepl( "kinase", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                                 grepl("p38", ignore.case=TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]

endoG <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("endonuclease", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "G", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

PI3 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("phosphatidylinositol", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "3-kinase", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

#only AKT interacting protein
AKTB <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("protein", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("kinase", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("B",resRODTran_05_FULL_ENTREZGENE_DATA$description),]
#no protein kinase B!
akt <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("akt", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]

JNK <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("stress-activated", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "JNK", resRODTran_05_FULL_ENTREZGENE_DATA$description),]


p35 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("cyclin", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "dependent", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                           grepl( "5", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl("activator",resRODTran_05_FULL_ENTREZGENE_DATA$description),]
#none
cJUN <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("proto", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "oncogene", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("jun", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
#none
cJUN <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("c-jun", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
#1 
cJUNa <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("AP-1", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                              grepl("transcription",resRODTran_05_FULL_ENTREZGENE_DATA$description),]
CD151 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("CD151", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]

#none
BTG1 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("B", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl("translocation", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

#1 
BTG2 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("BTG", resRODTran_05_FULL_ENTREZGENE_DATA$description) ,]

IAP <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("inhibitor", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl("apoptosis", ignore.case=TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                            !grepl("death-associated", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                            !grepl("caspase", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                            !grepl("5-like", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
IAP2 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("baculoviral", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("IAP", ignore.case = TRUE,resRODTran_05_FULL_ENTREZGENE_DATA$description),]
API <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("inhibitor", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl("apoptosis", ignore.case=TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                            !grepl("death-associated", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                            !grepl("caspase", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl("5-like", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

DIAP <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("inhibitor", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("apoptosis", ignore.case=TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("death-associated", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description),]
CAAP <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("inhibitor", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                                    grepl("apoptosis", ignore.case=TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                                    !grepl("death-associated", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                                    grepl("caspase", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

PCD <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("programmed", resRODTran_05_FULL_ENTREZGENE_DATA$description) & 
                                            grepl("death", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl("cell", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                            !grepl("interacting", resRODTran_05_FULL_ENTREZGENE_DATA$description),]
NR13 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("NR13", resRODTran_05_FULL_ENTREZGENE_DATA$description),]


PAWR <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("WT1", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl("regulator", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

IFI44 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("interferon", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl("44", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

#no IFI6
IFI6 <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("interferon", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl("6", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

TLR <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("toll", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl("receptor", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

MyD88 <-resRODTran_05_FULL_ENTREZGENE_DATA[grepl("myeloid", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl("differentiation", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

  
CCAR <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("cell", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl("division", resRODTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl("apoptosis",resRODTran_05_FULL_ENTREZGENE_DATA$description),]

ASPP <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("apoptosis-stimulating", resRODTran_05_FULL_ENTREZGENE_DATA$description) & grepl("p53", resRODTran_05_FULL_ENTREZGENE_DATA$description),]

#none
GIMAP <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("GIMAP", resRODTran_05_FULL_ENTREZGENE_DATA$description) ,]
GIMAP <- resRODTran_05_FULL_ENTREZGENE_DATA[grepl("GTPase", ignore.case = TRUE, resRODTran_05_FULL_ENTREZGENE_DATA$description) &
  grepl("IMAP",resRODTran_05_FULL_ENTREZGENE_DATA$description),]


#Combine all lists from above into one new big table
#only include those that have hits and are confirmed to be correct
resRODTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED <- rbind(IP3R, Rho, Rac1, cdc42, PKCiota, PKCdelta, IkB,
                                                                         NFkB, GPCRsevenpass, GPCR, RIP,AC,CREB,
                                                                         Fas, TNF, FAIM, FADD, caspase, TRAF,
                                                                         ICAD, CAD, PDRP, cytoc, AIF, bcl, BAX, BAG,
                                                                         ceramide, MEKK_MAPK,endoG, PI3, JNK, p35,
                                                                         cJUNa, CD151, BTG2, IAP, IAP2, API, DIAP,
                                                                         CAAP, PCD, NR13, PAWR, PPP1R13B, IFI44, TLR,
                                                                         MyD88, CCAR, ASPP, GIMAP)

#subset only those that are significant
resRODTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED_significant <- resRODTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED %>%
  filter(Significance == "significant")
#add column with apoptosis pathwway and label genes with same name with extra (#)
write.csv(resRODTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED_significant, file="resRODTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED_significant.csv")
resRODTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED_significant_path <- read.csv(file="resRODTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED_significant.csv", header=TRUE)

####Graph Significantly differentially expressed proteins from ROD challenge, with log2 fold change assessed across time points ####
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#color based on their known pathway

#padj less than 10-5
label.ROD.df <- data.frame(description= c("G-protein coupled receptor 161-like", "caspase-8-like",
                                          "baculoviral IAP repeat-containing protein 7-like", "toll-like receptor 6 (1)",
                                          "toll-like receptor 1 (2)", "toll-like receptor 6 (2)",
"GTPase IMAP family member 4-like (1)", "GTPase IMAP family member 4-like (3)"), log2FoldChange=c(-22.5, -12.0,-3.0,-4.0,-5.0, -3.0,-5.5, 4.0 ))
ROD_Sig_PLOT <- ggplot(resRODTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED_significant_path) + 
  geom_col(aes(x=description, y=log2FoldChange, fill=Pathway)) + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + geom_hline(yintercept=0.0) +
  scale_y_continuous(name ="log2 Fold Change", breaks = scales::pretty_breaks(n = 20)) + 
  ggtitle("Significantly Differentially Expressed Apoptosis Genes in C. virginica after ROD Challenge") + 
  theme(plot.title = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5))

ROD_Sig_PLOT2 <- ROD_Sig_PLOT + geom_text(data=label.ROD.df, aes(x=description, y=log2FoldChange), label = c("***")) +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  scale_x_discrete(name="Description", limits=c("protein kinase C iota type-like",
"G-protein coupled receptor 161-like","receptor-interacting serine/threonine-protein kinase 4-like (1)","receptor-interacting serine/threonine-protein kinase 4-like (2)",
 "caspase-8-like","caspase-2-like", "TNF receptor-associated factor family protein DDB_G0272098-like","ceramide synthase 2-like",
"mitogen-activated protein kinase kinase kinase 7-like", "baculoviral IAP repeat-containing protein 7-like",
"caspase activity and apoptosis inhibitor 1-like","interferon-induced protein 44-like (1)",
"interferon-induced protein 44-like (2)","interferon-induced protein 44-like (3)","toll-like receptor 13",
 "toll-like receptor 10", "toll-like receptor 1 (1)", "toll-like receptor 6 (1)","toll-like receptor 1 (2)", 
"toll-like receptor 4 (1)","toll-like receptor 4 (2)","toll-like receptor 6 (2)","GTPase IMAP family member 4-like (1)",
"GTPase IMAP family member 4-like (2)", "GTPase IMAP family member 4-like (3)", "GTPase IMAP family member 4-like (4)",
"GTPase IMAP family member 7-like (1)","GTPase IMAP family member 7-like (2)"), 
label=c("protein kinase C iota type-like"="PKC iota type-like","G-protein coupled receptor 161-like" ="GPCR 161-like",
"receptor-interacting serine/threonine-protein kinase 4-like (1)"="RIPK 4-like","receptor-interacting serine/threonine-protein kinase 4-like (2)"="RIPK 4-like",
  "caspase-8-like"="caspase-8-like","caspase-2-like"="caspase-2-like","TNF receptor-associated factor family protein DDB_G0272098-like"="TRAF protein DDB_G0272098-like",
"ceramide synthase 2-like"="ceramide synthase 2-like","mitogen-activated protein kinase kinase kinase 7-like"="MEKK/MAPK 7-like",
 "baculoviral IAP repeat-containing protein 7-like"="BIR IAP repeat-containing protein 7-like","caspase activity and apoptosis inhibitor 1-like"="CAAP 1-like",
"interferon-induced protein 44-like (1)"="IFI44","interferon-induced protein 44-like (2)"="IFI44","interferon-induced protein 44-like (3)"="IFI44",
"toll-like receptor 13"="TLR13","toll-like receptor 10"="TLR10", "toll-like receptor 1 (1)"="TLR1", "toll-like receptor 6 (1)"="TLR6","toll-like receptor 1 (2)"="TLR1", 
"toll-like receptor 4 (1)"="TLR4","toll-like receptor 4 (2)"="TLR4","toll-like receptor 6 (2)"="TLR6","GTPase IMAP family member 4-like (1)"="GIMAP4-like",
"GTPase IMAP family member 4-like (2)"="GIMAP4-like", "GTPase IMAP family member 4-like (3)"="GIMAP4-like", "GTPase IMAP family member 4-like (4)"="GIMAP4-like",
"GTPase IMAP family member 7-like (1)"="GIMAP7-like","GTPase IMAP family member 7-like (2)"="GIMAP7-like"))


#### RI PROBIOTIC CHALLENGED Differential Gene Expression Analysis ####

#Subset Data for RIF Transcriptomes 
#Extract columns you want from the C_vir_TranscriptCountData, based on which column the correct SRA data

RIF_C_vir_TranColData <- C_vir_TranColData[16:21,]
RIF_C_vir_TranscriptCountData <- C_vir_TranscriptCountData[,16:21]
head(RIF_C_vir_TranscriptCountData)

#give the treatment column levels
RIF_C_vir_TranColData$treatment <- factor(RIF_C_vir_TranColData$treatment)
levels(RIF_C_vir_TranColData$treatment)

#### Create RIF DESeq Data Set from Matrix ####
# DESeqDataSet from count matrix and labels, separate into resistant and susceptible 
#add an interaction term to compare treatment between two conditions 
#layout used for interactions: https://support.bioconductor.org/p/58162/

ddsRIFTran <- DESeqDataSetFromMatrix(countData = RIF_C_vir_TranscriptCountData, 
                                     colData = RIF_C_vir_TranColData, 
                                     design =  ~ treatment)

ddsRIFTran<- ddsRIFTran[ rowSums(counts(ddsRIFTran)) > 1, ]

# review how the data set looks
head(ddsRIFTran)

#Relevel each to make sure that control is the first level in the treatment factor for each
ddsRIFTran$treatment <- relevel( ddsRIFTran$treatment, "control")

#Check we're looking at the right samples
as.data.frame( colData(ddsRIFTran) )

#Running the DEG pipeline
ddsRIFTran<- DESeq(ddsRIFTran) #for designs with interactions, recommends setting betaPrior=FALSE

#Inspect results
#extract contrasts between control and treatment values for interaction
resRIFTran<- results(ddsRIFTran)
head(resRIFTran)

#summary is just printing a table for you, you need to tell it what threshold you want
help("summary",package="DESeq2")
alpha <- 0.05 #set alpha to 0.05, this will control FDR
summary(resRIFTran) #default FDR is still 0.1
summary(resRIFTran, alpha) #now showing all genes with FRD < 0.05

#To get the significant genes
#The independent filtering in results() has an argument 'alpha'
#which is used to optimize a cutoff on mean normalized count
#to maximize the number of genes with padj < alpha
resRIFTran_05 <- results(ddsRIFTran, alpha= alpha) #set FDR to 0.05 now
resRIFTran_05_Sig <- resRIFTran[which(resRIFTran$padj < alpha),]
summary(resRIFTran_05) #this is all the genes
summary(resRIFTran_05_Sig) #this is the significant ones!
sum(resRIFTran_05$padj < 0.05, na.rm=TRUE) #1179 tells you how many genes have expected FDR ≤ 0.05
sum(resRIFTran_05_Sig$padj < 0.05, na.rm=TRUE) #1179
resRIFTran_05_Sig$Significance <- sig
resRIFTran_05_nonSig <- resRIFTran[which(resRIFTran$padj > alpha),] #create list of nonsig
nonsig <- "non-significant"
resRIFTran_05_nonSig$Significance <- nonsig
head(resRIFTran_05_Sig)
head(resRIFTran_05_nonSig)

#add ID column with the rownames so a merge can happen later
resRIFTran_05_Sig["ID"] <- rownames(resRIFTran_05_Sig) #add a new column with the rownames for match
resRIFTran_05_nonSig["ID"] <- rownames(resRIFTran_05_nonSig)
resRIFTran_05["ID"] <- rownames(resRIFTran_05)

#metadata on meaning of the columns
mcols(resRIFTran_05_Sig, use.names = TRUE)
mcols(resRIFTran_05_Sig)$description
#shows treatment vs. control with baseMean as the mean of normalized counts for all samples

#Order by Log2FC
head( resRIFTran_05[ order( resRIFTran_05$log2FoldChange ), ] ) #head for strongest downregulation
tail( resRIFTran_05[ order( resRIFTran_05$log2FoldChange ), ] ) #tail for strongest up regulation

####p-value correction for both sets of results ####
#First Visualize histograms
#histogram of P- values to visualize any "hills" or "U shape"
# hill means variance of the null distribution too high, U shape means variance assumed too low
hist(resRIFTran_05$pvalue, breaks= 20, col = "grey")
hist(resRIFTran_05_Sig$pvalue, breaks = 20, col = "grey") #hill

#looks good, don't need to do the correction

#Visualize Results with Diagnostic Plots#
#MA plot, useful overview for experiment with two-group comparison. Plots log2FC over mean of normalized counts
#genes with adjusted p value 
plotMA(resRIFTran_05)
plotMA(resRIFTran_05_Sig)

#Export Results to CSV
write.csv( as.data.frame(resRIFTran_05), file="RIF_resRIFTran_05.csv")
write.csv( as.data.frame(resRIFTran_05_Sig), file="RIF_resRIFTran_05_Sig.csv")
write.csv( as.data.frame(resRIFTran_05_nonSig), file = "RIF_resRIFTran_05_nonSig.csv")

####Match lines for only those that have a "LOC" in the C_vir_stringtie_transcripts_rna_LOC_separated$transcript_id ####
# %in% returns true for every value in the first argument that matches a value in the second argument.
# the order of arguments is important
#load file with the spaces changed to ; so that the column can be separated
C_vir_stringtie_transcripts_rna_LOC_separated_clean2 <- read.csv(file="C_vir_stringtie_transcripts_rna_LOC_separated_clean2.csv", header=TRUE)

#separate the transcript id column with the semicolon
C_vir_stringtie_transcripts_rna_LOC_separated_clean_transcript_id <- C_vir_stringtie_transcripts_rna_LOC_separated_clean2 %>% separate(transcript_id, c("transcript_id", "ID"), ";", extra="merge")
#check duplicate rownames
RIF_dupliate_names<- duplicated(resRIFTran_05_Sig$rownames)
grep("TRUE", RIF_dupliate_names) #0 duplicates
#must comvert to data frame
resRIFTran_05_Sig_df <- data.frame(resRIFTran_05_Sig)
#Merge columns together based on match in the "ID" column 
resRIFTran_05_Sig_df_FULL <- merge(resRIFTran_05_Sig_df, C_vir_stringtie_transcripts_rna_LOC_separated_clean_transcript_id[,c("ID", "gene_name")], by="ID")
nrow(resRIFTran_05_Sig_df_FULL) #517
#put LOC info in new column
resRIFTran_05_Sig_df_FULL <- resRIFTran_05_Sig_df_FULL %>% separate(gene_name, c("gene_name", "gene_ID"), ";")

#Repeat for non significant genes
resRIFTran_05_nonSig_df <- data.frame(resRIFTran_05_nonSig)
resRIFTran_05_df_nonSig_FULL <- merge(resRIFTran_05_nonSig_df, C_vir_stringtie_transcripts_rna_LOC_separated_clean_transcript_id[,c("ID", "gene_name")], by="ID")

#put LOC info in new column
resRIFTran_05_df_nonSig_FULL <- resRIFTran_05_df_nonSig_FULL %>% separate(gene_name, c("gene_name", "gene_ID"), ";")

####Lookup LOC values using Batch Entrez ####
#Sig
write.csv(resRIFTran_05_Sig_df_FULL$gene_ID, file="resRIFTran_05_Sig_df_FULL_gene_ID.csv")
#put column with LOC info into text file
#perform batch Entrez lookup to the gene database to retrieve IDs

#non_Sig
write.csv(resRIFTran_05_df_nonSig_FULL$gene_ID, file="resRIFTran_05_df_nonSig_FULL_gene_ID.csv")
#put column wit LOC into text file
#perform batch Entrez lookup with the gene database to retrieve IDs
#batch Entrez gets rid of duplicates

#####Merge the LOC values with gene name####
#SIG
#convert file to csv using excel 
resRIFTran_05_Sig_ENTREZGENE <- read.csv(file="resRIFTran_05_df_Sig_ENTREZ_RESULTS_FULL.csv", header = TRUE)
#change symbol to gene_ID
colnames(resRIFTran_05_Sig_ENTREZGENE)[6] <- "gene_ID"
#merge file based on the LOC column
resRIFTran_05_Sig_ENTREZGENE_DATA <-  merge(resRIFTran_05_Sig_ENTREZGENE, resRIFTran_05_Sig_df_FULL[,c("gene_ID", "log2FoldChange","pvalue","padj")], by="gene_ID")
head(resRIFTran_05_Sig_ENTREZGENE_DATA)

#NON SIG
resRIFTran_05_NON_Sig_ENTREZGENE <- read.csv(file="RIF_non_sig_ENTREZ_all.csv", header = TRUE)
#change symbol to gene_ID
colnames(resRIFTran_05_NON_Sig_ENTREZGENE)[6] <- "gene_ID"
#merge file based on the LOC column
resRIFTran_05_NON_Sig_ENTREZGENE_DATA <-  merge(resRIFTran_05_NON_Sig_ENTREZGENE, resRIFTran_05_df_nonSig_FULL[,c("gene_ID", "log2FoldChange","pvalue","padj")], by="gene_ID")
head(resRIFTran_05_Sig_ENTREZGENE_DATA)

#add columns for significance and merge the two files so only one search needs to be done
resRIFTran_05_Sig_ENTREZGENE_DATA$Significance <- sig
resRIFTran_05_NON_Sig_ENTREZGENE_DATA$Significance <- nonsig

resRIFTran_05_FULL_ENTREZGENE_DATA <- rbind(resRIFTran_05_Sig_ENTREZGENE_DATA,resRIFTran_05_NON_Sig_ENTREZGENE_DATA)

####Add full list of apoptosis related genes and aliases to grep for
#https://stackoverflow.com/questions/17130129/find-matches-of-a-vector-of-strings-in-another-vector-of-strings
#https://stackoverflow.com/questions/37684179/how-to-search-for-multiple-words-in-the-same-string-using-r

# full data frame
#make the description column a character
resRIFTran_05_FULL_ENTREZGENE_DATA$description <- as.character(resRIFTran_05_FULL_ENTREZGENE_DATA$description)

# Sample vector of keywords or phrases
#Search terms come from the strings that found hits in the C_virginica genome annotation view in NCBI:
#all terms with no hits in the annotation on NCBI were INCLUDED IN THIS LIST..just in case, and all aliases searched for

apoptosis_keywords <- as.character(c("IP3R","inositol","RhoA","RhoB","RhoC","rho","RhoE","ras","Rac1","Cdc42",
                                     "PKC","IkB","kappa","NFKB","NF-kappa-B","GPCR","Apo3","TWEAK","DR3LG","TNFSF12","Apo","TNFSF10","TNFSF","TRAIL",
                                     "Apo2L","c-FLIP","FLICE", "Casper","FLAME","FLAME-1","CASH","CLARP","MRIT","FADD-like","IL-1β","I-FLICE",
                                     "RIP","RIPK1","RIPK","TNFR-interacting","TNFR","AMP","adenosine","cyclic","kinase","PKC","adenylate","CREB",
                                     "TNFRSF12","Apo3","TRAMP","DDR3","LARD","WSL","TNFRSF","DR4","TNFRSF10A","TRAILR1",
                                     "TRAILR","APO2","TNFRS10B","TRAIL-R2","TRICK2","ZTNFR9","KILLER","Fas","FASL","TNFSF6","FasR",
                                     "TNF-alpha","alpha","Fas","FAIM","TNF","TRADD","TNFRSF1A","Fas","casp","casp2","casp8","caspase",
                                     "procasp8","caspase","caspase","procaspase","procasp","siglec","immunoglobulin","sialic",
                                     "Cgsiglec-1","IFNLP","IFN","interferon","interferon-like","TNF","TNFRSF1A","TRAF",
                                     "mitotic","ICAD","DFF-45","Dnase","CAD","DFF-40","CPAN","caspase","casp6","caspase","caspase",
                                     "casp3/7","caspase","casp7","caspase","casp3","caspase","casp14","caspase",
                                     "caspase","casp9","caspase","APAF","APIP","apoptotic","APAF-1","protease","SAC",
                                     "soluble","AC","BI-1", "bcl-2", "bcl","PDRP","p53","cyto","apoptosis","APAF-1","APAF",
                                     "procaspase","p55-BLK","Blk","MGC10442","myc","oncogene","aven","NY-CO-13","trp53","p53",
                                     "smac","DIABLO","IAP","IAP-binding","xL","bcl2","bcl-2", "xS","bcl-w","myeloid","mcl-1",
                                     "mcl","leukemia","A1","Bak","BCL2L7","antagonist", "Bcl-XL","BAD","BAG","BH3","Bid",
                                     "tBid","BIM","Bik", "JFY1","puma","ceramide", "MEK","MAPK","p38","phosphatidyl","cyclin",
                                     "JUN","p39","c-jun","AP-1","proto-oncogene", "AP1","JUN","CD151","BTG","caspase",
                                     "IAP","apoptosis","CAAP","caspase","NR13","bcl-2","PRKC","WT1","PPP1R13B",
                                     "p53","IFI44","interferon","IFI6","TLR","MyD88","myeloid","CCAR", "GIMAP", "ovarian", "bok"))

apoptosis_multiple_keywords <- as.character(c("protein kinase c","NF-kappa-B inhibitor", "nuclear factor","G protein","Apo3 ligand",
                                              "Apo2 ligand","Receptor interacting","threonine kinase", "serine threonine kinase","Protein kinase A","adenylate cyclase","cyclic responsive",
                                              "cAMP responsive", "death receptor","death receptor 4", "death receptor", "death receptor 5","fatty acid synthase","tumor necrosis factor",
                                              "tumor necrosis factor","tumor necrosis factor","apoptotic inhibitory","death domain","receptor associated","Fas associated death domain","mitotic apparatus",
                                              "nuclear mitotic","SP-H antigen","adenylyl cyclase","bax inhibitor","DNA damage","cytochrome c","apoptosis inducing","peptidase activity","protease activity",
                                              "apoptotic protease","apoptotic pedtidase","B lymphoid", "tyrosine kinase","tumor antigen","tumor protein","low PI","mitochondrial activator",
                                              "activator caspase","cell differentiation","cell death","antagonist killer","bax","bcl-2","ovarian","bok","ovarian killer","associated death promoter",
                                              "binding component","tBid","interacting killer","modulator of apoptosis", "MEK Kinase", "endonuclease G", "akt-kinase","protein kinase B",
                                              "kinase B","N-terminal kinase","cyclin dependent cyclase 5 activator","cyclase 5","proto-oncogene","transcription factor","B cell translocation",
                                              "programmed cell death","WT1 regulator", "apoptosis regulator", "apoptosis-stimulating","cell division","apoptosis regulator","death agonist"))


#remove duplicates
apoptosis_keywords_unique <- unique(unlist(strsplit(apoptosis_keywords, " ")))
apoptosis_multiple_keywords_no_space_unique <- unique(unlist(strsplit(apoptosis_multiple_keywords_no_space, " ")))

#searching for each one individually:
RIF_IP3R <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("inositol", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & 
                                             grepl("phosphate", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &  grepl("receptor", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_IP3Ra <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("IP3R", resRIFTran_05_FULL_ENTREZGENE_DATA$description),] #8 
RIF_Rho <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("ras", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) & 
                                            grepl("GTP", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) &  grepl("rho", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_Rac1 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("ras", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) & 
                                             grepl("botulinum", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) &  grepl("C3", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_cdc42 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("cdc42", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                              !grepl("activated", ignore.case=TRUE,resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_PKCiota <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("protein", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & 
                                                grepl("kinase", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl("iota", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_PKCdelta <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("protein", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & 
                                                 grepl("kinase", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl("delta", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                                 !grepl("calmodulin", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & !grepl("ribosomal", resRIFTran_05_FULL_ENTREZGENE_DATA$description) ,]

RIF_IkB <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("NF-kappa-B", ignore.case= TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) & 
                                            grepl( "inhibitor", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) & !grepl("ras-like", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_NFkB <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("NF-kappa-B", ignore.case= TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) & 
                                             grepl( "nuclear", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) & !grepl("inhibitor", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_GPCRsevenpass <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("G-type", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl("seven-pass", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_GPCR <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("G-protein" , resRIFTran_05_FULL_ENTREZGENE_DATA$description) & 
                                             grepl("coupled", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &  grepl("receptor", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
# 0 hits
RIF_Apo3 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("Apo3",ignore.case= TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
# 0 hits
RIF_Apo2 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("Apo2",ignore.case= TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
# 0 hits
RIF_ApoL <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("Apo" , ignore.case= TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl("ligand", ignore.case= TRUE,resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
# 0 hits
RIF_TWEAK <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("TWEAK",ignore.case= TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
# 0 hits
RIF_TRAIL <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("TRAIL",ignore.case= TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
# 0 hits
RIF_TRAILb <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("TNF-related", ignore.case= TRUE,  resRIFTran_05_FULL_ENTREZGENE_DATA$description) & 
                                               grepl("apoptosis",ignore.case= TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) &  grepl("ligand",ignore.case= TRUE,  resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
#0 hits
RIF_TRAILc <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("TNFSF12",ignore.case= TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
# 0 hits
RIF_cFLIP <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("c-FLIP",ignore.case= TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
# 0 hits
RIF_cFLIP2 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("FLICE",ignore.case= TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_cFLIP3 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("CASPER",ignore.case= TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_cFLIP4 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("FLAME",ignore.case= TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_cFLIP5 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("CASH",ignore.case= TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_cFLIP6 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("CLARP",ignore.case= TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_cFLIP7 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("MRIT",ignore.case= TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_cFLIP8 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("FADD-like",ignore.case= TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                               grepl("inhibitory",ignore.case= TRUE,  resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_RIP <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("receptor", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl("interacting", resRIFTran_05_FULL_ENTREZGENE_DATA$description)
                                          & !grepl("thyroid", resRIFTran_05_FULL_ENTREZGENE_DATA$description)
                                          & !grepl("glutamate", resRIFTran_05_FULL_ENTREZGENE_DATA$description)
                                          & !grepl("cannabinoid", resRIFTran_05_FULL_ENTREZGENE_DATA$description)
                                          & !grepl("AH", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
#non are cAMP specifically
RIF_cAMP <-  resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("cAMP",ignore.case= TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#no hits
RIF_cAMP2 <-  resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("cyclic",ignore.case= TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) 
                                             & grepl("adenosine", ignore.case = TRUE,resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
#none are PKA
RIF_PKA <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("protein", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & 
                                            grepl("kinase", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl("A", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_AC <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("adenylate", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl("cyclase", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_CREB <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("cAMP-responsive", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#0 hits
RIF_DR <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("death",ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl("receptor", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_DR1 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("TRAMP", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_DR2 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("WSL",ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_DR3 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("LARD",ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_DR4 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("DDR3",ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_Fas <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("fatty", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl( "acid", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl( "synthase", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                            !grepl("acyl",resRIFTran_05_FULL_ENTREZGENE_DATA$description),]


RIF_FasL <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("FasL", ignore.case = TRUE ,resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_FASL2 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("TNF", ignore.case = TRUE ,resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                              grepl("superfamily", ignore.case = TRUE ,resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                              grepl("6", ignore.case = TRUE ,resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
#none are what I want
RIF_Fas2 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("Fas", ignore.case = TRUE ,resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#no DR5
RIF_DR5 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("TNFRS10B", ignore.case = TRUE ,resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_DR5a <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("TRAIL-R2", ignore.case = TRUE ,resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_DR5b <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("TRICK2", ignore.case = TRUE ,resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_DR5c <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("KILLER", ignore.case = TRUE ,resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_DR5d <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("ZTNFR9", ignore.case = TRUE ,resRIFTran_05_FULL_ENTREZGENE_DATA$description),]


RIF_TNF <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("tumor", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & 
                                            grepl( "necrosis", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl( "factor", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                            !grepl("C1q", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                            !grepl("alpha-induced",resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                              !grepl("endothelial",resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_FAIM <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("fas", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & 
                                             grepl("apoptotic", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("inhibitory", resRIFTran_05_FULL_ENTREZGENE_DATA$description) ,]

#non are TRAD
RIF_TRADD <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("death", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl(  "domain", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                              grepl(  "receptor", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_FADD <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("FAS-associated", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("death", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_caspase <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("caspase", ignore.case = TRUE ,resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                                !grepl("activity",resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_procaspase <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("procaspase", ignore.case = TRUE ,resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#no sialic acid
RIF_siglec <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("siglec", ignore.case = TRUE ,resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_sialic <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("sialic", ignore.case = TRUE ,resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                               grepl("immunoglobulin", ignore.case = TRUE ,resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
#none are IFNLP
RIF_IFNLP <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("IFNLP", ignore.case = TRUE ,resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_IFNLPa <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("interferon", ignore.case = TRUE ,resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_TRAF <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("TNF", ignore.case = TRUE ,resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("receptor-associated", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) ,]

#no NuMA
RIF_NuMa <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("mitotic", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "apparatus", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_NuMAa <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("nuclear", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "mitotic", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_NuMAb<- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("SP-H", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_ICAD <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("DNA", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl("fragmentation", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("alpha", resRIFTran_05_FULL_ENTREZGENE_DATA$description) ,]
RIF_CAD <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("DNA", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl("fragmentation", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl("beta", resRIFTran_05_FULL_ENTREZGENE_DATA$description) ,]
RIF_CADb <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("caspase", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl("activated", resRIFTran_05_FULL_ENTREZGENE_DATA$description) ,]

#no APIP
RIF_APIP <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("apoptotic", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "peptidase", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_APIPb <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("apoptotic", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "protease", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_APIPc <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("APAF", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) ,]
RIF_APIPc <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("APIP", ignore.case=TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) ,]

#no sAC
RIF_sAC <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("adenylyl", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "cyclase", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#no BI-1
RIF_BI <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("bax", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "inhibitor", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_PDRP <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("p53", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "damage", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_cytoc <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("cytochrome", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "c-like", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_AIF <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("apoptosis", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & 
                                            grepl( "inducing", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl( "factor", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
#no Blk
RIF_Blk <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("blk", ignore.case=TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_Blkb <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("B", ignore.case=TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "lymphoid", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_Blkc <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("p55", ignore.case=TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#no myc
RIF_myc <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("myc", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_myc <-  resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("oncogene", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#no aven
RIF_aven <-   resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("aven", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                               !grepl("scavenger", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
#no p53..
RIF_p53 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("tumor", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "antigen", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_p53 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("p53", resRIFTran_05_FULL_ENTREZGENE_DATA$description) ,]
RIF_p53 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("tumor", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "protein", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#no smac
RIF_smac <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("low", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "PI", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_smac <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("mitochondrial", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "activator", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_smac1 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("IAP", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "direct", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_bcl <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("bcl",ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                                !grepl("interacting",resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#no mcl
RIF_mcl <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("myeloid",ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl("leukemia", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_mcl <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("cell", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "differentiation", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#no A1
RIF_A1 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("bcl-2",ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#no Bak
RIF_Bak <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("antagonist", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( " killer", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_Bak2 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("bcl2l7",ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_BAX <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("bax-like",ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#no Bok
RIF_Bok <-  resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("ovarian", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "killer", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_Bok <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("ovarian", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#no BAD
RIF_BAD <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("death", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "promoter", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_BAD2 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("binding", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "component", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_BAD3 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("BAD", ignore.case=TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_BAG <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("BAG",ignore.case=TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#no BH3
RIF_BH3 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("BH3",ignore.case=TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_Bid <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("Bid",ignore.case=TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#no BIM
RIF_BIM <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("BIM",ignore.case=TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

# no Bik
RIF_Bik <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("interacting", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "killer", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#no Puma
RIF_Puma <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("modulator", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "apoptosis", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_Puma <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("JFY", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_ceramide <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("ceramide", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                                 !grepl("galactosylceramide",resRIFTran_05_FULL_ENTREZGENE_DATA$description)&
                                                 !grepl("kinase",resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                                 grepl("synthase",resRIFTran_05_FULL_ENTREZGENE_DATA$description ) &
                                                 !grepl("phospho", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_MEKK_MAPK <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("mitogen-activated", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                                  grepl( "kinase", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                                  !grepl("kinase-binding", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                                    !grepl("interacting", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#no p38 MAPK
RIF_p38_MAPK <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("mitogen-activated", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                                 grepl( "kinase", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                                 grepl("p38", ignore.case=TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_endoG <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("endonuclease", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                                  grepl( "G", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                                  !grepl("flap", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_PI3 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("phosphatidylinositol", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "3-kinase", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#only AKT interacting protein
RIF_AKTB <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("protein", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("kinase", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("B",resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
#no protein kinase B!
RIF_akt <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("akt", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_JNK <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("stress-activated", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "JNK", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]


RIF_p35 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("cyclin", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "dependent", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl( "5", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl("activator",resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
#none
RIF_cJUN <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("proto", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl( "oncogene", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("jun", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
#none
RIF_cJUN <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("c-jun", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
#1 
RIF_cJUNa <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("AP-1", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                              grepl("transcription",resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_CD151 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("CD151", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#none
RIF_BTG1 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("B", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl("translocation", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#1 
RIF_BTG2 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("BTG", resRIFTran_05_FULL_ENTREZGENE_DATA$description) ,]

RIF_IAP <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("inhibitor", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl("apoptosis", ignore.case=TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                            !grepl("death-associated", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                            !grepl("caspase", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                            !grepl("5-like", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_IAP2 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("baculoviral", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("IAP", ignore.case = TRUE,resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
RIF_API <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("inhibitor", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl("apoptosis", ignore.case=TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                            !grepl("death-associated", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                            !grepl("caspase", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl("5-like", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_DIAP <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("inhibitor", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("apoptosis", ignore.case=TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("death-associated", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
#no CAAP
RIF_CAAP <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("inhibitor", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("apoptosis", ignore.case=TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                             !grepl("death-associated", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("caspase", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_PCD <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("programmed", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & 
                                            grepl("death", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                            grepl("cell", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                            !grepl("interacting", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_NR13 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("NR13", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]


RIF_PAWR <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("WT1", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl("regulator", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_IFI44 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("interferon", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl("44", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#no IFI6
RIF_IFI6 <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("interferon", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl("6", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_TLR <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("toll", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl("receptor", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

RIF_MyD88 <-resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("myeloid", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl("differentiation", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#no CCAR
RIF_CCAR <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("cell", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl("division", resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                             grepl("apoptosis",resRIFTran_05_FULL_ENTREZGENE_DATA$description),]
#no ASPP
RIF_ASPP <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("apoptosis-stimulating", resRIFTran_05_FULL_ENTREZGENE_DATA$description) & grepl("p53", resRIFTran_05_FULL_ENTREZGENE_DATA$description),]

#
RIF_GIMAP <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("GIMAP", resRIFTran_05_FULL_ENTREZGENE_DATA$description) ,]
RIF_GIMAP <- resRIFTran_05_FULL_ENTREZGENE_DATA[grepl("GTPase", ignore.case = TRUE, resRIFTran_05_FULL_ENTREZGENE_DATA$description) &
                                              grepl("IMAP",resRIFTran_05_FULL_ENTREZGENE_DATA$description),]


#Combine all lists from above into one new big table
#only include those that have hits and are confirmed to be correct
resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED <- rbind(RIF_IP3R, RIF_Rho, RIF_Rac1, RIF_cdc42, RIF_PKCiota, RIF_PKCdelta, RIF_IkB,
                                                                         RIF_NFkB, RIF_GPCRsevenpass, RIF_GPCR, RIF_RIP, RIF_AC, RIF_CREB, RIF_Fas,
                                                                         RIF_TNF, RIF_FAIM,RIF_FADD,RIF_caspase, RIF_TRAF, RIF_ICAD, RIF_CAD, RIF_cytoc,
                                                                         RIF_AIF, RIF_bcl, RIF_BAX, RIF_BAG, RIF_ceramide, RIF_MEKK_MAPK, RIF_endoG,
                                                                         RIF_PI3, RIF_p35, RIF_CD151, RIF_BTG2, RIF_IAP, RIF_IAP2, RIF_API, RIF_DIAP,
                                                                         RIF_PCD, RIF_PAWR, RIF_TLR, RIF_IFI44, RIF_MyD88, RIF_GIMAP)

#subset only those that are significant
resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED_significant <- resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED %>%
  filter(Significance == "significant")
#add column with apoptosis pathwway and label genes with same name with extra (#)
write.csv(resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED_significant, file="resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED_significant.csv")
resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED_significant_path <- read.csv(file="resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED_significant.csv", header=TRUE)

####Graph Significantly differentially expressed proteins from ROD challenge, with log2 fold change assessed across time points ####
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#color based on their known pathway

#padj less than 10^-5
label.RIF.df <- data.frame(description= c("mitogen-activated protein kinase kinase kinase 7-like",
                                          "GTPase IMAP family member 4-like"), log2FoldChange=c(-13.5, 11.0))
RIF_Sig_PLOT <- ggplot(resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED_significant_path) + 
  geom_col(aes(x=description, y=log2FoldChange, fill=Pathway)) + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + geom_hline(yintercept=0.0) +
  scale_y_continuous(name ="log2 Fold Change", breaks = scales::pretty_breaks(n = 20)) + 
  ggtitle("Significantly Differentially Expressed Apoptosis Genes in C. virginica after RIF Challenge") + 
  theme(plot.title = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5))

RIF_Sig_PLOT2 <- RIF_Sig_PLOT + geom_text(data=label.RIF.df, aes(x=description, y=log2FoldChange), label = c("***")) +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  scale_x_discrete(name="Description", limits=c("leucine-rich repeat-containing G-protein coupled receptor 4-like",
                                                "receptor-interacting serine/threonine-protein kinase 4-like",
                                                "fatty acid synthase-like",
                                                 "mitogen-activated protein kinase kinase kinase 7-like",
                                                "GTPase IMAP family member 4-like"), 
                   label=c("leucine-rich repeat-containing G-protein coupled receptor 4-like"="GPCR 4-like",
                           "receptor-interacting serine/threonine-protein kinase 4-like"="RIPK 4-like",
                           "fatty acid synthase-like"="Fas-like",
                           "mitogen-activated protein kinase kinase kinase 7-like"="MEKK/MAPK 7-like",
                           "GTPase IMAP family member 4-like"="GIMAP4-like"))



#references: 
#https://github.com/sr320/LabDocs/tree/master/code/DESeq
#StringTie manual : http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4743157/
