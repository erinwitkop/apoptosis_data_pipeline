#05_C_gig_Bac_Viral_DESeq2_transcript_analysis.R

####INPUT GATHERED FROM CRASSOSTREA gigas ####


#This script takes as input the output Bac_Viral_transcript_count_matrix.csv data prepared from prepDE.py and performs
#differential transcript expression analysis, and subsets out isoforms of GIMAPs and CgIAPs and graphs their
#relative abundance.

#call the DESeq2 library 
#source("http://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade") 
#biocLite("DESeq2")
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
#install.packages("tm")
library(tm)
library(stringr)
# Construct Bac_Viral_PHENO_DATA.csv that contains SRA run information, such as which contrast, tissue, etc.

####Match "rna#" value with the Gene LOC name in the stringtie file
#Read in stringtie merged file
C_gig_stringtie <- read.csv(file= "Bac_viral_stringtie_merged.gtf", sep="\t", header = FALSE)
#set colnames for the attributes
colnames(C_gig_stringtie) <- c('seqname', 'source', 'feature','start','end','score','strand','frame', 'attribute')
#subset for just transcript lines
C_gig_stringtie_transcripts <- C_gig_stringtie %>% filter(feature =="transcript")
#grep out only the lines with "rna" in the line
C_gig_stringtie_transcripts_rna <-filter(C_gig_stringtie_transcripts, grepl("rna", attribute))
#grep out only the lines with "rna" and a "LOC" value
C_gig_stringtie_transcripts_rna_LOC <- filter(C_gig_stringtie_transcripts_rna, grepl("LOC", attribute))
#make into df
C_gig_stringtie_transcripts_rna_LOC <- as.data.frame(C_gig_stringtie_transcripts_rna_LOC)
#separate attributes column with the ; separator in tidyr
#use the extra  argument to control what happens when every row doesn't split into the same number of pieces
C_gig_stringtie_transcripts_rna_LOC_separated <- C_gig_stringtie_transcripts_rna_LOC %>% separate(attribute, c("gene_id","transcript_id","gene_name","ref_gene_id"), ";", extra="merge")   
#clean up white space in the data
C_gig_stringtie_transcripts_rna_LOC_separated_clean <- data.frame(lapply(C_gig_stringtie_transcripts_rna_LOC_separated, trimws),stringsAsFactors = FALSE)
#remove spaces in terminal
# tr ' ' ';'<C_gig_stringtie_transcripts_rna_LOC_separated_clean.csv >C_gig_stringtie_transcripts_rna_LOC_separated_clean2.csv
write.csv(C_gig_stringtie_transcripts_rna_LOC_separated_clean, file ="C_gig_stringtie_transcripts_rna_LOC_separated_clean.csv")

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
oshv1TranColData <- read.csv("OsHV1_PHENO_DATA.csv", header=TRUE)
print(oshv1TranColData)
rownames(oshv1TranColData) <- oshv1TranColData$sampleID
colnames(oshv1TranCountData) <- oshv1TranColData$sampleID
head(oshv1TranCountData)
head(oshv1TranColData)

#Give the condition column levels
oshv1TranColData$condition <- factor(oshv1TranColData$condition)
levels(oshv1TranColData$condition) #check to see that it has levels 

#Give the treatment column levels
oshv1TranColData$treatment <- factor(oshv1TranColData$treatment)
levels(oshv1TranColData$treatment)

# Check all sample IDs in oshv1ColData are also in oshv1CountData and match their orders
all(rownames(oshv1TranColData) %in% colnames(oshv1TranCountData))  #Should return TRUE
# returns TRUE
all(rownames(oshv1TranColData) == colnames(oshv1TranCountData))    # should return TRUE
#returns TRUE

#### Create oshv1 DESeq Data Set from Matrix ####
# DESeqDataSet from count matrix and labels, separate into resistant and susceptible 
#add an interaction term to compare treatment between two conditions 
#layout used for interactions: https://support.bioconductor.org/p/58162/

#if I include the interactions here that means I am analyzing which genes are significantly
#different when you look at control versus treatment in each group and then over time.. 
# so if I look at the apoptosis genes that are important over that whole time frame then I 
# am really looking at those genes that are the most different across all time points 
# which really is what I care about
ddsoshvTran <- DESeqDataSetFromMatrix(countData = oshv1TranCountData, 
                                     colData = oshv1TranColData, 
                                     design =  ~ condition + treatment + condition:treatment)

ddsoshvTran<- ddsoshvTran[ rowSums(counts(ddsoshvTran)) > 1, ]

# review how the data set looks
head(ddsoshvTran)

#Relevel each to make sure that control is the first level in the treatment factor for each
ddsoshvTran$condition <- relevel(ddsoshvTran$condition, "A")
ddsoshvTran$treatment <- relevel( ddsoshvTran$treatment, "control")

#Check we're looking at the right samples
as.data.frame( colData(ddsoshvTran) )

#Running the DEG pipeline
ddsoshvTran<- DESeq(ddsoshvTran, betaPrior = FALSE) #for designs with interactions, recommends setting betaPrior=FALSE

#Inspect results
#extract contrasts between control and treatment values for interaction
resoshvTran<- results(ddsoshvTran)
head(resoshvTran)

# this last line is all that is needed because the interaction term is last
# in the design formula (said my Michael Love)
#Small p-values for the interaction term indicate that the log fold change
# due to treatment is significantly different for the two conditions.

#summary is just printing a table for you, you need to tell it what threshold you want
help("summary",package="DESeq2")
alpha <- 0.05 #set alpha to 0.05, this will control FDR
summary(resoshvTran) #default FDR is still 0.1
summary(resoshvTran, alpha) #now showing all genes with FRD < 0.05

#To get the significant genes
#The independent filtering in results() has an argument 'alpha'
#which is used to optimize a cutoff on mean normalized count
#to maximize the number of genes with padj < alpha
resoshvTran_05 <- results(ddsoshvTran, alpha= alpha) #set FDR to 0.05 now
resoshvTran_05_Sig <- resoshvTran[which(resoshvTran$padj < alpha),]
summary(resoshvTran_05) #this is all the genes
summary(resoshvTran_05_Sig) #this is the significant ones!
sum(resoshvTran_05$padj < 0.05, na.rm=TRUE) #1148 tells you how many genes have expected FDR ≤ 0.05
sum(resoshvTran_05_Sig$padj < 0.05, na.rm=TRUE) #1168
sig="significant"
resoshvTran_05_Sig$Significance <- sig
resoshvTran_05_nonSig <- resoshvTran[which(resoshvTran$padj > alpha),] #create list of nonsig
nonsig <- "non-significant"
resoshvTran_05_nonSig$Significance <- nonsig

#add ID column with the rownames so a merge can happen later
resoshvTran_05_Sig["ID"] <- rownames(resoshvTran_05_Sig) #add a new column with the rownames for match
resoshvTran_05_nonSig["ID"] <- rownames(resoshvTran_05_nonSig)
resoshvTran_05["ID"] <- rownames(resoshvTran_05)
resoshvTran_05_sig_non_sig <- rbind(resoshvTran_05_Sig, resoshvTran_05_nonSig)

#metadata on meaning of the columns
mcols(resoshvTran_05_Sig, use.names = TRUE)
mcols(resoshvTran_05_Sig)$description
#shows treatment vs. control with baseMean as the mean of normalized counts for all samples

#Order by Log2FC
head( resoshvTran_05_sig_non_sig[ order( resoshvTran_05_sig_non_sig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resoshvTran_05_sig_non_sig[ order( resoshvTran_05_sig_non_sig$log2FoldChange ), ] ) #tail for strongest up regulation

####p-value correction for both sets of results ####
#First Visualize histograms
#histogram of P- values to visualize any "hills" or "U shape"
# hill means variance of the null distribution too high, U shape means variance assumed too low
hist(resoshvTran_05$pvalue, breaks= 20, col = "grey")
hist(resoshvTran_05_Sig$pvalue, breaks = 20, col = "grey") #hill

#looks okay, not doing the correction because we want to avoid extra false positives

#Visualize Results with Diagnostic Plots#
#MA plot, useful overview for experiment with two-group comparison. Plots log2FC over mean of normalized counts
#genes with adjusted p value 
plotMA(resoshvTran_05)
plotMA(resoshvTran_05_Sig)

#Export Results to CSV
write.csv( as.data.frame(resoshvTran_05), file="OsHV_resoshvTran_05_NEW.csv")
write.csv( as.data.frame(resoshvTran_05_Sig), file="OsHV_resoshvTran_05_Sig_NEW.csv")
write.csv( as.data.frame(resoshvTran_05_sig_non_sig), file = "OsHV_resoshvTran_05_sig_non_sig_NEW.csv")
#only do the Uniprot search for OsHV_resoshvTran_05_sig_non_sig_NEW.csv

####Subset files for only those that have Transcript IDs####
OsHV_withID_subset_resoshvTran_05_df <- 
  resoshvTran_05_sig_non_sig[grep("transcript:", rownames(resoshvTran_05_sig_non_sig)), ]
head(OsHV_withID_subset_resoshvTran_05_df)
#split the ID column by :
#install.packages("stringr")
library(stringr)
OsHV_withID_subset_resoshvTran_05_df_split <- str_split_fixed(OsHV_withID_subset_resoshvTran_05_df$ID, ":", 2)
head(OsHV_withID_subset_resoshvTran_05_df_split)
OsHV_withID_subset_resoshvTran_05_df_split <- toString(OsHV_withID_subset_resoshvTran_05_df_split[,2], sep=',')
write(OsHV_withID_subset_resoshvTran_05_df_split, "OsHV_withID_subset_resoshvTran_05_df_split", sep = ",")
#Add ID column to OsHV_withID_subset_resoshvTran_05_df
OsHV_withID_subset_resoshvTran_05_df["ID"] <- OsHV_withID_subset_resoshvTran_05_df_split[,2]
head(OsHV_withID_subset_resoshvTran_05_df)
#perform lookup with Unirpot using the Enemble Genomes Transcript option

####Upload transcripts from the UniProt.ws website ####
oshv1_transcriptIDs_UniProt_SIG <- read.csv("OsHV1_resoshv1Tran_dfSig_transcriptIDstring.csv", header=TRUE)
#make sure to upload as csv version so that the protein names load correctly
head(oshv1_transcriptIDs_UniProt_SIG)
oshv1_transcriptIDs_UniProt_SIG["Challenge"] <- "OsHV1"


####Match lines for only those that have a shared ID in resoshv_non_sig_sig_UNIPROT_NEW.csv
# %in% returns true for every value in the first argument that matches a value in the second argument.
# the order of arguments is important
#load file with the spaces changed to ; so that the column can be separated
resoshv_non_sig_sig_UNIPROT_NEW <- read.csv(file="resoshv_non_sig_sig_UNIPROT_NEW.csv", header=TRUE)

#Merge columns together based on match in the "ID" column 
resoshvTran_05_df_FULL <- merge(OsHV_withID_subset_resoshvTran_05_df, resoshv_non_sig_sig_UNIPROT_NEW[c("ID", "Protein.names","Gene.ontology..GO.","Gene.ontology.IDs")], by="ID")
nrow(resoshvTran_05_df_FULL) #1854
#put LOC info in new column

#### Apoptosis terms in OsHV1 challenge ####
#full list of terms
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

#make sure its a dataframe
resoshvTran_05_df_FULL <- as.data.frame(resoshvTran_05_df_FULL)

#searching for each one individually:
IP3R <- resoshvTran_05_df_FULL[grepl("inositol", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & 
                                             grepl("phosphate", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl("receptor", resoshvTran_05_df_FULL$Protein.names),]
#IP3Ra <- resoshvTran_05_df_FULL[grepl("IP3R", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]  
Rho <- resoshvTran_05_df_FULL[grepl("rho", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                grepl("GTP", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                !grepl("activating",ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                !grepl("mitochondrial", ignore.case=TRUE,resoshvTran_05_df_FULL$Protein.names),]
Rac1 <- resoshvTran_05_df_FULL[grepl("ras", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names) & 
                                             grepl("botulinum", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names) &  grepl("C3", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names),]
cdc42 <- resoshvTran_05_df_FULL[grepl("cdc42", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                              !grepl("activated", ignore.case=TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                              !grepl("RICS",ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]
PKCiota <- resoshvTran_05_df_FULL[grepl("protein",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) & 
                                                grepl("kinase",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                                grepl("EC 2.7.11.13", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names),]
#none
PKCdelta <- resoshvTran_05_df_FULL[grepl("protein", resoshvTran_05_df_FULL$Protein.names) & 
                                                 grepl("kinase",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) & grepl("delta", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) 
                                   ,]

IkB <- resoshvTran_05_df_FULL[grepl("kappa", ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) & 
                                            grepl( "inhibitor", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names) & grepl("nuclear", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names),]

#none are NFkB
NFkB <- resoshvTran_05_df_FULL[grepl("NF-kappa-B", ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) & 
                                            !grepl("inhibitor", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names),]
NFkB <- resoshvTran_05_df_FULL[grepl("nuclear", ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) & 
                                 grepl("kappa", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names),]


GPCRsevenpass <- resoshvTran_05_df_FULL[grepl("G-type", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl("seven-pass", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]
GPCR <- resoshvTran_05_df_FULL[grepl("G-protein" , ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & 
                                             grepl("coupled", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &  grepl("receptor", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]

# 0 hits
Apo3 <- resoshvTran_05_df_FULL[grepl("Apo3",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]
# 0 hits
Apo2 <- resoshvTran_05_df_FULL[grepl("Apo2",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]
# 0 hits
ApoL <- resoshvTran_05_df_FULL[grepl("Apo" , ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) & grepl("ligand", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]
# 0 hits
TWEAK <- resoshvTran_05_df_FULL[grepl("TWEAK",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]
# 0 hits
TRAIL <- resoshvTran_05_df_FULL[grepl("TRAIL",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]
# 0 hits
TRAILb <- resoshvTran_05_df_FULL[grepl("TNF", ignore.case= TRUE,  resoshvTran_05_df_FULL$Protein.names) & 
                                               grepl("apoptosis",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]
#0 hits
TRAILc <- resoshvTran_05_df_FULL[grepl("TNFSF12",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]
# 0 hits
cFLIP <- resoshvTran_05_df_FULL[grepl("c-FLIP",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]
# 0 hits
cFLIP2 <- resoshvTran_05_df_FULL[grepl("FLICE",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]
cFLIP3 <- resoshvTran_05_df_FULL[grepl("CASPER",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]
cFLIP4 <- resoshvTran_05_df_FULL[grepl("FLAME",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]
cFLIP5 <- resoshvTran_05_df_FULL[grepl("CASH",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]
cFLIP6 <- resoshvTran_05_df_FULL[grepl("CLARP",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]
cFLIP7 <- resoshvTran_05_df_FULL[grepl("MRIT",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]
cFLIP8 <- resoshvTran_05_df_FULL[grepl("Fas",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                               grepl("inhibitory",ignore.case= TRUE,  resoshvTran_05_df_FULL$Protein.names),]

RIP <- resoshvTran_05_df_FULL[grepl("receptor",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) & grepl("interacting", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names)
                                          & !grepl("thyroid",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names)
                                          & !grepl("glutamate",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names)
                                          & !grepl("cannabinoid", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names)
                                          & !grepl("AH", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                          !grepl('inositol', ignore.case=TRUE,resoshvTran_05_df_FULL$Protein.names ),]
#non are cAMP specifically
cAMP <-  resoshvTran_05_df_FULL[grepl("cAMP",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]

#no hits
cAMP2 <-  resoshvTran_05_df_FULL[grepl("cyclic",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) 
                                             & grepl("adenosine", ignore.case = TRUE,resoshvTran_05_df_FULL$Protein.names),]

#PKA <- resoshvTran_05_df_FULL[grepl("protein", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & 
                                           # grepl("kinase",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) & grepl("A",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]

PKA <- resoshvTran_05_df_FULL[grepl("cAMP",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                grepl("kinase",resoshvTran_05_df_FULL$Protein.names),]

AC <- resoshvTran_05_df_FULL[grepl("adenylate", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl("cyclase", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                               !grepl("cyclase-associated",ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]

CREB <- resoshvTran_05_df_FULL[grepl("cAMP-responsive", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names),]

#0 hits
DR <- resoshvTran_05_df_FULL[grepl("death",ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names) & grepl("receptor", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]
DR1 <- resoshvTran_05_df_FULL[grepl("TRAMP", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names),]
DR2 <- resoshvTran_05_df_FULL[grepl("WSL",ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names),]
DR3 <-resoshvTran_05_df_FULL[grepl("LARD",ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names),]
DR4 <- resoshvTran_05_df_FULL[grepl("DDR3",ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names),]

Fas <- resoshvTran_05_df_FULL[grepl("fatty", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                            grepl( "acid", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                            grepl( "synthase", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                            !grepl("elongation", ignore.case=TRUE, resoshvTran_05_df_FULL$Protein.names),]

#Fasb <-  resoshvTran_05_df_FULL[grepl("fas", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]
#FasL <- resoshvTran_05_df_FULL[grepl("FasL", ignore.case = TRUE ,resoshvTran_05_df_FULL$Protein.names),]
#FASL2 <- resoshvTran_05_df_FULL[grepl("TNF", ignore.case = TRUE ,resoshvTran_05_df_FULL$Protein.names) &
                                              grepl("superfamily", ignore.case = TRUE ,resoshvTran_05_df_FULL$Protein.names) &
                                              grepl("6", ignore.case = TRUE ,resoshvTran_05_df_FULL$Protein.names),]
#none are what I want
Fas2 <- resoshvTran_05_df_FULL[grepl("Fas", ignore.case = TRUE ,resoshvTran_05_df_FULL$Protein.names),]

#no DR5
DR5 <- resoshvTran_05_df_FULL[grepl("TNFRS10B", ignore.case = TRUE ,resoshvTran_05_df_FULL$Protein.names),]
DR5a <- resoshvTran_05_df_FULL[grepl("TRAIL-R2", ignore.case = TRUE ,resoshvTran_05_df_FULL$Protein.names),]
DR5b <- resoshvTran_05_df_FULL[grepl("TRICK2", ignore.case = TRUE ,resoshvTran_05_df_FULL$Protein.names),]
DR5c <- resoshvTran_05_df_FULL[grepl("KILLER", ignore.case = TRUE ,resoshvTran_05_df_FULL$Protein.names),]
DR5d <- resoshvTran_05_df_FULL[grepl("ZTNFR9", ignore.case = TRUE ,resoshvTran_05_df_FULL$Protein.names),]


TNF <- resoshvTran_05_df_FULL[grepl("tumor", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & 
                                            grepl( "necrosis", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                            grepl( "factor", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & 
                                !grepl("related",resoshvTran_05_df_FULL$Protein.names) &
                                !grepl("ligand", resoshvTran_05_df_FULL$Protein.names) &
                                !grepl("receptor", resoshvTran_05_df_FULL$Protein.names),]
TNFL <- resoshvTran_05_df_FULL[grepl("tumor", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & 
                                grepl( "necrosis", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                grepl( "factor", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & 
                                 !grepl("related",resoshvTran_05_df_FULL$Protein.names) &
                                 grepl("ligand", resoshvTran_05_df_FULL$Protein.names),]
TNFR <- resoshvTran_05_df_FULL[grepl("tumor", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & 
                              grepl( "necrosis", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                              grepl( "factor", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & 
                              !grepl("related",resoshvTran_05_df_FULL$Protein.names) &
                              !grepl("ligand", resoshvTran_05_df_FULL$Protein.names) &
                              grepl("receptor", ignore.case = TRUE,resoshvTran_05_df_FULL$Protein.names ),]
#none are FAIM
FAIM <- resoshvTran_05_df_FULL[grepl("fas", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & 
                                             grepl("apoptotic", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                             grepl("inhibitory", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) ,]

#none are TRAD
TRADD <- resoshvTran_05_df_FULL[grepl("death", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl(  "domain", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                              grepl(  "receptor", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]

#no FADD
FADD <- resoshvTran_05_df_FULL[grepl("Fas", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                             grepl("death", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names),]

caspase <- resoshvTran_05_df_FULL[grepl("caspase", ignore.case = TRUE ,resoshvTran_05_df_FULL$Protein.names) &
                                                !grepl("activity",ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                      !grepl("adapter", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]

#no procapases
procaspase <- resoshvTran_05_df_FULL[grepl("procasp", ignore.case = TRUE ,resoshvTran_05_df_FULL$Protein.names),]

#no sialic acid
siglec <- resoshvTran_05_df_FULL[grepl("siglec", ignore.case = TRUE ,resoshvTran_05_df_FULL$Protein.names),]
sialic <- resoshvTran_05_df_FULL[grepl("sialic", ignore.case = TRUE ,resoshvTran_05_df_FULL$Protein.names),]

#none are IFNLP
IFNLP <- resoshvTran_05_df_FULL[grepl("IFNLP", ignore.case = TRUE ,resoshvTran_05_df_FULL$Protein.names),]
IFNLPa <- resoshvTran_05_df_FULL[grepl("interferon", ignore.case = TRUE ,resoshvTran_05_df_FULL$Protein.names),]

TRAF <- resoshvTran_05_df_FULL[grepl("TNF", ignore.case = TRUE ,resoshvTran_05_df_FULL$Protein.names) &
                                             grepl("receptor-associated", ignore.case = TRUE ,resoshvTran_05_df_FULL$Protein.names) ,]

NuMa <- resoshvTran_05_df_FULL[grepl("mitotic", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl( "apparatus",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]
#NuMAa <- resoshvTran_05_df_FULL[grepl("nuclear", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl( "mitotic", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]
#NuMAb<- resoshvTran_05_df_FULL[grepl("SP-H", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]

#no CAD or ICAD
ICAD <- resoshvTran_05_df_FULL[grepl("DNA", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl("fragmentation", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]
CAD <- resoshvTran_05_df_FULL[grepl("DNA", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl("fragmentation", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                            grepl("beta", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) ,]
CADb <-resoshvTran_05_df_FULL[grepl("caspase", ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) & grepl("activated", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) ,]

#no APIP
APIP <- resoshvTran_05_df_FULL[grepl("apoptosis",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) & grepl( "peptidase", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]
APIPb <- resoshvTran_05_df_FULL[grepl("apoptosis", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl( "protease", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]
APIPc <- resoshvTran_05_df_FULL[grepl("APAF", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names) ,]
APIPc <- resoshvTran_05_df_FULL[grepl("APIP", ignore.case=TRUE, resoshvTran_05_df_FULL$Protein.names) ,]

#no sAC
sAC <- resoshvTran_05_df_FULL[grepl("adenylyl", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl( "cyclase", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]

BI <- resoshvTran_05_df_FULL[grepl("bax",ignore.case=TRUE, resoshvTran_05_df_FULL$Protein.names),]

#no PDRP
PDRP <- resoshvTran_05_df_FULL[grepl("p53", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]

cytoc <- resoshvTran_05_df_FULL[grepl("cytochrome", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl( "c", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                  !grepl("P450", resoshvTran_05_df_FULL$Protein.names) &
                                  !grepl("oxidase", resoshvTran_05_df_FULL$Protein.names ) &
                                  !grepl("reductase", resoshvTran_05_df_FULL$Protein.names ) &
                                  !grepl("lyase", resoshvTran_05_df_FULL$Protein.names ) &
                                  !grepl("b5", resoshvTran_05_df_FULL$Protein.names ) &
                                  !grepl("dehydrogenase", resoshvTran_05_df_FULL$Protein.names ) &
                                  grepl("heme", resoshvTran_05_df_FULL$Protein.names)
                                  ,]

AIF <- resoshvTran_05_df_FULL[grepl("apoptosis", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & 
                                            grepl( "inducing", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                            grepl( "factor", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]
#no Blk
Blk <- resoshvTran_05_df_FULL[grepl("blk", ignore.case=TRUE, resoshvTran_05_df_FULL$Protein.names),]
Blkb <- resoshvTran_05_df_FULL[grepl("B", ignore.case=TRUE, resoshvTran_05_df_FULL$Protein.names) & grepl( "lymphoid", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]
Blkc <- resoshvTran_05_df_FULL[grepl("p55", ignore.case=TRUE, resoshvTran_05_df_FULL$Protein.names),]


myc <- resoshvTran_05_df_FULL[grepl("C-myc", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]
#myc <-  resoshvTran_05_df_FULL[grepl("oncogene", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]

#no aven (it is aka Apoptosis And Caspase Activation Inhibitor)
aven <-  resoshvTran_05_df_FULL[grepl("caspase",ignore.case=TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                               !grepl("activation", ignore.case=TRUE, resoshvTran_05_df_FULL$Protein.names),]

#no p53..
p53 <- resoshvTran_05_df_FULL[grepl("tumor", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl( "antigen",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]
p53 <- resoshvTran_05_df_FULL[grepl("p53", ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) ,]
p53 <- resoshvTran_05_df_FULL[grepl("tumor", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl( "protein", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]

#no smac
smac <- resoshvTran_05_df_FULL[grepl("low", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl( "PI", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]
smac <- resoshvTran_05_df_FULL[grepl("mitochondrial", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl( "activator", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]
smac1 <- resoshvTran_05_df_FULL[grepl("IAP", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names) & grepl( "direct", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]

bcl <- resoshvTran_05_df_FULL[grepl("bcl",ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                grepl("apoptosis", ignore.case=TRUE,resoshvTran_05_df_FULL$Protein.names),]
bclw <-  resoshvTran_05_df_FULL[grepl("bcl-2",ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                  grepl("like",ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names),]
#no mcl
mcl <- resoshvTran_05_df_FULL[grepl("myeloid",ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                            grepl("leukemia", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names),]
mcl <- resoshvTran_05_df_FULL[grepl("cell", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl( "differentiation", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]

#no A1
A1 <- resoshvTran_05_df_FULL[grepl("bcl-2",ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names),]

#no Bak
Bak <- resoshvTran_05_df_FULL[grepl("antagonist", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl( "killer", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]
Bak2 <- resoshvTran_05_df_FULL[grepl("bcl2l7",ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names),]

#none are Bcl-2 associated protein 
BAX <- resoshvTran_05_df_FULL[grepl("bcl",ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names),]

#no Bok
Bok <-  resoshvTran_05_df_FULL[grepl("ovarian",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) & grepl( "killer", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]
Bok <- resoshvTran_05_df_FULL[grepl("ovarian",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]

#no BAD
BAD <- resoshvTran_05_df_FULL[grepl("death",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) & grepl( "promoter",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]
BAD2 <- resoshvTran_05_df_FULL[grepl("binding", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl( "component", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]
BAD3 <- resoshvTran_05_df_FULL[grepl("BAD", ignore.case=TRUE, resoshvTran_05_df_FULL$Protein.names),]

BAG <- resoshvTran_05_df_FULL[grepl("BAG",ignore.case=TRUE, resoshvTran_05_df_FULL$Protein.names),]

#no BH3
BH3 <- resoshvTran_05_df_FULL[grepl("BH3",ignore.case=TRUE, resoshvTran_05_df_FULL$Protein.names),]
Bid <- resoshvTran_05_df_FULL[grepl("Bid",ignore.case=TRUE, resoshvTran_05_df_FULL$Protein.names),]

#no BIM or Bcl-2-like protein 11
BIM <- resoshvTran_05_df_FULL[grepl("BIM",ignore.case=TRUE, resoshvTran_05_df_FULL$Protein.names),]
BIMb <- resoshvTran_05_df_FULL[grepl("Bcl-2",ignore.case=TRUE, resoshvTran_05_df_FULL$Protein.names),]

# no Bik
Bik <- resoshvTran_05_df_FULL[grepl("interacting", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl( "killer",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]

#no Puma
Puma <- resoshvTran_05_df_FULL[grepl("modulator",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) & grepl( "apoptosis",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]
Puma <- resoshvTran_05_df_FULL[grepl("JFY", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]

#no ceramide
ceramide <- resoshvTran_05_df_FULL[grepl("ceramide",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) ,]

MEKK_MAPK <- resoshvTran_05_df_FULL[grepl("mitogen-activated", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                                  grepl( "kinase", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                                  !grepl("kinase-binding", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]

#no p38 MAPK
p38_MAPK <-resoshvTran_05_df_FULL[grepl("mitogen-activated",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                                 grepl( "kinase",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                                 grepl("p38", ignore.case=TRUE, resoshvTran_05_df_FULL$Protein.names),]

#none
endoG <- resoshvTran_05_df_FULL[grepl("endonuclease",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) & grepl( "G",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]

PI3 <- resoshvTran_05_df_FULL[grepl("phosphatidylinositol", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl( "3-kinase",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]

#no AKT
AKTB <- resoshvTran_05_df_FULL[grepl("protein", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                             grepl("kinase",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                             grepl("B",ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]
#no protein kinase B!
akt <- resoshvTran_05_df_FULL[grepl("akt", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names),]

#no JNK
JNK <- resoshvTran_05_df_FULL[grepl("stress-activated", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl( "JNK", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]


#none
p35 <- resoshvTran_05_df_FULL[grepl("cyclin", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl( "dependent",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                            grepl( "5",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                            grepl("activator",ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]
p35 <- resoshvTran_05_df_FULL[grepl("p35", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]
                              

#cJUN <- resoshvTran_05_df_FULL[grepl("jun",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) ,]
#cJUN <- resoshvTran_05_df_FULL[grepl("c-jun",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]


cJUNa <- resoshvTran_05_df_FULL[grepl("AP-1", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                              grepl("transcription",ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]
#no CD151
CD151 <- resoshvTran_05_df_FULL[grepl("CD151", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names),]

#none
BTG1 <- resoshvTran_05_df_FULL[grepl("B",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) & grepl("translocation", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]

#1 
BTG2 <- resoshvTran_05_df_FULL[grepl("BTG", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) ,]

IAP <- resoshvTran_05_df_FULL[grepl("inhibitor", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                            grepl("apoptosis", ignore.case=TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                            !grepl("death-associated", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                            !grepl("caspase", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                            !grepl("5-like", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                            !grepl("TP53",resoshvTran_05_df_FULL$Protein.names),]

IAP2 <- resoshvTran_05_df_FULL[grepl("baculoviral", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                             grepl("IAP", ignore.case = TRUE,resoshvTran_05_df_FULL$Protein.names),]

#none
API <- resoshvTran_05_df_FULL[grepl("inhibitor",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                            grepl("apoptosis", ignore.case=TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                            !grepl("death", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                            !grepl("caspase",ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                            grepl("5", ignore.case = TRUE,resoshvTran_05_df_FULL$Protein.names),]
#none
DIAP <- resoshvTran_05_df_FULL[grepl("inhibitor", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                             grepl("apoptosis", ignore.case=TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                             grepl("death", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names),]
#none
CAAP <- resoshvTran_05_df_FULL[grepl("inhibitor", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                             grepl("apoptosis", ignore.case=TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                             !grepl("death", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                             grepl("caspase", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]

PCD <- resoshvTran_05_df_FULL[grepl("programmed",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) & 
                                            grepl("death",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                            grepl("cell", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                            !grepl("interacting",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]
#none
NR13 <- resoshvTran_05_df_FULL[grepl("NR13", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]


PAWR <- resoshvTran_05_df_FULL[grepl("WT1", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl("regulator", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]

IFI44 <- resoshvTran_05_df_FULL[grepl("interferon", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl("44",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]

#no IFI6
IFI6 <- resoshvTran_05_df_FULL[grepl("interferon",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) & grepl("6",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]

TLR <- resoshvTran_05_df_FULL[grepl("toll",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) & grepl("receptor", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]

MyD88 <-resoshvTran_05_df_FULL[grepl("myeloid", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl("differentiation",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]


CCAR <- resoshvTran_05_df_FULL[grepl("cell", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & grepl("division", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                             grepl("apoptosis",ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]

#none
ASPP <- resoshvTran_05_df_FULL[grepl("apoptosis-stimulating",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) & grepl("p53",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names),]


#GIMAP <- resoshvTran_05_df_FULL[grepl("GIMAP",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) ,]
GIMAP <- resoshvTran_05_df_FULL[grepl("GTPase", ignore.case = TRUE, resoshvTran_05_df_FULL$Protein.names) &
                                              grepl("IMAP",ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]


#Combine all lists from above into one new big table
#only include those that have hits and are confirmed to be correct
resoshvTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED <- rbind(IP3R, Rho, Rac1, cdc42, IkB, GPCRsevenpass, GPCR, RIP, PKA, AC, CREB, Fas,
                                                                         TNF, TNFL,TNFR, caspase, TRAF, NuMa, BI, cytoc, AIF, myc, bcl, bclw, 
                                                                         BAG, MEKK_MAPK, PI3, cJUNa, BTG2, IAP, IAP2, PCD, PAWR, IFI44, TLR, MyD88,
                                                                         CCAR, GIMAP, PKCiota)

#subset only those that are significant
resoshvTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_significant <- resoshvTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED %>%
  filter(Significance == "significant")
#add column with apoptosis pathwway and label genes with same name with extra (#)
write.csv(resoshvTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_significant, file="resoshvTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_significant.csv")
resoshvTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_significant_path <- read.csv(file="resoshvTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_significant.csv", header=TRUE)

#only 1 was significant without the P-value correction and p <0.05...the IP3R receptor

####Graphing OsHV1 data ####

#How many have a p <0.1 ?
resoshvTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_01 <- resoshvTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED %>% filter(padj < 0.1)
#add pathway type to this
resoshvTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_01["Type"] <- "extrinsic" 
resoshvTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_01[4, 12] = "intrinsic"


#labels for signifiance with *= p <0.05
label.oshv.sig.df <- data.frame(Protein.names= c("Inositol 1,4,5-trisphosphate receptor type 1"), log2FoldChange=c(6.5))
oshv_sig_1 <- ggplot(resoshvTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_01) + 
  geom_col(aes(x=Protein.names, y=log2FoldChange, fill=Type)) + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + geom_hline(yintercept=0.0) +
  scale_y_continuous(name ="log2 Fold Change", breaks = scales::pretty_breaks(n = 20)) + 
  ggtitle("Apoptosis Genes in M. gigas after OsHV-1 Challenge") + 
  theme(plot.title = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5))

Oshv_sig_2 <- oshv_sig_1 +
  geom_text(data=label.oshv.sig.df, aes(x=Protein.names, y=log2FoldChange), label = c("*")) +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  scale_x_discrete(name="Transcript Annotation")
                                                          
####GIMAP IAP graphing for OsHV-1 ####
# graphing all of them regardless of significance
oshv_GIMAP_IAP <- rbind(GIMAP, IAP, IAP2)
ncol(oshv_GIMAP_IAP) #11
write.csv(oshv_GIMAP_IAP, file="oshv_GIMAP_IAP.csv") # add Type column
oshv_GIMAP_IAP_updated <- read.csv(file="oshv_GIMAP_IAP.csv",header =TRUE)

oshv_GIMAP_IAP_1 <- ggplot(oshv_GIMAP_IAP_updated) + 
  geom_col(aes(x=Protein.names, y=log2FoldChange, fill=Type)) + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + geom_hline(yintercept=0.0) +
  scale_y_continuous(name ="log2 Fold Change", breaks = scales::pretty_breaks(n = 20)) + 
  ggtitle("GIMAP and IAP Apoptosis Genes in M. gigas after OsHV-1 Challenge") + 
  theme(plot.title = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5))

Oshv_GIMAP_IAP_2 <- oshv_GIMAP_IAP_1  +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  scale_x_discrete(name="Transcript Annotation", limits=c("Apoptosis 1 inhibitor (1)",
                                                          "Apoptosis 1 inhibitor (2)",
                                                          "Apoptosis 1 inhibitor (3)",
                                                          "Apoptosis 2 inhibitor (1)",
                                                          "Apoptosis 2 inhibitor (2)",
                                                          "Baculoviral IAP repeat-containing protein 2",
                                                          "Baculoviral IAP repeat-containing protein 3",
                                                          "Baculoviral IAP repeat-containing protein 3 (Fragment)",
                                                          "Baculoviral IAP repeat-containing protein 6",
                                                          "Baculoviral IAP repeat-containing protein 7",
                                                          "Baculoviral IAP repeat-containing protein 7-A",
                                                          "Baculoviral IAP repeat-containing protein 7-B (1)",
                                                          "Baculoviral IAP repeat-containing protein 7-B (20",
                                                          "Baculoviral IAP repeat-containing protein 7-B (3)",
                                                          "GTPase IMAP family member 1",
                                                          "GTPase IMAP family member 4 (1)",
                                                          "GTPase IMAP family member 4 (2)",
                                                          "GTPase IMAP family member 4 (3)",
                                                          "GTPase IMAP family member 4 (4)",
                                                          "GTPase IMAP family member 4 (5)",
                                                          "GTPase IMAP family member 4 (6)",
                                                          "GTPase IMAP family member 4 (7)",
                                                          "GTPase IMAP family member 4 (8)",
                                                          "GTPase IMAP family member 4 (9)",
                                                          "GTPase IMAP family member 4 (10)",
                                                          "GTPase IMAP family member 4 (11)",
                                                          "GTPase IMAP family member 7 (12)",
                                                          "GTPase IMAP family member 7 (13)",
                                                          "Inhibitor of apoptosis protein (1)",
                                                          "Putative apoptosis inhibitor ORF42",
                                                          "Putative inhibitor of apoptosis (2)",
                                                          "Putative inhibitor of apoptosis (3)",
                                                          "Putative inhibitor of apoptosis (4)"), labels= c("Apoptosis 1 inhibitor (1)"= "IAP1",
                                                                                                            "Apoptosis 1 inhibitor (2)"="IAP1",
                                                                                                            "Apoptosis 1 inhibitor (3)"="IAP1",
                                                                                                            "Apoptosis 2 inhibitor (1)"="IAP2",
                                                                                                            "Apoptosis 2 inhibitor (2)"="IAP2",
                                                                                                            "Baculoviral IAP repeat-containing protein 2"="BIR IAP repeat-containing 2",
                                                                                                            "Baculoviral IAP repeat-containing protein 3"="BIR IAP repeat-containing 3",
                                                                                                            "Baculoviral IAP repeat-containing protein 3 (Fragment)"="BIR IAP repeat-containing 3 (fragment)",
                                                                                                            "Baculoviral IAP repeat-containing protein 6"="BIR IAP repeat-containing 6",
                                                                                                            "Baculoviral IAP repeat-containing protein 7"="BIR IAP repeat-containing 7",
                                                                                                            "Baculoviral IAP repeat-containing protein 7-A"="BIR IAP repeat-containing 7-A",
                                                                                                            "Baculoviral IAP repeat-containing protein 7-B (1)"="BIR IAP repeat-containing 7-B",
                                                                                                            "Baculoviral IAP repeat-containing protein 7-B (20"="BIR IAP repeat-containing 7-B",
                                                                                                            "Baculoviral IAP repeat-containing protein 7-B (3)"="BIR IAP repeat-containing 7-B",
                                                                                                            "GTPase IMAP family member 1"="GIMAP1",
                                                                                                            "GTPase IMAP family member 4 (1)"="GIMAP4",
                                                                                                            "GTPase IMAP family member 4 (2)"="GIMAP4",
                                                                                                            "GTPase IMAP family member 4 (3)"="GIMAP4",
                                                                                                            "GTPase IMAP family member 4 (4)"="GIMAP4",
                                                                                                            "GTPase IMAP family member 4 (5)"="GIMAP4",
                                                                                                            "GTPase IMAP family member 4 (6)"="GIMAP4",
                                                                                                            "GTPase IMAP family member 4 (7)"="GIMAP4",
                                                                                                            "GTPase IMAP family member 4 (8)"="GIMAP4",
                                                                                                            "GTPase IMAP family member 4 (9)"="GIMAP4",
                                                                                                            "GTPase IMAP family member 4 (10)"="GIMAP4",
                                                                                                            "GTPase IMAP family member 4 (11)"="GIMAP4",
                                                                                                            "GTPase IMAP family member 7 (12)"="GIMAP7",
                                                                                                            "GTPase IMAP family member 7 (13)"="GIMAP7",
                                                                                                            "Inhibitor of apoptosis protein (1)"= "IAP",
                                                                                                            "Putative apoptosis inhibitor ORF42"= "putative IAP ORF42",
                                                                                                            "Putative inhibitor of apoptosis (2)"="putative IAP",
                                                                                                            "Putative inhibitor of apoptosis (3)"="putative IAP",
                                                                                                            "Putative inhibitor of apoptosis (4)"="putative IAP"))


#### Vibrio PROBIOTIC CHALLENGED Differential Gene Expression Analysis ####
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
#here I only care about treatment versus control
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
resBacTran <- results(ddsBacTran)

#summary is just printing a table for you, you need to tell it what threshold you want
help("summary",package="DESeq2")
alpha <- 0.05 #set alpha to 0.05, this will control FDR
summary(resBacTran) #default FDR is still 0.1
summary(resBacTran, alpha) #now showing all genes with FRD < 0.05

#To get the significant genes
#The independent filtering in results() has an argument 'alpha'
#which is used to optimize a cutoff on mean normalized count
#to maximize the number of genes with padj < alpha
resBacTran_05 <- results(ddsBacTran, alpha= alpha) #set FDR to 0.05 now
resBacTran_05_Sig <- resBacTran[which(resBacTran$padj < alpha),]
summary(resBacTran_05) #this is all the genes
summary(resBacTran_05_Sig) #this is the significant ones!
sum(resBacTran_05$padj < 0.05, na.rm=TRUE) #1437 tells you how many genes have expected FDR ≤ 0.05
sum(resBacTran_05_Sig$padj < 0.05, na.rm=TRUE) #1428
resBacTran_05_Sig$Significance <- sig
resBacTran_05_nonSig <- resBacTran[which(resBacTran$padj > alpha),] #create list of nonsig
nonsig <- "non-significant"
resBacTran_05_nonSig$Significance <- nonsig
head(resBacTran_05_Sig)
head(resBacTran_05_nonSig)

#add ID column with the rownames so a merge can happen later
resBacTran_05_Sig["ID"] <- rownames(resBacTran_05_Sig) #add a new column with the rownames for match
resBacTran_05_nonSig["ID"] <- rownames(resBacTran_05_nonSig)
resBacTran_05["ID"] <- rownames(resBacTran_05)
resBacTran_05_sig_non_sig <- rbind(resBacTran_05_Sig, resBacTran_05_nonSig)

#metadata on meaning of the columns
mcols(resBacTran_05_Sig, use.names = TRUE)
mcols(resBacTran_05_Sig)$description
#shows treatment vs. control with baseMean as the mean of normalized counts for all samples

#Order by Log2FC
head( resBacTran_05_sig_non_sig[ order( resBacTran_05_sig_non_sig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resBacTran_05_sig_non_sig[ order( resBacTran_05_sig_non_sig$log2FoldChange ), ] ) #tail for strongest up regulation

####p-value correction for both sets of results ####
#First Visualize histograms
#histogram of P- values to visualize any "hills" or "U shape"
# hill means variance of the null distribution too high, U shape means variance assumed too low
hist(resBacTran_05$pvalue, breaks= 20, col = "grey")
hist(resBacTran_05_Sig$pvalue, breaks = 20, col = "grey") #hill

#looks okay, not doing the correction because we want to avoid extra false positives
#Michael Love says he usually doesn't have a need to do one

#Visualize Results with Diagnostic Plots#
#MA plot, useful overview for experiment with two-group comparison. Plots log2FC over mean of normalized counts
#genes with adjusted p value 
plotMA(resBacTran_05)
plotMA(resBacTran_05_Sig)

#Export Results to CSV
write.csv( as.data.frame(resBacTran_05), file="Bac_resBacTran_05_NEW.csv")
write.csv( as.data.frame(resBacTran_05_Sig), file="Bac_resBacTran_05_Sig_NEW.csv")
write.csv( as.data.frame(resBacTran_05_sig_non_sig), file = "Bac_resBacTran_05_sig_non_sig_NEW.csv")
#only do the Uniprot search for OsHV_resoshvTran_05_sig_non_sig_NEW.csv

####Subset files for only those that have Transcript IDs####
Bac_withID_subset_resoshvTran_05_df <- 
  resBacTran_05_sig_non_sig[grep("transcript:", rownames(resBacTran_05_sig_non_sig)), ]
head(Bac_withID_subset_resoshvTran_05_df)
#split the ID column by :
Bac_withID_subset_resoshvTran_05_df_split <- str_split_fixed(Bac_withID_subset_resoshvTran_05_df$ID, ":", 2)
head(Bac_withID_subset_resoshvTran_05_df_split)
Bac_withID_subset_resoshvTran_05_df_split_2 <- toString(Bac_withID_subset_resoshvTran_05_df_split[,2], sep=',')
write(Bac_withID_subset_resoshvTran_05_df_split_2, "Bac_withID_subset_resoshvTran_05_df_split", sep = ",")
#Add ID column to OsHV_withID_subset_resoshvTran_05_df
Bac_withID_subset_resoshvTran_05_df["ID"] <- Bac_withID_subset_resoshvTran_05_df_split[,2]
head(Bac_withID_subset_resoshvTran_05_df)
#perform lookup with Uniprot using the Enemble Genomes Transcript option

####Upload transcripts from the UniProt.ws website ####
Bac_transcriptIDs_UniProt_sig_nonsig <- read.csv("Bac_Tran_sig_non_sig_Uniprot.csv", header=TRUE)
#make sure to upload as csv version so that the protein names load correctly
head(Bac_transcriptIDs_UniProt_sig_nonsig)
Bac_transcriptIDs_UniProt_sig_nonsig["Challenge"] <- "Bac"
head(Bac_transcriptIDs_UniProt_sig_nonsig)

####Match lines for only those that have a shared ID in Bac_Tran_sig_non_sig_Uniprot.csv
# %in% returns true for every value in the first argument that matches a value in the second argument.
# the order of arguments is important

#Merge columns together based on match in the "ID" column 
resBacTran_05_df_FULL <- merge(Bac_withID_subset_resoshvTran_05_df, Bac_transcriptIDs_UniProt_sig_nonsig[c("ID", "Protein.names","Gene.ontology..GO.","Gene.ontology.IDs", "Challenge")], by="ID")
nrow(resBacTran_05_df_FULL) #8333

#### Apoptosis terms in OsHV1 challenge ####
#full list of terms
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

#make sure its a dataframe
resBacTran_05_df_FULL <- as.data.frame(resBacTran_05_df_FULL)

#searching for each one individually:
oshv_IP3R <- resBacTran_05_df_FULL[grepl("inositol", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & 
                                 grepl("phosphate", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & 
                                grepl("receptor", resBacTran_05_df_FULL$Protein.names) &
                                !grepl("interacting", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names) &
                                !grepl("protein", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names),]
#IP3Ra <-resBacTran_05_df_FULL[grepl("IP3R", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]  
oshv_Rho <- resBacTran_05_df_FULL[grepl("rho", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names) &
                                grepl("GTP", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                !grepl("activating",ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                !grepl("mitochondrial", ignore.case=TRUE,resBacTran_05_df_FULL$Protein.names),]
oshv_Rac1 <- resBacTran_05_df_FULL[grepl("ras", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names) & 
                                 grepl("botulinum", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names) &  grepl("C3", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names),]
oshv_cdc42 <- resBacTran_05_df_FULL[grepl("cdc42", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names) &
                                  !grepl("activated", ignore.case=TRUE,resBacTran_05_df_FULL$Protein.names) &
                                  !grepl("RICS",ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]
oshv_PKC <- resBacTran_05_df_FULL[grepl("protein",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) & 
                                    grepl("kinase",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) &
                                    grepl('C', ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names) &
                                  grepl("EC 2.7.11.13", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names),]
oshv_PKCdelta <- resBacTran_05_df_FULL[grepl("protein", resBacTran_05_df_FULL$Protein.names) & 
                                     grepl("kinase",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) & grepl("delta", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]
                               

oshv_IkB <- resBacTran_05_df_FULL[grepl("kappa", ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) & 
                                grepl( "inhibitor", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names) & grepl("nuclear", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names),]

#none are NFkB
oshv_NFkB <- resBacTran_05_df_FULL[grepl("NF-kappa-B", ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) & 
                                 !grepl("inhibitor", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names),]
oshv_NFkB <- resBacTran_05_df_FULL[grepl("nuclear", ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) & 
                                 grepl("kappa", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names),]


oshv_GPCRsevenpass <- resBacTran_05_df_FULL[grepl("G-type", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl("seven-pass", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]
oshv_GPCR <- resBacTran_05_df_FULL[grepl("G-protein" , ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & 
                                 grepl("coupled", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &  grepl("receptor", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]

# 0 hits
oshv_Apo3 <- resBacTran_05_df_FULL[grepl("Apo3",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]
# 0 hits
oshv_Apo2 <- resBacTran_05_df_FULL[grepl("Apo2",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]
# 0 hits
oshv_ApoL <- resBacTran_05_df_FULL[grepl("Apo" , ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) & grepl("ligand", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]
# 0 hits
oshv_TWEAK <- resBacTran_05_df_FULL[grepl("TWEAK",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]
# 0 hits
oshv_TRAIL <- resBacTran_05_df_FULL[grepl("TRAIL",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]
# 0 hits
oshv_TRAILb <- resBacTran_05_df_FULL[grepl("Tumor", ignore.case= TRUE,  resBacTran_05_df_FULL$Protein.names) & 
                                   grepl("apoptosis",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]
#0 hits
oshv_TRAILc <- resBacTran_05_df_FULL[grepl("TNFSF12",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]
# 0 hits
oshv_cFLIP <- resBacTran_05_df_FULL[grepl("FLIP",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]
# 0 hits
oshv_cFLIP2 <- resBacTran_05_df_FULL[grepl("FLICE",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]
oshv_cFLIP3 <- resBacTran_05_df_FULL[grepl("CASPER",ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]
oshv_cFLIP4 <- resBacTran_05_df_FULL[grepl("FLAME",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]
oshv_cFLIP5 <- resBacTran_05_df_FULL[grepl("CASH",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]
oshv_cFLIP6 <- resBacTran_05_df_FULL[grepl("CLARP",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]
oshv_cFLIP7 <- resBacTran_05_df_FULL[grepl("MRIT",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]
oshv_cFLIP8 <- resBacTran_05_df_FULL[grepl("Fas",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) &
                                   grepl("inhibitory",ignore.case= TRUE,  resBacTran_05_df_FULL$Protein.names),]

oshv_RIPK <- resBacTran_05_df_FULL[grepl("receptor",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) & grepl("interacting", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names)
                              & !grepl("thyroid",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names)
                              & !grepl("glutamate",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names)
                              & !grepl("cannabinoid", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names)
                              & !grepl("AH", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                !grepl('inositol', ignore.case=TRUE,resBacTran_05_df_FULL$Protein.names) &
                                !grepl("nogo", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names),]
#non are cAMP specifically
oshv_cAMP <-  resBacTran_05_df_FULL[grepl("cAMP",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]

#no hits
oshv_cAMP2 <-  resBacTran_05_df_FULL[grepl("cyclic",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) 
                                 & grepl("adenosine", ignore.case = TRUE,resBacTran_05_df_FULL$Protein.names),]

#oshv_PKA <- resBacTran_05_df_FULL[grepl("protein", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & 
# grepl("kinase",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) & grepl("A",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]

oshv_PKA <- resBacTran_05_df_FULL[grepl("cAMP",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) &
                                grepl("kinase",resBacTran_05_df_FULL$Protein.names),]

oshv_AC <- resBacTran_05_df_FULL[grepl("adenylate", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl("cyclase", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                               !grepl("cyclase-associated",ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]

oshv_CREB <- resBacTran_05_df_FULL[grepl("cAMP-responsive", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names) &
                                     !grepl("modulator", resBacTran_05_df_FULL$Protein.names),]

#0 hits
oshv_DR <- resBacTran_05_df_FULL[grepl("death",ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names) & grepl("receptor", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]
oshv_DR1 <- resBacTran_05_df_FULL[grepl("TRAMP", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names),]
oshv_DR2 <- resBacTran_05_df_FULL[grepl("WSL",ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names),]
oshv_DR3 <-resBacTran_05_df_FULL[grepl("LARD",ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names),]
oshv_DR4 <- resBacTran_05_df_FULL[grepl("DDR3",ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names),]

oshv_Fas <- resBacTran_05_df_FULL[grepl("fatty", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                grepl( "acid", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                grepl( "synthase", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                !grepl("elongation", ignore.case=TRUE, resBacTran_05_df_FULL$Protein.names),]

#none are what I want
oshv_Fasb <-  resBacTran_05_df_FULL[grepl("fas", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]
oshv_FasL <- resBacTran_05_df_FULL[grepl("FasL", ignore.case = TRUE ,resBacTran_05_df_FULL$Protein.names),]
oshv_FASL2 <- resBacTran_05_df_FULL[grepl("TNF", ignore.case = TRUE ,resBacTran_05_df_FULL$Protein.names) &
  grepl("superfamily", ignore.case = TRUE ,resBacTran_05_df_FULL$Protein.names) &
  grepl("6", ignore.case = TRUE ,resBacTran_05_df_FULL$Protein.names),]
oshv_Fas2 <- resBacTran_05_df_FULL[grepl("Fas", ignore.case = TRUE ,resBacTran_05_df_FULL$Protein.names),]

#no DR5
oshv_DR5 <-resBacTran_05_df_FULL[grepl("TNFRS10B", ignore.case = TRUE ,resBacTran_05_df_FULL$Protein.names),]
oshv_DR5a <- resBacTran_05_df_FULL[grepl("TRAIL-R2", ignore.case = TRUE ,resBacTran_05_df_FULL$Protein.names),]
oshv_DR5b <- resBacTran_05_df_FULL[grepl("TRICK2", ignore.case = TRUE ,resBacTran_05_df_FULL$Protein.names),]
oshv_DR5c <- resBacTran_05_df_FULL[grepl("KILLER", ignore.case = TRUE ,resBacTran_05_df_FULL$Protein.names),]
oshv_DR5d <- resBacTran_05_df_FULL[grepl("ZTNFR9", ignore.case = TRUE ,resBacTran_05_df_FULL$Protein.names),]


oshv_TNF <- resBacTran_05_df_FULL[grepl("tumor", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & 
                                grepl( "necrosis", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                grepl( "factor", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & 
                                !grepl("related",resBacTran_05_df_FULL$Protein.names) &
                                !grepl("ligand", resBacTran_05_df_FULL$Protein.names) &
                                !grepl("receptor", resBacTran_05_df_FULL$Protein.names),]
oshv_TNFL <- resBacTran_05_df_FULL[grepl("tumor", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & 
                                 grepl( "necrosis", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                 grepl( "factor", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & 
                                 !grepl("related",resBacTran_05_df_FULL$Protein.names) &
                                 grepl("ligand", resBacTran_05_df_FULL$Protein.names),]
oshv_TNFR <- resBacTran_05_df_FULL[grepl("tumor", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & 
                                 grepl( "necrosis", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                 grepl( "factor", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & 
                                 !grepl("related",resBacTran_05_df_FULL$Protein.names) &
                                 !grepl("ligand", resBacTran_05_df_FULL$Protein.names) &
                                 grepl("receptor", ignore.case = TRUE,resBacTran_05_df_FULL$Protein.names ),]
#none are FAIM
oshv_FAIM <- resBacTran_05_df_FULL[grepl("fas", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & 
                                 grepl("inhibitory", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) ,]

#none are TRAD
oshv_TRADD <- resBacTran_05_df_FULL[grepl("death", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl(  "domain", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                  grepl(  "receptor", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]

#no FADD
oshv_FADD <- resBacTran_05_df_FULL[grepl("Fas", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names) &
                                 grepl("death", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names),]

oshv_caspase <- resBacTran_05_df_FULL[grepl("caspase", ignore.case = TRUE ,resBacTran_05_df_FULL$Protein.names) &
                                    !grepl("activity",ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                    !grepl("adapter", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]

#no procapases
oshv_procaspase <- resBacTran_05_df_FULL[grepl("procasp", ignore.case = TRUE ,resBacTran_05_df_FULL$Protein.names),]

#no sialic acid
oshv_siglec <- resBacTran_05_df_FULL[grepl("siglec", ignore.case = TRUE ,resBacTran_05_df_FULL$Protein.names),]
oshv_sialic <- resBacTran_05_df_FULL[grepl("sialic", ignore.case = TRUE ,resBacTran_05_df_FULL$Protein.names),]

#none are IFNLP
oshv_IFNLP <- resBacTran_05_df_FULL[grepl("IFNLP", ignore.case = TRUE ,resBacTran_05_df_FULL$Protein.names),]
oshv_IFNLPa <- resBacTran_05_df_FULL[grepl("interferon", ignore.case = TRUE ,resBacTran_05_df_FULL$Protein.names),]

oshv_TRAF <- resBacTran_05_df_FULL[grepl("TNF", ignore.case = TRUE ,resBacTran_05_df_FULL$Protein.names) &
                                 grepl("receptor-associated", ignore.case = TRUE ,resBacTran_05_df_FULL$Protein.names) ,]

oshv_NuMa <- resBacTran_05_df_FULL[grepl("mitotic", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl( "apparatus",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) &
                                        !grepl("p62", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names),]
#oshv_NuMAa <- resoshvTran_05_df_FULL[grepl("nuclear", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl( "mitotic", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]
#oshv_NuMAb<- resoshvTran_05_df_FULL[grepl("SP-H", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]

#no CAD or ICAD
oshv_ICAD <- resBacTran_05_df_FULL[grepl("DNA", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl("fragmentation", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]
oshv_CAD <- resBacTran_05_df_FULL[grepl("DNA", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl("fragmentation", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                grepl("beta", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) ,]
oshv_CADb <-resoshvTran_05_df_FULL[grepl("caspase", ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) & grepl("activated", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) ,]

#no APIP
oshv_APIP <- resBacTran_05_df_FULL[grepl("apoptosis",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) & grepl( "peptidase", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]
oshv_APIPb <- resBacTran_05_df_FULL[grepl("apoptosis", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl( "protease", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]
oshv_APIPc <- resBacTran_05_df_FULL[grepl("APAF", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names) ,]
oshv_APIPc <- resBacTran_05_df_FULL[grepl("APIP", ignore.case=TRUE, resBacTran_05_df_FULL$Protein.names) ,]

#no sAC
oshv_sAC <- resBacTran_05_df_FULL[grepl("adenylyl", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl( "soluble", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]

oshv_BI <- resBacTran_05_df_FULL[grepl("bax",ignore.case=TRUE, resBacTran_05_df_FULL$Protein.names),]

#no PDRP
oshv_PDRP <-resBacTran_05_df_FULL[grepl("p53", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]

oshv_cytoc <- resBacTran_05_df_FULL[grepl("cytochrome", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl( "c", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                  !grepl("P450", resBacTran_05_df_FULL$Protein.names) &
                                  !grepl("oxidase", resBacTran_05_df_FULL$Protein.names ) &
                                  !grepl("reductase", resBacTran_05_df_FULL$Protein.names ) &
                                  !grepl("lyase", resBacTran_05_df_FULL$Protein.names ) &
                                  !grepl("b5", resBacTran_05_df_FULL$Protein.names ) &
                                  !grepl("dehydrogenase", resBacTran_05_df_FULL$Protein.names ) &
                                  grepl("heme", resBacTran_05_df_FULL$Protein.names)
                                ,]

oshv_AIF <- resBacTran_05_df_FULL[grepl("apoptosis", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & 
                                grepl( "inducing", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                grepl( "factor", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]
#no Blk
oshv_Blk <- resBacTran_05_df_FULL[grepl("bik", ignore.case=TRUE, resBacTran_05_df_FULL$Protein.names),]
oshv_Blkb <- resBacTran_05_df_FULL[grepl("B", ignore.case=TRUE, resBacTran_05_df_FULL$Protein.names) & grepl( "lymphoid", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]
oshv_Blkc <- resBacTran_05_df_FULL[grepl("p55", ignore.case=TRUE, resBacTran_05_df_FULL$Protein.names),]


oshv_myc <- resBacTran_05_df_FULL[grepl("C-myc", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]
#oshv_myc <-  resBacTran_05_df_FULL[grepl("oncogene", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]

#no aven (it is aka Apoptosis And Caspase Activation Inhibitor)
oshv_aven <-  resBacTran_05_df_FULL[grepl("caspase",ignore.case=TRUE, resBacTran_05_df_FULL$Protein.names) &
                                  !grepl("activation", ignore.case=TRUE, resBacTran_05_df_FULL$Protein.names),]

#no p53..
oshv_p53 <- resBacTran_05_df_FULL[grepl("tumor", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl( "antigen",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]
oshv_p53 <- resBacTran_05_df_FULL[grepl("p53", ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) ,]
oshv_p53 <- resBacTran_05_df_FULL[grepl("tumor", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl( "protein", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]

#no smac
oshv_smac <- resBacTran_05_df_FULL[grepl("low", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl( "PI", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]
oshv_smac <- resBacTran_05_df_FULL[grepl("mitochondrial", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl( "activator", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]
oshv_smac1 <- resBacTran_05_df_FULL[grepl("IAP", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names) & grepl( "direct", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]

oshv_bcl <- resBacTran_05_df_FULL[grepl("bcl",ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names) &
                                grepl("apoptosis", ignore.case=TRUE,resBacTran_05_df_FULL$Protein.names),]
oshv_bclw <-  resBacTran_05_df_FULL[grepl("bcl-2",ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names) &
                                  grepl("like",ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names),]
#no mcl
oshv_mcl <- resBacTran_05_df_FULL[grepl("myeloid",ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names) &
                                grepl("leukemia", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names),]
oshv_mcl <- resBacTran_05_df_FULL[grepl("cell", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl( "differentiation", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]

#no A1
oshv_A1 <- resBacTran_05_df_FULL[grepl("bcl-2",ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names),]

#no Bak
oshv_Bak <- resBacTran_05_df_FULL[grepl("antagonist", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl( "killer", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]
oshv_Bak2 <- resBacTran_05_df_FULL[grepl("bcl2l7",ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names),]

#none are Bcl-2 associated protein 
oshv_BAX <- resBacTran_05_df_FULL[grepl("bax",ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names),]

#no Bok
oshv_Bok <-  resBacTran_05_df_FULL[grepl("ovarian",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) & grepl( "killer", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]
oshv_Bok <- resBacTran_05_df_FULL[grepl("ovarian",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]

#no BAD
oshv_BAD <- resBacTran_05_df_FULL[grepl("death",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) & grepl( "promoter",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]
oshv_BAD2 <- resBacTran_05_df_FULL[grepl("binding", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl( "component", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]
oshv_BAD3 <- resBacTran_05_df_FULL[grepl("BAD", ignore.case=TRUE, resBacTran_05_df_FULL$Protein.names),]

oshv_BAG <- resBacTran_05_df_FULL[grepl("BAG",ignore.case=TRUE, resBacTran_05_df_FULL$Protein.names),]

#no BH3
oshv_BH3 <- resBacTran_05_df_FULL[grepl("BH3",ignore.case=TRUE, resBacTran_05_df_FULL$Protein.names),]
oshv_Bid <- resBacTran_05_df_FULL[grepl("Bid",ignore.case=TRUE, resBacTran_05_df_FULL$Protein.names),]

#no BIM or Bcl-2-like protein 11
oshv_BIM <- resBacTran_05_df_FULL[grepl("BIM",ignore.case=TRUE, resBacTran_05_df_FULL$Protein.names),]
oshv_BIMb <- resBacTran_05_df_FULL[grepl("Bcl-2",ignore.case=TRUE, resBacTran_05_df_FULL$Protein.names),]

# no Bik
oshv_Bik <- resBacTran_05_df_FULL[grepl("interacting", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl( "killer",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]

#no Puma
oshv_Puma <- resBacTran_05_df_FULL[grepl("modulator",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) & grepl( "apoptosis",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]
oshv_Puma <- resBacTran_05_df_FULL[grepl("JFY", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]

#no ceramide
oshv_ceramide <- resBacTran_05_df_FULL[grepl("ceramide",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) ,]

oshv_MEKK_MAPK <- resBacTran_05_df_FULL[grepl("mitogen-activated", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                      grepl( "kinase", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                      !grepl("kinase-binding", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]

#no p38 MAPK
oshv_p38_MAPK <-resBacTran_05_df_FULL[grepl("mitogen-activated",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) &
                                    grepl( "kinase",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) &
                                    grepl("p38", ignore.case=TRUE, resBacTran_05_df_FULL$Protein.names),]

#none
oshv_endoG <-resBacTran_05_df_FULL[grepl("endonuclease",ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl( "G",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]

oshv_PI3 <- resBacTran_05_df_FULL[grepl("phosphatidylinositol", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl( "3-kinase",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]

#no AKT
oshv_AKTB <-resBacTran_05_df_FULL[grepl("protein", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names) &
                                 grepl("kinase",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) &
                                 grepl("B",ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]
#no protein kinase B!
oshv_akt <- resBacTran_05_df_FULL[grepl("akt", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names),]

#no JNK
oshv_JNK <- resBacTran_05_df_FULL[grepl("stress-activated", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl( "JNK", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]


#none
oshv_p35 <- resBacTran_05_df_FULL[grepl("cyclin", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl( "dependent",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) &
                                grepl( "5",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) &
                                grepl("activator",ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]
oshv_p35 <- resBacTran_05_df_FULL[grepl("p35", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]


#oshv_cJUN <- resBacTran_05_df_FULL[grepl("jun",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) ,]
oshv_cJUN <- resBacTran_05_df_FULL[grepl("c-jun",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]


oshv_AP1 <- resBacTran_05_df_FULL[grepl("AP-1", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                  grepl("transcription",ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]
#no CD151
oshv_CD151 <- resBacTran_05_df_FULL[grepl("CD151", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names),]

#oshv_BTG1 <- resBacTran_05_df_FULL[grepl("B",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) & grepl("translocation", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]

oshv_BTG2 <- resBacTran_05_df_FULL[grepl("BTG", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) ,]

oshv_IAP <- resBacTran_05_df_FULL[grepl("inhibitor", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                grepl("apoptosis", ignore.case=TRUE, resBacTran_05_df_FULL$Protein.names) &
                                !grepl("death-associated", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names) &
                                !grepl("caspase", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                !grepl("5-like", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                !grepl("TP53",resBacTran_05_df_FULL$Protein.names),]

oshv_IAP2 <- resBacTran_05_df_FULL[grepl("baculoviral", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names) &
                                 grepl("IAP", ignore.case = TRUE,resBacTran_05_df_FULL$Protein.names),]

#none
oshv_API <- resBacTran_05_df_FULL[grepl("inhibitor",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) &
                                grepl("apoptosis", ignore.case=TRUE, resBacTran_05_df_FULL$Protein.names) &
                                  !grepl("death", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names) &
                                !grepl("caspase",ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names) &
                                grepl("5", ignore.case = TRUE,resBacTran_05_df_FULL$Protein.names),]
#none
oshv_DIAP <- resBacTran_05_df_FULL[grepl("inhibitor", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                 grepl("apoptosis", ignore.case=TRUE, resBacTran_05_df_FULL$Protein.names) &
                                 grepl("death", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names),]
#none
oshv_CAAP <- resBacTran_05_df_FULL[grepl("inhibitor", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                 grepl("apoptosis", ignore.case=TRUE, resBacTran_05_df_FULL$Protein.names) &
                                 !grepl("death", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names) &
                                 grepl("caspase", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]

oshv_PCD <- resBacTran_05_df_FULL[grepl("programmed",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) & 
                                grepl("death",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) &
                                grepl("cell", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                !grepl("interacting",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]
#none
oshv_NR13 <- resBacTran_05_df_FULL[grepl("NR13", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]


oshv_PAWR <- resBacTran_05_df_FULL[grepl("WT1", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl("regulator", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]

oshv_IFI44 <- resBacTran_05_df_FULL[grepl("interferon", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl("44",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]

#no IFI6
oshv_IFI6 <- resBacTran_05_df_FULL[grepl("interferon",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) & grepl("6",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]

oshv_TLR <- resBacTran_05_df_FULL[grepl("toll",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) & grepl("receptor", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]

oshv_MyD88 <-resBacTran_05_df_FULL[grepl("myeloid", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl("differentiation",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]


oshv_CCAR <- resBacTran_05_df_FULL[grepl("cell", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) & grepl("division", ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names) &
                                 grepl("apoptosis",ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]

#none
oshv_ASPP <- resBacTran_05_df_FULL[grepl("apoptosis-stimulating",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) & grepl("p53",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names),]


#oshv_GIMAP <- resBacTran_05_df_FULL[grepl("GIMAP",ignore.case= TRUE, resBacTran_05_df_FULL$Protein.names) ,]
oshv_GIMAP <- resBacTran_05_df_FULL[grepl("GTPase", ignore.case = TRUE, resBacTran_05_df_FULL$Protein.names) &
                                  grepl("IMAP",ignore.case= TRUE,resBacTran_05_df_FULL$Protein.names),]


#Combine all lists from above into one new big table
#only include those that have hits and are confirmed to be correct
resBacTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED <- rbind(oshv_IP3R, oshv_Rho, oshv_Rac1, oshv_cdc42, oshv_PKC, 
                                                          oshv_PKCdelta, oshv_IkB, oshv_GPCRsevenpass, oshv_GPCR,
                                                          oshv_RIPK, oshv_PKA, oshv_AC, oshv_CREB, oshv_Fas, 
                                                          oshv_TNF, oshv_TNFL, oshv_TNFR, oshv_caspase, oshv_TRAF, 
                                                          oshv_NuMa, oshv_BI, oshv_cytoc, oshv_AIF, oshv_myc, 
                                                          oshv_bcl, oshv_bclw, oshv_BAG, oshv_MEKK_MAPK, 
                                                          oshv_PI3, oshv_cJUN, oshv_AP1, oshv_BTG2, oshv_IAP, 
                                                          oshv_IAP2, oshv_PCD, oshv_PAWR, oshv_IFI44, oshv_TLR, 
                                                          oshv_MyD88, oshv_CCAR, oshv_GIMAP)

#subset only those that are significant
resBacTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_significant <- resBacTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED %>%
  filter(Significance == "significant")
#add column with apoptosis pathwway and label genes with same name with extra (#)
write.csv(resBacTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_significant, file="resBacTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_significant.csv")
resBacTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_significant_path <- read.csv(file="resBacTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_significant.csv", header=TRUE)

####Graphing Significant BAC apoptosis genes p <0.05 ####
#labels for signifiance with *= p <0.05
label.Bac.sig.df <- data.frame(Protein.names= c("Rho-related GTP-binding protein RhoU",
                                                "Rho-related GTP-binding protein RhoJ",
                                                "Transmembrane BAX inhibitor motif-containing protein 4",
                                                "Phosphatidylinositol-4,5-bisphosphate 3-kinase catalytic subunit alpha isoform",
                                                "Apoptosis 2 inhibitor",
                                                "Baculoviral IAP repeat-containing protein 7-A",
                                                "Baculoviral IAP repeat-containing protein 7-B"), log2FoldChange=c(2.5,4.5,23,23.5,22,19.5,9.5))

Bac_sig_05_1 <- ggplot(resBacTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_significant_path) + 
  geom_col(aes(x=Protein.names, y=log2FoldChange, fill=Type)) + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + geom_hline(yintercept=0.0) +
  scale_y_continuous(name ="log2 Fold Change", breaks = scales::pretty_breaks(n = 20)) + 
  ggtitle("Apoptosis Genes in M. gigas after Bacterial Challenge") + 
  theme(plot.title = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5))

Bac_sig_05_2 <- Bac_sig_05_1 +
  geom_text(data=label.Bac.sig.df, aes(x=Protein.names, y=log2FoldChange), label = c("*")) +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  scale_x_discrete(name="Transcript Annotation", limits=c("Rho-related GTP-binding protein RhoU",
                                                          "Rho-related GTP-binding protein RhoJ",
                                                          "Transmembrane BAX inhibitor motif-containing protein 4",
                                                          "Phosphatidylinositol-4,5-bisphosphate 3-kinase catalytic subunit alpha isoform",
                                                          "Apoptosis 2 inhibitor",
                                                          "Baculoviral IAP repeat-containing protein 7-A",
                                                          "Baculoviral IAP repeat-containing protein 7-B"), labels=c("Rho-related GTP-binding protein RhoU"="RhoU",
                                                                                                                     "Rho-related GTP-binding protein RhoJ"="RhoJ",
                                                                                                                     "Transmembrane BAX inhibitor motif-containing protein 4"="BAX motif-containing protein 4",
                                                                                                                     "Phosphatidylinositol-4,5-bisphosphate 3-kinase catalytic subunit alpha isoform"="PI3K catalytic subunit alpha isoform",
                                                                                                                     "Apoptosis 2 inhibitor"="IAP2",
                                                                                                                     "Baculoviral IAP repeat-containing protein 7-A"="BIR IAP repeat-containing protein 7-A",
                                                                                                                     "Baculoviral IAP repeat-containing protein 7-B"="BIR IAP repeat-containing protein 7-B"))

#Apoptosis genes with p <0.1
resBacTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_01 <- resBacTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED  %>% filter(padj < 0.1)
#add pathway type to this
resBacTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_01
#add column with apoptosis pathwway and label genes with same name with extra (#)
write.csv(resBacTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_01, file="resBacTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_01.csv")
resBacTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_01_path <- read.csv(file="resBacTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_01.csv", header=TRUE)

#labels for signifiance with *= p <0.05
label.Bac.01.df <- data.frame(Protein.names= c("Rho-related GTP-binding protein RhoU",
                                               "Rho-related GTP-binding protein RhoJ",
                                               "Transmembrane BAX inhibitor motif-containing protein 4",
                                               "Phosphatidylinositol-4,5-bisphosphate 3-kinase catalytic subunit alpha isoform",
                                               "Apoptosis 2 inhibitor",
                                               "Baculoviral IAP repeat-containing protein 7-A",
                                               "Baculoviral IAP repeat-containing protein 7-B"), log2FoldChange=c(2.5,4.5,23,23.5,22,19.5,9.5))
Bac_sig_01 <- ggplot(resBacTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_01_path) + 
  geom_col(aes(x=Protein.names, y=log2FoldChange, fill=Type)) + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + geom_hline(yintercept=0.0) +
  scale_y_continuous(name ="log2 Fold Change", breaks = scales::pretty_breaks(n = 20)) + 
  ggtitle("Apoptosis Genes in M. gigas after Bacterial Challenge") + 
  theme(plot.title = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5))

Bac_sig01_2 <- Bac_sig_01 +
  geom_text(data=label.Bac.01.df, aes(x=Protein.names, y=log2FoldChange), label = c("*")) +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  scale_x_discrete(name="Transcript Annotation", limits=c("Inositol 1,4,5-trisphosphate receptor type 1",
                                                          "Rho-related GTP-binding protein RhoU",
                                                          "Rho-related GTP-binding protein RhoJ",
                                                          "Cadherin EGF LAG seven-pass G-type receptor 3",
                                                          "Fatty acid synthase",
                                                          "Lipopolysaccharide-induced tumor necrosis factor-alpha factor-like protein",
                                                          "Transmembrane BAX inhibitor motif-containing protein 4",
                                                          "Bcl-2-like protein 13",
                                                          "Phosphatidylinositol-4,5-bisphosphate 3-kinase catalytic subunit alpha isoform",
                                                          "Apoptosis 2 inhibitor",
                                                          "Baculoviral IAP repeat-containing protein 7-A",
                                                          "Baculoviral IAP repeat-containing protein 7-B",
                                                          "Myeloid differentiation primary response protein MyD88"), labels=c("Inositol 1,4,5-trisphosphate receptor type 1"="IP3R",
                                                                                                                              "Rho-related GTP-binding protein RhoU"="RhoU",
                                                                                                                              "Rho-related GTP-binding protein RhoJ"="RhoJ",
                                                                                                                              "Cadherin EGF LAG seven-pass G-type receptor 3"="Cadherin EGF LAG seven-pass GPCR3",
                                                                                                                              "Fatty acid synthase"="Fas",
                                                                                                                              "Lipopolysaccharide-induced tumor necrosis factor-alpha factor-like protein"="LPS-induced TNF-alpha factor-like",
                                                                                                                              "Transmembrane BAX inhibitor motif-containing protein 4"="BAX motif-containing protein 4",
                                                                                                                              "Bcl-2-like protein 13"="Bcl-2-like protein 13",
                                                                                                                              "Phosphatidylinositol-4,5-bisphosphate 3-kinase catalytic subunit alpha isoform"="PI3K catalytic subunit alpha isoform",
                                                                                                                              "Apoptosis 2 inhibitor"="IAP2",
                                                                                                                              "Baculoviral IAP repeat-containing protein 7-A"="BIR IAP repeat-containing protein 7-A",
                                                                                                                              "Baculoviral IAP repeat-containing protein 7-B"="BIR IAP repeat-containing protein 7-A",
                                                                                                                              "Myeloid differentiation primary response protein MyD88"="MyD88"))

####GIMAP IAP graphing for Bacterial challenge ####
# graphing all of them regardless of significance
Bac_GIMAP_IAP <- rbind(oshv_GIMAP, oshv_IAP, oshv_IAP2)
ncol(Bac_GIMAP_IAP) 
write.csv(Bac_GIMAP_IAP, file="Bac_GIMAP_IAP.csv") # add Type column
Bac_GIMAP_IAP_updated <- read.csv(file="Bac_GIMAP_IAP_type.csv",header =TRUE)
label.Bac.GIMAP.IAP.01.df <- data.frame(Protein.names= c("Apoptosis 2 inhibitor (3)",
                                                         "Baculoviral IAP repeat-containing protein 7-A (2)",
                                                         "Baculoviral IAP repeat-containing protein 7-B (2)"), log2FoldChange=c(22,19.5,9.5))
Bac_GIMAP_IAP_1 <- ggplot(Bac_GIMAP_IAP_updated) + 
  geom_col(aes(x=Protein.names, y=log2FoldChange, fill=Type)) + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + geom_hline(yintercept=0.0) +
  scale_y_continuous(name ="log2 Fold Change", breaks = scales::pretty_breaks(n = 20)) + 
  ggtitle("GIMAP and IAP Apoptosis Genes in M. gigas after Bacterial Challenge") + 
  theme(plot.title = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5))

#2 and 6 are too small to show up 
Bac_GIMAP_IAP_2 <- Bac_GIMAP_IAP_1  +
  geom_text(data=label.Bac.GIMAP.IAP.01.df, aes(x=Protein.names, y=log2FoldChange), label = c("*")) +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  scale_x_discrete(name="Transcript Annotation", limits=c("Apoptosis 1 inhibitor (1)",
                                                          "Apoptosis 1 inhibitor (2)",
                                                          "Apoptosis 2 inhibitor (1)",
                                                          "Apoptosis 2 inhibitor (2)",
                                                          "Apoptosis 2 inhibitor (3)",
                                                          "Baculoviral IAP repeat-containing protein 2",
                                                          "Baculoviral IAP repeat-containing protein 3 (1)",
                                                          "Baculoviral IAP repeat-containing protein 3 (2)",
                                                          "Baculoviral IAP repeat-containing protein 3 (Fragment)",
                                                          "Baculoviral IAP repeat-containing protein 6",
                                                          "Baculoviral IAP repeat-containing protein 7",
                                                          "Baculoviral IAP repeat-containing protein 7-A (1)",
                                                          "Baculoviral IAP repeat-containing protein 7-A (2)",
                                                          "Baculoviral IAP repeat-containing protein 7-B (1)",
                                                          "Baculoviral IAP repeat-containing protein 7-B (2)",
                                                          "GTPase IMAP family member 4 (1)",
                                                          "GTPase IMAP family member 4 (2)",
                                                          "GTPase IMAP family member 4 (3)",
                                                          "GTPase IMAP family member 4 (4)",
                                                          "GTPase IMAP family member 4 (5)",
                                                          "GTPase IMAP family member 4 (6)",
                                                          "GTPase IMAP family member 4 (7)",
                                                          "GTPase IMAP family member 4 (8)",
                                                          "GTPase IMAP family member 4  (9)",
                                                          "GTPase IMAP family member 7 (1)",
                                                          "GTPase IMAP family member 7 (2)",
                                                          "Putative apoptosis inhibitor ORF42",
                                                          "Putative inhibitor of apoptosis (1)",
                                                          "Putative inhibitor of apoptosis (2)",
                                                          "Putative inhibitor of apoptosis (3)"), labels= c("Apoptosis 1 inhibitor (1)" ="IAP1",
                                                                                                            "Apoptosis 1 inhibitor (2)"="IAP1",
                                                                                                            "Apoptosis 2 inhibitor (1)"="IAP2",
                                                                                                            "Apoptosis 2 inhibitor (2)"="IAP2",
                                                                                                            "Apoptosis 2 inhibitor (3)"="IAP2",
                                                                                                            "Baculoviral IAP repeat-containing protein 2"="BIR IAP repeat-containing 2",
                                                                                                            "Baculoviral IAP repeat-containing protein 3 (1)"="BIR IAP repeat-containing 3",
                                                                                                            "Baculoviral IAP repeat-containing protein 3 (2)"="BIR IAP repeat-containing 3",
                                                                                                            "Baculoviral IAP repeat-containing protein 3 (Fragment)"="BIR IAP repeat-containing 3 (fragment)",
                                                                                                            "Baculoviral IAP repeat-containing protein 6"="BIR IAP repeat-containing 6",
                                                                                                            "Baculoviral IAP repeat-containing protein 7"="BIR IAP repeat-containing 7",
                                                                                                            "Baculoviral IAP repeat-containing protein 7-A (1)"="BIR IAP repeat-containing 7-A",
                                                                                                            "Baculoviral IAP repeat-containing protein 7-A (2)"="BIR IAP repeat-containing 7-A",
                                                                                                            "Baculoviral IAP repeat-containing protein 7-B (1)"="BIR IAP repeat-containing 7-B",
                                                                                                            "Baculoviral IAP repeat-containing protein 7-B (2)"="BIR IAP repeat-containing 7-B",
                                                                                                            "GTPase IMAP family member 4 (1)"="GIMAP4",
                                                                                                            "GTPase IMAP family member 4 (2)"="GIMAP4",
                                                                                                            "GTPase IMAP family member 4 (3)"="GIMAP4",
                                                                                                            "GTPase IMAP family member 4 (4)"="GIMAP4",
                                                                                                            "GTPase IMAP family member 4 (5)"="GIMAP4",
                                                                                                            "GTPase IMAP family member 4 (6)"="GIMAP4",
                                                                                                            "GTPase IMAP family member 4 (7)"="GIMAP4",
                                                                                                            "GTPase IMAP family member 4 (8)"="GIMAP4",
                                                                                                            "GTPase IMAP family member 4  (9)"="GIMAP4",
                                                                                                            "GTPase IMAP family member 7 (1)"="GIMAP7",
                                                                                                            "GTPase IMAP family member 7 (2)"="GIMAP7",
                                                                                                            "Putative apoptosis inhibitor ORF42"="putative IAP ORF42",
                                                                                                            "Putative inhibitor of apoptosis (1)"="putative IAP",
                                                                                                            "Putative inhibitor of apoptosis (2)"="putative IAP",
                                                                                                      "Putative inhibitor of apoptosis (3)"="putative IAP"))
#### GIMAP and IAP across challenges #### FIXED ERROR FROM HERE DOWN
Bac_GIMAP_IAP <- droplevels(Bac_GIMAP_IAP)
as.data.frame(table(Bac_GIMAP_IAP$Protein.names))
oshv_GIMAP_IAP <- droplevels(oshv_GIMAP_IAP)
as.data.frame(table(oshv_GIMAP_IAP$Protein.names))


####Gather GIMAP, TLR, and IAP gene family information 
head(GIMAP)
GIMAP["Challenge_Type"] <- "Bac"
TLR["Challenge_Type"] <- "Bac"
Bac_IAP_combined <- rbind(oshv_IAP, oshv_IAP2)
Bac_IAP_combined["Challenge_Type"] <- "Bac"

oshv_GIMAP["Challenge_Type"] <- "oshv"
oshv_TLR["Challenge_Type"] <- "oshv"
oshv_IAP_combined <- rbind(IAP, IAP2)
oshv_IAP_combined["Challenge_Type"] <- "oshv" 

TLR_total <- rbind(TLR, oshv_TLR[,c(1,2,3,4,5,6,7,8,9,10,11,13)])
GIMAP_total <- rbind(GIMAP, oshv_GIMAP[,c(1,2,3,4,5,6,7,8,9,10,11,13)])
IAP_total <- rbind(Bac_IAP_combined, oshv_IAP_combined[,c(1,2,3,4,5,6,7,8,9,10,11,13)])

write.csv(TLR_total, file="M_gig_TLR_total.csv")
write.csv(GIMAP_total, file="M_gig_GIMAP_total.csv")
write.csv(IAP_total, file="M_gig_IAP_total.csv")

#references: 
#https://github.com/sr320/LabDocs/tree/master/code/DESeq
#StringTie manual : http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4743157/
