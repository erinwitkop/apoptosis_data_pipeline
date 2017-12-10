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
#none
PKCiota <- resoshvTran_05_df_FULL[grepl("protein",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) & 
                                                grepl("kinase",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) & grepl("iota", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names),]
#none
PKCdelta <- resoshvTran_05_df_FULL[grepl("protein", resoshvTran_05_df_FULL$Protein.names) & 
                                                 grepl("kinase",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) & grepl("delta", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) &
                                                 !grepl("calmodulin", ignore.case= TRUE,resoshvTran_05_df_FULL$Protein.names) & !grepl("ribosomal",ignore.case= TRUE, resoshvTran_05_df_FULL$Protein.names) ,]

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
                                                                         CCAR, GIMAP)

#subset only those that are significant
resoshvTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_significant <- resoshvTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED %>%
  filter(Significance == "significant")
#add column with apoptosis pathwway and label genes with same name with extra (#)
write.csv(resoshvTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_significant, file="resoshvTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_significant.csv")
resoshvTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_significant_path <- read.csv(file="resoshvTran_05_FULL_APOPTOSIS_ALL_TERMS_COMBINED_significant.csv", header=TRUE)

#only 1 was significant without the P-value correction...


#### Vibrio PROBIOTIC CHALLENGED Differential Gene Expression Analysis ####
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
write.csv(resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED, file="resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED.csv")
#subset only those that are significant
resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED_significant <- resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED[which(resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED$padj < 0.05),]

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


####Add to graph those with p value of < 0.1 rather than 0.05 ####

resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED_01 <- resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED[which(resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED$padj < 0.1),]
write.csv(resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED_01, file="resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED_01.csv")
resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED_01_path <- read.csv("resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED_01.csv", header=TRUE)

#Significance level of <=0.1 = *, less than 0.05= **, less than 10^-5= ***
label.RIF.df <- data.frame(description= c("mitogen-activated protein kinase kinase kinase 7-like",
                                          "GTPase IMAP family member 4-like"), log2FoldChange=c(-13.5, 11.0)) #less than 10^-5
label.RIF.df2 <- data.frame(description=c("apoptosis regulator BAX-like", "toll-like receptor 6"), log2FoldChange=c(2.5, 3.5)) #less than 0.1
label.RIF.df3 <- data.frame(description=c("leucine-rich repeat-containing G-protein coupled receptor 4-like",
                                          "fatty acid synthase-like",
                                          "receptor-interacting serine/threonine-protein kinase 4-like"), log2FoldChange=c(7.5, 3.0, 6.0)) #less than 0.05

RIF_Sig_PLOT_01 <- ggplot(resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED_01_path) + 
  geom_col(aes(x=description, y=log2FoldChange, fill=Pathway)) + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + geom_hline(yintercept=0.0) +
  scale_y_continuous(name ="log2 Fold Change", breaks = scales::pretty_breaks(n = 20)) + 
  ggtitle("Significantly Differentially Expressed Apoptosis Genes in C. virginica after RIF Challenge (p < 0.1)") + 
  theme(plot.title = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5))

RIF_Sig_PLOT2_01<- RIF_Sig_PLOT_01 + geom_text(data=label.RIF.df, aes(x=description, y=log2FoldChange), label = c("***")) +
  geom_text(data=label.RIF.df2, aes(x=description, y=log2FoldChange), label = c("*")) +
  geom_text(data=label.RIF.df3, aes(x=description, y=log2FoldChange), label = c("**")) +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  scale_x_discrete(name="Description", limits=c("leucine-rich repeat-containing G-protein coupled receptor 4-like",
                                                "receptor-interacting serine/threonine-protein kinase 4-like",
                                                "fatty acid synthase-like",
                                                "apoptosis regulator BAX-like",
                                                "mitogen-activated protein kinase kinase kinase 7-like",
                                                "toll-like receptor 6",
                                                "GTPase IMAP family member 4-like"), 
                   label=c("leucine-rich repeat-containing G-protein coupled receptor 4-like"="GPCR 4-like",
                           "receptor-interacting serine/threonine-protein kinase 4-like"="RIPK 4-like",
                           "fatty acid synthase-like"="Fas-like",
                           "apoptosis regulator BAX-like"="BAX-like",
                           "mitogen-activated protein kinase kinase kinase 7-like"="MEKK/MAPK 7-like",
                           "toll-like receptor 6"="TLR6",
                           "GTPase IMAP family member 4-like"="GIMAP4-like"))

####ROD GIMAP and IAP transcripts ####

RIF_GIMAP_IAP_combined <- rbind(RIF_GIMAP, RIF_IAP, RIF_IAP2)
#put two columns together to make the description column unique for the graphing process, use tidyr unite()
RIF_GIMAP_IAP_combined_2 <- unite(RIF_GIMAP_IAP_combined, description2, c(description, GeneID), remove=FALSE)
write.csv(RIF_GIMAP_IAP_combined_2, file= "RIF_GIMAP_IAP_combined_2.csv")
RIF_GIMAP_IAP_combined_3<- read.csv(file="RIF_GIMAP_IAP_combined_2.csv", header=TRUE) #make a new type column

#labels for signifiance with *= p <0.05
label.RIF.GIMAP_IAP.df <- data.frame(description2= c("GTPase IMAP family member 4-like_111105744"), log2FoldChange=c(10.5))
RIF_GIMAP_IAP_1 <- ggplot(RIF_GIMAP_IAP_combined_3) + 
  geom_col(aes(x=description2, y=log2FoldChange, fill=Type)) + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + geom_hline(yintercept=0.0) +
  scale_y_continuous(name ="log2 Fold Change", breaks = scales::pretty_breaks(n = 20)) + 
  ggtitle("GIMAP and IAP Apoptosis Genes in C. virginica after RIF Challenge") + 
  theme(plot.title = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5))

RIF_GIMAP_IAP_2 <- RIF_GIMAP_IAP_1 +
  geom_text(data=label.RIF.GIMAP_IAP.df, aes(x=description2, y=log2FoldChange), label = c("*")) +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  scale_x_discrete(name="Transcript Annotation", limits=c("baculoviral IAP repeat-containing protein 2-like_111100408 (1)",
                                                          "baculoviral IAP repeat-containing protein 2-like_111100408 (2)",
                                                          "baculoviral IAP repeat-containing protein 2-like_111100443",
                                                          "baculoviral IAP repeat-containing protein 2-like_111101035",
                                                          "baculoviral IAP repeat-containing protein 2-like_111101689 (1)",
                                                          "baculoviral IAP repeat-containing protein 2-like_111101689 (2)",
                                                          "baculoviral IAP repeat-containing protein 2-like_111101689 (2)",
                                                          "baculoviral IAP repeat-containing protein 2-like_111102770 (1)",
                                                          "baculoviral IAP repeat-containing protein 2-like_111102770 (2)",
                                                          "baculoviral IAP repeat-containing protein 2-like_111102858",
                                                          "baculoviral IAP repeat-containing protein 2-like_111102964",
                                                          "baculoviral IAP repeat-containing protein 2-like_111103826 (1)",
                                                          "baculoviral IAP repeat-containing protein 2-like_111103826 (2)",
                                                          "baculoviral IAP repeat-containing protein 2-like_111105503",
                                                          "baculoviral IAP repeat-containing protein 2-like_111106726",
                                                          "baculoviral IAP repeat-containing protein 2-like_111123894",
                                                          "baculoviral IAP repeat-containing protein 3-like_111100470 (1)",
                                                          "baculoviral IAP repeat-containing protein 3-like_111100470 (2)",
                                                          "baculoviral IAP repeat-containing protein 3-like_111100471 (1)",
                                                          "baculoviral IAP repeat-containing protein 3-like_111100471 (2)",
                                                          "baculoviral IAP repeat-containing protein 3-like_111100471 (3)",
                                                          "baculoviral IAP repeat-containing protein 3-like_111101018",
                                                          "baculoviral IAP repeat-containing protein 3-like_111101864",
                                                          "baculoviral IAP repeat-containing protein 3-like_111103155",
                                                          "baculoviral IAP repeat-containing protein 3-like_111103270 (1)",
                                                          "baculoviral IAP repeat-containing protein 3-like_111103270 (2)",
                                                          "baculoviral IAP repeat-containing protein 3-like_111103270 (3)",
                                                          "baculoviral IAP repeat-containing protein 3-like_111103270 (4)",
                                                          "baculoviral IAP repeat-containing protein 3-like_111103392 (1)",
                                                          "baculoviral IAP repeat-containing protein 3-like_111103392 (2)",
                                                          "baculoviral IAP repeat-containing protein 3-like_111103427 (1)",
                                                          "baculoviral IAP repeat-containing protein 3-like_111103427 (2)",
                                                          "baculoviral IAP repeat-containing protein 3-like_111103816",
                                                          "baculoviral IAP repeat-containing protein 3-like_111104430 (1)",
                                                          "baculoviral IAP repeat-containing protein 3-like_111104430 (2)",
                                                          "baculoviral IAP repeat-containing protein 3-like_111104430 (3)",
                                                          "baculoviral IAP repeat-containing protein 3-like_111105157",
                                                          "baculoviral IAP repeat-containing protein 6-like_111129365 (1)",
                                                          "baculoviral IAP repeat-containing protein 6-like_111129365 (2)",
                                                          "baculoviral IAP repeat-containing protein 6-like_111129365 (3)",
                                                          "baculoviral IAP repeat-containing protein 6-like_111129365 (4)",
                                                          "baculoviral IAP repeat-containing protein 6-like_111130310 (1)",
                                                          "baculoviral IAP repeat-containing protein 6-like_111130310 (2)",
                                                          "baculoviral IAP repeat-containing protein 6-like_111130310 (3)",
                                                          "baculoviral IAP repeat-containing protein 6-like_111130310 (4)",
                                                          "baculoviral IAP repeat-containing protein 7-A-like_111100432",
                                                          "baculoviral IAP repeat-containing protein 7-A-like_111103982",
                                                          "baculoviral IAP repeat-containing protein 7-A-like_111104279",
                                                          "baculoviral IAP repeat-containing protein 7-A-like_111104280",
                                                          "baculoviral IAP repeat-containing protein 7-B-like_111105148",
                                                          "baculoviral IAP repeat-containing protein 7-like_111100416",
                                                          "baculoviral IAP repeat-containing protein 7-like_111100417",
                                                          "baculoviral IAP repeat-containing protein 7-like_111100633 (1)",
                                                          "baculoviral IAP repeat-containing protein 7-like_111100633 (2)",
                                                          "baculoviral IAP repeat-containing protein 7-like_111105494 (1) ",
                                                          "baculoviral IAP repeat-containing protein 7-like_111105494 (2)",
                                                          "baculoviral IAP repeat-containing protein 7-like_111105494 (3)",
                                                          "baculoviral IAP repeat-containing protein 8-like_111103158",
                                                          "baculoviral IAP repeat-containing protein 8-like_111104229",
                                                          "baculoviral IAP repeat-containing protein 8-like_111109770 (1)",
                                                          "baculoviral IAP repeat-containing protein 8-like_111109770 (2)",
                                                          "GTPase IMAP family member 4-like_111100020 (1)",
                                                          "GTPase IMAP family member 4-like_111100020 (2)",
                                                          "GTPase IMAP family member 4-like_111105195",
                                                          "GTPase IMAP family member 4-like_111105333 (1)",
                                                          "GTPase IMAP family member 4-like_111105333 (2)",
                                                          "GTPase IMAP family member 4-like_111105333 (3)",
                                                          "GTPase IMAP family member 4-like_111105335",
                                                          "GTPase IMAP family member 4-like_111105336",
                                                          "GTPase IMAP family member 4-like_111105339",
                                                          "GTPase IMAP family member 4-like_111105744",
                                                          "GTPase IMAP family member 4-like_111105930",
                                                          "GTPase IMAP family member 4-like_111106079",
                                                          "GTPase IMAP family member 4-like_111106328",
                                                          "GTPase IMAP family member 4-like_111106343",
                                                          "GTPase IMAP family member 4-like_111107002",
                                                          "GTPase IMAP family member 4-like_111108253",
                                                          "GTPase IMAP family member 4-like_111109344",
                                                          "GTPase IMAP family member 4-like_111109667",
                                                          "GTPase IMAP family member 4-like_111109737",
                                                          "GTPase IMAP family member 4-like_111111775",
                                                          "GTPase IMAP family member 4-like_111119582 (1)",
                                                          "GTPase IMAP family member 4-like_111119582 (2)",
                                                          "GTPase IMAP family member 4-like_111119582 (3)",
                                                          "GTPase IMAP family member 4-like_111129930 (1)",
                                                          "GTPase IMAP family member 4-like_111129930 (2)",
                                                          "GTPase IMAP family member 4-like_111129932 (1)",
                                                          "GTPase IMAP family member 4-like_111129932 (2)",
                                                          "GTPase IMAP family member 4-like_111130153 (1)",
                                                          "GTPase IMAP family member 4-like_111130153 (2)",
                                                          "GTPase IMAP family member 4-like_111130153 (3)",
                                                          "GTPase IMAP family member 4-like_111130153 (4)",
                                                          "GTPase IMAP family member 4-like_111130153 (5)",
                                                          "GTPase IMAP family member 4-like_111130153 (6)",
                                                          "GTPase IMAP family member 4-like_111130155 (1)",
                                                          "GTPase IMAP family member 4-like_111130155 (2)",
                                                          "GTPase IMAP family member 7-like_111108121",
                                                          "GTPase IMAP family member 7-like_111108220",
                                                          "GTPase IMAP family member 7-like_111108559",
                                                          "GTPase IMAP family member 7-like_111109557 (1)",
                                                          "GTPase IMAP family member 7-like_111109557 (2)",
                                                          "GTPase IMAP family member 7-like_111110321",
                                                          "GTPase IMAP family member 8-like_111119581",
                                                          "GTPase IMAP family member 8-like_111120314",
                                                          "inhibitor of apoptosis protein-like_111133238",
                                                          "putative inhibitor of apoptosis_111100893",
                                                          "putative inhibitor of apoptosis_111100894",
                                                          "putative inhibitor of apoptosis_111101016 (1)",
                                                          "putative inhibitor of apoptosis_111101016 (2)",
                                                          "putative inhibitor of apoptosis_111102106",
                                                          "putative inhibitor of apoptosis_111102451",
                                                          "putative inhibitor of apoptosis_111103790",
                                                          "putative inhibitor of apoptosis_111104637",
                                                          "putative inhibitor of apoptosis_111132489"), label=c("baculoviral IAP repeat-containing protein 2-like_111100408 (1)"= "BIR IAP 2-like",
                                                                                                                "baculoviral IAP repeat-containing protein 2-like_111100408 (2)"= "BIR IAP 2-like",
                                                                                                                "baculoviral IAP repeat-containing protein 2-like_111100443"= "BIR IAP 2-like",
                                                                                                                "baculoviral IAP repeat-containing protein 2-like_111101035"="BIR IAP 2-like",
                                                                                                                "baculoviral IAP repeat-containing protein 2-like_111101689 (1)"="BIR IAP 2-like",
                                                                                                                "baculoviral IAP repeat-containing protein 2-like_111101689 (2)"="BIR IAP 2-like",
                                                                                                                "baculoviral IAP repeat-containing protein 2-like_111101689 (2)"="BIR IAP 2-like",
                                                                                                                "baculoviral IAP repeat-containing protein 2-like_111102770 (1)"="BIR IAP 2-like",
                                                                                                                "baculoviral IAP repeat-containing protein 2-like_111102770 (2)"="BIR IAP 2-like",
                                                                                                                "baculoviral IAP repeat-containing protein 2-like_111102858"="BIR IAP 2-like",
                                                                                                                "baculoviral IAP repeat-containing protein 2-like_111102964"="BIR IAP 2-like",
                                                                                                                "baculoviral IAP repeat-containing protein 2-like_111103826 (1)"="BIR IAP 2-like",
                                                                                                                "baculoviral IAP repeat-containing protein 2-like_111103826 (2)"="BIR IAP 2-like",
                                                                                                                "baculoviral IAP repeat-containing protein 2-like_111105503"="BIR IAP 2-like",
                                                                                                                "baculoviral IAP repeat-containing protein 2-like_111106726"="BIR IAP 2-like",
                                                                                                                "baculoviral IAP repeat-containing protein 2-like_111123894"="BIR IAP 2-like",
                                                                                                                "baculoviral IAP repeat-containing protein 3-like_111100470 (1)"="BIR IAP 3-like",
                                                                                                                "baculoviral IAP repeat-containing protein 3-like_111100470 (2)"="BIR IAP 3-like",
                                                                                                                "baculoviral IAP repeat-containing protein 3-like_111100471 (1)"="BIR IAP 3-like",
                                                                                                                "baculoviral IAP repeat-containing protein 3-like_111100471 (2)"="BIR IAP 3-like",
                                                                                                                "baculoviral IAP repeat-containing protein 3-like_111100471 (3)"="BIR IAP 3-like",
                                                                                                                "baculoviral IAP repeat-containing protein 3-like_111101018"="BIR IAP 3-like",
                                                                                                                "baculoviral IAP repeat-containing protein 3-like_111101864"="BIR IAP 3-like",
                                                                                                                "baculoviral IAP repeat-containing protein 3-like_111103155"="BIR IAP 3-like",
                                                                                                                "baculoviral IAP repeat-containing protein 3-like_111103270 (1)"="BIR IAP 3-like",
                                                                                                                "baculoviral IAP repeat-containing protein 3-like_111103270 (2)"="BIR IAP 3-like",
                                                                                                                "baculoviral IAP repeat-containing protein 3-like_111103270 (3)"="BIR IAP 3-like",
                                                                                                                "baculoviral IAP repeat-containing protein 3-like_111103270 (4)"="BIR IAP 3-like",
                                                                                                                "baculoviral IAP repeat-containing protein 3-like_111103392 (1)"="BIR IAP 3-like",
                                                                                                                "baculoviral IAP repeat-containing protein 3-like_111103392 (2)"="BIR IAP 3-like",
                                                                                                                "baculoviral IAP repeat-containing protein 3-like_111103427 (1)"="BIR IAP 3-like",
                                                                                                                "baculoviral IAP repeat-containing protein 3-like_111103427 (2)"="BIR IAP 3-like",
                                                                                                                "baculoviral IAP repeat-containing protein 3-like_111103816"="BIR IAP 3-like",
                                                                                                                "baculoviral IAP repeat-containing protein 3-like_111104430 (1)"="BIR IAP 3-like",
                                                                                                                "baculoviral IAP repeat-containing protein 3-like_111104430 (2)"="BIR IAP 3-like",
                                                                                                                "baculoviral IAP repeat-containing protein 3-like_111104430 (3)"="BIR IAP 3-like",
                                                                                                                "baculoviral IAP repeat-containing protein 3-like_111105157"="BIR IAP 3-like",
                                                                                                                "baculoviral IAP repeat-containing protein 6-like_111129365 (1)"="BIR IAP 6-like",
                                                                                                                "baculoviral IAP repeat-containing protein 6-like_111129365 (2)"="BIR IAP 6-like",
                                                                                                                "baculoviral IAP repeat-containing protein 6-like_111129365 (3)"="BIR IAP 6-like",
                                                                                                                "baculoviral IAP repeat-containing protein 6-like_111129365 (4)"="BIR IAP 6-like",
                                                                                                                "baculoviral IAP repeat-containing protein 6-like_111130310 (1)"="BIR IAP 6-like",
                                                                                                                "baculoviral IAP repeat-containing protein 6-like_111130310 (2)"="BIR IAP 6-like",
                                                                                                                "baculoviral IAP repeat-containing protein 6-like_111130310 (3)"="BIR IAP 6-like",
                                                                                                                "baculoviral IAP repeat-containing protein 6-like_111130310 (4)"="BIR IAP 6-like",
                                                                                                                "baculoviral IAP repeat-containing protein 7-A-like_111100432"="BIR IAP 7-A-like",
                                                                                                                "baculoviral IAP repeat-containing protein 7-A-like_111103982"="BIR IAP 7-A-like",
                                                                                                                "baculoviral IAP repeat-containing protein 7-A-like_111104279"="BIR IAP 7-A-like",
                                                                                                                "baculoviral IAP repeat-containing protein 7-A-like_111104280"="BIR IAP 7-A-like",
                                                                                                                "baculoviral IAP repeat-containing protein 7-B-like_111105148"="BIR IAP 7-B-like",
                                                                                                                "baculoviral IAP repeat-containing protein 7-like_111100416"="BIR IAP 7-like",
                                                                                                                "baculoviral IAP repeat-containing protein 7-like_111100417"="BIR IAP 7-like",
                                                                                                                "baculoviral IAP repeat-containing protein 7-like_111100633 (1)"="BIR IAP 7-like",
                                                                                                                "baculoviral IAP repeat-containing protein 7-like_111100633 (2)"="BIR IAP 7-like",
                                                                                                                "baculoviral IAP repeat-containing protein 7-like_111105494 (1) "="BIR IAP 7-like",
                                                                                                                "baculoviral IAP repeat-containing protein 7-like_111105494 (2)"="BIR IAP 7-like",
                                                                                                                "baculoviral IAP repeat-containing protein 7-like_111105494 (3)"="BIR IAP 7-like",
                                                                                                                "baculoviral IAP repeat-containing protein 8-like_111103158"="BIR IAP 8-like",
                                                                                                                "baculoviral IAP repeat-containing protein 8-like_111104229"="BIR IAP 8-like",
                                                                                                                "baculoviral IAP repeat-containing protein 8-like_111109770 (1)"="BIR IAP 8-like",
                                                                                                                "baculoviral IAP repeat-containing protein 8-like_111109770 (2)"="BIR IAP 8-like",
                                                                                                                "GTPase IMAP family member 4-like_111100020 (1)"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111100020 (2)"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111105195"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111105333 (1)"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111105333 (2)"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111105333 (3)"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111105335"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111105336"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111105339"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111105744"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111105930"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111106079"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111106328"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111106343"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111107002"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111108253"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111109344"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111109667"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111109737"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111111775"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111119582 (1)"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111119582 (2)"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111119582 (3)"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111129930 (1)"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111129930 (2)"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111129932 (1)"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111129932 (2)"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111130153 (1)"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111130153 (2)"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111130153 (3)"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111130153 (4)"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111130153 (5)"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111130153 (6)"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111130155 (1)"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 4-like_111130155 (2)"="GIMAP 4-like",
                                                                                                                "GTPase IMAP family member 7-like_111108121"="GIMAP 7-like",
                                                                                                                "GTPase IMAP family member 7-like_111108220"="GIMAP 7-like",
                                                                                                                "GTPase IMAP family member 7-like_111108559"="GIMAP 7-like",
                                                                                                                "GTPase IMAP family member 7-like_111109557 (1)"="GIMAP 7-like",
                                                                                                                "GTPase IMAP family member 7-like_111109557 (2)"="GIMAP 7-like",
                                                                                                                "GTPase IMAP family member 7-like_111110321"="GIMAP 7-like",
                                                                                                                "GTPase IMAP family member 8-like_111119581"="GIMAP 8-like",
                                                                                                                "GTPase IMAP family member 8-like_111120314"="GIMAP 8-like",
                                                                                                                "inhibitor of apoptosis protein-like_111133238"="IAP-like",
                                                                                                                "putative inhibitor of apoptosis_111100893"="putative IAP",
                                                                                                                "putative inhibitor of apoptosis_111100894"="putative IAP",
                                                                                                                "putative inhibitor of apoptosis_111101016 (1)"="putative IAP",
                                                                                                                "putative inhibitor of apoptosis_111101016 (2)"="putative IAP",
                                                                                                                "putative inhibitor of apoptosis_111102106"="putative IAP",
                                                                                                                "putative inhibitor of apoptosis_111102451"="putative IAP",
                                                                                                                "putative inhibitor of apoptosis_111103790"="putative IAP",
                                                                                                                "putative inhibitor of apoptosis_111104637"="putative IAP",
                                                                                                                "putative inhibitor of apoptosis_111132489"="putative IAP"))



#### Side by side graphing of the shared significantly differentially expressed apoptosis genes from both ####
# this is wrong
Shared_apoptosis <- merge(resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED, resRODTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED, by="description")  
resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED$source <- "X"
resRODTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED$source <- "Y"
C_Vir_Merged <- merge(x = resRIFTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED, y = resRODTran_05_FULL_ENTREZGENE_DATA_APOPTOSIS_ALL_TERMS_COMBINED,
                      all = TRUE, by = c("description"))
C_Vir_Merged$rowSource <- apply(C_Vir_Merged[c("source.x", "source.y")], 1, 
                                function(x) paste(na.omit(x), collapse = ""))
C_Vir_Merged


#### Compare the number of transcripts of each one of interest


#references: 
#https://github.com/sr320/LabDocs/tree/master/code/DESeq
#StringTie manual : http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4743157/
