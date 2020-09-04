# SCRIPT TO RUN SIGNED WGCNA ON PARASITE/PATHOGEN CHALLENGES 
# Erin Roberts, PhD Candidate University of Rhode Island 
# 4/9/2020

# This script perform WGCNA analysis on C. virginica and C. gigas experiments separately 

### LOAD PACKAGES ####
library(tidyverse)
library(limma)
library(WGCNA) # v WGCNA_1.68
options(stringsAsFactors = FALSE) # run every time
allowWGCNAThreads()
library(cluster)
library(anRichment)
library(anRichmentMethods)
library(plyr)
library(dplyr)
library(magicfor)
cor <- WGCNA::cor # run every time
library(UpSetR)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(VennDiagram)
library(ComplexHeatmap)
library(pheatmap)
library(gt)
# Using R version 3.6.1


# source("https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/installAnRichment.R"); installAnRichment(); 

#### LOAD APOPTOSIS DATA FRAMES AND PACKAGES ####

#### REMEMBER TO RELOAD AND REDO THIS WITH THE NEW NOW HAPLOTIG COLLAPSED C.VIRGINICA DATA (DONE August 5th, 2020)
Apoptosis_frames <- load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gig_C_vir_apoptosis_products.RData")
annotations <- load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gig_C_vir_annotations.RData")
Apop_LFC <- load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/apop_LFC.RData")
# IAP lists with domain type 
load(file = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/IAP_domain_structure_XM_CG.RData")
load(file = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/Apoptosis_Pathway_Annotation_Comparative_Genomics/Comparative_Analysis_Apoptosis_Gene_Families_Data/IAP_domain_structure_XM_CV.RData")

#### Helpful tutorials for meta-analysis
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/ModulePreservation/
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/JMiller/Tutorial%20document.pdf
# WGCNA Main tutorial: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
# anRichment: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/#manualInstall
# Package FAQs with some quidelines for running: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html
# tutorials on module preservation: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/ModulePreservation/Tutorials/
# Code for Differential network analysis: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/DifferentialNetworkAnalysis/

#### WGCNA C_VIRGINICA ####
 # WGCNA analysis resources
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/JMiller/
# Input data is variance stabilizing transformation counts data 
# see answers to general questions: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html

# Do I filter my gene set first? https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html
# We do not recommend filtering genes by differential expression. 
# WGCNA is designed to be an unsupervised analysis method that clusters 
# genes based on their expression profiles. Filtering genes by 
# differential expression will lead to a set of correlated genes that 
# will essentially form a single (or a few highly correlated) modules. 
# It also completely invalidates the scale-free topology assumption, 
# so choosing soft thresholding power by scale-free topology fit will fail. 

# What should I use for input data? https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html
# We then recommend a variance-stabilizing transformation.
# For example, package DESeq2 implements the function 
# varianceStabilizingTransformation which we have found useful, 
# but one could also start with normalized counts (or RPKM/FPKM data)
# and log-transform them using log2(x+1). For highly expressed features,
# the differences between full variance stabilization and a simple log transformation are small.
# Counts datatable needs to be in same format as DESeq2, each row is a transcript and each column is a sample

## See tutorial by Jeremy Miller that describes his process for meta-analysis:https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/JMiller/Tutorial%20document.pdf

## Remeber: batch effects represent the systematic technical differences when samples are processed and measured in different batches and which are unrelated to any biological variation

#### LOAD SAVED DATA ####
C_vir_expression = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_virginica_normalized_exprssion_WGCNA_input.RData")
C_vir_coldata = load(file='/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_virginica_coldata_WGCNA_input.RData')
C_gig_expression = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gigas_normalized_exprssion_WGCNA_input.RData")
C_gig_coldata = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gigas_coldata_WGCNA_input.RData")

#####  DATA FORMATTING, BATCH EFFECT REMOVAL ####
#Dermo_counts: Keep the following as two separate networks, correct for library prep date 
#Dermo_Tolerant_dds_vst
#Dermo_Susceptible_dds_vst
# formula ~lib.prep.time + time + condition

Dermo_Tolerant_dds_vst_limma <- Dermo_Tolerant_dds_vst
Dermo_Susceptible_dds_vst_limma <- Dermo_Susceptible_dds_vst

plotPCA(Dermo_Tolerant_dds_vst_limma, "Time") # large clustering by time
plotPCA(Dermo_Tolerant_dds_vst_limma, "Lib_prep_date") # large technical batch effects here 
plotPCA(Dermo_Susceptible_dds_vst_limma, "Lib_prep_date") # also technical batch effects here
mat_C_vir_dermo_Tol <- assay(Dermo_Tolerant_dds_vst_limma)
mat_C_vir_dermo_Sus <- assay(Dermo_Susceptible_dds_vst_limma)
mat_C_vir_dermo_Tol <- limma::removeBatchEffect(mat_C_vir_dermo_Tol, Dermo_Tolerant_dds_vst_limma$Lib_prep_date)
mat_C_vir_dermo_Sus <- limma::removeBatchEffect(mat_C_vir_dermo_Sus, Dermo_Susceptible_dds_vst_limma$Lib_prep_date)

assay(Dermo_Tolerant_dds_vst_limma) <- mat_C_vir_dermo_Tol
assay(Dermo_Susceptible_dds_vst_limma) <- mat_C_vir_dermo_Sus
plotPCA(Dermo_Tolerant_dds_vst_limma, "Lib_prep_date") # corrected
plotPCA(Dermo_Susceptible_dds_vst_limma, "Lib_prep_date") # corrected

# save as matrix with assay and transform
Dermo_Tolerant_dds_vst_matrix <- assay(Dermo_Tolerant_dds_vst_limma)
Dermo_Susceptible_dds_vst_matrix <-assay(Dermo_Susceptible_dds_vst_limma) 
class(Dermo_Tolerant_dds_vst_matrix)
class(Dermo_Susceptible_dds_vst_matrix)

Dermo_Tolerant_dds_vst_matrix  <- t(Dermo_Tolerant_dds_vst_matrix) 
Dermo_Susceptible_dds_vst_matrix <- t(Dermo_Susceptible_dds_vst_matrix)

#Probiotic_counts: no batch effect correction 
# Probiotic_dds_rlog , ~Time + condition
plotPCA(Probiotic_dds_rlog, "Time") # large clustering by time
Probiotic_dds_rlog_matrix <- assay(Probiotic_dds_rlog)
class(Probiotic_dds_rlog_matrix)
Probiotic_dds_rlog_matrix <- t(Probiotic_dds_rlog_matrix)
head(Probiotic_dds_rlog_matrix)
ncol(Probiotic_dds_rlog_matrix) 

#ROD_counts: no batch effect correction 
# ROD_Resistant_dds_rlog ~time+condition
# ROD_Susceptible_dds_rlog ~condition (no time)
plotPCA(ROD_Resistant_dds_rlog, "Time") # clustering by time
plotPCA(ROD_Susceptible_dds_rlog, "Time") 

ROD_Resistant_dds_rlog_matrix   <- assay(ROD_Resistant_dds_rlog)
ROD_Susceptible_dds_rlog_matrix <-  assay(ROD_Susceptible_dds_rlog)
ROD_Resistant_dds_rlog_matrix  <- t(ROD_Resistant_dds_rlog_matrix)
ROD_Susceptible_dds_rlog_matrix <- t(ROD_Susceptible_dds_rlog_matrix)

#Pro_RE22_counts : batch effect correction for "Family", which is the different larval sources used on different days
# Pro_RE22_dds_rlog  ~Time+ condition
plotPCA(Pro_RE22_dds_rlog, "Family") 
Pro_RE22_dds_rlog_limma <- Pro_RE22_dds_rlog

mat_C_vir_dermo_Pro_RE22 <- assay(Pro_RE22_dds_rlog_limma)
mat_C_vir_dermo_Pro_RE22 <- limma::removeBatchEffect(mat_C_vir_dermo_Pro_RE22, Pro_RE22_dds_rlog_limma$Family)
assay(Pro_RE22_dds_rlog_limma) <- mat_C_vir_dermo_Pro_RE22
plotPCA(Pro_RE22_dds_rlog_limma, "Family") # corrected

Pro_RE22_dds_rlog_matrix <- assay(Pro_RE22_dds_rlog_limma)
Pro_RE22_dds_rlog_matrix <- t(Pro_RE22_dds_rlog_matrix)

## LIMIT ANALYSIS TO TRANSCRIPTS EXPRESSED IN ALL EXPERIMENTS
ncol(Dermo_Tolerant_dds_vst_matrix) #49995
ncol(Dermo_Susceptible_dds_vst_matrix) # 49634
ncol(Probiotic_dds_rlog_matrix) # 51052
ncol(ROD_Resistant_dds_rlog_matrix) # 47205
ncol(ROD_Susceptible_dds_rlog_matrix)
ncol(Pro_RE22_dds_rlog_matrix) # 44221

Dermo_Tolerant_dds_vst_matrix_colnames <- colnames(Dermo_Tolerant_dds_vst_matrix)
Dermo_Susceptible_dds_vst_matrix_colnames <- colnames(Dermo_Susceptible_dds_vst_matrix)
Probiotic_dds_rlog_matrix_colnames <- colnames(Probiotic_dds_rlog_matrix)
ROD_Resistant_dds_rlog_matrix_colnames <- colnames(ROD_Resistant_dds_rlog_matrix)
ROD_Susceptible_dds_rlog_matrix_colnames <- colnames(ROD_Susceptible_dds_rlog_matrix)
Pro_RE22_dds_rlog_matrix_colnames  <- colnames(Pro_RE22_dds_rlog_matrix)

C_vir_common_vst_transcripts <- Reduce(intersect, list(Dermo_Tolerant_dds_vst_matrix_colnames, Dermo_Susceptible_dds_vst_matrix_colnames,
      Probiotic_dds_rlog_matrix_colnames,ROD_Resistant_dds_rlog_matrix_colnames,ROD_Susceptible_dds_rlog_matrix_colnames, Pro_RE22_dds_rlog_matrix_colnames))
head(C_vir_common_vst_transcripts)
length(C_vir_common_vst_transcripts) # 29440

Dermo_Tolerant_dds_vst_matrix_common <- Dermo_Tolerant_dds_vst_matrix[, C_vir_common_vst_transcripts]
ncol(Dermo_Tolerant_dds_vst_matrix_common)

Dermo_Susceptible_dds_vst_matrix_common <- Dermo_Susceptible_dds_vst_matrix[,C_vir_common_vst_transcripts]
ncol(Dermo_Susceptible_dds_vst_matrix_common)

Probiotic_dds_rlog_matrix_common <- Probiotic_dds_rlog_matrix[,C_vir_common_vst_transcripts]
ncol(Probiotic_dds_rlog_matrix_common)

ROD_Resistant_dds_rlog_matrix_common <- ROD_Resistant_dds_rlog_matrix[,C_vir_common_vst_transcripts]
ncol(ROD_Resistant_dds_rlog_matrix_common)

ROD_Susceptible_dds_rlog_matrix_common <- ROD_Susceptible_dds_rlog_matrix[,C_vir_common_vst_transcripts]
ncol(ROD_Susceptible_dds_rlog_matrix_common)

Pro_RE22_dds_rlog_matrix_common <- Pro_RE22_dds_rlog_matrix[,C_vir_common_vst_transcripts]
ncol(Pro_RE22_dds_rlog_matrix_common )

# Split Pro_RE22 also into  Probiotic and RE22 so both can be run
# subset the Pro_RE22 experiment
print(Pro_RE22_coldata %>% rownames_to_column("sample") %>% filter(Condition == "Control_no_treatment" | Condition == "Vibrio_coralliilyticus_RE22_exposure_6h"))

Pro_RE22_coldata_RE22 <- Pro_RE22_coldata %>% rownames_to_column("sample") %>% filter(Condition == "Control_no_treatment" | Condition == "Vibrio_coralliilyticus_RE22_exposure_6h")
Pro_RE22_coldata_Pro <-  Pro_RE22_coldata %>% rownames_to_column("sample") %>% filter(Condition == "Control_no_treatment" | Condition == "Phaeobacter_inhibens_S4_exposure_6h" |
                    Condition == "Phaeobacter_inhibens_S4_exposure_24h" | Condition == "Bacillus_pumilus_RI06_95_exposure_6h" | Condition == "Bacillus_pumilus_RI06_95_exposure_24h")     

Pro_RE22_dds_rlog_matrix_common_Pro <- Pro_RE22_dds_rlog_matrix_common[Pro_RE22_coldata_Pro$sample,]
Pro_RE22_dds_rlog_matrix_common_RE22 <- Pro_RE22_dds_rlog_matrix_common[Pro_RE22_coldata_RE22$sample,]

Pro_RE22_dds_rlog_matrix_Pro <- Pro_RE22_dds_rlog_matrix[Pro_RE22_coldata_Pro$sample,]
Pro_RE22_dds_rlog_matrix_RE22 <- Pro_RE22_dds_rlog_matrix[Pro_RE22_coldata_RE22$sample,]

# Save these
write.table(Pro_RE22_dds_rlog_matrix_Pro , file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_dds_rlog_matrix_Pro.table")
write.table(Pro_RE22_dds_rlog_matrix_RE22, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_dds_rlog_matrix_RE22.table")

# Do column names agree between all?
all(colnames(Dermo_Tolerant_dds_vst_matrix_common ) %in% colnames(Probiotic_dds_rlog_matrix_common)) # TRUE
all(colnames(Dermo_Tolerant_dds_vst_matrix_common ) == colnames(Probiotic_dds_rlog_matrix_common)) # TRUE 

#### SELECT SOFT THRESHOLDING POWER FOR EACH EXPERIMENT ####
# Choose a set of soft-thresholding powers
 powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function, following general recommendations to set network type to "signed hybrid" and using the "bicor" correlation
# Dermo_Tolerant_dds_vst_matrix_common_sft <- pickSoftThreshold(Dermo_Tolerant_dds_vst_matrix_common, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
# Dermo_Susceptible_dds_vst_matrix_common_sft <- pickSoftThreshold(Dermo_Susceptible_dds_vst_matrix_common, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
# Probiotic_dds_rlog_matrix_common_sft <- pickSoftThreshold(Probiotic_dds_rlog_matrix_common, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
# # Warning message:
# # executing %dopar% sequentially: no parallel backend registered 
# #From Peter Langfelder: https://bioinformatics.stackexchange.com/questions/10555/r-wgcna-error-code:
# #What you see is a warning, not an error. 
# #Your calculation will run fine, just slower. Unless you see other errors, you should be able to complete all steps of the analysis.
# 
# ROD_Resistant_dds_rlog_matrix_common_sft <- pickSoftThreshold(ROD_Resistant_dds_rlog_matrix_common, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
# ROD_Susceptible_dds_rlog_matrix_common_sft <- pickSoftThreshold(ROD_Susceptible_dds_rlog_matrix_common, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
# Pro_RE22_dds_rlog_matrix_common_sft <- pickSoftThreshold(Pro_RE22_dds_rlog_matrix_common, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
# Pro_RE22_dds_rlog_matrix_common_Pro_sft <- pickSoftThreshold(Pro_RE22_dds_rlog_matrix_common_Pro, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
# Pro_RE22_dds_rlog_matrix_common_RE22_sft <-  pickSoftThreshold(Pro_RE22_dds_rlog_matrix_common_RE22, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
# # export pick softthresholding so it is saved 
# save(Dermo_Tolerant_dds_vst_matrix_common_sft, Dermo_Susceptible_dds_vst_matrix_common_sft, Probiotic_dds_rlog_matrix_common_sft, ROD_Resistant_dds_rlog_matrix_common_sft,
# ROD_Susceptible_dds_rlog_matrix_common_sft, Pro_RE22_dds_rlog_matrix_common_sft, Pro_RE22_dds_rlog_matrix_common_Pro_sft, Pro_RE22_dds_rlog_matrix_common_RE22_sft,
# file= "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_virginica_individual_pickSoftThreshold_WGCNA_input.RData")

# Load saved data
C_vir_picksoftthreshold = load(file= "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_virginica_individual_pickSoftThreshold_WGCNA_input.RData")

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

#Dermo
# Scale-free topology fit index as a function of the soft-thresholding power
plot(Dermo_Tolerant_dds_vst_matrix_common_sft$fitIndices[,1], -sign(Dermo_Tolerant_dds_vst_matrix_common_sft$fitIndices[,3])*Dermo_Tolerant_dds_vst_matrix_common_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(Dermo_Tolerant_dds_vst_matrix_common_sft$fitIndices[,1], -sign(Dermo_Tolerant_dds_vst_matrix_common_sft$fitIndices[,3])*Dermo_Tolerant_dds_vst_matrix_common_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(Dermo_Tolerant_dds_vst_matrix_common_sft$fitIndices[,1], Dermo_Tolerant_dds_vst_matrix_common_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(Dermo_Tolerant_dds_vst_matrix_common_sft$fitIndices[,1], Dermo_Tolerant_dds_vst_matrix_common_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Softthreshold of 3 since this is lowest value past 0.9 we start to see flattening 

# Scale-free topology fit index as a function of the soft-thresholding power
plot(Dermo_Susceptible_dds_vst_matrix_common_sft$fitIndices[,1], -sign(Dermo_Susceptible_dds_vst_matrix_common_sft$fitIndices[,3])*Dermo_Susceptible_dds_vst_matrix_common_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(Dermo_Susceptible_dds_vst_matrix_common_sft$fitIndices[,1], -sign(Dermo_Susceptible_dds_vst_matrix_common_sft$fitIndices[,3])*Dermo_Susceptible_dds_vst_matrix_common_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(Dermo_Susceptible_dds_vst_matrix_common_sft$fitIndices[,1], Dermo_Susceptible_dds_vst_matrix_common_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(Dermo_Susceptible_dds_vst_matrix_common_sft$fitIndices[,1], Dermo_Susceptible_dds_vst_matrix_common_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Softthreshold of 2 

# Probiotic, Probiotic_counts_apop_dds_rlog_matrix_sft: Does not fit scale free topology!
# Scale-free topology fit index as a function of the soft-thresholding power
plot(Probiotic_dds_rlog_matrix_common_sft$fitIndices[,1], -sign(Probiotic_dds_rlog_matrix_common_sft$fitIndices[,3])*Probiotic_dds_rlog_matrix_common_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(Probiotic_dds_rlog_matrix_common_sft$fitIndices[,1], -sign(Probiotic_dds_rlog_matrix_common_sft$fitIndices[,3])*Probiotic_dds_rlog_matrix_common_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(Probiotic_dds_rlog_matrix_common_sft$fitIndices[,1],Probiotic_dds_rlog_matrix_common_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(Probiotic_dds_rlog_matrix_common_sft$fitIndices[,1], Probiotic_dds_rlog_matrix_common_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Poor fit of scale free topology: see recommendations here: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html
# This is likely caused by strong clustering by Time

# If the scale-free topology fit index fails to reach values above 0.8 for reasonable powers
# (less than 15 for unsigned or signed hybrid networks, and less than 30 for signed networks) 
# and the mean connectivity remains relatively high (in the hundreds or above), 
# chances are that the data exhibit a strong driver that makes a subset of the 
# samples globally different from the rest.The difference causes high correlation among large 
# groups of genes which invalidates the assumption of the scale-free topology approximation.
# If the lack of scale-free topology fit turns out to be caused by an interesting biological variable 
# that one does not want to remove (i.e., adjust the data for), the appropriate soft-thresholding power 
# can be chosen based on the number of samples as in the table below. 

# Because less than 20 samples, selecting soft threshold of 9 for signed hybrid network 

# ROD pick soft threshold
# Scale-free topology fit index as a function of the soft-thresholding power
plot(ROD_Resistant_dds_rlog_matrix_common_sft $fitIndices[,1], -sign(ROD_Resistant_dds_rlog_matrix_common_sft $fitIndices[,3])*ROD_Resistant_dds_rlog_matrix_common_sft $fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(ROD_Resistant_dds_rlog_matrix_common_sft $fitIndices[,1], -sign(ROD_Resistant_dds_rlog_matrix_common_sft $fitIndices[,3])*ROD_Resistant_dds_rlog_matrix_common_sft $fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(ROD_Resistant_dds_rlog_matrix_common_sft $fitIndices[,1],ROD_Resistant_dds_rlog_matrix_common_sft $fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(ROD_Resistant_dds_rlog_matrix_common_sft $fitIndices[,1], ROD_Resistant_dds_rlog_matrix_common_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Same story as above, selecting soft threshold of 9 because less than 20 samples

# Scale-free topology fit index as a function of the soft-thresholding power
plot(  ROD_Susceptible_dds_rlog_matrix_common_sft$fitIndices[,1], -sign(  ROD_Susceptible_dds_rlog_matrix_common_sft$fitIndices[,3])*  ROD_Susceptible_dds_rlog_matrix_common_sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
text( ROD_Susceptible_dds_rlog_matrix_common_sft$fitIndices[,1], -sign(ROD_Susceptible_dds_rlog_matrix_common_sft$fitIndices[,3])*ROD_Susceptible_dds_rlog_matrix_common_sft$fitIndices[,2],
      labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(ROD_Susceptible_dds_rlog_matrix_common_sft$fitIndices[,1],ROD_Susceptible_dds_rlog_matrix_common_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(ROD_Susceptible_dds_rlog_matrix_common_sft$fitIndices[,1], ROD_Susceptible_dds_rlog_matrix_common_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Pro_RE22
# Scale-free topology fit index as a function of the soft-thresholding power
plot(Pro_RE22_dds_rlog_matrix_common_sft$fitIndices[,1], -sign(Pro_RE22_dds_rlog_matrix_common_sft$fitIndices[,3])*Pro_RE22_dds_rlog_matrix_common_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(Pro_RE22_dds_rlog_matrix_common_sft$fitIndices[,1], -sign(Pro_RE22_dds_rlog_matrix_common_sft$fitIndices[,3])*Pro_RE22_dds_rlog_matrix_common_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(Pro_RE22_dds_rlog_matrix_common_sft$fitIndices[,1],Pro_RE22_dds_rlog_matrix_common_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(Pro_RE22_dds_rlog_matrix_common_sft$fitIndices[,1], Pro_RE22_dds_rlog_matrix_common_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Set soft thresholding to 8

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Pro_RE22_Pro
# Scale-free topology fit index as a function of the soft-thresholding power
plot(Pro_RE22_dds_rlog_matrix_common_Pro_sft$fitIndices[,1], -sign(Pro_RE22_dds_rlog_matrix_common_Pro_sft$fitIndices[,3])*Pro_RE22_dds_rlog_matrix_common_Pro_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(Pro_RE22_dds_rlog_matrix_common_Pro_sft$fitIndices[,1], -sign(Pro_RE22_dds_rlog_matrix_common_Pro_sft$fitIndices[,3])*Pro_RE22_dds_rlog_matrix_common_Pro_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(Pro_RE22_dds_rlog_matrix_common_Pro_sft$fitIndices[,1],Pro_RE22_dds_rlog_matrix_common_Pro_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(Pro_RE22_dds_rlog_matrix_common_Pro_sft$fitIndices[,1], Pro_RE22_dds_rlog_matrix_common_Pro_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# deos not meet network free topology, using soft thresholding power of 9 due to low sample size

# Pro_RE22_RE22
# Scale-free topology fit index as a function of the soft-thresholding power
plot(Pro_RE22_dds_rlog_matrix_common_RE22_sft$fitIndices[,1], -sign(Pro_RE22_dds_rlog_matrix_common_RE22_sft$fitIndices[,3])*Pro_RE22_dds_rlog_matrix_common_RE22_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(Pro_RE22_dds_rlog_matrix_common_RE22_sft$fitIndices[,1], -sign(Pro_RE22_dds_rlog_matrix_common_RE22_sft$fitIndices[,3])*Pro_RE22_dds_rlog_matrix_common_RE22_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(Pro_RE22_dds_rlog_matrix_common_RE22_sft$fitIndices[,1],Pro_RE22_dds_rlog_matrix_common_RE22_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(Pro_RE22_dds_rlog_matrix_common_RE22_sft$fitIndices[,1], Pro_RE22_dds_rlog_matrix_common_RE22_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Set soft thresholding to 9 

#### ONE STEP NETWORK CONSTRUCTION, MODULE DETECTION, MODULE DENDROGRAM INSPECTION #### 
# Dermo_Tol_net = blockwiseModules(Dermo_Tolerant_dds_vst_matrix_common, power = 3, # picked suitable power in the code above 
#                                  TOMType = "signed", # use signed TOM type
#                                  networkType= "signed hybrid", # use signed hybrid network type
#                                  corType = "bicor", # use suggested bicor
#                                  TminModuleSize = 30, # recommended default
#                                  reassignThreshold = 0, # recommended default
#                                  mergeCutHeight = 0.25, # recommended default
#                                  numericLabels = TRUE, # recommended default
#                                  pamRespectsDendro = FALSE,# recommended default
#                                  verbose = 3, 
#                                  maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
# save(Dermo_Tol_net)
# # How many modules identified
# table(Dermo_Tol_net$colors)
# # Plot dendrogram with colors
# # open a graphics window
# sizeGrWindow(12, 9)
# # Convert labels to colors for plotting
# Dermo_Tol_mergedColors = labels2colors(Dermo_Tol_net$colors)
# # Plot the dendrogram and the module colors underneath
# plotDendroAndColors(Dermo_Tol_net$dendrograms[[1]], Dermo_Tol_mergedColors[Dermo_Tol_net$blockGenes[[1]]],
#                     "Module colors",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# Dermo_Tol_moduleLabels = Dermo_Tol_net$colors
# Dermo_Tol_moduleColors = labels2colors(Dermo_Tol_net$colors)
# Dermo_Tol_MEs = Dermo_Tol_net$MEs
# Dermo_Tol_geneTree = Dermo_Tol_net$dendrograms[[1]]
# save network
#save(Dermo_Tol_net, Dermo_Tol_mergedColors, Dermo_Tol_moduleLabels, Dermo_Tol_moduleColors, Dermo_Tol_MEs, Dermo_Tol_geneTree, 
#     file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Tol_network.RData")
# Load network
Dermo_Tol_network = load(file= "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Tol_network.RData")

# Dermo Susceptible repeat 
# Dermo_Sus_net = blockwiseModules(Dermo_Susceptible_dds_vst_matrix_common, power = 2, # picked suitable power in the code above 
#                                  TOMType = "signed", # use signed TOM type
#                                  networkType= "signed hybrid", # use signed hybrid network type
#                                  corType = "bicor", # use suggested bicor
#                                  TminModuleSize = 30, # recommended default
#                                  reassignThreshold = 0, # recommended default
#                                  mergeCutHeight = 0.25, # recommended default
#                                  numericLabels = TRUE, # recommended default
#                                  pamRespectsDendro = FALSE,# recommended default
#                                  verbose = 3, 
#                                  maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
# # How many modules identified
# table(Dermo_Sus_net$colors)
# # Plot dendrogram with colors
# # open a graphics window
# sizeGrWindow(12, 9)
# # Convert labels to colors for plotting
# Dermo_Sus_mergedColors = labels2colors(Dermo_Sus_net$colors)
# # Plot the dendrogram and the module colors underneath
# plotDendroAndColors(Dermo_Sus_net$dendrograms[[1]], Dermo_Sus_mergedColors[Dermo_Sus_net$blockGenes[[1]]],
#                     "Module colors",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# Dermo_Sus_moduleLabels = Dermo_Sus_net$colors
# Dermo_Sus_moduleColors = labels2colors(Dermo_Sus_net$colors)
# Dermo_Sus_MEs = Dermo_Sus_net$MEs
# Dermo_Sus_geneTree = Dermo_Sus_net$dendrograms[[1]]
# save(Dermo_Sus_net, Dermo_Sus_mergedColors, Dermo_Sus_moduleLabels, Dermo_Sus_moduleColors, Dermo_Sus_MEs, Dermo_Sus_geneTree, 
#      file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Sus_network.RData")

# Load network
Dermo_Sus_network = load(file= "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Sus_network.RData")

# Probiotic
#Probiotic_net = blockwiseModules(Probiotic_dds_rlog_matrix_common, power = 9, # picked because less than 20 samples and isn't scale free
#                                 TOMType = "signed", # use signed TOM type
#                                 networkType= "signed hybrid", # use signed hybrid network type
#                                 corType = "bicor", # use suggested bicor
#                                 TminModuleSize = 30, # recommended default
#                                 reassignThreshold = 0, # recommended default
#                                 mergeCutHeight = 0.25, # recommended default
#                                 numericLabels = TRUE, # recommended default
#                                 pamRespectsDendro = FALSE,# recommended default
#                                 verbose = 3, 
#                                 maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
## How many modules identified
#table(Probiotic_net$colors)
## Plot dendrogram with colors
## open a graphics window
#sizeGrWindow(12, 9)
## Convert labels to colors for plotting
#Probiotic_mergedColors = labels2colors(Probiotic_net$colors)
## Plot the dendrogram and the module colors underneath
#plotDendroAndColors(Probiotic_net$dendrograms[[1]], Probiotic_mergedColors[Probiotic_net$blockGenes[[1]]],
#                    "Module colors",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05)
#Probiotic_moduleLabels = Probiotic_net$colors
#Probiotic_moduleColors = labels2colors(Probiotic_net$colors)
#Probiotic_MEs = Probiotic_net$MEs
#Probiotic_geneTree = Probiotic_net$dendrograms[[1]]
## save network externally
#save(Probiotic_net, Probiotic_mergedColors, Probiotic_moduleLabels, Probiotic_moduleColors, Probiotic_MEs, Probiotic_geneTree, file=
#       "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Probiotic_network.RData")

# Load network
Probiotic_network = load(file= "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Probiotic_network.RData")

# ROD Res
#ROD_Res_net = blockwiseModules(ROD_Resistant_dds_rlog_matrix_common, power = 9, # picked because less than 20 samples and isn't scale free
#                               TOMType = "signed", # use signed TOM type
#                               networkType= "signed hybrid", # use signed hybrid network type
#                               corType = "bicor", # use suggested bicor
#                               TminModuleSize = 30, # recommended default
#                               reassignThreshold = 0, # recommended default
#                               mergeCutHeight = 0.25, # recommended default
#                               numericLabels = TRUE, # recommended default
#                               pamRespectsDendro = FALSE,# recommended default
#                               verbose = 3, 
#                               maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
## How many modules identified
#table(ROD_Res_net$colors)
## Plot dendrogram with colors
## open a graphics window
#sizeGrWindow(12, 9)
## Convert labels to colors for plotting
#ROD_Res_mergedColors = labels2colors(ROD_Res_net$colors)
## Plot the dendrogram and the module colors underneath
#plotDendroAndColors(ROD_Res_net$dendrograms[[1]], ROD_Res_mergedColors[ROD_Res_net$blockGenes[[1]]],
#                    "Module colors",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05)
#ROD_Res_moduleLabels = ROD_Res_net$colors
#ROD_Res_moduleColors = labels2colors(ROD_Res_net$colors)
#ROD_Res_MEs = ROD_Res_net$MEs
#ROD_Res_geneTree = ROD_Res_net$dendrograms[[1]]
#save(ROD_Res_net, ROD_Res_mergedColors, ROD_Res_moduleLabels, ROD_Res_moduleColors, ROD_Res_MEs, ROD_Res_geneTree, file=
#       "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/ROD_Res_network.RData")

# Load network
ROD_Res_network = load(file= "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/ROD_Res_network.RData")

# ROD Sus
#ROD_Sus_net = blockwiseModules(ROD_Susceptible_dds_rlog_matrix_common, power = 9, # picked because less than 20 samples and isn't scale free
#                               TOMType = "signed", # use signed TOM type
#                               networkType= "signed hybrid", # use signed hybrid network type
#                               corType = "bicor", # use suggested bicor
#                               TminModuleSize = 30, # recommended default
#                               reassignThreshold = 0, # recommended default
#                               mergeCutHeight = 0.25, # recommended default
#                               numericLabels = TRUE, # recommended default
#                               pamRespectsDendro = FALSE,# recommended default
#                               verbose = 3, 
#                               maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
## How many modules identified
#table(ROD_Sus_net$colors)
## Plot dendrogram with colors
## open a graphics window
#sizeGrWindow(12, 9)
## Convert labels to colors for plotting
#ROD_Sus_mergedColors = labels2colors(ROD_Sus_net$colors)
## Plot the dendrogram and the module colors underneath
#plotDendroAndColors(ROD_Sus_net$dendrograms[[1]], ROD_Sus_mergedColors[ROD_Sus_net$blockGenes[[1]]],
#                    "Module colors",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05)
#ROD_Sus_moduleLabels = ROD_Sus_net$colors
#ROD_Sus_moduleColors = labels2colors(ROD_Sus_net$colors)
#ROD_Sus_MEs = ROD_Sus_net$MEs
#ROD_Sus_geneTree = ROD_Sus_net$dendrograms[[1]]
#
#save(ROD_Sus_net, ROD_Sus_mergedColors, ROD_Sus_moduleLabels, ROD_Sus_moduleColors, ROD_Sus_MEs, ROD_Sus_geneTree, file=
#       "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/ROD_Sus_network.RData")

# Load network
ROD_Sus_network = load(file= "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/ROD_Sus_network.RData")

# Pro_RE22
#Pro_RE22_net = blockwiseModules(Pro_RE22_dds_rlog_matrix_common, power = 8, 
#                                TOMType = "signed", # use signed TOM type
#                                networkType= "signed hybrid", # use signed hybrid network type
#                                corType = "bicor", # use suggested bicor
#                                TminModuleSize = 30, # recommended default
#                                reassignThreshold = 0, # recommended default
#                                mergeCutHeight = 0.25, # recommended default
#                                numericLabels = TRUE, # recommended default
#                                pamRespectsDendro = FALSE,# recommended default
#                                verbose = 3, 
#                                maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
## How many modules identified
#table(Pro_RE22_net$colors)
## Plot dendrogram with colors
## open a graphics window
#sizeGrWindow(12, 9)
## Convert labels to colors for plotting
#Pro_RE22_mergedColors = labels2colors(Pro_RE22_net$colors)
## Plot the dendrogram and the module colors underneath
#plotDendroAndColors(Pro_RE22_net$dendrograms[[1]], Pro_RE22_mergedColors[Pro_RE22_net$blockGenes[[1]]],
#                    "Module colors",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05)
#Pro_RE22_moduleLabels = Pro_RE22_net$colors
#Pro_RE22_moduleColors = labels2colors(Pro_RE22_net$colors)
#Pro_RE22_MEs = Pro_RE22_net$MEs
#Pro_RE22_geneTree = Pro_RE22_net$dendrograms[[1]]
#save(Pro_RE22_net, Pro_RE22_mergedColors, Pro_RE22_moduleLabels, Pro_RE22_moduleColors, Pro_RE22_MEs, Pro_RE22_geneTree, file=
#       "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_network.RData")

# Load network
Pro_RE22_network = load(file= "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_network.RData")

# Pro_RE22_Pro
Pro_RE22_Pro_net = blockwiseModules(Pro_RE22_dds_rlog_matrix_common_Pro, power = 9, 
                                TOMType = "signed", # use signed TOM type
                                networkType= "signed hybrid", # use signed hybrid network type
                                corType = "bicor", # use suggested bicor
                                TminModuleSize = 30, # recommended default
                                reassignThreshold = 0, # recommended default
                                mergeCutHeight = 0.25, # recommended default
                                numericLabels = TRUE, # recommended default
                                pamRespectsDendro = FALSE,# recommended default
                                verbose = 3, 
                                maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
# How many modules identified
table(Pro_RE22_Pro_net$colors)
# Plot dendrogram with colors
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
Pro_RE22_Pro_mergedColors = labels2colors(Pro_RE22_Pro_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(Pro_RE22_Pro_net$dendrograms[[1]], Pro_RE22_Pro_mergedColors[Pro_RE22_Pro_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
Pro_RE22_Pro_moduleLabels = Pro_RE22_Pro_net$colors
Pro_RE22_Pro_moduleColors = labels2colors(Pro_RE22_Pro_net$colors)
Pro_RE22_Pro_MEs = Pro_RE22_Pro_net$MEs
Pro_RE22_Pro_geneTree = Pro_RE22_Pro_net$dendrograms[[1]]

save(Pro_RE22_Pro_net, Pro_RE22_Pro_mergedColors, Pro_RE22_Pro_moduleLabels, Pro_RE22_Pro_moduleColors, Pro_RE22_Pro_MEs, Pro_RE22_Pro_geneTree, file=
       "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_Pro_network.RData")

# Pro_RE22_RE22
Pro_RE22_RE22_net = blockwiseModules(Pro_RE22_dds_rlog_matrix_common_RE22, power = 9, 
                                    TOMType = "signed", # use signed TOM type
                                    networkType= "signed hybrid", # use signed hybrid network type
                                    corType = "bicor", # use suggested bicor
                                    TminModuleSize = 30, # recommended default
                                    reassignThreshold = 0, # recommended default
                                    mergeCutHeight = 0.25, # recommended default
                                    numericLabels = TRUE, # recommended default
                                    pamRespectsDendro = FALSE,# recommended default
                                    verbose = 3, 
                                    maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
# How many modules identified
table(Pro_RE22_RE22_net$colors)
# Plot dendrogram with colors
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
Pro_RE22_RE22_mergedColors = labels2colors(Pro_RE22_RE22_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(Pro_RE22_RE22_net$dendrograms[[1]], Pro_RE22_RE22_mergedColors[Pro_RE22_RE22_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
Pro_RE22_RE22_moduleLabels = Pro_RE22_RE22_net$colors
Pro_RE22_RE22_moduleColors = labels2colors(Pro_RE22_RE22_net$colors)
Pro_RE22_RE22_MEs = Pro_RE22_RE22_net$MEs
Pro_RE22_RE22_geneTree = Pro_RE22_RE22_net$dendrograms[[1]]

save(Pro_RE22_RE22_net, Pro_RE22_RE22_mergedColors, Pro_RE22_RE22_moduleLabels, Pro_RE22_RE22_moduleColors, Pro_RE22_RE22_MEs, Pro_RE22_RE22_geneTree, file=
       "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_RE22_network.RData")

#### BINARIZE ALL CATEGORICAL VARIABLES TO TEST DISEASE CHALLENGE ASSOCIATIONS ####


## Dermo
Dermo_Tolerant_coldata_collapsed <- Dermo_Tolerant_coldata[,c("Sample_ID", "Condition","Time")]
Dermo_Tolerant_coldata_collapsed <- Dermo_Tolerant_coldata_collapsed[!duplicated(Dermo_Tolerant_coldata_collapsed$Sample_ID),]
row.names(Dermo_Tolerant_coldata_collapsed ) <- Dermo_Tolerant_coldata_collapsed$Sample_ID
Dermo_Tolerant_coldata_collapsed <- Dermo_Tolerant_coldata_collapsed[,c("Condition","Time")]
nrow(Dermo_Tolerant_coldata_collapsed) # 30
length(row.names(Dermo_Tolerant_dds_vst_matrix_common)) # 30

# check order
all(row.names(Dermo_Tolerant_coldata_collapsed) %in% row.names(Dermo_Tolerant_dds_vst_matrix_common)) # TRUE
all(row.names(Dermo_Tolerant_dds_vst_matrix_common) %in% row.names(Dermo_Tolerant_coldata_collapsed) ) # TRUE
all(row.names(Dermo_Tolerant_coldata_collapsed) == row.names(Dermo_Tolerant_dds_vst_matrix_common)) # FALSE
all(row.names(Dermo_Tolerant_dds_vst_matrix_common) == row.names(Dermo_Tolerant_coldata_collapsed) ) # FALSE

# fix order
Dermo_Tolerant_coldata_collapsed<- Dermo_Tolerant_coldata_collapsed[row.names(Dermo_Tolerant_dds_vst_matrix_common),]
all(row.names(Dermo_Tolerant_coldata_collapsed) == row.names(Dermo_Tolerant_dds_vst_matrix_common)) # TRUE
all(row.names(Dermo_Tolerant_dds_vst_matrix_common) == row.names(Dermo_Tolerant_coldata_collapsed) ) # TRUE

# Dermo sus
Dermo_Susceptible_coldata_collapsed<- Dermo_Susceptible_coldata[,c("Sample_ID", "Condition","Time")]
Dermo_Susceptible_coldata_collapsed <- Dermo_Susceptible_coldata_collapsed[!duplicated(Dermo_Susceptible_coldata_collapsed$Sample_ID),]
row.names(Dermo_Susceptible_coldata_collapsed ) <- Dermo_Susceptible_coldata_collapsed$Sample_ID
Dermo_Susceptible_coldata_collapsed <- Dermo_Susceptible_coldata_collapsed[,c("Condition","Time")]
nrow(Dermo_Susceptible_coldata_collapsed) # 32
length(row.names(Dermo_Susceptible_dds_vst_matrix_common)) # 32

all(row.names(Dermo_Susceptible_coldata_collapsed) %in% row.names(Dermo_Susceptible_dds_vst_matrix_common)) # TRUE
all(row.names(Dermo_Susceptible_dds_vst_matrix_common) %in% row.names(Dermo_Susceptible_coldata_collapsed) )  # TRUE
all(row.names(Dermo_Susceptible_coldata_collapsed) == row.names(Dermo_Susceptible_dds_vst_matrix_common)) # FALSE
all(row.names(Dermo_Susceptible_dds_vst_matrix_common) == row.names(Dermo_Susceptible_coldata_collapsed) )  # FALSE

# Fix order 
Dermo_Susceptible_coldata_collapsed <- Dermo_Susceptible_coldata_collapsed[row.names(Dermo_Susceptible_dds_vst_matrix_common),]
all(row.names(Dermo_Susceptible_coldata_collapsed) == row.names(Dermo_Susceptible_dds_vst_matrix_common)) # TRUE
all(row.names(Dermo_Susceptible_dds_vst_matrix_common) == row.names(Dermo_Susceptible_coldata_collapsed) ) # TRUE

# Binarize the data table
Dermo_Tolerant_coldata_collapsed_binarize <- binarizeCategoricalColumns.pairwise(Dermo_Tolerant_coldata_collapsed)
Dermo_Susceptible_coldata_collapsed_binarize <- binarizeCategoricalColumns.pairwise(Dermo_Susceptible_coldata_collapsed)
row.names(Dermo_Tolerant_coldata_collapsed_binarize) <- row.names(Dermo_Tolerant_coldata_collapsed)
row.names(Dermo_Susceptible_coldata_collapsed_binarize) <- row.names(Dermo_Susceptible_coldata_collapsed)

## Probiotic
#Probiotic_dds_rlog_matrix_common 
#Probiotic_coldata
Probiotic_coldata_collapsed <- Probiotic_coldata[,c("Condition","Time")]
length(row.names(Probiotic_coldata_collapsed )) # 6

# check order
all(row.names(Probiotic_coldata_collapsed) %in% row.names(Probiotic_dds_rlog_matrix_common )) # TRUE
all(row.names(Probiotic_dds_rlog_matrix_common) %in% row.names(Probiotic_coldata_collapsed) ) # TRUE
all(row.names(Probiotic_coldata_collapsed) == row.names(Probiotic_dds_rlog_matrix_common )) # TRUE
all(row.names(Probiotic_dds_rlog_matrix_common) == row.names(Probiotic_coldata_collapsed) ) # TRUE

Probiotic_coldata_collapsed_binarize <- binarizeCategoricalColumns.pairwise(Probiotic_coldata_collapsed)
row.names(Probiotic_coldata_collapsed_binarize) <- row.names(Probiotic_coldata_collapsed)

## ROD RES
ROD_Resistant_coldata
ROD_Susceptible_coldata
ROD_Resistant_dds_rlog_matrix_common
ROD_Susceptible_dds_rlog_matrix_common

ROD_Resistant_coldata_collapsed <-ROD_Resistant_coldata[,c("Condition","Time")]
nrow(ROD_Resistant_coldata_collapsed) # 8
length(row.names(ROD_Resistant_dds_rlog_matrix_common)) # 8

# check order
all(row.names(ROD_Resistant_coldata_collapsed) %in% row.names(ROD_Resistant_dds_rlog_matrix_common)) # TRUE
all(row.names(ROD_Resistant_dds_rlog_matrix_common) %in% row.names(ROD_Resistant_coldata_collapsed) ) # TRUE
all(row.names(ROD_Resistant_coldata_collapsed) == row.names(ROD_Resistant_dds_rlog_matrix_common)) # TRUE
all(row.names(ROD_Resistant_dds_rlog_matrix_common) == row.names(ROD_Resistant_coldata_collapsed) ) # TRUE

ROD_Susceptible_coldata_collapsed <-ROD_Susceptible_coldata[,c("Condition","Time")]
nrow(ROD_Susceptible_coldata_collapsed) # 4
length(row.names(ROD_Susceptible_dds_rlog_matrix_common)) # 4

# check order
all(row.names(ROD_Susceptible_coldata_collapsed) %in% row.names(ROD_Susceptible_dds_rlog_matrix_common)) # TRUE
all(row.names(ROD_Susceptible_dds_rlog_matrix_common) %in% row.names(ROD_Susceptible_coldata_collapsed) ) # TRUE
all(row.names(ROD_Susceptible_coldata_collapsed) == row.names(ROD_Susceptible_dds_rlog_matrix_common)) # TRUE
all(row.names(ROD_Susceptible_dds_rlog_matrix_common) == row.names(ROD_Susceptible_coldata_collapsed) ) # TRUE

# Binarize the data table
ROD_Resistant_coldata_collapsed_binarize <- binarizeCategoricalColumns.pairwise(ROD_Resistant_coldata_collapsed)
ROD_Susceptible_coldata_collapsed_binarize <- binarizeCategoricalColumns.pairwise(ROD_Susceptible_coldata_collapsed)
row.names(ROD_Resistant_coldata_collapsed_binarize) <- row.names(ROD_Resistant_coldata_collapsed)
row.names(ROD_Susceptible_coldata_collapsed_binarize) <- row.names(ROD_Susceptible_coldata_collapsed)

## Pro_RE22
Pro_RE22_coldata_collapsed <- Pro_RE22_coldata[,c("Condition","Time")]
length(row.names(Pro_RE22_coldata_collapsed )) # 18

# Gsub to just look at RE22, vs. RI0695 at any time, and S4 at any time
Bacillus_pumilus_RI06_95_exposure_24h
Bacillus_pumilus_RI06_95_exposure_6h
Phaeobacter_inhibens_S4_exposure_6h
Phaeobacter_inhibens_S4_exposure_24h

Pro_RE22_coldata_collapsed$Condition <- gsub("Bacillus_pumilus_RI06_95_exposure_24h", "Bacillus_pumilus",Pro_RE22_coldata_collapsed$Condition)
Pro_RE22_coldata_collapsed$Condition <- gsub("Bacillus_pumilus_RI06_95_exposure_6h", "Bacillus_pumilus",Pro_RE22_coldata_collapsed$Condition)
Pro_RE22_coldata_collapsed$Condition <- gsub("Phaeobacter_inhibens_S4_exposure_6h", "Phaeobacter_inhibens",Pro_RE22_coldata_collapsed$Condition)
Pro_RE22_coldata_collapsed$Condition <- gsub("Phaeobacter_inhibens_S4_exposure_24h", "Phaeobacter_inhibens",Pro_RE22_coldata_collapsed$Condition)

# check order
all(row.names(Pro_RE22_coldata_collapsed ) %in% row.names(Pro_RE22_dds_rlog_matrix_common)) # TRUE
all(row.names(Pro_RE22_dds_rlog_matrix_common) %in% row.names(Pro_RE22_coldata_collapsed)) # TRUE
all(row.names(Pro_RE22_coldata_collapsed) == row.names(Pro_RE22_dds_rlog_matrix_common)) # TRUE
all(row.names(Pro_RE22_dds_rlog_matrix_common) == row.names(Pro_RE22_coldata_collapsed)) # TRUE

# binarize the data table
Pro_RE22_coldata_collapsed_binarize <- binarizeCategoricalColumns.pairwise(Pro_RE22_coldata_collapsed)
row.names(Pro_RE22_coldata_collapsed_binarize) <- row.names(Pro_RE22_coldata_collapsed)

# binarize the split Pro_RE22 experiment
Pro_RE22_coldata_RE22_collapsed <- Pro_RE22_coldata_RE22[,c("Condition","Time")] 
Pro_RE22_coldata_Pro_collapsed <- Pro_RE22_coldata_Pro[,c("Condition","Time")]
row.names(Pro_RE22_coldata_RE22_collapsed) <- Pro_RE22_coldata_RE22$sample
row.names(Pro_RE22_coldata_Pro_collapsed) <- Pro_RE22_coldata_Pro$sample
Pro_RE22_coldata_Pro_collapsed$Condition <- gsub("Bacillus_pumilus_RI06_95_exposure_24h", "Bacillus_pumilus",   Pro_RE22_coldata_Pro_collapsed$Condition)
Pro_RE22_coldata_Pro_collapsed$Condition <- gsub("Bacillus_pumilus_RI06_95_exposure_6h", "Bacillus_pumilus",    Pro_RE22_coldata_Pro_collapsed$Condition)
Pro_RE22_coldata_Pro_collapsed$Condition <- gsub("Phaeobacter_inhibens_S4_exposure_6h", "Phaeobacter_inhibens", Pro_RE22_coldata_Pro_collapsed$Condition)
Pro_RE22_coldata_Pro_collapsed$Condition <- gsub("Phaeobacter_inhibens_S4_exposure_24h", "Phaeobacter_inhibens",Pro_RE22_coldata_Pro_collapsed$Condition)

Pro_RE22_coldata_RE22_collapsed_binarize <- binarizeCategoricalColumns.pairwise(Pro_RE22_coldata_RE22_collapsed)
row.names(Pro_RE22_coldata_RE22_collapsed_binarize) <- row.names(Pro_RE22_coldata_RE22_collapsed)
Pro_RE22_coldata_Pro_collapsed_binarize <- binarizeCategoricalColumns.pairwise(Pro_RE22_coldata_Pro_collapsed)
row.names(Pro_RE22_coldata_Pro_collapsed_binarize) <- row.names(Pro_RE22_coldata_Pro_collapsed)

#### QUANTIFYING MODULE TRAIT ASSOCIATIONS ####

#### DERMO TOL ####
# Define numbers of genes and samples
Dermo_Tol_nGenes = ncol(Dermo_Tolerant_dds_vst_matrix_common)
Dermo_Tol_nSamples = nrow(Dermo_Tolerant_dds_vst_matrix_common)

# Recalculate MEs with color labels
Dermo_Tol_MEs0 = moduleEigengenes(Dermo_Tolerant_dds_vst_matrix_common, Dermo_Tol_moduleColors)$eigengenes
Dermo_Tol_MEs = orderMEs(Dermo_Tol_MEs0)
Dermo_Tol_moduleTraitCor = cor(Dermo_Tol_MEs, Dermo_Tolerant_coldata_collapsed_binarize, use = "p");
Dermo_Tol_moduleTraitPvalue = corPvalueStudent(Dermo_Tol_moduleTraitCor, Dermo_Tol_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trai
sizeGrWindow(10,6)
# Will display correlations and their p-values
Dermo_Tol_textMatrix = paste(signif(Dermo_Tol_moduleTraitCor, 2), "\n(",
                   signif(Dermo_Tol_moduleTraitPvalue, 1), ")", sep = "");
dim(Dermo_Tol_textMatrix) = dim(Dermo_Tol_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
  # green is more negatively correlated)
labeledHeatmap(Matrix = Dermo_Tol_moduleTraitCor,
               xLabels = names(Dermo_Tolerant_coldata_collapsed_binarize),
               yLabels = names(Dermo_Tol_MEs),
               ySymbols = names(Dermo_Tol_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = Dermo_Tol_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
Dermo_Tol_moduleTraitCor_df <- as.data.frame(Dermo_Tol_moduleTraitCor)
Dermo_Tol_moduleTraitCor_df$mod_names <- row.names(Dermo_Tol_moduleTraitCor_df)
Dermo_Tol_moduleTraitCor_df <- Dermo_Tol_moduleTraitCor_df[,c("mod_names","Condition.Injected.vs.Control")]
Dermo_Tol_moduleTraitPvalue_df <- as.data.frame(Dermo_Tol_moduleTraitPvalue)
Dermo_Tol_moduleTraitPvalue_df$mod_names <- row.names(Dermo_Tol_moduleTraitPvalue_df)
Dermo_Tol_moduleTraitPvalue_df <- Dermo_Tol_moduleTraitPvalue_df[,c("mod_names","Condition.Injected.vs.Control")]
colnames(Dermo_Tol_moduleTraitPvalue_df)[2] <- "pvalue"

Dermo_Tol_moduleTraitCor_Pval_df <- join(Dermo_Tol_moduleTraitCor_df, Dermo_Tol_moduleTraitPvalue_df, by = "mod_names")

# Significantly correlated modules
Dermo_Tol_moduleTraitCor_Pval_df[order(Dermo_Tol_moduleTraitCor_Pval_df$pvalue),]
class(Dermo_Tol_moduleTraitCor_Pval_df$pvalue) # numeric
Dermo_Tol_moduleTraitCor_Pval_df_sig <- Dermo_Tol_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
Dermo_Tol_moduleTraitCor_Pval_df_sig

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
Dermo_Tol_moduleTraitCor_Pval_df_sig_list <- Dermo_Tol_moduleTraitCor_Pval_df_sig$mod_names
Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Dermo_Tol_moduleTraitCor_Pval_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
matrix_common= Dermo_Tolerant_dds_vst_matrix_common
moduleColors=Dermo_Tol_moduleColors
lookup =   C_vir_rtracklayer_apop_product_final

lookup_mod_apop <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$ID %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
names(Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm) <- c( "thistle1"  ,     "lightgreen"  ,   "pink"        ,   "royalblue"  ,    "honeydew1"  ,    "darkseagreen3",  "tan"     ,       "darkviolet" ,
                                                          "orangered4"  ,   "yellow"   ,      "antiquewhite2" , "lightcyan1"    ,
                                                          "plum2"     ,     "turquoise"   ,   "lightsteelblue")
Dermo_Tol_module_apop <- lapply(Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop)
Dermo_Tol_module_apop_df <- do.call(rbind,Dermo_Tol_module_apop)
Dermo_Tol_module_apop_df$mod_names <- gsub("\\..*","",row.names(Dermo_Tol_module_apop_df))
Dermo_Tol_module_apop_df$mod_names <- gsub("^","ME",Dermo_Tol_module_apop_df$mod_names)
# add module significance
Dermo_Tol_module_apop_df <- left_join(Dermo_Tol_module_apop_df,Dermo_Tol_moduleTraitCor_Pval_df_sig)
Dermo_Tol_module_apop_df$exp <- "Dermo_Tol"

## Gene relationship to trait and important modules: Gene Significance and Module Membership
# We quantify associations of individual genes with our trait of interest (weight) by defining Gene
# Significance GS as (the absolute value of) the correlation between the gene and the trait. 
# For each module, we also define a quantitative measure of module membership MM as the correlation
# of the module eigengene and the gene expression profile. This allows us to quantify 
# the similarity of all genes on the array to every module.
# higher the absolute value of GS, the more biologically relevant

# Define variable injected 
Dermo_Tol_injection = as.data.frame(Dermo_Tolerant_coldata_collapsed_binarize$Condition.Injected.vs.Control);
names(Dermo_Tol_injection) = "injection"
# names (colors) of the modules
Dermo_Tol_modNames = substring(names(Dermo_Tol_MEs), 3)
Dermo_Tol_geneModuleMembership = as.data.frame(cor(Dermo_Tolerant_dds_vst_matrix_common, Dermo_Tol_MEs, use = "p"))
Dermo_Tol_MMPvalue = as.data.frame(corPvalueStudent(as.matrix(Dermo_Tol_geneModuleMembership), Dermo_Tol_nSamples))

names(Dermo_Tol_geneModuleMembership) = paste("MM", Dermo_Tol_modNames, sep="")
names(Dermo_Tol_MMPvalue) = paste("p.MM", Dermo_Tol_modNames, sep="")

Dermo_Tol_geneTraitSignificance = as.data.frame(cor(Dermo_Tolerant_dds_vst_matrix_common,Dermo_Tol_injection, use = "p"))
Dermo_Tol_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(Dermo_Tol_geneTraitSignificance), Dermo_Tol_nSamples))

names(Dermo_Tol_geneTraitSignificance) = paste("GS.", names(Dermo_Tol_injection), sep="")
names(Dermo_Tol_GSPvalue) = paste("p.GS.", names(Dermo_Tol_injection), sep="")

## Intramodular analysis: identifying genes with high GS and MM
# Using the GS and MM measures, we can identify genes that have a high significance for disease challenge
# as well as high module membership in interesting modules. As an example, we look at the brown module 
# that has the highest association with weight. We plot a scatterplot of Gene Significance vs. Module Membership in the brown module:
Dermo_Tol_module = "lightgreen" # strong correlation
Dermo_Tol_column = match(Dermo_Tol_module, Dermo_Tol_modNames)
Dermo_Tol_moduleGenes = Dermo_Tol_moduleColors==Dermo_Tol_module
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(Dermo_Tol_geneModuleMembership[Dermo_Tol_moduleGenes, Dermo_Tol_column]),
                   abs(Dermo_Tol_geneTraitSignificance[Dermo_Tol_moduleGenes, 1]),
                   xlab = paste("Module Membership in", Dermo_Tol_module, "module"),
                   ylab = "Gene significance for challenge",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = Dermo_Tol_module)

## IDENTIFY HUB GENES IN EACH SIG MODULE ##
Dermo_Tol_colorh = c("thistle1","lightgreen","pink","royalblue","honeydew1","darkseagreen3","tan","darkviolet",
"orangered4","yellow","antiquewhite2","lightcyan1","plum2","turquoise","lightsteelblue")

#Dermo_Tol_Module_hub_genes <- chooseTopHubInEachModule(
#  Dermo_Tolerant_dds_vst_matrix_common, 
#  Dermo_Tol_colorh, 
#  power = 3,  # power used for the adjacency network
#  type = "signed hybrid", 
#  corFnc = "bicor"
#  )
class(Dermo_Tol_Module_hub_genes)
Dermo_Tol_Module_hub_genes_df <- as.data.frame(Dermo_Tol_Module_hub_genes)
colnames(Dermo_Tol_Module_hub_genes_df)[1] <- "ID"
Dermo_Tol_Module_hub_genes_apop <- merge(Dermo_Tol_Module_hub_genes_df, C_vir_rtracklayer, by = "ID")
nrow(Dermo_Tol_Module_hub_genes_apop) # 15

#### DERMO SUS ####
# Define numbers of genes and samples
Dermo_Sus_nGenes = ncol(Dermo_Susceptible_dds_vst_matrix_common)
Dermo_Sus_nSamples = nrow(Dermo_Susceptible_dds_vst_matrix_common)

# Recalculate MEs with color labels
Dermo_Sus_MEs0 = moduleEigengenes(Dermo_Susceptible_dds_vst_matrix_common, Dermo_Sus_moduleColors)$eigengenes
Dermo_Sus_MEs = orderMEs(Dermo_Sus_MEs0)
Dermo_Sus_moduleTraitCor = cor(Dermo_Sus_MEs, Dermo_Susceptible_coldata_collapsed_binarize, use = "p");
Dermo_Sus_moduleTraitPvalue = corPvalueStudent(Dermo_Sus_moduleTraitCor, Dermo_Sus_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trai
sizeGrWindow(10,6)
# Will display correlations and their p-values
Dermo_Sus_textMatrix = paste(signif(Dermo_Sus_moduleTraitCor, 2), "\n(",
                             signif(Dermo_Sus_moduleTraitPvalue, 1), ")", sep = "");
dim(Dermo_Sus_textMatrix) = dim(Dermo_Sus_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = Dermo_Sus_moduleTraitCor,
               xLabels = names(Dermo_Susceptible_coldata_collapsed_binarize),
               yLabels = names(Dermo_Sus_MEs),
               ySymbols = names(Dermo_Sus_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = Dermo_Sus_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
Dermo_Sus_moduleTraitCor_df <- as.data.frame(Dermo_Sus_moduleTraitCor)
Dermo_Sus_moduleTraitCor_df$mod_names <- row.names(Dermo_Sus_moduleTraitCor_df)
Dermo_Sus_moduleTraitCor_df <- Dermo_Sus_moduleTraitCor_df[,c("mod_names","Condition.Injected.vs.Control")]
Dermo_Sus_moduleTraitPvalue_df <- as.data.frame(Dermo_Sus_moduleTraitPvalue)
Dermo_Sus_moduleTraitPvalue_df$mod_names <- row.names(Dermo_Sus_moduleTraitPvalue_df)
Dermo_Sus_moduleTraitPvalue_df <- Dermo_Sus_moduleTraitPvalue_df[,c("mod_names","Condition.Injected.vs.Control")]
colnames(Dermo_Sus_moduleTraitPvalue_df)[2] <- "pvalue"

Dermo_Sus_moduleTraitCor_Pval_df <- join(Dermo_Sus_moduleTraitCor_df, Dermo_Sus_moduleTraitPvalue_df, by = "mod_names")

# Significantly correlated modules
Dermo_Sus_moduleTraitCor_Pval_df[order(Dermo_Sus_moduleTraitCor_Pval_df$pvalue),]
class(Dermo_Sus_moduleTraitCor_Pval_df$pvalue) # numeric
Dermo_Sus_moduleTraitCor_Pval_df_sig <- Dermo_Sus_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
Dermo_Sus_moduleTraitCor_Pval_df_sig

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
Dermo_Sus_moduleTraitCor_Pval_df_sig_list <- Dermo_Sus_moduleTraitCor_Pval_df_sig$mod_names
#"MElightgreen"  "MEgreenyellow"
Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Dermo_Sus_moduleTraitCor_Pval_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
matrix_common= Dermo_Susceptible_dds_vst_matrix_common
moduleColors=Dermo_Sus_moduleColors
lookup =   C_vir_rtracklayer_apop_product_final

lookup_mod_apop <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$ID %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
names(Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm) <- c("lightgreen",  "greenyellow")
Dermo_Sus_module_apop <- lapply(Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop)
Dermo_Sus_module_apop_df <- do.call(rbind,Dermo_Sus_module_apop)
Dermo_Sus_module_apop_df$mod_names <- gsub("\\..*","",row.names(Dermo_Sus_module_apop_df))
Dermo_Sus_module_apop_df$mod_names <- gsub("^","ME",Dermo_Sus_module_apop_df$mod_names)
# add module significance
Dermo_Sus_module_apop_df <- left_join(Dermo_Sus_module_apop_df,Dermo_Sus_moduleTraitCor_Pval_df_sig)
Dermo_Sus_module_apop_df$exp <- "Dermo_Sus"

## Gene relationship to trait and important modules: Gene Significance and Module Membership
# We quantify associations of individual genes with our trait of interest (weight) by defining Gene
# Significance GS as (the absolute value of) the correlation between the gene and the trait. 
# For each module, we also define a quantitative measure of module membership MM as the correlation
# of the module eigengene and the gene expression profile. This allows us to quantify 
# the similarity of all genes on the array to every module.
# higher the absolute value of GS, the more biologically relevant

# Define variable injected 
Dermo_Sus_injection = as.data.frame(Dermo_Susceptible_coldata_collapsed_binarize$Condition.Injected.vs.Control);
names(Dermo_Sus_injection) = "injection"
# names (colors) of the modules
Dermo_Sus_modNames = substring(names(Dermo_Sus_MEs), 3)
Dermo_Sus_geneModuleMembership = as.data.frame(cor(Dermo_Susceptible_dds_vst_matrix_common, Dermo_Sus_MEs, use = "p"))
Dermo_Sus_MMPvalue = as.data.frame(corPvalueStudent(as.matrix(Dermo_Sus_geneModuleMembership), Dermo_Sus_nSamples))

names(Dermo_Sus_geneModuleMembership) = paste("MM", Dermo_Sus_modNames, sep="")
names(Dermo_Sus_MMPvalue) = paste("p.MM", Dermo_Sus_modNames, sep="")

Dermo_Sus_geneTraitSignificance = as.data.frame(cor(Dermo_Susceptible_dds_vst_matrix_common,Dermo_Sus_injection, use = "p"))
Dermo_Sus_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(Dermo_Sus_geneTraitSignificance), Dermo_Sus_nSamples))

names(Dermo_Sus_geneTraitSignificance) = paste("GS.", names(Dermo_Sus_injection), sep="")
names(Dermo_Sus_GSPvalue) = paste("p.GS.", names(Dermo_Sus_injection), sep="")

## Intramodular analysis: identifying genes with high GS and MM
# Using the GS and MM measures, we can identify genes that have a high significance for disease challenge
# as well as high module membership in interesting modules. As an example, we look at the brown module 
# that has the highest association with weight. We plot a scatterplot of Gene Significance vs. Module Membership in the brown module:
Dermo_Sus_module = "lightgreen" # strong correlation
Dermo_Sus_column = match(Dermo_Sus_module, Dermo_Sus_modNames)
Dermo_Sus_moduleGenes = Dermo_Sus_moduleColors==Dermo_Sus_module
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(Dermo_Sus_geneModuleMembership[Dermo_Sus_moduleGenes, Dermo_Sus_column]),
                   abs(Dermo_Sus_geneTraitSignificance[Dermo_Sus_moduleGenes, 1]),
                   xlab = paste("Module Membership in", Dermo_Sus_module, "module"),
                   ylab = "Gene significance for challenge",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = Dermo_Sus_module)

## IDENTIFY HUB GENES IN EACH SIG MODULE ##
Dermo_Sus_colorh = c("MElightgreen","MEgreenyellow")

#Dermo_Sus_Module_hub_genes <- chooseTopHubInEachModule(
#  Dermo_Susceptible_dds_vst_matrix_common, 
#  Dermo_Sus_colorh, 
#  power = 2,  # power used for the adjacency network
#  type = "signed hybrid", 
#  corFnc = "bicor"
#)
class(Dermo_Sus_Module_hub_genes)
Dermo_Sus_Module_hub_genes_df <- as.data.frame(Dermo_Sus_Module_hub_genes)
colnames(Dermo_Sus_Module_hub_genes_df)[1] <- "ID"
Dermo_Sus_Module_hub_genes_apop <- merge(Dermo_Sus_Module_hub_genes_df, C_vir_rtracklayer, by = "ID")
nrow(Dermo_Sus_Module_hub_genes_apop) # 2 uncharacterized loci

#### PROBIOTIC ####
Probiotic_coldata_collapsed_binarize 
Probiotic_dds_rlog_matrix_common

# Define numbers of genes and samples
Probiotic_nGenes = ncol(Probiotic_dds_rlog_matrix_common)
Probiotic_nSamples = nrow(Probiotic_dds_rlog_matrix_common)

# Recalculate MEs with color labels
Probiotic_MEs0 = moduleEigengenes(Probiotic_dds_rlog_matrix_common, Probiotic_moduleColors)$eigengenes
Probiotic_MEs = orderMEs(Probiotic_MEs0)
Probiotic_moduleTraitCor = cor(Probiotic_MEs, Probiotic_coldata_collapsed_binarize , use = "p");
Probiotic_moduleTraitPvalue = corPvalueStudent(Probiotic_moduleTraitCor, Probiotic_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trait
sizeGrWindow(10,6)
# Will display correlations and their p-values
Probiotic_textMatrix = paste(signif(Probiotic_moduleTraitCor, 2), "\n(",
                             signif(Probiotic_moduleTraitPvalue, 1), ")", sep = "");
dim(Probiotic_textMatrix) = dim(Probiotic_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = Probiotic_moduleTraitCor,
               xLabels = names(Probiotic_coldata_collapsed_binarize),
               yLabels = names(Probiotic_MEs),
               ySymbols = names(Probiotic_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = Probiotic_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
Probiotic_moduleTraitCor_df <- as.data.frame(Probiotic_moduleTraitCor)
Probiotic_moduleTraitCor_df$mod_names <- row.names(Probiotic_moduleTraitCor_df)
Probiotic_moduleTraitCor_df <- Probiotic_moduleTraitCor_df[,c("mod_names","Condition.Bacillus_pumilus_RI0695.vs.Untreated_control")]
Probiotic_moduleTraitPvalue_df <- as.data.frame(Probiotic_moduleTraitPvalue)
Probiotic_moduleTraitPvalue_df$mod_names <- row.names(Probiotic_moduleTraitPvalue_df)
Probiotic_moduleTraitPvalue_df <- Probiotic_moduleTraitPvalue_df[,c("mod_names","Condition.Bacillus_pumilus_RI0695.vs.Untreated_control")]
colnames(Probiotic_moduleTraitPvalue_df)[2] <- "pvalue"
Probiotic_moduleTraitCor_Pval_df <- join(Probiotic_moduleTraitCor_df, Probiotic_moduleTraitPvalue_df, by = "mod_names")

# Significantly correlated modules
Probiotic_moduleTraitCor_Pval_df[order(Probiotic_moduleTraitCor_Pval_df$pvalue),]
class(Probiotic_moduleTraitCor_Pval_df$pvalue) # numeric
Probiotic_moduleTraitCor_Pval_df_sig <- Probiotic_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
Probiotic_moduleTraitCor_Pval_df_sig

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
Probiotic_moduleTraitCor_Pval_df_sig_list <- Probiotic_moduleTraitCor_Pval_df_sig$mod_names
Probiotic_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Probiotic_moduleTraitCor_Pval_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
matrix_common= Probiotic_dds_rlog_matrix_common
moduleColors= Probiotic_moduleColors
lookup =   C_vir_rtracklayer_apop_product_final

length(colnames(matrix_common)[moduleColors == "grey60"]) # 368
length(colnames(matrix_common)[moduleColors == "lightgreen"]) #348 

lookup_mod_apop <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$ID %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
names(Probiotic_moduleTraitCor_Pval_df_sig_list_rm) <- c("grey60","lightgreen")
Probiotic_module_apop <- lapply(Probiotic_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop)
Probiotic_module_apop_df <- do.call(rbind,Probiotic_module_apop)
Probiotic_module_apop_df$mod_names <- gsub("\\..*","",row.names(Probiotic_module_apop_df))
Probiotic_module_apop_df$mod_names <- gsub("^","ME",Probiotic_module_apop_df$mod_names)
# add module significance
Probiotic_module_apop_df <- left_join(Probiotic_module_apop_df,Probiotic_moduleTraitCor_Pval_df_sig)
Probiotic_module_apop_df$exp <- "Probiotic"

## Gene relationship to trait and important modules: Gene Significance and Module Membership
# Define variable injected 
Probiotic_injection = as.data.frame(Probiotic_coldata_collapsed_binarize$Condition.Bacillus_pumilus_RI0695.vs.Untreated_control);
names(Probiotic_injection) = "challenge"
# names (colors) of the modules
Probiotic_modNames = substring(names(Probiotic_MEs), 3)
Probiotic_geneModuleMembership = as.data.frame(cor(Probiotic_dds_rlog_matrix_common, Probiotic_MEs, use = "p"))
Probiotic_MMPvalue = as.data.frame(corPvalueStudent(as.matrix(Probiotic_geneModuleMembership), Probiotic_nSamples))

names(Probiotic_geneModuleMembership) = paste("MM", Probiotic_modNames, sep="")
names(Probiotic_MMPvalue) = paste("p.MM", Probiotic_modNames, sep="")

Probiotic_geneTraitSignificance = as.data.frame(cor(Probiotic_dds_rlog_matrix_common,Probiotic_injection, use = "p"))
Probiotic_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(Probiotic_geneTraitSignificance), Probiotic_nSamples))

names(Probiotic_geneTraitSignificance) = paste("GS.", names(Probiotic_injection), sep="")
names(Probiotic_GSPvalue) = paste("p.GS.", names(Probiotic_injection), sep="")

## Intramodular analysis: identifying genes with high GS and MM
Probiotic_module = "lightgreen"  # extremely strong correlation 
Probiotic_column = match(Probiotic_module, Probiotic_modNames)
Probiotic_moduleGenes = Probiotic_moduleColors==Probiotic_module
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(Probiotic_geneModuleMembership[Probiotic_moduleGenes, Probiotic_column]),
                   abs(Probiotic_geneTraitSignificance[Probiotic_moduleGenes, 1]),
                   xlab = paste("Module Membership in", Probiotic_module, "module"),
                   ylab = "Gene significance for challenge",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = Probiotic_module)

## IDENTIFY HUB GENES IN EACH SIG MODULE ##
Probiotic_colorh = c("MEgrey60"  ,   "MElightgreen")

#Probiotic_Module_hub_genes <- chooseTopHubInEachModule(
#  Probiotic_dds_rlog_matrix_common, 
#  Probiotic_colorh, 
#  power = 2,  # power used for the adjacency network
#  type = "signed hybrid", 
#  corFnc = "bicor"
#)
class(Probiotic_Module_hub_genes)
Probiotic_Module_hub_genes_df <- as.data.frame(Probiotic_Module_hub_genes)
colnames(Probiotic_Module_hub_genes_df)[1] <- "ID"
Probiotic_Module_hub_genes_apop <- merge(Probiotic_Module_hub_genes_df, C_vir_rtracklayer, by = "ID")
nrow(Probiotic_Module_hub_genes_apop) 

#### ROD RES ####
# Define numbers of genes and samples
ROD_Res_nGenes = ncol(ROD_Resistant_dds_rlog_matrix_common)
ROD_Res_nSamples = nrow(ROD_Resistant_dds_rlog_matrix_common)

# Recalculate MEs with color labels
ROD_Res_MEs0 = moduleEigengenes(ROD_Resistant_dds_rlog_matrix_common, ROD_Res_moduleColors)$eigengenes
ROD_Res_MEs = orderMEs(ROD_Res_MEs0)
ROD_Res_moduleTraitCor = cor(ROD_Res_MEs, ROD_Resistant_coldata_collapsed_binarize , use = "p");
ROD_Res_moduleTraitPvalue = corPvalueStudent(ROD_Res_moduleTraitCor, ROD_Res_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trait
sizeGrWindow(10,6)
# Will display correlations and their p-values
ROD_Res_textMatrix = paste(signif(ROD_Res_moduleTraitCor, 2), "\n(",
                             signif(ROD_Res_moduleTraitPvalue, 1), ")", sep = "");
dim(ROD_Res_textMatrix) = dim(ROD_Res_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = ROD_Res_moduleTraitCor,
               xLabels = names(ROD_Resistant_coldata_collapsed_binarize),
               yLabels = names(ROD_Res_MEs),
               ySymbols =names(ROD_Res_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = ROD_Res_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
ROD_Res_moduleTraitCor_df <- as.data.frame(ROD_Res_moduleTraitCor)
ROD_Res_moduleTraitCor_df$mod_names <- row.names(ROD_Res_moduleTraitCor_df)
ROD_Res_moduleTraitCor_df <- ROD_Res_moduleTraitCor_df[,c("mod_names","Condition.Resistant_Challenge.vs.Control_Resistant")]
ROD_Res_moduleTraitPvalue_df <- as.data.frame(ROD_Res_moduleTraitPvalue)
ROD_Res_moduleTraitPvalue_df$mod_names <- row.names(ROD_Res_moduleTraitPvalue_df)
ROD_Res_moduleTraitPvalue_df <- ROD_Res_moduleTraitPvalue_df[,c("mod_names","Condition.Resistant_Challenge.vs.Control_Resistant")]
colnames(ROD_Res_moduleTraitPvalue_df)[2] <- "pvalue"
ROD_Res_moduleTraitCor_Pval_df <- join(ROD_Res_moduleTraitCor_df, ROD_Res_moduleTraitPvalue_df, by = "mod_names")

# Significantly correlated modules
ROD_Res_moduleTraitCor_Pval_df[order(ROD_Res_moduleTraitCor_Pval_df$pvalue),]
class(ROD_Res_moduleTraitCor_Pval_df$pvalue) # numeric
ROD_Res_moduleTraitCor_Pval_df_sig <- ROD_Res_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
ROD_Res_moduleTraitCor_Pval_df_sig

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
ROD_Res_moduleTraitCor_Pval_df_sig_list <- ROD_Res_moduleTraitCor_Pval_df_sig$mod_names
ROD_Res_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(ROD_Res_moduleTraitCor_Pval_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
matrix_common= ROD_Resistant_dds_rlog_matrix_common
moduleColors= ROD_Res_moduleColors
lookup =   C_vir_rtracklayer_apop_product_final

lookup_mod_apop <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$ID %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
names(ROD_Res_moduleTraitCor_Pval_df_sig_list_rm) <- c("plum1",  "orange")
ROD_Res_module_apop <- lapply(ROD_Res_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop)
ROD_Res_module_apop_df <- do.call(rbind,ROD_Res_module_apop)
ROD_Res_module_apop_df$mod_names <- gsub("\\..*","",row.names(ROD_Res_module_apop_df))
ROD_Res_module_apop_df$mod_names <- gsub("^","ME",ROD_Res_module_apop_df$mod_names)
# add module significance
ROD_Res_module_apop_df <- left_join(ROD_Res_module_apop_df,ROD_Res_moduleTraitCor_Pval_df_sig)
ROD_Res_module_apop_df$exp <- "ROD_Res"

## Gene relationship to trait and important modules: Gene Significance and Module Membership
# Define variable injected 
ROD_Res_injection = as.data.frame(ROD_Resistant_coldata_collapsed_binarize$Condition.Resistant_Challenge.vs.Control_Resistant);
names(ROD_Res_injection) = "challenge"
# names (colors) of the modules
ROD_Res_modNames = substring(names(ROD_Res_MEs), 3)
ROD_Res_geneModuleMembership = as.data.frame(cor(ROD_Resistant_dds_rlog_matrix_common, ROD_Res_MEs, use = "p"))
ROD_Res_MMPvalue = as.data.frame(corPvalueStudent(as.matrix(ROD_Res_geneModuleMembership), ROD_Res_nSamples))

names(ROD_Res_geneModuleMembership) = paste("MM", ROD_Res_modNames, sep="")
names(ROD_Res_MMPvalue) = paste("p.MM", ROD_Res_modNames, sep="")

ROD_Res_geneTraitSignificance = as.data.frame(cor(ROD_Resistant_dds_rlog_matrix_common,ROD_Res_injection, use = "p"))
ROD_Res_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(ROD_Res_geneTraitSignificance), ROD_Res_nSamples))

names(ROD_Res_geneTraitSignificance) = paste("GS.", names(ROD_Res_injection), sep="")
names(ROD_Res_GSPvalue) = paste("p.GS.", names(ROD_Res_injection), sep="")

## Intramodular analysis: identifying genes with high GS and MM
ROD_Res_module = "orange"  
ROD_Res_column = match(ROD_Res_module, ROD_Res_modNames)
ROD_Res_moduleGenes = ROD_Res_moduleColors==ROD_Res_module
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(ROD_Res_geneModuleMembership[ROD_Res_moduleGenes, ROD_Res_column]),
                   abs(ROD_Res_geneTraitSignificance[ROD_Res_moduleGenes, 1]),
                   xlab = paste("Module Membership in", ROD_Res_module, "module"),
                   ylab = "Gene significance for challenge",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = ROD_Res_module)

## IDENTIFY HUB GENES IN EACH SIG MODULE ##
ROD_Res_colorh = c("MEplum1",  "MEorange")

#ROD_Res_Module_hub_genes <- chooseTopHubInEachModule(
#  ROD_Resistant_dds_rlog_matrix_common, 
#  ROD_Res_colorh, 
#  power = 2,  # power used for the adjacency network
#  type = "signed hybrid", 
#  corFnc = "bicor"
#)
class(ROD_Res_Module_hub_genes)
ROD_Res_Module_hub_genes_df <- as.data.frame(ROD_Res_Module_hub_genes)
colnames(ROD_Res_Module_hub_genes_df)[1] <- "ID"
ROD_Res_Module_hub_genes <- merge(ROD_Res_Module_hub_genes_df, C_vir_rtracklayer, by = "ID")
nrow(ROD_Res_Module_hub_genes) 

#### ROD SUS ####
# Define numbers of genes and samples
ROD_Sus_nGenes = ncol(ROD_Susceptible_dds_rlog_matrix_common)
ROD_Sus_nSamples = nrow(ROD_Susceptible_dds_rlog_matrix_common)

# Recalculate MEs with color labels
ROD_Sus_MEs0 = moduleEigengenes(ROD_Susceptible_dds_rlog_matrix_common, ROD_Sus_moduleColors)$eigengenes
ROD_Sus_MEs = orderMEs(ROD_Sus_MEs0)
ROD_Sus_moduleTraitCor = cor(ROD_Sus_MEs, ROD_Susceptible_coldata_collapsed_binarize , use = "p");
ROD_Sus_moduleTraitPvalue = corPvalueStudent(ROD_Sus_moduleTraitCor, ROD_Sus_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trait
sizeGrWindow(10,6)
# Will display correlations and their p-values
ROD_Sus_textMatrix = paste(signif(ROD_Sus_moduleTraitCor, 2), "\n(",
                           signif(ROD_Sus_moduleTraitPvalue, 1), ")", sep = "");
dim(ROD_Sus_textMatrix) = dim(ROD_Sus_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = ROD_Sus_moduleTraitCor,
               xLabels = names(ROD_Susceptible_coldata_collapsed_binarize),
               yLabels = names(ROD_Sus_MEs),
               ySymbols =names(ROD_Sus_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = ROD_Sus_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
ROD_Sus_moduleTraitCor_df <- as.data.frame(ROD_Sus_moduleTraitCor)
ROD_Sus_moduleTraitCor_df$mod_names <- row.names(ROD_Sus_moduleTraitCor_df)
ROD_Sus_moduleTraitCor_df <- ROD_Sus_moduleTraitCor_df[,c("mod_names","Condition.Late_Susecptible.vs.Early_Susceptible")]
ROD_Sus_moduleTraitPvalue_df <- as.data.frame(ROD_Sus_moduleTraitPvalue)
ROD_Sus_moduleTraitPvalue_df$mod_names <- row.names(ROD_Sus_moduleTraitPvalue_df)
ROD_Sus_moduleTraitPvalue_df <- ROD_Sus_moduleTraitPvalue_df[,c("mod_names","Condition.Late_Susecptible.vs.Early_Susceptible")]
colnames(ROD_Sus_moduleTraitPvalue_df)[2] <- "pvalue"
ROD_Sus_moduleTraitCor_Pval_df <- join(ROD_Sus_moduleTraitCor_df, ROD_Sus_moduleTraitPvalue_df, by = "mod_names")

# Significantly correlated modules
ROD_Sus_moduleTraitCor_Pval_df[order(ROD_Sus_moduleTraitCor_Pval_df$pvalue),]
class(ROD_Sus_moduleTraitCor_Pval_df$pvalue) # numeric
ROD_Sus_moduleTraitCor_Pval_df_sig <- ROD_Sus_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
ROD_Sus_moduleTraitCor_Pval_df_sig

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
ROD_Sus_moduleTraitCor_Pval_df_sig_list <- ROD_Sus_moduleTraitCor_Pval_df_sig$mod_names
ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(ROD_Sus_moduleTraitCor_Pval_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
matrix_common= ROD_Susceptible_dds_rlog_matrix_common
moduleColors= ROD_Sus_moduleColors
lookup =   C_vir_rtracklayer_apop_product_final

lookup_mod_apop <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$ID %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
names(ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm) <- c("turquoise", "brown" )
ROD_Sus_module_apop <- lapply(ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop)
ROD_Sus_module_apop_df <- do.call(rbind,ROD_Sus_module_apop)
ROD_Sus_module_apop_df$mod_names <- gsub("\\..*","",row.names(ROD_Sus_module_apop_df))
ROD_Sus_module_apop_df$mod_names <- gsub("^","ME",ROD_Sus_module_apop_df$mod_names)
# add module significance
ROD_Sus_module_apop_df <- left_join(ROD_Sus_module_apop_df,ROD_Sus_moduleTraitCor_Pval_df_sig)
ROD_Sus_module_apop_df$exp <- "ROD_Sus"

## Gene relationship to trait and important modules: Gene Significance and Module Membership
# Define variable injected 
ROD_Sus_injection = as.data.frame(ROD_Susceptible_coldata_collapsed_binarize$Condition.Late_Susecptible.vs.Early_Susceptible);
names(ROD_Sus_injection) = "challenge"
# names (colors) of the modules
ROD_Sus_modNames = substring(names(ROD_Sus_MEs), 3)
ROD_Sus_geneModuleMembership = as.data.frame(cor(ROD_Susceptible_dds_rlog_matrix_common, ROD_Sus_MEs, use = "p"))
ROD_Sus_MMPvalue = as.data.frame(corPvalueStudent(as.matrix(ROD_Sus_geneModuleMembership), ROD_Sus_nSamples))

names(ROD_Sus_geneModuleMembership) = paste("MM", ROD_Sus_modNames, sep="")
names(ROD_Sus_MMPvalue) = paste("p.MM", ROD_Sus_modNames, sep="")

ROD_Sus_geneTraitSignificance = as.data.frame(cor(ROD_Susceptible_dds_rlog_matrix_common,ROD_Sus_injection, use = "p"))
ROD_Sus_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(ROD_Sus_geneTraitSignificance), ROD_Sus_nSamples))

names(ROD_Sus_geneTraitSignificance) = paste("GS.", names(ROD_Sus_injection), sep="")
names(ROD_Sus_GSPvalue) = paste("p.GS.", names(ROD_Sus_injection), sep="")

## Intramodular analysis: identifying genes with high GS and MM
ROD_Sus_module = "brown"  
ROD_Sus_column = match(ROD_Sus_module, ROD_Sus_modNames)
ROD_Sus_moduleGenes = ROD_Sus_moduleColors==ROD_Sus_module
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(ROD_Sus_geneModuleMembership[ROD_Sus_moduleGenes, ROD_Sus_column]),
                   abs(ROD_Sus_geneTraitSignificance[ROD_Sus_moduleGenes, 1]),
                   xlab = paste("Module Membership in", ROD_Sus_module, "module"),
                   ylab = "Gene significance for challenge",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = ROD_Sus_module) # perfect correlation...something has probably gone wrong here!

## IDENTIFY HUB GENES IN EACH SIG MODULE ##
ROD_Sus_colorh = c("MEturquoise",  "MEbrown")

#ROD_Sus_Module_hub_genes <- chooseTopHubInEachModule(
#  ROD_Susceptible_dds_rlog_matrix_common, 
#  ROD_Sus_colorh, 
#  power = 2,  # power used for the adjacency network
#  type = "signed hybrid", 
#  corFnc = "bicor"
#)
class(ROD_Sus_Module_hub_genes)
ROD_Sus_Module_hub_genes_df <- as.data.frame(ROD_Sus_Module_hub_genes)
colnames(ROD_Sus_Module_hub_genes_df)[1] <- "ID"
ROD_Sus_Module_hub_genes <- merge(ROD_Sus_Module_hub_genes_df, C_vir_rtracklayer, by = "ID")
nrow(ROD_Res_Module_hub_genes) 

#### PRO_RE22 ####

# Define numbers of genes and samples
Pro_RE22_nGenes = ncol(Pro_RE22_dds_rlog_matrix_common)
Pro_RE22_nSamples = nrow(Pro_RE22_dds_rlog_matrix_common)

# Recalculate MEs with color labels
Pro_RE22_MEs0 = moduleEigengenes(Pro_RE22_dds_rlog_matrix_common, Pro_RE22_moduleColors)$eigengenes
Pro_RE22_MEs = orderMEs(Pro_RE22_MEs0)
Pro_RE22_moduleTraitCor = cor(Pro_RE22_MEs, Pro_RE22_coldata_collapsed_binarize , use = "p");
Pro_RE22_moduleTraitPvalue = corPvalueStudent(Pro_RE22_moduleTraitCor, Pro_RE22_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trait
sizeGrWindow(10,6)
# Will display correlations and their p-values
Pro_RE22_textMatrix = paste(signif(Pro_RE22_moduleTraitCor, 2), "\n(",
                           signif(Pro_RE22_moduleTraitPvalue, 1), ")", sep = "");
dim(Pro_RE22_textMatrix) = dim(Pro_RE22_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = Pro_RE22_moduleTraitCor,
               xLabels = names(Pro_RE22_coldata_collapsed_binarize),
               yLabels = names(Pro_RE22_MEs),
               ySymbols =names(Pro_RE22_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = Pro_RE22_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with  RI0695 (high correlation and low P value)?
Pro_RE22_RI_moduleTraitCor_df <- as.data.frame(Pro_RE22_moduleTraitCor)
Pro_RE22_RI_moduleTraitCor_df$mod_names <- row.names(Pro_RE22_RI_moduleTraitCor_df)
Pro_RE22_RI_moduleTraitCor_df <- Pro_RE22_RI_moduleTraitCor_df[,c("mod_names","Condition.Control_no_treatment.vs.Bacillus_pumilus")]
Pro_RE22_RI_moduleTraitPvalue_df <- as.data.frame(Pro_RE22_moduleTraitPvalue)
Pro_RE22_RI_moduleTraitPvalue_df$mod_names <- row.names(Pro_RE22_RI_moduleTraitPvalue_df)
Pro_RE22_RI_moduleTraitPvalue_df <- Pro_RE22_RI_moduleTraitPvalue_df[,c("mod_names","Condition.Control_no_treatment.vs.Bacillus_pumilus")]
colnames(Pro_RE22_RI_moduleTraitPvalue_df)[2] <- "pvalue"
Pro_RE22_RI_moduleTraitCor_Pval_df <- join(Pro_RE22_RI_moduleTraitCor_df, Pro_RE22_RI_moduleTraitPvalue_df, by = "mod_names")

Pro_RE22_S4_moduleTraitCor_df <- as.data.frame(Pro_RE22_moduleTraitCor)
Pro_RE22_S4_moduleTraitCor_df$mod_names <- row.names(Pro_RE22_S4_moduleTraitCor_df)
Pro_RE22_S4_moduleTraitCor_df <- Pro_RE22_S4_moduleTraitCor_df[,c("mod_names","Condition.Phaeobacter_inhibens.vs.Control_no_treatment" )]
Pro_RE22_S4_moduleTraitPvalue_df <- as.data.frame(Pro_RE22_moduleTraitPvalue)
Pro_RE22_S4_moduleTraitPvalue_df$mod_names <- row.names(Pro_RE22_S4_moduleTraitPvalue_df)
Pro_RE22_S4_moduleTraitPvalue_df <- Pro_RE22_S4_moduleTraitPvalue_df[,c("mod_names","Condition.Phaeobacter_inhibens.vs.Control_no_treatment" )]
colnames(Pro_RE22_S4_moduleTraitPvalue_df)[2] <- "pvalue"
Pro_RE22_S4_moduleTraitCor_Pval_df <- join(Pro_RE22_S4_moduleTraitCor_df, Pro_RE22_S4_moduleTraitPvalue_df, by = "mod_names")

Pro_RE22_moduleTraitCor_df <- as.data.frame(Pro_RE22_moduleTraitCor)
Pro_RE22_moduleTraitCor_df$mod_names <- row.names(Pro_RE22_moduleTraitCor_df)
Pro_RE22_moduleTraitCor_df <- Pro_RE22_moduleTraitCor_df[,c("mod_names","Condition.Vibrio_coralliilyticus_RE22_exposure_6h.vs.Control_no_treatment")]
Pro_RE22_moduleTraitPvalue_df <- as.data.frame(Pro_RE22_moduleTraitPvalue)
Pro_RE22_moduleTraitPvalue_df$mod_names <- row.names(Pro_RE22_moduleTraitPvalue_df)
Pro_RE22_moduleTraitPvalue_df <- Pro_RE22_moduleTraitPvalue_df[,c("mod_names","Condition.Vibrio_coralliilyticus_RE22_exposure_6h.vs.Control_no_treatment")]
colnames(Pro_RE22_moduleTraitPvalue_df)[2] <- "pvalue"
Pro_RE22_moduleTraitCor_Pval_df <- join(Pro_RE22_moduleTraitCor_df, Pro_RE22_moduleTraitPvalue_df, by = "mod_names")

# Significantly correlated modules
Pro_RE22_RI_moduleTraitCor_Pval_df[order(Pro_RE22_RI_moduleTraitCor_Pval_df$pvalue),]
class(Pro_RE22_RI_moduleTraitCor_Pval_df$pvalue) # numeric
# subset just for positive associations
Pro_RE22_RI_moduleTraitCor_Pval_df_sig <- Pro_RE22_RI_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
Pro_RE22_RI_moduleTraitCor_Pval_df_sig # 21

Pro_RE22_S4_moduleTraitCor_Pval_df[order(Pro_RE22_S4_moduleTraitCor_Pval_df$pvalue),]
class(Pro_RE22_S4_moduleTraitCor_Pval_df$pvalue) # numeric
Pro_RE22_S4_moduleTraitCor_Pval_df_sig <- Pro_RE22_S4_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
Pro_RE22_S4_moduleTraitCor_Pval_df_sig # 17

Pro_RE22_moduleTraitCor_Pval_df[order(Pro_RE22_moduleTraitCor_Pval_df$pvalue),]
class(Pro_RE22_moduleTraitCor_Pval_df$pvalue) # numeric
Pro_RE22_moduleTraitCor_Pval_df_sig <- Pro_RE22_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
Pro_RE22_moduleTraitCor_Pval_df_sig #31

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list <- Pro_RE22_RI_moduleTraitCor_Pval_df_sig$mod_names
Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list <- Pro_RE22_S4_moduleTraitCor_Pval_df_sig$mod_names
Pro_RE22_moduleTraitCor_Pval_df_sig_list <- Pro_RE22_moduleTraitCor_Pval_df_sig$mod_names

Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list, "ME")
Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list, "ME")
Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Pro_RE22_moduleTraitCor_Pval_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
names(Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm) <- c("mediumpurple3",  "black"    ,      "steelblue"  ,    "purple"   ,      "salmon"    ,     "midnightblue" ,  "plum1"      ,    "darkolivegreen", "cyan"   ,        "magenta"  ,      "green"  ,        "lightgreen",    
                                                         "saddlebrown" ,  "coral1"  ,       "salmon4"    ,    "blue"     ,      "royalblue" ,     "red"          ,  "darkgreen"  ,    "yellow"        , "orange"  )
names(Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm) <- c("mediumpurple3",  "midnightblue" ,  "plum1"   ,       "darkolivegreen", "darkorange"  ,   "cyan"  ,"greenyellow" ,   "magenta","green"  ,"lightgreen",     "saddlebrown" ,   "coral1" ,       
                                                           "paleturquoise",  "honeydew1"    ,  "salmon4" ,       "blue"          , "red" )
names(Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm) <- c( "darkseagreen4"  ,"white"     ,     "mediumpurple3",  "black"       ,   "steelblue"     , "brown4"    ,     "purple" , "salmon"         , "midnightblue" ,  "plum1"         , "darkolivegreen" ,"cyan"       ,   
                                                         "yellowgreen"   , "magenta"  ,      "green"       ,   "lightgreen" ,    "saddlebrown"  ,  "coral1"   ,      "blue"  ,  "lavenderblush3", "navajowhite2" ,  "darkslateblue" , "pink"           ,"lightpink4" ,   
                                                         "coral2"        , "maroon"   ,      "red"         ,   "darkred"    ,    "darkgreen"    ,  "yellow"   ,      "grey"  )
matrix_common= Pro_RE22_dds_rlog_matrix_common
moduleColors= Pro_RE22_moduleColors
lookup =   C_vir_rtracklayer_apop_product_final

lookup_mod_apop <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$ID %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
Pro_RE22_RI_module_apop <- lapply(Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop)
Pro_RE22_RI_module_apop_df <- do.call(rbind,Pro_RE22_RI_module_apop)
Pro_RE22_RI_module_apop_df$mod_names <- gsub("\\..*","",row.names(Pro_RE22_RI_module_apop_df))
Pro_RE22_RI_module_apop_df$mod_names <- gsub("^","ME",Pro_RE22_RI_module_apop_df$mod_names)
# add module significance
Pro_RE22_RI_module_apop_df <- left_join(Pro_RE22_RI_module_apop_df,Pro_RE22_RI_moduleTraitCor_Pval_df_sig)
Pro_RE22_RI_module_apop_df$exp <- "Pro_RE22_RI"

Pro_RE22_S4_module_apop <- lapply(Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop)
Pro_RE22_S4_module_apop_df <- do.call(rbind,Pro_RE22_S4_module_apop)
Pro_RE22_S4_module_apop_df$mod_names <- gsub("\\..*","",row.names(Pro_RE22_S4_module_apop_df))
Pro_RE22_S4_module_apop_df$mod_names <- gsub("^","ME",Pro_RE22_S4_module_apop_df$mod_names)
# add module significance
Pro_RE22_S4_module_apop_df <- left_join(Pro_RE22_S4_module_apop_df,Pro_RE22_S4_moduleTraitCor_Pval_df_sig)
Pro_RE22_S4_module_apop_df$exp <- "Pro_RE22_S4" 

Pro_RE22_module_apop <- lapply(Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop)
Pro_RE22_module_apop_df <- do.call(rbind,Pro_RE22_module_apop)
Pro_RE22_module_apop_df$mod_names <- gsub("\\..*","",row.names(Pro_RE22_module_apop_df))
Pro_RE22_module_apop_df$mod_names <- gsub("^","ME",Pro_RE22_module_apop_df$mod_names)
# add module significance
Pro_RE22_module_apop_df <- left_join(Pro_RE22_module_apop_df,Pro_RE22_moduleTraitCor_Pval_df_sig)
Pro_RE22_module_apop_df$exp <- "Pro_RE22_RE22"

## Gene relationship to trait and important modules: Gene Significance and Module Membership
# Define variable injected 
Pro_RE22_injection = as.data.frame(Pro_RE22_coldata_collapsed_binarize$Condition.Vibrio_coralliilyticus_RE22_exposure_6h.vs.Control_no_treatment);
names(Pro_RE22_injection) = "challenge"
# names (colors) of the modules
Pro_RE22_modNames = substring(names(Pro_RE22_MEs), 3)
Pro_RE22_geneModuleMembership = as.data.frame(cor(Pro_RE22_dds_rlog_matrix_common, Pro_RE22_MEs, use = "p"))
Pro_RE22_MMPvalue = as.data.frame(corPvalueStudent(as.matrix(Pro_RE22_geneModuleMembership), Pro_RE22_nSamples))

names(Pro_RE22_geneModuleMembership) = paste("MM", Pro_RE22_modNames, sep="")
names(Pro_RE22_MMPvalue) = paste("p.MM", Pro_RE22_modNames, sep="")

Pro_RE22_geneTraitSignificance = as.data.frame(cor(Pro_RE22_dds_rlog_matrix_common,Pro_RE22_injection, use = "p"))
Pro_RE22_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(Pro_RE22_geneTraitSignificance), Pro_RE22_nSamples))

names(Pro_RE22_geneTraitSignificance) = paste("GS.", names(Pro_RE22_injection), sep="")
names(Pro_RE22_GSPvalue) = paste("p.GS.", names(Pro_RE22_injection), sep="")

## Intramodular analysis: identifying genes with high GS and MM
Pro_RE22_module = "grey"  # poor correlation
Pro_RE22_column = match(Pro_RE22_module, Pro_RE22_modNames)
Pro_RE22_moduleGenes = Pro_RE22_moduleColors==Pro_RE22_module
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(Pro_RE22_geneModuleMembership[Pro_RE22_moduleGenes, Pro_RE22_column]),
                   abs(Pro_RE22_geneTraitSignificance[Pro_RE22_moduleGenes, 1]),
                   xlab = paste("Module Membership in", Pro_RE22_module, "module"),
                   ylab = "Gene significance for challenge",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = Pro_RE22_module) # perfect correlation...something has probably gone wrong here!


## IDENTIFY HUB GENES IN EACH SIG MODULE ##
Pro_RE22_colorh = c("MEmediumpurple3_apop","MEblack_apop",
"MEsteelblue_apop","MEpurple_apop","MEsalmon_apop","MEmidnightblue","MEplum1 ","MEcoral1",
"MEsalmon4","MEblue","MEred","MEorange","MEdarkolivegreen","MEdarkorange","MEcyan","MEgreenyellow",
"MEmagenta","MEgreen","MElightgreen","MEsaddlebrown","MEdarkseagreen4","MElavenderblush3","MEnavajowhite",
"MEpink","MElightpink4","MEcoral2","MEmaroon","MEdarkred","MEdarkgreen","MEyellow","MEgrey")
          
#Pro_RE22_Module_hub_genes <- chooseTopHubInEachModule(
#  Pro_RE22_dds_rlog_matrix_common, 
#  Pro_RE22_colorh, 
#  power = 2,  # power used for the adjacency network
#  type = "signed hybrid", 
#  corFnc = "bicor"
#)
class(Pro_RE22_Module_hub_genes)
Pro_RE22_Module_hub_genes_df <- as.data.frame(Pro_RE22_Module_hub_genes)
colnames(Pro_RE22_Module_hub_genes_df)[1] <- "ID"
Pro_RE22_Module_hub_genes <- merge(Pro_RE22_Module_hub_genes_df, C_vir_rtracklayer, by = "ID")
nrow(Pro_RE22_Module_hub_genes) 
# includes programmed cell death protein 6 as a hub gene, apoptosis regulatory protein Siva 

#### Pro_RE22 Pro ####
# define numbers of genes and sample
Pro_RE22_Pro_nGenes = ncol(Pro_RE22_dds_rlog_matrix_common_Pro)
Pro_RE22_Pro_nSamples = nrow(Pro_RE22_dds_rlog_matrix_common_Pro)

# Recalculate MEs with color labels
Pro_RE22_Pro_MEs0 = moduleEigengenes(Pro_RE22_dds_rlog_matrix_common_Pro, Pro_RE22_Pro_moduleColors)$eigengenes
Pro_RE22_Pro_MEs = orderMEs(Pro_RE22_Pro_MEs0)
Pro_RE22_Pro_moduleTraitCor = cor(Pro_RE22_Pro_MEs, Pro_RE22_coldata_Pro_collapsed_binarize , use = "p");
Pro_RE22_Pro_moduleTraitPvalue = corPvalueStudent(Pro_RE22_Pro_moduleTraitCor, Pro_RE22_Pro_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trait
sizeGrWindow(10,6)
# Will display correlations and their p-values
Pro_RE22_Pro_textMatrix = paste(signif(Pro_RE22_Pro_moduleTraitCor, 2), "\n(",
                            signif(Pro_RE22_Pro_moduleTraitPvalue, 1), ")", sep = "");
dim(Pro_RE22_Pro_textMatrix) = dim(Pro_RE22_Pro_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = Pro_RE22_Pro_moduleTraitCor,
               xLabels = names(Pro_RE22_coldata_Pro_collapsed_binarize),
               yLabels = names(Pro_RE22_Pro_MEs),
               ySymbols =names(Pro_RE22_Pro_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = Pro_RE22_Pro_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with probiotic treatment (high correlation and low P value)?
Pro_RE22_Pro_moduleTraitCor_df <- as.data.frame(Pro_RE22_Pro_moduleTraitCor)
Pro_RE22_Pro_moduleTraitCor_df$mod_names <- row.names(Pro_RE22_Pro_moduleTraitCor_df)
Pro_RE22_Pro_moduleTraitCor_df <- Pro_RE22_Pro_moduleTraitCor_df[,c("mod_names","Condition.Control_no_treatment.vs.Bacillus_pumilus")]
Pro_RE22_Pro_moduleTraitPvalue_df <- as.data.frame(Pro_RE22_Pro_moduleTraitPvalue)
Pro_RE22_Pro_moduleTraitPvalue_df$mod_names <- row.names(Pro_RE22_Pro_moduleTraitPvalue_df)
Pro_RE22_Pro_moduleTraitPvalue_df <- Pro_RE22_Pro_moduleTraitPvalue_df[,c("mod_names","Condition.Control_no_treatment.vs.Bacillus_pumilus")]
colnames(Pro_RE22_Pro_moduleTraitPvalue_df)[2] <- "pvalue"
Pro_RE22_Pro_moduleTraitCor_Pval_df <- join(Pro_RE22_Pro_moduleTraitCor_df, Pro_RE22_Pro_moduleTraitPvalue_df, by = "mod_names")

Pro_RE22_Pro_RI_moduleTraitCor_df <- as.data.frame(Pro_RE22_Pro_moduleTraitCor)
Pro_RE22_Pro_RI_moduleTraitCor_df$mod_names <- row.names(Pro_RE22_Pro_RI_moduleTraitCor_df)
Pro_RE22_Pro_RI_moduleTraitCor_df <- Pro_RE22_Pro_RI_moduleTraitCor_df[,c("mod_names","Condition.Phaeobacter_inhibens.vs.Control_no_treatment")]
Pro_RE22_Pro_RI_moduleTraitPvalue_df <- as.data.frame(Pro_RE22_Pro_moduleTraitPvalue)
Pro_RE22_Pro_RI_moduleTraitPvalue_df$mod_names <- row.names(Pro_RE22_Pro_RI_moduleTraitPvalue_df)
Pro_RE22_Pro_RI_moduleTraitPvalue_df <- Pro_RE22_Pro_RI_moduleTraitPvalue_df[,c("mod_names","Condition.Phaeobacter_inhibens.vs.Control_no_treatment")]
colnames(Pro_RE22_Pro_RI_moduleTraitPvalue_df)[2] <- "pvalue"
Pro_RE22_Pro_RI_moduleTraitCor_Pval_df <- join(Pro_RE22_Pro_RI_moduleTraitCor_df, Pro_RE22_Pro_RI_moduleTraitPvalue_df, by = "mod_names")

# Significantly correlated modules
Pro_RE22_Pro_moduleTraitCor_Pval_df[order(Pro_RE22_Pro_moduleTraitCor_Pval_df$pvalue),]
class(Pro_RE22_Pro_moduleTraitCor_Pval_df$pvalue) # numeric
# subset just for positive associations
Pro_RE22_Pro_moduleTraitCor_Pval_df_sig <- Pro_RE22_Pro_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
Pro_RE22_Pro_moduleTraitCor_Pval_df_sig # 27

# Significantly correlated modules
Pro_RE22_Pro_RI_moduleTraitCor_Pval_df[order(Pro_RE22_Pro_RI_moduleTraitCor_Pval_df$pvalue),]
class(Pro_RE22_Pro_RI_moduleTraitCor_Pval_df$pvalue) # numeric
# subset just for positive associations
Pro_RE22_Pro_RI_moduleTraitCor_Pval_df_sig <- Pro_RE22_Pro_RI_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
Pro_RE22_Pro_RI_moduleTraitCor_Pval_df_sig # 25

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
Pro_RE22_Pro_moduleTraitCor_Pval_df_sig_list <- Pro_RE22_Pro_moduleTraitCor_Pval_df_sig$mod_names
Pro_RE22_Pro_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Pro_RE22_Pro_moduleTraitCor_Pval_df_sig_list, "ME")

Pro_RE22_Pro_RI_moduleTraitCor_Pval_df_sig_list <- Pro_RE22_Pro_RI_moduleTraitCor_Pval_df_sig$mod_names
Pro_RE22_Pro_RI_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Pro_RE22_Pro_RI_moduleTraitCor_Pval_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
names(Pro_RE22_Pro_moduleTraitCor_Pval_df_sig_list_rm ) <- c("lavenderblush3", "mediumorchid"  , "brown"  ,"purple"        , "cyan"          , "plum1"         , "yellowgreen"  ,  "navajowhite2"  ,
                                                             "royalblue"     , "grey60"        , "yellow" ,"brown4"        , "antiquewhite2" , "palevioletred2", "skyblue3"     ,  "palevioletred3",
                                                              "black"        ,  "mediumpurple1",  "tan"   , "midnightblue" ,  "orange"       ,  "blue"         ,  "sienna4"     ,   "yellow3"       ,
                                                              "coral3"       ,  "pink"         ,  "red"  )
names(Pro_RE22_Pro_RI_moduleTraitCor_Pval_df_sig_list_rm ) <- c( "lavenderblush3" ,"brown"         , "salmon4"    ,    "purple"      ,   "yellowgreen" ,   "navajowhite2" ,  "royalblue"     , "yellow"        ,
                                                                 "darkslateblue"  ,"coral1"        , "skyblue1"   ,    "pink4"       ,   "black"       ,   "plum"         ,  "mediumpurple1" , "darkseagreen4" ,
                                                                  "lavenderblush2", "midnightblue" ,  "turquoise" ,     "lightgreen" ,    "orange"     ,    "coral2"      ,   "sienna4"      ,  "yellow3"       ,
                                                                  "darkolivegreen")
matrix_common= Pro_RE22_dds_rlog_matrix_common_Pro
moduleColors= Pro_RE22_Pro_moduleColors
lookup =   C_vir_rtracklayer_apop_product_final

lookup_mod_apop <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$ID %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
Pro_RE22_Pro_module_apop <- lapply(Pro_RE22_Pro_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop)
Pro_RE22_Pro_module_apop_df <- do.call(rbind,Pro_RE22_Pro_module_apop)
Pro_RE22_Pro_module_apop_df$mod_names <- gsub("\\..*","",row.names(Pro_RE22_Pro_module_apop_df))
Pro_RE22_Pro_module_apop_df$mod_names <- gsub("^","ME",Pro_RE22_Pro_module_apop_df$mod_names)
# add module significance
Pro_RE22_Pro_module_apop_df <- left_join(Pro_RE22_Pro_module_apop_df,Pro_RE22_Pro_moduleTraitCor_Pval_df_sig)
Pro_RE22_Pro_module_apop_df$exp <- "Pro_RE22_Pro_S4"

Pro_RE22_Pro_RI_module_apop <- lapply(Pro_RE22_Pro_RI_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop)
Pro_RE22_Pro_RI_module_apop_df <- do.call(rbind,Pro_RE22_Pro_RI_module_apop)
Pro_RE22_Pro_RI_module_apop_df$mod_names <- gsub("\\..*","",row.names(Pro_RE22_Pro_RI_module_apop_df))
Pro_RE22_Pro_RI_module_apop_df$mod_names <- gsub("^","ME",Pro_RE22_Pro_RI_module_apop_df$mod_names)
# add module significance
Pro_RE22_Pro_RI_module_apop_df <- left_join(Pro_RE22_Pro_RI_module_apop_df,Pro_RE22_Pro_RI_moduleTraitCor_Pval_df_sig)
Pro_RE22_Pro_RI_module_apop_df$exp <- "Pro_RE22_Pro_RI"

## Gene relationship to trait and important modules: Gene Significance and Module Membership
# Define variable injected 
# do later if necessary

## Intramodular analysis: identifying genes with high GS and MM
# later if necessary 

## IDENTIFY HUB GENES IN EACH SIG MODULE ##
Pro_RE22_Pro_colorh = c("MElavenderblush3", "MEmediumorchid"   ,"MEbrown"    , "MEpurple"        ,"MEcyan"  ,"MEplum1"         ,"MEyellowgreen"   ,
                        "MEnavajowhite2"  , "MEroyalblue"      ,"MEgrey60"   , "MEyellow"        ,"MEbrown4","MEantiquewhite2" ,"MEpalevioletred2",
                         "MEskyblue3"     ,  "MEpalevioletred3", "MEblack"   ,  "MEmediumpurple1", "MEtan"  , "MEmidnightblue" , "MEorange"        ,
                         "MEblue"         ,  "MEsienna4"       , "MEyellow3" ,  "MEcoral3"       , "MEpink" , "MEred" )

#Pro_RE22_Pro_Module_hub_genes <- chooseTopHubInEachModule(
#  Pro_RE22_dds_rlog_matrix_common_Pro, 
#  Pro_RE22_Pro_colorh, 
#  power = 2,  # power used for the adjacency network
#  type = "signed hybrid", 
#  corFnc = "bicor"
#)
class(Pro_RE22_Pro_Module_hub_genes)
Pro_RE22_Pro_Module_hub_genes_df <- as.data.frame(Pro_RE22_Pro_Module_hub_genes)
colnames(Pro_RE22_Pro_Module_hub_genes_df)[1] <- "ID"
Pro_RE22_Pro_Module_hub_genes <- merge(Pro_RE22_Pro_Module_hub_genes_df, C_vir_rtracklayer, by = "ID")
nrow(Pro_RE22_Pro_Module_hub_genes) 
# no apoptosis hub genes 

Pro_RE22_Pro_RI_colorh = c("MElavenderblush3", "MEbrown"         , "MEsalmon4"       , "MEpurple"       ,  "MEyellowgreen",    "MEnavajowhite2",   "MEroyalblue"  ,   
                          "MEyellow"         ,"MEdarkslateblue"  ,"MEcoral1"         ,"MEskyblue1"      , "MEpink4"       ,   "MEblack"        ,  "MEplum"        ,  
                           "MEmediumpurple1" , "MEdarkseagreen4" , "MElavenderblush2", "MEmidnightblue" ,  "MEturquoise"  ,    "MElightgreen"  ,   "MEorange"     ,   
                           "MEcoral2"        , "MEsienna4"       , "MEyellow3"       , "MEdarkolivegreen")

Pro_RE22_Pro_RI_Module_hub_genes <- chooseTopHubInEachModule(
  Pro_RE22_dds_rlog_matrix_common_Pro, 
  Pro_RE22_Pro_RI_colorh, 
  power = 2,  # power used for the adjacency network
  type = "signed hybrid", 
  corFnc = "bicor"
)
class(Pro_RE22_Pro_RI_Module_hub_genes)
Pro_RE22_Pro_RI_Module_hub_genes_df <- as.data.frame(Pro_RE22_Pro_RI_Module_hub_genes)
colnames(Pro_RE22_Pro_RI_Module_hub_genes_df)[1] <- "ID"
Pro_RE22_Pro_RI_Module_hub_genes <- merge(Pro_RE22_Pro_RI_Module_hub_genes_df, C_vir_rtracklayer, by = "ID")
nrow(Pro_RE22_Pro_RI_Module_hub_genes) 
#baculoviral IAP repeat-containing protein 7-like, p53 and DNA damage-regulated protein 1-like are both hub genes 

#### Pro_RE22_RE22 ####

# define numbers of genes and sample
Pro_RE22_RE22_nGenes = ncol(Pro_RE22_dds_rlog_matrix_common_RE22)
Pro_RE22_RE22_nSamples = nrow(Pro_RE22_dds_rlog_matrix_common_RE22)

# Recalculate MEs with color labels
Pro_RE22_RE22_MEs0 = moduleEigengenes(Pro_RE22_dds_rlog_matrix_common_RE22, Pro_RE22_RE22_moduleColors)$eigengenes
Pro_RE22_RE22_MEs = orderMEs(Pro_RE22_RE22_MEs0)
Pro_RE22_RE22_moduleTraitCor = cor(Pro_RE22_RE22_MEs, Pro_RE22_coldata_RE22_collapsed_binarize , use = "p");
Pro_RE22_RE22_moduleTraitPvalue = corPvalueStudent(Pro_RE22_RE22_moduleTraitCor, Pro_RE22_RE22_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trait
sizeGrWindow(10,6)
# Will display correlations and their p-values
Pro_RE22_RE22_textMatrix = paste(signif(Pro_RE22_RE22_moduleTraitCor, 2), "\n(",
                                signif(Pro_RE22_RE22_moduleTraitPvalue, 1), ")", sep = "");
dim(Pro_RE22_RE22_textMatrix) = dim(Pro_RE22_RE22_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = Pro_RE22_RE22_moduleTraitCor,
               xLabels = names(Pro_RE22_coldata_RE22_collapsed_binarize),
               yLabels = names(Pro_RE22_RE22_MEs),
               ySymbols =names(Pro_RE22_RE22_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = Pro_RE22_RE22_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with probiotic treatment (high correlation and low P value)?
Pro_RE22_RE22_moduleTraitCor_df <- as.data.frame(Pro_RE22_RE22_moduleTraitCor)
Pro_RE22_RE22_moduleTraitCor_df$mod_names <- row.names(Pro_RE22_RE22_moduleTraitCor_df)
Pro_RE22_RE22_moduleTraitCor_df <- Pro_RE22_RE22_moduleTraitCor_df[,c("mod_names","Condition.Vibrio_coralliilyticus_RE22_exposure_6h.vs.Control_no_treatment")]
Pro_RE22_RE22_moduleTraitPvalue_df <- as.data.frame(Pro_RE22_RE22_moduleTraitPvalue)
Pro_RE22_RE22_moduleTraitPvalue_df$mod_names <- row.names(Pro_RE22_RE22_moduleTraitPvalue_df)
Pro_RE22_RE22_moduleTraitPvalue_df <- Pro_RE22_RE22_moduleTraitPvalue_df[,c("mod_names","Condition.Vibrio_coralliilyticus_RE22_exposure_6h.vs.Control_no_treatment")]
colnames(Pro_RE22_RE22_moduleTraitPvalue_df)[2] <- "pvalue"
Pro_RE22_RE22_moduleTraitCor_Pval_df <- join(Pro_RE22_RE22_moduleTraitCor_df, Pro_RE22_RE22_moduleTraitPvalue_df, by = "mod_names")

# Significantly correlated modules
Pro_RE22_RE22_moduleTraitCor_Pval_df[order(Pro_RE22_RE22_moduleTraitCor_Pval_df$pvalue),]
class(Pro_RE22_RE22_moduleTraitCor_Pval_df$pvalue) # numeric
# subset just for positive associations
Pro_RE22_RE22_moduleTraitCor_Pval_df_sig <- Pro_RE22_RE22_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
Pro_RE22_RE22_moduleTraitCor_Pval_df_sig # 3

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
Pro_RE22_RE22_moduleTraitCor_Pval_df_sig_list <- Pro_RE22_RE22_moduleTraitCor_Pval_df_sig$mod_names
Pro_RE22_RE22_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Pro_RE22_RE22_moduleTraitCor_Pval_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
names(Pro_RE22_RE22_moduleTraitCor_Pval_df_sig_list_rm) <- c("turquoise",           "pink",      "lightcyan", "darkolivegreen")
matrix_common= Pro_RE22_dds_rlog_matrix_common_RE22
moduleColors= Pro_RE22_RE22_moduleColors
lookup =   C_vir_rtracklayer_apop_product_final

lookup_mod_apop <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$ID %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
Pro_RE22_RE22_module_apop <- lapply(Pro_RE22_RE22_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop)
Pro_RE22_RE22_module_apop_df <- do.call(rbind,Pro_RE22_RE22_module_apop)
Pro_RE22_RE22_module_apop_df$mod_names <- gsub("\\..*","",row.names(Pro_RE22_RE22_module_apop_df))
Pro_RE22_RE22_module_apop_df$mod_names <- gsub("^","ME",Pro_RE22_RE22_module_apop_df$mod_names)
# add module significance
Pro_RE22_RE22_module_apop_df <- left_join(Pro_RE22_RE22_module_apop_df,Pro_RE22_RE22_moduleTraitCor_Pval_df_sig)
Pro_RE22_RE22_module_apop_df$exp <- "Pro_RE22_RE22"

## Gene relationship to trait and important modules: Gene Significance and Module Membership
# Define variable injected 
# do later if necessary

## Intramodular analysis: identifying genes with high GS and MM
# later if necessary 

## IDENTIFY HUB GENES IN EACH SIG MODULE ##
Pro_RE22_RE22_colorh = c("MEturquoise", "MEcyan",      "MEred" )

#Pro_RE22_RE22_Module_hub_genes <- chooseTopHubInEachModule(
#  Pro_RE22_dds_rlog_matrix_common_RE22, 
#  Pro_RE22_RE22_colorh, 
#  power = 2,  # power used for the adjacency network
#  type = "signed hybrid", 
#  corFnc = "bicor"
#)
class(Pro_RE22_RE22_Module_hub_genes)
Pro_RE22_RE22_Module_hub_genes_df <- as.data.frame(Pro_RE22_RE22_Module_hub_genes)
colnames(Pro_RE22_RE22_Module_hub_genes_df)[1] <- "ID"
Pro_RE22_RE22_Module_hub_genes <- merge(Pro_RE22_RE22_Module_hub_genes_df, C_vir_rtracklayer, by = "ID")
nrow(Pro_RE22_RE22_Module_hub_genes) 
# no apoptosis hub genes 

#### MASTER C_VIR MODULE-TRAIT-SIG APOP GENES ####
colnames(Dermo_Tol_module_apop_df)[5] <- "mod_signif"
colnames(Dermo_Sus_module_apop_df)[5] <- "mod_signif"
colnames(ROD_Res_module_apop_df)[5] <- "mod_signif"
colnames(ROD_Sus_module_apop_df)[5] <- "mod_signif"
colnames(Probiotic_module_apop_df)[5] <- "mod_signif"
colnames(Pro_RE22_Pro_module_apop_df )[5] <- "mod_signif"    # use dataframes from the separated networks rather than the combined network, this top one is S4
colnames(Pro_RE22_Pro_RI_module_apop_df )[5] <- "mod_signif" # use dataframes from the separated networks rather than the combined network
colnames(Pro_RE22_RE22_module_apop_df)[5] <- "mod_signif" # use dataframes from the separated networks rather than the combined network
C_vir_all_exp_mod_sig_apop <- rbind(Dermo_Tol_module_apop_df,Dermo_Sus_module_apop_df,ROD_Res_module_apop_df,ROD_Sus_module_apop_df,
                                    Probiotic_module_apop_df,Pro_RE22_Pro_module_apop_df ,Pro_RE22_Pro_RI_module_apop_df ,Pro_RE22_RE22_module_apop_df)
nrow(C_vir_all_exp_mod_sig_apop) # 934 # updated August 21st, 2020
C_vir_all_exp_mod_sig_apop_positive <- C_vir_all_exp_mod_sig_apop %>% filter(mod_signif >0)

#### ONE WAY MEASURE MODULE PRESERVATION BETWEEN C. VIRGINICA DIFFERENT EXPERIMENTS ####

# Column names agree between the datasets, so we can go ahead and measure module preservation between all 
## 1. MEASURE OVERLAP BETWEEN RE22 AND RODs

# PRO_RE22 (all, not separated) VS. ROD RES ##
Pro_RE22_ROD_multiExpr = list(Pro_RE22=list(data=Pro_RE22_dds_rlog_matrix_common),ROD_Res=list(data=ROD_Resistant_dds_rlog_matrix_common))
Pro_RE22_multiColor = list(Pro_RE22 = Pro_RE22_moduleColors) 

#Load module preservation calculations
#Pro_RE22_ROD_Res_module_preservation = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_ROD_Res_module_preservation.RData")

#Pro_RE22_mp =modulePreservation(Pro_RE22_ROD_multiExpr, Pro_RE22_multiColor,
#                      referenceNetworks=1,
#                      verbose=3,
#                      networkType="signed hybrid", # use same signed hybrid as before
#                      corFnc = "bicor", # use recommended bicor as before 
#                      nPermutations=100,
#                      randomSeed = 1, # recommended in langfelder tutorial
#                      quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Pro_RE22_stats = Pro_RE22_mp$preservation$Z$ref.Pro_RE22$inColumnsAlsoPresentIn.ROD_Res
Pro_RE22_stats_order <- Pro_RE22_stats[order(-Pro_RE22_stats[,2]),c(1:2)]
Pro_RE22_stats_preserved <- Pro_RE22_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
Pro_RE22_stats_MR <- Pro_RE22_mp$preservation$observed$ref.Pro_RE22$inColumnsAlsoPresentIn.ROD_Res
Pro_RE22_stats_MR <- Pro_RE22_stats_MR [,c("moduleSize","medianRank.pres")]
Pro_RE22_stats_MR_order <- Pro_RE22_stats_MR[order(Pro_RE22_stats_MR[,2]),]
## Moderate preservation between Pro_RE22 and ROD Res
# brown                 1000     7.8915583 (does okay on median rank)
# red                   1000     7.1531975 (not in top twenty rank of module preservation)
# steelblue              175     6.5939398 (steelblue is the top module in the ranking)
# black                 1000     6.3774745 (not in top twenty)
# midnightblue           527     5.7505155 (in top twenty)
# purple                 731     5.7203439 (is in the top twenty)
# blue                  1000     5.1007921 (not in top twenty)
# pink                   787     5.0031835 (pink not in top twenty)

# Save module preservation calculation
# save(Pro_RE22_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_ROD_Res_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(Pro_RE22_stats_MR_order_less20 , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))
# mod_names moduleSize medianRank.pres
# 1    steelblue        175               1
# 2       coral1         37               5
# 3    honeydew1         38               5
# 4 navajowhite2         49               6
# 5       brown4         95               9

#Pro_RE22_stats_preserved_Zsum_medianRank <- merge(Pro_RE22_stats_MR_order_less20, Pro_RE22_stats_preserved)

# Were any of these preserved modules significant in Pro_RE22 or ROD_Res?
Pro_RE22_RI_SIG_ROD_Res_preserved_all <- Pro_RE22_stats_preserved[Pro_RE22_stats_preserved$mod_name %in% Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm,]
Pro_RE22_S4_SIG_ROD_Res_preserved_all <- Pro_RE22_stats_preserved[Pro_RE22_stats_preserved$mod_name %in% Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm,]
# midnightblue, blue, red 
Pro_RE22_RE22_SIG_ROD_Res_preserved_all <- Pro_RE22_stats_preserved[Pro_RE22_stats_preserved$mod_name %in% Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm,]
#  black       1000      6.377475
# 2         blue       1000      5.100792
# 3 midnightblue        527      5.750516
# 4         pink        787      5.003184
# 5       purple        731      5.720344
# 6          red       1000      7.153197
# 7    steelblue        175      6.593940
Pro_RE22_ROD_RES_SIG_preserved_all <- Pro_RE22_stats_preserved[Pro_RE22_stats_preserved$mod_name %in% ROD_Res_moduleTraitCor_Pval_df_sig_list_rm,]
# 0 

# Conclusions: 4 conserved between ROD_Res challenge and Pro_RE22, but are only significant with the RE22 and probiotic challenges. 
  # No overlap in ones significant in both experiments with disease challenge 

## Pro_RE22_RE22 VS. ROD_Res ##
Pro_RE22_RE22_ROD_multiExpr = list(Pro_RE22_RE22=list(data=Pro_RE22_dds_rlog_matrix_common_RE22),ROD_Res=list(data=ROD_Resistant_dds_rlog_matrix_common))
Pro_RE22_RE22_multiColor = list(Pro_RE22_RE22 = Pro_RE22_RE22_moduleColors) 
# Pro_RE22_ROD_mp =modulePreservation(Pro_RE22_RE22_ROD_multiExpr, Pro_RE22_RE22_multiColor,
#                                 referenceNetworks=1,
#                                 verbose=3,
#                                 networkType="signed hybrid", # use same signed hybrid as before
#                                 corFnc = "bicor", # use recommended bicor as before 
#                                 nPermutations=100,
#                                 randomSeed = 1, # recommended in langfelder tutorial
#                                 quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Pro_RE22_RE22_ROD_Res_stats = Pro_RE22_ROD_mp$preservation$Z$ref.Pro_RE22_RE22$inColumnsAlsoPresentIn.ROD_Res
Pro_RE22_RE22_ROD_Res_stats_order <- Pro_RE22_RE22_ROD_Res_stats[order(-Pro_RE22_RE22_ROD_Res_stats[,2]),c(1:2)]
Pro_RE22_RE22_ROD_Res_stats_preserved <- Pro_RE22_RE22_ROD_Res_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
# green module is preserved
Pro_RE22_RE22_ROD_Res_stats_MR <- Pro_RE22_ROD_mp$preservation$observed$ref.Pro_RE22_RE22$inColumnsAlsoPresentIn.ROD_Res
Pro_RE22_RE22_ROD_Res_stats_MR <- Pro_RE22_RE22_ROD_Res_stats_MR [,c("moduleSize","medianRank.pres")]
Pro_RE22_RE22_ROD_Res_stats_MR_order <-Pro_RE22_RE22_ROD_Res_stats_MR[order(-Pro_RE22_RE22_ROD_Res_stats_MR [,2]),]

# Save module preservation calculation
# save(Pro_RE22_ROD_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_RE22_ROD_Res_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(Pro_RE22_RE22_ROD_Res_stats_MR_order, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

# were any of the preserved modules also significant?
Pro_RE22_RE22_only_SIG_ROD_Res_preserved <-Pro_RE22_RE22_ROD_Res_stats_preserved [Pro_RE22_RE22_ROD_Res_stats_preserved $mod_name %in% Pro_RE22_RE22_moduleTraitCor_Pval_df_sig_list_rm,]
# 0
Pro_RE22_RE22_only_ROD_RES_SIG_preserved <- Pro_RE22_RE22_ROD_Res_stats_preserved [Pro_RE22_RE22_ROD_Res_stats_preserved $mod_name %in% ROD_Res_moduleTraitCor_Pval_df_sig_list_rm,]
# 0



##  ROD_Res VS. Pro_RE22_RE22 REVERSED 
Pro_RE22_RE22_ROD_RES_multiExpr = list(Pro_RE22_RE22=list(data=Pro_RE22_dds_rlog_matrix_common_RE22),ROD_Res=list(data=ROD_Resistant_dds_rlog_matrix_common))
# Names for the two sets
setLabels = c("Pro_RE22_RE22", "ROD_Res")
# Important: components of multiExpr must carry identificating names
names(Pro_RE22_RE22_ROD_RES_multiExpr) = setLabels
# Create an object (list) holding the module labels for each set:
colorList = list(Pro_RE22_RE22_moduleColors, ROD_Res_moduleColors)
# Components of the list must be named so that the names can be matched to the names of multiExpr
names(colorList) = setLabels
Pro_RE22_RE22_multiColor = list(Pro_RE22_RE22 = Pro_RE22_RE22_moduleColors) 
Pro_RE22_ROD_Res_BOTH_mp = modulePreservation(Pro_RE22_RE22_ROD_RES_multiExpr, colorList ,
                                              referenceNetworks=c(1:2),
                                              verbose=3,
                                              networkType="signed hybrid", # use same signed hybrid as before
                                              corFnc = "bicor", # use recommended bicor as before 
                                              nPermutations=100,
                                              randomSeed = 1, # recommended in langfelder tutorial
                                              quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
# save results
save(Pro_RE22_ROD_Res_BOTH_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_ROD_Res_BOTH_mp_network.RData")

Pro_RE22_RE22_ROD_Res_BOTHstats = Pro_RE22_ROD_mp$preservation$Z$ref.Pro_RE22_RE22$inColumnsAlsoPresentIn.ROD_Res
Pro_RE22_RE22_ROD_Res_BOTHstats_order <- Pro_RE22_RE22_ROD_Res_stats[order(-Pro_RE22_RE22_ROD_Res_stats[,2]),c(1:2)]
Pro_RE22_RE22_ROD_Res_BOTHstats_preserved <- Pro_RE22_RE22_ROD_Res_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
# green module is preserved
Pro_RE22_RE22_ROD_Res_stats_MR <- Pro_RE22_ROD_mp$preservation$observed$ref.Pro_RE22_RE22$inColumnsAlsoPresentIn.ROD_Res
Pro_RE22_RE22_ROD_Res_stats_MR <- Pro_RE22_RE22_ROD_Res_stats_MR [,c("moduleSize","medianRank.pres")]
Pro_RE22_RE22_ROD_Res_stats_MR_order <-Pro_RE22_RE22_ROD_Res_stats_MR[order(-Pro_RE22_RE22_ROD_Res_stats_MR [,2]),]

# Save module preservation calculation
# save(Pro_RE22_ROD_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_RE22_ROD_Res_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(Pro_RE22_RE22_ROD_Res_stats_MR_order, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

# were any of the preserved modules also significant?
Pro_RE22_RE22_only_SIG_ROD_Res_preserved <-Pro_RE22_RE22_ROD_Res_stats_preserved [Pro_RE22_RE22_ROD_Res_stats_preserved $mod_name %in% Pro_RE22_RE22_moduleTraitCor_Pval_df_sig_list_rm,]
# 0
Pro_RE22_RE22_only_ROD_RES_SIG_preserved <- Pro_RE22_RE22_ROD_Res_stats_preserved [Pro_RE22_RE22_ROD_Res_stats_preserved $mod_name %in% ROD_Res_moduleTraitCor_Pval_df_sig_list_rm,]
# 0

## Pro_RE22 (all not separated) vs. ROD_Sus modules ##
 Pro_RE22_Sus_multiExpr = list(Pro_RE22=list(data=Pro_RE22_dds_rlog_matrix_common),ROD_Sus=list(data=ROD_Susceptible_dds_rlog_matrix_common))

#Load module preservation calculations
#Pro_RE22_ROD_Sus_module_preservation = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_ROD_Sus_module_preservation.RData")

# Pro_RE22_Sus_mp =modulePreservation(Pro_RE22_Sus_multiExpr, Pro_RE22_multiColor,
#                                 referenceNetworks=1,
#                                 verbose=3,
#                                 networkType="signed hybrid", # use same signed hybrid as before
#                                 corFnc = "bicor", # use recommended bicor as before 
#                                 nPermutations=100,
#                                 randomSeed = 1, # recommended in langfelder tutorial
#                                 quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Pro_RE22_Sus_stats = Pro_RE22_Sus_mp$preservation$Z$ref.Pro_RE22$inColumnsAlsoPresentIn.ROD_Sus 
Pro_RE22_Sus_stats_order <- Pro_RE22_Sus_stats[order(-Pro_RE22_Sus_stats[,2]),c(1:2)]
Pro_RE22_Sus_stats_preserved <- Pro_RE22_Sus_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
#mod_name moduleSize Zsummary.pres
#1      red       1000      7.327336 (good median rank)
#2     blue       1000      5.772638 (good median rank)
#3    brown       1000      5.070500 (bad median rank)

Pro_RE22_Sus_stats_MR <- Pro_RE22_Sus_mp$preservation$observed$ref.Pro_RE22$inColumnsAlsoPresentIn.ROD_Sus
Pro_RE22_Sus_stats_MR <- Pro_RE22_Sus_stats_MR [,c("moduleSize","medianRank.pres")]
Pro_RE22_Sus_stats_MR_order <- Pro_RE22_Sus_stats_MR[order(Pro_RE22_Sus_stats_MR[,2]),]
Pro_RE22_Sus_stats_MR_order_less20 <- Pro_RE22_Sus_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save module preservation calculation
#save(Pro_RE22_Sus_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_ROD_Sus_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(Pro_RE22_Sus_stats_MR_order_less20 , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
 geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

Pro_RE22_Sus_stats_preserved_Zsum_medianRank <- merge(Pro_RE22_Sus_stats_MR_order_less20, Pro_RE22_Sus_stats_preserved)
# blue, red (brown )

Pro_RE22_RI_SIG_ROD_Sus_preserved_all <- Pro_RE22_Sus_stats_preserved[Pro_RE22_Sus_stats_preserved$mod_name %in% Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm,]
# 1      red       1000      7.327336
# 2     blue       1000      5.772638
Pro_RE22_S4_SIG_ROD_Sus_preserved_all <- Pro_RE22_Sus_stats_preserved[Pro_RE22_Sus_stats_preserved$mod_name %in%  Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm,]
# 1      red       1000      7.327336
# 2     blue       1000      5.772638
Pro_RE22_RE22_SIG_ROD_Sus_preserved_all <- Pro_RE22_Sus_stats_preserved[Pro_RE22_Sus_stats_preserved$mod_name %in%  Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm,]
# 1      red       1000      7.327336
# 2     blue       1000      5.772638
Pro_RE22_ROD_SUS_SIG_preserved_all <- Pro_RE22_Sus_stats_preserved[Pro_RE22_Sus_stats_preserved$mod_name %in%  ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
# brown

### Assess whether Pro_RE22_RE22 vs. ROD_Sus modules are preserved ##
Pro_RE22_RE22_ROD_Sus_multiExpr = list(Pro_RE22_RE22=list(data=Pro_RE22_dds_rlog_matrix_common_RE22),ROD_Sus=list(data=ROD_Susceptible_dds_rlog_matrix_common))
Pro_RE22_RE22_multiColor = list(Pro_RE22_RE22 = Pro_RE22_RE22_moduleColors) 
# Pro_RE22_ROD_Sus_mp =modulePreservation(Pro_RE22_RE22_ROD_Sus_multiExpr, Pro_RE22_RE22_multiColor,
#                                    referenceNetworks=1,
#                                    verbose=3,
#                                    networkType="signed hybrid", # use same signed hybrid as before
#                                    corFnc = "bicor", # use recommended bicor as before 
#                                    nPermutations=100,
#                                    randomSeed = 1, # recommended in langfelder tutorial
#                                    quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Pro_RE22_ROD_Sus_stats = Pro_RE22_ROD_Sus_mp$preservation$Z$ref.Pro_RE22_RE22$inColumnsAlsoPresentIn.ROD_Sus
Pro_RE22_ROD_Sus_stats_order <- Pro_RE22_ROD_Sus_stats[order(-Pro_RE22_ROD_Sus_stats[,2]),c(1:2)]
Pro_RE22_ROD_Sus_stats_preserved <- Pro_RE22_ROD_Sus_stats_order  %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
# 0 modules preserved between Pro_RE22_RE22 and the ROD Susceptible 
Pro_RE22_ROD_Sus_stats_MR <- Pro_RE22_ROD_Sus_mp$preservation$observed$ref.Pro_RE22_RE22$inColumnsAlsoPresentIn.ROD_Sus
Pro_RE22_ROD_Sus_stats_MR <- Pro_RE22_ROD_Sus_stats_MR  [,c("moduleSize","medianRank.pres")]
Pro_RE22_ROD_Sus_stats_MR_order <- Pro_RE22_ROD_Sus_stats_MR [order(Pro_RE22_ROD_Sus_stats_MR [,2]),]

# Save network
# save(Pro_RE22_ROD_Sus_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_RE22_ROD_Sus_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(Pro_RE22_ROD_Sus_stats_MR_order, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

# were any of the preserved modules also significant?
Pro_RE22_RE22_only_SIG_ROD_Sus_preserved <- Pro_RE22_ROD_Sus_stats_preserved[Pro_RE22_ROD_Sus_stats_preserved$mod_name %in% Pro_RE22_RE22_moduleTraitCor_Pval_df_sig_list_rm,]
#0 
Pro_RE22_RE22_only_ROD_Sus_SIG_preserved <- Pro_RE22_ROD_Sus_stats_preserved[Pro_RE22_ROD_Sus_stats_preserved$mod_name %in% ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
#0

## ROD_Sus vs. Pro_RE22_RE22 ##
Pro_RE22_RE22_ROD_Sus_rev_multiExpr = list(ROD_Sus=list(data=ROD_Susceptible_dds_rlog_matrix_common), Pro_RE22_RE22=list(data=Pro_RE22_dds_rlog_matrix_common_RE22))
ROD_Sus_rev_multiColor = list(ROD_Sus = ROD_Sus_moduleColors) 
#Pro_RE22_ROD_Sus_rev_mp =modulePreservation(Pro_RE22_RE22_ROD_Sus_rev_multiExpr, ROD_Sus_rev_multiColor,
#                                  referenceNetworks=1,
#                                  verbose=3,
#                                  networkType="signed hybrid", # use same signed hybrid as before
#                                  corFnc = "bicor", # use recommended bicor as before 
#                                  nPermutations=100,
#                                  randomSeed = 1, # recommended in langfelder tutorial
#                                  quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Pro_RE22_ROD_Sus_rev_stats = Pro_RE22_ROD_Sus_rev_mp$preservation$Z$ref.ROD_Sus$inColumnsAlsoPresentIn.Pro_RE22_RE22
Pro_RE22_ROD_Sus_rev_stats_order <- Pro_RE22_ROD_Sus_rev_stats[order(-Pro_RE22_ROD_Sus_rev_stats[,2]),c(1:2)]
Pro_RE22_ROD_Sus_rev_stats_preserved <- Pro_RE22_ROD_Sus_rev_stats_order  %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
#0 preserved
Pro_RE22_ROD_Sus_rev_stats_MR <- Pro_RE22_ROD_Sus_rev_mp$preservation$observed$ref.ROD_Sus$inColumnsAlsoPresentIn.Pro_RE22_RE22
Pro_RE22_ROD_Sus_rev_stats_MR <- Pro_RE22_ROD_Sus_rev_stats_MR  [,c("moduleSize","medianRank.pres")]
Pro_RE22_ROD_Sus_rev_stats_MR_order <- Pro_RE22_ROD_Sus_rev_stats_MR [order(Pro_RE22_ROD_Sus_rev_stats_MR[,2]),]

# Save network
#save(Pro_RE22_ROD_Sus_rev_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_RE22_ROD_Sus_rev_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(Pro_RE22_ROD_Sus_rev_stats_MR_order, aes(x= moduleSize, y = medianRank.pres, color = row.names(Pro_RE22_ROD_Sus_rev_stats_MR_order))) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank")

# were any of the preserved modules also significant?
Pro_RE22_RE22_only_SIG_ROD_Sus_rev_preserved <- Pro_RE22_ROD_Sus_rev_stats_preserved[Pro_RE22_ROD_Sus_rev_stats_preserved$mod_name %in% Pro_RE22_RE22_moduleTraitCor_Pval_df_sig_list_rm,]
#0
Pro_RE22_RE22_only_ROD_Sus_SIG_rev_preserved <- Pro_RE22_ROD_Sus_rev_stats_preserved[Pro_RE22_ROD_Sus_rev_stats_preserved$mod_name %in% ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
#0 

## Assess whether Pro_RE22 (all not separated) and Dermo_Tol modules are preserved ##
Pro_RE22_Dermo_multiExpr = list(Pro_RE22=list(data=Pro_RE22_dds_rlog_matrix_common),Dermo_Tol=list(data=Dermo_Tolerant_dds_vst_matrix_common))
Pro_RE22_multiColor = list(Pro_RE22 = Pro_RE22_moduleColors) 
# Load module preservation calculations
#Pro_RE22_Dermo_Tol = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_Dermo_Tol_module_preservation.RData")
# Pro_RE22_Dermo_mp =modulePreservation(Pro_RE22_Dermo_multiExpr, Pro_RE22_multiColor,
#                                 referenceNetworks=1,
#                                 verbose=3,
#                                 networkType="signed hybrid", # use same signed hybrid as before
#                                 corFnc = "bicor", # use recommended bicor as before 
#                                 nPermutations=100,
#                                 randomSeed = 1, # recommended in langfelder tutorial
#                                 quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Pro_RE22_Dermo_stats = Pro_RE22_Dermo_mp$preservation$Z$ref.Pro_RE22$inColumnsAlsoPresentIn.Dermo_Tol
Pro_RE22_Dermo_stats_order <- Pro_RE22_Dermo_stats[order(-Pro_RE22_Dermo_stats[,2]),c(1:2)]
Pro_RE22_Dermo_stats_preserved <- Pro_RE22_Dermo_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
# mod_name moduleSize Zsummary.pres
# 1        salmon        615     10.598431 (high median rank)
# 2        purple        731      9.399914 (decent median rank)
# 3     steelblue        175      9.376713 (high median rank)
# 4         brown       1000      9.373572 (decent median rank)
# 5           red       1000      8.038084 (lower median rank)
# 6     turquoise       1000      7.472569 (lower median rank)
# 7        yellow       1000      7.446067 (lower median rank)
# 8          blue       1000      6.243046 (lower median rank)
# 9         green       1000      5.794057 (lower median rank)
# 10    lightcyan        474      5.753286 (lower median rank)
# 11 midnightblue        527      5.580933 (lower median rank)
# 12        black       1000      5.354414 (lowest median rank in top twentry)

Pro_RE22_Dermo_stats_MR <- Pro_RE22_Dermo_mp$preservation$observed$ref.Pro_RE22$inColumnsAlsoPresentIn.Dermo_Tol
Pro_RE22_Dermo_stats_MR <- Pro_RE22_Dermo_stats_MR [,c("moduleSize","medianRank.pres")]
Pro_RE22_Dermo_stats_MR_order <- Pro_RE22_Dermo_stats_MR[order(Pro_RE22_Dermo_stats_MR[,2]),]
Pro_RE22_Dermo_stats_MR_order_less20 <- Pro_RE22_Dermo_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save network
#save(Pro_RE22_Dermo_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_Dermo_Tol_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(Pro_RE22_Dermo_stats_MR_order_less20 , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
 geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

Pro_RE22_Dermo_stats_preserved_Zsum_medianRank <- merge(Pro_RE22_Dermo_stats_MR_order_less20, Pro_RE22_Dermo_stats_preserved)
# 12 modules

Pro_RE22_RI_SIG_Dermo_Tol_preserved_all <- Pro_RE22_Dermo_stats_preserved[Pro_RE22_Dermo_stats_preserved$mod_names %in% Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm,]
#1        black       1000              19      5.354414
#2         blue       1000              13      6.243046
#3        green       1000              18      5.794057
#4 midnightblue        527              13      5.580933
#5       purple        731              11      9.399914
#6          red       1000              15      8.038084
#7       salmon        615               4     10.598431
#8    steelblue        175               1      9.376713
#9       yellow       1000              14      7.446067
Pro_RE22_S4_SIG_Dermo_Tol_preserved_all <- Pro_RE22_Dermo_stats_preserved[Pro_RE22_Dermo_stats_preserved$mod_names %in% Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm,]
#mod_name moduleSize medianRank.pres Zsummary.pres
#1         blue       1000              13      6.243046
#2        green       1000              18      5.794057
#3 midnightblue        527              13      5.580933
#4          red       1000              15      8.038084
Pro_RE22_RE22_SIG_RDermo_Tol_preserved_all <- Pro_RE22_Dermo_stats_preserved[Pro_RE22_Dermo_stats_preserved$mod_names %in% Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm,]
#       mod_name moduleSize medianRank.pres Zsummary.pres
# 1        black       1000              19      5.354414
# 2         blue       1000              13      6.243046
# 3        green       1000              18      5.794057
# 4 midnightblue        527              13      5.580933
# 5       purple        731              11      9.399914
# 6          red       1000              15      8.038084
# 7       salmon        615               4     10.598431
# 8    steelblue        175               1      9.376713
# 9       yellow       1000              14      7.446067
Pro_RE22_Dermo_Tol_SIG_preserved_all <- Pro_RE22_Dermo_stats_preserved[Pro_RE22_Dermo_stats_preserved$mod_names %in% Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm,]
# mod_name moduleSize medianRank.pres Zsummary.pres
# 1 turquoise       1000              14      7.472569
# 2    yellow       1000              14      7.446067

# Conclusion: probiotic and Dermo tolerant share two modules significantly correlated with disease challenge 

##  Pro_RE22_RE22 vs Dermo_Tol ##
Pro_RE22_RE22_Dermo_multiExpr = list(Pro_RE22_RE22=list(data=Pro_RE22_dds_rlog_matrix_common_RE22), Dermo_Tol=list(data=Dermo_Tolerant_dds_vst_matrix_common))
Pro_RE22_multiColor = list(Dermo_Tol = Pro_RE22_RE22_moduleColors) 
#Pro_RE22_RE22_Dermo_mp =modulePreservation(Pro_RE22_RE22_Dermo_rev_multiExpr, Pro_RE22_multiColor,
#                                               referenceNetworks=1,
#                                               verbose=3,
#                                               networkType="signed hybrid", # use same signed hybrid as before
#                                               corFnc = "bicor", # use recommended bicor as before 
#                                               nPermutations=100,
#                                               randomSeed = 1, # recommended in langfelder tutorial
#                                               quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Pro_RE22_RE22_Dermo_stats = Pro_RE22_RE22_Dermo_mp$preservation$Z$ref.Pro_RE22_RE22$inColumnsAlsoPresentIn.Dermo_Tol
Pro_RE22_RE22_Dermo_stats_order <- Pro_RE22_RE22_Dermo_stats[order(-Pro_RE22_RE22_Dermo_rev_stats[,2]),c(1:2)]
Pro_RE22_RE22_Dermo_stats_preserved <- Pro_RE22_RE22_Dermo_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
# mod_name moduleSize Zsummary.pres
# 1    green       1000      6.645777
# 2    black       1000      5.979723
Pro_RE22_RE22_Dermo_stats_MR <- Pro_RE22_RE22_Dermo_mp$preservation$observed$ref.Pro_RE22_RE22$inColumnsAlsoPresentIn.Dermo_Tol
Pro_RE22_RE22_Dermo_stats_MR <- Pro_RE22_RE22_Dermo_stats_MR[,c("moduleSize","medianRank.pres")]
Pro_RE22_RE22_Dermo_stats_MR_order <- Pro_RE22_RE22_Dermo_stats_MR[order(-Pro_RE22_RE22_Dermo_rev_stats_MR[,2]),]
Pro_RE22_RE22_Dermo_stats_MR_order_less20 <- Pro_RE22_RE22_Dermo_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save network
#save(Pro_RE22_RE22_Dermo_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_RE22_Dermo_Tol_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(Pro_RE22_RE22_Dermo_stats_MR_order_less20 , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

Pro_RE22_RE22_Dermo_Tol_SIG_RE22_preserved_all <- Pro_RE22_RE22_Dermo_stats_preserved[Pro_RE22_RE22_Dermo_stats_preserved$mod_name %in% Pro_RE22_RE22_moduleTraitCor_Pval_df_sig_list_rm,]
#0
Pro_RE22_RE22_Dermo_Tol_SIG_preserved_all <- Pro_RE22_RE22_Dermo_stats_preserved[Pro_RE22_RE22_Dermo_stats_preserved$mod_name %in% Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm,]
#0

# Conclusion: 

## Dermo_Tol vs Pro_RE22_RE22 ##
Pro_RE22_RE22_Dermo_rev_multiExpr = list(Dermo_Tol=list(data=Dermo_Tolerant_dds_vst_matrix_common), Pro_RE22_RE22=list(data=Pro_RE22_dds_rlog_matrix_common_RE22))
Dermo_Tol_multiColor = list(Dermo_Tol = Pro_RE22_RE22_moduleColors) 
#Pro_RE22_RE22_Dermo_rev_mp =modulePreservation(Pro_RE22_RE22_Dermo_rev_multiExpr, Dermo_Tol_multiColor,
#                                      referenceNetworks=1,
#                                      verbose=3,
#                                      networkType="signed hybrid", # use same signed hybrid as before
#                                      corFnc = "bicor", # use recommended bicor as before 
#                                      nPermutations=100,
#                                      randomSeed = 1, # recommended in langfelder tutorial
#                                      quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Pro_RE22_RE22_Dermo_rev_stats = Pro_RE22_RE22_Dermo_rev_mp$preservation$Z$ref.Dermo_Tol$inColumnsAlsoPresentIn.Pro_RE22_RE22
Pro_RE22_RE22_Dermo_rev_stats_order <- Pro_RE22_RE22_Dermo_rev_stats[order(-Pro_RE22_RE22_Dermo_rev_stats[,2]),c(1:2)]
Pro_RE22_RE22_Dermo_rev_stats_preserved <- Pro_RE22_RE22_Dermo_rev_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
# 39 preserved modules..27 as highly preserved
# green and black modules both highly preserved
Pro_RE22_RE22_Dermo_rev_stats_MR <- Pro_RE22_RE22_Dermo_rev_mp$preservation$observed$ref.Dermo_Tol$inColumnsAlsoPresentIn.Pro_RE22_RE22
Pro_RE22_RE22_Dermo_rev_stats_MR <- Pro_RE22_RE22_Dermo_rev_stats_MR[,c("moduleSize","medianRank.pres")]
Pro_RE22_RE22_Dermo_rev_stats_MR_order <- Pro_RE22_RE22_Dermo_rev_stats_MR[order(-Pro_RE22_RE22_Dermo_rev_stats_MR[,2]),]
Pro_RE22_RE22_Dermo_rev_stats_MR_order_less20 <- Pro_RE22_RE22_Dermo_rev_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save network
#save(Pro_RE22_RE22_Dermo_rev_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_RE22_Dermo_Tol_rev_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(Pro_RE22_RE22_Dermo_rev_stats_MR_order_less20 , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

Pro_RE22_RE22_Dermo_Tol_SIG_RE22_rev_preserved_all <- Pro_RE22_RE22_Dermo_rev_stats_preserved[Pro_RE22_RE22_Dermo_rev_stats_preserved$mod_name %in% Pro_RE22_RE22_moduleTraitCor_Pval_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#2       lightcyan        802      29.57142
#3            pink       1000      27.87111
#4       turquoise       1000      26.66202
#31 darkolivegreen        244      10.19548
Pro_RE22_RE22_Dermo_Tol_SIG_rev_preserved_all <- Pro_RE22_RE22_Dermo_rev_stats_preserved[Pro_RE22_RE22_Dermo_rev_stats_preserved$mod_name %in% Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#3        pink       1000     27.871115
#4   turquoise       1000     26.662016
#11     yellow       1000     17.154307
#13        tan        871     16.070944
#20 lightgreen        750     13.307754
#23  royalblue        711     12.804810
#39 orangered4         65      5.341695

# Assess whether Pro_RE22 (all, not separated) and Dermo_Sus modules are preserved ##
Pro_RE22_Dermo_Sus_multiExpr = list(Pro_RE22=list(data=Pro_RE22_dds_rlog_matrix_common),Dermo_Sus=list(data=Dermo_Susceptible_dds_vst_matrix_common))
Pro_RE22_multiColor = list(Pro_RE22 = Pro_RE22_moduleColors) 
# Load module preservation calculations
#Pro_RE22_Dermo_Sus = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_Dermo_Sus_module_preservation.RData")
# Pro_RE22_Dermo_Sus_mp =modulePreservation(Pro_RE22_Dermo_Sus_multiExpr, Pro_RE22_multiColor,
#                                       referenceNetworks=1,
#                                       verbose=3,
#                                       networkType="signed hybrid", # use same signed hybrid as before
#                                       corFnc = "bicor", # use recommended bicor as before 
#                                       nPermutations=100,
#                                       randomSeed = 1, # recommended in langfelder tutorial
#                                       quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Pro_RE22_Dermo_Sus_stats = Pro_RE22_Dermo_Sus_mp$preservation$Z$ref.Pro_RE22$inColumnsAlsoPresentIn.Dermo_Sus
Pro_RE22_Dermo_Sus_stats_order <- Pro_RE22_Dermo_Sus_stats[order(-Pro_RE22_Dermo_Sus_stats[,2]),c(1:2)]
Pro_RE22_Dermo_Sus_stats_preserved <- Pro_RE22_Dermo_Sus_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
# mod_name moduleSize Zsummary.pres
# 1    salmon        615      9.461935 (high median rank)
# 2 steelblue        175      7.608015 (highest median rank)
# 3    purple        731      6.821283 (good median rank)
# 4      blue       1000      6.176316 (medium median rank)
# 5       red       1000      6.060967 (medium median rank)
# 6 lightcyan        474      6.010655 (lower median rank)
# 7     brown       1000      5.165484 (not in top twentry)
Pro_RE22_Dermo_Sus_stats_MR <- Pro_RE22_Dermo_Sus_mp$preservation$observed$ref.Pro_RE22$inColumnsAlsoPresentIn.Dermo_Sus
Pro_RE22_Dermo_Sus_stats_MR <- Pro_RE22_Dermo_Sus_stats_MR [,c("moduleSize","medianRank.pres")]
Pro_RE22_Dermo_Sus_stats_MR_order <- Pro_RE22_Dermo_Sus_stats_MR[order(Pro_RE22_Dermo_Sus_stats_MR[,2]),]
Pro_RE22_Dermo_Sus_stats_MR_order_less20 <- Pro_RE22_Dermo_Sus_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save network
#save(Pro_RE22_Dermo_Sus_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_Dermo_Sus_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(Pro_RE22_Dermo_Sus_stats_MR_order_less20 , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

Pro_RE22_Dermo_Sus_stats_preserved_Zsum_medianRank <- Pro_RE22_Dermo_Sus_stats_MR_order_less20[Pro_RE22_Dermo_Sus_stats_MR_order_less20$mod_names %in% Pro_RE22_Dermo_Sus_stats_preserved,]
# 6 modules

# Were any of these preserved modules significant in Pro_RE22 or Dermo_Sus

Pro_RE22_RI_SIG_Dermo_Sus_preserved_all <- Pro_RE22_Dermo_Sus_stats_preserved[Pro_RE22_Dermo_Sus_stats_preserved$mod_name %in% Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm,]
# mod_name moduleSize Zsummary.pres
# 1      blue       1000      6.176316
# 2    purple        731      6.821283
# 3       red       1000      6.060967
# 4    salmon        615      9.461935
# 5 steelblue        175      7.608015
Pro_RE22_S4_SIG_Dermo_Sus_preserved_all <-Pro_RE22_Dermo_Sus_stats_preserved[Pro_RE22_Dermo_Sus_stats_preserved$mod_name %in% Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm,]
# blue, red
Pro_RE22_RE22_SIG_Dermo_Sus_preserved_all <-Pro_RE22_Dermo_Sus_stats_preserved[Pro_RE22_Dermo_Sus_stats_preserved$mod_name %in% Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm,]
# mod_name moduleSize Zsummary.pres
# 1      blue       1000      6.176316
# 2    purple        731      6.821283
# 3       red       1000      6.060967
# 4    salmon        615      9.461935
# 5 steelblue        175      7.608015
Pro_RE22_Dermo_Sus_SIG_preserved_all <-Pro_RE22_Dermo_Sus_stats_preserved[Pro_RE22_Dermo_Sus_stats_preserved$mod_name %in% Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
# 0 

# Assess whether Pro_RE22_RE22 and Dermo_Sus modules are preserved ##
Pro_RE22_RE22_Dermo_Sus_multiExpr = list(Pro_RE22_RE22=list(data=Pro_RE22_dds_rlog_matrix_common_RE22),Dermo_Sus=list(data=Dermo_Susceptible_dds_vst_matrix_common))
Pro_RE22_RE22_multiColor = list(Pro_RE22_RE22 = Pro_RE22_RE22_moduleColors) 
#Pro_RE22_RE22_Dermo_Sus_mp =modulePreservation(Pro_RE22_RE22_Dermo_Sus_multiExpr , Pro_RE22_RE22_multiColor,
#                                          referenceNetworks=1,
#                                          verbose=3,
#                                          networkType="signed hybrid", # use same signed hybrid as before
#                                          corFnc = "bicor", # use recommended bicor as before 
#                                          nPermutations=100,
#                                          randomSeed = 1, # recommended in langfelder tutorial
#                                          quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Pro_RE22_RE22_Dermo_Sus_stats = Pro_RE22_RE22_Dermo_Sus_mp$preservation$Z$ref.Pro_RE22_RE22$inColumnsAlsoPresentIn.Dermo_Sus
Pro_RE22_RE22_Dermo_Sus_stats_order <- Pro_RE22_RE22_Dermo_Sus_stats[order(-Pro_RE22_RE22_Dermo_Sus_stats[,2]),c(1:2)]
Pro_RE22_RE22_Dermo_Sus_stats_preserved <- Pro_RE22_RE22_Dermo_Sus_stats_order  %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
#mod_name moduleSize Zsummary.pres
#1    black       1000      5.182457
Pro_RE22_RE22_Dermo_Sus_stats_MR <- Pro_RE22_RE22_Dermo_Sus_mp$preservation$observed$ref.Pro_RE22_RE22$inColumnsAlsoPresentIn.Dermo_Sus
Pro_RE22_RE22_Dermo_Sus_stats_MR <- Pro_RE22_RE22_Dermo_Sus_stats_MR [,c("moduleSize","medianRank.pres")]
Pro_RE22_RE22_Dermo_Sus_stats_MR_order <- Pro_RE22_RE22_Dermo_Sus_stats_MR[order(-Pro_RE22_RE22_Dermo_Sus_stats_MR[,2]),]
Pro_RE22_RE22_Dermo_Sus_stats_MR_order_less20 <- Pro_RE22_RE22_Dermo_Sus_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save network
#save(Pro_RE22_RE22_Dermo_Sus_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_RE22_Dermo_Sus_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(Pro_RE22_RE22_Dermo_Sus_stats_MR_order_less20, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

# Were any of these preserved modules significant in Pro_RE22 or Dermo_Sus
Pro_RE22_RE22_SIG_Dermo_Sus_preserved_all <-Pro_RE22_RE22_Dermo_Sus_stats_preserved[Pro_RE22_RE22_Dermo_Sus_stats_preserved$mod_name %in% Pro_RE22_RE22_moduleTraitCor_Pval_df_sig_list_rm,]
# 0
Pro_RE22_RE22_Dermo_Sus_SIG_preserved_all <-Pro_RE22_RE22_Dermo_Sus_stats_preserved[Pro_RE22_RE22_Dermo_Sus_stats_preserved$mod_name %in% Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
#0

## Dermo_Sus vs.Pro_RE22_RE22 ##
Pro_RE22_RE22_Dermo_Sus_rev_multiExpr = list(Dermo_Sus=list(data=Dermo_Susceptible_dds_vst_matrix_common),Pro_RE22_RE22=list(data=Pro_RE22_dds_rlog_matrix_common_RE22))
Dermo_Sus_multiColor = list(Dermo_Sus = Dermo_Sus_moduleColors) 
#Pro_RE22_RE22_Dermo_Sus_rev_mp =modulePreservation(Pro_RE22_RE22_Dermo_Sus_rev_multiExpr, Dermo_Sus_multiColor,
#                                          referenceNetworks=1,
#                                          verbose=3,
#                                          networkType="signed hybrid", # use same signed hybrid as before
#                                          corFnc = "bicor", # use recommended bicor as before 
#                                          nPermutations=100,
#                                          randomSeed = 1, # recommended in langfelder tutorial
#                                          quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Pro_RE22_RE22_Dermo_Sus_rev_stats = Pro_RE22_RE22_Dermo_Sus_rev_mp$preservation$Z$ref.Dermo_Sus$inColumnsAlsoPresentIn.Pro_RE22_RE22
Pro_RE22_RE22_Dermo_Sus_rev_stats_order <- Pro_RE22_RE22_Dermo_Sus_rev_stats[order(-Pro_RE22_RE22_Dermo_Sus_rev_stats[,2]),c(1:2)]
Pro_RE22_RE22_Dermo_Sus_rev_stats_preserved <- Pro_RE22_RE22_Dermo_Sus_rev_stats_order  %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
Pro_RE22_RE22_Dermo_Sus_rev_stats_MR <- Pro_RE22_RE22_Dermo_Sus_rev_mp$preservation$observed$ref.Dermo_Sus$inColumnsAlsoPresentIn.Pro_RE22_RE22
Pro_RE22_RE22_Dermo_Sus_rev_stats_MR <- Pro_RE22_RE22_Dermo_Sus_rev_stats_MR[,c("moduleSize","medianRank.pres")]
Pro_RE22_RE22_Dermo_Sus_rev_stats_MR_order <- Pro_RE22_RE22_Dermo_Sus_rev_stats_MR[order(-Pro_RE22_RE22_Dermo_Sus_rev_stats_MR[,2]),]
Pro_RE22_RE22_Dermo_Sus_rev_stats_MR_order_less20 <- Pro_RE22_RE22_Dermo_Sus_rev_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save network
#save(Pro_RE22_RE22_Dermo_Sus_rev_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_RE22_Dermo_Sus_module_rev_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(Pro_RE22_RE22_Dermo_Sus_rev_stats_MR_order_less20, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

# Were any of these preserved modules significant in Pro_RE22 or Dermo_Sus
Pro_RE22_RE22_SIG_Dermo_Sus_rev_preserved_all <-Pro_RE22_RE22_Dermo_Sus_rev_stats_preserved[Pro_RE22_RE22_Dermo_Sus_rev_stats_preserved$mod_name %in% Pro_RE22_RE22_moduleTraitCor_Pval_df_sig_list_rm,]
#   mod_name moduleSize Zsummary.pres
# 3     pink        557      15.67606
Pro_RE22_RE22_Dermo_Sus_SIG_rev_preserved_all <-Pro_RE22_RE22_Dermo_Sus_rev_stats_preserved[Pro_RE22_RE22_Dermo_Sus_rev_stats_preserved$mod_name %in% Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
# 0 

## Assess whether Dermo_Tol and ROD_Res modules are preserved ## 
ROD_Dermo_multiExpr = list(Dermo_Tol=list(data=Dermo_Tolerant_dds_vst_matrix_common), ROD_Res=list(data=ROD_Resistant_dds_rlog_matrix_common))
ROD_Dermo_multiColor = list(Dermo_Tol = Dermo_Tol_moduleColors)
# Load module preservation
#Dermo_Tol_ROD_Res = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Tol_ROD_Res_module_preservation.RData")

#ROD_Dermo_mp =modulePreservation(ROD_Dermo_multiExpr, ROD_Dermo_multiColor,
#                                 referenceNetworks=1,
#                                 verbose=3,
#                                 networkType="signed hybrid", # use same signed hybrid as before
#                                 corFnc = "bicor", # use recommended bicor as before 
#                                 nPermutations=100,
#                                 randomSeed = 1, # recommended in langfelder tutorial
#                                 quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
ROD_Dermo_stats = ROD_Dermo_mp$preservation$Z$ref.Dermo_Tol$inColumnsAlsoPresentIn.ROD_Res
ROD_Dermo_stats_order <- ROD_Dermo_stats[order(-ROD_Dermo_stats[,2]),c(1:2)]
ROD_Dermo_stats_preserved <- ROD_Dermo_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
# mod_name moduleSize Zsummary.pres
# 1        pink        552     21.124021 (high median rank)
# 2     magenta        458     14.461520 (medium median rank)
# 3         tan        300     14.281723 (high median rank)
# 4       green        762     14.061856 (lower median rank)
# 5   turquoise       1000     13.976882 (lower median rank)
# 6   darkgreen        218     11.000436 (high median rank)
# 7        blue       1000      9.858004 (not in top 20)
# 8        cyan        261      8.939884 (high median rank)
# 9         red        628      8.528091 (not in top 20)
# 10      black        618      7.999981 (not in top 20)
# 11      plum2        127      6.854030 (good rank)
# 12 lightcyan1        146      6.714021 (medium rank)
# 13     grey60        246      6.669582 (not in top 20)
# 14     violet        174      6.052549 (medium rank)
# 15      white        187      5.852156 (low medium rank)
# 16      brown       1000      5.646930 (not in top twenty)
# 17 darkviolet         70      5.146672 (higest median rank)

ROD_Dermo_stats_MR <- ROD_Dermo_mp$preservation$observed$ref.Dermo_Tol$inColumnsAlsoPresentIn.ROD_Res
ROD_Dermo_stats_MR <-ROD_Dermo_stats_MR[,c("moduleSize","medianRank.pres")]
ROD_Dermo_stats_MR_order <- ROD_Dermo_stats_MR[order(ROD_Dermo_stats_MR[,2]),]
ROD_Dermo_stats_MR_order_less20 <- ROD_Dermo_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save network
#save(ROD_Dermo_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Tol_ROD_Res_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(ROD_Dermo_stats_MR_order_less20  , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

ROR_Dermo_stats_preserved_Zsum_medianRank <- ROD_Dermo_stats_MR_order_less20[ROD_Dermo_stats_MR_order_less20$mod_name , ROD_Dermo_stats_preserved,]
# 6 modules

# Were any of these preserved modules significant in Pro_RE22 or Dermo_Tol?

Dermo_Tol_ROD_Res <- ROD_Dermo_stats_preserved[ROD_Dermo_stats_preserved$mod_name %in% Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm,]
# mod_name moduleSize Zsummary.pres
# 1        pink        552     21.124021
# 3         tan        300     14.281723
# 5   turquoise       1000     13.976882
# 11      plum2        127      6.854030
# 12 lightcyan1        146      6.714021
# 17 darkviolet         70      5.146672
Dermo_Tol_ROD_Res_ROD <- ROD_Dermo_stats_preserved[ROD_Dermo_stats_preserved$mod_name %in% ROD_Res_moduleTraitCor_Pval_df_sig_list_rm,]
#0

## ROD_Res vs. Dermo_Tol ## 
ROD_Dermo_rev_multiExpr = list(ROD_Res=list(data=ROD_Resistant_dds_rlog_matrix_common), Dermo_Tol=list(data=Dermo_Tolerant_dds_vst_matrix_common))
ROD_Res_multiColor = list(ROD_Res = ROD_Res_moduleColors)
#ROD_Dermo_rev_mp =modulePreservation(ROD_Dermo_rev_multiExpr, ROD_Res_multiColor,
#                                 referenceNetworks=1,
#                                 verbose=3,
#                                 networkType="signed hybrid", # use same signed hybrid as before
#                                 corFnc = "bicor", # use recommended bicor as before 
#                                 nPermutations=100,
#                                 randomSeed = 1, # recommended in langfelder tutorial
#                                 quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
ROD_Dermo_rev_stats = ROD_Dermo_rev_mp$preservation$Z$ref.ROD_Res$inColumnsAlsoPresentIn.Dermo_Tol
ROD_Dermo_rev_stats_order <- ROD_Dermo_rev_stats[order(-ROD_Dermo_rev_stats[,2]),c(1:2)]
ROD_Dermo_rev_stats_preserved <- ROD_Dermo_rev_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
# mod_name moduleSize Zsummary.pres
# 1    brown       1000     21.729957
# 2    plum2        131      7.874752
# 3 darkgrey        456      7.201296
# 4   orange        450      6.209411

ROD_Dermo_rev_stats_MR <- ROD_Dermo_rev_mp$preservation$observed$ref.ROD_Res$inColumnsAlsoPresentIn.Dermo_Tol
ROD_Dermo_rev_stats_MR <-ROD_Dermo_rev_stats_MR[,c("moduleSize","medianRank.pres")]
ROD_Dermo_rev_stats_MR_order <- ROD_Dermo_rev_stats_MR[order(ROD_Dermo_rev_stats_MR[,2]),]
ROD_Dermo_rev_stats_MR_order_less20 <- ROD_Dermo_rev_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save network
#save(ROD_Dermo_rev_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Tol_ROD_Res_rev_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(ROD_Dermo_stats_MR_order_less20  , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

ROD_Dermo_stats_rev_preserved_Zsum_medianRank <- ROD_Dermo_rev_stats_MR_order_less20 [ROD_Dermo_rev_stats_MR_order_less20 $mod_name %in% ROD_Dermo_rev_stats_preserved$mod_name,]

# Were any of these preserved modules significant in Pro_RE22 or Dermo_Tol?
Dermo_Tol_ROD_rev_Res <- ROD_Dermo_rev_stats_preserved[ROD_Dermo_rev_stats_preserved$mod_name %in% Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#2    plum2        131      7.874752
Dermo_Tol_ROD_rev_Res_ROD <- ROD_Dermo_rev_stats_preserved[ROD_Dermo_rev_stats_preserved$mod_name %in% ROD_Res_moduleTraitCor_Pval_df_sig_list_rm,]
# mod_name moduleSize Zsummary.pres
# 4   orange        450      6.209411

## Assess whether Dermo_Sus and ROD_Sus modules are preserved ##
ROD_Dermo_Sus_multiExpr = list(Dermo_Sus=list(data=Dermo_Susceptible_dds_vst_matrix_common), ROD_Sus=list(data=ROD_Susceptible_dds_rlog_matrix_common))
ROD_Dermo_Sus_multiColor = list(Dermo_Sus = Dermo_Sus_moduleColors)
# Load module preservation
#Dermo_Sus_ROD_Sus = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Sus_ROD_Sus_module_preservation.RData")

#ROD_Dermo_Sus_mp =modulePreservation(ROD_Dermo_Sus_multiExpr, ROD_Dermo_Sus_multiColor,
#                                     referenceNetworks=1,
#                                     verbose=3,
#                                     networkType="signed hybrid", # use same signed hybrid as before
#                                     corFnc = "bicor", # use recommended bicor as before 
#                                     nPermutations=100,
#                                     randomSeed = 1, # recommended in langfelder tutorial
#                                     quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
ROD_Dermo_Sus_stats = ROD_Dermo_Sus_mp$preservation$Z$ref.Dermo_Sus$inColumnsAlsoPresentIn.ROD_Sus 
ROD_Dermo_Sus_stats_order <- ROD_Dermo_Sus_stats[order(-ROD_Dermo_Sus_stats[,2]),c(1:2)] 
ROD_Dermo_Sus_stats_preserved <- ROD_Dermo_Sus_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
# mod_name moduleSize Zsummary.pres
# 1       black        558     14.052762
# 2   turquoise       1000     10.807381 (only one in the top 20 median rank)
# 3 lightyellow        253      8.239593
# 4      purple        446      6.069997
# 5      salmon        402      5.345025

ROD_Dermo_Sus_stats_MR <- ROD_Dermo_Sus_mp$preservation$observed$ref.Dermo_Sus$inColumnsAlsoPresentIn.ROD_Sus
ROD_Dermo_Sus_stats_MR <-ROD_Dermo_stats_MR[,c("moduleSize","medianRank.pres")]
ROD_Dermo_Sus_stats_MR_order <- ROD_Dermo_Sus_stats_MR[order(ROD_Dermo_stats_MR[,2]),]
ROD_Dermo_Sus_stats_MR_order_less20 <- ROD_Dermo_Sus_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save network
#save(ROD_Dermo_Sus_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Sus_ROD_Sus_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(ROD_Dermo_Sus_stats_MR_order_less20  , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

ROD_Dermo_Sus_stats_preserved_Zsum_medianRank <- ROD_Dermo_Sus_stats_MR_order_less20[ROD_Dermo_Sus_stats_MR_order_less20$mod_name %in% ROD_Dermo_Sus_stats_preserved,]
# 1 modules, turqoise

Dermo_Sus_ROD_Sus <-ROD_Dermo_Sus_stats_preserved[ROD_Dermo_Sus_stats_preserved$mod_name %in%Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
#0
Dermo_Sus_ROD_Sus_ROD <-ROD_Dermo_Sus_stats_preserved[ROD_Dermo_Sus_stats_preserved$mod_name %in% ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
# turqoise

## ROD_Sus vs. Dermo_Sus modules are preserved ##
ROD_Dermo_Sus_rev_multiExpr = list(ROD_Sus=list(data=ROD_Susceptible_dds_rlog_matrix_common), Dermo_Sus=list(data=Dermo_Susceptible_dds_vst_matrix_common))
ROD_Sus_multiColor = list(ROD_Sus = ROD_Sus_moduleColors)
#ROD_Dermo_Sus_rev_mp =modulePreservation(ROD_Dermo_Sus_rev_multiExpr, ROD_Sus_multiColor,
#                                     referenceNetworks=1,
#                                     verbose=3,
#                                     networkType="signed hybrid", # use same signed hybrid as before
#                                     corFnc = "bicor", # use recommended bicor as before 
#                                     nPermutations=100,
#                                     randomSeed = 1, # recommended in langfelder tutorial
#                                     quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
ROD_Dermo_Sus_rev_stats = ROD_Dermo_Sus_rev_mp$preservation$Z$ref.ROD_Sus$inColumnsAlsoPresentIn.Dermo_Sus
ROD_Dermo_Sus_rev_stats_order <- ROD_Dermo_Sus_rev_stats[order(-ROD_Dermo_Sus_rev_stats[,2]),c(1:2)] 
ROD_Dermo_Sus_rev_stats_preserved <- ROD_Dermo_Sus_rev_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
# mod_name moduleSize Zsummary.pres
# 1 turquoise       1000      6.083756
ROD_Dermo_Sus_rev_stats_MR <- ROD_Dermo_Sus_rev_mp$preservation$observed$ref.ROD_Sus$inColumnsAlsoPresentIn.Dermo_Sus
ROD_Dermo_Sus_rev_stats_MR <-ROD_Dermo_Sus_rev_stats_MR [,c("moduleSize","medianRank.pres")]
ROD_Dermo_Sus_rev_stats_MR_order <- ROD_Dermo_Sus_rev_stats_MR [order(ROD_Dermo_Sus_rev_stats_MR [,2]),]
ROD_Dermo_Sus_rev_stats_MR_order_less20 <-ROD_Dermo_Sus_rev_stats_MR_order  %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save network
#save(ROD_Dermo_Sus_rev_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Sus_ROD_Sus_rev_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(ROD_Dermo_Sus_rev_stats_MR_order_less20  , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

ROD_Dermo_Sus_rev_stats_preserved_Zsum_medianRank <- ROD_Dermo_Sus_rev_stats_MR_order_less20 [ROD_Dermo_Sus_rev_stats_MR_order_less20 $mod_name %in% ROD_Dermo_Sus_rev_stats_preserved$mod_name ,]
Dermo_Sus_ROD_Sus_rev <-ROD_Dermo_Sus_rev_stats_preserved [ROD_Dermo_Sus_rev_stats_preserved $mod_name %in%Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
#0
Dermo_Sus_ROD_Sus_ROD_rev <-ROD_Dermo_Sus_rev_stats_preserved [ROD_Dermo_Sus_rev_stats_preserved $mod_name %in% ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
# mod_name moduleSize Zsummary.pres
# 1 turquoise       1000      6.083756

## Assess whether Dermo_Tol and ROD_Sus modules are preserved ##
# NEED TO RUN OVERNIGHT
ROD_Sus_Dermo_Tol_multiExpr = list(Dermo_Tol=list(data=Dermo_Tolerant_dds_vst_matrix_common), ROD_Sus=list(data=ROD_Susceptible_dds_rlog_matrix_common))
ROD_Sus_Dermo_Tol_multiColor = list(Dermo_Tol = Dermo_Tol_moduleColors)
#ROD_Sus_Dermo_Tol_mp =modulePreservation(ROD_Sus_Dermo_Tol_multiExpr, ROD_Sus_Dermo_Tol_multiColor,
#                                     referenceNetworks=1,
#                                     verbose=3,
#                                     networkType="signed hybrid", # use same signed hybrid as before
#                                     corFnc = "bicor", # use recommended bicor as before 
#                                     nPermutations=100,
#                                     randomSeed = 1, # recommended in langfelder tutorial
#                                     quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
ROD_Sus_Dermo_Tol_stats = ROD_Sus_Dermo_Tol_mp$preservation$Z$ref.Dermo_Tol$inColumnsAlsoPresentIn.ROD_Sus 
ROD_Sus_Dermo_Tol_stats_order <-ROD_Sus_Dermo_Tol_stats[order(-ROD_Sus_Dermo_Tol_stats[,2]),c(1:2)] 
ROD_Sus_Dermo_Tol_stats_preserved <- ROD_Sus_Dermo_Tol_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
#mod_name moduleSize Zsummary.pres
#1   turquoise       1000     14.802524
#2        pink        552     14.542658
#3       green        762     14.051725
#4         red        628     13.584773
#5       black        618      9.372605
#6         tan        300      8.762566
#7     magenta        458      7.512044
#8        gold       1000      6.783474
#9      yellow        806      6.067202
#10      white        187      5.489735
#11 lightgreen        246      5.225878
#12       blue       1000      5.126448

ROD_Sus_Dermo_Tol_stats_MR <- ROD_Sus_Dermo_Tol_mp$preservation$observed$ref.Dermo_Tol$inColumnsAlsoPresentIn.ROD_Sus
ROD_Sus_Dermo_Tol_stats_MR <- ROD_Sus_Dermo_Tol_stats_MR [,c("moduleSize","medianRank.pres")]
ROD_Sus_Dermo_Tol_stats_MR_order <- ROD_Sus_Dermo_Tol_stats_MR[order(ROD_Sus_Dermo_Tol_stats_MR[,2]),]
ROD_Sus_Dermo_Tol_stats_MR_order_less20 <- ROD_Sus_Dermo_Tol_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save network
#save(ROD_Sus_Dermo_Tol_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/ROD_Sus_Dermo_Tol_mp_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(ROD_Sus_Dermo_Tol_stats_MR_order_less20  , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

ROD_Sus_Dermo_Tol_stats_preserved_Zsum_medianRank <- ROD_Sus_Dermo_Tol_stats_MR_order_less20[ROD_Sus_Dermo_Tol_stats_MR_order_less20$mod_name %in% ROD_Sus_Dermo_Tol_stats_preserved$mod_name,]

ROD_Sus_Dermo_Tol_Dermo <- ROD_Sus_Dermo_Tol_stats_preserved[ROD_Sus_Dermo_Tol_stats_preserved$mod_name %in% Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#1   turquoise       1000     14.802524
#2        pink        552     14.542658
#6         tan        300      8.762566
#9      yellow        806      6.067202
#11 lightgreen        246      5.225878

ROD_Sus_Dermo_Tol_ROD <- ROD_Sus_Dermo_Tol_stats_preserved[ROD_Sus_Dermo_Tol_stats_preserved$mod_name %in%ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
# mod_name moduleSize Zsummary.pres
# 1 turquoise       1000      14.80252

## Assess whether Dermo_Sus and ROD_Res modules are preserved ##
# NEED TO RUN OVERNIGHT
ROD_Res_Dermo_Sus_multiExpr = list(Dermo_Sus=list(data=Dermo_Susceptible_dds_vst_matrix_common), ROD_Res=list(data=ROD_Resistant_dds_rlog_matrix_common))
ROD_Res_Dermo_Sus_multiColor = list(Dermo_Sus = Dermo_Sus_moduleColors)
#ROD_Res_Dermo_Sus_mp =modulePreservation(ROD_Res_Dermo_Sus_multiExpr, ROD_Res_Dermo_Sus_multiColor,
#                                         referenceNetworks=1,
#                                         verbose=3,
#                                         networkType="signed hybrid", # use same signed hybrid as before
#                                         corFnc = "bicor", # use recommended bicor as before 
#                                         nPermutations=100,
#                                         randomSeed = 1, # recommended in langfelder tutorial
#                                         quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
ROD_Res_Dermo_Sus_stats = ROD_Res_Dermo_Sus_mp$preservation$Z$ref.Dermo_Sus$inColumnsAlsoPresentIn.ROD_Res
ROD_Res_Dermo_Sus_stats_order <-ROD_Res_Dermo_Sus_stats[order(-ROD_Res_Dermo_Sus_stats[,2]),c(1:2)] 
ROD_Res_Dermo_Sus_stats_preserved <- ROD_Res_Dermo_Sus_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
#mod_name moduleSize Zsummary.pres
#1       black        558     20.105956
#2 lightyellow        253     10.880360
#3         tan        410     10.075431
#4   turquoise       1000      9.773883
#5        blue       1000      9.663325
#6     magenta        527      7.065072
#7      yellow        848      5.736910
ROD_Res_Dermo_Sus_stats_MR <- ROD_Res_Dermo_Sus_mp$preservation$observed$ref.Dermo_Sus$inColumnsAlsoPresentIn.ROD_Res
ROD_Res_Dermo_Sus_stats_MR <- ROD_Res_Dermo_Sus_stats_MR [,c("moduleSize","medianRank.pres")]
ROD_Res_Dermo_Sus_stats_MR_order <- ROD_Res_Dermo_Sus_stats_MR [order(ROD_Res_Dermo_Sus_stats_MR [,2]),]
ROD_Res_Dermo_Sus_stats_MR_order_less20 <- ROD_Res_Dermo_Sus_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save network
#save(ROD_Res_Dermo_Sus_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/ROD_Res_Dermo_Sus_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(ROD_Res_Dermo_Sus_stats_MR_order_less20  , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

ROD_Res_Dermo_Sus_stats_preserved_Zsum_medianRank <- ROD_Res_Dermo_Sus_stats_MR_order_less20[ROD_Res_Dermo_Sus_stats_MR_order_less20$mod_name %in% ROD_Res_Dermo_Sus_stats_preserved$mod_name ,]

ROD_Res_Dermo_Sus_Dermo <- ROD_Res_Dermo_Sus_stats_preserved[ROD_Res_Dermo_Sus_stats_preserved$mod_name %in% Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
#0

ROD_Res_Dermo_Sus_ROD <- ROD_Res_Dermo_Sus_stats_preserved[ROD_Res_Dermo_Sus_stats_preserved$mod_name %in% ROD_Res_moduleTraitCor_Pval_df_sig_list_rm,]
#0

## ROD_Res vs Dermo_sus ##
# NEED TO RUN OVERNIGHT 
ROD_Res_Dermo_Sus_multiExpr = list(ROD_Res=list(data=ROD_Resistant_dds_rlog_matrix_common), Dermo_Sus=list(data=Dermo_Susceptible_dds_vst_matrix_common))
ROD_Res_Dermo_Sus_multiColor = list(ROD = Dermo_Sus_moduleColors)
ROD_Res_Dermo_Sus_mp =modulePreservation(ROD_Res_Dermo_Sus_multiExpr, ROD_Res_Dermo_Sus_multiColor,
                                        referenceNetworks=1,
                                        verbose=3,
                                        networkType="signed hybrid", # use same signed hybrid as before
                                        corFnc = "bicor", # use recommended bicor as before 
                                        nPermutations=100,
                                        randomSeed = 1, # recommended in langfelder tutorial
                                        quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
ROD_Res_Dermo_Sus_stats = ROD_Res_Dermo_Sus_mp$preservation$Z$ref.Dermo_Sus$inColumnsAlsoPresentIn.ROD_Res
ROD_Res_Dermo_Sus_stats_order <-ROD_Res_Dermo_Sus_stats[order(-ROD_Res_Dermo_Sus_stats[,2]),c(1:2)] 
ROD_Res_Dermo_Sus_stats_preserved <- ROD_Res_Dermo_Sus_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
ROD_Res_Dermo_Sus_stats_MR <- ROD_Res_Dermo_Sus_mp$preservation$observed$ref.Dermo_Sus$inColumnsAlsoPresentIn.ROD_Res
ROD_Res_Dermo_Sus_stats_MR <- ROD_Res_Dermo_Sus_stats_MR [,c("moduleSize","medianRank.pres")]
ROD_Res_Dermo_Sus_stats_MR_order <- ROD_Res_Dermo_Sus_stats_MR [order(ROD_Res_Dermo_Sus_stats_MR [,2]),]
ROD_Res_Dermo_Sus_stats_MR_order_less20 <- ROD_Res_Dermo_Sus_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save network
#save(ROD_Res_Dermo_Sus_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/ROD_Res_Dermo_Sus_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(ROD_Res_Dermo_Sus_stats_MR_order_less20  , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

ROD_Res_Dermo_Sus_stats_preserved_Zsum_medianRank <- ROD_Res_Dermo_Sus_stats_MR_order_less20[ROD_Res_Dermo_Sus_stats_MR_order_less20$mod_name %in% ROD_Res_Dermo_Sus_stats_preserved$mod_name ,]
ROD_Res_Dermo_Sus_Dermo <- ROD_Res_Dermo_Sus_stats_preserved[ROD_Res_Dermo_Sus_stats_preserved$mod_name %in% Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
ROD_Res_Dermo_Sus_ROD <- ROD_Res_Dermo_Sus_stats_preserved[ROD_Res_Dermo_Sus_stats_preserved$mod_name %in% ROD_Res_moduleTraitCor_Pval_df_sig_list_rm,]


### Assess module preservation between Probiotic and Pro_RE22 (just the probiotic network, not the full network)
Pro_Pro_RE22_multiExpr = list(Pro=list(data=Probiotic_dds_rlog_matrix_common), Pro_RE22=list(data=Pro_RE22_dds_rlog_matrix_common_Pro))
Pro_Pro_RE22_multiColor = list(Pro = Probiotic_moduleColors)
#Pro_Pro_RE22 = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Probiotic_Pro_module_preservation.RData")
#Pro_Pro_RE22_mp =modulePreservation(Pro_Pro_RE22_multiExpr, Pro_Pro_RE22_multiColor,
#                                     referenceNetworks=1,
#                                     verbose=3,
#                                     networkType="signed hybrid", # use same signed hybrid as before
#                                     corFnc = "bicor", # use recommended bicor as before 
#                                     nPermutations=100,
#                                     randomSeed = 1, # recommended in langfelder tutorial
#                                     quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Pro_Pro_RE22_stats = Pro_Pro_RE22_mp$preservation$Z$ref.Pro$inColumnsAlsoPresentIn.Pro_RE22
Pro_Pro_RE22_stats_order <- Pro_Pro_RE22_stats[order(-Pro_Pro_RE22_stats[,2]),c(1:2)] 
Pro_Pro_RE22_stats_preserved <- Pro_Pro_RE22_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
# mod_name moduleSize Zsummary.pres
# 1 turquoise       1000      7.814080
# 2     brown       1000      5.279559
Pro_Pro_RE22_stats_MR <- Pro_Pro_RE22_mp$preservation$observed$ref.Pro$inColumnsAlsoPresentIn.Pro_RE22
Pro_Pro_RE22_stats_MR <- Pro_Pro_RE22_stats_MR[,c("moduleSize","medianRank.pres")]
Pro_Pro_RE22_stats_MR_order <- Pro_Pro_RE22_stats_MR[order(Pro_Pro_RE22_stats_MR[,2]),]
Pro_Pro_RE22_stats_MR_order_less20 <- Pro_Pro_RE22_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save network
#save(Pro_Pro_RE22_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Probiotic_Pro_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(Pro_Pro_RE22_stats_MR_order_less20  , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

Pro_Pro_RE22_pro <- Pro_Pro_RE22_stats_preserved[Pro_Pro_RE22_stats_preserved$mod_name %in% Probiotic_moduleTraitCor_Pval_df_sig_list_rm,]
# 0
Pro_Pro_RE22_Pro_S4 <- Pro_Pro_RE22_stats_preserved[Pro_Pro_RE22_stats_preserved$mod_name %in% Pro_RE22_Pro_moduleTraitCor_Pval_df_sig_list_rm,]
# mod_name moduleSize Zsummary.pres
# 2    brown       1000      5.279559
Pro_Pro_RE22_Pro_RI <-Pro_Pro_RE22_stats_preserved[Pro_Pro_RE22_stats_preserved$mod_name %in% Pro_RE22_Pro_RI_moduleTraitCor_Pval_df_sig_list_rm,]
# mod_name moduleSize Zsummary.pres
# 1 turquoise       1000      7.814080
# 2     brown       1000      5.279559

## Conclusions: May be some probiotic specific 

##  Pro_RE22 (just the probiotic network, not the full network) vs Probiotic 
Pro_Pro_RE22_rev_multiExpr = list(Pro_RE22=list(data=Pro_RE22_dds_rlog_matrix_common_Pro), Pro=list(data=Probiotic_dds_rlog_matrix_common))
Pro_Pro_RE22_rev_multiColor = list(Pro_RE22 = Pro_RE22_Pro_moduleColors)
#Pro_Pro_RE22_rev_mp =modulePreservation(Pro_Pro_RE22_rev_multiExpr, Pro_Pro_RE22_rev_multiColor,
#                                     referenceNetworks=1,
#                                     verbose=3,
#                                     networkType="signed hybrid", # use same signed hybrid as before
#                                     corFnc = "bicor", # use recommended bicor as before 
#                                     nPermutations=100,
#                                     randomSeed = 1, # recommended in langfelder tutorial
#                                     quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Pro_Pro_RE22_rev_stats = Pro_Pro_RE22_rev_mp$preservation$Z$ref.Pro_RE22$inColumnsAlsoPresentIn.Pro
Pro_Pro_RE22_rev_stats_order <- Pro_Pro_RE22_rev_stats[order(-Pro_Pro_RE22_rev_stats[,2]),c(1:2)] 
Pro_Pro_RE22_rev_stats_preserved <- Pro_Pro_RE22_rev_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
#mod_name moduleSize Zsummary.pres
#1          pink        843     25.807084
#2         green       1000     18.512317
#3           red        930     17.734625
#4        yellow       1000     10.464383
#5         black        907      8.796149
#6  mediumorchid         85      7.578267
#7  midnightblue        436      7.034171
#8   lightyellow        366      7.025025
#9          cyan        507      6.833761
#10       salmon        528      5.966004
#11       orange        258      5.517729
#12  greenyellow        548      5.492732
#13       purple        561      5.456145
Pro_Pro_RE22_rev_stats_MR <- Pro_Pro_RE22_rev_mp$preservation$observed$ref.Pro_RE22$inColumnsAlsoPresentIn.Pro
Pro_Pro_RE22_rev_stats_MR <- Pro_Pro_RE22_rev_stats_MR[,c("moduleSize","medianRank.pres")]
Pro_Pro_RE22_rev_stats_MR_order <- Pro_Pro_RE22_rev_stats_MR[order(Pro_Pro_RE22_rev_stats_MR[,2]),]
Pro_Pro_RE22_rev_stats_MR_order_less20 <- Pro_Pro_RE22_rev_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save network
#save(Pro_Pro_RE22_rev_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Probiotic_Pro_rev_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(Pro_Pro_RE22_rev_stats_MR_order_less20 , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

Pro_Pro_RE22_pro_Rev <- Pro_Pro_RE22_rev_stats_preserved[Pro_Pro_RE22_rev_stats_preserved$mod_name %in% Probiotic_moduleTraitCor_Pval_df_sig_list_rm,]
#0
Pro_Pro_RE22_Pro_S4_rev <- Pro_Pro_RE22_rev_stats_preserved[Pro_Pro_RE22_rev_stats_preserved$mod_name %in% Pro_RE22_Pro_moduleTraitCor_Pval_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#1          pink        843     25.807084
#3           red        930     17.734625
#4        yellow       1000     10.464383
#5         black        907      8.796149
#6  mediumorchid         85      7.578267
#7  midnightblue        436      7.034171
#9          cyan        507      6.833761
#11       orange        258      5.517729
#13       purple        561      5.456145
Pro_Pro_RE22_Pro_RI_rev <-Pro_Pro_RE22_rev_stats_preserved[Pro_Pro_RE22_rev_stats_preserved$mod_name %in% Pro_RE22_Pro_RI_moduleTraitCor_Pval_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#4        yellow       1000     10.464383
#5         black        907      8.796149
#7  midnightblue        436      7.034171
#11       orange        258      5.517729
#13       purple        561      5.456145

# Conclusions: 

#### IDENTIFY KEY MODULES OF INTEREST ####

# Potential Probiotic Challenge specific modules not conserved in other comparisons 
# Overall - clearly extrinsic pathway is at work here
# turqoise and brown conserved between both Probiotic and Pro_RE22 experiments, but only significant for trait in the Pro_RE22 challenge
# brown module = POSITIVE IN RI, NEGATIVE CORRELATION IN S4, not sig in the S4 only challenge
# in the Probiotic Challenge = cor: 0.14472289, pval: 0.7844312482
    # brown module genes: 
    View(C_vir_all_exp_mod_sig_apop %>% filter(mod_names == "MEbrown" & exp == "Pro_RE22_Pro_RI"))
    View(C_vir_all_exp_mod_sig_apop %>% filter(mod_names == "MEbrown" & exp == "Pro_RE22_Pro_S4"))
    # indicate TLRs
    View(C_vir_apop_LFC %>% filter(experiment == "Pro_RE22" & group_by_sim == "RI_6h" | group_by_sim == "RI_24h"))
    View(C_vir_apop_LFC %>% filter(experiment == "Pro_RE22" & group_by_sim == "S4_6h" | group_by_sim == "S4_24h"))

      # PATHWAY Observations from WGCNA: Shows clear pathway all the way from outer cell receptors to transcriptional regulators 
        # Inflammatory pathway responding directly to the bacterial challenge, some molecules inhibiting apoptosis, others promoting
        # Receptor: TLR 1/2
        # Signaling molecules: CRADD (involved in PIDDosome which can trigger apoptosis)
        # transcriptional regulators: transcriptional regulator Myc-A, 
        # Effectors: Macrophage migration inhibitor factor (MIF) is an inflammatory cytokine
    # turqoise: NEGATIVE IN RI 
        # Both Inflammatory and apoptotic pathway
        # PATHWAY Observations from WGCNA: TNF AIP, TRAFs, IP3R, MAPK, IAPs, caspase 3, caspase 6, BAG molecular chaperones, NF-kB, stimulator of interferon genes
        # Interferon pathway: interferon regulatory factor 2-binding protein, stimulator of interferon genes
        # TLR pathway/NFkB activation pathway clear: IRAK4..forms complexes with MAP3K7/TAK1-TRAF6 ....this activated IKK leadering to NFkB translocation, transcription factor p65 is part of NFkB complex
          # so we have two parts of the NFkB complex here. Also TNFAIP3 which can terminate NFkB activity..helps ensure transient nature of inflammatory signaling pathway. TNFAIP promotes dissassembly
          # of IL-IR and TNFR-1 pathways and affects E3 ligases TRAF6, TRAF2, BIRC2 
        # HTRA2 binds IAPs which interact with caspases and can promote cell death or inhibit cell death. BIRC5 inhibits Caspase3, BIRC6 inhibits caspase3,7,9 specifically inhibited by HTRA2 
        # Tyrosine-protein phosphatase non-receptor type 13 - negatively regaulates FAS-induced apoptosis

# Pathogenic bacteria observations
# RE22 VS ROD Res have green module conserved, BUT NOT SIGNIFICANT WITH CHALLENGE
# RE22 VS ROD Sus have 0 modules preserved

# Parasite Challenge observations 

# Pathogenic bacteria vs. Parasite
  # RE22 VS Dermo Tol: green and black modules preserved, BUT NOT SIGNIFICANT WITH CHALLENGE
  # RE22 VS Dermo Sus: black is preserved,  BUT NOT SIGNIFICANT WITH CHALLENGE
  # ROD Sus and Dermo Sus: the following are preserved, turqoise is significant for ROD Sus none sig for Dermo Sus
      # mod_name moduleSize Zsummary.pres
      # 1       black        558     14.052762
      # 2   turquoise       1000     10.807381 (only one in the top 20 median rank)
      # 3 lightyellow        253      8.239593
      # 4      purple        446      6.069997
      # 5      salmon        402      5.345025
  # ROD Tol and Dermo Tol: None significant with trait in ROD Res, but for Dermo Tol pink, tan, turquoise, plum2, lighcyan1, and darkviolet are significant 
      # mod_name moduleSize Zsummary.pres
      # 1        pink        552     21.124021 (high median rank)
      # 2     magenta        458     14.461520 (medium median rank)
      # 3         tan        300     14.281723 (high median rank)
      # 4       green        762     14.061856 (lower median rank)
      # 5   turquoise       1000     13.976882 (lower median rank)
      # 6   darkgreen        218     11.000436 (high median rank)
      # 7        blue       1000      9.858004 (not in top 20)
      # 8        cyan        261      8.939884 (high median rank)
      # 9         red        628      8.528091 (not in top 20)
      # 10      black        618      7.999981 (not in top 20)
      # 11      plum2        127      6.854030 (good rank)
      # 12 lightcyan1        146      6.714021 (medium rank)
      # 13     grey60        246      6.669582 (not in top 20)
      # 14     violet        174      6.052549 (medium rank)
      # 15      white        187      5.852156 (low medium rank)
      # 16      brown       1000      5.646930 (not in top twenty)
      # 17 darkviolet         70      5.146672 (higest median rank)
  # ROD Sus and Dermo Tol
  # ROD Res and Dermo Sus



#### HEATMAP OF MODULES OF INTEREST ACROSS EXPERIMENTS ####


#### PERFORM CONSENSUS NETWORK ANALYSIS ACROSS C.VIR PATHOGENIC NETWORKS BACTERIA ####
## Create MultiExpr from pathogenic bacterial challenges

# We work with two sets:
nSets = 3
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("Probiotic", "ROD_Resistant","ROD_Susceptible","Pro_RE22")
shortLabels = c("Probiotic", "ROD_Resistant","ROD_Susceptible","Pro_RE22")
# Form multi-set expression data: columns starting from 9 contain actual expression data.
C_vir_Bac_multiExpr = vector(mode = "list", length = nSets)

C_vir_Bac_multiExpr[[1]] = list(data = as.data.frame(ROD_Resistant_dds_rlog_matrix_common))
names(C_vir_Bac_multiExpr[[1]]$data) = row.names(ROD_Resistant_dds_rlog_matrix_common)
rownames(C_vir_Bac_multiExpr[[1]]$data) = row.names(ROD_Resistant_dds_rlog_matrix_common)

C_vir_Bac_multiExpr[[2]] = list(data = as.data.frame(ROD_Susceptible_dds_rlog_matrix_common))
names(C_vir_Bac_multiExpr[[2]]$data) = row.names(ROD_Susceptible_dds_rlog_matrix_common)
rownames(C_vir_Bac_multiExpr[[2]]$data) = row.names(ROD_Susceptible_dds_rlog_matrix_common)

C_vir_Bac_multiExpr[[3]] = list(data = as.data.frame(Pro_RE22_dds_rlog_matrix_common_RE22))
names(C_vir_Bac_multiExpr[[3]]$data) = row.names(Pro_RE22_dds_rlog_matrix_common_RE22)
rownames(C_vir_Bac_multiExpr[[3]]$data) = row.names(Pro_RE22_dds_rlog_matrix_common_RE22)
                      
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize_Cvir_Bac = checkSets(C_vir_Bac_multiExpr)

# Check that all genes and samples have sufficiently low numbers of missing values.
#gsg = goodSamplesGenesMS(C_vir_Bac_multiExpr, verbose = 3);
gsg$allOK # [1] TRUE

# Assign clinical trait data
# Create data frame with the sample information 
colnames(ROD_Resistant_coldata_collapsed)
colnames(ROD_Susceptible_coldata_collapsed)
colnames(Pro_RE22_coldata_RE22_collapsed)
# all the same

# Recode colnames to be just the sample name and the colname
C_vir_Bac_coldata <- rbind(ROD_Resistant_coldata_collapsed,
                           ROD_Susceptible_coldata_collapsed,
                           Pro_RE22_coldata_RE22_collapsed)
levels(C_vir_Bac_coldata$Condition)
# [1] "Untreated_control"                       "Bacillus_pumilus_RI0695"                 "Control_Resistant"                       "Resistant_Challenge"                    
# [5] "Early_Susceptible"                       "Late_Susecptible"                        "Control_no_treatment"                    "Bacillus_pumilus"                       
# [9] "Phaeobacter_inhibens"                    "Vibrio_coralliilyticus_RE22_exposure_6h"

C_vir_Bac_coldata$Condition <- str_replace(C_vir_Bac_coldata$Condition, "Control_Resistant", "Control") # ROD Res
C_vir_Bac_coldata$Condition <- str_replace(C_vir_Bac_coldata$Condition, "Resistant_Challenge", "Challenge") # ROD Res
C_vir_Bac_coldata$Condition <- str_replace(C_vir_Bac_coldata$Condition, "Early_Susceptible", "Control") # ROD Sus
C_vir_Bac_coldata$Condition <- str_replace(C_vir_Bac_coldata$Condition, "Late_Susecptible", "Challenge") # ROD Sus
C_vir_Bac_coldata$Condition <- str_replace(C_vir_Bac_coldata$Condition, "Control_no_treatment", "Control") # only recode the control for Pro_RE22, we want to keep 
C_vir_Bac_coldata$Condition <- str_replace(C_vir_Bac_coldata$Condition, "Vibrio_coralliilyticus_RE22_exposure_6h", "Challenge") # only recode the control for Pro_RE22, we want to keep 
C_vir_Bac_coldata <- C_vir_Bac_coldata %>% rownames_to_column("sample")
C_vir_Bac_coldata <- C_vir_Bac_coldata[,c("sample","Condition")]

C_vir_Bac_coldata_binarize <- binarizeCategoricalColumns.pairwise(C_vir_Bac_coldata[,2])
row.names(C_vir_Bac_coldata_binarize) <- C_vir_Bac_coldata$sample

  # See how big the traits are and what are the trait and sample names
dim(C_vir_Bac_coldata)
names(C_vir_Bac_coldata)
C_vir_Bac_coldata$sample

# Form a multi-set structure that will hold the clinical traits.
C_vir_Bac_coldata_all = vector(mode="list", length = nSets)
setSamples = rownames(C_vir_Bac_multiExpr[[1]]$data)
traitRows = match(setSamples, row.names(C_vir_Bac_coldata_binarize))
C_vir_Bac_coldata_all[[1]] = list(data = as.data.frame(C_vir_Bac_coldata_binarize[traitRows,]))
rownames(C_vir_Bac_coldata_all[[1]]$data) = rownames(C_vir_Bac_multiExpr[[1]]$data)
colnames(C_vir_Bac_coldata_all[[1]]$data) = "Condition"

setSamples = rownames(C_vir_Bac_multiExpr[[2]]$data)
traitRows = match(setSamples, row.names(C_vir_Bac_coldata_binarize))
C_vir_Bac_coldata_all[[2]] = list(data = as.data.frame(C_vir_Bac_coldata_binarize[traitRows,]))
rownames(C_vir_Bac_coldata_all[[2]]$data) = rownames(C_vir_Bac_multiExpr[[2]]$data)
colnames(C_vir_Bac_coldata_all[[2]]$data) = "Condition"

setSamples = rownames(C_vir_Bac_multiExpr[[3]]$data)
traitRows = match(setSamples, row.names(C_vir_Bac_coldata_binarize))
C_vir_Bac_coldata_all[[3]] = list(data = as.data.frame(C_vir_Bac_coldata_binarize[traitRows,]))
rownames(C_vir_Bac_coldata_all[[3]]$data) = rownames(C_vir_Bac_multiExpr[[3]]$data)
colnames(C_vir_Bac_coldata_all[[3]]$data) = "Condition"

# Define data set dimensions
nGenes = exprSize_Cvir_Bac$nGenes
nSamples = exprSize_Cvir_Bac$nSamples

# Choose soft-thresholing powers
# Choose set based on where all three datasets level off - picking 9 since none fit the topology free network criterion

# Network construction 
#C_vir_Bac_net = blockwiseConsensusModules(
#  C_vir_Bac_multiExpr,
#  power = 9,
#  minModuleSize = 30, 
#  deepSplit = 2,
#  pamRespectsDendro = FALSE,
#  mergeCutHeight = 0.25, 
#  numericLabels = TRUE,
#  minKMEtoStay = 0,
#  verbose = 5,
#  TOMType = "signed", # use signed TOM type
#  networkType= "signed hybrid", # use signed hybrid network type
#  corType = "bicor", # use suggested bicor
#  maxBlockSize = 20000)
names(C_vir_Bac_net)
C_vir_Bac_consMEs = C_vir_Bac_net$multiMEs
C_vir_Bac_moduleLabels = C_vir_Bac_net$colors
# Convert the numeric labels to color labels
C_vir_Bac_moduleColors = labels2colors(C_vir_Bac_moduleLabels)
C_vir_Bac_consTree = C_vir_Bac_net$dendrograms[[1]]

save(C_vir_Bac_net, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_vir_Bac_consensus_modules_net.RData")


# plot dendrogram
sizeGrWindow(8,6)
plotDendroAndColors(C_vir_Bac_consTree, C_vir_Bac_moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()

# Relate consensus modules to external traits
# Set up variables to contain the module-trait correlations
C_vir_Bac_moduleTraitCor = list()
C_vir_Bac_moduleTraitPvalue = list()
# Calculate the correlations
C_vir_Bac_moduleTraitCor[[1]] = cor(C_vir_Bac_consMEs[[1]]$data, C_vir_Bac_coldata_all[[1]]$data, use = "p")
C_vir_Bac_moduleTraitPvalue[[1]] = corPvalueFisher(C_vir_Bac_moduleTraitCor[[1]], exprSize_Cvir_Bac$nSamples[1])
C_vir_Bac_moduleTraitCor[[2]] = cor(C_vir_Bac_consMEs[[2]]$data, C_vir_Bac_coldata_all[[2]]$data, use = "p")
C_vir_Bac_moduleTraitPvalue[[2]] = corPvalueFisher(C_vir_Bac_moduleTraitCor[[2]], exprSize_Cvir_Bac$nSamples[2])
C_vir_Bac_moduleTraitCor[[3]] = cor(C_vir_Bac_consMEs[[3]]$data, C_vir_Bac_coldata_all[[3]]$data, use = "p")
C_vir_Bac_moduleTraitPvalue[[3]] = corPvalueFisher(C_vir_Bac_moduleTraitCor[[3]], exprSize_Cvir_Bac$nSamples[3])


## DISPLAY MODULE TRAIT CORRELATIONS FOR CONSENSUS MODULES
# Initialize matrices to hold the consensus correlation and p-value
C_vir_Bac_consensusCor = matrix(NA, nrow(C_vir_Bac_moduleTraitCor[[1]]), ncol(C_vir_Bac_moduleTraitCor[[1]]))
C_vir_Bac_consensusPvalue = matrix(NA, nrow(C_vir_Bac_moduleTraitCor[[1]]), ncol(C_vir_Bac_moduleTraitCor[[1]]))
# Find consensus negative correlations
C_vir_Bac_negative = C_vir_Bac_moduleTraitCor[[1]] < 0 & C_vir_Bac_moduleTraitCor[[2]] < 0
C_vir_Bac_consensusCor[C_vir_Bac_negative] = pmax(C_vir_Bac_moduleTraitCor[[1]][C_vir_Bac_negative], C_vir_Bac_moduleTraitCor[[2]][C_vir_Bac_negative])
C_vir_Bac_consensusPvalue[C_vir_Bac_negative] = pmax(C_vir_Bac_moduleTraitPvalue[[1]][C_vir_Bac_negative], C_vir_Bac_moduleTraitPvalue[[2]][C_vir_Bac_negative])
# Find consensus positive correlations
C_vir_Bac_positive = C_vir_Bac_moduleTraitCor[[1]] > 0 & C_vir_Bac_moduleTraitCor[[2]] > 0
C_vir_Bac_consensusCor[C_vir_Bac_positive] = pmin(C_vir_Bac_moduleTraitCor[[1]][C_vir_Bac_positive], C_vir_Bac_moduleTraitCor[[2]][C_vir_Bac_positive])
C_vir_Bac_consensusPvalue[C_vir_Bac_positive] = pmax(C_vir_Bac_moduleTraitPvalue[[1]][C_vir_Bac_positive], C_vir_Bac_moduleTraitPvalue[[2]][C_vir_Bac_positive])

# Convert numerical lables to colors for labeling of modules in the plot
C_vir_Bac_MEColors = labels2colors(as.numeric(substring(names(C_vir_Bac_consMEs[[1]]$data), 3)));
C_vir_Bac_MEColorNames = paste("ME", C_vir_Bac_MEColors, sep="");

# Plot the consensus modules for a set, but is same across sets
C_vir_Bac_textMatrix = paste(signif(C_vir_Bac_consensusCor, 2), "\n(",
                   signif(C_vir_Bac_consensusPvalue, 1), ")", sep = "")
dim(C_vir_Bac_textMatrix) = dim(C_vir_Bac_moduleTraitCor[[1]])
sizeGrWindow(10,7)
#pdf(file = "Plots/ModuleTraitRelationships-consensus.pdf", wi = 10, he = 7);
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = C_vir_Bac_consensusCor,
               xLabels = names(C_vir_Bac_coldata_all[[1]]$data),
               yLabels = C_vir_Bac_MEColorNames,
               ySymbols =C_vir_Bac_MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = C_vir_Bac_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Consensus module--trait relationships across\n",
                            paste(setLabels, collapse = " and ")))

# Conclusion: no significant consensus modules across all bacterial challenges

#### PERFORM CONSENSUS NETWORK ANALYSIS ACROSS C.VIR PROBIOTIC ####
## Create MultiExpr from pathogenic bacterial challenges
## Running consensus analysis on this doesn't make too much sense considering it's only two experiments
# Code has errors, check before re-running

# We work with two sets:
#nSets = 2
## For easier labeling of plots, create a vector holding descriptive names of the two sets.
#setLabels = c("Probiotic", "Pro_RE22")
#shortLabels = c("Probiotic", "Pro_RE22")
## Form multi-set expression data: columns starting from 9 contain actual expression data.
#C_vir_Pro_multiExpr = vector(mode = "list", length = nSets)
#
#C_vir_Pro_multiExpr[[1]] = list(data = as.data.frame(Probiotic_dds_rlog_matrix_common))
#names(C_vir_Pro_multiExpr[[1]]$data) = row.names(Probiotic_dds_rlog_matrix_common)
#rownames(C_vir_Pro_multiExpr[[1]]$data) = row.names(Probiotic_dds_rlog_matrix_common)
#
#C_vir_Pro_multiExpr[[2]] = list(data = as.data.frame(Pro_RE22_dds_rlog_matrix_common_Pro))
#names(C_vir_Pro_multiExpr[[2]]$data) = row.names(Pro_RE22_dds_rlog_matrix_common_Pro)
#rownames(C_vir_Pro_multiExpr[[2]]$data) = row.names(Pro_RE22_dds_rlog_matrix_common_Pro)
#
## Check that the data has the correct format for many functions operating on multiple sets:
#exprSize_Cvir_Pro = checkSets(C_vir_Pro_multiExpr)
#
## Check that all genes and samples have sufficiently low numbers of missing values.
##gsg = goodSamplesGenesMS(C_vir_Pro_multiExpr, verbose = 3)
#gsg$allOK # [1] TRUE
#
## Assign clinical trait data
## Create data frame with the sample information 
#colnames(Probiotic_coldata_collapsed)
#colnames(Pro_RE22_coldata_Pro_collapsed)
## all the same
#
## Recode colnames to be just the sample name and the colname
#C_vir_Pro_coldata <- rbind(Probiotic_coldata_collapsed,
#                           Pro_RE22_coldata_Pro_collapsed)
#levels(C_vir_Pro_coldata$Condition)
##[1] "Untreated_control"    "Bacillus_pumilus_RI0695" "Control_no_treatment"    "Bacillus_pumilus"        "Phaeobacter_inhibens"   
#
#C_vir_Pro_coldata$Condition <- str_replace(C_vir_Pro_coldata$Condition, "Untreated_control", "Control") # ROD Res
#C_vir_Pro_coldata$Condition <- str_replace(C_vir_Pro_coldata$Condition, "Bacillus_pumilus_RI0695", "Bacillus_pumilus") # ROD Res
#C_vir_Pro_coldata$Condition <- str_replace(C_vir_Pro_coldata$Condition, "Control_no_treatment", "Control") # only recode the control for Pro_RE22, we want to keep 
#C_vir_Pro_coldata <- C_vir_Pro_coldata %>% rownames_to_column("sample")
#C_vir_Pro_coldata <- C_vir_Pro_coldata[,c("sample","Condition")]
#
#C_vir_Pro_coldata_binarize <- binarizeCategoricalColumns.pairwise(C_vir_Pro_coldata[,2])
#row.names(C_vir_Pro_coldata_binarize) <- C_vir_Pro_coldata$sample
#
## See how big the traits are and what are the trait and sample names
#dim(C_vir_Pro_coldata)
#names(C_vir_Pro_coldata)
#C_vir_Pro_coldata$sample
#
## Form a multi-set structure that will hold the clinical traits.
#C_vir_Pro_coldata_all = vector(mode="list", length = nSets)
#setSamples = rownames(C_vir_Pro_multiExpr[[1]]$data)
#traitRows = match(setSamples, row.names(C_vir_Pro_coldata_binarize))
#C_vir_Pro_coldata_all[[1]] = list(data = as.data.frame(C_vir_Pro_coldata_binarize[traitRows,]))
#rownames(C_vir_Pro_coldata_all[[1]]$data) = rownames(C_vir_Pro_multiExpr[[1]]$data)
#colnames(C_vir_Pro_coldata_all[[1]]$data) = c("data.Control.vs.Bacillus_pumilus", "data.Phaeobacter_inhibens.vs.Bacillus_pumilus", "data.Phaeobacter_inhibens.vs.Control")
#
#setSamples = rownames(C_vir_Pro_multiExpr[[2]]$data)
#traitRows = match(setSamples, row.names(C_vir_Pro_coldata_binarize))
#C_vir_Pro_coldata_all[[2]] = list(data = as.data.frame(C_vir_Pro_coldata_binarize[traitRows,]))
#
## Define data set dimensions
#nGenes = exprSize_Cvir_Pro$nGenes
#nSamples = exprSize_Cvir_Pro$nSamples
#
## Choose soft-thresholing powers
## Choose set based on where all three datasets level off - picking 9 since none fit the topology free network criterion
#
## Network construction 
##C_vir_Pro_net = blockwiseConsensusModules(
##  C_vir_Pro_multiExpr,
##  power = 9,
##  minModuleSize = 30, 
##  deepSplit = 2,
##  pamRespectsDendro = FALSE,
##  mergeCutHeight = 0.25, 
##  numericLabels = TRUE,
##  minKMEtoStay = 0,
##  verbose = 5,
##  TOMType = "signed", # use signed TOM type
##  networkType= "signed hybrid", # use signed hybrid network type
##  corType = "bicor", # use suggested bicor
##  maxBlockSize = 20000)
#names(C_vir_Pro_net)
#C_vir_Pro_consMEs = C_vir_Pro_net$multiMEs
#C_vir_Pro_moduleLabels = C_vir_Pro_net$colors
## Convert the numeric labels to color labels
#C_vir_Pro_moduleColors = labels2colors(C_vir_Pro_moduleLabels)
#C_vir_Pro_consTree = C_vir_Pro_net$dendrograms[[1]]
#
#save(C_vir_Pro_net, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_vir_Pro_consensus_modules_net.RData")
#
## plot dendrogram
#sizeGrWindow(8,6)
#plotDendroAndColors(C_vir_Pro_consTree, C_vir_Pro_moduleColors,
#                    "Module colors",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05,
#                    main = "Consensus gene dendrogram and module colors")
#dev.off()
#
## Relate consensus modules to external traits
## Set up variables to contain the module-trait correlations
#C_vir_Pro_moduleTraitCor = list()
#C_vir_Pro_moduleTraitPvalue = list()
## Calculate the correlations
#C_vir_Pro_moduleTraitCor[[1]] = cor(C_vir_Pro_consMEs[[1]]$data, C_vir_Pro_coldata_all[[1]]$data, use = "p")
#C_vir_Pro_moduleTraitPvalue[[1]] = corPvalueFisher(C_vir_Pro_moduleTraitCor[[1]], exprSize_Cvir_Pro$nSamples[1])
#C_vir_Pro_moduleTraitCor[[2]] = cor(C_vir_Pro_consMEs[[2]]$data, C_vir_Pro_coldata_all[[2]]$data, use = "p")
#C_vir_Pro_moduleTraitPvalue[[2]] = corPvalueFisher(C_vir_Pro_moduleTraitCor[[2]], exprSize_Cvir_Pro$nSamples[2])
#
### DISPLAY MODULE TRAIT CORRELATIONS FOR CONSENSUS MODULES, I care about columns 1 and 3 
## Need to fix and invesigate this!
## Initialize matrices to hold the consensus correlation and p-value
#C_vir_Pro_consensusCor = matrix(NA, nrow(C_vir_Pro_moduleTraitCor[[1]]), ncol(C_vir_Pro_moduleTraitCor[[1]]))
#C_vir_Pro_consensusPvalue = matrix(NA, nrow(C_vir_Pro_moduleTraitCor[[1]]), ncol(C_vir_Pro_moduleTraitCor[[1]]))
## Find consensus negative correlations
#C_vir_Pro_negative = C_vir_Pro_moduleTraitCor[[1]] < 0 & C_vir_Pro_moduleTraitCor[[2]] < 0
#C_vir_Pro_consensusCor[C_vir_Pro_negative] = pmax(C_vir_Pro_moduleTraitCor[[1]][C_vir_Pro_negative], C_vir_Pro_moduleTraitCor[[2]][C_vir_Pro_negative])
#C_vir_Pro_consensusPvalue[C_vir_Pro_negative] = pmax(C_vir_Pro_moduleTraitPvalue[[1]][C_vir_Pro_negative], C_vir_Pro_moduleTraitPvalue[[2]][C_vir_Pro_negative])
## Find consensus positive correlations
#C_vir_Pro_positive = C_vir_Pro_moduleTraitCor[[1]] > 0 & C_vir_Pro_moduleTraitCor[[2]] > 0
#C_vir_Pro_consensusCor[C_vir_Pro_positive] = pmin(C_vir_Pro_moduleTraitCor[[1]][C_vir_Pro_positive], C_vir_Pro_moduleTraitCor[[2]][C_vir_Pro_positive])
#C_vir_Pro_consensusPvalue[C_vir_Pro_positive] = pmax(C_vir_Pro_moduleTraitPvalue[[1]][C_vir_Pro_positive], C_vir_Pro_moduleTraitPvalue[[2]][C_vir_Pro_positive])
#
## Convert numerical lables to colors for labeling of modules in the plot
#C_vir_Pro_MEColors = labels2colors(as.numeric(substring(names(C_vir_Pro_consMEs[[1]]$data), 3)));
#C_vir_Pro_MEColorNames = paste("ME", C_vir_Pro_MEColors, sep="");
#
## Plot the consensus modules for a set, but is same across sets
#C_vir_Pro_textMatrix = paste(signif(C_vir_Pro_consensusCor, 2), "\n(",
#                             signif(C_vir_Pro_consensusPvalue, 1), ")", sep = "")
#dim(C_vir_Pro_textMatrix) = dim(C_vir_Pro_moduleTraitCor[[1]])
#sizeGrWindow(10,7)
##pdf(file = "Plots/ModuleTraitRelationships-consensus.pdf", wi = 10, he = 7);
#par(mar = c(6, 8.8, 3, 2.2));
#labeledHeatmap(Matrix = C_vir_Pro_consensusCor,
#               xLabels = names(C_vir_Pro_coldata_all[[1]]$data),
#               yLabels = C_vir_Pro_MEColorNames,
#               ySymbols =C_vir_Pro_MEColorNames,
#               colorLabels = FALSE,
#               colors = greenWhiteRed(50),
#               textMatrix = C_vir_Pro_textMatrix,
#               setStdMargins = FALSE,
#               cex.text = 0.5,
#               zlim = c(-1,1),
#               main = paste("Consensus module--trait relationships across\n",
#                            paste(setLabels, collapse = " and ")))
#

#### RUNNING DERMO TOLERANT FULL DATASET ####

# Pick soft threshold
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# save for running TOM from cluster
write.table(Dermo_Tolerant_dds_vst_matrix, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Tolerant_dds_vst_matrix.table")

# Call the network topology analysis function, following general recommendations to set network type to "signed hybrid" and using the "bicor" correlation
#Dermo_Tolerant_dds_vst_matrix_sft <- pickSoftThreshold(Dermo_Tolerant_dds_vst_matrix, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
save(Dermo_Tolerant_dds_vst_matrix_sft, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Tolerant_dds_vst_matrix_sft")
#Dermo
# Scale-free topology fit index as a function of the soft-thresholding power
plot(Dermo_Tolerant_dds_vst_matrix_sft$fitIndices[,1], -sign(Dermo_Tolerant_dds_vst_matrix_sft$fitIndices[,3])*Dermo_Tolerant_dds_vst_matrix_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(Dermo_Tolerant_dds_vst_matrix_sft$fitIndices[,1], -sign(Dermo_Tolerant_dds_vst_matrix_sft$fitIndices[,3])*Dermo_Tolerant_dds_vst_matrix_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(Dermo_Tolerant_dds_vst_matrix_sft$fitIndices[,1], Dermo_Tolerant_dds_vst_matrix_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(Dermo_Tolerant_dds_vst_matrix_sft$fitIndices[,1], Dermo_Tolerant_dds_vst_matrix_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Softthreshold of 3 since this is lowest value past 0.9 we start to see flattening 

## ONE STEP NETWORK CONSTRUCTION, MODULE DETECTION, MODULE DENDROGRAM INSPECTION ##
Dermo_Tol_full_net <- load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Tol_full_network.RData")
# Dermo_Tol_full_net = blockwiseModules(Dermo_Tolerant_dds_vst_matrix, power = 3, # picked suitable power in the code above 
#                                 TOMType = "signed", # use signed TOM type
#                                 networkType= "signed hybrid", # use signed hybrid network type
#                                 corType = "bicor", # use suggested bicor
#                                 TminModuleSize = 30, # recommended default
#                                 reassignThreshold = 0, # recommended default
#                                 mergeCutHeight = 0.25, # recommended default
#                                 numericLabels = TRUE, # recommended default
#                                 pamRespectsDendro = FALSE,# recommended default
#                                 verbose = 3, 
#                                 maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
# 
# How many modules identified
table(Dermo_Tol_full_net$colors)
# Plot dendrogram with colors
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
#Dermo_Tol_full_mergedColors = labels2colors(Dermo_Tol_full_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(Dermo_Tol_full_net$dendrograms[[1]], Dermo_Tol_full_mergedColors[Dermo_Tol_full_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#Dermo_Tol_full_moduleLabels = Dermo_Tol_full_net$colors
#Dermo_Tol_full_moduleColors = labels2colors(Dermo_Tol_full_net$colors)
#Dermo_Tol_full_MEs = Dermo_Tol_full_net$MEs
#Dermo_Tol_full_geneTree = Dermo_Tol_full_net$dendrograms[[1]]
# save network
#save(Dermo_Tol_full_net, Dermo_Tol_full_mergedColors, Dermo_Tol_full_moduleLabels, Dermo_Tol_full_moduleColors, Dermo_Tol_full_MEs, Dermo_Tol_full_geneTree, 
#     file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Tol_full_network.RData")

# Quantify module trait associations 
# Define numbers of genes and samples
Dermo_Tol_full_nGenes = ncol(Dermo_Tolerant_dds_vst_matrix)
Dermo_Tol_full_nSamples = nrow(Dermo_Tolerant_dds_vst_matrix)

# Recalculate MEs with color labels
Dermo_Tol_full_MEs0 = moduleEigengenes(Dermo_Tolerant_dds_vst_matrix, Dermo_Tol_full_moduleColors)$eigengenes
Dermo_Tol_full_MEs = orderMEs(Dermo_Tol_full_MEs0)
Dermo_Tol_full_moduleTraitCor = cor(Dermo_Tol_full_MEs, Dermo_Tolerant_coldata_collapsed_binarize, use = "p");
Dermo_Tol_full_moduleTraitPvalue = corPvalueStudent(Dermo_Tol_full_moduleTraitCor, Dermo_Tol_full_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trai
sizeGrWindow(10,6)
# Will display correlations and their p-values
Dermo_Tol_full_textMatrix = paste(signif(Dermo_Tol_full_moduleTraitCor, 2), "\n(",
                             signif(Dermo_Tol_full_moduleTraitPvalue, 1), ")", sep = "");
dim(Dermo_Tol_full_textMatrix) = dim(Dermo_Tol_full_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = Dermo_Tol_full_moduleTraitCor,
               xLabels = names(Dermo_Tolerant_coldata_collapsed_binarize),
               yLabels = names(Dermo_Tol_full_MEs),
               ySymbols = names(Dermo_Tol_full_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = Dermo_Tol_full_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
Dermo_Tol_full_moduleTraitCor_df <- as.data.frame(Dermo_Tol_full_moduleTraitCor)
Dermo_Tol_full_moduleTraitCor_df$mod_names <- row.names(Dermo_Tol_full_moduleTraitCor_df)
Dermo_Tol_full_moduleTraitCor_df <- Dermo_Tol_full_moduleTraitCor_df[,c("mod_names","Condition.Injected.vs.Control")]
Dermo_Tol_full_moduleTraitPvalue_df <- as.data.frame(Dermo_Tol_full_moduleTraitPvalue)
Dermo_Tol_full_moduleTraitPvalue_df$mod_names <- row.names(Dermo_Tol_full_moduleTraitPvalue_df)
Dermo_Tol_full_moduleTraitPvalue_df <- Dermo_Tol_full_moduleTraitPvalue_df[,c("mod_names","Condition.Injected.vs.Control")]
colnames(Dermo_Tol_full_moduleTraitPvalue_df)[2] <- "pvalue"

Dermo_Tol_full_moduleTraitCor_Pval_df <- join(Dermo_Tol_full_moduleTraitCor_df, Dermo_Tol_full_moduleTraitPvalue_df, by = "mod_names")

# Significantly correlated modules
Dermo_Tol_full_moduleTraitCor_Pval_df[order(Dermo_Tol_full_moduleTraitCor_Pval_df$pvalue),]
class(Dermo_Tol_full_moduleTraitCor_Pval_df$pvalue) # numeric
Dermo_Tol_full_moduleTraitCor_Pval_df_sig <- Dermo_Tol_full_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.055)
Dermo_Tol_full_moduleTraitCor_Pval_df_sig #7 

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
Dermo_Tol_full_moduleTraitCor_Pval_df_sig_list <- Dermo_Tol_full_moduleTraitCor_Pval_df_sig$mod_names
Dermo_Tol_full_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Dermo_Tol_full_moduleTraitCor_Pval_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
matrix_common= Dermo_Tolerant_dds_vst_matrix_common
moduleColors=Dermo_Tol_full_moduleColors
lookup =   C_vir_rtracklayer_apop_product_final

lookup_mod_apop <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$ID %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
names(Dermo_Tol_full_moduleTraitCor_Pval_df_sig_list_rm) <- c("darkslateblue", "turquoise",     "greenyellow",   "skyblue3" ,     "cyan"  ,        "red"  ,         "tan" )

Dermo_Tol_full_module_apop <- lapply(Dermo_Tol_full_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop)
Dermo_Tol_full_module_apop_df <- do.call(rbind,Dermo_Tol_full_module_apop)
Dermo_Tol_full_module_apop_df$mod_names <- gsub("\\..*","",row.names(Dermo_Tol_full_module_apop_df))
Dermo_Tol_full_module_apop_df$mod_names <- gsub("^","ME",Dermo_Tol_full_module_apop_df$mod_names)
# add module significance
Dermo_Tol_full_module_apop_df <- left_join(Dermo_Tol_full_module_apop_df,Dermo_Tol_full_moduleTraitCor_Pval_df_sig)
Dermo_Tol_full_module_apop_df$exp <- "Dermo_Tol"

## Gene relationship to trait and important modules: Gene Significance and Module Membership
# Define variable injected 
Dermo_Tol_full_injection = as.data.frame(Dermo_Tolerant_coldata_collapsed_binarize$Condition.Injected.vs.Control);
names(Dermo_Tol_full_injection) = "injection"
# names (colors) of the modules
Dermo_Tol_full_modNames = substring(names(Dermo_Tol_full_MEs), 3)
Dermo_Tol_full_geneModuleMembership = as.data.frame(cor(Dermo_Tolerant_dds_vst_matrix, Dermo_Tol_full_MEs, use = "p"))
Dermo_Tol_full_MMPvalue = as.data.frame(corPvalueStudent(as.matrix(Dermo_Tol_full_geneModuleMembership), Dermo_Tol_full_nSamples))

names(Dermo_Tol_full_geneModuleMembership) = paste("MM", Dermo_Tol_full_modNames, sep="")
names(Dermo_Tol_full_MMPvalue) = paste("p.MM", Dermo_Tol_full_modNames, sep="")

Dermo_Tol_full_geneTraitSignificance = as.data.frame(cor(Dermo_Tolerant_dds_vst_matrix,Dermo_Tol_full_injection, use = "p"))
Dermo_Tol_full_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(Dermo_Tol_full_geneTraitSignificance), Dermo_Tol_full_nSamples))

names(Dermo_Tol_full_geneTraitSignificance) = paste("GS.", names(Dermo_Tol_full_injection), sep="")
names(Dermo_Tol_full_GSPvalue) = paste("p.GS.", names(Dermo_Tol_full_injection), sep="")

## Intramodular analysis: identifying genes with high GS and MM
# Using the GS and MM measures, we can identify genes that have a high significance for disease challenge
# as well as high module membership in interesting modules. As an example, we look at the brown module 
# that has the highest association with weight. We plot a scatterplot of Gene Significance vs. Module Membership in the brown module:
#Dermo_Tol_full_module = "skyblue3" # not that strong of an association
Dermo_Tol_full_module = "turquoise" # cor = 0.37
#Dermo_Tol_full_module = "greenyellow" # cor = 0.33
Dermo_Tol_full_column = match(Dermo_Tol_full_module, Dermo_Tol_full_modNames)
Dermo_Tol_full_moduleGenes = Dermo_Tol_full_moduleColors==Dermo_Tol_full_module
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(Dermo_Tol_full_geneModuleMembership[Dermo_Tol_full_moduleGenes, Dermo_Tol_full_column]),
                   abs(Dermo_Tol_full_geneTraitSignificance[Dermo_Tol_full_moduleGenes, 1]),
                   xlab = paste("Module Membership in", Dermo_Tol_full_module, "module"),
                   ylab = "Gene significance for challenge",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = Dermo_Tol_full_module)

## IDENTIFY HUB GENES IN EACH SIG MODULE ##
Dermo_Tol_full_colorh = c("darkslateblue", "turquoise",     "greenyellow",   "skyblue3" ,     "cyan"  ,        "red"  ,         "tan" )

Dermo_Tol_full_Module_hub_genes <- chooseTopHubInEachModule(
  Dermo_Tolerant_dds_vst_matrix_common, 
  Dermo_Tol_full_colorh, 
  power = 3,  # power used for the adjacency network
  type = "signed hybrid", 
  corFnc = "bicor"
  )
class(Dermo_Tol_full_Module_hub_genes)
Dermo_Tol_full_Module_hub_genes_df <- as.data.frame(Dermo_Tol_full_Module_hub_genes)
colnames(Dermo_Tol_full_Module_hub_genes_df)[1] <- "ID"
Dermo_Tol_full_Module_hub_genes_apop <- merge(Dermo_Tol_full_Module_hub_genes_df, C_vir_rtracklayer, by = "ID")
nrow(Dermo_Tol_full_Module_hub_genes_apop) # 7, none involed in apoptosis

## Compare turquoise module from consensus network to the full network...doesn't make sense to do this since I haven't confirmed that they are preserved yet 

Dermo_Tol_full_module_apop_df_turq <- Dermo_Tol_full_module_apop_df %>% filter(mod_names == "MEturquoise")
Dermo_Tol_full_module_apop_df_turq$type <- "full"
Dermo_Tol_module_apop_df_turq <- Dermo_Tol_module_apop_df %>% filter(mod_names == "MEturquoise")
Dermo_Tol_module_apop_df_turq$type <- "consensus"

Dermo_Tol_turq_comparison <- full_join(Dermo_Tol_full_module_apop_df_turq[,c("product","transcript_id","type")], Dermo_Tol_module_apop_df_turq[,c("product","transcript_id","type")], by ="transcript_id")
# few shared genes 

## Export modules to cytoscape for visualization ###
# Recalculate topological overlap if needed
#Dermo_Tol_full_TOM = TOMsimilarityFromExpr(Dermo_Tolerant_dds_vst_matrix,
#                                           power = 3, # picked suitable power in the code above 
#                                           TOMType = "signed", # use signed TOM type
#                                           networkType= "signed hybrid", # use signed hybrid network type
#                                           corType = "bicor") # use suggested bicor
#
#save(Dermo_Tol_full_TOM, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Tol_full_TOM.RData" )

# Select modules
Dermo_Tol_full_modules = c("darkslateblue", "turquoise", "greenyellow",   "skyblue3" , "cyan"  ,"red"  , "tan" )
# Select module probes
Dermo_Tol_full_probes = colnames(Dermo_Tolerant_dds_vst_matrix)
# export moduleColors file for use in cluster
write.table(Dermo_Tol_full_moduleColors, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Tol_full_moduleColors.table")

Dermo_Tol_full_inModule = is.finite(match(Dermo_Tol_full_moduleColors, Dermo_Tol_full_modules))
Dermo_Tol_full_modProbes = Dermo_Tol_full_probes[Dermo_Tol_full_inModule]
Dermo_Tol_full_modGenes = C_vir_rtracklayer$ID[match(Dermo_Tol_full_modProbes, C_vir_rtracklayer$ID)]
# Select the corresponding Topological Overlap
Dermo_Tol_full_modTOM = Dermo_Tol_full_TOM[Dermo_Tol_full_inModule, Dermo_Tol_full_inModule]

dimnames(Dermo_Tol_full_modTOM) = list(Dermo_Tol_full_modProbes, Dermo_Tol_full_modProbes)
# Export the network into edge and node list files Cytoscape can read
#Dermo_Tol_full_cyt = exportNetworkToCytoscape(Dermo_Tol_full_modTOM,
#                               edgeFile = paste("CytoscapeInput-edges-", paste(Dermo_Tol_full_modules, collapse="-"), ".txt", sep=""),
#                               nodeFile = paste("CytoscapeInput-nodes-", paste(Dermo_Tol_full_modules, collapse="-"), ".txt", sep=""),
#                               weighted = TRUE,
#                               threshold = 0.02,
#                               nodeNames = Dermo_Tol_full_modProbes,
#                               altNodeNames = Dermo_Tol_full_modGenes,
#                               nodeAttr = Dermo_Tol_full_moduleColors[Dermo_Tol_full_inModule])
#

## Upload finished cytoscape network node file so that annotation information can be added
# The below code may not be necessary because I can add the apoptosis data table as a separate data table cytoscape will automatically attach 
  # it to the network 
#Dermo_Tol_full_cyt_node <- read.table(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/CytoscapeInput-nodes-Dermo_Tull_fulldarkslateblue-turquoise-greenyellow-skyblue3-cyan-red-tan.txt",
#                                      sep="\t", skip=1) # skip the first line because it has the original column names
#
## original col names were : nodeName  altName nodeAttr[nodesPresent, ]
#
#colnames(Dermo_Tol_full_cyt_node)[c(1:3)] <- c("ID","altName","nodesPresent")
#Dermo_Tol_full_cyt_node <- left_join(Dermo_Tol_full_cyt_node, C_vir_rtracklayer_apop_product_final[,c("product","gene","ID","transcript_id")], by ="ID")
#
#write.table(Dermo_Tol_full_cyt_node, sep = " ", quote= FALSE, row.names=FALSE, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/CytoscapeInput-nodes-Dermo_Tull_fulldarkslateblue-turquoise-greenyellow-skyblue3-cyan-red-tan_ANNOT.txt")
## Compare IAPs and GIMAPs between Consensus set and Full set

# Consensus set in significant modules
Dermo_Tol_module_apop_df_IAP <- Dermo_Tol_module_apop_df[grepl("IAP", Dermo_Tol_module_apop_df$product, ignore.case = TRUE),]
Dermo_Tol_module_apop_df_IAP <- Dermo_Tol_module_apop_df_IAP[order(Dermo_Tol_module_apop_df_IAP$product),]
Dermo_Tol_module_apop_df_GIMAP <- Dermo_Tol_module_apop_df[grepl("IMAP", Dermo_Tol_module_apop_df$product, ignore.case = TRUE),]
Dermo_Tol_module_apop_df_GIMAP <- Dermo_Tol_module_apop_df_GIMAP[order(Dermo_Tol_module_apop_df_GIMAP$product),]

# full set in significant modules
Dermo_Tol_full_module_apop_df_IAP <-   Dermo_Tol_full_module_apop_df[grepl("IAP", Dermo_Tol_full_module_apop_df$product, ignore.case = TRUE),]
Dermo_Tol_full_module_apop_df_IAP <-   Dermo_Tol_full_module_apop_df_IAP[order(Dermo_Tol_full_module_apop_df_IAP$product),]
Dermo_Tol_full_module_apop_df_GIMAP <- Dermo_Tol_full_module_apop_df[grepl("IMAP", Dermo_Tol_full_module_apop_df$product, ignore.case = TRUE),]
Dermo_Tol_full_module_apop_df_GIMAP <- Dermo_Tol_full_module_apop_df_GIMAP[order(Dermo_Tol_full_module_apop_df_GIMAP$product),]

setdiff(Dermo_Tol_module_apop_df_IAP$transcript_id,Dermo_Tol_full_module_apop_df_IAP$transcript_id)
setdiff(Dermo_Tol_full_module_apop_df_IAP$transcript_id,Dermo_Tol_module_apop_df_IAP$transcript_id)


#### RUNNING PRO_RE22 FULL DATASET ####

# Pick soft threshold
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function, following general recommendations to set network type to "signed hybrid" and using the "bicor" correlation
#Pro_RE22_dds_rlog_matrix_RE22_sft <- pickSoftThreshold(Pro_RE22_dds_rlog_matrix_RE22, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
#save(Pro_RE22_dds_rlog_matrix_RE22_sft, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_dds_rlog_matrix_RE22_sft")
#Dermo
# Scale-free topology fit index as a function of the soft-thresholding power
plot(Pro_RE22_dds_rlog_matrix_RE22_sft$fitIndices[,1], -sign(Pro_RE22_dds_rlog_matrix_RE22_sft$fitIndices[,3])*Pro_RE22_dds_rlog_matrix_RE22_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(Pro_RE22_dds_rlog_matrix_RE22_sft$fitIndices[,1], -sign(Pro_RE22_dds_rlog_matrix_RE22_sft$fitIndices[,3])*Pro_RE22_dds_rlog_matrix_RE22_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(Pro_RE22_dds_rlog_matrix_RE22_sft$fitIndices[,1], Pro_RE22_dds_rlog_matrix_RE22_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(Pro_RE22_dds_rlog_matrix_RE22_sft$fitIndices[,1], Pro_RE22_dds_rlog_matrix_RE22_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Softthreshold of 9 since low number of samples and doesnt fit scale free topology 

## ONE STEP NETWORK CONSTRUCTION, MODULE DETECTION, MODULE DENDROGRAM INSPECTION ##
# Pro_RE22_RE22_full_net <- load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_RE22_full_net.RData")
#Pro_RE22_RE22_full_net = blockwiseModules(Pro_RE22_dds_rlog_matrix_RE22, power = 9, # picked suitable power in the code above 
#                                 TOMType = "signed", # use signed TOM type
#                                 networkType= "signed hybrid", # use signed hybrid network type
#                                 corType = "bicor", # use suggested bicor
#                                 TminModuleSize = 30, # recommended default
#                                 reassignThreshold = 0, # recommended default
#                                 mergeCutHeight = 0.25, # recommended default
#                                 numericLabels = TRUE, # recommended default
#                                 pamRespectsDendro = FALSE,# recommended default
#                                 verbose = 3, 
#                                 maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
# 
# How many modules identified
table(Pro_RE22_RE22_full_net$colors)
# Plot dendrogram with colors
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
Pro_RE22_RE22_full_net_mergedColors = labels2colors(Pro_RE22_RE22_full_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(Pro_RE22_RE22_full_net$dendrograms[[1]], Pro_RE22_RE22_full_net_mergedColors[Pro_RE22_RE22_full_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
Pro_RE22_RE22_full_net_moduleLabels = Pro_RE22_RE22_full_net$colors
Pro_RE22_RE22_full_net_moduleColors = labels2colors(Pro_RE22_RE22_full_net$colors)
Pro_RE22_RE22_full_net_MEs = Pro_RE22_RE22_full_net$MEs
Pro_RE22_RE22_full_net_geneTree = Pro_RE22_RE22_full_net$dendrograms[[1]]
# save network
#save(Pro_RE22_RE22_full_net, Pro_RE22_RE22_full_net_mergedColors, Pro_RE22_RE22_full_net_moduleLabels, Pro_RE22_RE22_full_net_moduleColors, Pro_RE22_RE22_full_net_MEs, Pro_RE22_RE22_full_net_geneTree, 
#     file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_RE22_full_net.RData")

# Quantify module trait associations 
# Define numbers of genes and samples
Pro_RE22_RE22_full_nGenes = ncol(Pro_RE22_dds_rlog_matrix_RE22)
Pro_RE22_RE22_full_nSamples = nrow(Pro_RE22_dds_rlog_matrix_RE22)

# Recalculate MEs with color labels
Pro_RE22_RE22_full_MEs0 = moduleEigengenes(Pro_RE22_dds_rlog_matrix_RE22, Pro_RE22_RE22_full_net_moduleColors)$eigengenes
Pro_RE22_RE22_full_MEs = orderMEs(Pro_RE22_RE22_full_MEs0)
Pro_RE22_RE22_full_moduleTraitCor = cor(Pro_RE22_RE22_full_MEs, Pro_RE22_coldata_RE22_collapsed_binarize, use = "p");
Pro_RE22_RE22_full_moduleTraitPvalue = corPvalueStudent(Pro_RE22_RE22_full_moduleTraitCor, Pro_RE22_RE22_full_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trai
sizeGrWindow(10,6)
# Will display correlations and their p-values
Pro_RE22_RE22_full_textMatrix = paste(signif(Pro_RE22_RE22_full_moduleTraitCor, 2), "\n(",
                                  signif(Pro_RE22_RE22_full_moduleTraitPvalue, 1), ")", sep = "");
dim(Pro_RE22_RE22_full_textMatrix) = dim(Pro_RE22_RE22_full_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = Pro_RE22_RE22_full_moduleTraitCor,
               xLabels = names(Pro_RE22_coldata_RE22_collapsed_binarize),
               yLabels = names(Pro_RE22_RE22_full_MEs),
               ySymbols = names(Pro_RE22_RE22_full_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = Pro_RE22_RE22_full_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
Pro_RE22_RE22_full_moduleTraitCor_df <- as.data.frame(Pro_RE22_RE22_full_moduleTraitCor)
Pro_RE22_RE22_full_moduleTraitCor_df$mod_names <- row.names(Pro_RE22_RE22_full_moduleTraitCor_df)
Pro_RE22_RE22_full_moduleTraitCor_df <- Pro_RE22_RE22_full_moduleTraitCor_df[,c("mod_names","Condition.Vibrio_coralliilyticus_RE22_exposure_6h.vs.Control_no_treatment")]
Pro_RE22_RE22_full_moduleTraitPvalue_df <- as.data.frame(Pro_RE22_RE22_full_moduleTraitPvalue)
Pro_RE22_RE22_full_moduleTraitPvalue_df$mod_names <- row.names(Pro_RE22_RE22_full_moduleTraitPvalue_df)
Pro_RE22_RE22_full_moduleTraitPvalue_df <- Pro_RE22_RE22_full_moduleTraitPvalue_df[,c("mod_names","Condition.Vibrio_coralliilyticus_RE22_exposure_6h.vs.Control_no_treatment")]
colnames(Pro_RE22_RE22_full_moduleTraitPvalue_df)[2] <- "pvalue"

Pro_RE22_RE22_full_moduleTraitCor_Pval_df <- join(Pro_RE22_RE22_full_moduleTraitCor_df,Pro_RE22_RE22_full_moduleTraitPvalue_df, by = "mod_names")

# Significantly correlated modules
Pro_RE22_RE22_full_moduleTraitCor_Pval_df[order(Pro_RE22_RE22_full_moduleTraitCor_Pval_df$pvalue),]
class(Pro_RE22_RE22_full_moduleTraitCor_Pval_df$pvalue) # numeric
Pro_RE22_RE22_full_moduleTraitCor_Pval_df_sig <- Pro_RE22_RE22_full_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.055)
Pro_RE22_RE22_full_moduleTraitCor_Pval_df_sig # 4

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
Pro_RE22_RE22_full_moduleTraitCor_Pval_df_sig_list <- Pro_RE22_RE22_full_moduleTraitCor_Pval_df_sig$mod_names
Pro_RE22_RE22_full_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Pro_RE22_RE22_full_moduleTraitCor_Pval_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
matrix_common= Pro_RE22_dds_rlog_matrix_RE22
moduleColors=Pro_RE22_RE22_full_net_moduleColors
lookup =   C_vir_rtracklayer_apop_product_final

lookup_mod_apop <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$ID %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
names(Pro_RE22_RE22_full_moduleTraitCor_Pval_df_sig_list_rm) <- c("darkgreen",  "tan" ,       "turquoise",  "darkorange" )

Pro_RE22_RE22_full_module_apop <- lapply(Pro_RE22_RE22_full_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop)
Pro_RE22_RE22_full_module_apop_df <- do.call(rbind,Pro_RE22_RE22_full_module_apop)
Pro_RE22_RE22_full_module_apop_df$mod_names <- gsub("\\..*","",row.names(Pro_RE22_RE22_full_module_apop_df))
Pro_RE22_RE22_full_module_apop_df$mod_names <- gsub("^","ME",Pro_RE22_RE22_full_module_apop_df$mod_names)
# add module significance
Pro_RE22_RE22_full_module_apop_df <- left_join(Pro_RE22_RE22_full_module_apop_df,Pro_RE22_RE22_full_moduleTraitCor_Pval_df_sig)
Pro_RE22_RE22_full_module_apop_df$exp <- "Pro_RE22_RE22_full"

## Gene relationship to trait and important modules: Gene Significance and Module Membership
# Define variable injected 
Pro_RE22_RE22_full_treatment = as.data.frame(Pro_RE22_coldata_RE22_collapsed_binarize$Condition.Vibrio_coralliilyticus_RE22_exposure_6h.vs.Control_no_treatment);
names(Pro_RE22_RE22_full_treatment) = "RE22"
# names (colors) of the modules
Pro_RE22_RE22_full_modNames = substring(names(Pro_RE22_RE22_full_MEs), 3)
Pro_RE22_RE22_full_geneModuleMembership = as.data.frame(cor(Pro_RE22_dds_rlog_matrix_RE22, Pro_RE22_RE22_full_MEs, use = "p"))
Pro_RE22_RE22_full_MMPvalue = as.data.frame(corPvalueStudent(as.matrix(Pro_RE22_RE22_full_geneModuleMembership),Pro_RE22_RE22_full_nSamples))

names(Pro_RE22_RE22_full_geneModuleMembership) = paste("MM", Pro_RE22_RE22_full_modNames, sep="")
names(Pro_RE22_RE22_full_MMPvalue) = paste("p.MM", Pro_RE22_RE22_full_modNames, sep="")

Pro_RE22_RE22_full_geneTraitSignificance = as.data.frame(cor(Pro_RE22_dds_rlog_matrix_RE22,Pro_RE22_RE22_full_treatment, use = "p"))
Pro_RE22_RE22_full_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(Pro_RE22_RE22_full_geneTraitSignificance), Pro_RE22_RE22_full_nSamples))

names(Pro_RE22_RE22_full_geneTraitSignificance) = paste("GS.", names(Pro_RE22_RE22_full_treatment), sep="")
names(Pro_RE22_RE22_full_GSPvalue) = paste("p.GS.", names(Pro_RE22_RE22_full_treatment), sep="")

## Intramodular analysis: identifying genes with high GS and MM
# Using the GS and MM measures, we can identify genes that have a high significance for disease challenge
# as well as high module membership in interesting modules. As an example, we look at the brown module 
# that has the highest association with weight. We plot a scatterplot of Gene Significance vs. Module Membership in the brown module:
Pro_RE22_RE22_full_module = "turquoise" #  0.88 very high correlation!
Pro_RE22_RE22_full_column = match(Pro_RE22_RE22_full_module, Pro_RE22_RE22_full_modNames)
Pro_RE22_RE22_full_moduleGenes = Pro_RE22_RE22_full_net_moduleColors==Pro_RE22_RE22_full_module
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(Pro_RE22_RE22_full_geneModuleMembership[Pro_RE22_RE22_full_moduleGenes, Pro_RE22_RE22_full_column]),
                   abs(Pro_RE22_RE22_full_geneTraitSignificance[Pro_RE22_RE22_full_moduleGenes, 1]),
                   xlab = paste("Module Membership in", Pro_RE22_RE22_full_module, "module"),
                   ylab = "Gene significance for challenge",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = Dermo_Tol_full_module)

## IDENTIFY HUB GENES IN EACH SIG MODULE ##
Pro_RE22_RE22_full_colorh = c("darkgreen",  "tan" ,       "turquoise",  "darkorange")

Pro_RE22_RE22_full_Module_hub_genes <- chooseTopHubInEachModule(
  Pro_RE22_dds_rlog_matrix_RE22, 
  Pro_RE22_RE22_full_colorh, 
  power = 9,  # power used for the adjacency network
  type = "signed hybrid", 
  corFnc = "bicor"
)
class(Pro_RE22_RE22_full_Module_hub_genes)
Pro_RE22_RE22_full_Module_hub_genes_df <- as.data.frame(Pro_RE22_RE22_full_Module_hub_genes)
colnames(Pro_RE22_RE22_full_Module_hub_genes_df)[1] <- "ID"
Pro_RE22_RE22_full_Module_hub_genes_apop <- merge(Pro_RE22_RE22_full_Module_hub_genes_df, C_vir_rtracklayer, by = "ID")
nrow(Pro_RE22_RE22_full_Module_hub_genes_apop) # 4, none involed in apoptosis

## Export modules to cytoscape for visualization ###
# Recalculate topological overlap if needed
#Pro_RE22_RE22_full_TOM = TOMsimilarityFromExpr(Pro_RE22_dds_rlog_matrix_RE22,
#                                           power = 9, # picked suitable power in the code above 
#                                           TOMType = "signed", # use signed TOM type
#                                           networkType= "signed hybrid", # use signed hybrid network type
#                                           corType = "bicor") # use suggested bicor
#
#save(Pro_RE22_RE22_full_TOM, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_RE22_full_TOM.RData" )

# Select modules
Pro_RE22_RE22_full_modules = c("darkgreen",  "tan" ,       "turquoise",  "darkorange" )
# Select module probes
Pro_RE22_RE22_full_probes = colnames(Pro_RE22_dds_rlog_matrix_RE22)
# export moduleColors file for use in cluster
write.table(Pro_RE22_RE22_full_net_moduleColors, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_RE22_full_net_moduleColors.table")

Pro_RE22_RE22_full_inModule = is.finite(match(Pro_RE22_RE22_full_net_moduleColors, Pro_RE22_RE22_full_modules))
Pro_RE22_RE22_full_modProbes = Pro_RE22_RE22_full_probes[Pro_RE22_RE22_full_inModule]
Pro_RE22_RE22_full_modGenes = C_vir_rtracklayer$ID[match(Pro_RE22_RE22_full_modProbes, C_vir_rtracklayer$ID)]
# Select the corresponding Topological Overlap
Pro_RE22_RE22_full_modTOM = Pro_RE22_RE22_full_TOM[Pro_RE22_RE22_full_inModule, Pro_RE22_RE22_full_inModule]

dimnames(Pro_RE22_RE22_full_modTOM ) = list(Pro_RE22_RE22_full_modProbes, Pro_RE22_RE22_full_modProbes)
# Export the network into edge and node list files Cytoscape can read
Pro_RE22_RE22_full_cyt = exportNetworkToCytoscape(Pro_RE22_RE22_full_modTOM,
                                              edgeFile = paste("CytoscapeInput-edges-Pro_RE22_RE22_full", paste(Pro_RE22_RE22_full_modules, collapse="-"), ".txt", sep=""),
                                              nodeFile = paste("CytoscapeInput-nodes-Pro_RE22_RE22_full", paste(Pro_RE22_RE22_full_modules, collapse="-"), ".txt", sep=""),
                                              weighted = TRUE,
                                              threshold = 0.02,
                                              nodeNames = Pro_RE22_RE22_full_modProbes,
                                              altNodeNames = Pro_RE22_RE22_full_modGenes,
                                              nodeAttr = Pro_RE22_RE22_full_net_moduleColors[Pro_RE22_RE22_full_inModule])


#### RUNNING ALL CVIRGINICA ON "FULL" SET OF TRANSCRIPTS ####

### SELECT SOFT THRESHOLDING POWER FOR EACH EXPERIMENT ###
# Choose a set of soft-thresholding powers
#powers = c(c(1:10), seq(from = 12, to=20, by=2))
#
## Call the network topology analysis function, following general recommendations to set network type to "signed hybrid" and using the "bicor" correlation
#Dermo_Susceptible_dds_vst_matrix_sft <- pickSoftThreshold(Dermo_Susceptible_dds_vst_matrix, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
#Probiotic_dds_rlog_matrix_sft <- pickSoftThreshold(Probiotic_dds_rlog_matrix, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
## # Warning message:
## # executing %dopar% sequentially: no parallel backend registered 
## #From Peter Langfelder: https://bioinformatics.stackexchange.com/questions/10555/r-wgcna-error-code:
## #What you see is a warning, not an error. 
## #Your calculation will run fine, just slower. Unless you see other errors, you should be able to complete all steps of the analysis.
## 
#ROD_Resistant_dds_rlog_matrix_sft <- pickSoftThreshold(ROD_Resistant_dds_rlog_matrix, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
#ROD_Susceptible_dds_rlog_matrix_sft <- pickSoftThreshold(ROD_Susceptible_dds_rlog_matrix, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
#Pro_RE22_dds_rlog_matrix_Pro_sft <- pickSoftThreshold(Pro_RE22_dds_rlog_matrix_Pro, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
#
## export pick softthresholding so it is saved 
#save(Dermo_Susceptible_dds_vst_matrix_sft, Probiotic_dds_rlog_matrix_sft, ROD_Resistant_dds_rlog_matrix_sft,
# ROD_Susceptible_dds_rlog_matrix_sft, Pro_RE22_dds_rlog_matrix_Pro_sft, Pro_RE22_dds_rlog_matrix_RE22_sft, file= "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_virginica_FULL_individual_pickSoftThreshold_WGCNA_input.RData")
#
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

#Dermo
# Scale-free topology fit index as a function of the soft-thresholding power
plot(Dermo_Susceptible_dds_vst_matrix_sft$fitIndices[,1], -sign(Dermo_Susceptible_dds_vst_matrix_sft$fitIndices[,3])*Dermo_Susceptible_dds_vst_matrix_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(Dermo_Susceptible_dds_vst_matrix_sft$fitIndices[,1], -sign(Dermo_Susceptible_dds_vst_matrix_sft$fitIndices[,3])*Dermo_Susceptible_dds_vst_matrix_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(Dermo_Susceptible_dds_vst_matrix_sft$fitIndices[,1], Dermo_Susceptible_dds_vst_matrix_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(Dermo_Susceptible_dds_vst_matrix_sft$fitIndices[,1], Dermo_Susceptible_dds_vst_matrix_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Softthreshold of 6 is more reasonable

# Probiotic, Probiotic_counts_apop_dds_rlog_matrix_sft: Does not fit scale free topology!
# Scale-free topology fit index as a function of the soft-thresholding power
plot(Probiotic_dds_rlog_matrix_sft$fitIndices[,1], -sign(Probiotic_dds_rlog_matrix_sft$fitIndices[,3])*Probiotic_dds_rlog_matrix_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(Probiotic_dds_rlog_matrix_sft$fitIndices[,1], -sign(Probiotic_dds_rlog_matrix_sft$fitIndices[,3])*Probiotic_dds_rlog_matrix_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(Probiotic_dds_rlog_matrix_sft$fitIndices[,1],Probiotic_dds_rlog_matrix_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(Probiotic_dds_rlog_matrix_sft$fitIndices[,1], Probiotic_dds_rlog_matrix_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Poor fit of scale free topology: see recommendations here: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html
# This is likely caused by strong clustering by Time, using 9 
# Because less than 20 samples, selecting soft threshold of 9 for signed hybrid network 

# ROD pick soft threshold
# Scale-free topology fit index as a function of the soft-thresholding power
plot(ROD_Resistant_dds_rlog_matrix_sft $fitIndices[,1], -sign(ROD_Resistant_dds_rlog_matrix_sft $fitIndices[,3])*ROD_Resistant_dds_rlog_matrix_sft $fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(ROD_Resistant_dds_rlog_matrix_sft $fitIndices[,1], -sign(ROD_Resistant_dds_rlog_matrix_sft $fitIndices[,3])*ROD_Resistant_dds_rlog_matrix_sft $fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(ROD_Resistant_dds_rlog_matrix_sft $fitIndices[,1],ROD_Resistant_dds_rlog_matrix_sft $fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(ROD_Resistant_dds_rlog_matrix_sft $fitIndices[,1], ROD_Resistant_dds_rlog_matrix_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Same story as above, selecting soft threshold of 9 because less than 20 samples

# Scale-free topology fit index as a function of the soft-thresholding power
plot(  ROD_Susceptible_dds_rlog_matrix_sft$fitIndices[,1], -sign(  ROD_Susceptible_dds_rlog_matrix_sft$fitIndices[,3])*  ROD_Susceptible_dds_rlog_matrix_sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
text( ROD_Susceptible_dds_rlog_matrix_sft$fitIndices[,1], -sign(ROD_Susceptible_dds_rlog_matrix_sft$fitIndices[,3])*ROD_Susceptible_dds_rlog_matrix_sft$fitIndices[,2],
      labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(ROD_Susceptible_dds_rlog_matrix_sft$fitIndices[,1],ROD_Susceptible_dds_rlog_matrix_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(ROD_Susceptible_dds_rlog_matrix_sft$fitIndices[,1], ROD_Susceptible_dds_rlog_matrix_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# no scale free topology, selecting 9 

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Pro_RE22_Pro
# Scale-free topology fit index as a function of the soft-thresholding power
plot(Pro_RE22_dds_rlog_matrix_Pro_sft$fitIndices[,1], -sign(Pro_RE22_dds_rlog_matrix_Pro_sft$fitIndices[,3])*Pro_RE22_dds_rlog_matrix_Pro_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(Pro_RE22_dds_rlog_matrix_Pro_sft$fitIndices[,1], -sign(Pro_RE22_dds_rlog_matrix_Pro_sft$fitIndices[,3])*Pro_RE22_dds_rlog_matrix_Pro_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(Pro_RE22_dds_rlog_matrix_Pro_sft$fitIndices[,1],Pro_RE22_dds_rlog_matrix_Pro_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(Pro_RE22_dds_rlog_matrix_Pro_sft$fitIndices[,1], Pro_RE22_dds_rlog_matrix_Pro_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# selecting scale of 4 

# Pro_RE22_RE22
# Scale-free topology fit index as a function of the soft-thresholding power
plot(Pro_RE22_dds_rlog_matrix_RE22_sft$fitIndices[,1], -sign(Pro_RE22_dds_rlog_matrix_RE22_sft$fitIndices[,3])*Pro_RE22_dds_rlog_matrix_RE22_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(Pro_RE22_dds_rlog_matrix_RE22_sft$fitIndices[,1], -sign(Pro_RE22_dds_rlog_matrix_RE22_sft$fitIndices[,3])*Pro_RE22_dds_rlog_matrix_RE22_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(Pro_RE22_dds_rlog_matrix_RE22_sft$fitIndices[,1],Pro_RE22_dds_rlog_matrix_RE22_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(Pro_RE22_dds_rlog_matrix_RE22_sft$fitIndices[,1], Pro_RE22_dds_rlog_matrix_RE22_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# no scale free topology, selecting 9 because few samples

### ONE STEP NETWORK CONSTRUCTION, MODULE DETECTION, MODULE DENDROGRAM INSPECTION ###

# Dermo Susceptible FULL
 Dermo_Sus_full_net = blockwiseModules(Dermo_Susceptible_dds_vst_matrix, power = 6, # picked suitable power in the code above 
                                  TOMType = "signed", # use signed TOM type
                                  networkType= "signed hybrid", # use signed hybrid network type
                                  corType = "bicor", # use suggested bicor
                                  TminModuleSize = 30, # recommended default
                                  reassignThreshold = 0, # recommended default
                                  mergeCutHeight = 0.25, # recommended default
                                  numericLabels = TRUE, # recommended default
                                  pamRespectsDendro = FALSE,# recommended default
                                  verbose = 3, 
                                  maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
# # How many modules identified
table(Dermo_Sus_full_net$colors)
# Plot dendrogram with colors
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
Dermo_Sus_full_mergedColors = labels2colors(Dermo_Sus_full_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(Dermo_Sus_full_net$dendrograms[[1]], Dermo_Sus_full_mergedColors[Dermo_Sus_full_net$blockGenes[[1]]],
                   "Module colors",
                   dendroLabels = FALSE, hang = 0.03,
                   addGuide = TRUE, guideHang = 0.05)
Dermo_Sus_full_moduleLabels = Dermo_Sus_full_net$colors
Dermo_Sus_full_moduleColors = labels2colors(Dermo_Sus_full_net$colors)
Dermo_Sus_full_MEs = Dermo_Sus_full_net$MEs
Dermo_Sus_full_geneTree = Dermo_Sus_full_net$dendrograms[[1]]
save(Dermo_Sus_full_net, Dermo_Sus_full_mergedColors, Dermo_Sus_full_moduleLabels, Dermo_Sus_full_moduleColors, Dermo_Sus_full_MEs, Dermo_Sus_full_geneTree, 
    file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Sus_full_network.RData")

# Load network
#Dermo_Sus_full_network = load(file= "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Sus_full_network.RData")

# Probiotic
Probiotic_full_net = blockwiseModules(Probiotic_dds_rlog_matrix, power = 9, # picked because less than 20 samples and isn't scale free
                                 TOMType = "signed", # use signed TOM type
                                 networkType= "signed hybrid", # use signed hybrid network type
                                 corType = "bicor", # use suggested bicor
                                 TminModuleSize = 30, # recommended default
                                 reassignThreshold = 0, # recommended default
                                 mergeCutHeight = 0.25, # recommended default
                                 numericLabels = TRUE, # recommended default
                                 pamRespectsDendro = FALSE,# recommended default
                                 verbose = 3, 
                                 maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
# How many modules identified
table(Probiotic_full_net$colors)
# Plot dendrogram with colors
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
Probiotic_full_mergedColors = labels2colors(Probiotic_full_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(Probiotic_full_net$dendrograms[[1]], Probiotic_full_mergedColors[Probiotic_full_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
Probiotic_full_moduleLabels = Probiotic_full_net$colors
Probiotic_full_moduleColors = labels2colors(Probiotic_full_net$colors)
Probiotic_full_MEs = Probiotic_full_net$MEs
Probiotic_full_geneTree = Probiotic_full_net$dendrograms[[1]]
# save network externally
save(Probiotic_full_net, Probiotic_full_mergedColors, Probiotic_full_moduleLabels, Probiotic_full_moduleColors, Probiotic_full_MEs, Probiotic_full_geneTree, file=
       "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Probiotic_full_network.RData")

# ROD Res
ROD_Res_full_net = blockwiseModules(ROD_Resistant_dds_rlog_matrix, power = 9, # picked because less than 20 samples and isn't scale free
                               TOMType = "signed", # use signed TOM type
                               networkType= "signed hybrid", # use signed hybrid network type
                               corType = "bicor", # use suggested bicor
                               TminModuleSize = 30, # recommended default
                               reassignThreshold = 0, # recommended default
                               mergeCutHeight = 0.25, # recommended default
                               numericLabels = TRUE, # recommended default
                               pamRespectsDendro = FALSE,# recommended default
                               verbose = 3, 
                               maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
# How many modules identified
table(ROD_Res_full_net$colors)
# Plot dendrogram with colors
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
ROD_Res_full_mergedColors = labels2colors(ROD_Res_full_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(ROD_Res_full_net$dendrograms[[1]], ROD_Res_full_mergedColors[ROD_Res_full_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
ROD_Res_full_moduleLabels = ROD_Res_full_net$colors
ROD_Res_full_moduleColors = labels2colors(ROD_Res_full_net$colors)
ROD_Res_full_MEs = ROD_Res_full_net$MEs
ROD_Res_full_geneTree = ROD_Res_full_net$dendrograms[[1]]
save(ROD_Res_full_net, ROD_Res_full_mergedColors, ROD_Res_full_moduleLabels, ROD_Res_full_moduleColors, ROD_Res_full_MEs, ROD_Res_full_geneTree, file=
       "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/ROD_Res_full_network.RData")

# ROD Sus
ROD_Sus_full_net = blockwiseModules(ROD_Susceptible_dds_rlog_matrix, power = 9, # picked because less than 20 samples and isn't scale free
                               TOMType = "signed", # use signed TOM type
                               networkType= "signed hybrid", # use signed hybrid network type
                               corType = "bicor", # use suggested bicor
                               TminModuleSize = 30, # recommended default
                               reassignThreshold = 0, # recommended default
                               mergeCutHeight = 0.25, # recommended default
                               numericLabels = TRUE, # recommended default
                               pamRespectsDendro = FALSE,# recommended default
                               verbose = 3, 
                               maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
# How many modules identified
table(ROD_Sus_full_net$colors)
# Plot dendrogram with colors
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
ROD_Sus_full_mergedColors = labels2colors(ROD_Sus_full_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(ROD_Sus_full_net$dendrograms[[1]], ROD_Sus_full_mergedColors[ROD_Sus_full_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
ROD_Sus_full_moduleLabels = ROD_Sus_full_net$colors
ROD_Sus_full_moduleColors = labels2colors(ROD_Sus_full_net$colors)
ROD_Sus_full_MEs = ROD_Sus_full_net$MEs
ROD_Sus_full_geneTree = ROD_Sus_full_net$dendrograms[[1]]

save(ROD_Sus_full_net, ROD_Sus_full_mergedColors, ROD_Sus_full_moduleLabels, ROD_Sus_full_moduleColors, ROD_Sus_full_MEs, ROD_Sus_full_geneTree, file=
       "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/ROD_Sus_full_network.RData")

# Pro_RE22_Pro
Pro_RE22_Pro_full_net = blockwiseModules(Pro_RE22_dds_rlog_matrix_Pro, power = 4, 
                                    TOMType = "signed", # use signed TOM type
                                    networkType= "signed hybrid", # use signed hybrid network type
                                    corType = "bicor", # use suggested bicor
                                    TminModuleSize = 30, # recommended default
                                    reassignThreshold = 0, # recommended default
                                    mergeCutHeight = 0.25, # recommended default
                                    numericLabels = TRUE, # recommended default
                                    pamRespectsDendro = FALSE,# recommended default
                                    verbose = 3, 
                                    maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
# How many modules identified
table(Pro_RE22_Pro_full_net$colors)
# Plot dendrogram with colors
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
Pro_RE22_Pro_full_mergedColors = labels2colors(Pro_RE22_Pro_full_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(Pro_RE22_Pro_full_net$dendrograms[[1]], Pro_RE22_Pro_full_mergedColors[Pro_RE22_Pro_full_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
Pro_RE22_Pro_full_moduleLabels = Pro_RE22_Pro_full_net$colors
Pro_RE22_Pro_full_moduleColors = labels2colors(Pro_RE22_Pro_full_net$colors)
Pro_RE22_Pro_full_MEs = Pro_RE22_Pro_full_net$MEs
Pro_RE22_Pro_full_geneTree = Pro_RE22_Pro_full_net$dendrograms[[1]]

save(Pro_RE22_Pro_full_net, Pro_RE22_Pro_full_mergedColors, Pro_RE22_Pro_full_moduleLabels, Pro_RE22_Pro_full_moduleColors, Pro_RE22_Pro_full_MEs, Pro_RE22_Pro_full_geneTree, file=
       "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_Pro_full_network.RData")

### QUANTIFYING MODULE TRAIT ASSOCIATIONS ###
### DERMO TOL ALREADY DONE ABOVE ###

### DERMO SUS ###
# Define numbers of genes and samples
Dermo_Sus_full_nGenes = ncol(Dermo_Susceptible_dds_vst_matrix)
Dermo_Sus_full_nSamples = nrow(Dermo_Susceptible_dds_vst_matrix)

# Recalculate MEs with color labels
Dermo_Sus_full_MEs0 = moduleEigengenes(Dermo_Susceptible_dds_vst_matrix, Dermo_Sus_full_moduleColors)$eigengenes
Dermo_Sus_full_MEs = orderMEs(Dermo_Sus_full_MEs0)
Dermo_Sus_full_moduleTraitCor = cor(Dermo_Sus_full_MEs, Dermo_Susceptible_coldata_collapsed_binarize, use = "p");
Dermo_Sus_full_moduleTraitPvalue = corPvalueStudent(Dermo_Sus_full_moduleTraitCor, Dermo_Sus_full_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trai
sizeGrWindow(10,6)
# Will display correlations and their p-values
Dermo_Sus_full_textMatrix = paste(signif(Dermo_Sus_full_moduleTraitCor, 2), "\n(",
                             signif(Dermo_Sus_full_moduleTraitPvalue, 1), ")", sep = "");
dim(Dermo_Sus_full_textMatrix) = dim(Dermo_Sus_full_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = Dermo_Sus_full_moduleTraitCor,
               xLabels = names(Dermo_Susceptible_coldata_collapsed_binarize),
               yLabels = names(Dermo_Sus_full_MEs),
               ySymbols = names(Dermo_Sus_full_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = Dermo_Sus_full_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
Dermo_Sus_full_moduleTraitCor_df <- as.data.frame(Dermo_Sus_full_moduleTraitCor)
Dermo_Sus_full_moduleTraitCor_df$mod_names <- row.names(Dermo_Sus_full_moduleTraitCor_df)
Dermo_Sus_full_moduleTraitCor_df <- Dermo_Sus_full_moduleTraitCor_df[,c("mod_names","Condition.Injected.vs.Control")]
Dermo_Sus_full_moduleTraitPvalue_df <- as.data.frame(Dermo_Sus_full_moduleTraitPvalue)
Dermo_Sus_full_moduleTraitPvalue_df$mod_names <- row.names(Dermo_Sus_full_moduleTraitPvalue_df)
Dermo_Sus_full_moduleTraitPvalue_df <- Dermo_Sus_full_moduleTraitPvalue_df[,c("mod_names","Condition.Injected.vs.Control")]
colnames(Dermo_Sus_full_moduleTraitPvalue_df)[2] <- "pvalue"

Dermo_Sus_full_moduleTraitCor_Pval_df <- join(Dermo_Sus_full_moduleTraitCor_df, Dermo_Sus_full_moduleTraitPvalue_df, by = "mod_names")

# Significantly correlated modules
Dermo_Sus_full_moduleTraitCor_Pval_df[order(Dermo_Sus_full_moduleTraitCor_Pval_df$pvalue),]
class(Dermo_Sus_full_moduleTraitCor_Pval_df$pvalue) # numeric
Dermo_Sus_full_moduleTraitCor_Pval_df_sig <- Dermo_Sus_full_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.055)
Dermo_Sus_full_moduleTraitCor_Pval_df_sig # 4 significant models 

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
Dermo_Sus_full_moduleTraitCor_Pval_df_sig_list <- Dermo_Sus_full_moduleTraitCor_Pval_df_sig$mod_names
Dermo_Sus_full_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Dermo_Sus_full_moduleTraitCor_Pval_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
matrix_common= Dermo_Susceptible_dds_vst_matrix
moduleColors=Dermo_Sus_full_moduleColors
lookup =   C_vir_rtracklayer_apop_product_final

lookup_mod_apop <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$ID %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")] 
}
# specify names for list of lists
names(Dermo_Sus_full_moduleTraitCor_Pval_df_sig_list_rm) <- c("ivory",       "skyblue" ,    "lightpink4" , "lightyellow")
Dermo_Sus_full_module_apop <- lapply(Dermo_Sus_full_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop)
Dermo_Sus_full_module_apop_df <- do.call(rbind,Dermo_Sus_full_module_apop)
Dermo_Sus_full_module_apop_df$mod_names <- gsub("\\..*","",row.names(Dermo_Sus_full_module_apop_df))
Dermo_Sus_full_module_apop_df$mod_names <- gsub("^","ME",Dermo_Sus_full_module_apop_df$mod_names)
# add module significance
Dermo_Sus_full_module_apop_df <- left_join(Dermo_Sus_full_module_apop_df,Dermo_Sus_full_moduleTraitCor_Pval_df_sig)
Dermo_Sus_full_module_apop_df$exp <- "Dermo_Sus"

## Gene relationship to trait and important modules: Gene Significance and Module Membership
# Define variable injected 
Dermo_Sus_full_injection = as.data.frame(Dermo_Susceptible_coldata_collapsed_binarize$Condition.Injected.vs.Control);
names(Dermo_Sus_full_injection) = "injection"
# names (colors) of the modules
Dermo_Sus_full_modNames = substring(names(Dermo_Sus_full_MEs), 3)
Dermo_Sus_full_geneModuleMembership = as.data.frame(cor(Dermo_Susceptible_dds_vst_matrix, Dermo_Sus_full_MEs, use = "p"))
Dermo_Sus_full_MMPvalue = as.data.frame(corPvalueStudent(as.matrix(Dermo_Sus_full_geneModuleMembership), Dermo_Sus_full_nSamples))

names(Dermo_Sus_full_geneModuleMembership) = paste("MM", Dermo_Sus_full_modNames, sep="")
names(Dermo_Sus_full_MMPvalue) = paste("p.MM", Dermo_Sus_full_modNames, sep="")

Dermo_Sus_full_geneTraitSignificance = as.data.frame(cor(Dermo_Susceptible_dds_vst_matrix,Dermo_Sus_full_injection, use = "p"))
Dermo_Sus_full_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(Dermo_Sus_full_geneTraitSignificance), Dermo_Sus_full_nSamples))

names(Dermo_Sus_full_geneTraitSignificance) = paste("GS.", names(Dermo_Sus_full_injection), sep="")
names(Dermo_Sus_full_GSPvalue) = paste("p.GS.", names(Dermo_Sus_full_injection), sep="")

## Intramodular analysis: identifying genes with high GS and MM
# Using the GS and MM measures, we can identify genes that have a high significance for disease challenge
# as well as high module membership in interesting modules. As an example, we look at the brown module 
# that has the highest association with weight. We plot a scatterplot of Gene Significance vs. Module Membership in the brown module:
Dermo_Sus_full_module = "ivory" # strong correlation
Dermo_Sus_full_column = match(Dermo_Sus_full_module, Dermo_Sus_full_modNames)
Dermo_Sus_full_moduleGenes = Dermo_Sus_full_moduleColors==Dermo_Sus_full_module
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(Dermo_Sus_full_geneModuleMembership[Dermo_Sus_full_moduleGenes, Dermo_Sus_full_column]),
                   abs(Dermo_Sus_full_geneTraitSignificance[Dermo_Sus_full_moduleGenes, 1]),
                   xlab = paste("Module Membership in", Dermo_Sus_full_module, "module"),
                   ylab = "Gene significance for challenge",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = Dermo_Sus_full_module)

## IDENTIFY HUB GENES IN EACH SIG MODULE ##
Dermo_Sus_full_colorh = c("MEivory" ,      "MEskyblue"  ,   "MElightpink4"  ,"MElightyellow")

#Dermo_Sus_full_Module_hub_genes <- chooseTopHubInEachModule(
#  Dermo_Susceptible_dds_vst_matrix, 
#  Dermo_Sus_full_colorh, 
#  power = 6,  # power used for the adjacency network
#  type = "signed hybrid", 
#  corFnc = "bicor"
#)
class(Dermo_Sus_full_Module_hub_genes)
Dermo_Sus_full_Module_hub_genes_df <- as.data.frame(Dermo_Sus_full_Module_hub_genes)
colnames(Dermo_Sus_full_Module_hub_genes_df)[1] <- "ID"
Dermo_Sus_full_Module_hub_genes_apop <- merge(Dermo_Sus_full_Module_hub_genes_df, C_vir_rtracklayer, by = "ID")
nrow(Dermo_Sus_full_Module_hub_genes_apop) 

### PROBIOTIC ###
# Define numbers of genes and samples
Probiotic_full_nGenes = ncol(Probiotic_dds_rlog_matrix)
Probiotic_full_nSamples = nrow(Probiotic_dds_rlog_matrix)

# Recalculate MEs with color labels
Probiotic_full_MEs0 = moduleEigengenes(Probiotic_dds_rlog_matrix, Probiotic_full_moduleColors)$eigengenes
Probiotic_full_MEs = orderMEs(Probiotic_full_MEs0)
Probiotic_full_moduleTraitCor = cor(Probiotic_full_MEs, Probiotic_coldata_collapsed_binarize , use = "p");
Probiotic_full_moduleTraitPvalue = corPvalueStudent(Probiotic_full_moduleTraitCor, Probiotic_full_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trait
sizeGrWindow(10,6)
# Will display correlations and their p-values
Probiotic_full_textMatrix = paste(signif(Probiotic_full_moduleTraitCor, 2), "\n(",
                             signif(Probiotic_full_moduleTraitPvalue, 1), ")", sep = "");
dim(Probiotic_full_textMatrix) = dim(Probiotic_full_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = Probiotic_full_moduleTraitCor,
               xLabels = names(Probiotic_coldata_collapsed_binarize),
               yLabels = names(Probiotic_full_MEs),
               ySymbols = names(Probiotic_full_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = Probiotic_full_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
Probiotic_full_moduleTraitCor_df <- as.data.frame(Probiotic_full_moduleTraitCor)
Probiotic_full_moduleTraitCor_df$mod_names <- row.names(Probiotic_full_moduleTraitCor_df)
Probiotic_full_moduleTraitCor_df <- Probiotic_full_moduleTraitCor_df[,c("mod_names","Condition.Bacillus_pumilus_RI0695.vs.Untreated_control")]
Probiotic_full_moduleTraitPvalue_df <- as.data.frame(Probiotic_full_moduleTraitPvalue)
Probiotic_full_moduleTraitPvalue_df$mod_names <- row.names(Probiotic_full_moduleTraitPvalue_df)
Probiotic_full_moduleTraitPvalue_df <- Probiotic_full_moduleTraitPvalue_df[,c("mod_names","Condition.Bacillus_pumilus_RI0695.vs.Untreated_control")]
colnames(Probiotic_full_moduleTraitPvalue_df)[2] <- "pvalue"
Probiotic_full_moduleTraitCor_Pval_df <- join(Probiotic_full_moduleTraitCor_df, Probiotic_full_moduleTraitPvalue_df, by = "mod_names")

# Significantly correlated modules
Probiotic_full_moduleTraitCor_Pval_df[order(Probiotic_full_moduleTraitCor_Pval_df$pvalue),]
class(Probiotic_full_moduleTraitCor_Pval_df$pvalue) # numeric
Probiotic_full_moduleTraitCor_Pval_df_sig <- Probiotic_full_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
Probiotic_full_moduleTraitCor_Pval_df_sig # 3 significant modules

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
Probiotic_full_moduleTraitCor_Pval_df_sig_list <- Probiotic_full_moduleTraitCor_Pval_df_sig$mod_names
Probiotic_full_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Probiotic_full_moduleTraitCor_Pval_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
matrix_common= Probiotic_dds_rlog_matrix
moduleColors= Probiotic_full_moduleColors
lookup =   C_vir_rtracklayer_apop_product_final

lookup_mod_apop <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$ID %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
names(Probiotic_full_moduleTraitCor_Pval_df_sig_list_rm) <- c("purple",      "lightcyan" ,  "saddlebrown")
Probiotic_full_module_apop <- lapply(Probiotic_full_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop)
Probiotic_full_module_apop_df <- do.call(rbind,Probiotic_full_module_apop)
Probiotic_full_module_apop_df$mod_names <- gsub("\\..*","",row.names(Probiotic_full_module_apop_df))
Probiotic_full_module_apop_df$mod_names <- gsub("^","ME",Probiotic_full_module_apop_df$mod_names)
# add module significance
Probiotic_full_module_apop_df <- left_join(Probiotic_full_module_apop_df,Probiotic_full_moduleTraitCor_Pval_df_sig)
Probiotic_full_module_apop_df$exp <- "Probiotic"

## Gene relationship to trait and important modules: Gene Significance and Module Membership
# Define variable injected 
Probiotic_full_injection = as.data.frame(Probiotic_coldata_collapsed_binarize$Condition.Bacillus_pumilus_RI0695.vs.Untreated_control);
names(Probiotic_full_injection) = "challenge"
# names (colors) of the modules
Probiotic_full_modNames = substring(names(Probiotic_full_MEs), 3)
Probiotic_full_geneModuleMembership = as.data.frame(cor(Probiotic_dds_rlog_matrix, Probiotic_full_MEs, use = "p"))
Probiotic_full_MMPvalue = as.data.frame(corPvalueStudent(as.matrix(Probiotic_full_geneModuleMembership), Probiotic_full_nSamples))

names(Probiotic_full_geneModuleMembership) = paste("MM", Probiotic_full_modNames, sep="")
names(Probiotic_full_MMPvalue) = paste("p.MM", Probiotic_full_modNames, sep="")

Probiotic_full_geneTraitSignificance = as.data.frame(cor(Probiotic_dds_rlog_matrix,Probiotic_full_injection, use = "p"))
Probiotic_full_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(Probiotic_full_geneTraitSignificance), Probiotic_full_nSamples))

names(Probiotic_full_geneTraitSignificance) = paste("GS.", names(Probiotic_full_injection), sep="")
names(Probiotic_full_GSPvalue) = paste("p.GS.", names(Probiotic_full_injection), sep="")

## Intramodular analysis: identifying genes with high GS and MM
Probiotic_full_module = "lightcyan"  # cor=0.75, very strong correlation with trait
Probiotic_full_column = match(Probiotic_full_module, Probiotic_full_modNames)
Probiotic_full_moduleGenes = Probiotic_full_moduleColors==Probiotic_full_module
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(Probiotic_full_geneModuleMembership[Probiotic_full_moduleGenes, Probiotic_full_column]),
                   abs(Probiotic_full_geneTraitSignificance[Probiotic_full_moduleGenes, 1]),
                   xlab = paste("Module Membership in", Probiotic_full_module, "module"),
                   ylab = "Gene significance for challenge",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = Probiotic_full_module)

## IDENTIFY HUB GENES IN EACH SIG MODULE ##
#Probiotic_full_colorh = c("MEpurple",      "MElightcyan" ,  "MEsaddlebrown")
#
#Probiotic_full_Module_hub_genes <- chooseTopHubInEachModule(
#  Probiotic_dds_rlog_matrix, 
#  Probiotic_full_colorh, 
#  power = 9,  # power used for the adjacency network
#  type = "signed hybrid", 
#  corFnc = "bicor"
#)
#class(Probiotic_full_Module_hub_genes)
#Probiotic_full_Module_hub_genes_df <- as.data.frame(Probiotic_full_Module_hub_genes)
#colnames(Probiotic_full_Module_hub_genes_df)[1] <- "ID"
#Probiotic_full_Module_hub_genes_apop <- merge(Probiotic_full_Module_hub_genes_df, C_vir_rtracklayer, by = "ID")
#nrow(Probiotic_full_Module_hub_genes_apop) 

### ROD RES ###
# Define numbers of genes and samples
ROD_Res_full_nGenes = ncol(ROD_Resistant_dds_rlog_matrix)
ROD_Res_full_nSamples = nrow(ROD_Resistant_dds_rlog_matrix)

# Recalculate MEs with color labels
ROD_Res_full_MEs0 = moduleEigengenes(ROD_Resistant_dds_rlog_matrix, ROD_Res_full_moduleColors)$eigengenes
ROD_Res_full_MEs = orderMEs(ROD_Res_full_MEs0)
ROD_Res_full_moduleTraitCor = cor(ROD_Res_full_MEs, ROD_Resistant_coldata_collapsed_binarize , use = "p");
ROD_Res_full_moduleTraitPvalue = corPvalueStudent(ROD_Res_full_moduleTraitCor, ROD_Res_full_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trait
sizeGrWindow(10,6)
# Will display correlations and their p-values
ROD_Res_full_textMatrix = paste(signif(ROD_Res_full_moduleTraitCor, 2), "\n(",
                           signif(ROD_Res_full_moduleTraitPvalue, 1), ")", sep = "");
dim(ROD_Res_full_textMatrix) = dim(ROD_Res_full_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = ROD_Res_full_moduleTraitCor,
               xLabels = names(ROD_Resistant_coldata_collapsed_binarize),
               yLabels = names(ROD_Res_full_MEs),
               ySymbols =names(ROD_Res_full_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = ROD_Res_full_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
ROD_Res_full_moduleTraitCor_df <- as.data.frame(ROD_Res_full_moduleTraitCor)
ROD_Res_full_moduleTraitCor_df$mod_names <- row.names(ROD_Res_full_moduleTraitCor_df)
ROD_Res_full_moduleTraitCor_df <- ROD_Res_full_moduleTraitCor_df[,c("mod_names","Condition.Resistant_Challenge.vs.Control_Resistant")]
ROD_Res_full_moduleTraitPvalue_df <- as.data.frame(ROD_Res_full_moduleTraitPvalue)
ROD_Res_full_moduleTraitPvalue_df$mod_names <- row.names(ROD_Res_full_moduleTraitPvalue_df)
ROD_Res_full_moduleTraitPvalue_df <- ROD_Res_full_moduleTraitPvalue_df[,c("mod_names","Condition.Resistant_Challenge.vs.Control_Resistant")]
colnames(ROD_Res_full_moduleTraitPvalue_df)[2] <- "pvalue"
ROD_Res_full_moduleTraitCor_Pval_df <- join(ROD_Res_full_moduleTraitCor_df, ROD_Res_full_moduleTraitPvalue_df, by = "mod_names")

# Significantly correlated modules
ROD_Res_full_moduleTraitCor_Pval_df[order(ROD_Res_full_moduleTraitCor_Pval_df$pvalue),]
class(ROD_Res_full_moduleTraitCor_Pval_df$pvalue) # numeric
ROD_Res_full_moduleTraitCor_Pval_df_sig <- ROD_Res_full_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
ROD_Res_full_moduleTraitCor_Pval_df_sig # 3 significant modules

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
ROD_Res_full_moduleTraitCor_Pval_df_sig_list <- ROD_Res_full_moduleTraitCor_Pval_df_sig$mod_names
ROD_Res_full_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(ROD_Res_full_moduleTraitCor_Pval_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
matrix_common= ROD_Resistant_dds_rlog_matrix
moduleColors= ROD_Res_full_moduleColors
lookup =   C_vir_rtracklayer_apop_product_final

lookup_mod_apop <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$ID %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
names(ROD_Res_full_moduleTraitCor_Pval_df_sig_list_rm) <- c("bisque4"    ,    "brown4"      ,   "lightsteelblue1")
ROD_Res_full_module_apop <- lapply(ROD_Res_full_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop)
ROD_Res_full_module_apop_df <- do.call(rbind,ROD_Res_full_module_apop)
ROD_Res_full_module_apop_df$mod_names <- gsub("\\..*","",row.names(ROD_Res_full_module_apop_df))
ROD_Res_full_module_apop_df$mod_names <- gsub("^","ME",ROD_Res_full_module_apop_df$mod_names)
# add module significance
ROD_Res_full_module_apop_df <- left_join(ROD_Res_full_module_apop_df,ROD_Res_full_moduleTraitCor_Pval_df_sig)
ROD_Res_full_module_apop_df$exp <- "ROD_Res"

## Gene relationship to trait and important modules: Gene Significance and Module Membership
# Define variable injected 
ROD_Res_full_injection = as.data.frame(ROD_Resistant_coldata_collapsed_binarize$Condition.Resistant_Challenge.vs.Control_Resistant);
names(ROD_Res_full_injection) = "challenge"
# names (colors) of the modules
ROD_Res_full_modNames = substring(names(ROD_Res_full_MEs), 3)
ROD_Res_full_geneModuleMembership = as.data.frame(cor(ROD_Resistant_dds_rlog_matrix, ROD_Res_full_MEs, use = "p"))
ROD_Res_full_MMPvalue = as.data.frame(corPvalueStudent(as.matrix(ROD_Res_full_geneModuleMembership), ROD_Res_full_nSamples))

names(ROD_Res_full_geneModuleMembership) = paste("MM", ROD_Res_full_modNames, sep="")
names(ROD_Res_full_MMPvalue) = paste("p.MM", ROD_Res_full_modNames, sep="")

ROD_Res_full_geneTraitSignificance = as.data.frame(cor(ROD_Resistant_dds_rlog_matrix,ROD_Res_full_injection, use = "p"))
ROD_Res_full_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(ROD_Res_full_geneTraitSignificance), ROD_Res_full_nSamples))

names(ROD_Res_full_geneTraitSignificance) = paste("GS.", names(ROD_Res_full_injection), sep="")
names(ROD_Res_full_GSPvalue) = paste("p.GS.", names(ROD_Res_full_injection), sep="")

## Intramodular analysis: identifying genes with high GS and MM
ROD_Res_full_module = "bisque4"  # very strong correlation, 0.81 
ROD_Res_full_column = match(ROD_Res_full_module, ROD_Res_full_modNames)
ROD_Res_full_moduleGenes = ROD_Res_full_moduleColors==ROD_Res_full_module
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(ROD_Res_full_geneModuleMembership[ROD_Res_full_moduleGenes, ROD_Res_full_column]),
                   abs(ROD_Res_full_geneTraitSignificance[ROD_Res_full_moduleGenes, 1]),
                   xlab = paste("Module Membership in", ROD_Res_full_module, "module"),
                   ylab = "Gene significance for challenge",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = ROD_Res_full_module)

## IDENTIFY HUB GENES IN EACH SIG MODULE ##
ROD_Res_full_colorh = c("MEbisque4"  ,       "MEbrown4"    ,      "MElightsteelblue1")

#ROD_Res_full_Module_hub_genes <- chooseTopHubInEachModule(
#  ROD_Resistant_dds_rlog_matrix, 
#  ROD_Res_full_colorh, 
#  power = 9,  # power used for the adjacency network
#  type = "signed hybrid", 
#  corFnc = "bicor"
#)
class(ROD_Res_full_Module_hub_genes)
ROD_Res_full_Module_hub_genes_df <- as.data.frame(ROD_Res_full_Module_hub_genes)
colnames(ROD_Res_full_Module_hub_genes_df)[1] <- "ID"
ROD_Res_full_Module_hub_genes <- merge(ROD_Res_full_Module_hub_genes_df, C_vir_rtracklayer, by = "ID")
nrow(ROD_Res_full_Module_hub_genes) 

### ROD SUS ###
# Define numbers of genes and samples
ROD_Sus_full_nGenes = ncol(ROD_Susceptible_dds_rlog_matrix)
ROD_Sus_full_nSamples = nrow(ROD_Susceptible_dds_rlog_matrix)

# Recalculate MEs with color labels
ROD_Sus_full_MEs0 = moduleEigengenes(ROD_Susceptible_dds_rlog_matrix, ROD_Sus_full_moduleColors)$eigengenes
ROD_Sus_full_MEs = orderMEs(ROD_Sus_full_MEs0)
ROD_Sus_full_moduleTraitCor = cor(ROD_Sus_full_MEs, ROD_Susceptible_coldata_collapsed_binarize , use = "p");
ROD_Sus_full_moduleTraitPvalue = corPvalueStudent(ROD_Sus_full_moduleTraitCor, ROD_Sus_full_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trait
sizeGrWindow(10,6)
# Will display correlations and their p-values
ROD_Sus_full_textMatrix = paste(signif(ROD_Sus_full_moduleTraitCor, 2), "\n(",
                           signif(ROD_Sus_full_moduleTraitPvalue, 1), ")", sep = "");
dim(ROD_Sus_full_textMatrix) = dim(ROD_Sus_full_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = ROD_Sus_full_moduleTraitCor,
               xLabels = names(ROD_Susceptible_coldata_collapsed_binarize),
               yLabels = names(ROD_Sus_full_MEs),
               ySymbols =names(ROD_Sus_full_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = ROD_Sus_full_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
ROD_Sus_full_moduleTraitCor_df <- as.data.frame(ROD_Sus_full_moduleTraitCor)
ROD_Sus_full_moduleTraitCor_df$mod_names <- row.names(ROD_Sus_full_moduleTraitCor_df)
ROD_Sus_full_moduleTraitCor_df <- ROD_Sus_full_moduleTraitCor_df[,c("mod_names","Condition.Late_Susecptible.vs.Early_Susceptible")]
ROD_Sus_full_moduleTraitPvalue_df <- as.data.frame(ROD_Sus_full_moduleTraitPvalue)
ROD_Sus_full_moduleTraitPvalue_df$mod_names <- row.names(ROD_Sus_full_moduleTraitPvalue_df)
ROD_Sus_full_moduleTraitPvalue_df <- ROD_Sus_full_moduleTraitPvalue_df[,c("mod_names","Condition.Late_Susecptible.vs.Early_Susceptible")]
colnames(ROD_Sus_full_moduleTraitPvalue_df)[2] <- "pvalue"
ROD_Sus_full_moduleTraitCor_Pval_df <- join(ROD_Sus_full_moduleTraitCor_df, ROD_Sus_full_moduleTraitPvalue_df, by = "mod_names")

# Significantly correlated modules
ROD_Sus_full_moduleTraitCor_Pval_df[order(ROD_Sus_full_moduleTraitCor_Pval_df$pvalue),]
class(ROD_Sus_full_moduleTraitCor_Pval_df$pvalue) # numeric
ROD_Sus_full_moduleTraitCor_Pval_df_sig <- ROD_Sus_full_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
ROD_Sus_full_moduleTraitCor_Pval_df_sig # 0 are less than 0.05....the turqousie module is 0.051...

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
ROD_Sus_full_moduleTraitCor_Pval_df_sig_list <- ROD_Sus_full_moduleTraitCor_Pval_df_sig$mod_names
ROD_Sus_full_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(ROD_Sus_full_moduleTraitCor_Pval_df_sig_list, "ME")

## Use function to lookup all apop names for each significant module
matrix_common= ROD_Susceptible_dds_rlog_matrix
moduleColors= ROD_Sus_full_moduleColors
lookup =   C_vir_rtracklayer_apop_product_final

lookup_mod_apop <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$ID %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
names(ROD_Sus_full_moduleTraitCor_Pval_df_sig_list_rm) <- "turquoise"
ROD_Sus_full_module_apop <- lapply(ROD_Sus_full_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop)
ROD_Sus_full_module_apop_df <- do.call(rbind,ROD_Sus_full_module_apop)
ROD_Sus_full_module_apop_df$mod_names <- gsub("\\..*","",row.names(ROD_Sus_full_module_apop_df))
ROD_Sus_full_module_apop_df$mod_names <- gsub("^","ME",ROD_Sus_full_module_apop_df$mod_names)
# add module significance
ROD_Sus_full_module_apop_df <- left_join(ROD_Sus_full_module_apop_df,ROD_Sus_full_moduleTraitCor_Pval_df_sig)
ROD_Sus_full_module_apop_df$exp <- "ROD_Sus"

## Gene relationship to trait and important modules: Gene Significance and Module Membership
# Define variable injected 
ROD_Sus_full_injection = as.data.frame(ROD_Susceptible_coldata_collapsed_binarize$Condition.Late_Susecptible.vs.Early_Susceptible);
names(ROD_Sus_full_injection) = "challenge"
# names (colors) of the modules
ROD_Sus_full_modNames = substring(names(ROD_Sus_full_MEs), 3)
ROD_Sus_full_geneModuleMembership = as.data.frame(cor(ROD_Susceptible_dds_rlog_matrix, ROD_Sus_full_MEs, use = "p"))
ROD_Sus_full_MMPvalue = as.data.frame(corPvalueStudent(as.matrix(ROD_Sus_full_geneModuleMembership), ROD_Sus_full_nSamples))

names(ROD_Sus_full_geneModuleMembership) = paste("MM", ROD_Sus_full_modNames, sep="")
names(ROD_Sus_full_MMPvalue) = paste("p.MM", ROD_Sus_full_modNames, sep="")
#
ROD_Sus_full_geneTraitSignificance = as.data.frame(cor(ROD_Susceptible_dds_rlog_matrix,ROD_Sus_full_injection, use = "p"))
ROD_Sus_full_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(ROD_Sus_full_geneTraitSignificance), ROD_Sus_full_nSamples))
#
names(ROD_Sus_full_geneTraitSignificance) = paste("GS.", names(ROD_Sus_full_injection), sep="")
names(ROD_Sus_full_GSPvalue) = paste("p.GS.", names(ROD_Sus_full_injection), sep="")
#
### Intramodular analysis: identifying genes with high GS and MM
ROD_Sus_full_module = "turquoise"  # weird bell shape, but 0.66 correlation
#ROD_Sus_full_column = match(ROD_Sus_full_module, ROD_Sus_full_modNames)
ROD_Sus_full_moduleGenes = ROD_Sus_full_moduleColors==ROD_Sus_full_module
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
#verboseScatterplot(abs(ROD_Sus_full_geneModuleMembership [ROD_Sus_full_moduleGenes, ROD_Sus_full_column]),
#                   abs(ROD_Sus_full_geneTraitSignificance[ROD_Sus_full_moduleGenes, 1]),
#                   xlab = paste("Module Membership in", ROD_Sus_full_module, "module"),
#                   ylab = "Gene significance for challenge",
#                   main = paste("Module membership vs. gene significance\n"),
#                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = ROD_Sus_full_module) 
#
### IDENTIFY HUB GENES IN EACH SIG MODULE ##
##ROD_Sus_full_colorh = c("MEturquoise")
##
##ROD_Sus_full_Module_hub_genes <- chooseTopHubInEachModule(
##  ROD_Susceptible_dds_rlog_matrix, 
##  ROD_Sus_full_colorh, 
##  power = 9,  # power used for the adjacency network
##  type = "signed hybrid", 
##  corFnc = "bicor"
##)
##class(ROD_Sus_full_Module_hub_genes)
##ROD_Sus_full_Module_hub_genes_df <- as.data.frame(ROD_Sus_full_Module_hub_genes)
##colnames(ROD_Sus_full_Module_hub_genes_df)[1] <- "ID"
##ROD_Sus_full_Module_hub_genes <- merge(ROD_Sus_full_Module_hub_genes_df, C_vir_rtracklayer, by = "ID")
##nrow(ROD_Sus_Module_hub_genes) 
#
#
### Pro_RE22 Pro ###
# define numbers of genes and sample
Pro_RE22_Pro_full_nGenes = ncol(Pro_RE22_dds_rlog_matrix_Pro)
Pro_RE22_Pro_full_nSamples = nrow(Pro_RE22_dds_rlog_matrix_Pro)

# Recalculate MEs with color labels
Pro_RE22_Pro_full_MEs0 = moduleEigengenes(Pro_RE22_dds_rlog_matrix_Pro, Pro_RE22_Pro_full_moduleColors)$eigengenes
Pro_RE22_Pro_full_MEs = orderMEs(Pro_RE22_Pro_full_MEs0)
Pro_RE22_Pro_full_moduleTraitCor = cor(Pro_RE22_Pro_full_MEs, Pro_RE22_coldata_Pro_collapsed_binarize , use = "p");
Pro_RE22_Pro_full_moduleTraitPvalue = corPvalueStudent(Pro_RE22_Pro_full_moduleTraitCor, Pro_RE22_Pro_full_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trait
sizeGrWindow(10,6)
# Will display correlations and their p-values
Pro_RE22_Pro_full_textMatrix = paste(signif(Pro_RE22_Pro_full_moduleTraitCor, 2), "\n(",
                                signif(Pro_RE22_Pro_full_moduleTraitPvalue, 1), ")", sep = "");
dim(Pro_RE22_Pro_full_textMatrix) = dim(Pro_RE22_Pro_full_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = Pro_RE22_Pro_full_moduleTraitCor,
               xLabels = names(Pro_RE22_coldata_Pro_collapsed_binarize),
               yLabels = names(Pro_RE22_Pro_full_MEs),
               ySymbols =names(Pro_RE22_Pro_full_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = Pro_RE22_Pro_full_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with probiotic treatment (high correlation and low P value)?
Pro_RE22_Pro_full_moduleTraitCor_df <- as.data.frame(Pro_RE22_Pro_full_moduleTraitCor)
Pro_RE22_Pro_full_moduleTraitCor_df$mod_names <- row.names(Pro_RE22_Pro_full_moduleTraitCor_df)
Pro_RE22_Pro_full_moduleTraitCor_df <- Pro_RE22_Pro_full_moduleTraitCor_df[,c("mod_names","Condition.Control_no_treatment.vs.Bacillus_pumilus")]
Pro_RE22_Pro_full_moduleTraitPvalue_df <- as.data.frame(Pro_RE22_Pro_full_moduleTraitPvalue)
Pro_RE22_Pro_full_moduleTraitPvalue_df$mod_names <- row.names(Pro_RE22_Pro_full_moduleTraitPvalue_df)
Pro_RE22_Pro_full_moduleTraitPvalue_df <- Pro_RE22_Pro_full_moduleTraitPvalue_df[,c("mod_names","Condition.Control_no_treatment.vs.Bacillus_pumilus")]
colnames(Pro_RE22_Pro_full_moduleTraitPvalue_df)[2] <- "pvalue"
Pro_RE22_Pro_full_moduleTraitCor_Pval_df <- join(Pro_RE22_Pro_full_moduleTraitCor_df, Pro_RE22_Pro_full_moduleTraitPvalue_df, by = "mod_names")

Pro_RE22_Pro_RI_full_moduleTraitCor_df <- as.data.frame(Pro_RE22_Pro_full_moduleTraitCor)
Pro_RE22_Pro_RI_full_moduleTraitCor_df$mod_names <- row.names(Pro_RE22_Pro_RI_full_moduleTraitCor_df)
Pro_RE22_Pro_RI_full_moduleTraitCor_df <- Pro_RE22_Pro_RI_full_moduleTraitCor_df[,c("mod_names","Condition.Phaeobacter_inhibens.vs.Control_no_treatment")]
Pro_RE22_Pro_RI_full_moduleTraitPvalue_df <- as.data.frame(Pro_RE22_Pro_full_moduleTraitPvalue)
Pro_RE22_Pro_RI_full_moduleTraitPvalue_df$mod_names <- row.names(Pro_RE22_Pro_RI_full_moduleTraitPvalue_df)
Pro_RE22_Pro_RI_full_moduleTraitPvalue_df <- Pro_RE22_Pro_RI_full_moduleTraitPvalue_df[,c("mod_names","Condition.Phaeobacter_inhibens.vs.Control_no_treatment")]
colnames(Pro_RE22_Pro_RI_full_moduleTraitPvalue_df)[2] <- "pvalue"
Pro_RE22_Pro_RI_full_moduleTraitCor_Pval_df <- join(Pro_RE22_Pro_RI_full_moduleTraitCor_df, Pro_RE22_Pro_RI_full_moduleTraitPvalue_df, by = "mod_names")

# Significantly correlated modules
Pro_RE22_Pro_full_moduleTraitCor_Pval_df[order(Pro_RE22_Pro_full_moduleTraitCor_Pval_df$pvalue),]
class(Pro_RE22_Pro_full_moduleTraitCor_Pval_df$pvalue) # numeric
Pro_RE22_Pro_full_moduleTraitCor_Pval_df_sig <- Pro_RE22_Pro_full_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
Pro_RE22_Pro_full_moduleTraitCor_Pval_df_sig # 27 (0.051 adds 1 module)

# Significantly correlated modules
Pro_RE22_Pro_RI_full_moduleTraitCor_Pval_df[order(Pro_RE22_Pro_RI_full_moduleTraitCor_Pval_df$pvalue),]
class(Pro_RE22_Pro_RI_full_moduleTraitCor_Pval_df$pvalue) # numeric
Pro_RE22_Pro_RI_full_moduleTraitCor_Pval_df_sig <- Pro_RE22_Pro_RI_full_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
Pro_RE22_Pro_RI_full_moduleTraitCor_Pval_df_sig # 35 (three modules added if I make the threshold 0.51 rather than 0.05)

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
Pro_RE22_Pro_full_moduleTraitCor_Pval_df_sig_list <- Pro_RE22_Pro_full_moduleTraitCor_Pval_df_sig$mod_names
Pro_RE22_Pro_full_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Pro_RE22_Pro_full_moduleTraitCor_Pval_df_sig_list, "ME")

Pro_RE22_Pro_RI_full_moduleTraitCor_Pval_df_sig_list <- Pro_RE22_Pro_RI_full_moduleTraitCor_Pval_df_sig$mod_names
Pro_RE22_Pro_RI_full_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Pro_RE22_Pro_RI_full_moduleTraitCor_Pval_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
names(Pro_RE22_Pro_full_moduleTraitCor_Pval_df_sig_list_rm ) <- c( "floralwhite"    , "darkolivegreen1", "palevioletred3" , "royalblue"  ,     "magenta4"        ,"greenyellow" ,    "plum1"      ,     "chocolate4"  ,   
                                                                   "white"          , "lightslateblue" , "skyblue"        , "thistle"    ,     "sienna4"         ,"skyblue2"    ,    "turquoise"  ,     "lightpink4"  ,   
                                                                    "sienna3"       ,  "darkorange2"   ,  "darkred"       ,  "yellow2"   ,      "darkolivegreen" , "lightpink2" ,     "tan"       ,      "pink3"      ,    
                                                                    "steelblue"     ,  "lightpink3"    ,  "maroon"          )
names(Pro_RE22_Pro_RI_full_moduleTraitCor_Pval_df_sig_list_rm ) <- c("lightgreen"  ,    "darkolivegreen1", "palevioletred2" , "palevioletred3",  "mediumpurple" ,  "red"          ,   "green4"        ,  "lightblue3"     ,
                                                                    "skyblue3"     ,   "darkgrey"        ,"orangered"       ,"magenta4"       , "thistle4"      , "yellowgreen"   ,  "plum"           , "lightslateblue" ,
                                                                     "skyblue"     ,    "thistle"        , "coral"          , "turquoise"     ,  "powderblue"   ,  "darkred"      ,   "darkolivegreen",  "darkolivegreen2",
                                                                     "blue"        ,    "magenta"        , "indianred4"     , "darkslateblue" ,  "tan"          ,  "pink3"        ,   "steelblue"     ,  "lightyellow"    ,
                                                                     "lightpink3"  ,    "thistle3"       , "lavenderblush3"  )
matrix_common= Pro_RE22_dds_rlog_matrix_Pro
moduleColors= Pro_RE22_Pro_full_moduleColors
lookup =   C_vir_rtracklayer_apop_product_final

lookup_mod_apop <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$ID %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
Pro_RE22_Pro_full_module_apop <- lapply(Pro_RE22_Pro_full_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop)
Pro_RE22_Pro_full_module_apop_df <- do.call(rbind,Pro_RE22_Pro_full_module_apop)
Pro_RE22_Pro_full_module_apop_df$mod_names <- gsub("\\..*","",row.names(Pro_RE22_Pro_full_module_apop_df))
Pro_RE22_Pro_full_module_apop_df$mod_names <- gsub("^","ME",Pro_RE22_Pro_full_module_apop_df$mod_names)
# add module significance
Pro_RE22_Pro_full_module_apop_df <- left_join(Pro_RE22_Pro_full_module_apop_df,Pro_RE22_Pro_full_moduleTraitCor_Pval_df_sig)
Pro_RE22_Pro_full_module_apop_df$exp <- "Pro_RE22_Pro_S4"

Pro_RE22_Pro_RI_full_module_apop <- lapply(Pro_RE22_Pro_RI_full_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop)
Pro_RE22_Pro_RI_full_module_apop_df <- do.call(rbind,Pro_RE22_Pro_RI_full_module_apop)
Pro_RE22_Pro_RI_full_module_apop_df$mod_names <- gsub("\\..*","",row.names(Pro_RE22_Pro_RI_full_module_apop_df))
Pro_RE22_Pro_RI_full_module_apop_df$mod_names <- gsub("^","ME",Pro_RE22_Pro_RI_full_module_apop_df$mod_names)
# add module significance
Pro_RE22_Pro_RI_full_module_apop_df <- left_join(Pro_RE22_Pro_RI_full_module_apop_df,Pro_RE22_Pro_RI_full_moduleTraitCor_Pval_df_sig)
Pro_RE22_Pro_RI_full_module_apop_df$exp <- "Pro_RE22_Pro_RI"

## Gene relationship to trait and important modules: Gene Significance and Module Membership
# Define variable injected 
# do later if necessary

## Intramodular analysis: identifying genes with high GS and MM
# later if necessary 

## IDENTIFY HUB GENES IN EACH SIG MODULE ##
# do this later if necessary 
Pro_RE22_Pro_colorh = c()

#Pro_RE22_Pro_Module_hub_genes <- chooseTopHubInEachModule(
#  Pro_RE22_dds_rlog_matrix_common_Pro, 
#  Pro_RE22_Pro_colorh, 
#  power = 2,  # power used for the adjacency network
#  type = "signed hybrid", 
#  corFnc = "bicor"
#)
class(Pro_RE22_Pro_Module_hub_genes)
Pro_RE22_Pro_Module_hub_genes_df <- as.data.frame(Pro_RE22_Pro_Module_hub_genes)
colnames(Pro_RE22_Pro_Module_hub_genes_df)[1] <- "ID"
Pro_RE22_Pro_Module_hub_genes <- merge(Pro_RE22_Pro_Module_hub_genes_df, C_vir_rtracklayer, by = "ID")
nrow(Pro_RE22_Pro_Module_hub_genes) 
# no apoptosis hub genes 

Pro_RE22_Pro_RI_colorh = c()

#Pro_RE22_Pro_RI_Module_hub_genes <- chooseTopHubInEachModule(
#  Pro_RE22_dds_rlog_matrix_common_Pro, 
#  Pro_RE22_Pro_RI_colorh, 
#  power = 2,  # power used for the adjacency network
#  type = "signed hybrid", 
#  corFnc = "bicor"
#)
class(Pro_RE22_Pro_RI_Module_hub_genes)
Pro_RE22_Pro_RI_Module_hub_genes_df <- as.data.frame(Pro_RE22_Pro_RI_Module_hub_genes)
colnames(Pro_RE22_Pro_RI_Module_hub_genes_df)[1] <- "ID"
Pro_RE22_Pro_RI_Module_hub_genes <- merge(Pro_RE22_Pro_RI_Module_hub_genes_df, C_vir_rtracklayer, by = "ID")
nrow(Pro_RE22_Pro_RI_Module_hub_genes) 
#baculoviral IAP repeat-containing protein 7-like, p53 and DNA damage-regulated protein 1-like are both hub genes 


#### MASTER C_VIR FULL MODULE-TRAIT-SIG APOP GENES ####
colnames(Dermo_Tol_full_module_apop_df)[5] <- "mod_signif"
colnames(Dermo_Sus_full_module_apop_df)[5] <- "mod_signif"
colnames(ROD_Res_full_module_apop_df)[5] <- "mod_signif"
#colnames(ROD_Sus_full_module_apop_df)[5] <- "mod_signif" # significance was 0.051...think about whether this is valid to include with others at 0.05 level and if I
  # need to include all decided for the sake of consistency to take out ROD_Sus because it was 0.051 p-value
colnames(Probiotic_full_module_apop_df)[5] <- "mod_signif"
colnames(Pro_RE22_Pro_full_module_apop_df )[5] <- "mod_signif"    # use dataframes from the separated networks rather than the combined network, this top one is S4
colnames(Pro_RE22_Pro_RI_full_module_apop_df )[5] <- "mod_signif" # use dataframes from the separated networks rather than the combined network
colnames(Pro_RE22_RE22_full_module_apop_df)[5] <- "mod_signif" # use dataframes from the separated networks rather than the combined network
C_vir_full_all_exp_mod_sig_apop <- rbind(Dermo_Tol_full_module_apop_df,Dermo_Sus_full_module_apop_df,ROD_Res_full_module_apop_df,
                                         #ROD_Sus_full_module_apop_df,
                                    Probiotic_full_module_apop_df,Pro_RE22_Pro_full_module_apop_df ,Pro_RE22_Pro_RI_full_module_apop_df ,Pro_RE22_RE22_full_module_apop_df)
nrow(C_vir_full_all_exp_mod_sig_apop) # 727 (with ROD_Sus removed) Aug 24th, 2020
C_vir_full_all_exp_mod_sig_apop_positive <- C_vir_full_all_exp_mod_sig_apop %>% filter(mod_signif >0)

#### C. GIGAS WGCNA ####

#### DATA FORMATTING, BATCH EFFECT REMOVAL ####

## ZHANG, no batch effect removal  

# save as matrix with assay and transform
Zhang_dds_rlog_matrix <- assay(Zhang_dds_rlog)
class(Zhang_dds_rlog_matrix)
Zhang_dds_rlog_matrix  <- t(Zhang_dds_rlog_matrix) 

## Rubio, no batch effect removal
Rubio_dds_rlog_matrix <- assay(Rubio_dds_rlog)
class(Rubio_dds_rlog_matrix)
Rubio_dds_rlog_matrix <- t(Rubio_dds_rlog_matrix)

## DeLorgeril, no batch effect correction  
deLorgeril_Resistant_dds_vst_matrix   <- assay(deLorgeril_Resistant_dds_vst)
deLorgeril_Susceptible_dds_vst_matrix <-  assay(deLorgeril_Susceptible_dds_vst)
class(deLorgeril_Resistant_dds_vst_matrix  )
class(deLorgeril_Susceptible_dds_vst_matrix)

deLorgeril_Resistant_dds_vst_matrix  <- t(deLorgeril_Resistant_dds_vst_matrix  )
deLorgeril_Susceptible_dds_vst_matrix <- t(deLorgeril_Susceptible_dds_vst_matrix)

## He, no batch effect correction 
He_dds_vst_matrix <- assay(He_dds_vst)
class(He_dds_vst_matrix)
He_dds_vst_matrix <- t(He_dds_vst_matrix)

## LIMIT ANALYSIS TO TRANSCRIPTS EXPRESSED IN ALL EXPERIMENTS
ncol(Zhang_dds_rlog_matrix) # 33418
ncol(Rubio_dds_rlog_matrix)  # 44705
ncol(deLorgeril_Resistant_dds_vst_matrix) # 44415
ncol(deLorgeril_Susceptible_dds_vst_matrix) # 44817
ncol(He_dds_vst_matrix) # 38633

Zhang_dds_rlog_matrix_colnames <- colnames(Zhang_dds_rlog_matrix)
Rubio_dds_rlog_matrix_colnames <- colnames(Rubio_dds_rlog_matrix)
deLorgeril_Resistant_dds_vst_matrix_colnames <- colnames(deLorgeril_Resistant_dds_vst_matrix)
deLorgeril_Susceptible_dds_vst_matrix_colnames <- colnames(deLorgeril_Susceptible_dds_vst_matrix)
He_dds_vst_matrix_colnames <- colnames(He_dds_vst_matrix)

C_gig_common_vst_transcripts <- Reduce(intersect, list(Zhang_dds_rlog_matrix_colnames,
                                                       Rubio_dds_rlog_matrix_colnames,
                                                       deLorgeril_Resistant_dds_vst_matrix_colnames,
                                                       deLorgeril_Susceptible_dds_vst_matrix_colnames,
                                                       He_dds_vst_matrix_colnames))
head(C_gig_common_vst_transcripts)
length(C_gig_common_vst_transcripts) #27678

Zhang_dds_rlog_matrix_common <- Zhang_dds_rlog_matrix[, C_gig_common_vst_transcripts]
Rubio_dds_rlog_matrix_common <- Rubio_dds_rlog_matrix [, C_gig_common_vst_transcripts]
deLorgeril_Resistant_dds_vst_matrix_common <- deLorgeril_Resistant_dds_vst_matrix[, C_gig_common_vst_transcripts]
deLorgeril_Susceptible_dds_vst_matrix_common <- deLorgeril_Susceptible_dds_vst_matrix[, C_gig_common_vst_transcripts]
He_dds_vst_matrix_common <- He_dds_vst_matrix[, C_gig_common_vst_transcripts]
 
ncol(Zhang_dds_rlog_matrix_common) # 27678
ncol(Rubio_dds_rlog_matrix_common) # 27678
ncol(deLorgeril_Resistant_dds_vst_matrix_common) # 27678
ncol(deLorgeril_Susceptible_dds_vst_matrix_common) # 27678
ncol(He_dds_vst_matrix_common) # 27678

# Do column names agree between all?
all(colnames(Zhang_dds_rlog_matrix_common ) %in% colnames(He_dds_vst_matrix_common)) # TRUE
all(colnames(Zhang_dds_rlog_matrix_common) == colnames(He_dds_vst_matrix_common)) # TRUE 

#### SELECT SOFT THRESHOLDING POWER FOR EACH EXPERIMENT ####
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Load picksoftthreshold
C_gigas_picksoft = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gigas_individual_pickSoftThreshold_WGCNA_input.RData")

# Call the network topology analysis function, following general recommendations to set network type to "signed hybrid" and using the "bicor" correlation
#Zhang_dds_rlog_matrix_common_sft <- pickSoftThreshold(Zhang_dds_rlog_matrix_common, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
#Rubio_dds_rlog_matrix_common_sft <- pickSoftThreshold(Rubio_dds_rlog_matrix_common, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
#deLorgeril_Resistant_dds_vst_matrix_common_sft <- pickSoftThreshold(deLorgeril_Resistant_dds_vst_matrix_common, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
#deLorgeril_Susceptible_dds_vst_matrix_common_sft <- pickSoftThreshold(deLorgeril_Susceptible_dds_vst_matrix_common, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
#He_dds_vst_matrix_common_sft <- pickSoftThreshold(He_dds_vst_matrix_common, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
# save soft thresholding results
#save(Zhang_dds_rlog_matrix_common_sft,Rubio_dds_rlog_matrix_common_sft,deLorgeril_Resistant_dds_vst_matrix_common_sft,deLorgeril_Susceptible_dds_vst_matrix_common_sft,
#     He_dds_vst_matrix_common_sft, file= "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gigas_individual_pickSoftThreshold_WGCNA_input.RData")

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Zhang 
# Scale-free topology fit index as a function of the soft-thresholding power
plot(Zhang_dds_rlog_matrix_common_sft$fitIndices[,1], -sign(Zhang_dds_rlog_matrix_common_sft$fitIndices[,3])*Zhang_dds_rlog_matrix_common_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(Zhang_dds_rlog_matrix_common_sft$fitIndices[,1], -sign(Zhang_dds_rlog_matrix_common_sft$fitIndices[,3])*Zhang_dds_rlog_matrix_common_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(Zhang_dds_rlog_matrix_common_sft$fitIndices[,1], Zhang_dds_rlog_matrix_common_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(Zhang_dds_rlog_matrix_common_sft$fitIndices[,1], Zhang_dds_rlog_matrix_common_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Softthreshold of 4  since this is lowest value past 0.9 we start to see flattening 


# Rubio
# Scale-free topology fit index as a function of the soft-thresholding power
plot(Rubio_dds_rlog_matrix_common_sft$fitIndices[,1], -sign(Rubio_dds_rlog_matrix_common_sft$fitIndices[,3])*Rubio_dds_rlog_matrix_common_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(Rubio_dds_rlog_matrix_common_sft$fitIndices[,1], -sign(Rubio_dds_rlog_matrix_common_sft$fitIndices[,3])*Rubio_dds_rlog_matrix_common_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(Rubio_dds_rlog_matrix_common_sft$fitIndices[,1], Rubio_dds_rlog_matrix_common_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(Rubio_dds_rlog_matrix_common_sft$fitIndices[,1], Rubio_dds_rlog_matrix_common_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Softthreshold of 3 since this is lowest value past 0.9 we start to see flattening 

# deLorg Res
# Scale-free topology fit index as a function of the soft-thresholding power
plot(deLorgeril_Resistant_dds_vst_matrix_common_sft$fitIndices[,1], -sign(deLorgeril_Resistant_dds_vst_matrix_common_sft$fitIndices[,3])*deLorgeril_Resistant_dds_vst_matrix_common_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(deLorgeril_Resistant_dds_vst_matrix_common_sft$fitIndices[,1], -sign(deLorgeril_Resistant_dds_vst_matrix_common_sft$fitIndices[,3])*deLorgeril_Resistant_dds_vst_matrix_common_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(deLorgeril_Resistant_dds_vst_matrix_common_sft$fitIndices[,1], deLorgeril_Resistant_dds_vst_matrix_common_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(deLorgeril_Resistant_dds_vst_matrix_common_sft$fitIndices[,1], deLorgeril_Resistant_dds_vst_matrix_common_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Softthreshold of 4 since this is lowest value past 0.9 we start to see flattening 

# deLorg Sus = 21 samples
# Scale-free topology fit index as a function of the soft-thresholding power
plot(deLorgeril_Susceptible_dds_vst_matrix_common_sft$fitIndices[,1], -sign(deLorgeril_Susceptible_dds_vst_matrix_common_sft$fitIndices[,3])*deLorgeril_Susceptible_dds_vst_matrix_common_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(deLorgeril_Susceptible_dds_vst_matrix_common_sft$fitIndices[,1], -sign(deLorgeril_Susceptible_dds_vst_matrix_common_sft$fitIndices[,3])*deLorgeril_Susceptible_dds_vst_matrix_common_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(deLorgeril_Susceptible_dds_vst_matrix_common_sft$fitIndices[,1], deLorgeril_Susceptible_dds_vst_matrix_common_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(deLorgeril_Susceptible_dds_vst_matrix_common_sft$fitIndices[,1], deLorgeril_Susceptible_dds_vst_matrix_common_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Softthreshold of 8 since this is lowest value past 0.9 we start to see flattening and 21 samples

# He
# Scale-free topology fit index as a function of the soft-thresholding power
plot(He_dds_vst_matrix_common_sft$fitIndices[,1], -sign(He_dds_vst_matrix_common_sft$fitIndices[,3])*He_dds_vst_matrix_common_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(He_dds_vst_matrix_common_sft$fitIndices[,1], -sign(He_dds_vst_matrix_common_sft$fitIndices[,3])*He_dds_vst_matrix_common_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(He_dds_vst_matrix_common_sft$fitIndices[,1], He_dds_vst_matrix_common_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(He_dds_vst_matrix_common_sft$fitIndices[,1], He_dds_vst_matrix_common_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Softthreshold of 6 since this is lowest value past 0.9 we start to see flattening 

#### ONE STEP NETWORK CONSTRUCTION, MODULE DETECTION, MODULE DENDROGRAM INSPECTION #### 

# Zhang 
#Zhang_net = blockwiseModules(Zhang_dds_rlog_matrix_common, power = 4, # picked suitable power in the code above 
#                                 TOMType = "signed", # use signed TOM type
#                                 networkType= "signed hybrid", # use signed hybrid network type
#                                 corType = "bicor", # use suggested bicor
#                                 TminModuleSize = 30, # recommended default
#                                 reassignThreshold = 0, # recommended default
#                                 mergeCutHeight = 0.25, # recommended default
#                                 numericLabels = TRUE, # recommended default
#                                 pamRespectsDendro = FALSE,# recommended default
#                                 verbose = 3, 
#                                 maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
# How many modules identified
table(Zhang_net$colors)
# Plot dendrogram with colors
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
Zhang_mergedColors = labels2colors(Zhang_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(Zhang_net$dendrograms[[1]], Zhang_mergedColors[Zhang_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
Zhang_moduleLabels = Zhang_net$colors
Zhang_moduleColors = labels2colors(Zhang_net$colors)
Zhang_MEs = Zhang_net$MEs
Zhang_geneTree = Zhang_net$dendrograms[[1]]
#
## Save Zhang network
#save(Zhang_net, Zhang_mergedColors, Zhang_moduleLabels, Zhang_moduleColors, Zhang_MEs,Zhang_geneTree, 
#     file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Zhang_network.RData")

# Load Zhang network
#Zhang_network = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Zhang_network.RData")


# Rubio
# Rubio_net = blockwiseModules(Rubio_dds_rlog_matrix_common, power = 3, # picked suitable power in the code above 
#                              TOMType = "signed", # use signed TOM type
#                              networkType= "signed hybrid", # use signed hybrid network type
#                              corType = "bicor", # use suggested bicor
#                              TminModuleSize = 30, # recommended default
#                              reassignThreshold = 0, # recommended default
#                              mergeCutHeight = 0.25, # recommended default
#                              numericLabels = TRUE, # recommended default
#                              pamRespectsDendro = FALSE,# recommended default
#                              verbose = 3, 
#                              maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
# # How many modules identified
# table(Rubio_net$colors)
# # Plot dendrogram with colors
# # open a graphics window
# sizeGrWindow(12, 9)
# # Convert labels to colors for plotting
# Rubio_mergedColors = labels2colors(Rubio_net$colors)
# # Plot the dendrogram and the module colors underneath
# plotDendroAndColors(Rubio_net$dendrograms[[1]], Rubio_mergedColors[Rubio_net$blockGenes[[1]]],
#                     "Module colors",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# Rubio_moduleLabels = Rubio_net$colors
# Rubio_moduleColors = labels2colors(Rubio_net$colors)
# Rubio_MEs = Rubio_net$MEs
# Rubio_geneTree = Rubio_net$dendrograms[[1]]

# save network
#save(Rubio_net, Rubio_mergedColors, Rubio_moduleLabels, Rubio_moduleColors, Rubio_MEs,Rubio_geneTree, 
#         file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Rubio_network.RData")

# Load Rubio network
Rubio_network = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Rubio_network.RData")


# deLorg Res
#deLorg_Res_net = blockwiseModules(deLorgeril_Resistant_dds_vst_matrix_common, power = 4, # picked suitable power in the code above 
#                             TOMType = "signed", # use signed TOM type
#                             networkType= "signed hybrid", # use signed hybrid network type
#                             corType = "bicor", # use suggested bicor
#                             TminModuleSize = 30, # recommended default
#                             reassignThreshold = 0, # recommended default
#                             mergeCutHeight = 0.25, # recommended default
#                             numericLabels = TRUE, # recommended default
#                             pamRespectsDendro = FALSE,# recommended default
#                             verbose = 3, 
#                             maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
## How many modules identified
#table(deLorg_Res_net$colors)
## Plot dendrogram with colors
## open a graphics window
#sizeGrWindow(12, 9)
## Convert labels to colors for plotting
#deLorg_Res_mergedColors = labels2colors(deLorg_Res_net$colors)
## Plot the dendrogram and the module colors underneath
#plotDendroAndColors(deLorg_Res_net$dendrograms[[1]], deLorg_Res_mergedColors[deLorg_Res_net$blockGenes[[1]]],
#                    "Module colors",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05)
#deLorg_Res_moduleLabels = deLorg_Res_net$colors
#deLorg_Res_moduleColors = labels2colors(deLorg_Res_net$colors)
#deLorg_Res_MEs = deLorg_Res_net$MEs
#deLorg_Res_geneTree = deLorg_Res_net$dendrograms[[1]]
#
#save(deLorg_Res_net, deLorg_Res_mergedColors,deLorg_Res_moduleLabels, deLorg_Res_moduleColors, deLorg_Res_MEs,deLorg_Res_geneTree, 
#     file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/deLorg_Res_network.RData")

# load delorg_Res network
delorg_Res_network = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/deLorg_Res_network.RData")

# deLorg Sus
#deLorg_Sus_net = blockwiseModules(deLorgeril_Susceptible_dds_vst_matrix_common, power = 8, # picked suitable power in the code above 
#                                  TOMType = "signed", # use signed TOM type
#                                  networkType= "signed hybrid", # use signed hybrid network type
#                                  corType = "bicor", # use suggested bicor
#                                  TminModuleSize = 30, # recommended default
#                                  reassignThreshold = 0, # recommended default
#                                  mergeCutHeight = 0.25, # recommended default
#                                  numericLabels = TRUE, # recommended default
#                                  pamRespectsDendro = FALSE,# recommended default
#                                  verbose = 3, 
#                                  maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
## How many modules identified
#table(deLorg_Sus_net$colors)
## Plot dendrogram with colors
## open a graphics window
#sizeGrWindow(12, 9)
## Convert labels to colors for plotting
#deLorg_Sus_mergedColors = labels2colors(deLorg_Sus_net$colors)
## Plot the dendrogram and the module colors underneath
#plotDendroAndColors(deLorg_Sus_net$dendrograms[[1]], deLorg_Sus_mergedColors[deLorg_Sus_net$blockGenes[[1]]],
#                    "Module colors",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05)
#deLorg_Sus_moduleLabels = deLorg_Sus_net$colors
#deLorg_Sus_moduleColors = labels2colors(deLorg_Sus_net$colors)
#deLorg_Sus_MEs = deLorg_Sus_net$MEs
#deLorg_Sus_geneTree = deLorg_Sus_net$dendrograms[[1]]
#
#save(deLorg_Sus_net, deLorg_Sus_mergedColors,deLorg_Sus_moduleLabels, deLorg_Sus_moduleColors, deLorg_Sus_MEs,deLorg_Sus_geneTree, 
#     file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/deLorg_Sus_network.RData")

# load delorg_Res network
delorg_Sus_network = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/deLorg_Sus_network.RData")

# He 
#He_net = blockwiseModules(He_dds_vst_matrix_common, power = 6, # picked suitable power in the code above 
#                                  TOMType = "signed", # use signed TOM type
#                                  networkType= "signed hybrid", # use signed hybrid network type
#                                  corType = "bicor", # use suggested bicor
#                                  TminModuleSize = 30, # recommended default
#                                  reassignThreshold = 0, # recommended default
#                                  mergeCutHeight = 0.25, # recommended default
#                                  numericLabels = TRUE, # recommended default
#                                  pamRespectsDendro = FALSE,# recommended default
#                                  verbose = 3, 
#                                  maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
## How many modules identified
#table(He_net$colors)
## Plot dendrogram with colors
## open a graphics window
#sizeGrWindow(12, 9)
## Convert labels to colors for plotting
#He_mergedColors = labels2colors(He_net$colors)
## Plot the dendrogram and the module colors underneath
#plotDendroAndColors(He_net$dendrograms[[1]], He_mergedColors[He_net$blockGenes[[1]]],
#                    "Module colors",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05)
#He_moduleLabels = He_net$colors
#He_moduleColors = labels2colors(He_net$colors)
#He_MEs = He_net$MEs
#He_geneTree = He_net$dendrograms[[1]]
#
#save(He_net, He_mergedColors,He_moduleLabels, He_moduleColors, He_MEs,He_geneTree, 
#     file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/He_network.RData")

# load He network
He_network = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/He_network.RData")

#### BINARIZE ALL CATEGORICAL VARIABLES TO TEST DISEASE CHALLENGE ASSOCIATIONS ####
# out = binarizeCategoricalVariable(x,
# includePairwise = TRUE,
# includeLevelVsAll = FALSE);

## Zhang
Zhang_dds_rlog_matrix_common
Zhang_coldata_collapsed <- Zhang_coldata[,c("condition","group_by_sim")]
Zhang_coldata_collapsed$group_by_sim <- str_replace(Zhang_coldata_collapsed$group_by_sim, "V_aes_V_alg1_V_alg2", "Vibrio")
Zhang_coldata_collapsed$group_by_sim <- str_replace(Zhang_coldata_collapsed$group_by_sim, "V_tub_V_ang", "Vibrio")

all(row.names(Zhang_coldata_collapsed) %in% row.names(Zhang_dds_rlog_matrix_common)) # TRUE
all(row.names(Zhang_dds_rlog_matrix_common) %in% row.names(Zhang_coldata_collapsed) ) # TRUE
all(row.names(Zhang_coldata_collapsed) == row.names(Zhang_dds_rlog_matrix_common)) # TRUE
all(row.names(Zhang_dds_rlog_matrix_common) == row.names(Zhang_coldata_collapsed) ) # TRUE

# binarize
Zhang_coldata_collapsed_binarize <- binarizeCategoricalColumns.pairwise(Zhang_coldata_collapsed)
row.names(Zhang_coldata_collapsed_binarize) <- row.names(Zhang_coldata_collapsed)
colnames(Zhang_coldata_collapsed_binarize)
Zhang_coldata_collapsed_binarize <- Zhang_coldata_collapsed_binarize[,c(37:39)]

## Rubio
Rubio_dds_rlog_matrix_common
Rubio_coldata_collapsed <- as.data.frame(Rubio_coldata[,2])
row.names(Rubio_coldata_collapsed) <- row.names(Rubio_coldata)
colnames(Rubio_coldata_collapsed)[1] <- "Group"
nrow(Rubio_coldata_collapsed) # 18
nrow(Rubio_dds_rlog_matrix_common) # 18

all(row.names(Rubio_coldata_collapsed) %in% row.names(Rubio_dds_rlog_matrix_common)) # TRUE
all(row.names(Rubio_dds_rlog_matrix_common) %in% row.names(Rubio_coldata_collapsed) ) # TRUE
all(row.names(Rubio_coldata_collapsed) == row.names(Rubio_dds_rlog_matrix_common)) # TRUE
all(row.names(Rubio_dds_rlog_matrix_common) == row.names(Rubio_coldata_collapsed) ) # TRUE

# binarize
Rubio_coldata_collapsed_binarize <- binarizeCategoricalColumns.pairwise(Rubio_coldata_collapsed)
row.names(Rubio_coldata_collapsed_binarize) <- row.names(Rubio_coldata_collapsed)
colnames(Rubio_coldata_collapsed_binarize)

## deLorg Res

deLorgeril_Resistant_coldata_collapsed <- as.data.frame(deLorgeril_Resistant_coldata[,c(1)])
row.names(deLorgeril_Resistant_coldata_collapsed) <- row.names(deLorgeril_Resistant_coldata)
colnames(deLorgeril_Resistant_coldata_collapsed) <- "Group"
nrow(deLorgeril_Resistant_coldata_collapsed) # 21
nrow(deLorgeril_Resistant_dds_vst_matrix_common) # 21

all(row.names(deLorgeril_Resistant_coldata_collapsed) %in% row.names(deLorgeril_Resistant_dds_vst_matrix_common)) # TRUE
all(row.names(deLorgeril_Resistant_dds_vst_matrix_common) %in% row.names(deLorgeril_Resistant_coldata_collapsed) ) # TRUE
all(row.names(deLorgeril_Resistant_coldata_collapsed) == row.names(deLorgeril_Resistant_dds_vst_matrix_common)) # TRUE
all(row.names(deLorgeril_Resistant_dds_vst_matrix_common) == row.names(deLorgeril_Resistant_coldata_collapsed) ) # TRUE

# binarize
deLorgeril_Resistant_coldata_collapsed_binarize <- binarizeCategoricalColumns.pairwise(deLorgeril_Resistant_coldata_collapsed)
row.names(deLorgeril_Resistant_coldata_collapsed_binarize) <- row.names(deLorgeril_Resistant_coldata_collapsed)
colnames(deLorgeril_Resistant_coldata_collapsed_binarize)

## deLorg Sus
deLorgeril_Susceptible_coldata_collapsed <- as.data.frame(deLorgeril_Susceptible_coldata[,c(1)])
row.names(deLorgeril_Susceptible_coldata_collapsed) <- row.names(deLorgeril_Susceptible_coldata)
colnames(deLorgeril_Susceptible_coldata_collapsed) <- "Group"
nrow(deLorgeril_Susceptible_coldata_collapsed) # 21
nrow(deLorgeril_Susceptible_dds_vst_matrix_common) # 21

all(row.names(deLorgeril_Susceptible_coldata_collapsed) %in% row.names(deLorgeril_Susceptible_dds_vst_matrix_common)) # TRUE
all(row.names(deLorgeril_Susceptible_dds_vst_matrix_common) %in% row.names(deLorgeril_Susceptible_coldata_collapsed) ) # TRUE
all(row.names(deLorgeril_Susceptible_coldata_collapsed) == row.names(deLorgeril_Susceptible_dds_vst_matrix_common)) # TRUE
all(row.names(deLorgeril_Susceptible_dds_vst_matrix_common) == row.names(deLorgeril_Susceptible_coldata_collapsed) ) # TRUE

# binarize
deLorgeril_Susceptible_coldata_collapsed_binarize <- binarizeCategoricalColumns.pairwise(deLorgeril_Susceptible_coldata_collapsed)
row.names(deLorgeril_Susceptible_coldata_collapsed_binarize) <- row.names(deLorgeril_Susceptible_coldata_collapsed)
colnames(deLorgeril_Susceptible_coldata_collapsed_binarize)

## He
He_coldata_collapsed <- as.data.frame(He_coldata[,c(1)])
row.names(He_coldata_collapsed) <- row.names(He_coldata)
colnames(He_coldata_collapsed) <- "Group"
nrow(He_coldata_collapsed) # 32
nrow(He_dds_vst_matrix_common) # 32

all(row.names(He_coldata_collapsed) %in% row.names(He_dds_vst_matrix_common)) # TRUE
all(row.names(He_dds_vst_matrix_common) %in% row.names(He_coldata_collapsed) ) # TRUE
all(row.names(He_coldata_collapsed) == row.names(He_dds_vst_matrix_common)) # TRUE
all(row.names(He_dds_vst_matrix_common) == row.names(He_coldata_collapsed) ) # TRUE

# binarize
He_coldata_collapsed_binarize <- binarizeCategoricalColumns.pairwise(He_coldata_collapsed)
row.names(He_coldata_collapsed_binarize) <- row.names(He_coldata_collapsed)
colnames(He_coldata_collapsed_binarize)

#### QUANTIFYING MODULE TRAIT ASSOCIATIONS ####

#### Zhang ####
# Define numbers of genes and samples
Zhang_nGenes = ncol(Zhang_dds_rlog_matrix_common)
Zhang_nSamples = nrow(Zhang_dds_rlog_matrix_common)

# Recalculate MEs with color labels
Zhang_MEs0 = moduleEigengenes(Zhang_dds_rlog_matrix_common, Zhang_moduleColors)$eigengenes
Zhang_MEs = orderMEs(Zhang_MEs0)
Zhang_moduleTraitCor = cor(Zhang_MEs, Zhang_coldata_collapsed_binarize, use = "p");
Zhang_moduleTraitPvalue = corPvalueStudent(Zhang_moduleTraitCor, Zhang_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trai
sizeGrWindow(10,6)
# Will display correlations and their p-values
Zhang_textMatrix = paste(signif(Zhang_moduleTraitCor, 2), "\n(",
                             signif(Zhang_moduleTraitPvalue, 1), ")", sep = "");
dim(Zhang_textMatrix) = dim(Zhang_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = Zhang_moduleTraitCor,
               xLabels = names(Zhang_coldata_collapsed_binarize),
               yLabels = names(Zhang_MEs),
               ySymbols = names(Zhang_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = Zhang_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
Zhang_moduleTraitCor_df <- as.data.frame(Zhang_moduleTraitCor)
Zhang_moduleTraitCor_df$mod_names <- row.names(Zhang_moduleTraitCor_df)
Zhang_moduleTraitCor_LPS_df <- Zhang_moduleTraitCor_df[,c("mod_names","group_by_sim.LPS_M_lut.vs.control")]
Zhang_moduleTraitPvalue_df <- as.data.frame(Zhang_moduleTraitPvalue)
Zhang_moduleTraitPvalue_df$mod_names <- row.names(Zhang_moduleTraitPvalue_df )
Zhang_moduleTraitPvalue_LPS_df <- Zhang_moduleTraitPvalue_df[,c("mod_names","group_by_sim.LPS_M_lut.vs.control")]
colnames(Zhang_moduleTraitPvalue_LPS_df)[2] <- "pvalue"
Zhang_moduleTraitCor_Pval_LPS_df <- join(Zhang_moduleTraitCor_LPS_df, Zhang_moduleTraitPvalue_LPS_df, by = "mod_names")

# Significantly correlated modules
Zhang_moduleTraitCor_Pval_LPS_df[order(Zhang_moduleTraitCor_Pval_LPS_df$pvalue),]
class(Zhang_moduleTraitCor_Pval_LPS_df$pvalue) # numeric
Zhang_moduleTraitCor_Pval_LPS_df_sig <- Zhang_moduleTraitCor_Pval_LPS_df %>% filter(pvalue <= 0.05)
Zhang_moduleTraitCor_Pval_LPS_df_sig # 47

Zhang_moduleTraitCor_df <- as.data.frame(Zhang_moduleTraitCor)
Zhang_moduleTraitCor_df$mod_names <- row.names(Zhang_moduleTraitCor_df)
Zhang_moduleTraitCor_Vibrio_df <- Zhang_moduleTraitCor_df[,c("mod_names","group_by_sim.Vibrio.vs.control")]
Zhang_moduleTraitPvalue_Vibrio_df <- Zhang_moduleTraitPvalue_df[,c("mod_names","group_by_sim.Vibrio.vs.control")]
colnames(Zhang_moduleTraitPvalue_Vibrio_df)[2] <- "pvalue"
Zhang_moduleTraitCor_Pval_Vibrio_df <- join(Zhang_moduleTraitCor_Vibrio_df, Zhang_moduleTraitPvalue_Vibrio_df, by = "mod_names")

# Significantly correlated modules
Zhang_moduleTraitCor_Pval_Vibrio_df[order(Zhang_moduleTraitCor_Pval_Vibrio_df$pvalue),]
class(Zhang_moduleTraitCor_Pval_Vibrio_df$pvalue) # numeric
Zhang_moduleTraitCor_Pval_Vibrio_df_sig <- Zhang_moduleTraitCor_Pval_Vibrio_df %>% filter(pvalue <= 0.05)
Zhang_moduleTraitCor_Pval_Vibrio_df_sig # 16

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
Zhang_moduleTraitCor_Pval_LPS_df_sig_list <- Zhang_moduleTraitCor_Pval_LPS_df_sig$mod_names
Zhang_moduleTraitCor_Pval_LPS_df_sig_list_rm <- str_remove(Zhang_moduleTraitCor_Pval_LPS_df_sig_list, "ME")
Zhang_moduleTraitCor_Pval_Vibrio_df_sig_list <- Zhang_moduleTraitCor_Pval_Vibrio_df_sig$mod_names
Zhang_moduleTraitCor_Pval_Vibrio_df_sig_list_rm <- str_remove(Zhang_moduleTraitCor_Pval_Vibrio_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
names(Zhang_moduleTraitCor_Pval_Vibrio_df_sig_list_rm) <- c("blue"      ,"blue4"  ,"honeydew1"  ,"brown1"         , "yellow2"   ,"pink4"     ,      "salmon2","greenyellow" ,    "thistle1",       
                                                            "darkgreen" ,"purple" ,"lightblue4" ,"darkolivegreen2", "lightgreen","orangered" ,      "lightskyblue4"   )
matrix_common= Zhang_dds_rlog_matrix_common
moduleColors= Zhang_moduleColors
lookup =   C_gig_rtracklayer_apop_product_final

lookup_mod_apop_CG <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$transcript_id %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
Zhang_module_apop <- lapply(Zhang_moduleTraitCor_Pval_Vibrio_df_sig_list_rm,  lookup_mod_apop_CG)
Zhang_module_apop_df <- do.call(rbind,Zhang_module_apop)
Zhang_module_apop_df$mod_names <- gsub("\\..*","",row.names(Zhang_module_apop_df))
Zhang_module_apop_df$mod_names <- gsub("^","ME",Zhang_module_apop_df$mod_names)
# add module significance
Zhang_module_Vibrio_apop_df <- left_join(Zhang_module_apop_df,Zhang_moduleTraitCor_Pval_Vibrio_df)
Zhang_module_Vibrio_apop_df$exp <- "Zhang_Vibrio"

# specify names for list of lists
names(Zhang_moduleTraitCor_Pval_LPS_df_sig_list_rm) <- c("darkolivegreen4","chocolate4"      ,"darkorange2","honeydew"       ,"pink"          ,"sienna4"       ,"darkseagreen2"   ,"magenta4"      ,"grey60"     ,"tan3"           ,
                                                          "ivory"         , "black"          , "powderblue", "darkorange"    , "firebrick4"   , "lightpink2"   , "lightsteelblue" , "midnightblue" , "blue"      , "darkseagreen4"  ,
                                                          "navajowhite"   , "magenta"        , "lightpink3", "mediumpurple2" , "brown1"       , "yellow2"      , "pink4"          , "paleturquoise", "turquoise" , "thistle1"       ,
                                                          "cyan"          , "darkolivegreen2", "green4"    , "palevioletred2", "darkseagreen3", "mediumpurple3", "orangered3"     , "tan"          , "indianred3", "mistyrose"      ,
                                                          "orangered"     , "coral2"         , "plum3"     , "coral4"        , "lightskyblue4", "blue3"        , "navajowhite2"  )
Zhang_LPS_module_apop <- lapply(Zhang_moduleTraitCor_Pval_LPS_df_sig_list_rm,  lookup_mod_apop_CG)
Zhang_LPS_module_apop_df <- do.call(rbind,Zhang_LPS_module_apop)
Zhang_LPS_module_apop_df$mod_names <- gsub("\\..*","",row.names(Zhang_LPS_module_apop_df))
Zhang_LPS_module_apop_df$mod_names <- gsub("^","ME",Zhang_LPS_module_apop_df$mod_names)
# add module significance
Zhang_LPS_module_apop_df <- left_join(Zhang_LPS_module_apop_df,Zhang_moduleTraitCor_Pval_LPS_df_sig)
Zhang_LPS_module_apop_df$exp <- "Zhang_LPS"

## Gene relationship to trait and important modules: Gene Significance and Module Membership
# do later

## Intramodular analysis: identifying genes with high GS and MM
# Using the GS and MM measures, we can identify genes that have a high significance for disease challenge
# as well as high module membership in interesting modules. As an example, we look at the brown module 
# that has the highest association with weight. We plot a scatterplot of Gene Significance vs. Module Membership in the brown module:

## IDENTIFY HUB GENES IN EACH SIG MODULE ##

#### Rubio ####
# Define numbers of genes and samples
Rubio_nGenes = ncol(Rubio_dds_rlog_matrix_common)
Rubio_nSamples = nrow(Rubio_dds_rlog_matrix_common)

# Recalculate MEs with color labels
Rubio_MEs0 = moduleEigengenes(Rubio_dds_rlog_matrix_common, Rubio_moduleColors)$eigengenes
Rubio_MEs = orderMEs(Rubio_MEs0)
Rubio_moduleTraitCor = cor(Rubio_MEs, Rubio_coldata_collapsed_binarize, use = "p");
Rubio_moduleTraitPvalue = corPvalueStudent(Rubio_moduleTraitCor, Rubio_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trai
sizeGrWindow(10,6)
# Will display correlations and their p-values
Rubio_textMatrix = paste(signif(Rubio_moduleTraitCor, 2), "\n(",
                         signif(Rubio_moduleTraitPvalue, 1), ")", sep = "");
dim(Rubio_textMatrix) = dim(Rubio_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = Rubio_moduleTraitCor,
               xLabels = names(Rubio_coldata_collapsed_binarize),
               yLabels = names(Rubio_MEs),
               ySymbols = names(Rubio_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = Rubio_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
Rubio_moduleTraitCor_df <- as.data.frame(Rubio_moduleTraitCor)
Rubio_moduleTraitCor_df$mod_names <- row.names(Rubio_moduleTraitCor_df)
Rubio_moduleTraitCor_NV_df <- Rubio_moduleTraitCor_df[,c("mod_names","Group.Non_virulent.vs.Control")]
Rubio_moduleTraitPvalue_df <- as.data.frame(Rubio_moduleTraitPvalue)
Rubio_moduleTraitPvalue_df$mod_names <- row.names(Rubio_moduleTraitPvalue_df )
Rubio_moduleTraitPvalue_NV_df <- Rubio_moduleTraitPvalue_df[,c("mod_names","Group.Non_virulent.vs.Control")]
colnames(Rubio_moduleTraitPvalue_NV_df)[2] <- "pvalue"
Rubio_moduleTraitCor_Pval_NV_df <- join(Rubio_moduleTraitCor_NV_df, Rubio_moduleTraitPvalue_NV_df, by = "mod_names")

# Significantly correlated modules
Rubio_moduleTraitCor_Pval_NV_df[order(Rubio_moduleTraitCor_Pval_NV_df$pvalue),]
class(Rubio_moduleTraitCor_Pval_NV_df$pvalue) # numeric
Rubio_moduleTraitCor_Pval_NV_df_sig <- Rubio_moduleTraitCor_Pval_NV_df %>% filter(pvalue <= 0.05)
Rubio_moduleTraitCor_Pval_NV_df_sig # 8

Rubio_moduleTraitCor_V_df <- Rubio_moduleTraitCor_df[,c("mod_names","Group.Virulent.vs.Control")]
Rubio_moduleTraitPvalue_V_df <- Rubio_moduleTraitPvalue_df[,c("mod_names","Group.Virulent.vs.Control")]
colnames(Rubio_moduleTraitPvalue_V_df)[2] <- "pvalue"
Rubio_moduleTraitCor_Pval_V_df <- join(Rubio_moduleTraitCor_V_df, Rubio_moduleTraitPvalue_V_df, by = "mod_names")

# Significantly correlated modules
Rubio_moduleTraitCor_Pval_V_df[order(Rubio_moduleTraitCor_Pval_V_df$pvalue),]
class(Rubio_moduleTraitCor_Pval_V_df$pvalue) # numeric
Rubio_moduleTraitCor_Pval_V_df_sig <- Rubio_moduleTraitCor_Pval_V_df %>% filter(pvalue <= 0.05)
Rubio_moduleTraitCor_Pval_V_df_sig # 9

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
Rubio_moduleTraitCor_Pval_NV_df_sig_list <- Rubio_moduleTraitCor_Pval_NV_df_sig$mod_names
Rubio_moduleTraitCor_Pval_V_df_sig_list <- Rubio_moduleTraitCor_Pval_V_df_sig$mod_names
Rubio_moduleTraitCor_Pval_NV_df_sig_list_rm <- str_remove(Rubio_moduleTraitCor_Pval_NV_df_sig_list, "ME")
Rubio_moduleTraitCor_Pval_V_df_sig_list_rm <- str_remove(Rubio_moduleTraitCor_Pval_V_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
names(Rubio_moduleTraitCor_Pval_NV_df_sig_list_rm) <- c( "brown" ,"red" ,"violet","paleturquoise", "cyan","green","darkturquoise", "turquoise")
names(Rubio_moduleTraitCor_Pval_V_df_sig_list_rm) <- c("brown","red","blue","orange","cyan","green","darkturquoise", "turquoise","yellow")

matrix_common= Rubio_dds_rlog_matrix_common
moduleColors= Rubio_moduleColors
lookup =  C_gig_rtracklayer_apop_product_final

lookup_mod_apop_CG <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$transcript_id %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
Rubio_NV_module_apop <- lapply(Rubio_moduleTraitCor_Pval_NV_df_sig_list_rm,  lookup_mod_apop_CG)
Rubio_NV_module_apop_df <- do.call(rbind,Rubio_NV_module_apop)
Rubio_NV_module_apop_df$mod_names <- gsub("\\..*","",row.names(Rubio_NV_module_apop_df))
Rubio_NV_module_apop_df$mod_names <- gsub("^","ME",Rubio_NV_module_apop_df$mod_names)
# add module significance
Rubio_NV_module_apop_df <- left_join(Rubio_NV_module_apop_df,Rubio_moduleTraitCor_Pval_NV_df_sig)
Rubio_NV_module_apop_df$exp <- "Rubio_NV"

Rubio_V_module_apop <- lapply(Rubio_moduleTraitCor_Pval_V_df_sig_list_rm,  lookup_mod_apop_CG)
Rubio_V_module_apop_df <- do.call(rbind,Rubio_V_module_apop)
Rubio_V_module_apop_df$mod_names <- gsub("\\..*","",row.names(Rubio_V_module_apop_df))
Rubio_V_module_apop_df$mod_names <- gsub("^","ME",Rubio_V_module_apop_df$mod_names)
# add module significance
Rubio_V_module_apop_df <- left_join(Rubio_V_module_apop_df,Rubio_moduleTraitCor_Pval_V_df_sig)
Rubio_V_module_apop_df$exp <- "Rubio_V"


## Gene relationship to trait and important modules: Gene Significance and Module Membership
# do later

## Intramodular analysis: identifying genes with high GS and MM
# Using the GS and MM measures, we can identify genes that have a high significance for disease challenge
# as well as high module membership in interesting modules. As an example, we look at the brown module 
# that has the highest association with weight. We plot a scatterplot of Gene Significance vs. Module Membership in the brown module:

## IDENTIFY HUB GENES IN EACH SIG MODULE ##

# look at this later

#### deLorg Res ####
# Define numbers of genes and samples
deLorg_Resistant_nGenes = ncol(deLorgeril_Resistant_dds_vst_matrix_common)
deLorg_Resistant_nSamples = nrow(deLorgeril_Resistant_dds_vst_matrix_common)

# Recalculate MEs with color labels
deLorg_Res_MEs0 = moduleEigengenes(deLorgeril_Resistant_dds_vst_matrix_common, deLorg_Res_moduleColors)$eigengenes
deLorg_Res_MEs = orderMEs(deLorg_Res_MEs0)
deLorg_Res_moduleTraitCor = cor(deLorg_Res_MEs, deLorgeril_Resistant_coldata_collapsed_binarize, use = "p");
deLorg_Res_moduleTraitPvalue = corPvalueStudent(deLorg_Res_moduleTraitCor, deLorg_Resistant_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trai
sizeGrWindow(10,6)
# Will display correlations and their p-values
deLorg_Res_textMatrix = paste(signif(deLorg_Res_moduleTraitCor, 2), "\n(",
                         signif(deLorg_Res_moduleTraitPvalue, 1), ")", sep = "");
dim(deLorg_Res_textMatrix) = dim(deLorg_Res_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = deLorg_Res_moduleTraitCor,
               xLabels = names(deLorgeril_Resistant_coldata_collapsed_binarize),
               yLabels = names(deLorg_Res_MEs),
               ySymbols = names(deLorg_Res_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = deLorg_Res_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
deLorg_Res_moduleTraitCor_df <- as.data.frame(deLorg_Res_moduleTraitCor)
deLorg_Res_moduleTraitCor_df$mod_names <- row.names(deLorg_Res_moduleTraitCor_df)
deLorg_Res_moduleTraitCor_df <- deLorg_Res_moduleTraitCor_df[,c("mod_names","Group.AF21_Resistant.vs.AF21_Resistant_control")]
deLorg_Res_moduleTraitPvalue_df <- as.data.frame(deLorg_Res_moduleTraitPvalue)
deLorg_Res_moduleTraitPvalue_df$mod_names <- row.names(deLorg_Res_moduleTraitPvalue_df )
deLorg_Res_moduleTraitPvalue_df <- deLorg_Res_moduleTraitPvalue_df[,c("mod_names","Group.AF21_Resistant.vs.AF21_Resistant_control")]
colnames(deLorg_Res_moduleTraitPvalue_df)[2] <- "pvalue"
deLorg_Res_moduleTraitCor_Pval_df <- join(deLorg_Res_moduleTraitCor_df, deLorg_Res_moduleTraitPvalue_df, by = "mod_names")

# Significantly correlated modules
deLorg_Res_moduleTraitCor_Pval_df[order(deLorg_Res_moduleTraitCor_Pval_df$pvalue),]
class(deLorg_Res_moduleTraitCor_Pval_df$pvalue) # numeric
deLorg_Res_moduleTraitCor_Pval_df_sig <- deLorg_Res_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
deLorg_Res_moduleTraitCor_Pval_df_sig #25
deLorg_Res_moduleTraitCor_Pval_df_sig <- deLorg_Res_moduleTraitCor_Pval_df_sig %>% filter(mod_names != "MEgrey")

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
deLorg_Res_moduleTraitCor_Pval_df_sig_list <- deLorg_Res_moduleTraitCor_Pval_df_sig$mod_names
deLorg_Res_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(deLorg_Res_moduleTraitCor_Pval_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
names(deLorg_Res_moduleTraitCor_Pval_df_sig_list_rm) <- c("orangered4" ,     "skyblue3"        ,"indianred2" ,     "lightgreen"   ,   "floralwhite"     ,"lightskyblue4" ,  "coral4"    ,     "lightslateblue" , "green"          ,
                                                           "turquoise" ,      "antiquewhite2"  , "blue"      ,      "plum1"       ,    "darkolivegreen2", "pink"         ,   "skyblue2" ,      "coral1"        ,  "lightpink2"     ,
                                                           "deeppink1" ,      "salmon4"        , "white"     ,      "greenyellow" ,    "royalblue"      , "brown"  )
matrix_common= deLorgeril_Resistant_dds_vst_matrix_common
moduleColors= deLorg_Res_moduleColors
lookup =   C_gig_rtracklayer_apop_product_final

lookup_mod_apop_CG <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$transcript_id %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
deLorg_Res_module_apop <- lapply(deLorg_Res_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop_CG)
deLorg_Res_module_apop_df <- do.call(rbind,deLorg_Res_module_apop)
deLorg_Res_module_apop_df$mod_names <- gsub("\\..*","",row.names(deLorg_Res_module_apop_df))
deLorg_Res_module_apop_df$mod_names <- gsub("^","ME",deLorg_Res_module_apop_df$mod_names)
# add module significance
deLorg_Res_module_apop_df <- left_join(deLorg_Res_module_apop_df,deLorg_Res_moduleTraitCor_Pval_df_sig)
deLorg_Res_module_apop_df$exp <- "deLorg_Res"

#### deLorg Sus ####
# Define numbers of genes and samples
deLorg_Susceptible_nGenes = ncol(deLorgeril_Susceptible_dds_vst_matrix_common)
deLorg_Susceptible_nSamples = nrow(deLorgeril_Susceptible_dds_vst_matrix_common)

# Recalculate MEs with color labels
deLorg_Sus_MEs0 = moduleEigengenes(deLorgeril_Susceptible_dds_vst_matrix_common, deLorg_Sus_moduleColors)$eigengenes
deLorg_Sus_MEs = orderMEs(deLorg_Sus_MEs0)
deLorg_Sus_moduleTraitCor = cor(deLorg_Sus_MEs, deLorgeril_Susceptible_coldata_collapsed_binarize , use = "p");
deLorg_Sus_moduleTraitPvalue = corPvalueStudent(deLorg_Sus_moduleTraitCor, deLorg_Susceptible_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trai
sizeGrWindow(10,6)
# Will display correlations and their p-values
deLorg_Sus_textMatrix = paste(signif(deLorg_Sus_moduleTraitCor, 2), "\n(",
                              signif(deLorg_Sus_moduleTraitPvalue, 1), ")", sep = "");
dim(deLorg_Sus_textMatrix) = dim(deLorg_Sus_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = deLorg_Sus_moduleTraitCor,
               xLabels = names(deLorgeril_Susceptible_coldata_collapsed_binarize ),
               yLabels = names(deLorg_Sus_MEs),
               ySymbols = names(deLorg_Sus_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = deLorg_Sus_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
deLorg_Sus_moduleTraitCor_df <- as.data.frame(deLorg_Sus_moduleTraitCor)
deLorg_Sus_moduleTraitCor_df$mod_names <- row.names(deLorg_Sus_moduleTraitCor_df)
deLorg_Sus_moduleTraitCor_df <- deLorg_Sus_moduleTraitCor_df[,c("mod_names","Group.AF11_Susceptible.vs.AF11_Susceptible_control")]
deLorg_Sus_moduleTraitPvalue_df <- as.data.frame(deLorg_Sus_moduleTraitPvalue)
deLorg_Sus_moduleTraitPvalue_df$mod_names <- row.names(deLorg_Sus_moduleTraitPvalue_df )
deLorg_Sus_moduleTraitPvalue_df <- deLorg_Sus_moduleTraitPvalue_df[,c("mod_names","Group.AF11_Susceptible.vs.AF11_Susceptible_control")]
colnames(deLorg_Sus_moduleTraitPvalue_df)[2] <- "pvalue"
deLorg_Sus_moduleTraitCor_Pval_df <- join(deLorg_Sus_moduleTraitCor_df, deLorg_Sus_moduleTraitPvalue_df, by = "mod_names")

# Significantly correlated modules
deLorg_Sus_moduleTraitCor_Pval_df[order(deLorg_Sus_moduleTraitCor_Pval_df$pvalue),]
class(deLorg_Sus_moduleTraitCor_Pval_df$pvalue) # numeric
deLorg_Sus_moduleTraitCor_Pval_df_sig <- deLorg_Sus_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
deLorg_Sus_moduleTraitCor_Pval_df_sig # 4

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
deLorg_Sus_moduleTraitCor_Pval_df_sig_list <- deLorg_Sus_moduleTraitCor_Pval_df_sig$mod_names
deLorg_Sus_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(deLorg_Sus_moduleTraitCor_Pval_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
deLorg_Sus_moduleTraitCor_Pval_df_sig_list_rm
names(deLorg_Sus_moduleTraitCor_Pval_df_sig_list_rm) <- c( "purple",       "turquoise" ,   "midnightblue", "blue" )

matrix_common= deLorgeril_Susceptible_dds_vst_matrix_common
moduleColors= deLorg_Sus_moduleColors
lookup =   C_gig_rtracklayer_apop_product_final

lookup_mod_apop_CG <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$transcript_id %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
deLorg_Sus_module_apop <- lapply(deLorg_Sus_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop_CG)
deLorg_Sus_module_apop_df <- do.call(rbind,deLorg_Sus_module_apop)
deLorg_Sus_module_apop_df$mod_names <- gsub("\\..*","",row.names(deLorg_Sus_module_apop_df))
deLorg_Sus_module_apop_df$mod_names <- gsub("^","ME",deLorg_Sus_module_apop_df$mod_names)
# add module significance
deLorg_Sus_module_apop_df <- left_join(deLorg_Sus_module_apop_df,deLorg_Sus_moduleTraitCor_Pval_df_sig)
deLorg_Sus_module_apop_df$exp <- "deLorg_Sus"

#### He ####
# Define numbers of genes and samples
He_nGenes = ncol(He_dds_vst_matrix_common)
He_nSamples = nrow(He_dds_vst_matrix_common)

# Recalculate MEs with color labels
He_MEs0 = moduleEigengenes(He_dds_vst_matrix_common, He_moduleColors)$eigengenes
He_MEs = orderMEs(He_MEs0)
He_moduleTraitCor = cor(He_MEs, He_coldata_collapsed_binarize, use = "p");
He_moduleTraitPvalue = corPvalueStudent(He_moduleTraitCor, He_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trai
sizeGrWindow(10,6)
# Will display correlations and their p-values
He_textMatrix = paste(signif(He_moduleTraitCor, 2), "\n(",
                              signif(He_moduleTraitPvalue, 1), ")", sep = "");
dim(He_textMatrix) = dim(He_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = He_moduleTraitCor,
               xLabels = names(He_coldata_collapsed_binarize),
               yLabels = names(He_MEs),
               ySymbols = names(He_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = He_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
He_moduleTraitCor_df <- as.data.frame(He_moduleTraitCor)
He_moduleTraitCor_df$mod_names <- row.names(He_moduleTraitCor_df)
He_moduleTraitCor_df <- He_moduleTraitCor_df[,c("mod_names","Group.OsHV1.vs.control")]
He_moduleTraitPvalue_df <- as.data.frame(He_moduleTraitPvalue)
He_moduleTraitPvalue_df$mod_names <- row.names(He_moduleTraitPvalue_df )
He_moduleTraitPvalue_df <- He_moduleTraitPvalue_df[,c("mod_names","Group.OsHV1.vs.control")]
colnames(He_moduleTraitPvalue_df)[2] <- "pvalue"
He_moduleTraitCor_Pval_df <- join(He_moduleTraitCor_df, He_moduleTraitPvalue_df, by = "mod_names")

# Significantly correlated modules
He_moduleTraitCor_Pval_df[order(He_moduleTraitCor_Pval_df$pvalue),]
class(He_moduleTraitCor_Pval_df$pvalue) # numeric
He_moduleTraitCor_Pval_df_sig <- He_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
He_moduleTraitCor_Pval_df_sig # 22
# remove grey from list
He_moduleTraitCor_Pval_df_sig <- He_moduleTraitCor_Pval_df_sig %>% filter(mod_names !="MEgrey")

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
He_moduleTraitCor_Pval_df_sig_list <- He_moduleTraitCor_Pval_df_sig$mod_names
He_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(He_moduleTraitCor_Pval_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
names(He_moduleTraitCor_Pval_df_sig_list_rm) <- c("black"       ,"skyblue"  ,"darkolivegreen" , "magenta" , "white"   ,"darkgreen", "darkslateblue","salmon4",         "red" ,           
                                                   "tan"        , "yellow"  , "midnightblue"  ,  "salmon" ,  "lightsteelblue1", "brown4"  ,  "orangered4"  , "brown" ,          "blue",           
                                                   "darkorange2", "darkgrey", "turquoise"  )
matrix_common= He_dds_vst_matrix_common
moduleColors= He_moduleColors
lookup =   C_gig_rtracklayer_apop_product_final

lookup_mod_apop_CG <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$transcript_id %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
He_module_apop <- lapply(He_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop_CG)
He_module_apop_df <- do.call(rbind,He_module_apop)
He_module_apop_df$mod_names <- gsub("\\..*","",row.names(He_module_apop_df))
He_module_apop_df$mod_names <- gsub("^","ME",He_module_apop_df$mod_names)
# add module significance
He_module_apop_df <- left_join(He_module_apop_df,He_moduleTraitCor_Pval_df_sig)
He_module_apop_df$exp <- "He"

#### MASTER C_GIG MODULE-TRAIT-SIG APOP GENES ####
colnames(Zhang_LPS_module_apop_df)[5] <- "mod_signif"
colnames(Zhang_module_Vibrio_apop_df)[5] <- "mod_signif"
colnames(Rubio_NV_module_apop_df   )[5] <- "mod_signif"
colnames(Rubio_V_module_apop_df    )[5] <- "mod_signif"
colnames(deLorg_Res_module_apop_df )[5] <- "mod_signif"
colnames(deLorg_Sus_module_apop_df )[5] <- "mod_signif"   
colnames(He_module_apop_df)[5] <- "mod_signif" 

C_gig_all_exp_mod_sig_apop <- rbind(Zhang_LPS_module_apop_df,
                                    Zhang_module_Vibrio_apop_df,
                                    Rubio_NV_module_apop_df  ,
                                    Rubio_V_module_apop_df   ,
                                    deLorg_Res_module_apop_df,
                                    deLorg_Sus_module_apop_df,
                                    He_module_apop_df)
nrow(C_gig_all_exp_mod_sig_apop) # 1337
C_gig_all_exp_mod_sig_apop_positive <- C_gig_all_exp_mod_sig_apop %>% filter(mod_signif >0)

#### MEASURE MODULE PRESERVATION BETWEEN C. GIGAS DIFFERENT EXPERIMENTS ####

# Column names agree between the datasets, so we can go ahead and measure module preservation between all 
## 1. MEASURE OVERLAP BETWEEN Zhang and Rubio, and delorgeril and He, and then comparing bacteria and virus 

# compare more to Rubio for the bacterial challenge

# Assess whether Delorgeril_Sus and He have conservation (since He was also susceptible)
deLorg_He_multiExpr = list(deLorg_Sus=list(data=deLorgeril_Susceptible_dds_vst_matrix_common),He=list(data=He_dds_vst_matrix_common))
deLorg_He_multiColor = list(deLorg_Sus = deLorg_Sus_moduleColors) 
#deLorg_He_mp =modulePreservation(deLorg_He_multiExpr, deLorg_He_multiColor,
#                                   referenceNetworks=1,
#                                   verbose=3,
#                                   networkType="signed hybrid", # use same signed hybrid as before
#                                   corFnc = "bicor", # use recommended bicor as before 
#                                   nPermutations=100,
#                                   randomSeed = 1, # recommended in langfelder tutorial
#                                   quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
deLorg_He_stats = deLorg_He_mp$preservation$Z$ref.deLorg_Sus$inColumnsAlsoPresentIn.He
deLorg_He_stats_order <- deLorg_He_stats[order(-deLorg_He_stats[,2]),c(1:2)]
deLorg_He_stats_preserved <- deLorg_He_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
# mod_name moduleSize Zsummary.pres
# 1      green       1000     33.794665
# 2     yellow       1000     30.735406
# 3    magenta        356     24.497166
# 4        red        526     23.668887
# 5      brown       1000     21.318501
# 6       pink        406     16.401608
# 7       blue       1000     15.766808
# 8       gold       1000     11.828504
# 9  turquoise       1000     11.567648
# 10     black        484      5.603149
deLorg_He_stats_MR <- deLorg_He_mp$preservation$observed$ref.deLorg_Sus$inColumnsAlsoPresentIn.He
deLorg_He_stats_MR <- deLorg_He_stats_MR  [,c("moduleSize","medianRank.pres")]
deLorg_He_stats_MR_order <- deLorg_He_stats_MR [order(-deLorg_He_stats_MR [,2]),]
deLorg_He_stats_MR_order_less20 <- deLorg_He_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save module preservation results
#save(deLorg_He_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/deLorg_Sus_He_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(deLorg_He_stats_MR_order_less20, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

deLorg_He_stats_preserved_Zsum_medianRank <- deLorg_He_stats_MR_order_less20[deLorg_He_stats_MR_order_less20$mod_name %in% deLorg_He_stats_preserved$mod_name,]

# Were any of these preserved modules significant in deLorg_Sus and He
deLorg_Sus_preserved_all <- deLorg_He_stats_preserved[deLorg_He_stats_preserved$mod_name %in% deLorg_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#7      blue       1000      15.76681
#9 turquoise       1000      11.56765

He_preserved_all <- deLorg_He_stats_preserved[deLorg_He_stats_preserved$mod_name %in% He_moduleTraitCor_Pval_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#2     yellow       1000     30.735406
#3    magenta        356     24.497166
#4        red        526     23.668887
#5      brown       1000     21.318501
#7       blue       1000     15.766808
#9  turquoise       1000     11.567648
#10     black        484      5.603149

## He vs. Delorgeril_Sus 
He_deLorg_Sus_multiExpr = list(He=list(data=He_dds_vst_matrix_common),deLorg_Sus=list(data=deLorgeril_Susceptible_dds_vst_matrix_common))
He_deLorg_Sus_multiColor = list(He = He_moduleColors) 
#He_deLorg_Sus_mp =modulePreservation(He_deLorg_Sus_multiExpr, He_deLorg_Sus_multiColor,
#                                   referenceNetworks=1,
#                                   verbose=3,
#                                   networkType="signed hybrid", # use same signed hybrid as before
#                                   corFnc = "bicor", # use recommended bicor as before 
#                                   nPermutations=100,
#                                   randomSeed = 1, # recommended in langfelder tutorial
#                                   quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
He_deLorg_Sus_stats = He_deLorg_Sus_mp$preservation$Z$ref.He$inColumnsAlsoPresentIn.deLorg_Sus
He_deLorg_Sus_stats_order <- He_deLorg_Sus_stats[order(-He_deLorg_Sus_stats[,2]),c(1:2)]
He_deLorg_Sus_stats_preserved <- He_deLorg_Sus_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
#mod_name moduleSize Zsummary.pres
#1              red        450     40.851138
#2           yellow        794     36.369011
#3              tan        370     29.637540
#4          magenta        387     23.378484
#5    mediumpurple3         86     21.048026
#6           violet        132     19.980389
#7            brown        819     19.910803
#8     midnightblue        339     19.559415
#9        turquoise       1000     19.230363
#10       darkgreen        248     18.663438
#11            blue       1000     18.180394
#12         skyblue        163     17.047495
#13     saddlebrown        154     15.920245
#14            cyan        348     14.604253
#15            pink        398     13.196043
#16          salmon        365     12.450763
#17        skyblue3         93     11.335303
#18           white        178     11.253941
#19            gold       1000     11.210571
#20       lightcyan        339     11.128849
#21        darkgrey        209     10.404264
#22   darkslateblue         49     10.213339
#23        thistle2         46      9.137793
#24 lightsteelblue1         79      7.507454
#25   darkturquoise        219      6.854947
#26           black        413      6.690562
#27       steelblue        140      6.541944
#28         salmon4         38      6.498236
#29      orangered4         86      6.296053
#30          purple        386      5.698736
#31           green        682      5.510336
#32           ivory         72      5.070521
He_deLorg_Sus_stats_MR <- He_deLorg_Sus_mp$preservation$observed$ref.He$inColumnsAlsoPresentIn.deLorg_Sus
He_deLorg_Sus_stats_MR <- He_deLorg_Sus_stats_MR  [,c("moduleSize","medianRank.pres")]
He_deLorg_Sus_stats_MR_order <- He_deLorg_Sus_stats_MR [order(-He_deLorg_Sus_stats_MR [,2]),]
He_deLorg_Sus_stats_MR_order_less20 <- He_deLorg_Sus_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save module preservation results
# save(He_deLorg_Sus_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/He_deLorg_Sus_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(He_deLorg_Sus_stats_MR_order_less20, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

He_deLorg_Sus_stats_preserved_Zsum_medianRank <- He_deLorg_Sus_stats_MR_order_less20[He_deLorg_Sus_stats_MR_order_less20$mod_name %in% He_deLorg_Sus_stats_preserved$mod_name,]

# Were any of these preserved modules significant in deLorg_Sus and He
He_deLorg_Sus_preserved_all <- He_deLorg_Sus_stats_preserved[He_deLorg_Sus_stats_preserved$mod_name %in% deLorg_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#8  midnightblue        339     19.559415
#9     turquoise       1000     19.230363
#11         blue       1000     18.180394
#30       purple        386      5.698736
He_deLorg_Sus_he_preserved_all <- He_deLorg_Sus_stats_preserved[He_deLorg_Sus_stats_preserved$mod_name %in% He_moduleTraitCor_Pval_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#1              red        450     40.851138
#2           yellow        794     36.369011
#3              tan        370     29.637540
#4          magenta        387     23.378484
#7            brown        819     19.910803
#8     midnightblue        339     19.559415
#9        turquoise       1000     19.230363
#10       darkgreen        248     18.663438
#11            blue       1000     18.180394
#12         skyblue        163     17.047495
#16          salmon        365     12.450763
#18           white        178     11.253941
#21        darkgrey        209     10.404264
#22   darkslateblue         49     10.213339
#24 lightsteelblue1         79      7.507454
#26           black        413      6.690562
#28         salmon4         38      6.498236
#29      orangered4         86      6.296053

## He vs. Delorgeril_Tol 
He_deLorg_Tol_multiExpr = list(He=list(data=He_dds_vst_matrix_common),deLorg_Tol=list(data=deLorgeril_Resistant_dds_vst_matrix_common))
He_deLorg_Tol_multiColor = list(He = He_moduleColors) 
He_deLorg_Tol_mp =modulePreservation(He_deLorg_Tol_multiExpr, He_deLorg_Tol_multiColor,
                                   referenceNetworks=1,
                                   verbose=3,
                                   networkType="signed hybrid", # use same signed hybrid as before
                                   corFnc = "bicor", # use recommended bicor as before 
                                   nPermutations=100,
                                   randomSeed = 1, # recommended in langfelder tutorial
                                   quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
He_deLorg_Tol_stats = He_deLorg_Tol_mp$preservation$Z$ref.He$inColumnsAlsoPresentIn.deLorg_Tol
He_deLorg_Tol_stats_order <- He_deLorg_Tol_stats[order(-He_deLorg_Tol_stats[,2]),c(1:2)]
He_deLorg_Tol_stats_preserved <- He_deLorg_Tol_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
He_deLorg_Tol_stats_MR <- He_deLorg_Tol_mp$preservation$observed$ref.He$inColumnsAlsoPresentIn.deLorg_Tol
He_deLorg_Tol_stats_MR <- He_deLorg_Tol_stats_MR  [,c("moduleSize","medianRank.pres")]
He_deLorg_Tol_stats_MR_order <- He_deLorg_Tol_stats_MR [order(-He_deLorg_Sus_stats_MR [,2]),]
He_deLorg_Tol_stats_MR_order_less20 <- He_deLorg_Tol_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save module preservation results
save(He_deLorg_Tol_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/He_deLorg_Tol_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(He_deLorg_Tol_stats_MR_order_less20, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

He_deLorg_Tol_stats_preserved_Zsum_medianRank <- He_deLorg_Tol_stats_MR_order_less20[He_deLorg_Tol_stats_MR_order_less20$mod_name %in% He_deLorg_Tol_stats_preserved$mod_name,]

# Were any of these preserved modules significant in deLorg_Tol and He
He_deLorg_Tol_preserved_all <-    He_deLorg_Tol_stats_preserved[He_deLorg_Tol_stats_preserved$mod_name %in% deLorg_Res_moduleTraitCor_Pval_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#6         blue       1000     17.918548
#11   turquoise       1000     15.434665
#18    skyblue3         93     12.944139
#19 greenyellow        381     12.500035
#22        pink        398     10.192308
#25  lightgreen        309      9.194302
#29       white        178      7.437909
#30       brown        819      7.203497
#31 floralwhite         71      6.806136
#35     salmon4         38      5.885579
He_deLorg_Tol_he_preserved_all <- He_deLorg_Tol_stats_preserved[He_deLorg_Tol_stats_preserved$mod_name %in% He_moduleTraitCor_Pval_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#1              red        450     76.257779
#2              tan        370     41.447376
#3          skyblue        163     26.891464
#4          magenta        387     20.875964
#5           yellow        794     19.879381
#6             blue       1000     17.918548
#7     midnightblue        339     16.942338
#10       darkgreen        248     16.225612
#11       turquoise       1000     15.434665
#14          salmon        365     14.377787
#23           black        413     10.020731
#26        darkgrey        209      9.108802
#28 lightsteelblue1         79      7.565908
#29           white        178      7.437909
#30           brown        819      7.203497
#32   darkslateblue         49      6.563884
#35         salmon4         38      5.885579
#36  darkolivegreen        130      5.609281


## Delorgeril_Tol vs. He
deLorg_Tol_He_multiExpr = list(deLorg_Tol=list(data=deLorgeril_Resistant_dds_vst_matrix_common), He=list(data=He_dds_vst_matrix_common))
deLorg_Tol_He_multiColor = list(deLorg_Tol = deLorg_Res_moduleColors) 
deLorg_Tol_He_mp =modulePreservation(deLorg_Tol_He_multiExpr, deLorg_Tol_He_multiColor,
                                     referenceNetworks=1,
                                     verbose=3,
                                     networkType="signed hybrid", # use same signed hybrid as before
                                     corFnc = "bicor", # use recommended bicor as before 
                                     nPermutations=100,
                                     randomSeed = 1, # recommended in langfelder tutorial
                                     quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
deLorg_Tol_He_stats = deLorg_Tol_He_mp$preservation$Z$ref.deLorg_Tol$inColumnsAlsoPresentIn.He
deLorg_Tol_He_stats_order <- deLorg_Tol_He_stats[order(-deLorg_Tol_He_stats[,2]),c(1:2)]
deLorg_Tol_He_stats_preserved <- deLorg_Tol_He_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
deLorg_Tol_He_stats_MR <- deLorg_Tol_He_mp$preservation$observed$ref.deLorg_Tol$inColumnsAlsoPresentIn.He
deLorg_Tol_He_stats_MR <- deLorg_Tol_He_stats_MR [,c("moduleSize","medianRank.pres")]
deLorg_Tol_He_stats_MR_order <- deLorg_Tol_He_stats_MR [order(-deLorg_Tol_He_stats_MR [,2]),]
deLorg_Tol_He_stats_MR_order_less20 <- deLorg_Tol_He_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save module preservation results
save(deLorg_Tol_He_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/deLorg_Tol_He_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(deLorg_Tol_He_stats_MR_order_less20, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

deLorg_Tol_He_stats_preserved_Zsum_medianRank <- deLorg_Tol_He_stats_MR_order_less20[deLorg_Tol_He_stats_MR_order_less20$mod_name %in% deLorg_Tol_He_stats_preserved$mod_name,]

# Were any of these preserved modules significant in deLorg_Tol and He
deLorg_Tol_He_Tol_preserved_all <-    deLorg_Tol_He_stats_preserved[deLorg_Tol_He_stats_preserved$mod_name %in% deLorg_Res_moduleTraitCor_Pval_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#2         blue       1000      31.03493
#5    turquoise       1000      18.49194
#6         pink        525      17.33296
#9        brown       1000      14.28718
#11 greenyellow        475      13.79025
#12       white        110      12.59825
#15       green        690      10.34434

deLorg_Tol_He_He_preserved_all <- deLorg_Tol_He_stats_preserved[deLorg_Tol_He_stats_preserved$mod_name %in% He_moduleTraitCor_Pval_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#1         yellow        870     35.478220
#2           blue       1000     31.034927
#3          black        566     21.762627
#4            red        680     21.071764
#5      turquoise       1000     18.491943
#9          brown       1000     14.287180
#10           tan        475     14.277122
#12         white        110     12.598251
#16     darkgreen        165      9.701913
#17       magenta        510      9.167497
#22 darkslateblue         69      5.681363
#24  midnightblue        355      5.017819

# Assess whether deLorg_Sus and deLorgTol have conservation 
#deLorg_Sus_Tol_multiExpr = list(deLorg_Sus=list(data=deLorgeril_Susceptible_dds_vst_matrix_common),deLorg_Tol=list(data=deLorgeril_Resistant_dds_vst_matrix_common))
#deLorg_Sus_Tol_multiColor = list(deLorg_Sus = deLorg_Sus_moduleColors) 
##deLorg_Sus_Tol_mp =modulePreservation(deLorg_Sus_Tol_multiExpr, deLorg_Sus_Tol_multiColor,
##                                 referenceNetworks=1,
##                                 verbose=3,
##                                 networkType="signed hybrid", # use same signed hybrid as before
##                                 corFnc = "bicor", # use recommended bicor as before 
##                                 nPermutations=100,
##                                 randomSeed = 1, # recommended in langfelder tutorial
##                                 quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
#deLorg_Sus_Tol_stats = deLorg_Sus_Tol_mp$preservation$Z$ref.deLorg_Sus$inColumnsAlsoPresentIn.deLorg_Tol
#deLorg_Sus_Tol_stats_order <- deLorg_Sus_Tol_stats[order(-deLorg_Sus_Tol_stats[,2]),c(1:2)]
#deLorg_Sus_Tol_stats_preserved <- deLorg_Sus_Tol_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
##mod_name moduleSize Zsummary.pres
##1       yellow       1000     49.473174
##2        green       1000     34.124863
##3    turquoise       1000     20.453268
##4         blue       1000     17.552898
##5          red        526     16.676912
##6         pink        406     13.949355
##7         gold       1000     13.572896
##8      magenta        356     12.053217
##9        brown       1000      9.691178
##10 greenyellow        164      8.095585
##11      purple        208      6.741922
##12       black        484      6.005315
#
#deLorg_Sus_Tol_stats_MR <- deLorg_Sus_Tol_mp$preservation$observed$ref.deLorg_Sus$inColumnsAlsoPresentIn.deLorg_Tol
#deLorg_Sus_Tol_stats_MR <- deLorg_Sus_Tol_stats_MR[,c("moduleSize","medianRank.pres")]
#deLorg_Sus_Tol_stats_MR_order <- deLorg_Sus_Tol_stats_MR[order(-deLorg_Sus_Tol_stats_MR[,2]),]
#deLorg_Sus_Tol_stats_MR_order_less20 <- deLorg_Sus_Tol_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)
#
## Save module preservation results
##save(deLorg_Sus_Tol_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/deLorg_Sus_deLorg_Res_module_preservation.RData")
#
## Plot module preservation
## remember - low median rank means high preservation
#ggplot(deLorg_Sus_Tol_stats_MR_order_less20 , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
#  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))
#
#deLorg_Sus_Tol_stats_preserved_Zsum_medianRank <- deLorg_Sus_Tol_stats_MR_order_less20[deLorg_Sus_Tol_stats_MR_order_less20$mod_name %in% deLorg_Sus_Tol_stats_preserved$mod_name,]
#
## Were any of these preserved modules significant in deLorg_Sus and He
#deLorg_Sus_Tol_Suspreserved_all <- deLorg_Sus_Tol_stats_preserved[deLorg_Sus_Tol_stats_preserved$mod_name %in% deLorg_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
## mod_name moduleSize Zsummary.pres
## 3  turquoise       1000     20.453268
## 4       blue       1000     17.552898
## 11    purple        208      6.741922
#deLorg_Sus_Tol_Tol_preserved_all <- deLorg_Sus_Tol_stats_preserved[deLorg_Sus_Tol_stats_preserved$mod_name %in% deLorg_Res_moduleTraitCor_Pval_df_sig_list_rm,]
## mod_name moduleSize Zsummary.pres
## 2        green       1000     34.124863
## 3    turquoise       1000     20.453268
## 4         blue       1000     17.552898
## 6         pink        406     13.949355
## 9        brown       1000      9.691178
## 10 greenyellow        164      8.095585

## deLorg_Tol and deLorg_Sus ##  
deLorg_Tol_Sus_multiExpr = list(deLorg_Tol=list(data=deLorgeril_Resistant_dds_vst_matrix_common), deLorg_Sus=list(data=deLorgeril_Susceptible_dds_vst_matrix_common))
deLorg_Tol_Sus_multiColor = list(deLorg_Tol = deLorg_Res_moduleColors) 
deLorg_Tol_Sus_mp =modulePreservation(deLorg_Tol_Sus_multiExpr, deLorg_Tol_Sus_multiColor,
                                 referenceNetworks=1,
                                 verbose=3,
                                 networkType="signed hybrid", # use same signed hybrid as before
                                 corFnc = "bicor", # use recommended bicor as before 
                                 nPermutations=100,
                                 randomSeed = 1, # recommended in langfelder tutorial
                                 quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
deLorg_Tol_Sus_stats = deLorg_Tol_Sus_mp$preservation$Z$ref.deLorg_Tol$inColumnsAlsoPresentIn.deLorg_Sus
deLorg_Tol_Sus_stats_order <- deLorg_Tol_Sus_stats[order(-deLorg_Tol_Sus_stats[,2]),c(1:2)]
deLorg_Tol_Sus_stats_preserved <- deLorg_Tol_Sus_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
deLorg_Tol_Sus_stats_MR <- deLorg_Tol_Sus_mp$preservation$observed$ref.deLorg_Tol$inColumnsAlsoPresentIn.deLorg_Sus
deLorg_Tol_Sus_stats_MR <- deLorg_Tol_Sus_stats_MR[,c("moduleSize","medianRank.pres")]
deLorg_Tol_Sus_stats_MR_order <- deLorg_Tol_Sus_stats_MR[order(-deLorg_Tol_Sus_stats_MR[,2]),]
deLorg_Tol_Sus_stats_MR_order_less20 <- deLorg_Tol_Sus_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save module preservation results
save(deLorg_Tol_Sus_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/deLorg_Tol_Sus_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(deLorg_Tol_Sus_stats_MR_order_less20 , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

deLorg_Tol_Sus_stats_preserved_Zsum_medianRank <- deLorg_Tol_Sus_stats_MR_order_less20[deLorg_Tol_Sus_stats_MR_order_less20$mod_name %in% deLorg_Tol_Sus_stats_preserved$mod_name,]

# Were any of these preserved modules significant in deLorg_Sus and He
deLorg_Tol_Sus_Suspreserved_all <-  deLorg_Tol_Sus_stats_preserved[deLorg_Tol_Sus_stats_preserved$mod_name %in% deLorg_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#1       blue       1000     43.977505
#4  turquoise       1000     24.620722
#21    purple        481      7.691346
deLorg_Tol_Sus_Tol_preserved_all <- deLorg_Tol_Sus_stats_preserved[deLorg_Tol_Sus_stats_preserved$mod_name %in% deLorg_Res_moduleTraitCor_Pval_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#1           blue       1000     43.977505
#3          brown       1000     26.525211
#4      turquoise       1000     24.620722
#5          green        690     22.984064
#7    greenyellow        475     17.587972
#9          white        110     14.414030
#19          pink        525      8.806550
#23 antiquewhite2         42      6.659519

# Assess whether Rubio and He have conservation (since He was also susceptible)
Rubio_He_multiExpr = list(Rubio=list(data=Rubio_dds_rlog_matrix_common),He=list(data=He_dds_vst_matrix_common))
Rubio_He_multiColor = list(Rubio = Rubio_moduleColors) 
#Rubio_He_mp =modulePreservation(Rubio_He_multiExpr, Rubio_He_multiColor,
#                                      referenceNetworks=1,
#                                      verbose=3,
#                                      networkType="signed hybrid", # use same signed hybrid as before
#                                      corFnc = "bicor", # use recommended bicor as before 
#                                      nPermutations=100,
#                                      randomSeed = 1, # recommended in langfelder tutorial
#                                      quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Rubio_He_stats = Rubio_He_mp$preservation$Z$ref.Rubio$inColumnsAlsoPresentIn.He
Rubio_He_stats_order <- Rubio_He_stats[order(-Rubio_He_stats[,2]),c(1:2)]
Rubio_He_stats_preserved <- Rubio_He_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
#mod_name moduleSize Zsummary.pres
#1    greenyellow        857     30.926061
#2          brown       1000     20.578240
#3          green       1000     18.799122
#4            tan        642     17.368799
#5      turquoise       1000     16.208233
#6    lightyellow        322     14.832272
#7      lightcyan        367     12.060367
#8  darkturquoise        249     11.250846
#9           gold       1000     11.036436
#10          blue       1000     10.921887
#11    lightgreen        334      9.823507
#12          cyan        527      9.722070
#13           red       1000      8.040630
#14         black       1000      6.498954
#15       skyblue         93      5.785529
#16        salmon        633      5.172902
#17        orange        166      5.028439
#18       darkred        271      5.012814
Rubio_He_stats_MR <- Rubio_He_mp$preservation$observed$ref.Rubio$inColumnsAlsoPresentIn.He
Rubio_He_stats_MR <- Rubio_He_stats_MR[,c("moduleSize","medianRank.pres")]
Rubio_He_stats_MR_order <- Rubio_He_stats_MR[order(-Rubio_He_stats_MR [,2]),]
Rubio_He_stats_MR_order_less20 <- Rubio_He_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save module preservation results
# save(Rubio_He_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Rubio_He_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(Rubio_He_stats_MR_order_less20, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

Rubio_He_stats_preserved_Zsum_medianRank <- Rubio_He_stats_MR_order_less20[Rubio_He_stats_MR_order_less20$mod_name %in% Rubio_He_stats_preserved,]

# Were any of these preserved modules significant in deLorg_Sus and He
Rubio_He_NV_preserved_all <- Rubio_He_stats_preserved[Rubio_He_stats_preserved$mod_name %in% Rubio_moduleTraitCor_Pval_NV_df_sig_list_rm,]
# mod_name moduleSize Zsummary.pres
# 2          brown       1000      20.57824
# 3          green       1000      18.79912
# 5      turquoise       1000      16.20823
# 8  darkturquoise        249      11.25085
# 12          cyan        527       9.72207
# 13           red       1000       8.04063
Rubio_He_V_preserved_all <- Rubio_He_stats_preserved[Rubio_He_stats_preserved$mod_name %in% Rubio_moduleTraitCor_Pval_V_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#2          brown       1000     20.578240
#3          green       1000     18.799122
#5      turquoise       1000     16.208233
#8  darkturquoise        249     11.250846
#10          blue       1000     10.921887
#12          cyan        527      9.722070
#13           red       1000      8.040630
#17        orange        166      5.028439
He_preserved_all<- Rubio_He_stats_preserved[Rubio_He_stats_preserved$mod_name %in% He_moduleTraitCor_Pval_df_sig_list_rm,]
# mod_name moduleSize Zsummary.pres
# 2      brown       1000     20.578240
# 4        tan        642     17.368799
# 5  turquoise       1000     16.208233
# 10      blue       1000     10.921887
# 13       red       1000      8.040630
# 14     black       1000      6.498954
# 15   skyblue         93      5.785529
# 16    salmon        633      5.172902

# He vs Rubio 
He_Rubio_multiExpr = list(He=list(data=He_dds_vst_matrix_common), Rubio=list(data=Rubio_dds_rlog_matrix_common))
He_Rubio_multiColor = list(He = He_moduleColors) 
#He_Rubio_mp =modulePreservation(He_Rubio_multiExpr, He_Rubio_multiColor,
#                                      referenceNetworks=1,
#                                      verbose=3,
#                                      networkType="signed hybrid", # use same signed hybrid as before
#                                      corFnc = "bicor", # use recommended bicor as before 
#                                      nPermutations=100,
#                                      randomSeed = 1, # recommended in langfelder tutorial
#                                      quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
He_Rubio_stats = He_Rubio_mp$preservation$Z$ref.He$inColumnsAlsoPresentIn.Rubio
He_Rubio_stats_order <- He_Rubio_stats[order(-He_Rubio_stats[,2]),c(1:2)]
He_Rubio_stats_preserved <- He_Rubio_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
#mod_name moduleSize Zsummary.pres
#1          magenta        387     27.535801
#2           yellow        794     20.153655
#3           salmon        365     19.019356
#4             blue       1000     17.397045
#5     midnightblue        339     13.156157
#6      floralwhite         71     12.471832
#7             cyan        348     11.657512
#8           purple        386     11.358633
#9         skyblue3         93     11.356415
#10             tan        370     10.813415
#11       turquoise       1000     10.090613
#12           white        178      9.831406
#13           brown        819      9.601788
#14     saddlebrown        154      9.385677
#15         skyblue        163      9.035803
#16            gold       1000      9.007044
#17            pink        398      8.988689
#18          violet        132      8.772543
#19 lightsteelblue1         79      8.240611
#20       lightcyan        339      8.207468
#21     greenyellow        381      7.775783
#22   darkslateblue         49      6.350780
#23       steelblue        140      5.989767
#24         salmon4         38      5.839396
#25   paleturquoise        135      5.401110
#26          orange        195      5.114006
#27             red        450      5.082108
#28       darkgreen        248      5.054699

He_Rubio_stats_MR <- He_Rubio_mp$preservation$observed$ref.He$inColumnsAlsoPresentIn.Rubio
He_Rubio_stats_MR <-He_Rubio_stats_MR[,c("moduleSize","medianRank.pres")]
He_Rubio_stats_MR_order <- He_Rubio_stats_MR[order(-He_Rubio_stats_MR [,2]),]
He_Rubio_stats_MR_order_less20 <- He_Rubio_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save module preservation results
#save(He_Rubio_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/He_Rubio_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(He_Rubio_stats_MR_order_less20 , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

He_Rubio_stats_preserved_Zsum_medianRank <- Rubio_He_stats_MR_order_less20[Rubio_He_stats_MR_order_less20$mod_name %in% He_Rubio_stats_preserved$mod_name,]

# Were any of these preserved modules significant in deLorg_Sus and He
He_Rubio_NV_preserved_all <- He_Rubio_stats_preserved[He_Rubio_stats_preserved$mod_name %in% Rubio_moduleTraitCor_Pval_NV_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#7           cyan        348     11.657512
#11     turquoise       1000     10.090613
#13         brown        819      9.601788
#18        violet        132      8.772543
#25 paleturquoise        135      5.401110
#27           red        450      5.082108
He_Rubio_V_preserved_all <-  He_Rubio_stats_preserved[He_Rubio_stats_preserved$mod_name %in% Rubio_moduleTraitCor_Pval_V_df_sig_list_rm,]
# mod_name moduleSize Zsummary.pres
# 2     yellow        794     20.153655
# 4       blue       1000     17.397045
# 7       cyan        348     11.657512
# 11 turquoise       1000     10.090613
# 13     brown        819      9.601788
# 26    orange        195      5.114006
# 27       red        450      5.082108
He_Rubio_He_preserved_all<-  He_Rubio_stats_preserved[He_Rubio_stats_preserved$mod_name %in% He_moduleTraitCor_Pval_df_sig_list_rm,]
# mod_name moduleSize Zsummary.pres
# 1          magenta        387     27.535801
# 2           yellow        794     20.153655
# 3           salmon        365     19.019356
# 4             blue       1000     17.397045
# 5     midnightblue        339     13.156157
# 10             tan        370     10.813415
# 11       turquoise       1000     10.090613
# 12           white        178      9.831406
# 13           brown        819      9.601788
# 15         skyblue        163      9.035803
# 19 lightsteelblue1         79      8.240611
# 22   darkslateblue         49      6.350780
# 24         salmon4         38      5.839396
# 27             red        450      5.082108
# 28       darkgreen        248      5.054699


## Assess whether Zhang and He have conservation (since He was also susceptible) ##
Zhang_He_multiExpr = list(Zhang=list(data=Zhang_dds_rlog_matrix_common),He=list(data=He_dds_vst_matrix_common))
Zhang_He_multiColor = list(Zhang = Zhang_moduleColors) 
#Zhang_He_mp =modulePreservation(Zhang_He_multiExpr, Zhang_He_multiColor,
#                                referenceNetworks=1,
#                                verbose=3,
#                                networkType="signed hybrid", # use same signed hybrid as before
#                                corFnc = "bicor", # use recommended bicor as before 
#                                nPermutations=100,
#                                randomSeed = 1, # recommended in langfelder tutorial
#                                quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Zhang_He_stats = Zhang_He_mp$preservation$Z$ref.Zhang$inColumnsAlsoPresentIn.He
Zhang_He_stats_order <- Zhang_He_stats[order(-Zhang_He_stats[,2]),c(1:2)]
Zhang_He_stats_preserved <- Zhang_He_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
#mod_name moduleSize Zsummary.pres
#1          blue       1000     15.979158
#2       magenta        588     12.212596
#3        yellow        941      9.095257
#4     turquoise       1000      8.963291
#5 paleturquoise        219      7.572312
#6        purple        531      5.964199

Zhang_He_stats_MR <- Zhang_He_mp$preservation$observed$ref.Zhang$inColumnsAlsoPresentIn.He
Zhang_He_stats_MR <- Zhang_He_stats_MR[,c("moduleSize","medianRank.pres")]
Zhang_He_stats_MR_order <- Zhang_He_stats_MR[order(-Zhang_He_stats_MR[,2]),]
Zhang_He_stats_MR_order_less20 <- Zhang_He_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save module preservation results
#save(Zhang_He_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Zhang_He_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(Zhang_He_stats_MR_order_less20, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

Zhang_He_stats_preserved_Zsum_medianRank <- Zhang_He_stats_MR_order_less20[Zhang_He_stats_MR_order_less20$mod_name %in% Zhang_He_stats_preserved$mod_name, ]

# Were any of these preserved modules significant in deLorg_Sus and He
Zhang_He_Zhang_LPS_preserved_all <- Zhang_He_stats_preserved[Zhang_He_stats_preserved$mod_name %in% Zhang_moduleTraitCor_Pval_LPS_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#1          blue       1000     15.979158
#2       magenta        588     12.212596
#4     turquoise       1000      8.963291
#5 paleturquoise        219      7.572312

Zhang_He_Zhang_Vibrio_preserved_all <- Zhang_He_stats_preserved[Zhang_He_stats_preserved$mod_name %in% Zhang_moduleTraitCor_Pval_Vibrio_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#1     blue       1000     15.979158
#6   purple        531      5.964199

Zhang_He_preserved_He_all <- Zhang_He_stats_preserved[Zhang_He_stats_preserved$mod_name %in%  He_moduleTraitCor_Pval_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#1      blue       1000     15.979158
#2   magenta        588     12.212596
#3    yellow        941      9.095257
#4 turquoise       1000      8.963291

# He vs. Zhang
He_Zhang_multiExpr = list(He=list(data=He_dds_vst_matrix_common),Zhang=list(data=Zhang_dds_rlog_matrix_common))
He_Zhang_multiColor = list(He = He_moduleColors) 
#He_Zhang_mp =modulePreservation(He_Zhang_multiExpr , He_Zhang_multiColor,
#                                referenceNetworks=1,
#                                verbose=3,
#                                networkType="signed hybrid", # use same signed hybrid as before
#                                corFnc = "bicor", # use recommended bicor as before 
#                                nPermutations=100,
#                                randomSeed = 1, # recommended in langfelder tutorial
#                                quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
He_Zhang_stats = He_Zhang_mp$preservation$Z$ref.He$inColumnsAlsoPresentIn.Zhang
He_Zhang_stats_order <- He_Zhang_stats[order(-He_Zhang_stats[,2]),c(1:2)]
He_Zhang_stats_preserved <- He_Zhang_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
#mod_name moduleSize Zsummary.pres
#1    lightcyan        339     12.024891
#2      magenta        387     10.822159
#3         blue       1000      8.309949
#4       purple        386      7.213162
#5  floralwhite         71      7.018956
#6          tan        370      6.578664
#7         cyan        348      6.389133
#8      skyblue        163      5.914049
#9       salmon        365      5.075173
#10    skyblue3         93      5.072137
He_Zhang_stats_MR <- He_Zhang_mp$preservation$observed$ref.He$inColumnsAlsoPresentIn.Zhang
He_Zhang_stats_MR <- He_Zhang_stats_MR[,c("moduleSize","medianRank.pres")]
He_Zhang_stats_MR_order <- He_Zhang_stats_MR[order(-He_Zhang_stats_MR[,2]),]
He_Zhang_stats_MR_order_less20 <- He_Zhang_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save module preservation results
#save(He_Zhang_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/He_Zhang_mp_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(He_Zhang_stats_MR_order_less20, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

He_Zhang_stats_preserved_Zsum_medianRank <- He_Zhang_stats_MR_order_less20[He_Zhang_stats_MR_order_less20$mod_name %in% He_Zhang_stats_preserved$mod_name, ]

# Were any of these preserved modules significant in deLorg_Sus and He
He_Zhang_Zhang_LPS_preserved_all <- He_Zhang_stats_preserved[He_Zhang_stats_preserved$mod_name %in% Zhang_moduleTraitCor_Pval_LPS_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#2  magenta        387     10.822159
#3     blue       1000      8.309949
#6      tan        370      6.578664
#7     cyan        348      6.389133

He_Zhang_Zhang_Vibrio_preserved_all <- He_Zhang_stats_preserved[He_Zhang_stats_preserved$mod_name %in% Zhang_moduleTraitCor_Pval_Vibrio_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#3     blue       1000      8.309949
#4   purple        386      7.213162
He_Zhang_He_preserved_He_all <- He_Zhang_stats_preserved[He_Zhang_stats_preserved$mod_name %in%  He_moduleTraitCor_Pval_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#2  magenta        387     10.822159
#3     blue       1000      8.309949
#6      tan        370      6.578664
#8  skyblue        163      5.914049
#9   salmon        365      5.075173


## Assess whether Rubio and deLorg_Sus have conservation (both susceptible susceptible)
Rubio_deLorg_Sus_multiExpr = list(Rubio=list(data=Rubio_dds_rlog_matrix_common), deLorg_Sus=list(data=deLorgeril_Susceptible_dds_vst_matrix_common))
Rubio_deLorg_Sus_multiColor = list(Rubio = Rubio_moduleColors) 
#Rubio_deLorg_Sus_mp =modulePreservation(Rubio_deLorg_Sus_multiExpr, Rubio_deLorg_Sus_multiColor,
#                                referenceNetworks=1,
#                                verbose=3,
#                                networkType="signed hybrid", # use same signed hybrid as before
#                                corFnc = "bicor", # use recommended bicor as before 
#                                nPermutations=100,
#                                randomSeed = 1, # recommended in langfelder tutorial
#                                quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Rubio_deLorg_Sus_stats = Rubio_deLorg_Sus_mp$preservation$Z$ref.Rubio$inColumnsAlsoPresentIn.deLorg_Sus
Rubio_deLorg_Sus_stats_order <- Rubio_deLorg_Sus_stats[order(-Rubio_deLorg_Sus_stats[,2]),c(1:2)]
Rubio_deLorg_Sus_stats_preserved <- Rubio_deLorg_Sus_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
#mod_name moduleSize Zsummary.pres
#1    greenyellow        857     34.234560
#2          brown       1000     27.352964
#3      turquoise       1000     18.251003
#4           cyan        527     15.659252
#5            tan        642     14.426980
#6  darkturquoise        249     11.512318
#7           blue       1000     10.752557
#8    lightyellow        322     10.424208
#9           gold       1000      9.891569
#10         green       1000      9.854796
#11           red       1000      9.823248
#12      darkgrey        171      6.377679
#13         black       1000      6.203096
#14        salmon        633      5.168682

Rubio_deLorg_Sus_stats_MR <- Rubio_deLorg_Sus_mp$preservation$observed$ref.Rubio$inColumnsAlsoPresentIn.deLorg_Sus
Rubio_deLorg_Sus_stats_MR <- Rubio_deLorg_Sus_stats_MR [,c("moduleSize","medianRank.pres")]
Rubio_deLorg_Sus_stats_MR_order <- Rubio_deLorg_Sus_stats_MR [order(-Rubio_deLorg_Sus_stats_MR [,2]),]
Rubio_deLorg_Sus_stats_MR_order_less20 <- Rubio_deLorg_Sus_stats_MR_order  %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save module preservation results
#save(Rubio_deLorg_Sus_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Rubio_deLorg_Sus_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(Rubio_deLorg_Sus_stats_MR_order_less20, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

Rubio_deLorg_Sus_stats_preserved_Zsum_medianRank <- Rubio_deLorg_Sus_stats_MR_order_less20[Rubio_deLorg_Sus_stats_MR_order_less20$mod_name %in% Rubio_deLorg_Sus_stats_preserved $mod_name, ]

# Were any of these preserved modules significant in deLorg_Sus and He
Rubio_deLorg_Sus_delorg_Sus_preserved_all <- Rubio_deLorg_Sus_stats_preserved [Rubio_deLorg_Sus_stats_preserved $mod_name %in% deLorg_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#3 turquoise       1000      18.25100
#7      blue       1000      10.75256
Rubio_deLorg_Sus_Rubio_NV_preserved_He_all<- Rubio_deLorg_Sus_stats_preserved [Rubio_deLorg_Sus_stats_preserved $mod_name %in%  Rubio_moduleTraitCor_Pval_NV_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#2          brown       1000     27.352964
#3      turquoise       1000     18.251003
#4           cyan        527     15.659252
#6  darkturquoise        249     11.512318
#10         green       1000      9.854796
#11           red       1000      9.823248

Rubio_deLorg_Sus_Rubio_V_preserved_He_all<- Rubio_deLorg_Sus_stats_preserved [Rubio_deLorg_Sus_stats_preserved $mod_name %in%  Rubio_moduleTraitCor_Pval_V_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#2          brown       1000     27.352964
#3      turquoise       1000     18.251003
#4           cyan        527     15.659252
#6  darkturquoise        249     11.512318
#7           blue       1000     10.752557
#10         green       1000      9.854796
#11           red       1000      9.823248

# DeLorg Sus vs. Rubio 
deLorg_Sus_Rubio_multiExpr = list(deLorg_Sus=list(data=deLorgeril_Susceptible_dds_vst_matrix_common), Rubio=list(data=Rubio_dds_rlog_matrix_common))
deLorg_Sus_Rubio_multiColor = list(deLorg_Sus = deLorg_Sus_moduleColors) 
#deLorg_Sus_Rubio_mp = modulePreservation(deLorg_Sus_Rubio_multiExpr, deLorg_Sus_Rubio_multiColor,
#                                referenceNetworks=1,
#                                verbose=3,
#                                networkType="signed hybrid", # use same signed hybrid as before
#                                corFnc = "bicor", # use recommended bicor as before 
#                                nPermutations=100,
#                                randomSeed = 1, # recommended in langfelder tutorial
#                                quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
deLorg_Sus_Rubio_stats = deLorg_Sus_Rubio_mp$preservation$Z$ref.deLorg_Sus$inColumnsAlsoPresentIn.Rubio
deLorg_Sus_Rubio_stats_order <- deLorg_Sus_Rubio_stats[order(-deLorg_Sus_Rubio_stats[,2]),c(1:2)]
deLorg_Sus_Rubio_stats_preserved <- deLorg_Sus_Rubio_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
#mod_name moduleSize Zsummary.pres
#1       red        526     17.918318
#2    yellow       1000     13.070538
#3     green       1000     12.151436
#4     brown       1000     11.733508
#5      blue       1000      9.901168
#6      pink        406      9.555296
#7 turquoise       1000      8.585884
#8   magenta        356      8.546037
#9      gold       1000      8.315236

deLorg_Sus_Rubio_stats_MR <- deLorg_Sus_Rubio_mp$preservation$observed$ref.deLorg_Sus$inColumnsAlsoPresentIn.Rubio
deLorg_Sus_Rubio_stats_MR <- deLorg_Sus_Rubio_stats_MR[,c("moduleSize","medianRank.pres")]
deLorg_Sus_Rubio_stats_MR_order <- deLorg_Sus_Rubio_stats_MR[order(-deLorg_Sus_Rubio_stats_MR[,2]),]
deLorg_Sus_Rubio_stats_MR_order_less20 <- deLorg_Sus_Rubio_stats_MR_order  %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save module preservation results
# save(deLorg_Sus_Rubio_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/deLorg_Sus_Rubio_mp_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(deLorg_Sus_Rubio_stats_MR_order_less20, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

deLorg_Sus_Rubio_stats_preserved_Zsum_medianRank <-deLorg_Sus_Rubio_stats_MR_order_less20[deLorg_Sus_Rubio_stats_MR_order_less20$mod_name %in% deLorg_Sus_Rubio_stats_preserved$mod_name, ]

# Were any of these preserved modules significant in deLorg_Sus and He
deLorg_Sus_Rubio_delorg_Sus_preserved_all <- deLorg_Sus_Rubio_stats_preserved [deLorg_Sus_Rubio_stats_preserved$mod_name %in% deLorg_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
# mod_name moduleSize Zsummary.pres
# 5      blue       1000      9.901168
# 7 turquoise       1000      8.585884
deLorg_Sus_Rubio_Rubio_NV_preserved_He_all<- deLorg_Sus_Rubio_stats_preserved [deLorg_Sus_Rubio_stats_preserved$mod_name %in%  Rubio_moduleTraitCor_Pval_NV_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#1       red        526     17.918318
#3     green       1000     12.151436
#4     brown       1000     11.733508
#7 turquoise       1000      8.585884

deLorg_Sus_Rubio_Rubio_V_preserved_He_all<-  deLorg_Sus_Rubio_stats_preserved [deLorg_Sus_Rubio_stats_preserved$mod_name %in%  Rubio_moduleTraitCor_Pval_V_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#1       red        526     17.918318
#2    yellow       1000     13.070538
#3     green       1000     12.151436
#4     brown       1000     11.733508
#5      blue       1000      9.901168
#7 turquoise       1000      8.585884

## DeLorg Tol vs. Rubio 
deLorg_Tol_Rubio_multiExpr = list(deLorg_Tol=list(data=deLorgeril_Tolerant_dds_vst_matrix_common), Rubio=list(data=Rubio_dds_rlog_matrix_common))
deLorg_Tol_Rubio_multiColor = list(deLorg_Tol = deLorg_Tol_moduleColors) 
deLorg_Tol_Rubio_mp = modulePreservation(deLorg_Tol_Rubio_multiExpr, deLorg_Tol_Rubio_multiColor,
                                referenceNetworks=1,
                                verbose=3,
                                networkType="signed hybrid", # use same signed hybrid as before
                                corFnc = "bicor", # use recommended bicor as before 
                                nPermutations=100,
                                randomSeed = 1, # recommended in langfelder tutorial
                                quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
deLorg_Tol_Rubio_stats = deLorg_Tol_Rubio_mp$preservation$Z$ref.deLorg_Tol$inColumnsAlsoPresentIn.Rubio
deLorg_Tol_Rubio_stats_order <- deLorg_Tol_Rubio_stats[order(-deLorg_Tol_Rubio_stats[,2]),c(1:2)]
deLorg_Tol_Rubio_stats_preserved <- deLorg_Tol_Rubio_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
deLorg_Tol_Rubio_stats_MR <- deLorg_Tol_Rubio_mp$preservation$observed$ref.deLorg_Tol$inColumnsAlsoPresentIn.Rubio
deLorg_Tol_Rubio_stats_MR <- deLorg_Tol_Rubio_stats_MR[,c("moduleSize","medianRank.pres")]
deLorg_Tol_Rubio_stats_MR_order <- deLorg_Tol_Rubio_stats_MR[order(-deLorg_Tol_Rubio_stats_MR[,2]),]
deLorg_Tol_Rubio_stats_MR_order_less20 <- deLorg_Tol_Rubio_stats_MR_order  %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save module preservation results
save(deLorg_Tol_Rubio_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/deLorg_Tol_Rubio_mp_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(deLorg_Tol_Rubio_stats_MR_order_less20, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

deLorg_Tol_Rubio_stats_preserved_Zsum_medianRank <-deLorg_Tol_Rubio_stats_MR_order_less20[deLorg_Tol_Rubio_stats_MR_order_less20$mod_name %in% deLorg_Tol_Rubio_stats_preserved$mod_name, ]

# Were any of these preserved modules significant in deLorg_Sus and He
deLorg_Tol_Rubio_delorg_Tol_preserved_all <- deLorg_Tol_Rubio_stats_preserved [deLorg_Tol_Rubio_stats_preserved$mod_name %in% deLorg_Tol_moduleTraitCor_Pval_df_sig_list_rm,]
deLorg_Tol_Rubio_Rubio_NV_preserved_He_all<- deLorg_Tol_Rubio_stats_preserved [deLorg_Tol_Rubio_stats_preserved$mod_name %in%  Rubio_moduleTraitCor_Pval_NV_df_sig_list_rm,]
deLorg_Tol_Rubio_Rubio_V_preserved_He_all<-  deLorg_Tol_Rubio_stats_preserved [deLorg_Tol_Rubio_stats_preserved$mod_name %in%  Rubio_moduleTraitCor_Pval_V_df_sig_list_rm,]

## Rubio vs. deLorg Tol
Rubio_deLorg_Tol_multiExpr = list(Rubio=list(data=Rubio_dds_rlog_matrix_common), deLorg_Tol=list(data=deLorgeril_Tolerant_dds_vst_matrix_common))
Rubio_deLorg_Tol_multiColor = list(Rubio = Rubio_moduleColors) 
Rubio_deLorg_Tol_mp = modulePreservation(Rubio_deLorg_Tol_multiExpr, Rubio_deLorg_Tol_multiColor,
                                         referenceNetworks=1,
                                         verbose=3,
                                         networkType="signed hybrid", # use same signed hybrid as before
                                         corFnc = "bicor", # use recommended bicor as before 
                                         nPermutations=100,
                                         randomSeed = 1, # recommended in langfelder tutorial
                                         quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Rubio_deLorg_Tol_stats = Rubio_deLorg_Tol_mp$preservation$Z$ref.Rubio$inColumnsAlsoPresentIn.deLorg_Tol
Rubio_deLorg_Tol_stats_order <- Rubio_deLorg_Tol_stats[order(-Rubio_deLorg_Tol_stats[,2]),c(1:2)]
Rubio_deLorg_Tol_stats_preserved <- Rubio_deLorg_Tol_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
Rubio_deLorg_Tol_stats_MR <- Rubio_deLorg_Tol_mp$preservation$observed$ref.Rubio$inColumnsAlsoPresentIn.deLorg_Tol
Rubio_deLorg_Tol_stats_MR <- Rubio_deLorg_Tol_stats_MR[,c("moduleSize","medianRank.pres")]
Rubio_deLorg_Tol_stats_MR_order <- Rubio_deLorg_Tol_stats_MR[order(-Rubio_deLorg_Tol_stats_MR[,2]),]
Rubio_deLorg_Tol_stats_MR_order_less20 <- Rubio_deLorg_Tol_stats_MR_order  %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save module preservation results
save(Rubio_deLorg_Tol_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Rubio_deLorg_Tol_mp_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(Rubio_deLorg_Tol_stats_MR_order_less20, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

Rubio_deLorg_Tol_stats_preserved_Zsum_medianRank <- Rubio_deLorg_Tol_stats_MR_order_less20[Rubio_deLorg_Tol_stats_MR_order_less20$mod_name %in% Rubio_deLorg_Tol_stats_preserved$mod_name, ]

# Were any of these preserved modules significant in deLorg_Sus and He
Rubio_deLorg_Tol_delorg_Tol_preserved_all <- Rubio_deLorg_Tol_stats_preserved [Rubio_deLorg_Tol_stats_preserved$mod_name %in% deLorg_Tol_moduleTraitCor_Pval_df_sig_list_rm,]
Rubio_deLorg_Tol_Rubio_NV_preserved_He_all<- Rubio_deLorg_Tol_stats_preserved [Rubio_deLorg_Tol_stats_preserved$mod_name %in%  Rubio_moduleTraitCor_Pval_NV_df_sig_list_rm,]
Rubio_deLorg_Tol_Rubio_V_preserved_He_all<-  Rubio_deLorg_Tol_stats_preserved [Rubio_deLorg_Tol_stats_preserved$mod_name %in%  Rubio_moduleTraitCor_Pval_V_df_sig_list_rm,]

# Assess whether Zhang and Rubio modules are preserved
Zhang_Rubio_multiExpr = list(Zhang=list(data=Zhang_dds_rlog_matrix_common),Rubio=list(data=Rubio_dds_rlog_matrix_common))
Zhang_multiColor = list(Zhang = Zhang_moduleColors) 
#Zhang_Rubio_mp =modulePreservation(Zhang_Rubio_multiExpr , Zhang_multiColor,
#                                   referenceNetworks=1,
#                                   verbose=3,
#                                   networkType="signed hybrid", # use same signed hybrid as before
#                                   corFnc = "bicor", # use recommended bicor as before 
#                                   nPermutations=100, #need to change to 100
#                                   randomSeed = 1, # recommended in langfelder tutorial
#                                   quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Zhang_Rubio_stats = Zhang_Rubio_mp$preservation$Z$ref.Zhang$inColumnsAlsoPresentIn.Rubio
Zhang_Rubio_stats_order <- Zhang_Rubio_stats[order(-Zhang_Rubio_stats[,2]),c(1:2)]
Zhang_Rubio_stats_preserved <- Zhang_Rubio_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
# mod_name moduleSize Zsummary.pres
# 1     turquoise       1000     17.390204
# 2          gold       1000      9.128892 # dont count
# 3        yellow        941      8.219197
# 4          grey       1000      6.899878 # dont count
# 5         green        849      6.420218
# 6           red        751      5.746211
# 7 paleturquoise        219      5.481892
# 8     darkgreen        270      5.241935
# 9        purple        531      5.045024
Zhang_Rubio_stats_MR <- Zhang_Rubio_mp$preservation$observed$ref.Zhang$inColumnsAlsoPresentIn.Rubio
Zhang_Rubio_stats_MR <- Zhang_Rubio_stats_MR [,c("moduleSize","medianRank.pres")]
Zhang_Rubio_stats_MR_order <- Zhang_Rubio_stats_MR[order(Zhang_Rubio_stats_MR[,2]),]
Zhang_Rubio_stats_MR_order_less20 <- Zhang_Rubio_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save module preservation results
#save(Zhang_Rubio_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Zhang_Rubio_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(Zhang_Rubio_stats_MR_order_less20, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

Zhang_Rubio_stats_preserved_Zsum_medianRank <- Zhang_Rubio_stats_MR_order_less20[Zhang_Rubio_stats_MR_order_less20$mod_name %in% Zhang_Rubio_stats_preserved,]

# Were any of these preserved modules significant in Zhang or Rubio
Zhang_LPS_preserved_all <- Zhang_Rubio_stats_preserved[Zhang_Rubio_stats_preserved$mod_name %in% Zhang_moduleTraitCor_Pval_LPS_df_sig_list_rm,]
#       mod_name moduleSize Zsummary.pres
# 1     turquoise       1000     17.390204
# 7 paleturquoise        219      5.481892
Zhang_Vibrio_preserved_all <- Zhang_Rubio_stats_preserved[Zhang_Rubio_stats_preserved$mod_name %in%  Zhang_moduleTraitCor_Pval_Vibrio_df_sig_list_rm,]
#mod_name moduleSize Zsummary.pres
#8 darkgreen        270      5.434683
Rubio_NV_preserved_all <- Zhang_Rubio_stats_preserved[Zhang_Rubio_stats_preserved$mod_name %in%  Rubio_moduleTraitCor_Pval_NV_df_sig_list_rm,]
# mod_name moduleSize Zsummary.pres
# 1     turquoise       1000     17.660194
# 4         green        849      6.617005
# 6           red        751      5.755566
# 7 paleturquoise        219      5.601135
Rubio_V_preserved_all <- Zhang_Rubio_stats_preserved[Zhang_Rubio_stats_preserved$mod_name %in%  Rubio_moduleTraitCor_Pval_V_df_sig_list_rm,]
# mod_name moduleSize Zsummary.pres
# 1 turquoise       1000     17.660194
# 3    yellow        941      8.639548
# 4     green        849      6.617005
# 6       red        751      5.755566

# Rubio vs. Zhang 
Rubio_Zhang_multiExpr = list(Rubio=list(data=Rubio_dds_rlog_matrix_common), Zhang=list(data=Zhang_dds_rlog_matrix_common))
Rubio_multiColor = list(Rubio = Rubio_moduleColors) 
Rubio_Zhang_mp =modulePreservation(Rubio_Zhang_multiExpr , Rubio_multiColor,
                                   referenceNetworks=1,
                                   verbose=3,
                                   networkType="signed hybrid", # use same signed hybrid as before
                                   corFnc = "bicor", # use recommended bicor as before 
                                   nPermutations=100, #need to change to 100
                                   randomSeed = 1, # recommended in langfelder tutorial
                                   quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Rubio_Zhang_stats = Rubio_Zhang_mp$preservation$Z$ref.Rubio$inColumnsAlsoPresentIn.Zhang
Rubio_Zhang_stats_order <- Rubio_Zhang_stats[order(-Rubio_Zhang_stats[,2]),c(1:2)]
Rubio_Zhang_stats_preserved <- Rubio_Zhang_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
Rubio_Zhang_stats_MR <- Rubio_Zhang_mp$preservation$observed$ref.Rubio$inColumnsAlsoPresentIn.Zhang
Rubio_Zhang_stats_MR <- Rubio_Zhang_stats_MR[,c("moduleSize","medianRank.pres")]
Rubio_Zhang_stats_MR_order <- Rubio_Zhang_stats_MR[order(Rubio_Zhang_stats_MR[,2]),]
Rubio_Zhang_stats_MR_order_less20 <- Rubio_Zhang_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save module preservation results
save(Rubio_Zhang_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Rubio_Zhang_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(Rubio_Zhang_stats_MR_order_less20, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

Rubio_Zhang_stats_preserved_Zsum_medianRank <- Rubio_Zhang_stats_MR_order_less20[Rubio_Zhang_stats_MR_order_less20$mod_name %in% Rubio_Zhang_stats_preserved,]

# Were any of these preserved modules significant in Zhang or Rubio
Rubio_Zhang_Zhang_LPS_preserved_all <- Rubio_Zhang_stats_preserved[Rubio_Zhang_stats_preserved$mod_name %in% Zhang_moduleTraitCor_Pval_LPS_df_sig_list_rm,]
Rubio_Zhang_Zhang_Vibrio_preserved_all <- Rubio_Zhang_stats_preserved[Rubio_Zhang_stats_preserved$mod_name %in%  Zhang_moduleTraitCor_Pval_Vibrio_df_sig_list_rm,]
Rubio_Zhang_Rubio_NV_preserved_all <- Rubio_Zhang_stats_preserved[Rubio_Zhang_stats_preserved$mod_name %in%  Rubio_moduleTraitCor_Pval_NV_df_sig_list_rm,]
Rubio_Zhang_Rubio_V_preserved_all <- Rubio_Zhang_stats_preserved[Rubio_Zhang_stats_preserved$mod_name %in%  Rubio_moduleTraitCor_Pval_V_df_sig_list_rm,]

#### PERFORM CONSENSUS NETWORK ANALYSIS ACROSS C.GIGAS NETWORKS BACTERIA ####

# haven't performed yet for bacteria and virus 


#### Running all CGIGAS "FULL" TRANSCRIPT SETS" ####

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function, following general recommendations to set network type to "signed hybrid" and using the "bicor" correlation
#Zhang_dds_rlog_matrix_sft <- pickSoftThreshold(Zhang_dds_rlog_matrix, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
#Rubio_dds_rlog_matrix_sft <- pickSoftThreshold(Rubio_dds_rlog_matrix, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
#deLorgeril_Resistant_dds_vst_matrix_sft <- pickSoftThreshold(deLorgeril_Resistant_dds_vst_matrix, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
#deLorgeril_Susceptible_dds_vst_matrix_sft <- pickSoftThreshold(deLorgeril_Susceptible_dds_vst_matrix, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
#He_dds_vst_matrix_sft <- pickSoftThreshold(He_dds_vst_matrix, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
##save soft thresholding results
#save(Zhang_dds_rlog_matrix_sft,Rubio_dds_rlog_matrix_sft,deLorgeril_Resistant_dds_vst_matrix_sft,deLorgeril_Susceptible_dds_vst_matrix_sft,
#     He_dds_vst_matrix_sft, file= "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gigas_individual_FULL_pickSoftThreshold_WGCNA_input.RData")

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Zhang 
# Scale-free topology fit index as a function of the soft-thresholding power
plot(Zhang_dds_rlog_matrix_sft$fitIndices[,1], -sign(Zhang_dds_rlog_matrix_sft$fitIndices[,3])*Zhang_dds_rlog_matrix_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(Zhang_dds_rlog_matrix_sft$fitIndices[,1], -sign(Zhang_dds_rlog_matrix_sft$fitIndices[,3])*Zhang_dds_rlog_matrix_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(Zhang_dds_rlog_matrix_sft$fitIndices[,1], Zhang_dds_rlog_matrix_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(Zhang_dds_rlog_matrix_sft$fitIndices[,1], Zhang_dds_rlog_matrix_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Softthreshold of 4since this is lowest value past 0.9 we start to see flattening 


# Rubio
# Scale-free topology fit index as a function of the soft-thresholding power
plot(Rubio_dds_rlog_matrix_sft$fitIndices[,1], -sign(Rubio_dds_rlog_matrix_sft$fitIndices[,3])*Rubio_dds_rlog_matrix_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(Rubio_dds_rlog_matrix_sft$fitIndices[,1], -sign(Rubio_dds_rlog_matrix_sft$fitIndices[,3])*Rubio_dds_rlog_matrix_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(Rubio_dds_rlog_matrix_sft$fitIndices[,1], Rubio_dds_rlog_matrix_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(Rubio_dds_rlog_matrix_sft$fitIndices[,1], Rubio_dds_rlog_matrix_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Softthreshold of 7 since this is lowest value past 0.9 we start to see flattening 

# deLorg Res
# Scale-free topology fit index as a function of the soft-thresholding power
plot(deLorgeril_Resistant_dds_vst_matrix_sft$fitIndices[,1], -sign(deLorgeril_Resistant_dds_vst_matrix_sft$fitIndices[,3])*deLorgeril_Resistant_dds_vst_matrix_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(deLorgeril_Resistant_dds_vst_matrix_sft$fitIndices[,1], -sign(deLorgeril_Resistant_dds_vst_matrix_sft$fitIndices[,3])*deLorgeril_Resistant_dds_vst_matrix_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(deLorgeril_Resistant_dds_vst_matrix_sft$fitIndices[,1], deLorgeril_Resistant_dds_vst_matrix_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(deLorgeril_Resistant_dds_vst_matrix_sft$fitIndices[,1], deLorgeril_Resistant_dds_vst_matrix_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# doesnt quite fit scale free topology now, use 8 since there are 21 samples

# deLorg Sus = 21 samples
# Scale-free topology fit index as a function of the soft-thresholding power
plot(deLorgeril_Susceptible_dds_vst_matrix_sft$fitIndices[,1], -sign(deLorgeril_Susceptible_dds_vst_matrix_sft$fitIndices[,3])*deLorgeril_Susceptible_dds_vst_matrix_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(deLorgeril_Susceptible_dds_vst_matrix_sft$fitIndices[,1], -sign(deLorgeril_Susceptible_dds_vst_matrix_sft$fitIndices[,3])*deLorgeril_Susceptible_dds_vst_matrix_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(deLorgeril_Susceptible_dds_vst_matrix_sft$fitIndices[,1], deLorgeril_Susceptible_dds_vst_matrix_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(deLorgeril_Susceptible_dds_vst_matrix_sft$fitIndices[,1], deLorgeril_Susceptible_dds_vst_matrix_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Softthreshold of 3 since this is lowest value past 0.9 we start to see flattening and 21 samples

# He
# Scale-free topology fit index as a function of the soft-thresholding power
plot(He_dds_vst_matrix_sft$fitIndices[,1], -sign(He_dds_vst_matrix_sft$fitIndices[,3])*He_dds_vst_matrix_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(He_dds_vst_matrix_sft$fitIndices[,1], -sign(He_dds_vst_matrix_sft$fitIndices[,3])*He_dds_vst_matrix_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(He_dds_vst_matrix_sft$fitIndices[,1], He_dds_vst_matrix_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(He_dds_vst_matrix_sft$fitIndices[,1], He_dds_vst_matrix_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Softthreshold of 5 since this is lowest value past 0.9 we start to see flattening 

## ONE STEP NETWORK CONSTRUCTION, MODULE DETECTION, MODULE DENDROGRAM INSPECTION ##

# Zhang 
#Zhang_full_net = blockwiseModules(Zhang_dds_rlog_matrix, power = 4, # picked suitable power in the code above 
#                                 TOMType = "signed", # use signed TOM type
#                                 networkType= "signed hybrid", # use signed hybrid network type
#                                 corType = "bicor", # use suggested bicor
#                                 TminModuleSize = 30, # recommended default
#                                 reassignThreshold = 0, # recommended default
#                                 mergeCutHeight = 0.25, # recommended default
#                                 numericLabels = TRUE, # recommended default
#                                 pamRespectsDendro = FALSE,# recommended default
#                                 verbose = 3, 
#                                 maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
## How many modules identified
#table(Zhang_full_net$colors)
## Plot dendrogram with colors
## open a graphics window
#sizeGrWindow(12, 9)
## Convert labels to colors for plotting
#Zhang_full_mergedColors = labels2colors(Zhang_full_net$colors)
## Plot the dendrogram and the module colors underneath
#plotDendroAndColors(Zhang_full_net$dendrograms[[1]], Zhang_full_mergedColors[Zhang_full_net$blockGenes[[1]]],
#                    "Module colors",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05)
#Zhang_full_moduleLabels = Zhang_full_net$colors
#Zhang_full_moduleColors = labels2colors(Zhang_full_net$colors)
#Zhang_full_MEs = Zhang_full_net$MEs
#Zhang_full_geneTree = Zhang_full_net$dendrograms[[1]]
##
### Save Zhang full network
#save(Zhang_full_net, Zhang_full_mergedColors, Zhang_full_moduleLabels, Zhang_full_moduleColors, Zhang_full_MEs,Zhang_full_geneTree, 
#     file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Zhang_full_network.RData")
#
## Rubio
# Rubio_full_net = blockwiseModules(Rubio_dds_rlog_matrix, power = 7, # picked suitable power in the code above 
#                              TOMType = "signed", # use signed TOM type
#                              networkType= "signed hybrid", # use signed hybrid network type
#                              corType = "bicor", # use suggested bicor
#                              TminModuleSize = 30, # recommended default
#                              reassignThreshold = 0, # recommended default
#                              mergeCutHeight = 0.25, # recommended default
#                              numericLabels = TRUE, # recommended default
#                              pamRespectsDendro = FALSE,# recommended default
#                              verbose = 3, 
#                              maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
## # How many modules identified
#table(Rubio_full_net$colors)
## Plot dendrogram with colors
## open a graphics window
#sizeGrWindow(12, 9)
## Convert labels to colors for plotting
#Rubio_full_mergedColors = labels2colors(Rubio_full_net$colors)
## Plot the dendrogram and the module colors underneath
#plotDendroAndColors(Rubio_full_net$dendrograms[[1]], Rubio_full_mergedColors[Rubio_full_net$blockGenes[[1]]],
#                    "Module colors",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05)
#Rubio_full_moduleLabels = Rubio_full_net$colors
#Rubio_full_moduleColors = labels2colors(Rubio_full_net$colors)
#Rubio_full_MEs = Rubio_full_net$MEs
#Rubio_full_geneTree = Rubio_full_net$dendrograms[[1]]
#
## save network
#save(Rubio_full_net, Rubio_full_mergedColors, Rubio_full_moduleLabels, Rubio_full_moduleColors, Rubio_full_MEs,Rubio_full_geneTree, 
#         file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Rubio_full_network.RData")
#
#
## deLorg Res
#deLorg_Res_full_net = blockwiseModules(deLorgeril_Resistant_dds_vst_matrix, power = 8, # picked suitable power in the code above 
#                             TOMType = "signed", # use signed TOM type
#                             networkType= "signed hybrid", # use signed hybrid network type
#                             corType = "bicor", # use suggested bicor
#                             TminModuleSize = 30, # recommended default
#                             reassignThreshold = 0, # recommended default
#                             mergeCutHeight = 0.25, # recommended default
#                             numericLabels = TRUE, # recommended default
#                             pamRespectsDendro = FALSE,# recommended default
#                             verbose = 3, 
#                             maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
## How many modules identified
#table(deLorg_Res_full_net$colors)
## Plot dendrogram with colors
## open a graphics window
#sizeGrWindow(12, 9)
## Convert labels to colors for plotting
#deLorg_Res_full_mergedColors = labels2colors(deLorg_Res_full_net$colors)
## Plot the dendrogram and the module colors underneath
#plotDendroAndColors(deLorg_Res_full_net$dendrograms[[1]], deLorg_Res_full_mergedColors[deLorg_Res_full_net$blockGenes[[1]]],
#                    "Module colors",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05)
#deLorg_Res_full_moduleLabels = deLorg_Res_full_net$colors
#deLorg_Res_full_moduleColors = labels2colors(deLorg_Res_full_net$colors)
#deLorg_Res_full_MEs = deLorg_Res_full_net$MEs
#deLorg_Res_full_geneTree = deLorg_Res_full_net$dendrograms[[1]]
#
#save(deLorg_Res_full_net, deLorg_Res_full_mergedColors,deLorg_Res_full_moduleLabels, deLorg_Res_full_moduleColors, deLorg_Res_full_MEs,deLorg_Res_full_geneTree, 
#     file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/deLorg_Res_full_network.RData")
#
## deLorg Sus
#deLorg_Sus_full_net = blockwiseModules(deLorgeril_Susceptible_dds_vst_matrix, power = 3, # picked suitable power in the code above 
#                                  TOMType = "signed", # use signed TOM type
#                                  networkType= "signed hybrid", # use signed hybrid network type
#                                  corType = "bicor", # use suggested bicor
#                                  TminModuleSize = 30, # recommended default
#                                  reassignThreshold = 0, # recommended default
#                                  mergeCutHeight = 0.25, # recommended default
#                                  numericLabels = TRUE, # recommended default
#                                  pamRespectsDendro = FALSE,# recommended default
#                                  verbose = 3, 
#                                  maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
## How many modules identified
#table(deLorg_Sus_full_net$colors)
## Plot dendrogram with colors
## open a graphics window
#sizeGrWindow(12, 9)
## Convert labels to colors for plotting
#deLorg_Sus_full_mergedColors = labels2colors(deLorg_Sus_full_net$colors)
## Plot the dendrogram and the module colors underneath
#plotDendroAndColors(deLorg_Sus_full_net$dendrograms[[1]], deLorg_Sus_full_mergedColors[deLorg_Sus_full_net$blockGenes[[1]]],
#                    "Module colors",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05)
#deLorg_Sus_full_moduleLabels = deLorg_Sus_full_net$colors
#deLorg_Sus_full_moduleColors = labels2colors(deLorg_Sus_full_net$colors)
#deLorg_Sus_full_MEs = deLorg_Sus_full_net$MEs
#deLorg_Sus_full_geneTree = deLorg_Sus_full_net$dendrograms[[1]]
#
#save(deLorg_Sus_full_net, deLorg_Sus_full_mergedColors,deLorg_Sus_full_moduleLabels, deLorg_Sus_full_moduleColors, deLorg_Sus_full_MEs,deLorg_Sus_full_geneTree, 
#     file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/deLorg_Sus_full_network.RData")
#
## He 
#He_full_net = blockwiseModules(He_dds_vst_matrix, power = 5, # picked suitable power in the code above 
#                                  TOMType = "signed", # use signed TOM type
#                                  networkType= "signed hybrid", # use signed hybrid network type
#                                  corType = "bicor", # use suggested bicor
#                                  TminModuleSize = 30, # recommended default
#                                  reassignThreshold = 0, # recommended default
#                                  mergeCutHeight = 0.25, # recommended default
#                                  numericLabels = TRUE, # recommended default
#                                  pamRespectsDendro = FALSE,# recommended default
#                                  verbose = 3, 
#                                  maxBlockSize = 20000) # 20,000 should be okay because I have 16GB memory on my computer
## How many modules identified
#table(He_full_net$colors)
## Plot dendrogram with colors
## open a graphics window
#sizeGrWindow(12, 9)
## Convert labels to colors for plotting
#He_full_mergedColors = labels2colors(He_full_net$colors)
## Plot the dendrogram and the module colors underneath
#plotDendroAndColors(He_full_net$dendrograms[[1]], He_full_mergedColors[He_full_net$blockGenes[[1]]],
#                    "Module colors",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05)
#He_full_moduleLabels = He_full_net$colors
#He_full_moduleColors = labels2colors(He_full_net$colors)
#He_full_MEs = He_full_net$MEs
#He_full_geneTree = He_full_net$dendrograms[[1]]
#
#save(He_full_net, He_full_mergedColors,He_full_moduleLabels, He_full_moduleColors, He_full_MEs,He_full_geneTree, 
#     file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/He_full_network.RData")

### QUANTIFYING MODULE TRAIT ASSOCIATIONS ###

### Zhang ###
# Define numbers of genes and samples
Zhang_full_nGenes = ncol(Zhang_dds_rlog_matrix)
Zhang_full_nSamples = nrow(Zhang_dds_rlog_matrix)

# Recalculate MEs with color labels
Zhang_full_MEs0 = moduleEigengenes(Zhang_dds_rlog_matrix, Zhang_full_moduleColors)$eigengenes
Zhang_full_MEs = orderMEs(Zhang_MEs0)
Zhang_full_moduleTraitCor = cor(Zhang_full_MEs, Zhang_coldata_collapsed_binarize, use = "p");
Zhang_full_moduleTraitPvalue = corPvalueStudent(Zhang_full_moduleTraitCor, Zhang_full_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trai
sizeGrWindow(10,6)
# Will display correlations and their p-values
Zhang_full_textMatrix = paste(signif(Zhang_full_moduleTraitCor, 2), "\n(",
                         signif(Zhang_full_moduleTraitPvalue, 1), ")", sep = "");
dim(Zhang_full_textMatrix) = dim(Zhang_full_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = Zhang_full_moduleTraitCor,
               xLabels = names(Zhang_coldata_collapsed_binarize),
               yLabels =  names(Zhang_full_MEs),
               ySymbols = names(Zhang_full_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = Zhang_full_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
Zhang_full_moduleTraitCor_df <- as.data.frame(Zhang_full_moduleTraitCor)
Zhang_full_moduleTraitCor_df$mod_names <- row.names(Zhang_full_moduleTraitCor_df)
Zhang_full_moduleTraitCor_LPS_df <- Zhang_full_moduleTraitCor_df[,c("mod_names","group_by_sim.LPS_M_lut.vs.control")]
Zhang_full_moduleTraitPvalue_df <- as.data.frame(Zhang_full_moduleTraitPvalue)
Zhang_full_moduleTraitPvalue_df$mod_names <- row.names(Zhang_full_moduleTraitPvalue_df )
Zhang_full_moduleTraitPvalue_LPS_df <- Zhang_full_moduleTraitPvalue_df[,c("mod_names","group_by_sim.LPS_M_lut.vs.control")]
colnames(Zhang_full_moduleTraitPvalue_LPS_df)[2] <- "pvalue"
Zhang_full_moduleTraitCor_Pval_LPS_df <- join(Zhang_full_moduleTraitCor_LPS_df, Zhang_full_moduleTraitPvalue_LPS_df, by = "mod_names")

# Significantly correlated modules
Zhang_full_moduleTraitCor_Pval_LPS_df[order(Zhang_full_moduleTraitCor_Pval_LPS_df$pvalue),]
class(Zhang_full_moduleTraitCor_Pval_LPS_df$pvalue) # numeric
Zhang_full_moduleTraitCor_Pval_LPS_df_sig <- Zhang_full_moduleTraitCor_Pval_LPS_df %>% filter(pvalue <= 0.05)
Zhang_full_moduleTraitCor_Pval_LPS_df_sig # 47

Zhang_full_moduleTraitCor_df <- as.data.frame(Zhang_full_moduleTraitCor)
Zhang_full_moduleTraitCor_df$mod_names <- row.names(Zhang_full_moduleTraitCor_df)
Zhang_full_moduleTraitCor_Vibrio_df <- Zhang_full_moduleTraitCor_df[,c("mod_names","group_by_sim.Vibrio.vs.control")]
Zhang_full_moduleTraitPvalue_Vibrio_df <- Zhang_full_moduleTraitPvalue_df[,c("mod_names","group_by_sim.Vibrio.vs.control")]
colnames(Zhang_full_moduleTraitPvalue_Vibrio_df)[2] <- "pvalue"
Zhang_full_moduleTraitCor_Pval_Vibrio_df <- join(Zhang_full_moduleTraitCor_Vibrio_df, Zhang_full_moduleTraitPvalue_Vibrio_df, by = "mod_names")

# Significantly correlated modules
Zhang_full_moduleTraitCor_Pval_Vibrio_df[order(Zhang_full_moduleTraitCor_Pval_Vibrio_df$pvalue),]
class(Zhang_full_moduleTraitCor_Pval_Vibrio_df$pvalue) # numeric
Zhang_full_moduleTraitCor_Pval_Vibrio_df_sig <- Zhang_full_moduleTraitCor_Pval_Vibrio_df %>% filter(pvalue <= 0.05)
Zhang_full_moduleTraitCor_Pval_Vibrio_df_sig #16, making it 0.051 adds one more module

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
Zhang_full_moduleTraitCor_Pval_LPS_df_sig_list <- Zhang_full_moduleTraitCor_Pval_LPS_df_sig$mod_names
Zhang_full_moduleTraitCor_Pval_LPS_df_sig_list_rm <- str_remove(Zhang_full_moduleTraitCor_Pval_LPS_df_sig_list, "ME")
Zhang_full_moduleTraitCor_Pval_Vibrio_df_sig_list <- Zhang_full_moduleTraitCor_Pval_Vibrio_df_sig$mod_names
Zhang_full_moduleTraitCor_Pval_Vibrio_df_sig_list_rm <- str_remove(Zhang_full_moduleTraitCor_Pval_Vibrio_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
names(Zhang_full_moduleTraitCor_Pval_Vibrio_df_sig_list_rm) <- c( "blue"      ,      "blue4"       ,    "honeydew1"   ,    "brown1"      ,    "yellow2"        ,"pink4"      ,     "salmon2"   ,      "greenyellow"  ,  
                                                             "thistle1"  ,      "darkgreen"   ,    "purple"      ,    "lightblue4"  ,    "darkolivegreen2","lightgreen" ,     "orangered" ,      "lightskyblue4"  )
matrix_common= Zhang_dds_rlog_matrix
moduleColors= Zhang_full_moduleColors
lookup =   C_gig_rtracklayer_apop_product_final

lookup_mod_apop_CG <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$transcript_id %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
Zhang_full_module_apop <- lapply(Zhang_full_moduleTraitCor_Pval_Vibrio_df_sig_list_rm,  lookup_mod_apop_CG)
Zhang_full_module_apop_df <- do.call(rbind,Zhang_full_module_apop)
Zhang_full_module_apop_df$mod_names <- gsub("\\..*","",row.names(Zhang_full_module_apop_df))
Zhang_full_module_apop_df$mod_names <- gsub("^","ME",Zhang_full_module_apop_df$mod_names)
# add module significance
Zhang_full_module_Vibrio_apop_df <- left_join(Zhang_full_module_apop_df,Zhang_full_moduleTraitCor_Pval_Vibrio_df)
Zhang_full_module_Vibrio_apop_df$exp <- "Zhang_Vibrio"

# specify names for list of lists
names(Zhang_full_moduleTraitCor_Pval_LPS_df_sig_list_rm) <- c("darkolivegreen4","chocolate4"     , "darkorange2"    , "honeydew"       , "pink"          , "sienna4"      ,   "darkseagreen2" ,  "magenta4"        ,
                                                         "grey60"         ,"tan3"           , "ivory"          , "black"          , "powderblue"    , "darkorange"   ,   "firebrick4"    ,  "lightpink2"     ,
                                                          "lightsteelblue", "midnightblue"  ,  "blue"          ,  "darkseagreen4" ,  "navajowhite"  ,  "magenta"     ,    "lightpink3"   ,   "mediumpurple2"  ,
                                                          "brown1"        , "yellow2"       ,  "pink4"         ,  "paleturquoise" ,  "turquoise"    ,  "thistle1"    ,    "cyan"         ,   "darkolivegreen2",
                                                          "green4"        , "palevioletred2",  "darkseagreen3" ,  "mediumpurple3" ,  "orangered3"   ,  "tan"         ,    "indianred3"   ,   "mistyrose"      ,
                                                          "orangered"     , "coral2"        ,  "plum3"         ,  "coral4"        ,  "lightskyblue4",  "blue3"       ,    "navajowhite2"   )
Zhang_LPS_full_module_apop <- lapply(Zhang_full_moduleTraitCor_Pval_LPS_df_sig_list_rm,  lookup_mod_apop_CG)
Zhang_LPS_full_module_apop_df <- do.call(rbind,Zhang_LPS_full_module_apop)
Zhang_LPS_full_module_apop_df$mod_names <- gsub("\\..*","",row.names(Zhang_LPS_full_module_apop_df))
Zhang_LPS_full_module_apop_df$mod_names <- gsub("^","ME",Zhang_LPS_full_module_apop_df$mod_names)
# add module significance
Zhang_LPS_full_module_apop_df <- left_join(Zhang_LPS_full_module_apop_df,Zhang_full_moduleTraitCor_Pval_LPS_df_sig)
Zhang_LPS_full_module_apop_df$exp <- "Zhang_LPS"

## Gene relationship to trait and important modules: Gene Significance and Module Membership
# do later

## Intramodular analysis: identifying genes with high GS and MM
# Using the GS and MM measures, we can identify genes that have a high significance for disease challenge
# as well as high module membership in interesting modules. As an example, we look at the brown module 
# that has the highest association with weight. We plot a scatterplot of Gene Significance vs. Module Membership in the brown module:

## IDENTIFY HUB GENES IN EACH SIG MODULE ##

#### Rubio ####
# Define numbers of genes and samples
Rubio_full_nGenes = ncol(Rubio_dds_rlog_matrix)
Rubio_full_nSamples = nrow(Rubio_dds_rlog_matrix)

# Recalculate MEs with color labels
Rubio_full_MEs0 = moduleEigengenes(Rubio_dds_rlog_matrix, Rubio_full_moduleColors)$eigengenes
Rubio_full_MEs = orderMEs(Rubio_full_MEs0)
Rubio_full_moduleTraitCor = cor(Rubio_full_MEs, Rubio_coldata_collapsed_binarize, use = "p");
Rubio_full_moduleTraitPvalue = corPvalueStudent(Rubio_full_moduleTraitCor, Rubio_full_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trai
sizeGrWindow(10,6)
# Will display correlations and their p-values
Rubio_full_textMatrix = paste(signif(Rubio_full_moduleTraitCor, 2), "\n(",
                         signif(Rubio_full_moduleTraitPvalue, 1), ")", sep = "");
dim(Rubio_full_textMatrix) = dim(Rubio_full_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = Rubio_full_moduleTraitCor,
               xLabels = names(Rubio_coldata_collapsed_binarize),
               yLabels = names(Rubio_full_MEs),
               ySymbols = names(Rubio_full_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = Rubio_full_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
Rubio_full_moduleTraitCor_df <- as.data.frame(Rubio_full_moduleTraitCor)
Rubio_full_moduleTraitCor_df$mod_names <- row.names(Rubio_full_moduleTraitCor_df)
Rubio_full_moduleTraitCor_NV_df <- Rubio_full_moduleTraitCor_df[,c("mod_names","Group.Non_virulent.vs.Control")]
Rubio_full_moduleTraitPvalue_df <- as.data.frame(Rubio_full_moduleTraitPvalue)
Rubio_full_moduleTraitPvalue_df$mod_names <- row.names(Rubio_full_moduleTraitPvalue_df )
Rubio_full_moduleTraitPvalue_NV_df <- Rubio_full_moduleTraitPvalue_df[,c("mod_names","Group.Non_virulent.vs.Control")]
colnames(Rubio_full_moduleTraitPvalue_NV_df)[2] <- "pvalue"
Rubio_full_moduleTraitCor_Pval_NV_df <- join(Rubio_full_moduleTraitCor_NV_df, Rubio_full_moduleTraitPvalue_NV_df, by = "mod_names")

# Significantly correlated modules
Rubio_full_moduleTraitCor_Pval_NV_df[order(Rubio_full_moduleTraitCor_Pval_NV_df$pvalue),]
class(Rubio_full_moduleTraitCor_Pval_NV_df$pvalue) # numeric
Rubio_full_moduleTraitCor_Pval_NV_df_sig <- Rubio_full_moduleTraitCor_Pval_NV_df %>% filter(pvalue <= 0.05)
Rubio_full_moduleTraitCor_Pval_NV_df_sig # 9

Rubio_full_moduleTraitCor_V_df <- Rubio_full_moduleTraitCor_df[,c("mod_names","Group.Virulent.vs.Control")]
Rubio_full_moduleTraitPvalue_V_df <- Rubio_full_moduleTraitPvalue_df[,c("mod_names","Group.Virulent.vs.Control")]
colnames(Rubio_full_moduleTraitPvalue_V_df)[2] <- "pvalue"
Rubio_full_moduleTraitCor_Pval_V_df <- join(Rubio_full_moduleTraitCor_V_df, Rubio_full_moduleTraitPvalue_V_df, by = "mod_names")

# Significantly correlated modules
Rubio_full_moduleTraitCor_Pval_V_df[order(Rubio_full_moduleTraitCor_Pval_V_df$pvalue),]
class(Rubio_full_moduleTraitCor_Pval_V_df$pvalue) # numeric
Rubio_full_moduleTraitCor_Pval_V_df_sig <- Rubio_full_moduleTraitCor_Pval_V_df %>% filter(pvalue <= 0.05)
Rubio_full_moduleTraitCor_Pval_V_df_sig # 11

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
Rubio_full_moduleTraitCor_Pval_NV_df_sig_list <- Rubio_full_moduleTraitCor_Pval_NV_df_sig$mod_names
Rubio_full_moduleTraitCor_Pval_V_df_sig_list <-  Rubio_full_moduleTraitCor_Pval_V_df_sig$mod_names
Rubio_full_moduleTraitCor_Pval_NV_df_sig_list_rm <- str_remove(Rubio_full_moduleTraitCor_Pval_NV_df_sig_list, "ME")
Rubio_full_moduleTraitCor_Pval_V_df_sig_list_rm <-  str_remove(Rubio_full_moduleTraitCor_Pval_V_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
names(Rubio_full_moduleTraitCor_Pval_NV_df_sig_list_rm) <- c("darkgreen", "grey60" ,   "magenta" ,  "black"  ,   "turquoise", "blue"   ,   "brown" ,    "yellow",    "red" )
names(Rubio_full_moduleTraitCor_Pval_V_df_sig_list_rm) <- c( "darkgreen",   "grey60"  ,    "magenta"  ,   "black"  ,     "turquoise" ,  "blue"    ,    "darkorange" , "greenyellow", "brown"   ,    "yellow" ,    
                                                             "red")

matrix_common= Rubio_dds_rlog_matrix
moduleColors= Rubio_full_moduleColors
lookup =  C_gig_rtracklayer_apop_product_final

lookup_mod_apop_CG <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$transcript_id %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
Rubio_NV_full_module_apop <- lapply(Rubio_full_moduleTraitCor_Pval_NV_df_sig_list_rm,  lookup_mod_apop_CG)
Rubio_NV_full_module_apop_df <- do.call(rbind,Rubio_NV_full_module_apop)
Rubio_NV_full_module_apop_df$mod_names <- gsub("\\..*","",row.names(Rubio_NV_full_module_apop_df))
Rubio_NV_full_module_apop_df$mod_names <- gsub("^","ME",Rubio_NV_full_module_apop_df$mod_names)
# add module significance
Rubio_NV_full_module_apop_df <- left_join(Rubio_NV_full_module_apop_df,Rubio_full_moduleTraitCor_Pval_NV_df_sig)
Rubio_NV_full_module_apop_df$exp <- "Rubio_NV"

Rubio_V_full_module_apop <- lapply(Rubio_full_moduleTraitCor_Pval_V_df_sig_list_rm,  lookup_mod_apop_CG)
Rubio_V_full_module_apop_df <- do.call(rbind,Rubio_V_full_module_apop)
Rubio_V_full_module_apop_df$mod_names <- gsub("\\..*","",row.names(Rubio_V_full_module_apop_df))
Rubio_V_full_module_apop_df$mod_names <- gsub("^","ME",Rubio_V_full_module_apop_df$mod_names)
# add module significance
Rubio_V_full_module_apop_df <- left_join(Rubio_V_full_module_apop_df,Rubio_full_moduleTraitCor_Pval_V_df_sig)
Rubio_V_full_module_apop_df$exp <- "Rubio_V"

## Gene relationship to trait and important modules: Gene Significance and Module Membership
# do later

## Intramodular analysis: identifying genes with high GS and MM
# Using the GS and MM measures, we can identify genes that have a high significance for disease challenge
# as well as high module membership in interesting modules. As an example, we look at the brown module 
# that has the highest association with weight. We plot a scatterplot of Gene Significance vs. Module Membership in the brown module:

## IDENTIFY HUB GENES IN EACH SIG MODULE ##

# look at this later

#### deLorg Res ####
# Define numbers of genes and samples
deLorg_Res_full_nGenes = ncol(deLorgeril_Resistant_dds_vst_matrix)
deLorg_Res_full_nSamples = nrow(deLorgeril_Resistant_dds_vst_matrix)

# Recalculate MEs with color labels
deLorg_Res_full_MEs0 = moduleEigengenes(deLorgeril_Resistant_dds_vst_matrix, deLorg_Res_full_moduleColors)$eigengenes
deLorg_Res_full_MEs = orderMEs(deLorg_Res_full_MEs0)
deLorg_Res_full_moduleTraitCor = cor(deLorg_Res_full_MEs, deLorgeril_Resistant_coldata_collapsed_binarize, use = "p");
deLorg_Res_full_moduleTraitPvalue = corPvalueStudent(deLorg_Res_full_moduleTraitCor, deLorg_Res_full_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trai
sizeGrWindow(10,6)
# Will display correlations and their p-values
deLorg_Res_full_textMatrix = paste(signif(deLorg_Res_full_moduleTraitCor, 2), "\n(",
                              signif(deLorg_Res_full_moduleTraitPvalue, 1), ")", sep = "");
dim(deLorg_Res_full_textMatrix) = dim(deLorg_Res_full_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = deLorg_Res_full_moduleTraitCor,
               xLabels = names(deLorgeril_Resistant_coldata_collapsed_binarize),
               yLabels = names(deLorg_Res_full_MEs),
               ySymbols = names(deLorg_Res_full_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = deLorg_Res_full_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
deLorg_Res_full_moduleTraitCor_df <- as.data.frame(deLorg_Res_full_moduleTraitCor)
deLorg_Res_full_moduleTraitCor_df$mod_names <- row.names(deLorg_Res_full_moduleTraitCor_df)
deLorg_Res_full_moduleTraitCor_df <- deLorg_Res_full_moduleTraitCor_df[,c("mod_names","Group.AF21_Resistant.vs.AF21_Resistant_control")]
deLorg_Res_full_moduleTraitPvalue_df <- as.data.frame(deLorg_Res_full_moduleTraitPvalue)
deLorg_Res_full_moduleTraitPvalue_df$mod_names <- row.names(deLorg_Res_full_moduleTraitPvalue_df )
deLorg_Res_full_moduleTraitPvalue_df <- deLorg_Res_full_moduleTraitPvalue_df[,c("mod_names","Group.AF21_Resistant.vs.AF21_Resistant_control")]
colnames(deLorg_Res_full_moduleTraitPvalue_df)[2] <- "pvalue"
deLorg_Res_full_moduleTraitCor_Pval_df <- join(deLorg_Res_full_moduleTraitCor_df, deLorg_Res_full_moduleTraitPvalue_df, by = "mod_names")

# Significantly correlated modules
deLorg_Res_full_moduleTraitCor_Pval_df[order(deLorg_Res_full_moduleTraitCor_Pval_df$pvalue),]
class(deLorg_Res_full_moduleTraitCor_Pval_df$pvalue) # numeric
deLorg_Res_full_moduleTraitCor_Pval_df_sig <- deLorg_Res_full_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
deLorg_Res_full_moduleTraitCor_Pval_df_sig #8
deLorg_Res_full_moduleTraitCor_Pval_df_sig <- deLorg_Res_full_moduleTraitCor_Pval_df_sig %>% filter(mod_names != "MEgrey")

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
deLorg_Res_full_moduleTraitCor_Pval_df_sig_list <- deLorg_Res_full_moduleTraitCor_Pval_df_sig$mod_names
deLorg_Res_full_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(deLorg_Res_full_moduleTraitCor_Pval_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
names(deLorg_Res_full_moduleTraitCor_Pval_df_sig_list_rm) <- c("cyan" ,      "green",      "yellow" ,    "black" ,     "darkorange", "pink"   ,    "turquoise"   )

matrix_common= deLorgeril_Resistant_dds_vst_matrix
moduleColors= deLorg_Res_full_moduleColors
lookup =   C_gig_rtracklayer_apop_product_final

lookup_mod_apop_CG <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$transcript_id %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
deLorg_Res_full_module_apop <- lapply(deLorg_Res_full_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop_CG)
deLorg_Res_full_module_apop_df <- do.call(rbind,deLorg_Res_full_module_apop)
deLorg_Res_full_module_apop_df$mod_names <- gsub("\\..*","",row.names(deLorg_Res_full_module_apop_df))
deLorg_Res_full_module_apop_df$mod_names <- gsub("^","ME",deLorg_Res_full_module_apop_df$mod_names)
# add module significance
deLorg_Res_full_module_apop_df <- left_join(deLorg_Res_full_module_apop_df,deLorg_Res_full_moduleTraitCor_Pval_df_sig)
deLorg_Res_full_module_apop_df$exp <- "deLorg_Res"

#### deLorg Sus ####
# Define numbers of genes and samples
deLorg_Sus_full_nGenes = ncol(deLorgeril_Susceptible_dds_vst_matrix)
deLorg_Sus_full_nSamples = nrow(deLorgeril_Susceptible_dds_vst_matrix)

# Recalculate MEs with color labels
deLorg_Sus_full_MEs0 = moduleEigengenes(deLorgeril_Susceptible_dds_vst_matrix, deLorg_Sus_full_moduleColors)$eigengenes
deLorg_Sus_full_MEs = orderMEs(deLorg_Sus_full_MEs0)
deLorg_Sus_full_moduleTraitCor = cor(deLorg_Sus_full_MEs, deLorgeril_Susceptible_coldata_collapsed_binarize , use = "p");
deLorg_Sus_full_moduleTraitPvalue = corPvalueStudent(deLorg_Sus_full_moduleTraitCor, deLorg_Sus_full_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trai
sizeGrWindow(10,6)
# Will display correlations and their p-values
deLorg_Sus_full_textMatrix = paste(signif(deLorg_Sus_full_moduleTraitCor, 2), "\n(",
                              signif(deLorg_Sus_full_moduleTraitPvalue, 1), ")", sep = "");
dim(deLorg_Sus_full_textMatrix) = dim(deLorg_Sus_full_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = deLorg_Sus_full_moduleTraitCor,
               xLabels = names(deLorgeril_Susceptible_coldata_collapsed_binarize ),
               yLabels = names(deLorg_Sus_full_MEs),
               ySymbols = names(deLorg_Sus_full_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = deLorg_Sus_full_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
deLorg_Sus_full_moduleTraitCor_df <- as.data.frame(deLorg_Sus_full_moduleTraitCor)
deLorg_Sus_full_moduleTraitCor_df$mod_names <- row.names(deLorg_Sus_full_moduleTraitCor_df)
deLorg_Sus_full_moduleTraitCor_df <- deLorg_Sus_full_moduleTraitCor_df[,c("mod_names","Group.AF11_Susceptible.vs.AF11_Susceptible_control")]
deLorg_Sus_full_moduleTraitPvalue_df <- as.data.frame(deLorg_Sus_full_moduleTraitPvalue)
deLorg_Sus_full_moduleTraitPvalue_df$mod_names <- row.names(deLorg_Sus_full_moduleTraitPvalue_df )
deLorg_Sus_full_moduleTraitPvalue_df <- deLorg_Sus_full_moduleTraitPvalue_df[,c("mod_names","Group.AF11_Susceptible.vs.AF11_Susceptible_control")]
colnames(deLorg_Sus_full_moduleTraitPvalue_df)[2] <- "pvalue"
deLorg_Sus_full_moduleTraitCor_Pval_df <- join(deLorg_Sus_full_moduleTraitCor_df, deLorg_Sus_full_moduleTraitPvalue_df, by = "mod_names")

# Significantly correlated modules
deLorg_Sus_full_moduleTraitCor_Pval_df[order(deLorg_Sus_full_moduleTraitCor_Pval_df$pvalue),]
class(deLorg_Sus_full_moduleTraitCor_Pval_df$pvalue) # numeric
deLorg_Sus_full_moduleTraitCor_Pval_df_sig <- deLorg_Sus_full_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
deLorg_Sus_full_moduleTraitCor_Pval_df_sig # 6

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
deLorg_Sus_full_moduleTraitCor_Pval_df_sig_list <- deLorg_Sus_full_moduleTraitCor_Pval_df_sig$mod_names
deLorg_Sus_full_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(deLorg_Sus_full_moduleTraitCor_Pval_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
names(deLorg_Sus_full_moduleTraitCor_Pval_df_sig_list_rm) <- c( "turquoise",    "purple" ,      "magenta" ,     "midnightblue", "white"  ,      "blue"  )

matrix_common= deLorgeril_Susceptible_dds_vst_matrix
moduleColors= deLorg_Sus_full_moduleColors
lookup =   C_gig_rtracklayer_apop_product_final

lookup_mod_apop_CG <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$transcript_id %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
deLorg_Sus_full_module_apop <- lapply(deLorg_Sus_full_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop_CG)
deLorg_Sus_full_module_apop_df <- do.call(rbind,deLorg_Sus_full_module_apop)
deLorg_Sus_full_module_apop_df$mod_names <- gsub("\\..*","",row.names(deLorg_Sus_full_module_apop_df))
deLorg_Sus_full_module_apop_df$mod_names <- gsub("^","ME",deLorg_Sus_full_module_apop_df$mod_names)
# add module significance
deLorg_Sus_full_module_apop_df <- left_join(deLorg_Sus_full_module_apop_df,deLorg_Sus_full_moduleTraitCor_Pval_df_sig)
deLorg_Sus_full_module_apop_df$exp <- "deLorg_Sus"

#### He ####
# Define numbers of genes and samples
He_full_nGenes = ncol(He_dds_vst_matrix)
He_full_nSamples = nrow(He_dds_vst_matrix)

# Recalculate MEs with color labels
He_full_MEs0 = moduleEigengenes(He_dds_vst_matrix, He_full_moduleColors)$eigengenes
He_full_MEs = orderMEs(He_full_MEs0)
He_full_moduleTraitCor = cor(He_full_MEs, He_coldata_collapsed_binarize, use = "p");
He_full_moduleTraitPvalue = corPvalueStudent(He_full_moduleTraitCor, He_full_nSamples)

# Graph and color code each the strength of association (correlation) of module eigengenes and trai
sizeGrWindow(10,6)
# Will display correlations and their p-values
He_full_textMatrix = paste(signif(He_full_moduleTraitCor, 2), "\n(",
                      signif(He_full_moduleTraitPvalue, 1), ")", sep = "");
dim(He_full_textMatrix) = dim(He_full_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot, color coded by correlation value (red means more highly positively correlated,
# green is more negatively correlated)
labeledHeatmap(Matrix = He_full_moduleTraitCor,
               xLabels = names(He_coldata_collapsed_binarize),
               yLabels = names(He_full_MEs),
               ySymbols = names(He_full_MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = He_full_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

# Which modules have the highest associations with disease (high correlation and low P value)?
He_full_moduleTraitCor_df <- as.data.frame(He_full_moduleTraitCor)
He_full_moduleTraitCor_df$mod_names <- row.names(He_full_moduleTraitCor_df)
He_full_moduleTraitCor_df <- He_full_moduleTraitCor_df[,c("mod_names","Group.OsHV1.vs.control")]
He_full_moduleTraitPvalue_df <- as.data.frame(He_full_moduleTraitPvalue)
He_full_moduleTraitPvalue_df$mod_names <- row.names(He_full_moduleTraitPvalue_df )
He_full_moduleTraitPvalue_df <- He_full_moduleTraitPvalue_df[,c("mod_names","Group.OsHV1.vs.control")]
colnames(He_full_moduleTraitPvalue_df)[2] <- "pvalue"
He_full_moduleTraitCor_Pval_df <- join(He_full_moduleTraitCor_df, He_full_moduleTraitPvalue_df, by = "mod_names")

# Significantly correlated modules
He_full_moduleTraitCor_Pval_df[order(He_full_moduleTraitCor_Pval_df$pvalue),]
class(He_full_moduleTraitCor_Pval_df$pvalue) # numeric
He_full_moduleTraitCor_Pval_df_sig <- He_full_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
He_full_moduleTraitCor_Pval_df_sig # 9

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
He_full_moduleTraitCor_Pval_df_sig_list <- He_full_moduleTraitCor_Pval_df_sig$mod_names
He_full_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(He_full_moduleTraitCor_Pval_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
names(He_full_moduleTraitCor_Pval_df_sig_list_rm) <- c("greenyellow", "skyblue"  ,   "blue" ,       "turquoise",   "darkgrey"  ,  "purple" ,     "brown",       "red"   ,      "yellow" )

matrix_common= He_dds_vst_matrix
moduleColors= He_full_moduleColors
lookup =   C_gig_rtracklayer_apop_product_final

lookup_mod_apop_CG <- function(list) {
  list_vec <- colnames(matrix_common)[moduleColors == list]
  list_apop <- lookup[lookup$transcript_id %in% list_vec,]
  list_apop_short <- list_apop[,c("product","transcript_id","gene")]
}
# specify names for list of lists
He_full_module_apop <- lapply(He_full_moduleTraitCor_Pval_df_sig_list_rm,  lookup_mod_apop_CG)
He_full_module_apop_df <- do.call(rbind,He_full_module_apop)
He_full_module_apop_df$mod_names <- gsub("\\..*","",row.names(He_full_module_apop_df))
He_full_module_apop_df$mod_names <- gsub("^","ME",He_full_module_apop_df$mod_names)
# add module significance
He_full_module_apop_df <- left_join(He_full_module_apop_df,He_full_moduleTraitCor_Pval_df_sig)
He_full_module_apop_df$exp <- "He"

#### MASTER C_GIG FULL MODULE-TRAIT-SIG APOP GENES ####
colnames(Zhang_LPS_full_module_apop_df)[5] <- "mod_signif"
colnames(Zhang_full_module_Vibrio_apop_df)[5] <- "mod_signif"
colnames(Rubio_NV_full_module_apop_df   )[5] <- "mod_signif"
colnames(Rubio_V_full_module_apop_df    )[5] <- "mod_signif"
colnames(deLorg_Res_full_module_apop_df )[5] <- "mod_signif"
colnames(deLorg_Sus_full_module_apop_df )[5] <- "mod_signif"   
colnames(He_full_module_apop_df)[5] <- "mod_signif" 

C_gig_full_all_exp_mod_sig_apop <- rbind(Zhang_LPS_full_module_apop_df,
                                    Zhang_full_module_Vibrio_apop_df,
                                    Rubio_NV_full_module_apop_df  ,
                                    Rubio_V_full_module_apop_df   ,
                                    deLorg_Res_full_module_apop_df,
                                    deLorg_Sus_full_module_apop_df,
                                    He_full_module_apop_df)
nrow(C_gig_full_all_exp_mod_sig_apop) # 1235
C_gig_full_all_exp_mod_sig_apop_positive <- C_gig_full_all_exp_mod_sig_apop %>% filter(mod_signif >0)

#### ANALYSIS OF ALL IAP INTERACTIONS WITH FULL WGCNA DATA RUN FOR EACH EXPERIMENT SEPARATELY #### 
# create combined sig apop modules list 
C_vir_full_all_exp_mod_sig_apop$Species <- "Crassostrea_virginica"
C_gig_full_all_exp_mod_sig_apop$Species <- "Crassostrea_gigas"
C_vir_C_gig_full_all_exp_mod_sig_apop <- rbind(C_vir_full_all_exp_mod_sig_apop, C_gig_full_all_exp_mod_sig_apop)
# are all experiments found?
levels(factor(C_vir_C_gig_full_all_exp_mod_sig_apop$exp)) # all except ROD_Sus
#[1] "deLorg_Res"         "deLorg_Sus"         "Dermo_Sus"          "Dermo_Tol"          "He"                 "Pro_RE22_Pro_RI"    "Pro_RE22_Pro_S4"   
#[8] "Pro_RE22_RE22_full" "Probiotic"          "ROD_Res"                 "Rubio_NV"           "Rubio_V"            "Zhang_LPS"         
#[15] "Zhang_Vibrio"  

### Investigating apoptosis interaction partners for each domain structure type
# Use dataframes for domain structure loaded at top of whole script 
IAP_domain_structure_XM_CG
IAP_domain_structure_XM_CV

# join XM onto IAP_domain_structure_XM_CV using C_vir_rtracklayer_apop_product_final
IAP_domain_structure_XM_CV_XM <- left_join(IAP_domain_structure_XM_CV[,-26], C_vir_rtracklayer_apop_product_final[,c("ID","transcript_id")], by = "ID") 
      #get rid of IAP_domain_structure_XM_CV empty transcript_id column before joining

# combine data IAP_domain_name data frames, put columns in correct order and remove NA domain names
IAP_domain_structure_XM_filter <- rbind(IAP_domain_structure_XM_CG[,c("transcript_id","Domain_Name")], IAP_domain_structure_XM_CV_XM[,c("transcript_id","Domain_Name")]) %>%
  # change NA to be "not classified"
  mutate(Domain_Name = case_when(is.na(Domain_Name) ~ "not_classified",
                                       TRUE ~ Domain_Name))
  
# create search lists for each domain structure
IAP_domain_structure_df <- IAP_domain_structure_XM_filter %>%
  group_by(Domain_Name) %>%
  # put all proteins into a string
  dplyr::summarise(transcript_id = paste0(transcript_id, collapse = ",")) %>%
  filter(!is.na(Domain_Name)) %>% 
  # set column names
  column_to_rownames(., var= "Domain_Name")

# split strings into vector for grepping
IAP_domain_structure_df$transcript_id <- strsplit(IAP_domain_structure_df$transcript_id, ",")

# create list of lists
IAP_domain_structure_list <- as.list(as.data.frame(t(IAP_domain_structure_df)))

# loop through to get list of module hits for each group of transcripts
IAP_domain_structure_WGCNA_hits <- vector('list', length(IAP_domain_structure_list))
names(IAP_domain_structure_WGCNA_hits) <- names(IAP_domain_structure_list) # set names
for(i in seq_along(IAP_domain_structure_list)){
  for(j in seq_along(IAP_domain_structure_list[[i]])){
    IAP_domain_structure_WGCNA_hits[[i]][[j]]<- C_vir_C_gig_full_all_exp_mod_sig_apop %>% 
      group_by(mod_names, exp) %>%
      filter(any(grepl(paste(IAP_domain_structure_list[[i]][[j]], collapse="|"), transcript_id)))
  }
}
# put in dataframe (using purrr yay!)
IAP_domain_structure_WGCNA_hits_df <- map_df(IAP_domain_structure_WGCNA_hits, ~bind_rows(., .id="Domain_Name"), .id="Domain_Name") 
distinct(IAP_domain_structure_WGCNA_hits_df, Domain_Name, mod_names, exp) %>% View() # 140 module hits (some duplicated because include multiple domain types)
View(distinct(IAP_domain_structure_WGCNA_hits_df[c("mod_names", "exp")]))# 73 distinct module names
# are all experiments found- all except ROD_Sus
levels(factor(IAP_domain_structure_WGCNA_hits_df$exp))
# [1] "deLorg_Res"         "deLorg_Sus"         "Dermo_Sus"          "Dermo_Tol"          "He"                 "Pro_RE22_Pro_RI"    "Pro_RE22_Pro_S4"   
# [8] "Pro_RE22_RE22_full" "Probiotic"          "ROD_Res"                      "Rubio_NV"           "Rubio_V"            "Zhang_LPS"         
# [15] "Zhang_Vibrio"                           

# Make table with just the IAP hits for each one 
C_vir_rtracklayer_transcripts <- C_vir_rtracklayer %>% filter(type == "mRNA")
IAP_domain_structure_WGCNA_hits_df_IAP <- IAP_domain_structure_WGCNA_hits_df[IAP_domain_structure_WGCNA_hits_df$transcript_id %in% IAP_domain_structure_XM_filter$transcript_id,]

# remove those that only contain one transcript 
IAP_domain_structure_WGCNA_hits_df_modsize <- IAP_domain_structure_WGCNA_hits_df %>% dplyr::count(mod_names, exp) %>% arrange(desc(n)) %>% filter(n>1)
nrow(IAP_domain_structure_WGCNA_hits_df_modsize) # 69 modules with more than 1 transcript

## Remove modules with only 1 transcript (an IAP) using semi_join! # 
IAP_domain_structure_WGCNA_hits_df <- dplyr::semi_join(IAP_domain_structure_WGCNA_hits_df,  IAP_domain_structure_WGCNA_hits_df_modsize[,c("mod_names","exp")])
# check
distinct(IAP_domain_structure_WGCNA_hits_df, mod_names, exp) # 69 total without any module repeats
View(distinct(IAP_domain_structure_WGCNA_hits_df, Domain_Name, mod_names, exp))

## Count IAPs of each type in each experiment and make table
IAP_domain_structure_WGCNA_hits_df_IAP_count_exp_table <- IAP_domain_structure_WGCNA_hits_df_IAP %>%
  ungroup() %>% group_by(exp, Species, Domain_Name) %>% dplyr::summarise(IAP_count = n())

# Make table for paper (Table 4)
# Generate WGCNA IAP table across modules similar to the LFC table
IAP_domain_structure_WGCNA_hits_df_IAP_count_exp_TABLE <- IAP_domain_structure_WGCNA_hits_df_IAP_count_exp_table %>%
  # mutate experiment level names so I can change and add markdown formatting below
  # spread table into wide format
  ungroup() %>% dplyr::select(exp,IAP_count,Domain_Name) %>%
  spread(exp, IAP_count, fill = 0) %>%
  gt::gt(rowname_col = "Domain_Name") %>%
  tab_header(title = gt::md("**Domain Structure of IAPs in all Significant WGCNA Modules Per Experiment**")) %>%
  cols_label( Domain_Name = md("**Domain Structure Type**"),
    "deLorg_Res" = "de Lorgeril OsHv-1 Resistant",
    "deLorg_Sus" = "de Lorgeril OsHv-1 Susceptible",
    "Dermo_Sus"  = "Dermo Susceptible",
    "Dermo_Tol" = "Dermo Tolerant",
    "He"  = "He OsHV-1",
    "Pro_RE22_Pro_RI" = "Lab Pro. RI",
    "Pro_RE22_Pro_S4" = "Lab Pro. S4",
    "Pro_RE22_RE22_full" = "Lab RE22",
    "Probiotic"          = "Hatchery Pro. RI",
    "ROD_Res" = "ROD Resistant",
    "Rubio_NV" = md("Rubio. Non-virulent\n*Vibrio* spp."),
    "Rubio_V" =  md("Rubio. Virulent\n*Vibrio* spp."),
    "Zhang_LPS" = md("Zhang LPS, *M. Lut*"),
    "Zhang_Vibrio"   = md("Zhang *Vibrio* spp.")) %>%
  #add total column
  grand_summary_rows(fns = list(Total="sum"), formatter = fmt_number, decimals = 0) %>%
  # add spanner over experiments for each species
  tab_spanner(
    label = md("*Crassostrea gigas*"),
    columns = vars(deLorg_Res,deLorg_Sus,He, Rubio_NV, Rubio_V,Zhang_LPS, Zhang_Vibrio)) %>%
  tab_spanner(
    label = md("*Crassostrea virginica*"),
    columns = vars(Dermo_Sus, Dermo_Tol,Pro_RE22_Pro_RI,Pro_RE22_Pro_S4,Pro_RE22_RE22_full,Probiotic,ROD_Res)) %>%
  tab_source_note(source_note = md("\\* = *IAP Domain identified by Interproscan and not CDD search*")) %>%
  tab_options(table.font.color = "black")
# save as png
gtsave(IAP_domain_structure_WGCNA_hits_df_IAP_count_exp_TABLE, "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/TABLES/IAP_domain_structure_WGCNA_hits_df_IAP_count_exp_TABLE.png")

### FIND DOMAIN STRUCTURE SHARED ACROSS CHALLENGE TYPES 
# Goal is see if there is a response that might be specific to the IAP domain structure rather than a challenge specific response
levels(factor(IAP_domain_structure_WGCNA_hits_df$Domain_Name))
#[1] "BIR*"                  "BIR*-DD-RING"          "not_classified"        "NZBIR-TII-UBA-DD-RING" "TI-TII-DD-RING"        "TI-TII-RING"           "TI-TII-TII-UBA-RING"   "TII"                   "TII-BIR6-E2"          
#[10] "TII-DD"                "TII-DD-RING"           "TII-RING"              "TII-TII"               "TII-TII-RING"          "TX-TII"     
# ALL domain structure types were found

### Create Comb_domains for those modules that hit to multiple domains ###
## Condense data frame so that module members in modules that hit to multiple domain structures don't get counted twice
# Which domain structure types found in the same modules? Creating combined domain 
IAP_domain_structure_WGCNA_hits_df_modules_domain_hits <- IAP_domain_structure_WGCNA_hits_df %>% 
  distinct(Domain_Name, mod_names, exp) %>% 
  # count the total number of modules for each domain type and create comb domain if multiple for an experiment
  group_by(mod_names, exp) %>% dplyr::mutate(count = n(), comb_domain = paste(Domain_Name, collapse = ",")) %>% 
  distinct(mod_names, exp, count, comb_domain) %>% arrange(exp)
# join back with species for grouping in the code below
IAP_domain_structure_WGCNA_hits_df_modules_domain_hits <- left_join(IAP_domain_structure_WGCNA_hits_df_modules_domain_hits,   unique(IAP_domain_structure_WGCNA_hits_df[,c("exp","Species")]))

# how many domains are there with the comb_domains
comb_domain <-  as.data.frame(levels(factor(IAP_domain_structure_WGCNA_hits_df_modules_domain_hits$comb_domain)))
colnames(comb_domain )[1] <- "comb_domain"
nrow(comb_domain) # 36 total
# any modules with only 1 domain type?
comb_domain_unique <- comb_domain %>% filter(!grepl(",",comb_domain)) # 12 of the modules were found uniquely in a module
    #              comb_domain
    #1           BIR*-DD-RING
    #2         not_classified
    #3  NZBIR-TII-UBA-DD-RING
    #4         TI-TII-DD-RING
    #5            TI-TII-RING
    #6    TI-TII-TII-UBA-RING
    #7            TII-BIR6-E2
    #8                 TII-DD
    #9            TII-DD-RING
    #10               TII-TII
    #11          TII-TII-RING
    #12                TX-TII
# meaning 25 different domain combinations
# 25/36 = 69.4% of domains are combo domains 
comb_domain_unique$comb_domain_type <- "unique"

# are the domain combos common across experiments?
IAP_domain_structure_WGCNA_hits_df_modules_domain_hit_common_exp_species <-IAP_domain_structure_WGCNA_hits_df_modules_domain_hits %>% group_by(comb_domain, Species, exp) %>% 
  dplyr::count() %>% arrange(desc(n))
IAP_domain_structure_WGCNA_hits_df_modules_domain_hit_common_exp <-IAP_domain_structure_WGCNA_hits_df_modules_domain_hits %>% group_by(comb_domain) %>% 
  dplyr::count() %>% arrange(desc(n))
# which are only found once: 
IAP_domain_structure_WGCNA_hits_df_modules_domain_hit_common_exp %>% filter(n ==1)
# A tibble: 17 x 2
# Groups:   comb_domain [17]
#comb_domain                                                                    n
#<chr>                                                                      <int>
#  1 BIR*,BIR*-DD-RING                                                              1
#2 BIR*,BIR*-DD-RING,not_classified,TI-TII-DD-RING,TII-DD-RING,TII-TII,TX-TII     1
#3 BIR*,BIR*-DD-RING,not_classified,TI-TII-DD-RING,TII-TII,TX-TII                 1
#4 BIR*,TII-BIR6-E2                                                               1
#5 BIR*,TII-DD-RING                                                               1
#6 not_classified,TI-TII-DD-RING,TI-TII-TII-UBA-RING,TII-RING                     1
#7 not_classified,TI-TII-DD-RING,TII-DD,TII-DD-RING                               1
#8 not_classified,TI-TII-DD-RING,TII-TII-RING                                     1
#9 not_classified,TI-TII-DD-RING,TX-TII                                           1
#10 not_classified,TI-TII-TII-UBA-RING,TII-BIR6-E2                                 1
#11 not_classified,TII-DD-RING                                                     1
#12 not_classified,TII-RING,TX-TII                                                 1
#13 not_classified,TII-TII,TII-TII-RING                                            1
#14 TI-TII-DD-RING,TII-DD-RING                                                     1
#15 TII-BIR6-E2,TII-TII-RING                                                       1
#16 TII-TII-RING,TX-TII                                                            1
#17 TII,TII-DD-RING                                                                1

## How common are the individual domain name types across experiments if you split up the combo domains within the actual data 
IAP_domain_structure_WGCNA_hits_df_modules_domain_hit_common_exp_separate <- IAP_domain_structure_WGCNA_hits_df_modules_domain_hit_common_exp_species %>%
  # separate into rows
  separate_rows(comb_domain, sep = ",") %>% 
  # add together the numbers for each type 
  group_by(comb_domain, Species) %>% dplyr::summarize(total_times_across_exp_modules = sum(n))

# Export as formatted table 
IAP_domain_structure_WGCNA_hits_df_modules_domain_hit_common_exp_separate_TABLE <- IAP_domain_structure_WGCNA_hits_df_modules_domain_hit_common_exp_separate %>%
  # mutate experiment level names so I can change and add markdown formatting below
  # spread table into wide format
  ungroup() %>% 
  spread(Species, total_times_across_exp_modules, fill = 0) %>%
  gt::gt(rowname_col = "comb_domain") %>%
  tab_header(title = gt::md("**Number of Modules Including Domain Structure Per Experiment**")) %>%
  cols_label(Crassostrea_gigas = md("*Crassostrea gigas*"),
        Crassostrea_virginica = md("*Crassostrea virginica*")) %>%
    tab_source_note(source_note = md("\\* = *IAP Domain identified by Interproscan and not CDD search*")) %>%
  tab_options(table.font.color = "black")
# save as png
gtsave(IAP_domain_structure_WGCNA_hits_df_modules_domain_hit_common_exp_separate_TABLE, "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/ANNOTATION_DATA_FIGURES/TABLES/IAP_domain_structure_WGCNA_hits_df_modules_domain_hit_common_exp_separate_TABLE.png")

# Which experiments have the most different domain types overall 
IAP_domain_structure_WGCNA_hits_df_modules_domain_hit_common_separate_EXP <- IAP_domain_structure_WGCNA_hits_df_modules_domain_hit_common_exp_species %>%
  # separate into rows
  separate_rows(comb_domain, sep = ",") %>% 
  # add find all the different comb_domains associated with each experiment and add them up
   dplyr::distinct(Species,exp, comb_domain) %>% dplyr::count(exp) %>%
  arrange(desc(n))

# how common are the individual domain types within just the names (including both combo and unique)? Are some more specific than others?
comb_domain_freq <- as.data.frame(str_split(comb_domain$comb_domain, ",") %>% unlist(.))
colnames(comb_domain_freq)[1] <- "Domain_Name"
comb_domain_freq %>% dplyr::mutate(total = nrow(.)) %>% group_by(Domain_Name) %>% dplyr::mutate(count = n(), percent = count/total*100) %>% arrange(desc(count)) %>% distinct(Domain_Name, total, count, percent)

#Domain_Name           total count percent
#<chr>                 <int> <int>   <dbl>
#  1 not_classified           86    17   19.8 
#2 TI-TII-DD-RING           86    11   12.8 
#3 TII-DD-RING              86     9   10.5 
#4 TX-TII                   86     8    9.30
#5 BIR*                     86     6    6.98
#6 TII-DD                   86     6    6.98
#7 TII-TII-RING             86     6    6.98
#8 BIR*-DD-RING             86     5    5.81
#9 TII-BIR6-E2              86     5    5.81
#10 TII-TII                  86     4    4.65
#11 TI-TII-TII-UBA-RING      86     3    3.49
#12 TII-RING                 86     3    3.49
#13 NZBIR-TII-UBA-DD-RING    86     1    1.16
#14 TI-TII-RING              86     1    1.16
#15 TII                      86     1    1.16


## Create condensed data frame using the new comb_domain types
IAP_domain_structure_WGCNA_hits_df_condensed <- left_join(IAP_domain_structure_WGCNA_hits_df_modules_domain_hits, unique(IAP_domain_structure_WGCNA_hits_df[,c("mod_names","exp","product","transcript_id")]))

### Compare Usage of domains between experiments ###
# Calculate Count products per module in each experiment
IAP_domain_structure_WGCNA_hits_exp <- IAP_domain_structure_WGCNA_hits_df_condensed  %>%
  group_by(mod_names,exp, comb_domain) %>%
  dplyr::summarise(products_per_mod_exp = n()) %>% 
  distinct(mod_names, exp, comb_domain, products_per_mod_exp) %>%
  arrange(desc(products_per_mod_exp))

### Which experiments have the most module types involved? (to help assess whether a more broad or more specific response is elicited)
# Calculate number of differnt modules types for each experimet 
IAP_domain_structure_WGCNA_hits_df_condensed  %>% ungroup() %>%
  dplyr::distinct(exp, comb_domain) %>% # only want to count 1 per experiment
dplyr::count(exp) %>% arrange(desc(n))
# A tibble: 14 x 2
#exp                    n
#<chr>              <int>
#  1 Pro_RE22_Pro_RI       10
#2 Pro_RE22_Pro_S4        9
#3 Zhang_LPS              9
#4 Rubio_V                7
#5 Rubio_NV               6
#6 He                     5
#7 Dermo_Tol              4
#8 Zhang_Vibrio           3
#9 deLorg_Res             2
#10 deLorg_Sus             2
#11 Pro_RE22_RE22_full     2
#12 Probiotic              2
#13 Dermo_Sus              1
#14 ROD_Res                1

## Join experiments with challenge type: viral, bacterial, parasitic,
challenge_type <- data.frame(exp =c(
  "deLorg_Res",     
  "deLorg_Sus",     
  "Dermo_Tol", 
  "Dermo_Sus",
  "He",             
  "Pro_RE22_Pro_RI",
  "Pro_RE22_Pro_S4",
  "Pro_RE22_RE22_full"  ,
  "Probiotic",      
  "ROD_Res",        
  "ROD_Sus",        
  "Rubio_NV",       
  "Rubio_V",        
  "Zhang_LPS",      
  "Zhang_Vibrio"),  
  challenge_type = c(
    "viral" ,
    "viral",
    "parasite",
    "parasite",
    "viral",
    "bacterial",
    "bacterial",
    "bacterial",
    "bacterial",
    "bacterial",
    "bacterial",
    "bacterial",
    "bacterial", 
    "bacterial",
    "bacterial"))

IAP_domain_structure_WGCNA_hits_df_condensed_type <- left_join(IAP_domain_structure_WGCNA_hits_df_condensed, challenge_type)
# check
levels(factor(IAP_domain_structure_WGCNA_hits_df_condensed_type$exp)) # no NAs

# classify comb_domain as unique or combo
IAP_domain_structure_WGCNA_hits_df_condensed_type <- left_join(IAP_domain_structure_WGCNA_hits_df_condensed_type ,comb_domain_unique) %>% mutate(comb_domain_type = case_when(
  is.na(.$comb_domain_type) ~ "combo",
  TRUE ~ .$comb_domain_type))
    
### Upset (kinda) plot of domain types used across experiments
# see https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html#make-the-combination-matrix
IAP_domain_structure_WGCNA_hits_exp_upset <- IAP_domain_structure_WGCNA_hits_df_condensed_type %>% dplyr::distinct(mod_names, exp, comb_domain, challenge_type) %>%
  group_by(comb_domain,exp) %>% dplyr::mutate(count=n())

IAP_domain_structure_WGCNA_hits_exp_upset_plot <- ggplot(IAP_domain_structure_WGCNA_hits_exp_upset, aes(y=comb_domain, x=exp, fill= count)) + geom_tile() +
  facet_grid(.~challenge_type, scales = "free", space = "free") + 
  theme(axis.text.x = element_text(angle = 90, hjust =1, size = 14),
        axis.text.y = element_text(size=14),
        plot.title = element_text(size = 16))+
  labs(x = "Experiment", y = "Domain Name or Combination", title = "Occurence of Domain Structure Combinations Across Experiments", fill = "Modules\n Per Exp.") 

ggsave(plot = IAP_domain_structure_WGCNA_hits_exp_upset_plot, filename = "IAP_domain_structure_WGCNA_hits_exp_upset_plot.tiff", device = "tiff",
       width = 20, height = 10,
       path = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/WGCNA/")

### Upset (kinda) plot of domain types (split up) used across experiments
# Assessing if there are patterns in the experiment types where the domains are used
# this table has the number of times each domain type is found in each experiment
# How common are the individual domain name types across experiments if you split up the combo domains within the actual data 
IAP_domain_structure_WGCNA_hits_df_modules_domain_hit_common_exp_separate_EXP <- IAP_domain_structure_WGCNA_hits_df_modules_domain_hit_common_exp_species %>%
  # separate into rows
  separate_rows(comb_domain, sep = ",") %>% 
  # add together the numbers for each type 
  group_by(comb_domain, Species, exp) %>% dplyr::summarize(total_times_across_exp_modules = sum(n))

# Join with experiment type
IAP_domain_structure_WGCNA_hits_df_modules_domain_hit_common_exp_separate_EXP_upset <-  left_join(IAP_domain_structure_WGCNA_hits_df_modules_domain_hit_common_exp_separate_EXP, challenge_type)
# export table to compare with LFC domain usage 
save(IAP_domain_structure_WGCNA_hits_df_modules_domain_hit_common_exp_separate_EXP_upset, file = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/WGCNA/IAP_domain_structure_WGCNA_hits_df_modules_domain_hit_common_exp_separate_EXP_upset.RData")

IAP_domain_structure_WGCNA_hits_df_modules_domain_hit_common_exp_separate_EXP_upset_plot <- ggplot(IAP_domain_structure_WGCNA_hits_df_modules_domain_hit_common_exp_separate_EXP_upset , aes(y=comb_domain, x=exp, fill= total_times_across_exp_modules)) + 
  geom_tile() + facet_grid(.~challenge_type, scales = "free", space = "free") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 20),
        strip.text.x = element_text(size = 18)) +
  scale_fill_viridis_c(option="plasma") +
  labs(x = "Experiment", y = "Domain Name or Combination", title = "Occurence of Domain Types Across Experiments")

ggsave(plot = IAP_domain_structure_WGCNA_hits_df_modules_domain_hit_common_exp_separate_EXP_upset_plot, filename = "IAP_domain_structure_WGCNA_hits_df_modules_domain_hit_common_exp_separate_EXP_upset_plot.tiff", device = "tiff",
       width = 20, height = 10,
       path = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/WGCNA/")


#### Combined table with module stats, apoptosis transcript stats, distinct IAP types #### 



#### Investigation of Interaction Partners ####
## How many apoptosis transcripts are found total across all modules in each experiment
IAP_domain_structure_WGCNA_hits_df_condensed_type %>% ungroup() %>% dplyr::count(exp)
## A tibble: 14 x 2
#exp                    n
#<chr>              <int>
#  1 deLorg_Res           131
#2 deLorg_Sus            95
#3 Dermo_Sus              5
#4 Dermo_Tol             74
#5 He                   160
#6 Pro_RE22_Pro_RI      192
#7 Pro_RE22_Pro_S4      128
#8 Pro_RE22_RE22_full   100
#9 Probiotic             53
#10 ROD_Res                8
#11 Rubio_NV             210
#12 Rubio_V              219
#13 Zhang_LPS            109
#14 Zhang_Vibrio          23

## Number of shared transcripts across particular experiment combos and domain types in all domains - 
IAP_domain_structure_WGCNA_hits_PATHWAY_shared_between_exp_comb_domain <- IAP_domain_structure_WGCNA_hits_df_condensed_type %>% ungroup() %>% 
  dplyr::group_by(transcript_id, comb_domain) %>% dplyr::mutate(exp_combined = paste(exp, collapse = "_")) %>% 
  # get rid of any duplicate rows within experiments so that total count isn't double or triple counting transcripts
  distinct(transcript_id, exp,comb_domain, .keep_all = TRUE) %>%
  # count total shared transcripts in each exp_combined
  ungroup() %>% group_by(exp_combined, comb_domain )%>% dplyr::mutate(count = n()) %>% dplyr::distinct(exp_combined, count, comb_domain) %>%
  arrange(desc(count))

### Compare most common transcripts across modules for each domain type with each challenge type ###
# Split transcript variant info and find other proteins that show up most often as interaction partners ACROSS modules
IAP_domain_structure_WGCNA_hits_freq <- IAP_domain_structure_WGCNA_hits_df_condensed_type %>% 
  separate(product, into = c("product","transcript_variant"), sep = ", transcript") %>% 
  # add column for positive and negative association
  #mutate(module_sign = case_when(
  #  mod_signif >= 0 ~ "positive",   # decided to not separate by positive and negative 
  #  TRUE ~ "negative")) %>% 
  group_by(comb_domain, product, challenge_type) %>%
  dplyr::summarise(product_freq_domain_total = n()) %>% 
  arrange(desc(product_freq_domain_total))
# this gives us the number of times a product is showing up across different modules of a domain type and challenge type

# Create same frequency plot but for experiments rather than type (in case effect is experiment specific)
IAP_domain_structure_WGCNA_hits_freq_exp <- IAP_domain_structure_WGCNA_hits_df_condensed_type %>% 
  separate(product, into = c("product","transcript_variant"), sep = ", transcript") %>% 
  # add column for positive and negative association
  #mutate(module_sign = case_when(
  #  mod_signif >= 0 ~ "positive",   # decided to not separate by positive and negative 
  #  TRUE ~ "negative")) %>% 
  group_by(comb_domain, product, challenge_type, exp) %>%
  dplyr::summarise(product_freq_domain_total = n()) %>% 
  arrange(desc(product_freq_domain_total))

# Create same frequency plot but for domain rather than type 
IAP_domain_structure_WGCNA_hits_freq_domain <- IAP_domain_structure_WGCNA_hits_df_condensed_type %>% 
  separate(product, into = c("product","transcript_variant"), sep = ", transcript") %>% 
  # add column for positive and negative association
  #mutate(module_sign = case_when(
  #  mod_signif >= 0 ~ "positive",   # decided to not separate by positive and negative 
  #  TRUE ~ "negative")) %>% 
  group_by(comb_domain, product) %>%
  dplyr::summarise(product_freq_domain_total = n()) %>% 
  arrange(desc(product_freq_domain_total))

# Create same frequency plot but for transcript_id
IAP_domain_structure_WGCNA_hits_freq_transcript <- IAP_domain_structure_WGCNA_hits_df_condensed_type %>% 
  group_by(comb_domain, transcript_id) %>%
  dplyr::summarise(product_freq_trans_total = n()) %>% 
  arrange(desc(product_freq_trans_total))

# Create list of unique products
IAP_domain_structure_WGCNA_hits_df_condensed_type_product <- IAP_domain_structure_WGCNA_hits_freq %>% ungroup() %>% dplyr::distinct(product)
nrow(IAP_domain_structure_WGCNA_hits_df_condensed_type_product) # 323

# plot frequency of products across modules for particular domain structures 
# heatmap plot
IAP_domain_structure_WGCNA_hits_freq_plot <- IAP_domain_structure_WGCNA_hits_freq %>% 
  filter(product_freq_domain_total > 2) %>% 
  ggplot(aes(x = product, y = comb_domain, fill= product_freq_domain_total)) + geom_tile() + 
  scale_fill_viridis_c(option="plasma") + 
  theme(axis.text.x = element_text(angle = 90, hjust =1, size = 20),
        axis.text.y = element_text(size=20)) + coord_flip()

ggsave(plot = IAP_domain_structure_WGCNA_hits_freq_plot, filename = "IAP_domain_structure_WGCNA_hits_freq_plot.tiff", device = "tiff",
       width = 30, height = 30, limitsize = FALSE,
       path = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/WGCNA/")

## Cluster the above data by profile similarity using pheatmap 
  # Transform the data so that there is one row per transcript, and the columns are the experiment and comb_domain type 
length(unique(IAP_domain_structure_WGCNA_hits_freq_domain $product)) # 323 unique 
# spread the product freq_domain 
IAP_domain_structure_WGCNA_hits_freq_heatmap <- spread(IAP_domain_structure_WGCNA_hits_freq_domain, comb_domain, product_freq_domain_total, fill = 0)
nrow(IAP_domain_structure_WGCNA_hits_freq_heatmap ) # 323
IAP_domain_structure_WGCNA_hits_freq_heatmap <-  column_to_rownames(IAP_domain_structure_WGCNA_hits_freq_heatmap, var = "product") 
IAP_domain_structure_WGCNA_hits_freq_heatmap_mat <- as.matrix(IAP_domain_structure_WGCNA_hits_freq_heatmap)
IAP_domain_structure_WGCNA_hits_freq_heatmap_plot <- pheatmap(IAP_domain_structure_WGCNA_hits_freq_heatmap_mat)

ggsave(plot = IAP_domain_structure_WGCNA_hits_freq_heatmap_plot, filename = "IAP_domain_structure_WGCNA_hits_pheatmap.tiff", device = "tiff",
       width = 30, height = 60, limitsize = FALSE,
       path = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/WGCNA/")

# repeat but don't use frequency to cluster because this can be somewhat biased
# create table where presence of a product is a 1
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod <- IAP_domain_structure_WGCNA_hits_df_condensed_type %>% ungroup() %>%
  separate(product, into = c("product","transcript_variant"), sep = ", transcript") %>% dplyr::distinct(comb_domain, product) %>% mutate(count = 1)
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod <- spread(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod, comb_domain, count, fill = 0)
nrow(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod ) # 323
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod <-  column_to_rownames(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod, var = "product") 
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_mat <- as.matrix(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod)
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_plot <- pheatmap(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_mat)

ggsave(plot = IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_plot, filename = "IAP_domain_structure_WGCNA_hits_product_pheatmap.tiff", device = "tiff",
       width = 30, height = 60, limitsize = FALSE,
       path = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/WGCNA/")

## Make heatmap but make a row for each experiment, domain combo
# create table where presence of a product is a 1
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp <- IAP_domain_structure_WGCNA_hits_df_condensed_type %>% ungroup() %>%
  separate(product, into = c("product","transcript_variant"), sep = ", transcript") %>% dplyr::mutate(dom_exp = paste(exp,comb_domain, sep = ":")) %>%
  dplyr::distinct(dom_exp, product) %>% mutate(count = 1)
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp <- spread(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp, dom_exp, count, fill = 0)
nrow(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp ) # 323
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp <-  column_to_rownames(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp, var = "product") 
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat <- as.matrix(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp)
# make column annotation dataframe 
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot <- data.frame(colnames(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp)) 
colnames(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot)[1] <- "dom_exp"
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot$key <- IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot$dom_exp
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot <- separate(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot,key, into = c("exp", "comb_domain"),sep = ":") %>% 
  mutate(species = case_when(
    exp =="deLorg_Res"  | exp =="deLorg_Sus"|  exp =="He" | exp ==  "Rubio_NV"| exp =="Rubio_V"| exp =="Zhang_LPS"| exp =="Zhang_Vibrio"  ~ "C_gigas",    
    exp == "Dermo_Sus"| exp =="Dermo_Tol"| exp =="Pro_RE22_Pro_RI"| exp =="Pro_RE22_Pro_S4"| exp =="Pro_RE22_RE22_full" | exp =="Probiotic"| exp =="ROD_Res"  ~ "C_virginica",
      TRUE ~ NA_character_))
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot <- column_to_rownames(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot, var = "dom_exp")

# make row annotation dataframe to look at subpathway, join with subpathway
combined_gene_name_org_yes_no_table_unique_pathway_joined_edited <- combined_gene_name_org_yes_no_table_unique_pathway_joined
# change gene_name to product
colnames(combined_gene_name_org_yes_no_table_unique_pathway_joined_edited)[2] <- "product"
# join with frequence table by experiment
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_row <- as.data.frame(rownames(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp))
colnames(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_row)[1] <- "product"
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_row <- left_join(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_row, combined_gene_name_org_yes_no_table_unique_pathway_joined_edited[,c("product","Sub_pathway")])
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_row[is.na(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_row),] # check for NAs
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_row <- column_to_rownames(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_row, var = "product") # make product the rownames 

# generate heatmap
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_plot <- pheatmap(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat, 
                                                                       annotation_col = IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot,
                                                                       annotation_row = IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_row)

ggsave(plot = IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_plot, filename = "IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_pheatmap.tiff", device = "tiff",
       width = 30, height = 60, limitsize = FALSE,
       path = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/WGCNA/")

## Repeat heatmap above but remove "-like" from product names
# create table where presence of a product is a 1
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_like <- IAP_domain_structure_WGCNA_hits_df_condensed_type %>% ungroup() %>%
  filter(!grepl("toll-like", product)) %>% # filter out then add back in so it isn't accidentally removed
  separate(product, into = c("product","transcript_variant"), sep = ", transcript") %>% separate(product, into = c("product","like"), sep = "-like") %>% dplyr::mutate(dom_exp = paste(exp,comb_domain, sep = ":")) %>%
  dplyr::distinct(dom_exp, product) %>% mutate(count = 1)
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_like_toll <- IAP_domain_structure_WGCNA_hits_df_condensed_type %>% ungroup() %>%
  filter(grepl("toll-like", product)) %>% # filter out then add back in so it isn't accidentally removed
  separate(product, into = c("product","transcript_variant"), sep = ", transcript") %>% 
  dplyr::mutate(dom_exp = paste(exp,comb_domain, sep = ":")) %>%
  dplyr::distinct(dom_exp, product) %>% mutate(count = 1) 
#combine back
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_like <- rbind(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_like, IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_like_toll )

IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_like <- spread(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_like, dom_exp, count, fill = 0)
nrow(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_like) # 234
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_like <-  column_to_rownames(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_like, var = "product") 
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_like_mat <- as.matrix(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_like)
# make column annotation dataframe 
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_like <- data.frame(colnames(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_like)) 
colnames(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_like)[1] <- "dom_exp"
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_like$key <- IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_like$dom_exp
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_like <- separate(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_like,key, into = c("exp", "comb_domain"),sep = ":") %>% 
  mutate(species = case_when(
    exp =="deLorg_Res"  | exp =="deLorg_Sus"|  exp =="He" | exp ==  "Rubio_NV"| exp =="Rubio_V"| exp =="Zhang_LPS"| exp =="Zhang_Vibrio"  ~ "C_gigas",    
    exp == "Dermo_Sus"| exp =="Dermo_Tol"| exp =="Pro_RE22_Pro_RI"| exp =="Pro_RE22_Pro_S4"| exp =="Pro_RE22_RE22_full" | exp =="Probiotic"| exp =="ROD_Res"  ~ "C_virginica",
    TRUE ~ NA_character_))
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_like <- column_to_rownames(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_like, var = "dom_exp")

# make row annotation dataframe to look at subpathway, join with subpathway
combined_gene_name_org_yes_no_table_unique_pathway_joined_edited_like <- combined_gene_name_org_yes_no_table_unique_pathway_joined_edited %>% separate(product, into = c("product","like"), sep = "-like") %>%
  distinct(product, Sub_pathway)
# join with frequence table by experiment
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_row_like <- as.data.frame(rownames(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_like))
colnames(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_row_like)[1] <- "product"
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_row_like <- left_join(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_row_like, combined_gene_name_org_yes_no_table_unique_pathway_joined_edited_like[,c("product","Sub_pathway")])
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_row_like[is.na(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_row_like),] # check for NAs
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_row_like <- column_to_rownames(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_row_like, var = "product") # make product the rownames 

# generate heatmap
IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_like_plot <- pheatmap(IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_like_mat, 
                                                                       annotation_col = IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_like,
                                                                       annotation_row = IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_mat_annot_row_like,
                                                                       fontsize = 20)

ggsave(plot = IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_like_plot, filename = "IAP_domain_structure_WGCNA_hits_freq_heatmap_prod_exp_like_pheatmap.tiff", device = "tiff",
       width = 43, height = 60, limitsize = FALSE,
       path = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/WGCNA/")

## Repeat heatmap but with individual transcript_ids not the common names
# Transform the data so that there is one row per transcript, and the columns are the experiment and comb_domain type 
length(unique(IAP_domain_structure_WGCNA_hits_freq_transcript$transcript_id)) # 855 unique 
# create table where presence of a product is a 1
IAP_domain_structure_WGCNA_hits_freq_transcript_heatmap <- IAP_domain_structure_WGCNA_hits_freq_transcript %>% ungroup() %>% dplyr::distinct(transcript_id,comb_domain) %>% 
      dplyr::mutate(value = 1) %>% spread(., comb_domain, value, fill = 0)
nrow(IAP_domain_structure_WGCNA_hits_freq_transcript_heatmap ) # 855 
IAP_domain_structure_WGCNA_hits_freq_transcript_heatmap <-  column_to_rownames(IAP_domain_structure_WGCNA_hits_freq_transcript_heatmap, var = "transcript_id") 
IAP_domain_structure_WGCNA_hits_freq_transcript_heatmap_mat <- as.matrix(IAP_domain_structure_WGCNA_hits_freq_transcript_heatmap)
IAP_domain_structure_WGCNA_hits_freq_transcript_heatmap_plot <- pheatmap(IAP_domain_structure_WGCNA_hits_freq_transcript_heatmap_mat)

ggsave(plot = IAP_domain_structure_WGCNA_hits_freq_transcript_heatmap_plot, filename = "IAP_domain_structure_WGCNA_hits_freq_transcript_pheatmap.tiff", device = "tiff",
       width = 30, height = 60, limitsize = FALSE,
       path = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/WGCNA/")

# how many transcripts are shared between experiments 
# Create same frequency plot but for transcript_id
IAP_domain_structure_WGCNA_hits_df_condensed_type %>% 
  dplyr::distinct(comb_domain, exp, transcript_id) %>% 
  # nrow() # 1507
  group_by(comb_domain, transcript_id) %>%
  filter(n() >1) # 562 total rows
# what percent of total WGCNA transcripts are shared
562/1507*100 # 37.3% are shared, the rest are unique
100-37.3 # 62.7


# Typing my thoughts...
# So what I really want to isolate are those products that are doing different things across domain modules within a type
# products only in specific modules may help better understand specificity to particular stressors  
# the products that show up with all domain types speak more to a ubiquitous response to the stressor
# once I identify particular products and domain structures of interest and can use those to pull out specific domains 

# Goal now is to identify those products with high frequency in some types but not others...
# count how many domain types these products are shared across. Those with less sharing are unique to particular domain structures
IAP_domain_structure_WGCNA_hits_freq_by_domain_challenge_type <-
  IAP_domain_structure_WGCNA_hits_freq %>% 
  group_by(challenge_type, product) %>%
  dplyr::mutate(number_domain_groups_shared = n()) %>%
  distinct() %>%
  arrange(comb_domain, number_domain_groups_shared)

# Find products unique to a particular domain in each type
IAP_domain_structure_WGCNA_hits_freq_by_domain_challenge_type_unique <- IAP_domain_structure_WGCNA_hits_freq_by_domain_challenge_type %>%
  filter(number_domain_groups_shared ==1) %>% arrange(comb_domain)
unique(IAP_domain_structure_WGCNA_hits_freq_by_domain_challenge_type_unique$comb_domain) 

## only distinct transcripts
IAP_domain_structure_WGCNA_hits_freq_by_domain_challenge_type_unique_plot <- IAP_domain_structure_WGCNA_hits_freq_by_domain_challenge_type_unique %>%  # only look at the most common
  ggplot(aes(x = product, y = product_freq_domain_total, fill = comb_domain)) + geom_col() + 
  facet_grid(comb_domain~., scales="free", space ="free") + 
  theme(axis.text.x = element_text(angle = 90, hjust =1)) + coord_flip()

# export
ggsave(plot = IAP_domain_structure_WGCNA_hits_freq_by_domain_challenge_type_unique_plot, filename = "IAP_domain_structure_WGCNA_hits_unique_per_domain.tiff", device = "tiff",
       width = 30, height = 25,
       path = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/WGCNA")

### Are there any modules that only have IAP domain structures of one type? in other words: do any modules for same experiment have different domain ID
IAP_domain_structure_WGCNA_hits_exp_unique <- IAP_domain_structure_WGCNA_hits_df_condensed %>%
   filter(count == 1) %>% distinct(mod_names, exp, comb_domain)
IAP_domain_structure_WGCNA_hits_exp_unique_products_count <- left_join(IAP_domain_structure_WGCNA_hits_exp_unique, IAP_domain_structure_WGCNA_hits_exp)

# those modules in this list that only have one domain type in a module may reveal specific pathways triggered by the domain types
IAP_domain_structure_WGCNA_hits_exp_unique_cytoscape <- IAP_domain_structure_WGCNA_hits_exp_unique_products_count %>%
  filter(products_per_mod_exp >=10) # 11 of the modules have greater than 10 products and only 1 domain type, but none are viral or parasite! 

### ANALYSIS OF PATHWAYS ASSOCIATED WITH SPECIFIC DOMAIN STRUCTURES ###
# use sorted table with the product frequency with each domain type
IAP_domain_structure_WGCNA_hits_freq %>% arrange(product, product_freq_domain_total) %>% View()

## Investigate Caspases
IAP_domain_structure_WGCNA_hits_freq_caspase <- IAP_domain_structure_WGCNA_hits_freq %>% filter(grepl("caspase", product) & !(grepl("activity", product)))

# plot domain usage across caspase
ggplot(IAP_domain_structure_WGCNA_hits_freq_caspase, aes(x= product, y = product_freq_domain_total, fill=comb_domain)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust=1)) + coord_flip() + 
  facet_grid(.~challenge_type)
ggplot(IAP_domain_structure_WGCNA_hits_freq_caspase, aes(x= comb_domain, y = product_freq_domain_total, fill=product)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust=1)) + coord_flip()+ 
  facet_grid(.~challenge_type)

## PATHWAY SPECIFIC DOMAIN STRUCTURES 
# ADD IN PATHWAY FOR PLOTTING
IAP_domain_structure_WGCNA_hits_freq_caspase_pathway <- IAP_domain_structure_WGCNA_hits_freq_caspase %>% 
  mutate(pathway = case_when(
    product == "caspase-1-like"  ~ "inflammatory",
    product == "caspase-2"      ~ "DNA_damage",
    product == "caspase-2-like"  ~ "DNA_damage",
    product == "caspase-3-like"  ~ "execution",
    product == "caspase-6"       ~ "execution",
    product == "caspase-6-like"  ~ "execution",
    product == "caspase-7"       ~ "execution",
    product == "caspase-7-like"  ~ "execution",
    product == "caspase-8"       ~ "extrinsic",
    product == "caspase-8-like"  ~ "extrinsic",
    product == "caspase-9-like" ~ "intrinsic"
  ))
# plot 
IAP_domain_structure_WGCNA_hits_freq_caspase_pathway_plot <- ggplot(IAP_domain_structure_WGCNA_hits_freq_caspase_pathway, aes(x= comb_domain, y = product_freq_domain_total, fill=pathway)) + geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust =1, size = 20),
        axis.text.y = element_text(size=20),
        legend.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.x = element_text(size = 25)) + coord_flip() + 
  labs(y = "Caspase Transcript Number", x = "Domain Structure Type") + 
  facet_grid(.~challenge_type) # chan choose to split it up or not
# export
ggsave(plot = IAP_domain_structure_WGCNA_hits_freq_caspase_pathway_plot , filename = "IAP_domain_structure_WGCNA_hits_freq_caspase_pathway_plot.tiff", device = "tiff",
       width = 25, height = 20,
       path = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/WGCNA")

## OVERALL OBSERVATIONS: 

# PATHWAY SPECIFIC CASPASES:
# A. EXTRINSIC: Caspase 8 domain structure types:  involves ALL BUT "TII"          "TII-BIR6-E2"  "TII-TII-RING"
IAP_domain_structure_WGCNA_hits_freq_caspase8_dom_name <-  IAP_domain_structure_WGCNA_hits_freq_caspase %>% filter(grepl("caspase-8",product)) %>% distinct(Domain_Name)
setdiff(unique(IAP_domain_structure_WGCNA_hits_df_type$Domain_Name), IAP_domain_structure_WGCNA_hits_freq_caspase8_dom_name$Domain_Name) # 
#TX-TII              caspase-8-like
#2 TII-DD              caspase-8-like
#3 TI-TII-DD-RING      caspase-8     
#4 TI-TII-DD-RING      caspase-8-like
#5 TI-TII-TII-UBA-RING caspase-8-like
#6 TII-DD-RING         caspase-8     
#7 TII-TII             caspase-8     
#8 TX-TII              caspase-8     
#9 BIR*-DD-RING        caspase-8-like
#10 TII-DD              caspase-8     
#11 TII-DD-RING         caspase-8-like
#12 TII-RING            caspase-8-like 

# B. INTRINSIC: Caspase 9 domain structure types: NO UBA repeats involved, no TI repeats involved
IAP_domain_structure_WGCNA_hits_freq_caspase9 <-  IAP_domain_structure_WGCNA_hits_freq_caspase %>% filter(grepl("caspase-9",product)) %>% distinct(Domain_Name)
setdiff(unique(IAP_domain_structure_WGCNA_hits_df_type$Domain_Name), IAP_domain_structure_WGCNA_hits_freq_caspase9$Domain_Name) #
  # missing: "TI-TII-TII-UBA-RING" "TII"                 "TII-BIR6-E2"         "TII-DD"              "TII-TII"             "TII-TII-RING" 
#BIR*-DD-RING   caspase-9-like
#2 TII-DD-RING    caspase-9-like
#3 TI-TII-DD-RING caspase-9-like
#4 TII-RING       caspase-9-like
#5 TX-TII         caspase-9-like

# C. DNA damage: caspase 2 
IAP_domain_structure_WGCNA_hits_freq_caspase2 <-  IAP_domain_structure_WGCNA_hits_freq_caspase %>% filter(grepl("caspase-2",product)) %>% distinct(Domain_Name)
setdiff(unique(IAP_domain_structure_WGCNA_hits_df_type$Domain_Name), IAP_domain_structure_WGCNA_hits_freq_caspase2$Domain_Name)   
  # missing: "BIR*-DD-RING" "TII"          "TII-BIR6-E2"  "TII-DD"       "TII-TII-RING"
# 1 TI-TII-DD-RING      caspase-2-like
# 2 TI-TII-TII-UBA-RING caspase-2     
# 3 TII-DD-RING         caspase-2-like
# 4 TII-RING            caspase-2-like
# 5 TII-TII             caspase-2-like
# 6 TX-TII              caspase-2-like

# D. Inflammatory: caspase 1: all domain structures involved except TII-DD can ubiquitinate, interesting pattern
IAP_domain_structure_WGCNA_hits_freq_caspase1 <-  IAP_domain_structure_WGCNA_hits_freq_caspase %>% filter(grepl("caspase-1",product)) %>% distinct(Domain_Name)
setdiff(unique(IAP_domain_structure_WGCNA_hits_df_type$Domain_Name), IAP_domain_structure_WGCNA_hits_freq_caspase9$Domain_Name) #"NZBIR-TII-UBA-DD-RING", "TI-TII-DD-RING","TI-TII-TII-UBA-RING", "TII","TII-DD-RING","TII-RING","TII-TII-RING"  
#1 TI-TII-DD-RING      caspase-1-like
#2 TI-TII-TII-UBA-RING caspase-1-like     

### TLR Analysis ###
IAP_domain_structure_WGCNA_hits_freq_TLR <- IAP_domain_structure_WGCNA_hits_freq %>% filter(grepl("toll", product))

# plot domain usage across caspase
ggplot(IAP_domain_structure_WGCNA_hits_freq_TLR, aes(x= product, y = product_freq_domain_total, fill=comb_domain)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust=1)) + coord_flip() + 
  facet_grid(.~challenge_type)
IAP_domain_structure_WGCNA_hits_freq_TLR_all_plot <- ggplot(IAP_domain_structure_WGCNA_hits_freq_TLR, aes(x= comb_domain, y = product_freq_domain_total, fill=product)) + 
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust =1, size = 20),
    axis.text.y = element_text(size=20),
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    strip.text.x = element_text(size = 25)) + 
  coord_flip() + 
  labs(y = "Caspase Transcript Number", x = "Domain Structure Type") + 
  facet_grid(.~challenge_type)

# export
ggsave(plot = IAP_domain_structure_WGCNA_hits_freq_TLR_all_plot , filename = "IAP_domain_structure_WGCNA_hits_freq_TLR_all_plot.tiff", device = "tiff",
       width = 25, height = 20,
       path = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/WGCNA")

# split into intracellular and extracellular localization 
# (Takeda and Yamamoto 2010): "TLRs can be divided into extracellular and intracellular TLRs. TLR1, TLR2, TLR4, TLR5, TLR6, and TLR11 recognize their ligands on the cell surface. 
# On the other hand, TLR3, TLR7, TLR8, and TLR9 are intracellularly localized"

IAP_domain_structure_WGCNA_hits_freq_TLR_type <- IAP_domain_structure_WGCNA_hits_freq_TLR %>%
  mutate(type = case_when(
    product == "toll-like receptor 1"~ "extracellular",
    product == "toll-like receptor 10"~ "extracellular", #(uniprot says membrane bound)
    product == "toll-like receptor 13"~ "intracellular", # in the endosome membrane, binds RNA
    product == "toll-like receptor 2"~ "extracellular",
    product == "toll-like receptor 2 type-1"~ "extracellular",
    product == "toll-like receptor 2 type-2"~ "extracellular",
    product == "toll-like receptor 3"     ~ "intracellular",
    product == "toll-like receptor 4"~ "extracellular",
    product == "toll-like receptor 6"~ "extracellular",
    product == "toll-like receptor 7"~ "intracellular",
    product == "toll-like receptor 8"~ "intracellular",
    product == "toll-like receptor 9" ~ "intracellular",
    product == "toll-like receptor Tollo"   ~ "intracellular" # toll 8 is AKA tollo (Akhouayri et al., 2011)
  )) 

# fill with intracellular vs. extracellular
IAP_domain_structure_WGCNA_hits_freq_TLR_type_plot <- ggplot(IAP_domain_structure_WGCNA_hits_freq_TLR_type, aes(x= comb_domain, y = product_freq_domain_total, fill=type)) +
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust =1, size = 20),
        axis.text.y = element_text(size=20),
        legend.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.x = element_text(size = 25)) + 
  coord_flip() + 
  labs(y = "Caspase Transcript Number", x = "Domain Structure Type") + 
  facet_grid(.~challenge_type)

# export
ggsave(plot = IAP_domain_structure_WGCNA_hits_freq_TLR_type_plot, filename = "IAP_domain_structure_WGCNA_hits_freq_TLR_type_plot.tiff", device = "tiff",
       width = 25, height = 20,
       path = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/WGCNA")

# OBSERVATIONS
# 1. DOMAIN TYPE SPECIFC TLRS 
# - TLR 8 (intracellular): specific to TII and bacterial challenge only
# - TLR 10 (extracellular: specific to TII-TII
# 2. TLRs involved in response to all domain types
# 3. Split into intracellular and extracellular - *ALL DOMAIN STRUCTURE TYPES ASSOCIATED WITH BOTH INTRACELLULAR AND EXTRACELLULAR TLRS*
# - * ALL EXPERIMENT TYPES ASSOCIATED WITH BOTH INTRACELLULAR AND EXTRACELLULAR TLRS* 
# intracellular vs extracellular doesnt seem to make a difference

### Bcl-2 family analysis ###
IAP_domain_structure_WGCNA_hits_freq_bcl2 <- IAP_domain_structure_WGCNA_hits_freq %>% filter(grepl("BAX", product) |
                                                                                               grepl("bcl-2-related protein A1", product) | grepl("B-cell translocation gene 1",product) | 
                                                                                               grepl("bcl-2 homologous antagonist/killer", product) | grepl("bcl-2-like protein 1",product))

# plot domain usage across caspase
ggplot(IAP_domain_structure_WGCNA_hits_freq_bcl2 , aes(x= product, y = product_freq_domain_total, fill=Domain_Name)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust=1)) + coord_flip() + 
  facet_grid(.~challenge_type)
ggplot(IAP_domain_structure_WGCNA_hits_freq_bcl2 , aes(x= Domain_Name, y = product_freq_domain_total, fill=product)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust=1)) + coord_flip()+ 
  facet_grid(.~challenge_type)

# which domain structures missing
setdiff(unique(IAP_domain_structure_WGCNA_hits_df_type$Domain_Name), IAP_domain_structure_WGCNA_hits_freq_bcl2$Domain_Name) # "BIR*-DD-RING" "NZBIR-TII-UBA-DD-RING" "TI-TII-TII-UBA-RING" 

## Observations:
# Goal is to look for consistent patterns across experiment challenge in molecules associated 
# 1. TI-TII-DD-RING maintains the most diversity across experimental challenges and associated with bcl-1, A1, bcl-xl and bak across experiments
# main thing this tells me is that bcl2 family molecules are not associated with "BIR*-DD-RING" "NZBIR-TII-UBA-DD-RING" "TI-TII-TII-UBA-RING". Possible the UBA is only 
# associated with extracellular pathways?

### MAPK analysis ###
IAP_domain_structure_WGCNA_hits_freq_MAPK <- IAP_domain_structure_WGCNA_hits_freq %>% filter(grepl("mitogen-activated protein", product))

# plot domain usage across caspase
ggplot(IAP_domain_structure_WGCNA_hits_freq_MAPK , aes(x= product, y = product_freq_domain_total, fill=comb_domain)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust=1)) + coord_flip() + 
  facet_grid(.~challenge_type)
ggplot(IAP_domain_structure_WGCNA_hits_freq_MAPK , aes(x= comb_domain, y = product_freq_domain_total, fill=product)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust=1)) + coord_flip()+ 
  facet_grid(.~challenge_type)

# which domain structures missing
setdiff(unique(IAP_domain_structure_WGCNA_hits_df_type$comb_domain), IAP_domain_structure_WGCNA_hits_freq_MAPK$Domain_Name) # "TII-TII-RING"

## Observations:
# TII-TII-RING is only domain structure type not associated with MAPK

#### TNFR Analysis ###
IAP_domain_structure_WGCNA_hits_freq_TNFR <- IAP_domain_structure_WGCNA_hits_freq %>% filter(grepl("tumor necrosis", product) & !grepl("lipopolysaccharide",product) &
                                                                                               !grepl("alpha-induced",product))
# plot domain usage across caspase
ggplot(IAP_domain_structure_WGCNA_hits_freq_TNFR , aes(x= product, y = product_freq_domain_total, fill=comb_domain)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust=1)) + coord_flip() + 
  facet_grid(.~challenge_type)
ggplot(IAP_domain_structure_WGCNA_hits_freq_TNFR , aes(x= comb_domain, y = product_freq_domain_total, fill=product)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust=1)) + coord_flip()+ 
  facet_grid(.~challenge_type)

# which domain structures missing
setdiff(unique(IAP_domain_structure_WGCNA_hits_df_type$Domain_Name), IAP_domain_structure_WGCNA_hits_freq_TNFR$Domain_Name) # "TII-TII-RING"

## Observations 
# TNFSF5 unique to TII-DD, TII-BIR6-E2, TI-TII-DD-RING. 
# Domain specific observations:
# TII-BIR6-E2 only associated with TNFSF5
### this strategy of looking at particular molecules is not very helpful! 


### ANALYSIS OF IAP TRANSCRIPT USAGE ####
# Overall question is whether in shared pathways with similar domain structures, are the same or different IAP transcripts being used 

# How many IAPs hits in each of these modules
IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP <- IAP_domain_structure_WGCNA_hits_df_condensed_type[IAP_domain_structure_WGCNA_hits_df_condensed_type$transcript_id %in% IAP_domain_structure_XM_filter$transcript_id,]

# how many IAPs identified total 
nrow(IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP) # 156 

# how many unique IAPS
length(unique(IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP$transcript_id)) #96 total 

# plot IAP counts in each individual module instead of module counts above
IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP_count <- IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP %>% group_by(mod_names, exp, comb_domain) %>% dplyr::mutate(IAP_count = n()) %>%
  dplyr::distinct(mod_names, exp, comb_domain,challenge_type, IAP_count, comb_domain_type) 

IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP_plot <- 
  ggplot(IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP_count, aes(y=comb_domain, x=exp, fill= IAP_count)) + geom_tile() +
  facet_grid(.~challenge_type, scales = "free", space = "free") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 20),
        strip.text.x = element_text(size = 18)) +
  scale_fill_viridis_c(option="plasma") + 
  labs(x = "Experiment", y = "Module Classification", title = "IAP Counts Modules Across Experiments", fill = "IAP transcripts\nPer Exp.") 

ggsave(plot = IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP_plot, filename = "IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP_plot.tiff", device = "tiff",
       width = 20, height = 10,
       path = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/WGCNA/")

## What percent of IAPs are in unique vs. combo domains?
IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP_count %>% ungroup() %>% dplyr::group_by(comb_domain_type) %>% dplyr::summarize(count = sum(IAP_count)) %>%
  ungroup() %>% dplyr::mutate(percent_total = count/sum(count)*100)
#comb_domain_type count percent_total
#<chr>            <int>         <dbl>
#  1 combo              115          73.7
#2 unique              41          26.3

# Are any IAP transcripts shared within domain type and between experiments?
IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP_shared <- x   %>% ungroup() %>% 
  # remove rows with duplicated module name and shared_exp, these are those where the module membership was the same
  group_by(mod_names, shared_exp) %>% filter(n()==1) %>% ungroup() %>% distinct(comb_domain, transcript_id, shared_exp)
# 3 total shared IAPs between all WGCNA modules
#A tibble: 3 x 3
#comb_domain           transcript_id  shared_exp               
#<chr>                 <chr>          <chr>                    
#1 TII-DD-RING           XM_022438072.1 Dermo_Tol_Pro_RE22_Pro_S4
#2 NZBIR-TII-UBA-DD-RING XM_020067151.1 He_Zhang_LPS             
#3 not_classified        XM_011416128.2 Rubio_V_Zhang_LPS        

# Are IAP transcripts shared across experiments, regardless of domain type
# remove grouping by comb_domain and just look at sharing of experiments
IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP_shared_all <- IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP %>% group_by(transcript_id) %>% filter(n()>1) %>% 
  group_by(transcript_id) %>%  dplyr::mutate(shared_exp = paste(exp, collapse = "_")) %>% ungroup() %>% 
  # remove rows with duplicated module name and shared_exp, these are those where the module membership was the same
  group_by(mod_names, shared_exp) %>% filter(n()==1) %>% ungroup() %>% distinct(transcript_id, shared_exp)
# 20  transcripts are shared
# what is their domain type
IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP_shared_all_type <- left_join(IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP_shared_all, IAP_domain_structure_XM_filter)
IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP_shared_all_plot <- 
  IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP_shared_all_type %>% ungroup() %>% dplyr::mutate(count=1) %>% 
  ggplot(aes(x= shared_exp, y = count, fill = Domain_Name)) + geom_col() + coord_flip()

IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP_shared_all_type_domain_plot <-   IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP_shared_all_type %>% ungroup() %>% dplyr::mutate(count=1) %>% 
  ggplot(aes(x= Domain_Name, y = count, fill = shared_exp)) + geom_col() + coord_flip()

ggsave(plot = IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP_shared_all_plot, filename = "IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP_shared_all_plot.tiff", device = "tiff",
       width = 10, height = 10,
       path = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/WGCNA/")

ggsave(plot = IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP_shared_all_type_domain_plot, filename = "IAP_domain_structure_WGCNA_hits_df_condensed_type_IAP_shared_all_type_domain_plot.tiff", device = "tiff",
       width = 10, height = 10,
       path = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/WGCNA/")




#### PATHWAY COMPARISON BETWEEN SIGNIFICANT WGCNA MODULES ####

## note in above IAP_domain_structure_WGCNA_hits_freq and IAP_domain_structure_WGCNA_hits_freq_exp, grouping was done by name
## grouping by name gets rid of looking at transcript diversity

# Create same frequency plot but for experiments rather than type (in case effect is experiment specific)
IAP_domain_structure_WGCNA_hits_transcript_id_freq_exp <- IAP_domain_structure_WGCNA_hits_df_condensed_type %>% 
  separate(product, into = c("product","transcript_variant"), sep = ", transcript") %>% 
  # add column for positive and negative association
  #mutate(module_sign = case_when(
  #  mod_signif >= 0 ~ "positive",   # decided to not separate by positive and negative 
  #  TRUE ~ "negative")) %>% 
  group_by(comb_domain, transcript_id, product, challenge_type, exp) %>%
  dplyr::summarise(product_freq_domain_total = n()) %>% 
  arrange(desc(product_freq_domain_total))

## Join product frequency tables with products gene_name_pathway_key_merged loaded at the top of script
# Note that in my orignal script to create this table I got rid of "like" for products so I could better compare which might be truly missing from one or the other
combined_gene_name_org_yes_no_table_unique_pathway_joined_edited <- combined_gene_name_org_yes_no_table_unique_pathway_joined
# change gene_name to product
colnames(combined_gene_name_org_yes_no_table_unique_pathway_joined_edited)[2] <- "product"

# join with frequence table by experiment
IAP_domain_structure_WGCNA_hits_freq_exp_PATHWAY <- left_join(IAP_domain_structure_WGCNA_hits_transcript_id_freq_exp, combined_gene_name_org_yes_no_table_unique_pathway_joined_edited[,c("product","Sub_pathway")])

# join with original table, not condensed by frequency
IAP_domain_structure_WGCNA_hits_PATHWAY <- IAP_domain_structure_WGCNA_hits_df_condensed_type %>% 
  separate(product, into = c("product","transcript_variant"), sep = ", transcript") %>% 
  left_join(., combined_gene_name_org_yes_no_table_unique_pathway_joined_edited[,c("product","Sub_pathway")]) %>% dplyr::select(-transcript_variant)

# check for NAs
IAP_domain_structure_WGCNA_hits_freq_exp_PATHWAY %>% filter(is.na(Sub_pathway)) %>% View() # no NAs all have been assigned a pathway
IAP_domain_structure_WGCNA_hits_PATHWAY %>% filter(is.na(Sub_pathway)) %>% View() # no NAs all have been assigned a pathway


### QUESTION: WHICH PATHWAYS ARE UNIQUELY INVOLVED WITH INDIVIDUAL IAP DOMAIN TYPES: EXAMINE PATHWAYS INVOLVED IN MODULES WITH ONLY ONE IAP DOMAIN TYPE ###
IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain <- IAP_domain_structure_WGCNA_hits_PATHWAY %>% filter(comb_domain_type == "unique")
IAP_domain_structure_WGCNA_hits_PATHWAY_comb_domain <- IAP_domain_structure_WGCNA_hits_PATHWAY %>% filter(comb_domain_type == "combo")
unique(IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain$comb_domain) 
#[1] "TII-TII-RING"          "TII-BIR6-E2"           "BIR*"                  "not_classified"        "TII-DD-RING"           "NZBIR-TII-UBA-DD-RING" "TI-TII-RING"           "TX-TII"               
#[9] "TII-DD"                "BIR*-DD-RING"          "TI-TII-TII-UBA-RING"   "TI-TII-DD-RING"  

# View types one at a time
IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain %>% filter(comb_domain == "TII-TII-RING") %>% View()
IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain %>% filter(comb_domain == "TII-BIR6-E2") %>% View()
IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain %>% filter(comb_domain == "TI-TII-RING") %>% View()
IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain %>% filter(comb_domain == "BIR*" ) %>% View()
IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain %>% filter(comb_domain == "TII-DD-RING") %>% View()
IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain %>% filter(comb_domain == "NZBIR-TII-UBA-DD-RING") %>% View()
IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain %>% filter(comb_domain == "TI-TII-RING" ) %>% View()
IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain %>% filter(comb_domain == "TX-TII") %>% View()
IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain %>% filter(comb_domain == "TII-DD") %>% View()
IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain %>% filter(comb_domain == "BIR*-DD-RING") %>% View()
IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain %>% filter(comb_domain == "TI-TII-TII-UBA-RING" ) %>% View()
IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain %>% filter(comb_domain == "TI-TII-DD-RING") %>% View()

# Which transcripts shared between pathways in unique domains?
IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain_IAP_shared <- IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain %>% ungroup()  %>% dplyr::count(transcript_id, product, exp, challenge_type, Sub_pathway)
ggplot(IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain_IAP_shared, aes(x = transcript_id, y = n, fill = Sub_pathway)) + geom_col() +
  theme(axis.text.x = element_text(angle = 90, hjust  =1 )) + facet_grid(.~challenge_type) 
ggplot(IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain_IAP_shared, aes(x = transcript_id, y = n, fill = exp)) + geom_col() +
  theme(axis.text.x = element_text(angle = 90, hjust  =1 )) + facet_grid(.~challenge_type)

# Number of shared transcripts across particular experiment combos in unique domains- 
IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain_IAP_shared_comb <- IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain_IAP_shared %>% ungroup() %>% 
  dplyr::group_by(transcript_id) %>% dplyr::mutate(exp_combined = paste(exp, collapse = "_")) %>% 
  # get rid of any duplicate rows within experiments so that total count isn't double or triple counting transcripts
  distinct(transcript_id, exp, .keep_all = TRUE) %>%
  # count total shared transcripts in each exp_combined
  ungroup() %>% group_by(exp_combined )%>% dplyr::mutate(count = n()) %>% dplyr::distinct(exp_combined, count) %>%
  arrange(desc(count))

# exp_combined                    count
# <chr>                           <int>
#   1 Rubio_NV_Rubio_V                   94
# 2 Probiotic                          50
# 3 He_Rubio_NV_Rubio_V                45
# 4 He                                 39
# 5 Pro_RE22_Pro_S4                    25
# 6 Pro_RE22_Pro_RI                    22
# 7 Zhang_Vibrio                       20
# 8 Zhang_LPS                          19
# 9 Dermo_Tol                          18
# 10 Pro_RE22_Pro_RI_Pro_RE22_Pro_S4    10

## QUESTION: ARE THERE SPECIFIC INTERACTION PARTNERS ASSOCIATED ONLY WITH A PARTICULAR DOMAIN TYPE
# Which products are only found once in the unique domain groups?
# correct first for those modules with duplcate 
dplyr::distinct(IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain[,c("mod_names","exp")]) %>% View()
  # MEbrown, MEred, MEyellow Rubio_NV and V have same module membership, should only keep 1 of each of these 
  # MEsteelblue is the same for Pro_RI and Pro_S4

IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain_uniq_product <- IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain %>% 
  ungroup() %>% group_by(product) %>% filter(n()==1) # doing it this way ensures you're only getting things that occur once 
nrow(IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain_uniq_product)
### Need to fix this because for experiments like PRo_re22 and Rubio NV and V the modules have the exact same members, dont want to remove these - August 24th, 2020 decided to keep this in 

# Are these unique molecules for each domain type also found in their respective groups with combo domains? Use semi_join so I get all rows that a match
IAP_domain_structure_WGCNA_hits_PATHWAY_comb_domain_uniq_prod <- inner_join(IAP_domain_structure_WGCNA_hits_PATHWAY_unique_domain_uniq_product[,c("product","transcript_id", "comb_domain")], IAP_domain_structure_WGCNA_hits_PATHWAY_comb_domain, by = c("product","transcript_id")) 
# which of these found transcripts share the same comb_domain involvement: 7 transcripts
        # - tumor necrosis factor receptor superfamily member 5	XM_011458363.2	TI-TII-DD-RING	MEturquoise	deLorg_Sus
        # - mitogen-activated protein kinase kinase kinase 13-A	XM_011450976.2	TI-TII-DD-RING	MEturquoise	deLorg_Sus
        # - NF-kappa-B inhibitor epsilon-like	NM_001308876.1	TI-TII-DD-RING	MEturquoise	deLorg_Res  
        # - baculoviral IAP repeat-containing protein 6		XM_011438295.2	TII-BIR6-E2	MEpurple	He
        # - uncharacterized LOC111122723	XM_022464605.1	BIR*	MEroyalblue	Pro_RE22_Pro_S4
        # - dynamin-like 120 kDa protein, mitochondrial	XM_022486900.1	BIR*	MEroyalblue	Pro_RE22_Pro_S4
        # - programmed cell death protein 6-like	XM_022483983.1	BIR*	MEturquoise	Pro_RE22_RE22_full

## Are there any unique when considering all domain groups and combo groups
IAP_domain_structure_WGCNA_hits_PATHWAY_unique <- IAP_domain_structure_WGCNA_hits_PATHWAY %>% group_by(product) %>% filter(n() ==1)
nrow(IAP_domain_structure_WGCNA_hits_PATHWAY_unique) # 65 are unique 

## List of only the IAP transcripts in full table based on combined IAP list IAP_domain_structure_XM_CV_XM
IAP_domain_structure_WGCNA_hits_PATHWAY_IAP_only <-IAP_domain_structure_WGCNA_hits_PATHWAY[(IAP_domain_structure_WGCNA_hits_PATHWAY$transcript_id %in% IAP_domain_structure_XM_CV_XM$transcript_id),]

#### ANALYSIS OF SIGNIFICANT IAP INTERACTIONS WITH FULL WGCNA DATA RUN FOR EACH EXPERIMENT SEPARATELY #### 
C_vir_C_gig_full_all_exp_mod_sig_apop 
  
# Find modules for each experiment that are also significant in their DEG experiments
# Load DESeq2 dataframe that was recoded for collapsed haplotigs in the IAP_GIMAP_Gene_Expansion.R script
load(file = "/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/WGCNA/LFC_cont_comb_domain_type.RData")

LFC_comb_domain_type <- LFC_cont_comb_domain_type %>% filter(Data_Type == "LFC")

# join with transcript info 
LFC_comb_domain_type_CV_XM <- LFC_comb_domain_type %>% filter(Species == "Crassostrea_virginica") %>% left_join(., IAP_domain_structure_XM_CV_XM[,c("protein_id","transcript_id")])
LFC_comb_domain_type_CG_XM <- LFC_comb_domain_type %>% filter(Species == "Crassostrea_gigas") %>% left_join(., IAP_domain_structure_XM_CG[,c("protein_id","transcript_id")], by ="protein_id")
LFC_comb_domain_type_XM <- rbind(LFC_comb_domain_type_CV_XM, LFC_comb_domain_type_CG_XM) 

# check for non matched XMs
is.na(LFC_comb_domain_type_XM$transcript_id) # none
# create search lists for Experiment
LFC_comb_domain_type_XM_df <- LFC_comb_domain_type_XM %>%
  group_by(experiment) %>%
  # put all proteins into a string
  dplyr::summarise(transcript_id = paste0(transcript_id, collapse = ",")) %>%
  # set column names
  column_to_rownames(., var= "experiment")

# split strings into vector for grepping
LFC_comb_domain_type_XM_df$transcript_id <- strsplit(LFC_comb_domain_type_XM_df$transcript_id, ",")

# create list of lists
LFC_comb_domain_type_XM_df_list <- as.list(as.data.frame(t(LFC_comb_domain_type_XM_df)))

# loop through to get list of module hits for each group of transcripts
LFC_comb_domain_type_XM_exp_hits <- vector('list', length(LFC_comb_domain_type_XM_df_list))
names(LFC_comb_domain_type_XM_exp_hits) <- names(LFC_comb_domain_type_XM_df_list) # set names
for(i in seq_along(LFC_comb_domain_type_XM_df_list)){
  for(j in seq_along(LFC_comb_domain_type_XM_df_list[[i]])){
    LFC_comb_domain_type_XM_exp_hits[[i]][[j]]<- C_vir_C_gig_full_all_exp_mod_sig_apop %>% 
      group_by(mod_names, exp) %>%
      filter(any(grepl(paste(LFC_comb_domain_type_XM_df_list[[i]][[j]], collapse="|"), transcript_id)))
  }
}
# put in dataframe (using purrr yay!)
LFC_comb_domain_type_XM_exp_hits_df <- map_df(LFC_comb_domain_type_XM_exp_hits, ~bind_rows(., .id="DESeq_experiment"), .id="DESeq_experiment") 

# subset for each experiments
levels(factor(LFC_comb_domain_type_XM_exp_hits_df$DESeq_experiment))
#[1] "de Lorgeril<br>Res. OsHV-1" "de Lorgeril<br>Sus. OsHV-1" "Dermo"                      "Hatchery\nPro. RI"          "He OsHV-1"                  "Lab Pro. S4, RI\n or RE22" 
#[7] "ROD"                        "Rubio<br>*Vibrio* spp."     "Zhang<br>*Vibrio* spp." 

# Subset for rows where the DESeq experiment is the same as the WGCNA experiment - meaning shared IAPs
LFC_comb_domain_type_XM_exp_hits_df_deLorg_Res <- LFC_comb_domain_type_XM_exp_hits_df %>% filter(grepl("de Lorgeril<br>Res. OsHV-1", DESeq_experiment)) %>% filter(grepl("deLorg_Res", exp))
LFC_comb_domain_type_XM_exp_hits_df_deLorg_Sus <- LFC_comb_domain_type_XM_exp_hits_df %>% filter(grepl("de Lorgeril<br>Sus. OsHV-1", DESeq_experiment)) %>% filter(grepl("deLorg_Sus", exp))
LFC_comb_domain_type_XM_exp_hits_df_Dermo <- LFC_comb_domain_type_XM_exp_hits_df %>% filter(grepl("Dermo", DESeq_experiment)) %>% filter(grepl("Dermo", exp))
LFC_comb_domain_type_XM_exp_hits_df_Probiotic <- LFC_comb_domain_type_XM_exp_hits_df %>% filter(grepl("Hatchery\nPro. RI", DESeq_experiment)) %>% filter(grepl("Probiotic", exp))
LFC_comb_domain_type_XM_exp_hits_df_Pro_RE22 <- LFC_comb_domain_type_XM_exp_hits_df %>% filter(grepl("Lab Pro. S4, RI\n or RE22", DESeq_experiment)) %>% filter(grepl("Pro_", exp))
LFC_comb_domain_type_XM_exp_hits_df_He <- LFC_comb_domain_type_XM_exp_hits_df %>% filter(grepl("He OsHV-1", DESeq_experiment)) %>% filter(grepl("He", exp))
LFC_comb_domain_type_XM_exp_hits_df_ROD <- LFC_comb_domain_type_XM_exp_hits_df %>% filter(grepl("ROD", DESeq_experiment)) %>% filter(grepl("ROD", exp))
LFC_comb_domain_type_XM_exp_hits_df_Rubio <- LFC_comb_domain_type_XM_exp_hits_df %>% filter(grepl("Rubio", DESeq_experiment)) %>% filter(grepl("Rubio", exp))
LFC_comb_domain_type_XM_exp_hits_df_Zhang <- LFC_comb_domain_type_XM_exp_hits_df %>% filter(grepl("Zhang", DESeq_experiment)) %>% filter(grepl("Zhang", exp))

# Find Shared IAPs between the two experiments
LFC_comb_domain_type_XM_exp_hits_df_deLorg_Res_IAP <- left_join(LFC_comb_domain_type_XM_exp_hits_df_deLorg_Res, LFC_comb_domain_type_XM[,c("Domain_Name","transcript_id")]) %>% filter(!is.na(Domain_Name))
LFC_comb_domain_type_XM_exp_hits_df_deLorg_Sus_IAP <- left_join(LFC_comb_domain_type_XM_exp_hits_df_deLorg_Sus, LFC_comb_domain_type_XM[,c("Domain_Name","transcript_id")]) %>% filter(!is.na(Domain_Name))
LFC_comb_domain_type_XM_exp_hits_df_Dermo_IAP <- left_join(LFC_comb_domain_type_XM_exp_hits_df_Dermo, LFC_comb_domain_type_XM[,c("Domain_Name","transcript_id")]) %>% filter(!is.na(Domain_Name))
#LFC_comb_domain_type_XM_exp_hits_df_Probiotic_IAP  <- left_join(LFC_comb_domain_type_XM_exp_hits_df_Probiotic, LFC_comb_domain_type_XM[,c("Domain_Name","transcript_id")]) %>% filter(!is.na(Domain_Name))
LFC_comb_domain_type_XM_exp_hits_df_Pro_RE22_IAP <- left_join(LFC_comb_domain_type_XM_exp_hits_df_Pro_RE22, LFC_comb_domain_type_XM[,c("Domain_Name","transcript_id")]) %>% filter(!is.na(Domain_Name))
LFC_comb_domain_type_XM_exp_hits_df_He_IAP <- left_join(LFC_comb_domain_type_XM_exp_hits_df_He, LFC_comb_domain_type_XM[,c("Domain_Name","transcript_id")]) %>% filter(!is.na(Domain_Name))
#LFC_comb_domain_type_XM_exp_hits_df_ROD_IAP <- left_join(LFC_comb_domain_type_XM_exp_hits_df_ROD, LFC_comb_domain_type_XM[,c("Domain_Name","transcript_id")]) %>% filter(!is.na(Domain_Name))
LFC_comb_domain_type_XM_exp_hits_df_Rubio_IAP <- left_join(LFC_comb_domain_type_XM_exp_hits_df_Rubio, LFC_comb_domain_type_XM[,c("Domain_Name","transcript_id")]) %>% filter(!is.na(Domain_Name))
LFC_comb_domain_type_XM_exp_hits_df_Zhang_IAP <- left_join(LFC_comb_domain_type_XM_exp_hits_df_Zhang, LFC_comb_domain_type_XM[,c("Domain_Name","transcript_id")]) %>% filter(!is.na(Domain_Name))

# Find list of modules to export
LFC_comb_domain_type_XM_exp_hits_df_deLorg_Res_IAP %>% distinct(exp, mod_names) # MEturquoise deLorg_Res
LFC_comb_domain_type_XM_exp_hits_df_deLorg_Sus_IAP %>% distinct(exp, mod_names) # MEturquoise deLorg_Sus
LFC_comb_domain_type_XM_exp_hits_df_Dermo_IAP %>% distinct(exp, mod_names) # MEturquoise Dermo_Tol
LFC_comb_domain_type_XM_exp_hits_df_Pro_RE22_IAP %>% distinct(exp, mod_names) # MEroyalblue, MEturquoise, MEsteelblue,  MEdarkslateblue
# are the modules the same between experiments, yes so only need to export these four above
#mod_names       exp            
#<chr>           <chr>          
#1 MEroyalblue     Pro_RE22_Pro_S4
#2 MEturquoise     Pro_RE22_Pro_S4
#3 MEsteelblue     Pro_RE22_Pro_S4
#4 MEturquoise     Pro_RE22_Pro_RI
#5 MEdarkslateblue Pro_RE22_Pro_RI
#6 MEsteelblue     Pro_RE22_Pro_RI

LFC_comb_domain_type_XM_exp_hits_df_He_IAP %>% distinct(exp, mod_names)
#mod_names exp  
#<chr>     <chr>
#  1 MEpurple  He   
#2 MEyellow  He 

LFC_comb_domain_type_XM_exp_hits_df_Rubio_IAP %>% distinct(exp, mod_names) # only export 4 because they are the same
#mod_names   exp     
#<chr>       <chr>   
#1 MEmagenta   Rubio_NV
#2 MEturquoise Rubio_NV
#3 MEblue      Rubio_NV
#4 MEbrown     Rubio_NV

LFC_comb_domain_type_XM_exp_hits_df_Zhang_IAP %>% distinct(exp, mod_names) # MEblack  

#### EXPORT WGCNA MODULES TO CYTOSCAPE ####

# Which modules? 
# deLorg_Res: MEturquoise
# deLorg_Sus: MEturquoise
# Dermo: MEturquoise
# Pro_RE22_Pro_RI:MEturquoise    , MEdarkslateblue, MEsteelblue    
# Pro_RE22_Pro_S4: MEroyalblue, MEturquoise, MEsteelblue
# He:  MEpurple, MEyellow
# Rubio: MEmagenta, MEturquoise, MEblue, MEbrown    
# Zhang: MEblack

# Code is in the TOMsim_cluster.R
# Export full Matrices as tables for each experiment so that I can export the matrices in bluewaves 
# Repeating Dermo_tol_modules. Then move this data to bluewaves
write.table(Dermo_Tolerant_dds_vst_matrix, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Tolerant_dds_vst_matrix.table")
write.table(Pro_RE22_dds_rlog_matrix_Pro, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_dds_rlog_matrix_Pro.table")
write.table(Zhang_dds_rlog_matrix, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Zhang_dds_rlog_matrix.table")
write.table(Rubio_dds_rlog_matrix, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Rubio_dds_rlog_matrix.table")
write.table(deLorgeril_Resistant_dds_vst_matrix, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/deLorgeril_Resistant_dds_vst_matrix.table")
write.table(deLorgeril_Susceptible_dds_vst_matrix,file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/deLorgeril_Susceptible_dds_vst_matrix.table")
write.table(He_dds_vst_matrix,file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/He_dds_vst_matrix.table")

# Export moduleColors file for each 
save(Dermo_Tol_full_moduleColors, file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Tol_full_moduleColors.RData")
save(Pro_RE22_Pro_full_moduleColors, file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_Pro_full_moduleColors.RData")
save(deLorg_Res_full_moduleColors, file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/deLorg_Res_full_moduleColors.RData")
save(deLorg_Sus_full_moduleColors, file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/deLorg_Sus_full_moduleColors.RData")
save(He_full_moduleColors, file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/He_full_moduleColors.RData")
save(Rubio_full_moduleColors, file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Rubio_full_moduleColors.RData")
save(Zhang_full_moduleColors, file= "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Zhang_full_moduleColors.RData") 

#### ANALYZE WGCNA NETWORK DIRECT IAP INTERACTIONS ####

# First export lists of significant IAP transcripts in each experiment to search for in the edges files
LFC_comb_domain_type_XM_exp_hits_df_deLorg_Res_IAP_unique <- LFC_comb_domain_type_XM_exp_hits_df_deLorg_Res_IAP %>% ungroup() %>% distinct(transcript_id)
write.table(LFC_comb_domain_type_XM_exp_hits_df_deLorg_Res_IAP_unique, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/LFC_comb_domain_type_XM_exp_hits_df_deLorg_Res_IAP_unique.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
LFC_comb_domain_type_XM_exp_hits_df_deLorg_Sus_IAP_unique <- LFC_comb_domain_type_XM_exp_hits_df_deLorg_Sus_IAP %>% ungroup() %>% distinct(transcript_id) 
write.table(LFC_comb_domain_type_XM_exp_hits_df_deLorg_Sus_IAP_unique, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/LFC_comb_domain_type_XM_exp_hits_df_deLorg_Sus_IAP_unique.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

LFC_comb_domain_type_XM_exp_hits_df_He_IAP_unique <- LFC_comb_domain_type_XM_exp_hits_df_He_IAP %>% ungroup() %>% distinct(transcript_id)
write.table(LFC_comb_domain_type_XM_exp_hits_df_He_IAP_unique, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/LFC_comb_domain_type_XM_exp_hits_df_He_IAP_unique.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

LFC_comb_domain_type_XM_exp_hits_df_Rubio_IAP_unique <- LFC_comb_domain_type_XM_exp_hits_df_Rubio_IAP %>% ungroup() %>% distinct(transcript_id)
write.table(LFC_comb_domain_type_XM_exp_hits_df_Rubio_IAP_unique, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/LFC_comb_domain_type_XM_exp_hits_df_Rubio_IAP_unique.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

LFC_comb_domain_type_XM_exp_hits_df_Zhang_IAP_unique <- LFC_comb_domain_type_XM_exp_hits_df_Zhang_IAP %>% ungroup() %>% distinct(transcript_id)
write.table(LFC_comb_domain_type_XM_exp_hits_df_Zhang_IAP_unique, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/LFC_comb_domain_type_XM_exp_hits_df_Zhang_IAP_unique.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# also need to add in the rnaID so I can match it with the rnaID that is in the edge files 
LFC_comb_domain_type_XM_exp_hits_df_Dermo_IAP_unique <- LFC_comb_domain_type_XM_exp_hits_df_Dermo_IAP %>% ungroup() %>% distinct(transcript_id) %>% 
  left_join(.,IAP_domain_structure_XM_CV_XM[,c("transcript_id","ID")]) %>% distinct(ID) 
write.table(LFC_comb_domain_type_XM_exp_hits_df_Dermo_IAP_unique, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/LFC_comb_domain_type_XM_exp_hits_df_Dermo_IAP_unique.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

LFC_comb_domain_type_XM_exp_hits_df_Pro_RE22_IAP_unique <- LFC_comb_domain_type_XM_exp_hits_df_Pro_RE22_IAP %>% ungroup() %>% distinct(transcript_id) %>% 
  left_join(.,IAP_domain_structure_XM_CV_XM[,c("transcript_id","ID")]) %>% distinct(ID) 
write.table(LFC_comb_domain_type_XM_exp_hits_df_Pro_RE22_IAP_unique, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/LFC_comb_domain_type_XM_exp_hits_df_Pro_RE22_IAP_unique.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Search for these transcript interaction partners in the command line because the files are too large to load into R workspace (see September 4th code in notes)
# Dermo: MEturquoise
#$ grep -f LFC_comb_domain_type_XM_exp_hits_df_Dermo_IAP_unique.txt CytoscapeInput-edges-Dermo_Tol_fullturquoise.txt > Dermo_Tol_fullturquoise_IAP_hits.txt
# this transcript was not found in the edges file for some reason... not sure why. May need to rerun this to assess. Going to continue by checking the other sheets

# deLorg Res
#$ grep -f LFC_comb_domain_type_XM_exp_hits_df_deLorg_Res_IAP_unique.txt CytoscapeInput-edges-deLorg_Res_fullturquoise.txt > deLorg_Res_fullturquoise_IAP_hits.txt

# de Lorg sus
#grep -f LFC_comb_domain_type_XM_exp_hits_df_deLorg_Sus_IAP_unique.txt CytoscapeInput-edges-deLorg_Sus_fullturquoise.txt > deLorg_Sus_fullturquoise_IAP_hits.txt

# He:  MEpurple, MEyellow
#grep -f LFC_comb_domain_type_XM_exp_hits_df_He_IAP_unique.txt CytoscapeInput-edges-He_fullpurple.txt > He_fullpurple_IAP_hits.txt
#grep -f LFC_comb_domain_type_XM_exp_hits_df_He_IAP_unique.txt CytoscapeInput-edges-He_fullyellow.txt > He_fullyellow_IAP_hits.txt

# Zhang: MEblack
#grep -f LFC_comb_domain_type_XM_exp_hits_df_Zhang_IAP_unique.txt CytoscapeInput-edges-Zhang_fullblack.txt > Zhang_fullblack_IAP_hits.txt

# Rubio: MEmagenta, MEturquoise, MEblue, MEbrown
#grep -f LFC_comb_domain_type_XM_exp_hits_df_Rubio_IAP_unique.txt CytoscapeInput-edges-Rubio_fullblue.txt > Rubio_fullblue_IAP_hits.txt
#grep -f LFC_comb_domain_type_XM_exp_hits_df_Rubio_IAP_unique.txt CytoscapeInput-edges-Rubio_fullbrown.txt > Rubio_fullbrown_IAP_hits.txt
#grep -f LFC_comb_domain_type_XM_exp_hits_df_Rubio_IAP_unique.txt CytoscapeInput-edges-Rubio_fullmagenta.txt > Rubio_fullmagenta_IAP_hits.txt
#grep -f LFC_comb_domain_type_XM_exp_hits_df_Rubio_IAP_unique.txt CytoscapeInput-edges-Rubio_fullturquoise.txt > Rubio_fullturquoise_IAP_hits.txt

## Pro_RE22_Pro_RI:MEturquoise    , MEdarkslateblue, MEsteelblue , MEroyalblue
#grep -f LFC_comb_domain_type_XM_exp_hits_df_Pro_RE22_IAP_unique.txt CytoscapeInput-edges-Pro_RE22_Pro_fulldarkslateblue.txt > Pro_RE22_Pro_fulldarkslateblue_IAP_hits.txt
#grep -f LFC_comb_domain_type_XM_exp_hits_df_Pro_RE22_IAP_unique.txt CytoscapeInput-edges-Pro_RE22_Pro_fullroyalblue.txt > Pro_RE22_Pro_fullroyalblue_IAP_hits.txt
#grep -f LFC_comb_domain_type_XM_exp_hits_df_Pro_RE22_IAP_unique.txt CytoscapeInput-edges-Pro_RE22_Pro_fullsteelblue.txt > Pro_RE22_Pro_fullsteelblue_IAP_hits.txt
#grep -f LFC_comb_domain_type_XM_exp_hits_df_Pro_RE22_IAP_unique.txt CytoscapeInput-edges-Pro_RE22_Pro_fullturquoise.txt > Pro_RE22_Pro_fullturquoise_IAP_hits.txt

# Check that all the edge files were searched
# $ ls *IAP_hits.txt | wc -l # 14 this is the same as the the number of modules to go through

## Import IAP interaction hits into R
deLorg_Res_fullturquoise_IAP_hits <- read.table(file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/deLorg_Res_fullturquoise_IAP_hits.txt", col.names = c("fromNode",	"toNode",	"weight",	"direction",	"fromAltName",	"toAltName"))
deLorg_Sus_fullturquoise_IAP_hits <- read.table(file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/deLorg_Sus_fullturquoise_IAP_hits.txt", col.names = c("fromNode",	"toNode",	"weight",	"direction",	"fromAltName",	"toAltName"))
He_fullpurple_IAP_hits <- read.table(file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/He_fullpurple_IAP_hits.txt", col.names = c("fromNode",	"toNode",	"weight",	"direction",	"fromAltName",	"toAltName"))
He_fullyellow_IAP_hits <- read.table(file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/He_fullyellow_IAP_hits.txt", col.names = c("fromNode",	"toNode",	"weight",	"direction",	"fromAltName",	"toAltName"))
Zhang_fullblack_IAP_hits <- read.table(file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Zhang_fullblack_IAP_hits.txt", col.names = c("fromNode",	"toNode",	"weight",	"direction",	"fromAltName",	"toAltName"))
Rubio_fullblue_IAP_hits <- read.table(file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Rubio_fullblue_IAP_hits.txt", col.names = c("fromNode",	"toNode",	"weight",	"direction",	"fromAltName",	"toAltName"))
Rubio_fullmagenta_IAP_hits <- read.table(file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Rubio_fullbrown_IAP_hits.txt", col.names = c("fromNode",	"toNode",	"weight",	"direction",	"fromAltName",	"toAltName"))
Rubio_fullmagenta_IAP_hits <- read.table(file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Rubio_fullmagenta_IAP_hits.txt", col.names = c("fromNode",	"toNode",	"weight",	"direction",	"fromAltName",	"toAltName"))
Rubio_fullturquoise_IAP_hits <- read.table(file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Rubio_fullturquoise_IAP_hits.txt", col.names = c("fromNode",	"toNode",	"weight",	"direction",	"fromAltName",	"toAltName"))
Pro_RE22_Pro_fulldarkslateblue_IAP_hits <- read.table(file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_Pro_fulldarkslateblue_IAP_hits.txt", col.names = c("fromNode",	"toNode",	"weight",	"direction",	"fromAltName",	"toAltName"))
Pro_RE22_Pro_fullroyalblue_IAP_hits <- read.table(file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_Pro_fullroyalblue_IAP_hits.txt", col.names = c("fromNode",	"toNode",	"weight",	"direction",	"fromAltName",	"toAltName"))
Pro_RE22_Pro_fullsteelblue_IAP_hits <- read.table(file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_Pro_fullsteelblue_IAP_hits.txt", col.names = c("fromNode",	"toNode",	"weight",	"direction",	"fromAltName",	"toAltName"))
Pro_RE22_Pro_fullturquoise_IAP_hits <- read.table(file = "/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_Pro_fullturquoise_IAP_hits.txt", col.names = c("fromNode",	"toNode",	"weight",	"direction",	"fromAltName",	"toAltName"))

## Annotate fromNode and toNodes by joining with apoptosis data frame and annotate BIR domain type (we only care about IAP interactions with other apoptosis pathway proteins, NAs will be non apoptosis transcripts)
deLorg_Res_fullturquoise_IAP_hits_annot <- deLorg_Res_fullturquoise_IAP_hits %>% dplyr::select(fromNode, toNode, weight, direction) %>% dplyr::rename(transcript_id = fromNode) %>% left_join(., C_gig_rtracklayer_apop_product_final[,c("gene","product", "transcript_id")]) %>%
  left_join(., LFC_comb_domain_type_XM_exp_hits_df_deLorg_Res_IAP[,c("transcript_id","Domain_Name")]) %>% dplyr::rename (fromNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (fromNodetranscript = transcript_id) %>% dplyr::rename(fromNodeproduct = product) %>% dplyr::rename(fromNodegene = gene) %>% dplyr::rename(transcript_id = toNode) %>% left_join(., C_gig_rtracklayer_apop_product_final[,c("gene","product", "transcript_id")]) %>%
  left_join(., LFC_comb_domain_type_XM_exp_hits_df_deLorg_Res_IAP[,c("transcript_id","Domain_Name")]) %>% dplyr::rename (toNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (toNodetranscript = transcript_id) %>% dplyr::rename(toNodeproduct = product) %>% dplyr::rename(toNodegene = gene) %>% filter(!is.na(fromNodeproduct)) %>% filter(!is.na(toNodeproduct)) %>% distinct()

deLorg_Sus_fullturquoise_IAP_hits_annot <- deLorg_Sus_fullturquoise_IAP_hits  %>% dplyr::select(fromNode, toNode, weight, direction) %>% dplyr::rename(transcript_id = fromNode) %>% left_join(., C_gig_rtracklayer_apop_product_final[,c("gene","product", "transcript_id")]) %>%
  left_join(., LFC_comb_domain_type_XM_exp_hits_df_deLorg_Sus_IAP[,c("transcript_id","Domain_Name")]) %>% dplyr::rename (fromNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (fromNodetranscript = transcript_id) %>% dplyr::rename(fromNodeproduct = product) %>% dplyr::rename(fromNodegene = gene) %>% dplyr::rename(transcript_id = toNode) %>% left_join(., C_gig_rtracklayer_apop_product_final[,c("gene","product", "transcript_id")]) %>%
  left_join(., LFC_comb_domain_type_XM_exp_hits_df_deLorg_Sus_IAP[,c("transcript_id","Domain_Name")]) %>% dplyr::rename (toNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (toNodetranscript = transcript_id) %>% dplyr::rename(toNodeproduct = product) %>% dplyr::rename(toNodegene = gene) %>% filter(!is.na(fromNodeproduct)) %>% filter(!is.na(toNodeproduct)) %>% distinct()

He_fullpurple_IAP_hits_annot <- He_fullpurple_IAP_hits  %>% dplyr::select(fromNode, toNode, weight, direction) %>% dplyr::rename(transcript_id = fromNode) %>% left_join(., C_gig_rtracklayer_apop_product_final[,c("gene","product", "transcript_id")]) %>%
  left_join(., LFC_comb_domain_type_XM_exp_hits_df_He_IAP[,c("transcript_id","Domain_Name")]) %>% dplyr::rename (fromNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (fromNodetranscript = transcript_id) %>% dplyr::rename(fromNodeproduct = product) %>% dplyr::rename(fromNodegene = gene) %>% dplyr::rename(transcript_id = toNode) %>% left_join(., C_gig_rtracklayer_apop_product_final[,c("gene","product", "transcript_id")]) %>%
  left_join(., LFC_comb_domain_type_XM_exp_hits_df_He_IAP[,c("transcript_id","Domain_Name")]) %>% dplyr::rename (toNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (toNodetranscript = transcript_id) %>% dplyr::rename(toNodeproduct = product) %>% dplyr::rename(toNodegene = gene) %>% filter(!is.na(fromNodeproduct)) %>% filter(!is.na(toNodeproduct)) %>% distinct()

He_fullyellow_IAP_hits_annot <- He_fullyellow_IAP_hits  %>% dplyr::select(fromNode, toNode, weight, direction) %>% dplyr::rename(transcript_id = fromNode) %>% left_join(., C_gig_rtracklayer_apop_product_final[,c("gene","product", "transcript_id")]) %>%
  left_join(., LFC_comb_domain_type_XM_exp_hits_df_He_IAP[,c("transcript_id","Domain_Name")]) %>% dplyr::rename (fromNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (fromNodetranscript = transcript_id) %>% dplyr::rename(fromNodeproduct = product) %>% dplyr::rename(fromNodegene = gene) %>% dplyr::rename(transcript_id = toNode) %>% left_join(., C_gig_rtracklayer_apop_product_final[,c("gene","product", "transcript_id")]) %>%
  left_join(., LFC_comb_domain_type_XM_exp_hits_df_He_IAP[,c("transcript_id","Domain_Name")]) %>% dplyr::rename (toNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (toNodetranscript = transcript_id) %>% dplyr::rename(toNodeproduct = product) %>% dplyr::rename(toNodegene = gene) %>% filter(!is.na(fromNodeproduct)) %>% filter(!is.na(toNodeproduct)) %>% distinct()

Zhang_fullblack_IAP_hits_annot <- Zhang_fullblack_IAP_hits  %>% dplyr::select(fromNode, toNode, weight, direction) %>% dplyr::rename(transcript_id = fromNode) %>% left_join(., C_gig_rtracklayer_apop_product_final[,c("gene","product", "transcript_id")]) %>%
  left_join(., LFC_comb_domain_type_XM_exp_hits_df_Zhang_IAP[,c("transcript_id","Domain_Name")]) %>% dplyr::rename (fromNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (fromNodetranscript = transcript_id) %>% dplyr::rename(fromNodeproduct = product) %>% dplyr::rename(fromNodegene = gene) %>% dplyr::rename(transcript_id = toNode) %>% left_join(., C_gig_rtracklayer_apop_product_final[,c("gene","product", "transcript_id")]) %>%
  left_join(., LFC_comb_domain_type_XM_exp_hits_df_Zhang_IAP[,c("transcript_id","Domain_Name")]) %>% dplyr::rename (toNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (toNodetranscript = transcript_id) %>% dplyr::rename(toNodeproduct = product) %>% dplyr::rename(toNodegene = gene) %>% filter(!is.na(fromNodeproduct)) %>% filter(!is.na(toNodeproduct)) %>% distinct()

Rubio_fullblue_IAP_hits_annot <- Rubio_fullblue_IAP_hits  %>% dplyr::select(fromNode, toNode, weight, direction) %>% dplyr::rename(transcript_id = fromNode) %>% left_join(., C_gig_rtracklayer_apop_product_final[,c("gene","product", "transcript_id")]) %>%
  left_join(., LFC_comb_domain_type_XM_exp_hits_df_Rubio_IAP[,c("transcript_id","Domain_Name")]) %>% dplyr::rename (fromNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (fromNodetranscript = transcript_id) %>% dplyr::rename(fromNodeproduct = product) %>% dplyr::rename(fromNodegene = gene) %>% dplyr::rename(transcript_id = toNode) %>% left_join(., C_gig_rtracklayer_apop_product_final[,c("gene","product", "transcript_id")]) %>%
  left_join(., LFC_comb_domain_type_XM_exp_hits_df_Rubio_IAP[,c("transcript_id","Domain_Name")]) %>% dplyr::rename (toNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (toNodetranscript = transcript_id) %>% dplyr::rename(toNodeproduct = product) %>% dplyr::rename(toNodegene = gene) %>% filter(!is.na(fromNodeproduct)) %>% filter(!is.na(toNodeproduct)) %>% distinct()

Rubio_fullmagenta_IAP_hits_annot <- Rubio_fullmagenta_IAP_hits  %>% dplyr::select(fromNode, toNode, weight, direction) %>% dplyr::rename(transcript_id = fromNode) %>% left_join(., C_gig_rtracklayer_apop_product_final[,c("gene","product", "transcript_id")]) %>%
  left_join(., LFC_comb_domain_type_XM_exp_hits_df_Rubio_IAP[,c("transcript_id","Domain_Name")]) %>% dplyr::rename (fromNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (fromNodetranscript = transcript_id) %>% dplyr::rename(fromNodeproduct = product) %>% dplyr::rename(fromNodegene = gene) %>% dplyr::rename(transcript_id = toNode) %>% left_join(., C_gig_rtracklayer_apop_product_final[,c("gene","product", "transcript_id")]) %>%
  left_join(., LFC_comb_domain_type_XM_exp_hits_df_Rubio_IAP[,c("transcript_id","Domain_Name")]) %>% dplyr::rename (toNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (toNodetranscript = transcript_id) %>% dplyr::rename(toNodeproduct = product) %>% dplyr::rename(toNodegene = gene) %>% filter(!is.na(fromNodeproduct)) %>% filter(!is.na(toNodeproduct)) %>% distinct()

Rubio_fullmagenta_IAP_hits_annot <- Rubio_fullmagenta_IAP_hits %>% dplyr::select(fromNode, toNode, weight, direction) %>% dplyr::rename(transcript_id = fromNode) %>% left_join(., C_gig_rtracklayer_apop_product_final[,c("gene","product", "transcript_id")]) %>%
  left_join(., LFC_comb_domain_type_XM_exp_hits_df_Rubio_IAP[,c("transcript_id","Domain_Name")]) %>% dplyr::rename (fromNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (fromNodetranscript = transcript_id) %>% dplyr::rename(fromNodeproduct = product) %>% dplyr::rename(fromNodegene = gene) %>% dplyr::rename(transcript_id = toNode) %>% left_join(., C_gig_rtracklayer_apop_product_final[,c("gene","product", "transcript_id")]) %>%
  left_join(., LFC_comb_domain_type_XM_exp_hits_df_Rubio_IAP[,c("transcript_id","Domain_Name")]) %>% dplyr::rename (toNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (toNodetranscript = transcript_id) %>% dplyr::rename(toNodeproduct = product) %>% dplyr::rename(toNodegene = gene) %>% filter(!is.na(fromNodeproduct)) %>% filter(!is.na(toNodeproduct)) %>% distinct()

Rubio_fullturquoise_IAP_hits_annot <- Rubio_fullturquoise_IAP_hits  %>% dplyr::select(fromNode, toNode, weight, direction) %>% dplyr::rename(transcript_id = fromNode) %>% left_join(., C_gig_rtracklayer_apop_product_final[,c("gene","product", "transcript_id")]) %>%
  left_join(., LFC_comb_domain_type_XM_exp_hits_df_Rubio_IAP[,c("transcript_id","Domain_Name")]) %>% dplyr::rename (fromNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (fromNodetranscript = transcript_id) %>% dplyr::rename(fromNodeproduct = product) %>% dplyr::rename(fromNodegene = gene) %>% dplyr::rename(transcript_id = toNode) %>% left_join(., C_gig_rtracklayer_apop_product_final[,c("gene","product", "transcript_id")]) %>%
  left_join(., LFC_comb_domain_type_XM_exp_hits_df_Rubio_IAP[,c("transcript_id","Domain_Name")]) %>% dplyr::rename (toNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (toNodetranscript = transcript_id) %>% dplyr::rename(toNodeproduct = product) %>% dplyr::rename(toNodegene = gene) %>% filter(!is.na(fromNodeproduct)) %>% filter(!is.na(toNodeproduct)) %>% distinct()

# annot C. virginica modules 
Pro_RE22_Pro_fulldarkslateblue_IAP_hits_annot <- Pro_RE22_Pro_fulldarkslateblue_IAP_hits %>% dplyr::select(fromNode, toNode, weight, direction) %>% dplyr::rename(ID = fromNode) %>% left_join(., C_vir_rtracklayer_apop_product_final[,c("gene","product", "ID")]) %>%
  left_join(., IAP_domain_structure_XM_CV_XM[,c("ID","Domain_Name")]) %>% dplyr::rename (fromNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (fromNodetranscript = ID) %>% dplyr::rename(fromNodeproduct = product) %>% dplyr::rename(fromNodegene = gene) %>% dplyr::rename(ID = toNode) %>% left_join(., C_vir_rtracklayer_apop_product_final[,c("gene","product", "ID")]) %>%
  left_join(., IAP_domain_structure_XM_CV_XM[,c("ID","Domain_Name")]) %>% dplyr::rename (toNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (toNodetranscript = ID) %>% dplyr::rename(toNodeproduct = product) %>% dplyr::rename(toNodegene = gene) %>% filter(!is.na(fromNodeproduct)) %>% filter(!is.na(toNodeproduct)) %>% distinct()

Pro_RE22_Pro_fullroyalblue_IAP_hits_annot <- Pro_RE22_Pro_fullroyalblue_IAP_hits %>% dplyr::select(fromNode, toNode, weight, direction) %>% dplyr::rename(ID = fromNode) %>% left_join(., C_vir_rtracklayer_apop_product_final[,c("gene","product", "ID")]) %>%
  left_join(., IAP_domain_structure_XM_CV_XM[,c("ID","Domain_Name")]) %>% dplyr::rename (fromNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (fromNodetranscript = ID) %>% dplyr::rename(fromNodeproduct = product) %>% dplyr::rename(fromNodegene = gene) %>% dplyr::rename(ID = toNode) %>% left_join(., C_vir_rtracklayer_apop_product_final[,c("gene","product", "ID")]) %>%
  left_join(., IAP_domain_structure_XM_CV_XM[,c("ID","Domain_Name")]) %>% dplyr::rename (toNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (toNodetranscript = ID) %>% dplyr::rename(toNodeproduct = product) %>% dplyr::rename(toNodegene = gene) %>% filter(!is.na(fromNodeproduct)) %>% filter(!is.na(toNodeproduct)) %>% distinct()

Pro_RE22_Pro_fullsteelblue_IAP_hits_annot <- Pro_RE22_Pro_fullsteelblue_IAP_hits %>% dplyr::select(fromNode, toNode, weight, direction) %>% dplyr::rename(ID = fromNode) %>% left_join(., C_vir_rtracklayer_apop_product_final[,c("gene","product", "ID")]) %>%
  left_join(., IAP_domain_structure_XM_CV_XM[,c("ID","Domain_Name")]) %>% dplyr::rename (fromNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (fromNodetranscript = ID) %>% dplyr::rename(fromNodeproduct = product) %>% dplyr::rename(fromNodegene = gene) %>% dplyr::rename(ID = toNode) %>% left_join(., C_vir_rtracklayer_apop_product_final[,c("gene","product", "ID")]) %>%
  left_join(., IAP_domain_structure_XM_CV_XM[,c("ID","Domain_Name")]) %>% dplyr::rename (toNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (toNodetranscript = ID) %>% dplyr::rename(toNodeproduct = product) %>% dplyr::rename(toNodegene = gene) %>% filter(!is.na(fromNodeproduct)) %>% filter(!is.na(toNodeproduct)) %>% distinct()

Pro_RE22_Pro_fullturquoise_IAP_hits_annot <- Pro_RE22_Pro_fullturquoise_IAP_hits %>% dplyr::select(fromNode, toNode, weight, direction) %>% dplyr::rename(ID = fromNode) %>% left_join(., C_vir_rtracklayer_apop_product_final[,c("gene","product", "ID")]) %>%
  left_join(., IAP_domain_structure_XM_CV_XM[,c("ID","Domain_Name")]) %>% dplyr::rename (fromNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (fromNodetranscript = ID) %>% dplyr::rename(fromNodeproduct = product) %>% dplyr::rename(fromNodegene = gene) %>% dplyr::rename(ID = toNode) %>% left_join(., C_vir_rtracklayer_apop_product_final[,c("gene","product", "ID")]) %>%
  left_join(., IAP_domain_structure_XM_CV_XM[,c("ID","Domain_Name")]) %>% dplyr::rename (toNodeDomain_Name = Domain_Name) %>% 
  dplyr::rename (toNodetranscript = ID) %>% dplyr::rename(toNodeproduct = product) %>% dplyr::rename(toNodegene = gene) %>% filter(!is.na(fromNodeproduct)) %>% filter(!is.na(toNodeproduct)) %>% distinct()

### GOALS OF THIS SECTION 
#1. Identify any specific important interactions of IAPs with molecules of interes 
    #a. caspases
    #b. TLRs, TRAFs, RIPK
    #c. Bcl-2 family members 
#2. Identify patterns between which domain types are interacting with which group of molecules 

#### COMPARE CONSENSUS AND FULL IAP AND GIMAP ACROSS ALL DATASETS, NOT JUST THOSE IN SIGNIFICANT MODULES ####

# Number of GIMAPs and IAPs in consensus set 
C_vir_common_vst_transcripts_df <- as.data.frame(C_vir_common_vst_transcripts)
colnames(C_vir_common_vst_transcripts_df)[1] <- "ID"
C_vir_common_vst_transcripts_df_annot <- left_join(C_vir_common_vst_transcripts_df, C_vir_rtracklayer_apop_product_final[,c("ID","product","transcript_id","gene")])
C_vir_common_vst_transcripts_df_annot <- C_vir_common_vst_transcripts_df_annot %>% filter(!is.na(transcript_id ))
nrow(C_vir_common_vst_transcripts_df_annot ) # 598
nrow(C_vir_rtracklayer_apop_product_final) # 1245 ....I've lost more than half of my genes by doing it this way

# Number of IAPs in consensus set
C_vir_common_vst_transcripts_df_annot_IAP <- C_vir_common_vst_transcripts_df_annot[grepl("IAP", C_vir_common_vst_transcripts_df_annot$product, ignore.case = TRUE),]
nrow(C_vir_common_vst_transcripts_df_annot_IAP) # 48 unique IAPs remained in my consensus set 
# IAP genes in consensus 
C_vir_common_vst_transcripts_df_annot_IAP %>% group_by(gene) %>% dplyr::summarise(count=n()) %>% nrow() # 38 
C_vir_IAP_consensus_genes <- C_vir_common_vst_transcripts_df_annot_IAP %>% group_by(gene) %>% dplyr::summarise(product= first(product), n=n())
C_vir_IAP_consensus_genes$product <- str_remove(C_vir_IAP_consensus_genes$product, ", transcript variant X1")
C_vir_IAP_consensus_genes$product <- str_remove(C_vir_IAP_consensus_genes$product, ", transcript variant X2")
C_vir_IAP_consensus_genes$product <- str_remove(C_vir_IAP_consensus_genes$product, ", transcript variant X7")

View(C_vir_IAP_consensus_genes %>% group_by(product) %>% dplyr::summarise(n=n()))

# Number of GIMAPs in consensus set
C_vir_common_vst_transcripts_df_annot_GIMAP <- C_vir_common_vst_transcripts_df_annot[grepl("IMAP", C_vir_common_vst_transcripts_df_annot$product, ignore.case = TRUE),]
nrow(C_vir_common_vst_transcripts_df_annot_GIMAP) # 25 unique GIMAP transcripts remained in my consensus set 
# GIMAP genes in consensus 
C_vir_common_vst_transcripts_df_annot_GIMAP %>% group_by(gene) %>% dplyr::summarise(count=n()) %>% nrow() # 19
C_vir_GIMAP_consensus_genes <- C_vir_common_vst_transcripts_df_annot_GIMAP %>% group_by(gene) %>% dplyr::summarise(product= first(product), n=n())
C_vir_GIMAP_consensus_genes$product <- str_remove(C_vir_GIMAP_consensus_genes$product, ", transcript variant X1")
C_vir_GIMAP_consensus_genes$product <- str_remove(C_vir_GIMAP_consensus_genes$product, ", transcript variant X2")
View(C_vir_GIMAP_consensus_genes %>% group_by(product) %>% dplyr::summarise(n=n()))

# IAPs in each set 
C_vir_rtracklayer_apop_product_final_IAP_list <- C_vir_rtracklayer_apop_product_final[grepl("IAP",C_vir_rtracklayer_apop_product_final$product, ignore.case = TRUE),]
C_vir_rtracklayer_apop_product_final_IAP_list <- C_vir_rtracklayer_apop_product_final_IAP_list$ID
length(C_vir_rtracklayer_apop_product_final_IAP_list) # 105

Dermo_Tolerant_dds_vst_matrix_IAP <- colnames(Dermo_Tolerant_dds_vst_matrix)[colnames(Dermo_Tolerant_dds_vst_matrix) %in% C_vir_rtracklayer_apop_product_final_IAP_list]
Dermo_Susceptible_dds_vst_matrix_IAP <- colnames(Dermo_Susceptible_dds_vst_matrix)[colnames(Dermo_Susceptible_dds_vst_matrix) %in% C_vir_rtracklayer_apop_product_final_IAP_list]
Probiotic_dds_rlog_matrix_IAP <- colnames(Probiotic_dds_rlog_matrix)[colnames(Probiotic_dds_rlog_matrix) %in% C_vir_rtracklayer_apop_product_final_IAP_list]
ROD_Susceptible_dds_rlog_matrix_IAP <- colnames(ROD_Susceptible_dds_rlog_matrix)[colnames(ROD_Susceptible_dds_rlog_matrix) %in% C_vir_rtracklayer_apop_product_final_IAP_list]
ROD_Resistant_dds_rlog_matrix_IAP <- colnames(ROD_Resistant_dds_rlog_matrix)[colnames(ROD_Resistant_dds_rlog_matrix) %in% C_vir_rtracklayer_apop_product_final_IAP_list]
Pro_RE22_dds_rlog_matrix_IAP <- colnames(Pro_RE22_dds_rlog_matrix)[colnames(Pro_RE22_dds_rlog_matrix) %in% C_vir_rtracklayer_apop_product_final_IAP_list]

Dermo_Tolerant_dds_vst_matrix_IAP <- C_vir_rtracklayer_apop_product_final[C_vir_rtracklayer_apop_product_final$ID %in% Dermo_Tolerant_dds_vst_matrix_IAP,]
Dermo_Susceptible_dds_vst_matrix_IAP <- C_vir_rtracklayer_apop_product_final[C_vir_rtracklayer_apop_product_final$ID %in% Dermo_Susceptible_dds_vst_matrix_IAP,]
Probiotic_dds_rlog_matrix_IAP <- C_vir_rtracklayer_apop_product_final[C_vir_rtracklayer_apop_product_final$ID %in% Probiotic_dds_rlog_matrix_IAP,]
ROD_Susceptible_dds_rlog_matrix_IAP <- C_vir_rtracklayer_apop_product_final[C_vir_rtracklayer_apop_product_final$ID %in% ROD_Susceptible_dds_rlog_matrix_IAP,]
ROD_Resistant_dds_rlog_matrix_IAP <- C_vir_rtracklayer_apop_product_final[C_vir_rtracklayer_apop_product_final$ID %in% ROD_Resistant_dds_rlog_matrix_IAP,]
Pro_RE22_dds_rlog_matrix_IAP <- C_vir_rtracklayer_apop_product_final[C_vir_rtracklayer_apop_product_final$ID %in% Pro_RE22_dds_rlog_matrix_IAP,]
nrow(Pro_RE22_dds_rlog_matrix_IAP) # 76 
nrow(Dermo_Tolerant_dds_vst_matrix_IAP) # 77

Dermo_Tolerant_dds_vst_matrix_IAP$exp <- "Dermo_Tolerant"
Dermo_Susceptible_dds_vst_matrix_IAP$exp <- "Dermo_Susceptible"
Probiotic_dds_rlog_matrix_IAP$exp <- "Probiotic"
ROD_Susceptible_dds_rlog_matrix_IAP$exp <- "ROD_Susceptible"
ROD_Resistant_dds_rlog_matrix_IAP$exp <- "ROD_Resistant"
Pro_RE22_dds_rlog_matrix_IAP$exp <- "Pro_RE22"

C_vir_IAP <- rbind(Dermo_Tolerant_dds_vst_matrix_IAP,
                   Dermo_Susceptible_dds_vst_matrix_IAP,
                   Probiotic_dds_rlog_matrix_IAP,
                   ROD_Susceptible_dds_rlog_matrix_IAP,
                   ROD_Resistant_dds_rlog_matrix_IAP,
                   Pro_RE22_dds_rlog_matrix_IAP)

## of IAP genes from genome utilized across all the transcriptomes
C_vir_rtracklayer_apop_product_final_IAP_list[!(C_vir_rtracklayer_apop_product_final_IAP_list$transcript_id %in% C_vir_IAP$transcript_id),]
# 4 transcripts not present in the transcriptomes
C_vir_rtracklayer_apop_product_final_IAP_list[!(C_vir_rtracklayer_apop_product_final_IAP_list$gene %in% C_vir_IAP$gene),]
# 0 genes not present in the transcriptomes 

C_vir_IAP_number <- C_vir_IAP %>% group_by(exp) %>% dplyr::summarise(count = n())
C_vir_IAP_number$species <- "C_vir"
C_vir_IAP_number$type <- "IAP"

# GIMAPs in each set 
# IAPs in each set 
C_vir_rtracklayer_apop_product_final_IMAP_list <- C_vir_rtracklayer_apop_product_final[grepl("IMAP",C_vir_rtracklayer_apop_product_final$product, ignore.case = TRUE),]
C_vir_rtracklayer_apop_product_final_IMAP_list <- C_vir_rtracklayer_apop_product_final_IMAP_list$ID
length(C_vir_rtracklayer_apop_product_final_IMAP_list) # 105

Dermo_Tolerant_dds_vst_matrix_IMAP <- colnames(Dermo_Tolerant_dds_vst_matrix)[colnames(Dermo_Tolerant_dds_vst_matrix) %in% C_vir_rtracklayer_apop_product_final_IMAP_list]
Dermo_Susceptible_dds_vst_matrix_IMAP <- colnames(Dermo_Susceptible_dds_vst_matrix)[colnames(Dermo_Susceptible_dds_vst_matrix) %in% C_vir_rtracklayer_apop_product_final_IMAP_list]
Probiotic_dds_rlog_matrix_IMAP <- colnames(Probiotic_dds_rlog_matrix)[colnames(Probiotic_dds_rlog_matrix) %in% C_vir_rtracklayer_apop_product_final_IMAP_list]
ROD_Susceptible_dds_rlog_matrix_IMAP <- colnames(ROD_Susceptible_dds_rlog_matrix)[colnames(ROD_Susceptible_dds_rlog_matrix) %in% C_vir_rtracklayer_apop_product_final_IMAP_list]
ROD_Resistant_dds_rlog_matrix_IMAP <- colnames(ROD_Resistant_dds_rlog_matrix)[colnames(ROD_Resistant_dds_rlog_matrix) %in% C_vir_rtracklayer_apop_product_final_IMAP_list]
Pro_RE22_dds_rlog_matrix_IMAP <- colnames(Pro_RE22_dds_rlog_matrix)[colnames(Pro_RE22_dds_rlog_matrix) %in% C_vir_rtracklayer_apop_product_final_IMAP_list]

Dermo_Tolerant_dds_vst_matrix_IMAP <- C_vir_rtracklayer_apop_product_final[C_vir_rtracklayer_apop_product_final$ID %in% Dermo_Tolerant_dds_vst_matrix_IMAP,]
Dermo_Susceptible_dds_vst_matrix_IMAP <- C_vir_rtracklayer_apop_product_final[C_vir_rtracklayer_apop_product_final$ID %in% Dermo_Susceptible_dds_vst_matrix_IMAP,]
Probiotic_dds_rlog_matrix_IMAP <- C_vir_rtracklayer_apop_product_final[C_vir_rtracklayer_apop_product_final$ID %in% Probiotic_dds_rlog_matrix_IMAP,]
ROD_Susceptible_dds_rlog_matrix_IMAP <- C_vir_rtracklayer_apop_product_final[C_vir_rtracklayer_apop_product_final$ID %in% ROD_Susceptible_dds_rlog_matrix_IMAP,]
ROD_Resistant_dds_rlog_matrix_IMAP <- C_vir_rtracklayer_apop_product_final[C_vir_rtracklayer_apop_product_final$ID %in% ROD_Resistant_dds_rlog_matrix_IMAP,]
Pro_RE22_dds_rlog_matrix_IMAP <- C_vir_rtracklayer_apop_product_final[C_vir_rtracklayer_apop_product_final$ID %in% Pro_RE22_dds_rlog_matrix_IMAP,]
nrow(Pro_RE22_dds_rlog_matrix_IMAP) # 39
nrow(Dermo_Tolerant_dds_vst_matrix_IMAP ) # 48

Dermo_Tolerant_dds_vst_matrix_IMAP$exp <- "Dermo_Tolerant"
Dermo_Susceptible_dds_vst_matrix_IMAP$exp <- "Dermo_Susceptible"
Probiotic_dds_rlog_matrix_IMAP$exp <- "Probiotic"
ROD_Susceptible_dds_rlog_matrix_IMAP$exp <- "ROD_Susceptible"
ROD_Resistant_dds_rlog_matrix_IMAP$exp <- "ROD_Resistant"
Pro_RE22_dds_rlog_matrix_IMAP$exp <- "Pro_RE22"

C_vir_IMAP <- rbind(Dermo_Tolerant_dds_vst_matrix_IMAP,
                    Dermo_Susceptible_dds_vst_matrix_IMAP,
                    Probiotic_dds_rlog_matrix_IMAP,
                    ROD_Susceptible_dds_rlog_matrix_IMAP,
                    ROD_Resistant_dds_rlog_matrix_IMAP,
                    Pro_RE22_dds_rlog_matrix_IMAP)

## of IAP genes from genome utilized across all the transcriptomes
nrow(C_vir_rtracklayer_apop_product_final_IMAP_list[!(C_vir_rtracklayer_apop_product_final_IMAP_list$transcript_id %in% C_vir_IMAP$transcript_id),])
# 32 transcripts not present in the transcriptomes
nrow(C_vir_rtracklayer_apop_product_final_IMAP_list[!(C_vir_rtracklayer_apop_product_final_IMAP_list$gene %in% C_vir_IMAP$gene),])
# 10 genes not present in the transcriptomes 

C_vir_IMAP_number <- C_vir_IMAP %>% group_by(exp) %>% dplyr::summarise(count = n())
C_vir_IMAP_number$species <- "C_vir"
C_vir_IMAP_number$type <- "GIMAP"


#### COMPARE CGIGAS CONSENSUS AND FULL IAP AND GIMAP ####

# Number of GIMAPs and IAPs in consensus set 
C_gig_common_vst_transcripts_df <- as.data.frame(C_gig_common_vst_transcripts)
colnames(C_gig_common_vst_transcripts_df)[1] <- "transcript_id"
C_gig_common_vst_transcripts_df_annot <- left_join(C_gig_common_vst_transcripts_df, C_gig_rtracklayer_apop_product_final[,c("product","transcript_id","gene")])
C_gig_common_vst_transcripts_df_annot <- C_gig_common_vst_transcripts_df_annot %>% filter(!is.na(product))
nrow(C_gig_common_vst_transcripts_df_annot ) # 546
nrow(C_gig_rtracklayer_apop_product_final) # 833

# Number of IAPs in consensus set
C_gig_common_vst_transcripts_df_annot_IAP <- C_gig_common_vst_transcripts_df_annot[grepl("IAP", C_gig_common_vst_transcripts_df_annot$product, ignore.case = TRUE),]
nrow(C_gig_common_vst_transcripts_df_annot_IAP) # 24 unique IAPs remained in my consensus set 
# IAP genes in consensus 
C_gig_common_vst_transcripts_df_annot_IAP %>% group_by(gene) %>% dplyr::summarise(count=n()) %>% nrow() # 25 
C_gig_IAP_consensus_genes <- C_gig_common_vst_transcripts_df_annot_IAP %>% group_by(gene) %>% dplyr::summarise(product= first(product), n=n())
C_gig_IAP_consensus_genes$product <- str_remove(C_gig_IAP_consensus_genes$product, ", transcript variant X1")
C_gig_IAP_consensus_genes$product <- str_remove(C_gig_IAP_consensus_genes$product, ", transcript variant X2")
C_gig_IAP_consensus_genes$product <- str_remove(C_gig_IAP_consensus_genes$product, ", transcript variant X3")
View(C_gig_IAP_consensus_genes %>% group_by(product) %>% dplyr::summarise(n=n()))

# Number of GIMAPs in consensus set
C_gig_common_vst_transcripts_df_annot_GIMAP <- C_gig_common_vst_transcripts_df_annot[grepl("IMAP", C_gig_common_vst_transcripts_df_annot$product, ignore.case = TRUE),]
nrow(C_gig_common_vst_transcripts_df_annot_GIMAP) # 15  unique GIMAP transcripts remained in my consensus set 
# GIMAP genes in consensus 
C_gig_common_vst_transcripts_df_annot_GIMAP %>% group_by(gene) %>% dplyr::summarise(count=n()) %>% nrow() # 13 
C_gig_GIMAP_consensus_genes <- C_gig_common_vst_transcripts_df_annot_GIMAP %>% group_by(gene) %>% dplyr::summarise(product= first(product), n=n())
C_gig_GIMAP_consensus_genes$product <- str_remove(C_gig_GIMAP_consensus_genes$product, ", transcript variant X2")
View(C_gig_GIMAP_consensus_genes %>% group_by(product) %>% dplyr::summarise(n=n()))

# IAPs in each set 
C_gig_rtracklayer_apop_product_final_IAP_list <- C_gig_rtracklayer_apop_product_final[grepl("IAP",C_gig_rtracklayer_apop_product_final$product, ignore.case = TRUE),]
C_gig_rtracklayer_apop_product_final_IAP_list <- C_gig_rtracklayer_apop_product_final_IAP_list$transcript_id
length(C_gig_rtracklayer_apop_product_final_IAP_list) # 51

Zhang_dds_rlog_matrix_IAP <- colnames(Zhang_dds_rlog_matrix)[colnames(Zhang_dds_rlog_matrix) %in% C_gig_rtracklayer_apop_product_final_IAP_list]
Zhang_dds_rlog_matrix_IAP <- C_gig_rtracklayer_apop_product_final[C_gig_rtracklayer_apop_product_final$transcript_id %in% Zhang_dds_rlog_matrix_IAP,]
Zhang_dds_rlog_matrix_IAP$exp <- "Zhang"
Rubio_dds_rlog_matrix_IAP <- colnames(Rubio_dds_rlog_matrix)[colnames(Zhang_dds_rlog_matrix) %in% C_gig_rtracklayer_apop_product_final_IAP_list]
Rubio_dds_rlog_matrix_IAP <- C_gig_rtracklayer_apop_product_final[C_gig_rtracklayer_apop_product_final$transcript_id %in% Rubio_dds_rlog_matrix_IAP ,]
Rubio_dds_rlog_matrix_IAP$ exp <- "Rubio"
deLorgeril_Resistant_dds_vst_matrix_IAP <- colnames(deLorgeril_Resistant_dds_vst_matrix)[colnames(deLorgeril_Resistant_dds_vst_matrix) %in% C_gig_rtracklayer_apop_product_final_IAP_list]
deLorgeril_Resistant_dds_vst_matrix_IAP <- C_gig_rtracklayer_apop_product_final[C_gig_rtracklayer_apop_product_final$transcript_id %in% deLorgeril_Resistant_dds_vst_matrix_IAP,]
deLorgeril_Resistant_dds_vst_matrix_IAP$exp <-"deLorg_Res"
deLorgeril_Susceptible_dds_vst_matrix_IAP <- colnames(deLorgeril_Susceptible_dds_vst_matrix)[colnames(deLorgeril_Susceptible_dds_vst_matrix) %in% C_gig_rtracklayer_apop_product_final_IAP_list]
deLorgeril_Susceptible_dds_vst_matrix_IAP <- C_gig_rtracklayer_apop_product_final[C_gig_rtracklayer_apop_product_final$transcript_id %in% deLorgeril_Susceptible_dds_vst_matrix_IAP,]
deLorgeril_Susceptible_dds_vst_matrix_IAP$exp <- "deLorg_Sus"
He_dds_vst_matrix_IAP <- colnames(He_dds_vst_matrix)[colnames(He_dds_vst_matrix) %in% C_gig_rtracklayer_apop_product_final_IAP_list]
He_dds_vst_matrix_IAP <- C_gig_rtracklayer_apop_product_final[C_gig_rtracklayer_apop_product_final$transcript_id %in% He_dds_vst_matrix_IAP,]
He_dds_vst_matrix_IAP$ exp <- "He"

C_gig_IAP <- rbind(Zhang_dds_rlog_matrix_IAP, Rubio_dds_rlog_matrix_IAP, 
               deLorgeril_Resistant_dds_vst_matrix_IAP, deLorgeril_Susceptible_dds_vst_matrix_IAP, He_dds_vst_matrix_IAP)

## of IAP genes from genome utilized across all the transcriptomes
C_gig_rtracklayer_apop_product_final_IAP_list[!(C_gig_rtracklayer_apop_product_final_IAP_list$transcript_id %in% C_gig_IAP$transcript_id),]
# 0 not present in the transcriptomes
C_gig_rtracklayer_apop_product_final_IAP_list[!(C_gig_rtracklayer_apop_product_final_IAP_list$gene %in% C_gig_IAP$gene),]
# 0 gene not present in the transcriptomes 

C_gig_IAP_number <- C_gig_IAP %>% group_by(exp) %>% dplyr::summarise(count = n())
C_gig_IAP_number$species <- "C_gigas"
C_gig_IAP_number$type <- "IAP"

# GIMAPs in each set 
C_gig_rtracklayer_apop_product_final_IMAP_list <- C_gig_rtracklayer_apop_product_final[grepl("IMAP",C_gig_rtracklayer_apop_product_final$product, ignore.case = TRUE),]
C_gig_rtracklayer_apop_product_final_IMAP_list <- C_gig_rtracklayer_apop_product_final_IMAP_list$transcript_id
length(C_gig_rtracklayer_apop_product_final_IMAP_list) #36

Zhang_dds_rlog_matrix_IMAP <- colnames(Zhang_dds_rlog_matrix)[colnames(Zhang_dds_rlog_matrix) %in% C_gig_rtracklayer_apop_product_final_IMAP_list]
Zhang_dds_rlog_matrix_IMAP <- C_gig_rtracklayer_apop_product_final[C_gig_rtracklayer_apop_product_final$transcript_id %in% Zhang_dds_rlog_matrix_IMAP,]
Zhang_dds_rlog_matrix_IMAP$exp <- "Zhang"
Rubio_dds_rlog_matrix_IMAP <- colnames(Rubio_dds_rlog_matrix)[colnames(Zhang_dds_rlog_matrix) %in% C_gig_rtracklayer_apop_product_final_IMAP_list]
Rubio_dds_rlog_matrix_IMAP <- C_gig_rtracklayer_apop_product_final[C_gig_rtracklayer_apop_product_final$transcript_id %in% Rubio_dds_rlog_matrix_IMAP ,]
Rubio_dds_rlog_matrix_IMAP$exp <- "Rubio"

deLorgeril_Resistant_dds_vst_matrix_IMAP <- colnames(deLorgeril_Resistant_dds_vst_matrix)[colnames(deLorgeril_Resistant_dds_vst_matrix) %in% C_gig_rtracklayer_apop_product_final_IMAP_list]
deLorgeril_Resistant_dds_vst_matrix_IMAP <- C_gig_rtracklayer_apop_product_final[C_gig_rtracklayer_apop_product_final$transcript_id %in% deLorgeril_Resistant_dds_vst_matrix_IMAP,]
deLorgeril_Resistant_dds_vst_matrix_IMAP$exp <- "deLorg_Res"
deLorgeril_Susceptible_dds_vst_matrix_IMAP <- colnames(deLorgeril_Susceptible_dds_vst_matrix)[colnames(deLorgeril_Susceptible_dds_vst_matrix) %in% C_gig_rtracklayer_apop_product_final_IMAP_list]
deLorgeril_Susceptible_dds_vst_matrix_IMAP <- C_gig_rtracklayer_apop_product_final[C_gig_rtracklayer_apop_product_final$transcript_id %in% deLorgeril_Susceptible_dds_vst_matrix_IMAP,]
deLorgeril_Susceptible_dds_vst_matrix_IMAP$exp <- "deLorg_Sus"
He_dds_vst_matrix_IMAP <- colnames(He_dds_vst_matrix)[colnames(He_dds_vst_matrix) %in% C_gig_rtracklayer_apop_product_final_IMAP_list]
He_dds_vst_matrix_IMAP <- C_gig_rtracklayer_apop_product_final[C_gig_rtracklayer_apop_product_final$transcript_id %in% He_dds_vst_matrix_IMAP,]
He_dds_vst_matrix_IMAP$exp <- "He"

C_gig_IMAP <- rbind(Zhang_dds_rlog_matrix_IMAP, Rubio_dds_rlog_matrix_IMAP, 
                   deLorgeril_Resistant_dds_vst_matrix_IMAP, deLorgeril_Susceptible_dds_vst_matrix_IMAP, He_dds_vst_matrix_IMAP)

## of GIMAP genes from genome utilized across all the transcriptomes
C_gig_rtracklayer_apop_product_final_IMAP_list[!(C_gig_rtracklayer_apop_product_final_IMAP_list$transcript_id %in% C_gig_IMAP$transcript_id),]
# two transcripts not present in the transcriptomes
C_gig_rtracklayer_apop_product_final_IMAP_list[!(C_gig_rtracklayer_apop_product_final_IMAP_list$gene %in% C_gig_IMAP$gene),]
# only 1 gene not present in the transcriptomes 

C_gig_IMAP_number <- C_gig_IMAP %>% group_by(exp) %>% dplyr::summarise(count = n())
C_gig_IMAP_number$species <- "C_gigas"
C_gig_IMAP_number$type <- "GIMAP"

## Combine C_gig IAP and GIMAP with C_vir

genome_IAP_GIMAP_number = data.frame(exp = c("C_vir_IAP","C_gig_IAP","C_vir_GIMAP","C_gig_GIMAP"), count=c("105","51","107","36"),
                                     species=c("C_vir","C_gigas","C_vir","C_gigas"), type=c("IAP", "IAP","GIMAP","GIMAP"), ref=c("ref","ref","ref","ref"))
genome_IAP_GIMAP_number$count <- as.numeric(genome_IAP_GIMAP_number$count)

C_gig_IMAP_number$ref <- "GIMAP"
C_vir_IMAP_number$ref <- "GIMAP"
C_gig_IAP_number$ref <- "IAP"
C_vir_IAP_number$ref <- "IAP"

consensus_IAP_GIMAP_number <- rbind(C_gig_IMAP_number, C_vir_IMAP_number, C_gig_IAP_number, C_vir_IAP_number,genome_IAP_GIMAP_number)

# Plot
ggplot(consensus_IAP_GIMAP_number, aes(x=exp, y =count, fill=ref)) + geom_col(position="dodge") + coord_flip() + 
  ggtitle("GIMAP and IAP transcripts in each transcript data set") + xlab("Experiment") + ylab("Number of Transcripts") + facet_grid(.~type)

#### UPSET PLOT OF SIGNIFICANT IAP AND GIMAP TRANSCRIPTS ####
# helpful tutorial for doing this: http://genomespot.blogspot.com/2017/09/upset-plots-as-replacement-to-venn.html
# http://crazyhottommy.blogspot.com/2016/01/upset-plot-for-overlapping-chip-seq.html
# https://www.littlemissdata.com/blog/set-analysis
# UpsetR vignette: https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html

C_vir_all_exp_mod_sig_apop_IAP <- C_vir_all_exp_mod_sig_apop[grepl("IAP",C_vir_all_exp_mod_sig_apop$product, ignore.case = TRUE),]
C_gig_all_exp_mod_sig_apop_IAP <- C_gig_all_exp_mod_sig_apop[grepl("IAP",C_gig_all_exp_mod_sig_apop$product, ignore.case = TRUE),]
C_vir_all_exp_mod_sig_apop_IMAP <- C_vir_all_exp_mod_sig_apop[grepl("IMAP",C_vir_all_exp_mod_sig_apop$product, ignore.case = TRUE),]
C_gig_all_exp_mod_sig_apop_IMAP <- C_gig_all_exp_mod_sig_apop[grepl("IMAP",C_gig_all_exp_mod_sig_apop$product, ignore.case = TRUE),]

C_vir_all_exp_mod_sig_apop_IAP_upset <- C_vir_all_exp_mod_sig_apop_IAP [,c("transcript_id","exp")]
C_gig_all_exp_mod_sig_apop_IAP_upset <- C_gig_all_exp_mod_sig_apop_IAP [,c("transcript_id","exp")]
C_vir_all_exp_mod_sig_apop_IMAP_upset <- C_vir_all_exp_mod_sig_apop_IMAP[,c("transcript_id","exp")]
C_gig_all_exp_mod_sig_apop_IMAP_upset <- C_gig_all_exp_mod_sig_apop_IMAP[,c("transcript_id","exp")]

# Convert into wide format using reshape
C_vir_all_exp_mod_sig_apop_IAP_upset_wide  <- C_vir_all_exp_mod_sig_apop_IAP_upset %>%mutate(value=1) %>% spread(exp, value, fill =0 )
C_gig_all_exp_mod_sig_apop_IAP_upset_wide  <- C_gig_all_exp_mod_sig_apop_IAP_upset %>%mutate(value=1) %>% spread(exp, value, fill =0 )
C_vir_all_exp_mod_sig_apop_IMAP_upset_wide <- C_vir_all_exp_mod_sig_apop_IMAP_upset %>%mutate(value=1) %>% spread(exp, value, fill =0 )
C_gig_all_exp_mod_sig_apop_IMAP_upset_wide <- C_gig_all_exp_mod_sig_apop_IMAP_upset %>%mutate(value=1) %>% spread(exp, value, fill =0 )


# Make upset plots
C_vir_all_exp_mod_sig_apop_IAP_upset_wide_GROUP <- upset(C_vir_all_exp_mod_sig_apop_IAP_upset_wide,  mainbar.y.label = "Transcript id Intersections", 
                                         sets.x.label = "C_vir IAP upset ", text.scale = c(1.3, 1.3, 1.3, 1.3, 2, 1.0), order.by="freq")

C_gig_all_exp_mod_sig_apop_IAP_upset_wide_GROUP <- upset(C_gig_all_exp_mod_sig_apop_IAP_upset_wide,  mainbar.y.label = "Transcript id Intersections", 
                                                         sets.x.label = "C_gig IAP upset ", text.scale = c(1.3, 1.3, 1.3, 1.3, 2, 1.0), order.by="freq")

C_gig_all_exp_mod_sig_apop_IMAP_upset_wide_GROUP <- upset(C_gig_all_exp_mod_sig_apop_IMAP_upset_wide,  mainbar.y.label = "Transcript id Intersections", 
                                                         sets.x.label = "C_gig IAP upset ", text.scale = c(1.3, 1.3, 1.3, 1.3, 2, 1.0), order.by="freq")

C_vir_all_exp_mod_sig_apop_IMAP_upset_wide_GROUP <- upset(C_vir_all_exp_mod_sig_apop_IMAP_upset_wide,  mainbar.y.label = "Transcript id Intersections", 
                                                          sets.x.label = "C_gig IAP upset ", text.scale = c(1.3, 1.3, 1.3, 1.3, 2, 1.0), order.by="freq")
## Extract products that appear only once 
#C_vir_apop_LFC_notshared  <- C_vir_apop_LFC[!(duplicated(C_vir_apop_LFC$transcript_id) | duplicated(C_vir_apop_LFC$transcript_id, fromLast = TRUE)), ]
#C_gig_apop_LFC_notshared <- C_gig_apop_LFC[!(duplicated(C_gig_apop_LFC$Name) | duplicated(C_gig_apop_LFC$Name, fromLast = TRUE)), ]

#### SIGNIFICANT MODULES: GENE LEVEL FREQUENCY PLOT AND HEATMAP IAP AND GIMAP ####

C_vir_all_exp_mod_sig_apop_IAP_GENE_upset <- C_vir_all_exp_mod_sig_apop_IAP [,c("gene","exp", "product")]
C_gig_all_exp_mod_sig_apop_IAP_GENE_upset <- C_gig_all_exp_mod_sig_apop_IAP [,c("gene","exp", "product")]
C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset <- C_vir_all_exp_mod_sig_apop_IMAP[,c("gene","exp",  "product")]
C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset <- C_gig_all_exp_mod_sig_apop_IMAP[,c("gene","exp",  "product")]

C_vir_all_exp_mod_sig_apop_IAP_GENE_upset$value <- as.numeric(c("1"))
C_gig_all_exp_mod_sig_apop_IAP_GENE_upset$value <- as.numeric(c("1"))
C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset$value <- as.numeric(c("1"))
C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset$value <- as.numeric(c("1"))

key= data.frame(exp=c("Pro_RE22_RE22","ROD_Res","ROD_Sus","Pro_RE22_Pro_RI", "Pro_RE22_Pro_S4","Dermo_Tol","Dermo_Sus","Probiotic"), type = c("Path_Bac","Path_Bac","Path_Bac","Probiotic",
                                                                                                                                  "Probiotic", "Dermo","Dermo","Probiotic"))
C_vir_all_exp_mod_sig_apop_IAP_GENE_upset <- left_join(C_vir_all_exp_mod_sig_apop_IAP_GENE_upset   , key)
C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset <- left_join(C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset, key)

C_gig_key = data.frame(exp= c("Zhang_LPS",    "Zhang_Vibrio", "Rubio_NV",     "Rubio_V" ,     "deLorg_Res",   "deLorg_Sus",   "He"   ), type=c("Path_Bac","Path_Bac","Nonpath_Bac",
                                                                                                                                               "Path_Bac", "Viral","Viral","Viral"  ))
C_gig_all_exp_mod_sig_apop_IAP_GENE_upset <- left_join(C_gig_all_exp_mod_sig_apop_IAP_GENE_upset  , C_gig_key)
C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset <- left_join(C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset  , C_gig_key)

## Condense product name so I can plot by it
C_vir_all_exp_mod_sig_apop_IAP_GENE_upset $product <- str_remove(C_vir_all_exp_mod_sig_apop_IAP_GENE_upset $product, ", transcript variant X1")
C_vir_all_exp_mod_sig_apop_IAP_GENE_upset $product <- str_remove(C_vir_all_exp_mod_sig_apop_IAP_GENE_upset $product, ", transcript variant X2")
C_vir_all_exp_mod_sig_apop_IAP_GENE_upset $product <- str_remove(C_vir_all_exp_mod_sig_apop_IAP_GENE_upset $product, ", transcript variant X3")
C_vir_all_exp_mod_sig_apop_IAP_GENE_upset $product <- str_remove(C_vir_all_exp_mod_sig_apop_IAP_GENE_upset $product, ", transcript variant X4")
C_vir_all_exp_mod_sig_apop_IAP_GENE_upset $product <- str_remove(C_vir_all_exp_mod_sig_apop_IAP_GENE_upset $product, ", transcript variant X7")

C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset$product <- str_remove(C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset$product, ", transcript variant X1")
C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset$product <- str_remove(C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset$product, ", transcript variant X2")
C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset$product <- str_remove(C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset$product, ", transcript variant X3")
C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset$product <- str_remove(C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset$product, ", transcript variant X4")
C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset$product <- str_remove(C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset$product, ", transcript variant X8")

C_gig_all_exp_mod_sig_apop_IAP_GENE_upset $product <- str_remove(C_gig_all_exp_mod_sig_apop_IAP_GENE_upset $product, ", transcript variant X1")
C_gig_all_exp_mod_sig_apop_IAP_GENE_upset $product <- str_remove(C_gig_all_exp_mod_sig_apop_IAP_GENE_upset $product, ", transcript variant X2")
C_gig_all_exp_mod_sig_apop_IAP_GENE_upset $product <- str_remove(C_gig_all_exp_mod_sig_apop_IAP_GENE_upset $product, ", transcript variant X3")
C_gig_all_exp_mod_sig_apop_IAP_GENE_upset $product <- str_remove(C_gig_all_exp_mod_sig_apop_IAP_GENE_upset $product, ", transcript variant X4")
C_gig_all_exp_mod_sig_apop_IAP_GENE_upset $product <- str_remove(C_gig_all_exp_mod_sig_apop_IAP_GENE_upset $product, ", transcript variant X5")

C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset$product <- str_remove(C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset$product, ", transcript variant X1")
C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset$product <- str_remove(C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset$product, ", transcript variant X2")
C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset$product <- str_remove(C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset$product, ", transcript variant X3")
C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset$product <- str_remove(C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset$product, ", transcript variant X4")

# plot frequency across experiments
C_vir_all_exp_mod_sig_apop_IAP_GENE_upset_freq_plot <- ggplot( C_vir_all_exp_mod_sig_apop_IAP_GENE_upset, aes(x=product, y =value, fill=exp)) + geom_col() + coord_flip() + theme(axis.text.x = element_text(hjust=1.0, angle=90))  + ylab("number of genes") + xlab("gene ID") + ggtitle("C_vir IAP shared genes across transcriptomes")
#ggplot( C_vir_all_exp_mod_sig_apop_IAP_GENE_upset, aes(x=gene, y =value, fill=type)) + geom_col() + coord_flip() + theme(axis.text.x = element_text(hjust=1.0, angle=90))  + ylab("number of genes") + xlab("gene ID") + ggtitle("C_vir IAP shared genes across transcriptomes")

C_gig_all_exp_mod_sig_apop_IAP_GENE_upset_freq_plot <- ggplot( C_gig_all_exp_mod_sig_apop_IAP_GENE_upset, aes(x=product, y =value, fill=exp)) + geom_col() + coord_flip() + theme(axis.text.x = element_text(hjust=1.0, angle=90))  + ylab("number of genes") + xlab("gene ID") + ggtitle("C_gig IAP shared genes across transcriptomes")
#ggplot( C_gig_all_exp_mod_sig_apop_IAP_GENE_upset, aes(x=gene, y =value, fill=type)) + geom_col() + coord_flip() + theme(axis.text.x = element_text(hjust=1.0, angle=90))  + ylab("number of genes") + xlab("gene ID") + ggtitle("C_gig IAP shared genes across transcriptomes")

C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset_freq_plot <- ggplot( C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset, aes(x=product, y =value, fill=exp)) + geom_col() + coord_flip() + theme(axis.text.x = element_text(hjust=1.0, angle=90)) + ylab("number of genes") + xlab("gene ID") + ggtitle("C_vir GIMAP shared genes across transcriptomes")
#ggplot( C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset, aes(x=gene, y =value, fill=type)) + geom_col() + coord_flip() + theme(axis.text.x = element_text(hjust=1.0, angle=90)) + ylab("number of genes") + xlab("gene ID") + ggtitle("C_vir GIMAP shared genes across transcriptomes")

C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset_freq_plot <-ggplot( C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset, aes(x=product, y =value, fill=exp)) + geom_col() + coord_flip() + theme(axis.text.x = element_text(hjust=1.0, angle=90)) + ylab("number of genes") + xlab("gene ID") + ggtitle("C_gig GIMAP shared genes across transcriptomes")
#ggplot( C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset, aes(x=gene, y =value, fill=type)) + geom_col() + coord_flip() + theme(axis.text.x = element_text(hjust=1.0, angle=90)) + ylab("number of genes") + xlab("gene ID") + ggtitle("C_gig GIMAP shared genes across transcriptomes")

## Percentage plot 
C_vir_all_exp_mod_sig_apop_IAP_GENE_upset_percent_plot <- ggplot( C_vir_all_exp_mod_sig_apop_IAP_GENE_upset, aes(x=exp, y =value, fill=product)) + geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(hjust=1.0, angle=90))  + ylab("number of genes") + xlab("gene ID") + ggtitle("C_vir IAP shared genes across transcriptomes WGCNA")

C_gig_all_exp_mod_sig_apop_IAP_GENE_upset_percent_plot <- ggplot(C_gig_all_exp_mod_sig_apop_IAP_GENE_upset, aes(x=exp, y =value, fill=product)) + geom_bar(position="fill", stat="identity") + 
  theme(axis.text.x = element_text(hjust=1.0, angle=90))  + ylab("number of genes") + xlab("gene ID") + ggtitle("C_gig IAP shared genes across transcriptomes WGCNA")

C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset_percent_plot <- ggplot( C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset, aes(x=exp, y =value, fill=product)) + geom_bar(position="fill", stat="identity") + 
  theme(axis.text.x = element_text(hjust=1.0, angle=90))  + ylab("number of genes") + xlab("gene ID") + ggtitle("C_vir GIMAP shared genes across transcriptomes WGCNA")

C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset_percent_plot <- ggplot( C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset, aes(x=exp, y =value, fill=product)) + geom_bar(position="fill", stat="identity") + 
  theme(axis.text.x = element_text(hjust=1.0, angle=90))  + ylab("number of genes") + xlab("gene ID") + ggtitle("C_gig GIMAP shared genes across transcriptomes WGCNA")

### PLOT AS HEATMAP ###
C_vir_all_exp_mod_sig_apop_IAP_GENE_upset$gene <- factor(C_vir_all_exp_mod_sig_apop_IAP_GENE_upset$gene , levels=unique((C_vir_all_exp_mod_sig_apop_IAP_GENE_upset$gene )[order(C_vir_all_exp_mod_sig_apop_IAP_GENE_upset$product)]))
C_gig_all_exp_mod_sig_apop_IAP_GENE_upset$gene <- factor(C_gig_all_exp_mod_sig_apop_IAP_GENE_upset$gene, levels=unique((C_gig_all_exp_mod_sig_apop_IAP_GENE_upset$gene)[order(C_gig_all_exp_mod_sig_apop_IAP_GENE_upset$product)]))
C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset$gene <- factor(C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset$gene, levels=unique((C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset$gene)[order(C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset$product)]))
C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset$gene <- factor(C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset$gene, levels=unique((C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset$gene)[order(C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset$product)]))

C_vir_all_exp_mod_sig_apop_IAP_GENE_upset$exp <- factor(C_vir_all_exp_mod_sig_apop_IAP_GENE_upset$exp, levels=c("Probiotic","Pro_RE22_Pro_S4","Pro_RE22_Pro_RI","Pro_RE22_RE22","ROD_Res","ROD_Sus","Dermo_Tol","Dermo_Sus"))
C_gig_all_exp_mod_sig_apop_IAP_GENE_upset$exp <- factor(C_gig_all_exp_mod_sig_apop_IAP_GENE_upset$exp,   levels=c("Zhang_LPS","Zhang_Vibrio","Rubio_NV","Rubio_V","deLorg_Res","deLorg_Sus","He"))
C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset$exp <- factor(C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset$exp, levels=c("Probiotic","Pro_RE22_Pro_S4","Pro_RE22_Pro_RI","Pro_RE22_RE22","ROD_Res","ROD_Sus","Dermo_Tol","Dermo_Sus"))
C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset$exp <- factor(C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset$exp, levels=c("Zhang_LPS","Zhang_Vibrio","Rubio_NV","Rubio_V","deLorg_Res","deLorg_Sus","He"))

C_vir_all_exp_mod_sig_apop_IAP_GENE_upset_tile <- ggplot(C_vir_all_exp_mod_sig_apop_IAP_GENE_upset, aes(x=gene, y =exp, fill=product)) + 
  geom_tile()  + 
  coord_flip() + 
  scale_fill_brewer(palette = "PiYG")+
  theme_minimal() +
  labs(y="Experiment", x="Gene ID",
         title ="C. virginica IAP Gene Significance Across Experiments")
  
C_gig_all_exp_mod_sig_apop_IAP_GENE_upset_tile <- ggplot(C_gig_all_exp_mod_sig_apop_IAP_GENE_upset, aes(x=gene, y =exp, fill=product)) + 
  geom_tile()  + 
  coord_flip() + 
  scale_fill_brewer(palette = "PiYG")+
  theme_minimal() +
  labs(y="Experiment", x="Gene ID",
         title="C. gigas IAP Gene Significance Across Experiments")

C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset_tile <- ggplot(C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset, aes(x=gene, y =exp, fill=product)) + 
  geom_tile()  + 
  coord_flip() + 
  scale_fill_brewer(palette = "PiYG")+
  theme_minimal() +
  labs(y="Experiment", x="Gene ID",
         title="C. virginica GIMAP Gene Significance Across Experiments")

C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset_tile <- ggplot(C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset, aes(x=gene, y =exp, fill=product)) + 
  geom_tile()  + 
  coord_flip() + 
  scale_fill_brewer(palette = "PiYG")+
  theme_minimal() +
  labs(y="Experiment", x="Gene ID",
         title="C. gigas GIMAP Gene Significance Across Experiments")
# IAP compare
cowplot::plot_grid(C_vir_all_exp_mod_sig_apop_IAP_GENE_upset_tile, C_gig_all_exp_mod_sig_apop_IAP_GENE_upset_tile,
                   align="hv", nrow=2)
# GIMAP compare 
cowplot::plot_grid(C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset_tile,C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset_tile,
                   align="hv", nrow=2)

#### REPEAT ABOVE WITH THE DESEQ DATA 
C_vir_apop_LFC_IAP <- C_vir_apop_LFC[grepl("IAP",C_vir_apop_LFC$product, ignore.case = TRUE),]
C_gig_apop_LFC_IAP <- C_gig_apop_LFC[grepl("IAP",C_gig_apop_LFC$product, ignore.case = TRUE),]
C_vir_apop_LFC_IMAP <- C_vir_apop_LFC[grepl("IMAP",C_vir_apop_LFC$product, ignore.case = TRUE),]
C_gig_apop_LFC_IMAP <- C_gig_apop_LFC[grepl("IMAP",C_gig_apop_LFC$product, ignore.case = TRUE),]

C_vir_apop_LFC_IAP_upset <- C_vir_apop_LFC_IAP [,c("gene","experiment",  "product","group_by_sim")]
C_gig_apop_LFC_IAP_upset <- C_gig_apop_LFC_IAP [,c("gene","experiment",  "product","group_by_sim")]
C_vir_apop_LFC_IMAP_upset <-C_vir_apop_LFC_IMAP[,c("gene","experiment",  "product","group_by_sim")]
C_gig_apop_LFC_IMAP_upset <-C_gig_apop_LFC_IMAP[,c("gene","experiment",  "product","group_by_sim")]

C_vir_apop_LFC_IAP_upset $value <- as.numeric(c("1"))
C_gig_apop_LFC_IAP_upset $value <- as.numeric(c("1"))
C_vir_apop_LFC_IMAP_upset$value <- as.numeric(c("1"))
C_gig_apop_LFC_IMAP_upset$value <- as.numeric(c("1"))

## Condense product name so I can plot by it
C_vir_apop_LFC_IAP_upset $product <- gsub("\\,.*","", C_vir_apop_LFC_IAP_upset $product )
C_gig_apop_LFC_IAP_upset $product <- gsub("\\,.*","", C_gig_apop_LFC_IAP_upset $product )
C_vir_apop_LFC_IMAP_upset$product <- gsub("\\,.*","", C_vir_apop_LFC_IMAP_upset$product )
C_gig_apop_LFC_IMAP_upset$product <- gsub("\\,.*","", C_gig_apop_LFC_IMAP_upset$product )

# Plot as frequency table and compare with WGCNA
## General frequency patterns seem the same 
C_vir_apop_LFC_IAP_upset_freq_plot <- ggplot( C_vir_apop_LFC_IAP_upset, aes(x=product, y =value, fill=group_by_sim)) + geom_col() + coord_flip() + 
theme(axis.text.x = element_text(hjust=1.0, angle=90))  + ylab("number of genes") + xlab("gene ID") + ggtitle("C_vir IAP shared genes across transcriptomes DESeq")
cowplot::plot_grid(C_vir_all_exp_mod_sig_apop_IAP_GENE_upset_freq_plot , C_vir_apop_LFC_IAP_upset_freq_plot, align ="hv",nrow=2)

C_gig_apop_LFC_IAP_upset_freq_plot <- ggplot( C_gig_apop_LFC_IAP_upset, aes(x=product, y =value, fill=group_by_sim)) + geom_col() + coord_flip() + 
  theme(axis.text.x = element_text(hjust=1.0, angle=90))  + ylab("number of genes") + xlab("gene ID") + ggtitle("C_gig IAP shared genes across transcriptomes DESeq")
cowplot::plot_grid(C_gig_all_exp_mod_sig_apop_IAP_GENE_upset_freq_plot , C_gig_apop_LFC_IAP_upset_freq_plot, align ="hv",nrow=2)

C_gig_apop_LFC_GIMAP_upset_freq_plot <- ggplot( C_gig_apop_LFC_IMAP_upset, aes(x=product, y =value, fill=group_by_sim)) + geom_col() + coord_flip() + 
  theme(axis.text.x = element_text(hjust=1.0, angle=90))  + ylab("number of genes") + xlab("gene ID") + ggtitle("C_gig GIMAP shared genes across transcriptomes DESeq")
cowplot::plot_grid(C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset_freq_plot , C_gig_apop_LFC_GIMAP_upset_freq_plot, align ="hv",nrow=2)

C_vir_apop_LFC_GIMAP_upset_freq_plot <- ggplot( C_vir_apop_LFC_IMAP_upset, aes(x=product, y =value, fill=group_by_sim)) + geom_col() + coord_flip() + 
  theme(axis.text.x = element_text(hjust=1.0, angle=90))  + ylab("number of genes") + xlab("gene ID") + ggtitle("C_vir GIMAP shared genes across transcriptomes DESeq")
cowplot::plot_grid(C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset_freq_plot , C_vir_apop_LFC_GIMAP_upset_freq_plot, align ="hv",nrow=2)

## Percentage plot 
# Need to go back in and explicitly set colors 
C_vir_apop_LFC_IAP_upset_percent_plot <- ggplot( C_vir_apop_LFC_IAP_upset, aes(x=group_by_sim, y =value, fill=product)) + geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(hjust=1.0, angle=90))  + ylab("number of genes") + xlab("gene ID") + ggtitle("C_vir IAP shared genes across transcriptomes DESeq")

cowplot::plot_grid(C_vir_apop_LFC_IAP_upset_percent_plot , C_vir_all_exp_mod_sig_apop_IAP_GENE_upset_percent_plot , align ="hv",nrow=2)

C_gig_apop_LFC_IAP_upset_percent_plot <- ggplot( C_gig_apop_LFC_IAP_upset, aes(x=group_by_sim, y =value, fill=product)) + geom_bar(position="fill", stat="identity") + 
  theme(axis.text.x = element_text(hjust=1.0, angle=90))  + ylab("number of genes") + xlab("gene ID") + ggtitle("C_gig IAP shared genes across transcriptomes DESeq")
cowplot::plot_grid(C_gig_apop_LFC_IAP_upset_percent_plot , C_gig_all_exp_mod_sig_apop_IAP_GENE_upset_percent_plot , align ="hv",nrow=2)

C_gig_apop_LFC_GIMAP_upset_percent_plot <- ggplot( C_gig_apop_LFC_IMAP_upset, aes(x=group_by_sim, y =value, fill=product)) + geom_bar(position="fill", stat="identity") + 
  theme(axis.text.x = element_text(hjust=1.0, angle=90))  + ylab("number of genes") + xlab("gene ID") + ggtitle("C_gig GIMAP shared genes across transcriptomes DESeq")

cowplot::plot_grid(C_gig_apop_LFC_GIMAP_upset_percent_plot  ,  C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset_percent_plot, align ="hv",nrow=2)

C_vir_apop_LFC_GIMAP_upset_percent_plot <- ggplot( C_vir_apop_LFC_IMAP_upset, aes(x=group_by_sim, y =value, fill=product)) + geom_bar(position="fill", stat="identity") + 
  theme(axis.text.x = element_text(hjust=1.0, angle=90))  + ylab("number of genes") + xlab("gene ID") + ggtitle("C_vir GIMAP shared genes across transcriptomes DESeq")
cowplot::plot_grid(C_vir_apop_LFC_GIMAP_upset_percent_plot , C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset_percent_plot, align ="hv",nrow=2)

# refactor data so I can plot as heatmap
C_vir_apop_LFC_IAP_upset$gene <- factor(C_vir_apop_LFC_IAP_upset$gene , levels=unique((C_vir_apop_LFC_IAP_upset$gene )[order(C_vir_apop_LFC_IAP_upset$product)]))
C_gig_apop_LFC_IAP_upset $gene <- factor(C_gig_apop_LFC_IAP_upset $gene, levels=unique((C_gig_apop_LFC_IAP_upset $gene)[order(C_gig_apop_LFC_IAP_upset $product)]))
C_vir_apop_LFC_IMAP_upset$gene <- factor(C_vir_apop_LFC_IMAP_upset$gene, levels=unique((C_vir_apop_LFC_IMAP_upset$gene)[order(C_vir_apop_LFC_IMAP_upset$product)]))
C_gig_apop_LFC_IMAP_upset$gene <- factor(C_gig_apop_LFC_IMAP_upset$gene, levels=unique((C_gig_apop_LFC_IMAP_upset$gene)[order(C_gig_apop_LFC_IMAP_upset$product)]))

## need to reorder the experiments here

C_vir_apop_LFC_IAP_upset_tile <- ggplot(C_vir_apop_LFC_IAP_upset, aes(x=gene, y =group_by_sim, fill=product)) + 
  geom_tile()  + 
  coord_flip() + 
  scale_fill_brewer(palette = "PiYG")+
  theme_minimal() +
  labs(y="Experiment", x="Gene ID",
       title ="C. virginica IAP Gene Significance Across Experiments DESeq")+theme(axis.text.x = element_text(hjust=1.0, angle=90))

C_gig_apop_LFC_IAP_upset_tile <- ggplot(C_gig_apop_LFC_IAP_upset, aes(x=gene, y =group_by_sim, fill=product)) + 
  geom_tile()  + 
  coord_flip() + 
  scale_fill_brewer(palette = "PiYG")+
  theme_minimal() +
  labs(y="Experiment", x="Gene ID",
       title="C. gigas IAP Gene Significance Across Experiments DESeq") + theme(axis.text.x = element_text(hjust=1.0, angle=90))

C_vir_apop_LFC_IMAP_upset_tile <- ggplot(C_vir_apop_LFC_IMAP_upset, aes(x=gene, y =group_by_sim, fill=product)) + 
  geom_tile()  + 
  coord_flip() + 
  scale_fill_brewer(palette = "PiYG")+
  theme_minimal() +
  labs(y="Experiment", x="Gene ID",
       title="C. virginica GIMAP Gene Significance Across Experiments DESeq") +theme(axis.text.x = element_text(hjust=1.0, angle=90))

C_gig_apop_LFC_IMAP_upset_tile <- ggplot(C_gig_apop_LFC_IMAP_upset, aes(x=gene, y =group_by_sim, fill=product)) + 
  geom_tile()  + 
  coord_flip() + 
  scale_fill_brewer(palette = "PiYG")+
  theme_minimal() +
  labs(y="Experiment", x="Gene ID",
       title="C. gigas GIMAP Gene Significance Across Experiments DESeq") + theme(axis.text.x = element_text(hjust=1.0, angle=90))
# IAP compare
cowplot::plot_grid(C_vir_apop_LFC_IAP_upset_tile, C_gig_apop_LFC_IAP_upset_tile ,
                   align="hv", nrow=2)
# GIMAP compare 
cowplot::plot_grid(C_vir_apop_LFC_IMAP_upset_tile,C_gig_apop_LFC_IMAP_upset_tile,
                   align="hv", nrow=2)

## COMPARE WGCNA AND DESEQ2 ##
cowplot::plot_grid( C_vir_all_exp_mod_sig_apop_IAP_GENE_upset_tile,C_vir_apop_LFC_IAP_upset_tile,
                   align="hv", nrow=2)
# GIMAP compare 
cowplot::plot_grid( C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset_tile,C_vir_apop_LFC_IMAP_upset_tile,
                   align="hv", nrow=2)

cowplot::plot_grid(C_gig_all_exp_mod_sig_apop_IAP_GENE_upset_tile, C_gig_apop_LFC_IAP_upset_tile,
                   align="hv", nrow=2)
# GIMAP compare 
cowplot::plot_grid(C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset_tile,C_gig_apop_LFC_IMAP_upset_tile,
                   align="hv", nrow=2)

#### VENNDIAGRAM OF GENES BETWEEN WGCNA AND DESEQ2 DATA ####

## COMPARE SIGNIFICANT GENES BETWEEN WGNCA AND IAP
C_vir_apop_LFC_IAP_upset
C_gig_apop_LFC_IAP_upset
C_vir_apop_LFC_IMAP_upset
C_gig_apop_LFC_IMAP_upset
C_vir_all_exp_mod_sig_apop_IAP_GENE_upset
C_gig_all_exp_mod_sig_apop_IAP_GENE_upset
C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset
C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset

# Basic Venndiagram 
venn.diagram(
  x = list(unique(C_vir_apop_LFC_IAP$gene), unique(C_vir_all_exp_mod_sig_apop_IAP_GENE_upset$gene)),
  category.names = c("C_vir IAP DESeq" , "C_vir IAP WGCNA"),
  filename="C_vir_IAP_DESeq_WGCNA_comparison.png"
)

#### PERCENTAGE PLOT OF SIGNIFICANT GENES FOR EACH GENE MEMBER ####

C_vir_apop_LFC_IAP_upset
C_gig_apop_LFC_IAP_upset
C_vir_apop_LFC_IMAP_upset
C_gig_apop_LFC_IMAP_upset
C_vir_all_exp_mod_sig_apop_IAP_GENE_upset
C_gig_all_exp_mod_sig_apop_IAP_GENE_upset
C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset
C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset

### FRACTION OF IAP AND GIMAP GENES SIGNIFICANT WITH CHALLENGE ####

## For WGCNA data
C_gig_sig_apop_fraction_IAP <- C_gig_rtracklayer_apop_product_final_IAP_list[(C_gig_rtracklayer_apop_product_final_IAP_list$gene %in%  C_gig_all_exp_mod_sig_apop_IAP_GENE_upset$gene),]
length(unique(C_gig_sig_apop_fraction_IAP$gene)) # 23

C_vir_sig_apop_fraction_IAP <- C_vir_rtracklayer_apop_product_final_IAP_list[(C_vir_rtracklayer_apop_product_final_IAP_list$gene %in%  C_vir_all_exp_mod_sig_apop_IAP_GENE_upset$gene),]
length(unique(C_vir_sig_apop_fraction_IAP$gene)) # 28

C_gig_sig_apop_fraction_GIMAP <- C_gig_rtracklayer_apop_product_final_IMAP_list[(C_gig_rtracklayer_apop_product_final_IMAP_list$gene %in%  C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset$gene),]
length(unique(C_gig_sig_apop_fraction_GIMAP$gene)) # 11

C_vir_sig_apop_fraction_GIMAP <- C_vir_rtracklayer_apop_product_final_IMAP_list[(C_vir_rtracklayer_apop_product_final_IMAP_list$gene %in%  C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset$gene),]
length(unique(C_vir_sig_apop_fraction_GIMAP$gene)) # 16

## For DESeq2 data
C_gig_apop_LFC_IAP_fraction <- C_gig_rtracklayer_apop_product_final_IAP_list[(C_gig_rtracklayer_apop_product_final_IAP_list$gene %in%  C_gig_apop_LFC_IAP$gene),]
length(unique(C_gig_apop_LFC_IAP_fraction$gene)) #18

C_vir_apop_LFC_IAP_fraction <- C_vir_rtracklayer_apop_product_final_IAP_list[(C_vir_rtracklayer_apop_product_final_IAP_list$gene %in%  C_vir_apop_LFC_IAP$gene),]
length(unique(C_vir_apop_LFC_IAP_fraction$gene)) # 20

C_gig_apop_LFC_GIMAP_fraction <- C_gig_rtracklayer_apop_product_final_IMAP_list[(C_gig_rtracklayer_apop_product_final_IMAP_list$gene %in%  C_gig_apop_LFC_IMAP$gene),]
length(unique(C_gig_apop_LFC_GIMAP_fraction$gene)) # 12

C_vir_apop_LFC_GIMAP_fraction <- C_vir_rtracklayer_apop_product_final_IMAP_list[(C_vir_rtracklayer_apop_product_final_IMAP_list$gene %in%  C_vir_apop_LFC_IMAP$gene),]
length(unique(C_vir_apop_LFC_GIMAP_fraction$gene)) # 5

####  COMPARING SIGNIFICANT IAP AND GIMAP GENES BETWEEN EXPERIMENTS ####

### FIND C VIR IAP AND GIMAP SHARED BY ALL EXPERIMENTS
## SIG IAPS
C_vir_full_all_exp_mod_sig_apop_IAP <- C_vir_full_all_exp_mod_sig_apop[grepl("IAP", C_vir_full_all_exp_mod_sig_apop$product, ignore.case=TRUE),]

# Find IAPS shared by all experiments
C_vir_full_all_exp_mod_sig_apop_IAP_split <- split(C_vir_full_all_exp_mod_sig_apop_IAP$gene, C_vir_full_all_exp_mod_sig_apop_IAP$exp)
Reduce(intersect, C_vir_full_all_exp_mod_sig_apop_IAP_split)
# 0 

## SIG GIMAPS 
C_vir_full_all_exp_mod_sig_apop_GIMAP <- C_vir_full_all_exp_mod_sig_apop[grepl("IMAP", C_vir_full_all_exp_mod_sig_apop$product, ignore.case=TRUE),]

# Find GIMAPS shared by all experiments
C_vir_full_all_exp_mod_sig_apop_GIMAP_split <- split(C_vir_full_all_exp_mod_sig_apop_GIMAP$gene, C_vir_full_all_exp_mod_sig_apop_GIMAP$exp)
Reduce(intersect, C_vir_full_all_exp_mod_sig_apop_GIMAP_split)
# 0 

### FIND SHARED IAP AND GIMAP ACROSS EXPERIMENT TYPE
# Using code from the section above this which has it joined by type 
C_vir_all_exp_mod_sig_apop_IAP_GENE_upset_split <- split(C_vir_all_exp_mod_sig_apop_IAP_GENE_upset$gene, C_vir_all_exp_mod_sig_apop_IAP_GENE_upset$type)
Reduce(intersect, C_vir_all_exp_mod_sig_apop_IAP_GENE_upset_split)

C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset_split <- split(C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset$gene, C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset$type)
Reduce(intersect, C_vir_all_exp_mod_sig_apop_IMAP_GENE_upset_split)

### FIND UNIQUE IAP AND GIMAP IN EACH EXPERIMENT 
# table of total per experiment
View(C_vir_full_all_exp_mod_sig_apop_IAP %>% group_by(exp) %>% distinct(gene) %>%  dplyr::summarise(n=n()))
View(C_vir_full_all_exp_mod_sig_apop_GIMAP %>% group_by(exp) %>% distinct(gene) %>%  dplyr::summarise(n=n()))

C_vir_full_all_exp_mod_sig_apop_IAP_minus_Dermo_Tol <- C_vir_full_all_exp_mod_sig_apop_IAP %>% filter(exp != "Dermo_Tol")
C_vir_full_all_exp_mod_sig_apop_IAP_minus_Dermo_Sus <- C_vir_full_all_exp_mod_sig_apop_IAP %>% filter(exp != "Dermo_Sus")
C_vir_full_all_exp_mod_sig_apop_IAP_minus_ROD_Res <- C_vir_full_all_exp_mod_sig_apop_IAP %>% filter(exp != "ROD_Res")
C_vir_full_all_exp_mod_sig_apop_IAP_minus_Rod_Sus <- C_vir_full_all_exp_mod_sig_apop_IAP %>% filter(exp != "ROD_Sus")
C_vir_full_all_exp_mod_sig_apop_IAP_minus_Probiotic <- C_vir_full_all_exp_mod_sig_apop_IAP %>% filter(exp != "Probiotic")
C_vir_full_all_exp_mod_sig_apop_IAP_minus_Pro_RE22_S4 <- C_vir_full_all_exp_mod_sig_apop_IAP %>% filter(exp != "Pro_RE22_Pro_S4")
C_vir_full_all_exp_mod_sig_apop_IAP_minus_Pro_RE22_RI <- C_vir_full_all_exp_mod_sig_apop_IAP %>% filter(exp != "Pro_RE22_Pro_RI")
C_vir_full_all_exp_mod_sig_apop_IAP_minus_Pro_RE22_RE22 <- C_vir_full_all_exp_mod_sig_apop_IAP %>% filter(exp != "Pro_RE22_RE22_full")

C_vir_full_all_exp_mod_sig_apop_GIMAP_minus_Dermo_Tol <-     C_vir_full_all_exp_mod_sig_apop_GIMAP %>% filter(exp != "Dermo_Tol")
C_vir_full_all_exp_mod_sig_apop_GIMAP_minus_Dermo_Sus <-     C_vir_full_all_exp_mod_sig_apop_GIMAP %>% filter(exp != "Dermo_Sus")
C_vir_full_all_exp_mod_sig_apop_GIMAP_minus_ROD_Res <-       C_vir_full_all_exp_mod_sig_apop_GIMAP %>% filter(exp != "ROD_Res")
C_vir_full_all_exp_mod_sig_apop_GIMAP_minus_Rod_Sus <-       C_vir_full_all_exp_mod_sig_apop_GIMAP %>% filter(exp != "ROD_Sus")
C_vir_full_all_exp_mod_sig_apop_GIMAP_minus_Probiotic <-     C_vir_full_all_exp_mod_sig_apop_GIMAP %>% filter(exp != "Probiotic")
C_vir_full_all_exp_mod_sig_apop_GIMAP_minus_Pro_RE22_S4 <-   C_vir_full_all_exp_mod_sig_apop_GIMAP %>% filter(exp != "Pro_RE22_Pro_S4")
C_vir_full_all_exp_mod_sig_apop_GIMAP_minus_Pro_RE22_RI <-   C_vir_full_all_exp_mod_sig_apop_GIMAP %>% filter(exp != "Pro_RE22_Pro_RI")
C_vir_full_all_exp_mod_sig_apop_GIMAP_minus_Pro_RE22_RE22 <- C_vir_full_all_exp_mod_sig_apop_GIMAP %>% filter(exp != "Pro_RE22_RE22_full")

# select IAP and GIMAPs from apop dataframes
Dermo_Sus_full_module_apop_df_IAP <- Dermo_Sus_full_module_apop_df[grepl("IAP", Dermo_Sus_full_module_apop_df$product, ignore.case = TRUE),]
ROD_Res_full_module_apop_df_IAP <- ROD_Res_full_module_apop_df[grepl("IAP",ROD_Res_full_module_apop_df$product, ignore.case = TRUE),]
ROD_Sus_full_module_apop_df_IAP <- ROD_Sus_full_module_apop_df[grepl("IAP",ROD_Sus_full_module_apop_df$product, ignore.case = TRUE),]
Probiotic_full_module_apop_df_IAP <-       Probiotic_full_module_apop_df[grepl("IAP",Probiotic_full_module_apop_df$product, ignore.case = TRUE),]
Pro_RE22_Pro_full_module_apop_df_IAP <-    Pro_RE22_Pro_full_module_apop_df[grepl("IAP",Pro_RE22_Pro_full_module_apop_df$product, ignore.case = TRUE),]
Pro_RE22_Pro_RI_full_module_apop_df_IAP <- Pro_RE22_Pro_RI_full_module_apop_df[grepl("IAP", Pro_RE22_Pro_RI_full_module_apop_df$product, ignore.case = TRUE),]

Dermo_Sus_full_module_apop_df_GIMAP <- Dermo_Sus_full_module_apop_df[grepl("IMAP", Dermo_Sus_full_module_apop_df$product, ignore.case = TRUE),]
ROD_Res_full_module_apop_df_GIMAP <- ROD_Res_full_module_apop_df[grepl("IMAP",ROD_Res_full_module_apop_df$product, ignore.case = TRUE),]
ROD_Sus_full_module_apop_df_GIMAP <- ROD_Sus_full_module_apop_df[grepl("IMAP",ROD_Sus_full_module_apop_df$product, ignore.case = TRUE),]
Probiotic_full_module_apop_df_GIMAP <-       Probiotic_full_module_apop_df[grepl("IMAP",Probiotic_full_module_apop_df$product, ignore.case = TRUE),]
Pro_RE22_Pro_full_module_apop_df_GIMAP <-    Pro_RE22_Pro_full_module_apop_df[grepl("IMAP",Pro_RE22_Pro_full_module_apop_df$product, ignore.case = TRUE),]
Pro_RE22_Pro_RI_full_module_apop_df_GIMAP <- Pro_RE22_Pro_RI_full_module_apop_df[grepl("IMAP", Pro_RE22_Pro_RI_full_module_apop_df$product, ignore.case = TRUE),]

Dermo_Tol_full_module_apop_df_IAP[!(Dermo_Tol_full_module_apop_df_IAP$gene %in% C_vir_full_all_exp_mod_sig_apop_IAP_minus_Dermo_Tol$gene), ]
#64 baculoviral IAP repeat-containing protein 3-like, transcript variant X1 XM_022434497.1 LOC111101864      MEcyan                     0.4459765 0.01350638 Dermo_Tol
#34                        baculoviral IAP repeat-containing protein 7-like XM_022432289.1 LOC111100417 MEturquoise                    -0.3944922 0.03098507 Dermo_Tol

Dermo_Sus_full_module_apop_df_IAP[!(Dermo_Sus_full_module_apop_df_IAP$gene %in% C_vir_full_all_exp_mod_sig_apop_IAP_minus_Dermo_Sus$gene),]
#0
ROD_Res_full_module_apop_df_IAP[!(ROD_Res_full_module_apop_df_IAP$gene %in% C_vir_full_all_exp_mod_sig_apop_IAP_minus_ROD_Res$gene),]
#0
ROD_Sus_full_module_apop_df_IAP[!(ROD_Sus_full_module_apop_df_IAP$gene %in% C_vir_full_all_exp_mod_sig_apop_IAP_minus_Rod_Sus$gene),]
#product  transcript_id         gene   mod_names mod_signif     pvalue     exp
#124 baculoviral IAP repeat-containing protein 7-like, transcript variant X3 XM_022439819.1 LOC111105494 MEturquoise  0.9487006 0.05129943 ROD_Sus

Probiotic_full_module_apop_df_IAP[!(Probiotic_full_module_apop_df_IAP$gene %in% C_vir_full_all_exp_mod_sig_apop_IAP_minus_Probiotic$gene),]
#product  transcript_id         gene   mod_names mod_signif      pvalue       exp
#31 baculoviral IAP repeat-containing protein 2-like, transcript variant X2 XM_022432349.1 LOC111100443 MElightcyan  0.9500564 0.003679259 Probiotic

Pro_RE22_Pro_full_module_apop_df_IAP[!(Pro_RE22_Pro_full_module_apop_df_IAP$gene %in% C_vir_full_all_exp_mod_sig_apop_IAP_minus_Pro_RE22_S4$gene),]
#13  baculoviral IAP repeat-containing protein 2-like, transcript variant X1 XM_022466570.1 LOC111123894   MEroyalblue -0.6238093 0.012949352 Pro_RE22_Pro_S4
#38                         baculoviral IAP repeat-containing protein 3-like XM_022437326.1 LOC111103816 MEgreenyellow -0.7069681 0.003206386 Pro_RE22_Pro_S4
#115 baculoviral IAP repeat-containing protein 2-like, transcript variant X2 XM_022437360.1 LOC111103826  MElightpink4 -0.6601048 0.007404182 Pro_RE22_Pro_S4
#168                        baculoviral IAP repeat-containing protein 8-like XM_022436216.1 LOC111103158      MEmaroon  0.6594826 0.007479934 Pro_RE22_Pro_S4

Pro_RE22_Pro_RI_full_module_apop_df_IAP[!(Pro_RE22_Pro_RI_full_module_apop_df_IAP$gene %in% C_vir_full_all_exp_mod_sig_apop_IAP_minus_Pro_RE22_RI$gene),]
#7   baculoviral IAP repeat-containing protein 2-like, transcript variant X1 XM_022433267.1 LOC111101035  MElightgreen  0.5201691 0.046842984 Pro_RE22_Pro_RI
#65  baculoviral IAP repeat-containing protein 2-like, transcript variant X1 XM_022435920.1 LOC111102964    MEskyblue3  0.7315844 0.001936156 Pro_RE22_Pro_RI
#175                        baculoviral IAP repeat-containing protein 7-like XM_022432288.1 LOC111100416        MEblue -0.7299198 0.002006628 Pro_RE22_Pro_RI
#176 baculoviral IAP repeat-containing protein 3-like, transcript variant X1 XM_022433225.1 LOC111101018        MEblue -0.7299198 0.002006628 Pro_RE22_Pro_RI
#199 baculoviral IAP repeat-containing protein 7-like, transcript variant X1 XM_022432712.1 LOC111100633     MEmagenta -0.5713401 0.026091580 Pro_RE22_Pro_RI
#230 baculoviral IAP repeat-containing protein 6-like, transcript variant X1 XM_022477206.1 LOC111130310 MElightyellow -0.6420686 0.009859033 Pro_RE22_Pro_RI
#233 baculoviral IAP repeat-containing protein 2-like, transcript variant X2 XM_022433268.1 LOC111101035 MElightyellow -0.6420686 0.009859033 Pro_RE22_Pro_RI

Pro_RE22_RE22_full_module_apop_df_IAP[!(Pro_RE22_RE22_full_module_apop_df_IAP$gene %in% C_vir_full_all_exp_mod_sig_apop_IAP_minus_Pro_RE22_RE22$gene),]
#98    baculoviral IAP repeat-containing protein 2-like, transcript variant X1 XM_022435637.1 LOC111102770  MEturquoise
#125 baculoviral IAP repeat-containing protein 7-A-like, transcript variant X1 XM_022432323.1 LOC111100432 MEdarkorange
#93               E3 ubiquitin-protein ligase XIAP-like, transcript variant X1 XM_022438073.1 LOC111104230  MEturquoise
#92               E3 ubiquitin-protein ligase XIAP-like, transcript variant X2 XM_022438074.1 LOC111104230  MEturquoise

# run comparison for GIMAP
Dermo_Tol_full_module_apop_df_GIMAP[!(Dermo_Tol_full_module_apop_df_GIMAP$gene %in% C_vir_full_all_exp_mod_sig_apop_GIMAP_minus_Dermo_Tol$gene), ]
#0
Dermo_Sus_full_module_apop_df_GIMAP[!(Dermo_Sus_full_module_apop_df_GIMAP$gene %in% C_vir_full_all_exp_mod_sig_apop_GIMAP_minus_Dermo_Sus$gene),]
#0
ROD_Res_full_module_apop_df_GIMAP[!(ROD_Res_full_module_apop_df_GIMAP$gene %in% C_vir_full_all_exp_mod_sig_apop_GIMAP_minus_ROD_Res$gene),]
#0
ROD_Sus_full_module_apop_df_GIMAP[!(ROD_Sus_full_module_apop_df_GIMAP$gene %in% C_vir_full_all_exp_mod_sig_apop_GIMAP_minus_Rod_Sus$gene),]
#60  GTPase IMAP family member 4-like, transcript variant X1 XM_022476897.1 LOC111130155 MEturquoise  0.9487006 0.05129943 ROD_Sus
#61  GTPase IMAP family member 4-like, transcript variant X2 XM_022476898.1 LOC111130155 MEturquoise  0.9487006 0.05129943 ROD_Sus
#94                         GTPase IMAP family member 4-like XM_022479991.1 LOC111132212 MEturquoise  0.9487006 0.05129943 ROD_Sus
#115 GTPase IMAP family member 4-like, transcript variant X2 XM_022439579.1 LOC111105339 MEturquoise  0.9487006 0.05129943 ROD_Sus
#127                        GTPase IMAP family member 4-like XM_022444820.1 LOC111108760 MEturquoise  0.9487006 0.05129943 ROD_Sus
#129                        GTPase IMAP family member 4-like XM_022440119.1 LOC111105744 MEturquoise  0.9487006 0.05129943 ROD_Sus
#141                        GTPase IMAP family member 7-like XM_022444543.1 LOC111108559 MEturquoise  0.9487006 0.05129943 ROD_Sus

Probiotic_full_module_apop_df_GIMAP[!(Probiotic_full_module_apop_df_GIMAP$gene %in% C_vir_full_all_exp_mod_sig_apop_GIMAP_minus_Probiotic$gene),]
#0
Pro_RE22_Pro_full_module_apop_df_GIMAP[!(Pro_RE22_Pro_full_module_apop_df_GIMAP$gene %in% C_vir_full_all_exp_mod_sig_apop_GIMAP_minus_Pro_RE22_S4$gene),]
#product  transcript_id         gene    mod_names mod_signif      pvalue             exp
#117 GTPase IMAP family member 4-like XM_022440606.1 LOC111106079 MElightpink4 -0.6601048 0.007404182 Pro_RE22_Pro_S4

Pro_RE22_Pro_RI_full_module_apop_df_GIMAP[!(Pro_RE22_Pro_RI_full_module_apop_df_GIMAP$gene %in% C_vir_full_all_exp_mod_sig_apop_GIMAP_minus_Pro_RE22_RI$gene),]
#product  transcript_id         gene mod_names mod_signif      pvalue             exp
#37 GTPase IMAP family member 4-like, transcript variant X2 XM_022476479.1 LOC111129930     MEred   0.708263 0.003126275 Pro_RE22_Pro_RI

Pro_RE22_RE22_full_module_apop_df_GIMAP[!(Pro_RE22_RE22_full_module_apop_df_GIMAP$gene %in% C_vir_full_all_exp_mod_sig_apop_GIMAP_minus_Pro_RE22_RE22$gene),]
#product  transcript_id         gene   mod_names Condition.Vibrio_coralliilyticus_RE22_exposure_6h.vs.Control_no_treatment      pvalue
#102 GTPase IMAP family member 4-like XM_022445423.1 LOC111109343 MEturquoise                                                                 0.9581529 0.002590135
#21  GTPase IMAP family member 8-like XM_022461041.1 LOC111120314       MEtan                                                                -0.9311959 0.006938151


### FIND UNIQUE IAP AND GIMAP ACROSS EXPERIMENT TYPE

### FIND C GIG IAP AND GIMAP SHARED BY ALL EXPERIMENTS 
## SIG IAPS
C_gig_full_all_exp_mod_sig_apop_IAP <- C_gig_full_all_exp_mod_sig_apop[grepl("IAP", C_gig_full_all_exp_mod_sig_apop$product, ignore.case=TRUE),]

# Find IAPS shared by all experiments
C_gig_full_all_exp_mod_sig_apop_IAP_split <- split(C_gig_full_all_exp_mod_sig_apop_IAP$gene, C_gig_full_all_exp_mod_sig_apop_IAP$exp)
Reduce(intersect, C_gig_full_all_exp_mod_sig_apop_IAP_split)
# 0 

## SIG GIMAPS 
C_gig_full_all_exp_mod_sig_apop_GIMAP <- C_gig_full_all_exp_mod_sig_apop[grepl("IMAP", C_gig_full_all_exp_mod_sig_apop$product, ignore.case=TRUE),]

# Find GIMAPS shared by all experiments
C_gig_full_all_exp_mod_sig_apop_GIMAP_split <- split(C_gig_full_all_exp_mod_sig_apop_GIMAP$gene, C_gig_full_all_exp_mod_sig_apop_GIMAP$exp)
Reduce(intersect, C_gig_full_all_exp_mod_sig_apop_GIMAP_split)
# 0 

### FIND SHARED IAP AND GIMAP ACROSS EXPERIMENT TYPE
# Using code from the section above this which has it joined by type 
C_gig_all_exp_mod_sig_apop_IAP_GENE_upset_split <- split(C_gig_all_exp_mod_sig_apop_IAP_GENE_upset$gene, C_gig_all_exp_mod_sig_apop_IAP_GENE_upset$type)
C_gig_IAP_gene_upset_list <- Reduce(intersect, C_gig_all_exp_mod_sig_apop_IAP_GENE_upset_split)
View(C_gig_all_exp_mod_sig_apop[C_gig_all_exp_mod_sig_apop$gene %in% C_gig_IAP_gene_upset_list,])

C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset_split <- split(C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset$gene, C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset$type)
C_gig_IMAP_gene_upset_list <- Reduce(intersect, C_gig_all_exp_mod_sig_apop_IMAP_GENE_upset_split)
View(C_gig_all_exp_mod_sig_apop[C_gig_all_exp_mod_sig_apop$gene %in% C_gig_IMAP_gene_upset_list,])


### FIND UNIQUE IAP AND GIMAP IN EACH EXPERIMENT 
# table of total per experiment
View(C_gig_full_all_exp_mod_sig_apop_IAP %>% group_by(exp) %>% distinct(gene) %>%  dplyr::summarise(n=n()))
View(C_gig_full_all_exp_mod_sig_apop_GIMAP %>% group_by(exp) %>% distinct(gene) %>%  dplyr::summarise(n=n()))

C_gig_full_all_exp_mod_sig_apop_IAP_minus_Zhang_LPS <-     C_gig_full_all_exp_mod_sig_apop_IAP %>% filter(exp != "Zhang_LPS")
C_gig_full_all_exp_mod_sig_apop_IAP_minus_Dermo_Sus <-     C_gig_full_all_exp_mod_sig_apop_IAP %>% filter(exp != "Zhang_Vibrio")
C_gig_full_all_exp_mod_sig_apop_IAP_minus_Rubio_NV <-       C_gig_full_all_exp_mod_sig_apop_IAP %>% filter(exp != "Rubio_NV")
C_gig_full_all_exp_mod_sig_apop_IAP_minus_Rubio_V <-       C_gig_full_all_exp_mod_sig_apop_IAP %>% filter(exp != "Rubio_V")
C_gig_full_all_exp_mod_sig_apop_IAP_minus_deLorg_Res <-     C_gig_full_all_exp_mod_sig_apop_IAP %>% filter(exp != "deLorg_Res")
C_gig_full_all_exp_mod_sig_apop_IAP_minus_deLorg_Sus <-   C_gig_full_all_exp_mod_sig_apop_IAP %>% filter(exp != "deLorg_Sus")
C_gig_full_all_exp_mod_sig_apop_IAP_minus_He <-      C_gig_full_all_exp_mod_sig_apop_IAP %>% filter(exp != "He")

C_gig_full_all_exp_mod_sig_apop_GIMAP_minus_Zhang_LPS <-    C_gig_full_all_exp_mod_sig_apop_GIMAP %>% filter(exp != "Zhang_LPS")
C_gig_full_all_exp_mod_sig_apop_GIMAP_minus_Dermo_Sus <-    C_gig_full_all_exp_mod_sig_apop_GIMAP %>% filter(exp != "Zhang_Vibrio")
C_gig_full_all_exp_mod_sig_apop_GIMAP_minus_Rubio_NV <-     C_gig_full_all_exp_mod_sig_apop_GIMAP %>% filter(exp != "Rubio_NV")
C_gig_full_all_exp_mod_sig_apop_GIMAP_minus_Rubio_V <-      C_gig_full_all_exp_mod_sig_apop_GIMAP %>% filter(exp != "Rubio_V")
C_gig_full_all_exp_mod_sig_apop_GIMAP_minus_deLorg_Res <-   C_gig_full_all_exp_mod_sig_apop_GIMAP %>% filter(exp != "deLorg_Res")
C_gig_full_all_exp_mod_sig_apop_GIMAP_minus_deLorg_Sus <-   C_gig_full_all_exp_mod_sig_apop_GIMAP %>% filter(exp != "deLorg_Sus")
C_gig_full_all_exp_mod_sig_apop_GIMAP_minus_He <-           C_gig_full_all_exp_mod_sig_apop_GIMAP %>% filter(exp != "He")

# select IAP and GIMAPs from apop dataframes
Zhang_LPS_full_module_apop_df_IAP <- Zhang_LPS_full_module_apop_df[grepl("IAP",Zhang_LPS_full_module_apop_df$product, ignore.case = TRUE),]
Zhang_full_module_Vibrio_apop_df_IAP <- Zhang_full_module_Vibrio_apop_df[grepl("IAP",Zhang_full_module_Vibrio_apop_df$product, ignore.case = TRUE),]
Rubio_NV_full_module_apop_df_IAP <- Rubio_NV_full_module_apop_df[grepl("IAP",Rubio_NV_full_module_apop_df$product, ignore.case = TRUE),]
Rubio_V_full_module_apop_df_IAP <- Rubio_V_full_module_apop_df[grepl("IAP",Rubio_V_full_module_apop_df$product, ignore.case = TRUE),]
deLorg_Res_full_module_apop_df_IAP <- deLorg_Res_full_module_apop_df[grepl("IAP",deLorg_Res_full_module_apop_df$product, ignore.case = TRUE),]
deLorg_Sus_full_module_apop_df_IAP <- deLorg_Sus_full_module_apop_df[grepl("IAP",deLorg_Sus_full_module_apop_df$product, ignore.case = TRUE),]
He_full_module_apop_df_IAP <- He_full_module_apop_df[grepl("IAP",He_full_module_apop_df$product, ignore.case = TRUE),]

Zhang_LPS_full_module_apop_df_GIMAP <- Zhang_LPS_full_module_apop_df[grepl("IMAP",Zhang_LPS_full_module_apop_df$product, ignore.case = TRUE),]
Zhang_full_module_Vibrio_apop_df_GIMAP <- Zhang_full_module_Vibrio_apop_df[grepl("IMAP",Zhang_full_module_Vibrio_apop_df$product, ignore.case = TRUE),]
Rubio_NV_full_module_apop_df_GIMAP <- Rubio_NV_full_module_apop_df[grepl("IMAP",Rubio_NV_full_module_apop_df$product, ignore.case = TRUE),]
Rubio_V_full_module_apop_df_GIMAP <- Rubio_V_full_module_apop_df[grepl("IMAP",Rubio_V_full_module_apop_df$product, ignore.case = TRUE),]
deLorg_Res_full_module_apop_df_GIMAP <- deLorg_Res_full_module_apop_df[grepl("IMAP",deLorg_Res_full_module_apop_df$product, ignore.case = TRUE),]
deLorg_Sus_full_module_apop_df_GIMAP <- deLorg_Sus_full_module_apop_df[grepl("IMAP",deLorg_Sus_full_module_apop_df$product, ignore.case = TRUE),]
He_full_module_apop_df_GIMAP <- He_full_module_apop_df[grepl("IMAP",He_full_module_apop_df$product, ignore.case = TRUE),]

Zhang_LPS_full_module_apop_df_IAP[!(Zhang_LPS_full_module_apop_df_IAP$gene %in%       C_gig_full_all_exp_mod_sig_apop_IAP_minus_Zhang_LPS$gene),]
#product  transcript_id         gene      mod_names mod_signif       pvalue       exp
#79          E3 ubiquitin-protein ligase XIAP-like XM_011454001.1 LOC105345723 MEmidnightblue  0.8319345 5.414170e-03 Zhang_LPS
#175   baculoviral IAP repeat-containing protein 7 XM_020064340.1 LOC105320825    MEturquoise  0.7070316 3.317257e-02 Zhang_LPS
#229 baculoviral IAP repeat-containing protein 7-B XM_020073460.1 LOC105343529 MEnavajowhite2 -0.9648457 2.591327e-05 Zhang_LPS
Zhang_full_module_Vibrio_apop_df_IAP[!(Zhang_full_module_Vibrio_apop_df_IAP$gene %in% C_gig_full_all_exp_mod_sig_apop_IAP_minus_Dermo_Sus$gene),]
#product  transcript_id         gene mod_names mod_signif     pvalue          exp
#26 baculoviral IAP repeat-containing protein 7-A, transcript variant X1 XM_011456705.2 LOC105347559   MEblue4 -0.8012361 0.00943345 Zhang_Vibrio

Rubio_NV_full_module_apop_df_IAP[!(Rubio_NV_full_module_apop_df_IAP$gene %in%         C_gig_full_all_exp_mod_sig_apop_IAP_minus_Rubio_NV $gene),]
#0
Rubio_V_full_module_apop_df_IAP[!(Rubio_V_full_module_apop_df_IAP$gene %in%           C_gig_full_all_exp_mod_sig_apop_IAP_minus_Rubio_V$gene),]
#0
deLorg_Res_full_module_apop_df_IAP[!(deLorg_Res_full_module_apop_df_IAP$gene %in%     C_gig_full_all_exp_mod_sig_apop_IAP_minus_deLorg_Res$gene),]
#0
deLorg_Sus_full_module_apop_df_IAP[!(deLorg_Sus_full_module_apop_df_IAP$gene %in%     C_gig_full_all_exp_mod_sig_apop_IAP_minus_deLorg_Sus$gene),]
#product  transcript_id         gene mod_names mod_signif pvalue        exp
#140 baculoviral IAP repeat-containing protein 5 XM_011437323.2 LOC105334034   ME15647         NA     NA deLorg_Sus
He_full_module_apop_df_IAP[!(He_full_module_apop_df_IAP$gene %in%                     C_gig_full_all_exp_mod_sig_apop_IAP_minus_He   $gene),]
#0

Zhang_LPS_full_module_apop_df_GIMAP[!(Zhang_LPS_full_module_apop_df_GIMAP$gene %in%       C_gig_full_all_exp_mod_sig_apop_GIMAP_minus_Zhang_LPS$gene),]
#product  transcript_id         gene     mod_names mod_signif     pvalue       exp
#4        GTPase IMAP family member 4 XM_011443127.2 LOC105338133 MEdarkorange2  0.6676620 0.04940214 Zhang_LPS
#46       GTPase IMAP family member 4 XM_020063057.1 LOC105317921       MEblack  0.6835217 0.04236273 Zhang_LPS
#143 GTPase IMAP family member 4-like XM_011439295.2 LOC105335444   MEturquoise  0.7070316 0.03317257 Zhang_LPS

Zhang_full_module_Vibrio_apop_df_GIMAP[!(Zhang_full_module_Vibrio_apop_df_GIMAP$gene %in% C_gig_full_all_exp_mod_sig_apop_GIMAP_minus_Dermo_Sus$gene),]
#0
Rubio_NV_full_module_apop_df_GIMAP[!(Rubio_NV_full_module_apop_df_GIMAP$gene %in%         C_gig_full_all_exp_mod_sig_apop_GIMAP_minus_Rubio_NV $gene),]
#0
Rubio_V_full_module_apop_df_GIMAP[!(Rubio_V_full_module_apop_df_GIMAP$gene %in%           C_gig_full_all_exp_mod_sig_apop_GIMAP_minus_Rubio_V$gene),]
#product  transcript_id         gene     mod_names mod_signif     pvalue     exp
#147 GTPase IMAP family member 4 XM_011435754.2 LOC105332988 MEgreenyellow -0.5532173 0.01724299 Rubio_V

deLorg_Res_full_module_apop_df_GIMAP[!(deLorg_Res_full_module_apop_df_GIMAP$gene %in%     C_gig_full_all_exp_mod_sig_apop_GIMAP_minus_deLorg_Res$gene),]
#product  transcript_id         gene mod_names mod_signif     pvalue        exp
#29 GTPase IMAP family member 8 XM_011416333.2 LOC105318994   MEblack -0.5390338 0.01168657 deLorg_Res

deLorg_Sus_full_module_apop_df_GIMAP[!(deLorg_Sus_full_module_apop_df_GIMAP$gene %in%     C_gig_full_all_exp_mod_sig_apop_GIMAP_minus_deLorg_Sus$gene),]
#0
He_full_module_apop_df_GIMAP[!(He_full_module_apop_df_GIMAP$gene %in%                     C_gig_full_all_exp_mod_sig_apop_GIMAP_minus_He   $gene),]
#product  transcript_id         gene mod_names mod_signif pvalue exp
#34       GTPase IMAP family member 8 XM_011422327.2 LOC105323259   ME43722         NA     NA  He
#92  GTPase IMAP family member 4-like XM_011424091.2 LOC105324868   ME45311         NA     NA  He
#144 GTPase IMAP family member 4-like XM_011429379.2 LOC105328481     ME809         NA     NA  He

### FIND UNIQUE IAP AND GIMAP ACROSS EXPERIMENT TYPE
# can do later 


#### ANNOTATION DATA: GENE LEVEL IAP GIMAP SUMMARY PLOT####

## Load gene family annotation stats from annotation workspace
load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gig_annotation_gene_family_info.RData")
load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_vir_annotation_gene_family_info.RData")

C_vir_rtracklayer_apop_product_final_product_joined_by_gene_duplicates_IAP <- C_vir_rtracklayer_apop_product_final_product_joined_by_gene_duplicates[grepl("IAP",C_vir_rtracklayer_apop_product_final_product_joined_by_gene_duplicates$gene_name, ignore.case=TRUE),]
C_vir_rtracklayer_apop_product_final_product_joined_by_gene_duplicates_GIMAP <-  C_vir_rtracklayer_apop_product_final_product_joined_by_gene_duplicates[grepl("IMAP",C_vir_rtracklayer_apop_product_final_product_joined_by_gene_duplicates$gene_name, ignore.case=TRUE),]
C_vir_genes_transcripts_per_gene_family_IAP_GIMAP <- C_vir_genes_transcripts_per_gene_family[grepl("IAP", C_vir_genes_transcripts_per_gene_family$apoptosis_names_query, ignore.case = TRUE) | 
                                                                                               grepl("IMAP",C_vir_genes_transcripts_per_gene_family$apoptosis_names_query, ignore.case = TRUE), ]
C_vir_genes_transcripts_per_gene_family_IAP_GIMAP$species <- "C_virginica"
C_vir_rtracklayer_apop_product_final_product_joined_by_gene_duplicates_IAP $species <- "C_virginica"
C_vir_rtracklayer_apop_product_final_product_joined_by_gene_duplicates_GIMAP$species <- "C_virginica"

C_gig_rtracklayer_apop_product_final_product_joined_by_gene_duplicates_IAP <-  C_gig_rtracklayer_apop_product_final_product_joined_by_gene_duplicates[grepl("IAP",C_gig_rtracklayer_apop_product_final_product_joined_by_gene_duplicates$gene_name, ignore.case=TRUE),]
C_gig_rtracklayer_apop_product_final_product_joined_by_gene_duplicates_GIMAP <-  C_gig_rtracklayer_apop_product_final_product_joined_by_gene_duplicates[grepl("IMAP",C_gig_rtracklayer_apop_product_final_product_joined_by_gene_duplicates$gene_name, ignore.case=TRUE),]
C_gig_genes_transcripts_per_gene_family_IAP_GIMAP <- C_gig_genes_transcripts_per_gene_family[grepl("IAP", C_gig_genes_transcripts_per_gene_family$apoptosis_names_query, ignore.case = TRUE) | 
                                                                                              grepl("IMAP",C_gig_genes_transcripts_per_gene_family$apoptosis_names_query, ignore.case = TRUE), ]
C_gig_rtracklayer_apop_product_final_product_joined_by_gene_duplicates_IAP $species <- "C_gigas"
C_gig_rtracklayer_apop_product_final_product_joined_by_gene_duplicates_GIMAP$species <- "C_gigas"

C_gig_genes_transcripts_per_gene_family_IAP_GIMAP$species <- "C_gigas"
Gene_family_members_combined <- full_join(C_vir_genes_transcripts_per_gene_family_IAP_GIMAP, C_gig_genes_transcripts_per_gene_family_IAP_GIMAP )
Gene_family_members_combined <- Gene_family_members_combined[order(Gene_family_members_combined$apoptosis_names_query),] 
View(Gene_family_members_combined[,c(1,2,4)])




# Tutorials: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/

#### SESSION INFO ####
sessionInfo()
#R version 3.6.1 (2019-07-05)
#Platform: x86_64-apple-darwin15.6.0 (64-bit)
#Running under: macOS Mojave 10.14
#
#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
#LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
#
#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#
#attached base packages:
#  [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
#
#other attached packages:
#  [1] DESeq2_1.24.0               SummarizedExperiment_1.14.1 DelayedArray_0.10.0         BiocParallel_1.18.1         matrixStats_0.54.0         
#[6] Biobase_2.44.0              GenomicRanges_1.36.0        GenomeInfoDb_1.20.0         IRanges_2.18.2              S4Vectors_0.22.0           
#[11] BiocGenerics_0.30.0        
#
#loaded via a namespace (and not attached):
#  [1] bitops_1.0-6           robust_0.4-18.1        fit.models_0.5-14      bit64_0.9-7            doParallel_1.0.15      RColorBrewer_1.1-2     dynamicTreeCut_1.63-1 
#[8] tools_3.6.1            backports_1.1.5        R6_2.4.1               rpart_4.1-15           Hmisc_4.2-0            DBI_1.0.0              colorspace_1.4-1      
#[15] nnet_7.3-12            tidyselect_0.2.5       gridExtra_2.3          bit_1.1-14             compiler_3.6.1         preprocessCore_1.46.0  WGCNA_1.68            
#[22] htmlTable_1.13.1       scales_1.1.0           checkmate_1.9.4        DEoptimR_1.0-8         mvtnorm_1.0-11         robustbase_0.93-5      genefilter_1.66.0     
#[29] stringr_1.4.0          digest_0.6.23          foreign_0.8-72         XVector_0.24.0         rrcov_1.4-7            base64enc_0.1-3        pkgconfig_2.0.3       
#[36] htmltools_0.3.6        htmlwidgets_1.3        rlang_0.4.2            rstudioapi_0.10        RSQLite_2.1.2          impute_1.58.0          acepack_1.4.1         
#[43] dplyr_0.8.3            RCurl_1.95-4.12        magrittr_1.5           GO.db_3.8.2            GenomeInfoDbData_1.2.1 Formula_1.2-3          Matrix_1.2-17         
#[50] Rcpp_1.0.3             munsell_0.5.0          lifecycle_0.1.0        yaml_2.2.0             stringi_1.4.6          MASS_7.3-51.4          zlibbioc_1.30.0       
#[57] grid_3.6.1             blob_1.2.0             crayon_1.3.4           lattice_0.20-38        splines_3.6.1          annotate_1.62.0        locfit_1.5-9.1        
#[64] zeallot_0.1.0          knitr_1.24             pillar_1.4.3           fastcluster_1.1.25     geneplotter_1.62.0     codetools_0.2-16       XML_3.99-0.3          
#[71] glue_1.3.1             latticeExtra_0.6-28    data.table_1.12.8      vctrs_0.2.1            foreach_1.4.7          gtable_0.3.0           purrr_0.3.3           
#[78] assertthat_0.2.1       ggplot2_3.3.0          xfun_0.9               xtable_1.8-4           survival_2.44-1.1      pcaPP_1.9-73           tibble_2.1.3          
#[85] iterators_1.0.12       AnnotationDbi_1.46.1   memoise_1.1.0          cluster_2.1.0        
