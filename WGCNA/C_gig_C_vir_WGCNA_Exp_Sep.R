# SCRIPT TO RUN SIGNED WGCNA ON PARASITE/PATHOGEN CHALLENGES 
# Erin Roberts, PhD Candidate University of Rhode Island 
# 4/9/2020

# This script perform WGCNA analysis on C. virginica and C. gigas experiments separately 

### LOAD PACKAGES ####
library(tidyverse)
library(limma)
library(WGCNA)
options(stringsAsFactors = FALSE) # run every time
library(cluster)
library(anRichment)
library(anRichmentMethods)
library(dplyr)
library(plyr)
library(magicfor)
cor <- WGCNA::cor # run every time

# source("https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/installAnRichment.R"); installAnRichment(); 

#### LOAD APOPTOSIS DATA FRAMES AND PACKAGES ####
Apoptosis_frames <- load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gig_C_vir_apoptosis_products.RData")
annotations <- load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_gig_C_vir_annotations.RData")
Apop_LFC <- load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/apop_LFC.RData")
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

length(colnames(matrix_common)[moduleColors == "grey60"]) # 405
length(colnames(matrix_common)[moduleColors == "lightgreen"]) # 395

# Use function to lookup all apop names for each significant module
matrix_common= Probiotic_dds_rlog_matrix_common
moduleColors= Probiotic_moduleColors
lookup =   C_vir_rtracklayer_apop_product_final

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
names(Pro_RE22_RE22_moduleTraitCor_Pval_df_sig_list_rm) <- c("turquoise", "cyan" ,     "red",  "darkolivegreen"  )
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
C_vir_all_exp_mod_sig_apop_positive <- C_vir_all_exp_mod_sig_apop %>% filter(mod_signif >0)

#### ONE WAY MEASURE MODULE PRESERVATION BETWEEN C. VIRGINICA DIFFERENT EXPERIMENTS ####

# Column names agree between the datasets, so we can go ahead and measure module preservation between all 
## 1. MEASURE OVERLAP BETWEEN RE22 AND RODs

## PRO_RE22 (all, not separated) VS. ROD RES ##
#Pro_RE22_ROD_multiExpr = list(Pro_RE22=list(data=Pro_RE22_dds_rlog_matrix_common),ROD_Res=list(data=ROD_Resistant_dds_rlog_matrix_common))
#Pro_RE22_multiColor = list(Pro_RE22 = Pro_RE22_moduleColors) 
#
##Load module preservation calculations
##Pro_RE22_ROD_Res_module_preservation = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_ROD_Res_module_preservation.RData")
#
##Pro_RE22_mp =modulePreservation(Pro_RE22_ROD_multiExpr, Pro_RE22_multiColor,
##                      referenceNetworks=1,
##                      verbose=3,
##                      networkType="signed hybrid", # use same signed hybrid as before
##                      corFnc = "bicor", # use recommended bicor as before 
##                      nPermutations=100,
##                      randomSeed = 1, # recommended in langfelder tutorial
##                      quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
#Pro_RE22_stats = Pro_RE22_mp$preservation$Z$ref.Pro_RE22$inColumnsAlsoPresentIn.ROD_Res
#Pro_RE22_stats_order <- Pro_RE22_stats[order(-Pro_RE22_stats[,2]),c(1:2)]
#Pro_RE22_stats_preserved <- Pro_RE22_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
#Pro_RE22_stats_MR <- Pro_RE22_mp$preservation$observed$ref.Pro_RE22$inColumnsAlsoPresentIn.ROD_Res
#Pro_RE22_stats_MR <- Pro_RE22_stats_MR [,c("moduleSize","medianRank.pres")]
#Pro_RE22_stats_MR_order <- Pro_RE22_stats_MR[order(Pro_RE22_stats_MR[,2]),]
### Moderate preservation between Pro_RE22 and ROD Res
## brown                 1000     7.8915583 (does okay on median rank)
## red                   1000     7.1531975 (not in top twenty rank of module preservation)
## steelblue              175     6.5939398 (steelblue is the top module in the ranking)
## black                 1000     6.3774745 (not in top twenty)
## midnightblue           527     5.7505155 (in top twenty)
## purple                 731     5.7203439 (is in the top twenty)
## blue                  1000     5.1007921 (not in top twenty)
## pink                   787     5.0031835 (pink not in top twenty)
#
## Save module preservation calculation
## save(Pro_RE22_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_ROD_Res_module_preservation.RData")
#
## Plot module preservation
## remember - low median rank means high preservation
#ggplot(Pro_RE22_stats_MR_order_less20 , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
#  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))
## mod_names moduleSize medianRank.pres
## 1    steelblue        175               1
## 2       coral1         37               5
## 3    honeydew1         38               5
## 4 navajowhite2         49               6
## 5       brown4         95               9
#
#Pro_RE22_stats_preserved_Zsum_medianRank <- merge(Pro_RE22_stats_MR_order_less20, Pro_RE22_stats_preserved)
#
## Were any of these preserved modules significant in Pro_RE22 or ROD_Res?
#Pro_RE22_RI_SIG_ROD_Res_preserved_all <- Pro_RE22_stats_preserved[Pro_RE22_stats_preserved$mod_name %in% Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm,]
#Pro_RE22_S4_SIG_ROD_Res_preserved_all <- Pro_RE22_stats_preserved[Pro_RE22_stats_preserved$mod_name %in% Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm,]
## midnightblue, blue, red 
#Pro_RE22_RE22_SIG_ROD_Res_preserved_all <- Pro_RE22_stats_preserved[Pro_RE22_stats_preserved$mod_name %in% Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm,]
##  black       1000      6.377475
## 2         blue       1000      5.100792
## 3 midnightblue        527      5.750516
## 4         pink        787      5.003184
## 5       purple        731      5.720344
## 6          red       1000      7.153197
## 7    steelblue        175      6.593940
#Pro_RE22_ROD_RES_SIG_preserved_all <- Pro_RE22_stats_preserved[Pro_RE22_stats_preserved$mod_name %in% ROD_Res_moduleTraitCor_Pval_df_sig_list_rm,]
## 0 
#
## Conclusions: 4 conserved between ROD_Res challenge and Pro_RE22, but are only significant with the RE22 and probiotic challenges. 
#  # No overlap in ones significant in both experiments with disease challenge 
#
### Pro_RE22_RE22 VS. ROD_Res ##
#Pro_RE22_RE22_ROD_multiExpr = list(Pro_RE22_RE22=list(data=Pro_RE22_dds_rlog_matrix_common_RE22),ROD_Res=list(data=ROD_Resistant_dds_rlog_matrix_common))
#Pro_RE22_RE22_multiColor = list(Pro_RE22_RE22 = Pro_RE22_RE22_moduleColors) 
## Pro_RE22_ROD_mp =modulePreservation(Pro_RE22_RE22_ROD_multiExpr, Pro_RE22_RE22_multiColor,
##                                 referenceNetworks=1,
##                                 verbose=3,
##                                 networkType="signed hybrid", # use same signed hybrid as before
##                                 corFnc = "bicor", # use recommended bicor as before 
##                                 nPermutations=100,
##                                 randomSeed = 1, # recommended in langfelder tutorial
##                                 quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
#Pro_RE22_RE22_ROD_Res_stats = Pro_RE22_ROD_mp$preservation$Z$ref.Pro_RE22_RE22$inColumnsAlsoPresentIn.ROD_Res
#Pro_RE22_RE22_ROD_Res_stats_order <- Pro_RE22_RE22_ROD_Res_stats[order(-Pro_RE22_RE22_ROD_Res_stats[,2]),c(1:2)]
#Pro_RE22_RE22_ROD_Res_stats_preserved <- Pro_RE22_RE22_ROD_Res_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
## green module is preserved
#Pro_RE22_RE22_ROD_Res_stats_MR <- Pro_RE22_ROD_mp$preservation$observed$ref.Pro_RE22_RE22$inColumnsAlsoPresentIn.ROD_Res
#Pro_RE22_RE22_ROD_Res_stats_MR <- Pro_RE22_RE22_ROD_Res_stats_MR [,c("moduleSize","medianRank.pres")]
#Pro_RE22_RE22_ROD_Res_stats_MR_order <-Pro_RE22_RE22_ROD_Res_stats_MR[order(-Pro_RE22_RE22_ROD_Res_stats_MR [,2]),]
#
## Save module preservation calculation
## save(Pro_RE22_ROD_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_RE22_ROD_Res_module_preservation.RData")
#
## Plot module preservation
## remember - low median rank means high preservation
#ggplot(Pro_RE22_RE22_ROD_Res_stats_MR_order, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
#  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))
#
## were any of the preserved modules also significant?
#Pro_RE22_RE22_only_SIG_ROD_Res_preserved <-Pro_RE22_RE22_ROD_Res_stats_preserved [Pro_RE22_RE22_ROD_Res_stats_preserved $mod_name %in% Pro_RE22_RE22_moduleTraitCor_Pval_df_sig_list_rm,]
## 0
#Pro_RE22_RE22_only_ROD_RES_SIG_preserved <- Pro_RE22_RE22_ROD_Res_stats_preserved [Pro_RE22_RE22_ROD_Res_stats_preserved $mod_name %in% ROD_Res_moduleTraitCor_Pval_df_sig_list_rm,]
## 0
#
#
#
###  ROD_Res VS. Pro_RE22_RE22 REVERSED ## RUNNING NOW need to code for comparing to network 2
#Pro_RE22_RE22_ROD_RES_multiExpr = list(Pro_RE22_RE22=list(data=Pro_RE22_dds_rlog_matrix_common_RE22),ROD_Res=list(data=ROD_Resistant_dds_rlog_matrix_common))
## Names for the two sets
#setLabels = c("Pro_RE22_RE22", "ROD_Res")
## Important: components of multiExpr must carry identificating names
#names(Pro_RE22_RE22_ROD_RES_multiExpr) = setLabels
## Create an object (list) holding the module labels for each set:
#colorList = list(Pro_RE22_RE22_moduleColors, ROD_Res_moduleColors)
## Components of the list must be named so that the names can be matched to the names of multiExpr
#names(colorList) = setLabels
#Pro_RE22_RE22_multiColor = list(Pro_RE22_RE22 = Pro_RE22_RE22_moduleColors) 
#Pro_RE22_ROD_Res_BOTH_mp = modulePreservation(Pro_RE22_RE22_ROD_RES_multiExpr, colorList ,
#                                              referenceNetworks=c(1:2),
#                                              verbose=3,
#                                              networkType="signed hybrid", # use same signed hybrid as before
#                                              corFnc = "bicor", # use recommended bicor as before 
#                                              nPermutations=100,
#                                              randomSeed = 1, # recommended in langfelder tutorial
#                                              quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
## save results
#save(Pro_RE22_ROD_Res_BOTH_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_ROD_Res_BOTH_mp_network.RData")
#
#Pro_RE22_RE22_ROD_Res_BOTHstats = Pro_RE22_ROD_mp$preservation$Z$ref.Pro_RE22_RE22$inColumnsAlsoPresentIn.ROD_Res
#Pro_RE22_RE22_ROD_Res_BOTHstats_order <- Pro_RE22_RE22_ROD_Res_stats[order(-Pro_RE22_RE22_ROD_Res_stats[,2]),c(1:2)]
#Pro_RE22_RE22_ROD_Res_BOTHstats_preserved <- Pro_RE22_RE22_ROD_Res_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
## green module is preserved
#Pro_RE22_RE22_ROD_Res_stats_MR <- Pro_RE22_ROD_mp$preservation$observed$ref.Pro_RE22_RE22$inColumnsAlsoPresentIn.ROD_Res
#Pro_RE22_RE22_ROD_Res_stats_MR <- Pro_RE22_RE22_ROD_Res_stats_MR [,c("moduleSize","medianRank.pres")]
#Pro_RE22_RE22_ROD_Res_stats_MR_order <-Pro_RE22_RE22_ROD_Res_stats_MR[order(-Pro_RE22_RE22_ROD_Res_stats_MR [,2]),]
#
## Save module preservation calculation
## save(Pro_RE22_ROD_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_RE22_ROD_Res_module_preservation.RData")
#
## Plot module preservation
## remember - low median rank means high preservation
#ggplot(Pro_RE22_RE22_ROD_Res_stats_MR_order, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
#  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))
#
## were any of the preserved modules also significant?
#Pro_RE22_RE22_only_SIG_ROD_Res_preserved <-Pro_RE22_RE22_ROD_Res_stats_preserved [Pro_RE22_RE22_ROD_Res_stats_preserved $mod_name %in% Pro_RE22_RE22_moduleTraitCor_Pval_df_sig_list_rm,]
## 0
#Pro_RE22_RE22_only_ROD_RES_SIG_preserved <- Pro_RE22_RE22_ROD_Res_stats_preserved [Pro_RE22_RE22_ROD_Res_stats_preserved $mod_name %in% ROD_Res_moduleTraitCor_Pval_df_sig_list_rm,]
## 0
#
### Pro_RE22 (all not separated) vs. ROD_Sus modules ##
# Pro_RE22_Sus_multiExpr = list(Pro_RE22=list(data=Pro_RE22_dds_rlog_matrix_common),ROD_Sus=list(data=ROD_Susceptible_dds_rlog_matrix_common))
#
##Load module preservation calculations
##Pro_RE22_ROD_Sus_module_preservation = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_ROD_Sus_module_preservation.RData")
#
## Pro_RE22_Sus_mp =modulePreservation(Pro_RE22_Sus_multiExpr, Pro_RE22_multiColor,
##                                 referenceNetworks=1,
##                                 verbose=3,
##                                 networkType="signed hybrid", # use same signed hybrid as before
##                                 corFnc = "bicor", # use recommended bicor as before 
##                                 nPermutations=100,
##                                 randomSeed = 1, # recommended in langfelder tutorial
##                                 quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
#Pro_RE22_Sus_stats = Pro_RE22_Sus_mp$preservation$Z$ref.Pro_RE22$inColumnsAlsoPresentIn.ROD_Sus 
#Pro_RE22_Sus_stats_order <- Pro_RE22_Sus_stats[order(-Pro_RE22_Sus_stats[,2]),c(1:2)]
#Pro_RE22_Sus_stats_preserved <- Pro_RE22_Sus_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
##mod_name moduleSize Zsummary.pres
##1      red       1000      7.327336 (good median rank)
##2     blue       1000      5.772638 (good median rank)
##3    brown       1000      5.070500 (bad median rank)
#
#Pro_RE22_Sus_stats_MR <- Pro_RE22_Sus_mp$preservation$observed$ref.Pro_RE22$inColumnsAlsoPresentIn.ROD_Sus
#Pro_RE22_Sus_stats_MR <- Pro_RE22_Sus_stats_MR [,c("moduleSize","medianRank.pres")]
#Pro_RE22_Sus_stats_MR_order <- Pro_RE22_Sus_stats_MR[order(Pro_RE22_Sus_stats_MR[,2]),]
#Pro_RE22_Sus_stats_MR_order_less20 <- Pro_RE22_Sus_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)
#
## Save module preservation calculation
##save(Pro_RE22_Sus_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_ROD_Sus_module_preservation.RData")
#
## Plot module preservation
## remember - low median rank means high preservation
#ggplot(Pro_RE22_Sus_stats_MR_order_less20 , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
# geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))
#
#Pro_RE22_Sus_stats_preserved_Zsum_medianRank <- merge(Pro_RE22_Sus_stats_MR_order_less20, Pro_RE22_Sus_stats_preserved)
## blue, red (brown )
#
#Pro_RE22_RI_SIG_ROD_Sus_preserved_all <- Pro_RE22_Sus_stats_preserved[Pro_RE22_Sus_stats_preserved$mod_name %in% Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm,]
## 1      red       1000      7.327336
## 2     blue       1000      5.772638
#Pro_RE22_S4_SIG_ROD_Sus_preserved_all <- Pro_RE22_Sus_stats_preserved[Pro_RE22_Sus_stats_preserved$mod_name %in%  Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm,]
## 1      red       1000      7.327336
## 2     blue       1000      5.772638
#Pro_RE22_RE22_SIG_ROD_Sus_preserved_all <- Pro_RE22_Sus_stats_preserved[Pro_RE22_Sus_stats_preserved$mod_name %in%  Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm,]
## 1      red       1000      7.327336
## 2     blue       1000      5.772638
#Pro_RE22_ROD_SUS_SIG_preserved_all <- Pro_RE22_Sus_stats_preserved[Pro_RE22_Sus_stats_preserved$mod_name %in%  ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
## brown

### Assess whether Pro_RE22_RE22 vs. ROD_Sus modules are preserved ##
#Pro_RE22_RE22_ROD_Sus_multiExpr = list(Pro_RE22_RE22=list(data=Pro_RE22_dds_rlog_matrix_common_RE22),ROD_Sus=list(data=ROD_Susceptible_dds_rlog_matrix_common))
#Pro_RE22_RE22_multiColor = list(Pro_RE22_RE22 = Pro_RE22_RE22_moduleColors) 
## Pro_RE22_ROD_Sus_mp =modulePreservation(Pro_RE22_RE22_ROD_Sus_multiExpr, Pro_RE22_RE22_multiColor,
##                                    referenceNetworks=1,
##                                    verbose=3,
##                                    networkType="signed hybrid", # use same signed hybrid as before
##                                    corFnc = "bicor", # use recommended bicor as before 
##                                    nPermutations=100,
##                                    randomSeed = 1, # recommended in langfelder tutorial
##                                    quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
#Pro_RE22_ROD_Sus_stats = Pro_RE22_ROD_Sus_mp$preservation$Z$ref.Pro_RE22_RE22$inColumnsAlsoPresentIn.ROD_Sus
#Pro_RE22_ROD_Sus_stats_order <- Pro_RE22_ROD_Sus_stats[order(-Pro_RE22_ROD_Sus_stats[,2]),c(1:2)]
#Pro_RE22_ROD_Sus_stats_preserved <- Pro_RE22_ROD_Sus_stats_order  %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
## 0 modules preserved between Pro_RE22_RE22 and the ROD Susceptible 
#Pro_RE22_ROD_Sus_stats_MR <- Pro_RE22_ROD_Sus_mp$preservation$observed$ref.Pro_RE22_RE22$inColumnsAlsoPresentIn.ROD_Sus
#Pro_RE22_ROD_Sus_stats_MR <- Pro_RE22_ROD_Sus_stats_MR  [,c("moduleSize","medianRank.pres")]
#Pro_RE22_ROD_Sus_stats_MR_order <- Pro_RE22_ROD_Sus_stats_MR [order(Pro_RE22_ROD_Sus_stats_MR [,2]),]
#
## Save network
## save(Pro_RE22_ROD_Sus_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_RE22_ROD_Sus_module_preservation.RData")
#
## Plot module preservation
## remember - low median rank means high preservation
#ggplot(Pro_RE22_ROD_Sus_stats_MR_order, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
#  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))
#
## were any of the preserved modules also significant?
#Pro_RE22_RE22_only_SIG_ROD_Sus_preserved <- Pro_RE22_ROD_Sus_stats_preserved[Pro_RE22_ROD_Sus_stats_preserved$mod_name %in% Pro_RE22_RE22_moduleTraitCor_Pval_df_sig_list_rm,]
##0 
#Pro_RE22_RE22_only_ROD_Sus_SIG_preserved <- Pro_RE22_ROD_Sus_stats_preserved[Pro_RE22_ROD_Sus_stats_preserved$mod_name %in% ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
##0
#
#Pro_RE22_RE22_ROD_Sus_multiExpr = list(Pro_RE22_RE22=list(data=Pro_RE22_dds_rlog_matrix_common_RE22),ROD_Sus=list(data=ROD_Susceptible_dds_rlog_matrix_common))
#Pro_RE22_RE22_multiColor = list(Pro_RE22_RE22 = Pro_RE22_RE22_moduleColors) 
## Pro_RE22_ROD_Sus_mp =modulePreservation(Pro_RE22_RE22_ROD_Sus_multiExpr, Pro_RE22_RE22_multiColor,
##                                    referenceNetworks=1,
##                                    verbose=3,
##                                    networkType="signed hybrid", # use same signed hybrid as before
##                                    corFnc = "bicor", # use recommended bicor as before 
##                                    nPermutations=100,
##                                    randomSeed = 1, # recommended in langfelder tutorial
##                                    quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
#Pro_RE22_ROD_Sus_stats = Pro_RE22_ROD_Sus_mp$preservation$Z$ref.Pro_RE22_RE22$inColumnsAlsoPresentIn.ROD_Sus
#Pro_RE22_ROD_Sus_stats_order <- Pro_RE22_ROD_Sus_stats[order(-Pro_RE22_ROD_Sus_stats[,2]),c(1:2)]
#Pro_RE22_ROD_Sus_stats_preserved <- Pro_RE22_ROD_Sus_stats_order  %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
## 0 modules preserved between Pro_RE22_RE22 and the ROD Susceptible 
#Pro_RE22_ROD_Sus_stats_MR <- Pro_RE22_ROD_Sus_mp$preservation$observed$ref.Pro_RE22_RE22$inColumnsAlsoPresentIn.ROD_Sus
#Pro_RE22_ROD_Sus_stats_MR <- Pro_RE22_ROD_Sus_stats_MR  [,c("moduleSize","medianRank.pres")]
#Pro_RE22_ROD_Sus_stats_MR_order <- Pro_RE22_ROD_Sus_stats_MR [order(Pro_RE22_ROD_Sus_stats_MR [,2]),]
#
## Save network
## save(Pro_RE22_ROD_Sus_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_RE22_ROD_Sus_module_preservation.RData")
#
## Plot module preservation
## remember - low median rank means high preservation
#ggplot(Pro_RE22_ROD_Sus_stats_MR_order, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
#  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))
#
## were any of the preserved modules also significant?
#Pro_RE22_RE22_only_SIG_ROD_Sus_preserved <- Pro_RE22_ROD_Sus_stats_preserved[Pro_RE22_ROD_Sus_stats_preserved$mod_name %in% Pro_RE22_RE22_moduleTraitCor_Pval_df_sig_list_rm,]
##0 
#Pro_RE22_RE22_only_ROD_Sus_SIG_preserved <- Pro_RE22_ROD_Sus_stats_preserved[Pro_RE22_ROD_Sus_stats_preserved$mod_name %in% ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
##0

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

### Assess whether Pro_RE22 (all not separated) and Dermo_Tol modules are preserved ##
#Pro_RE22_Dermo_multiExpr = list(Pro_RE22=list(data=Pro_RE22_dds_rlog_matrix_common),Dermo_Tol=list(data=Dermo_Tolerant_dds_vst_matrix_common))
#Pro_RE22_multiColor = list(Pro_RE22 = Pro_RE22_moduleColors) 
## Load module preservation calculations
##Pro_RE22_Dermo_Tol = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_Dermo_Tol_module_preservation.RData")
## Pro_RE22_Dermo_mp =modulePreservation(Pro_RE22_Dermo_multiExpr, Pro_RE22_multiColor,
##                                 referenceNetworks=1,
##                                 verbose=3,
##                                 networkType="signed hybrid", # use same signed hybrid as before
##                                 corFnc = "bicor", # use recommended bicor as before 
##                                 nPermutations=100,
##                                 randomSeed = 1, # recommended in langfelder tutorial
##                                 quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
#Pro_RE22_Dermo_stats = Pro_RE22_Dermo_mp$preservation$Z$ref.Pro_RE22$inColumnsAlsoPresentIn.Dermo_Tol
#Pro_RE22_Dermo_stats_order <- Pro_RE22_Dermo_stats[order(-Pro_RE22_Dermo_stats[,2]),c(1:2)]
#Pro_RE22_Dermo_stats_preserved <- Pro_RE22_Dermo_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
## mod_name moduleSize Zsummary.pres
## 1        salmon        615     10.598431 (high median rank)
## 2        purple        731      9.399914 (decent median rank)
## 3     steelblue        175      9.376713 (high median rank)
## 4         brown       1000      9.373572 (decent median rank)
## 5           red       1000      8.038084 (lower median rank)
## 6     turquoise       1000      7.472569 (lower median rank)
## 7        yellow       1000      7.446067 (lower median rank)
## 8          blue       1000      6.243046 (lower median rank)
## 9         green       1000      5.794057 (lower median rank)
## 10    lightcyan        474      5.753286 (lower median rank)
## 11 midnightblue        527      5.580933 (lower median rank)
## 12        black       1000      5.354414 (lowest median rank in top twentry)
#
#Pro_RE22_Dermo_stats_MR <- Pro_RE22_Dermo_mp$preservation$observed$ref.Pro_RE22$inColumnsAlsoPresentIn.Dermo_Tol
#Pro_RE22_Dermo_stats_MR <- Pro_RE22_Dermo_stats_MR [,c("moduleSize","medianRank.pres")]
#Pro_RE22_Dermo_stats_MR_order <- Pro_RE22_Dermo_stats_MR[order(Pro_RE22_Dermo_stats_MR[,2]),]
#Pro_RE22_Dermo_stats_MR_order_less20 <- Pro_RE22_Dermo_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)
#
## Save network
##save(Pro_RE22_Dermo_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_Dermo_Tol_module_preservation.RData")
#
## Plot module preservation
## remember - low median rank means high preservation
#ggplot(Pro_RE22_Dermo_stats_MR_order_less20 , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
# geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))
#
#Pro_RE22_Dermo_stats_preserved_Zsum_medianRank <- merge(Pro_RE22_Dermo_stats_MR_order_less20, Pro_RE22_Dermo_stats_preserved)
## 12 modules
#
#Pro_RE22_RI_SIG_Dermo_Tol_preserved_all <- Pro_RE22_Dermo_stats_preserved[Pro_RE22_Dermo_stats_preserved$mod_names %in% Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm,]
##1        black       1000              19      5.354414
##2         blue       1000              13      6.243046
##3        green       1000              18      5.794057
##4 midnightblue        527              13      5.580933
##5       purple        731              11      9.399914
##6          red       1000              15      8.038084
##7       salmon        615               4     10.598431
##8    steelblue        175               1      9.376713
##9       yellow       1000              14      7.446067
#Pro_RE22_S4_SIG_Dermo_Tol_preserved_all <- Pro_RE22_Dermo_stats_preserved[Pro_RE22_Dermo_stats_preserved$mod_names %in% Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm,]
##mod_name moduleSize medianRank.pres Zsummary.pres
##1         blue       1000              13      6.243046
##2        green       1000              18      5.794057
##3 midnightblue        527              13      5.580933
##4          red       1000              15      8.038084
#Pro_RE22_RE22_SIG_RDermo_Tol_preserved_all <- Pro_RE22_Dermo_stats_preserved[Pro_RE22_Dermo_stats_preserved$mod_names %in% Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm,]
##       mod_name moduleSize medianRank.pres Zsummary.pres
## 1        black       1000              19      5.354414
## 2         blue       1000              13      6.243046
## 3        green       1000              18      5.794057
## 4 midnightblue        527              13      5.580933
## 5       purple        731              11      9.399914
## 6          red       1000              15      8.038084
## 7       salmon        615               4     10.598431
## 8    steelblue        175               1      9.376713
## 9       yellow       1000              14      7.446067
#Pro_RE22_Dermo_Tol_SIG_preserved_all <- Pro_RE22_Dermo_stats_preserved[Pro_RE22_Dermo_stats_preserved$mod_names %in% Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm,]
## mod_name moduleSize medianRank.pres Zsummary.pres
## 1 turquoise       1000              14      7.472569
## 2    yellow       1000              14      7.446067
#
## Conclusion: probiotic and Dermo tolerant share two modules significantly correlated with disease challenge 
#
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

# Conclusion: 

## Assess whether Pro_RE22 (all, not separated) and Dermo_Sus modules are preserved ##
#Pro_RE22_Dermo_Sus_multiExpr = list(Pro_RE22=list(data=Pro_RE22_dds_rlog_matrix_common),Dermo_Sus=list(data=Dermo_Susceptible_dds_vst_matrix_common))
#Pro_RE22_multiColor = list(Pro_RE22 = Pro_RE22_moduleColors) 
## Load module preservation calculations
##Pro_RE22_Dermo_Sus = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_Dermo_Sus_module_preservation.RData")
## Pro_RE22_Dermo_Sus_mp =modulePreservation(Pro_RE22_Dermo_Sus_multiExpr, Pro_RE22_multiColor,
##                                       referenceNetworks=1,
##                                       verbose=3,
##                                       networkType="signed hybrid", # use same signed hybrid as before
##                                       corFnc = "bicor", # use recommended bicor as before 
##                                       nPermutations=100,
##                                       randomSeed = 1, # recommended in langfelder tutorial
##                                       quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
#Pro_RE22_Dermo_Sus_stats = Pro_RE22_Dermo_Sus_mp$preservation$Z$ref.Pro_RE22$inColumnsAlsoPresentIn.Dermo_Sus
#Pro_RE22_Dermo_Sus_stats_order <- Pro_RE22_Dermo_Sus_stats[order(-Pro_RE22_Dermo_Sus_stats[,2]),c(1:2)]
#Pro_RE22_Dermo_Sus_stats_preserved <- Pro_RE22_Dermo_Sus_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
## mod_name moduleSize Zsummary.pres
## 1    salmon        615      9.461935 (high median rank)
## 2 steelblue        175      7.608015 (highest median rank)
## 3    purple        731      6.821283 (good median rank)
## 4      blue       1000      6.176316 (medium median rank)
## 5       red       1000      6.060967 (medium median rank)
## 6 lightcyan        474      6.010655 (lower median rank)
## 7     brown       1000      5.165484 (not in top twentry)
#Pro_RE22_Dermo_Sus_stats_MR <- Pro_RE22_Dermo_Sus_mp$preservation$observed$ref.Pro_RE22$inColumnsAlsoPresentIn.Dermo_Sus
#Pro_RE22_Dermo_Sus_stats_MR <- Pro_RE22_Dermo_Sus_stats_MR [,c("moduleSize","medianRank.pres")]
#Pro_RE22_Dermo_Sus_stats_MR_order <- Pro_RE22_Dermo_Sus_stats_MR[order(Pro_RE22_Dermo_Sus_stats_MR[,2]),]
#Pro_RE22_Dermo_Sus_stats_MR_order_less20 <- Pro_RE22_Dermo_Sus_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)
#
## Save network
##save(Pro_RE22_Dermo_Sus_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_Dermo_Sus_module_preservation.RData")
#
## Plot module preservation
## remember - low median rank means high preservation
#ggplot(Pro_RE22_Dermo_Sus_stats_MR_order_less20 , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
#  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))
#
#Pro_RE22_Dermo_Sus_stats_preserved_Zsum_medianRank <- Pro_RE22_Dermo_Sus_stats_MR_order_less20[Pro_RE22_Dermo_Sus_stats_MR_order_less20$mod_names %in% Pro_RE22_Dermo_Sus_stats_preserved,]
## 6 modules
#
## Were any of these preserved modules significant in Pro_RE22 or Dermo_Sus
#
#Pro_RE22_RI_SIG_Dermo_Sus_preserved_all <- Pro_RE22_Dermo_Sus_stats_preserved[Pro_RE22_Dermo_Sus_stats_preserved$mod_name %in% Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm,]
## mod_name moduleSize Zsummary.pres
## 1      blue       1000      6.176316
## 2    purple        731      6.821283
## 3       red       1000      6.060967
## 4    salmon        615      9.461935
## 5 steelblue        175      7.608015
#Pro_RE22_S4_SIG_Dermo_Sus_preserved_all <-Pro_RE22_Dermo_Sus_stats_preserved[Pro_RE22_Dermo_Sus_stats_preserved$mod_name %in% Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm,]
## blue, red
#Pro_RE22_RE22_SIG_Dermo_Sus_preserved_all <-Pro_RE22_Dermo_Sus_stats_preserved[Pro_RE22_Dermo_Sus_stats_preserved$mod_name %in% Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm,]
## mod_name moduleSize Zsummary.pres
## 1      blue       1000      6.176316
## 2    purple        731      6.821283
## 3       red       1000      6.060967
## 4    salmon        615      9.461935
## 5 steelblue        175      7.608015
#Pro_RE22_Dermo_Sus_SIG_preserved_all <-Pro_RE22_Dermo_Sus_stats_preserved[Pro_RE22_Dermo_Sus_stats_preserved$mod_name %in% Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
## 0 

## Assess whether Pro_RE22_RE22 and Dermo_Sus modules are preserved ##
#Pro_RE22_RE22_Dermo_Sus_multiExpr = list(Pro_RE22_RE22=list(data=Pro_RE22_dds_rlog_matrix_common_RE22),Dermo_Sus=list(data=Dermo_Susceptible_dds_vst_matrix_common))
#Pro_RE22_RE22_multiColor = list(Pro_RE22_RE22 = Pro_RE22_RE22_moduleColors) 
##Pro_RE22_RE22_Dermo_Sus_mp =modulePreservation(Pro_RE22_RE22_Dermo_Sus_multiExpr , Pro_RE22_RE22_multiColor,
##                                          referenceNetworks=1,
##                                          verbose=3,
##                                          networkType="signed hybrid", # use same signed hybrid as before
##                                          corFnc = "bicor", # use recommended bicor as before 
##                                          nPermutations=100,
##                                          randomSeed = 1, # recommended in langfelder tutorial
##                                          quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
#Pro_RE22_RE22_Dermo_Sus_stats = Pro_RE22_RE22_Dermo_Sus_mp$preservation$Z$ref.Pro_RE22_RE22$inColumnsAlsoPresentIn.Dermo_Sus
#Pro_RE22_RE22_Dermo_Sus_stats_order <- Pro_RE22_RE22_Dermo_Sus_stats[order(-Pro_RE22_RE22_Dermo_Sus_stats[,2]),c(1:2)]
#Pro_RE22_RE22_Dermo_Sus_stats_preserved <- Pro_RE22_RE22_Dermo_Sus_stats_order  %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
##mod_name moduleSize Zsummary.pres
##1    black       1000      5.182457
#Pro_RE22_RE22_Dermo_Sus_stats_MR <- Pro_RE22_RE22_Dermo_Sus_mp$preservation$observed$ref.Pro_RE22_RE22$inColumnsAlsoPresentIn.Dermo_Sus
#Pro_RE22_RE22_Dermo_Sus_stats_MR <- Pro_RE22_RE22_Dermo_Sus_stats_MR [,c("moduleSize","medianRank.pres")]
#Pro_RE22_RE22_Dermo_Sus_stats_MR_order <- Pro_RE22_RE22_Dermo_Sus_stats_MR[order(-Pro_RE22_RE22_Dermo_Sus_stats_MR[,2]),]
#Pro_RE22_RE22_Dermo_Sus_stats_MR_order_less20 <- Pro_RE22_RE22_Dermo_Sus_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)
#
## Save network
##save(Pro_RE22_RE22_Dermo_Sus_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Pro_RE22_RE22_Dermo_Sus_module_preservation.RData")
#
## Plot module preservation
## remember - low median rank means high preservation
#ggplot(Pro_RE22_RE22_Dermo_Sus_stats_MR_order_less20, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
#  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))
#
## Were any of these preserved modules significant in Pro_RE22 or Dermo_Sus
#Pro_RE22_RE22_SIG_Dermo_Sus_preserved_all <-Pro_RE22_RE22_Dermo_Sus_stats_preserved[Pro_RE22_RE22_Dermo_Sus_stats_preserved$mod_name %in% Pro_RE22_RE22_moduleTraitCor_Pval_df_sig_list_rm,]
## 0
#Pro_RE22_RE22_Dermo_Sus_SIG_preserved_all <-Pro_RE22_RE22_Dermo_Sus_stats_preserved[Pro_RE22_RE22_Dermo_Sus_stats_preserved$mod_name %in% Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
##0
#
## Dermo_Sus vs.  Pro_RE22_RE22 ##
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
#ROD_Dermo_multiExpr = list(Dermo_Tol=list(data=Dermo_Tolerant_dds_vst_matrix_common), ROD_Res=list(data=ROD_Resistant_dds_rlog_matrix_common))
#ROD_Dermo_multiColor = list(Dermo_Tol = Dermo_Tol_moduleColors)
## Load module preservation
##Dermo_Tol_ROD_Res = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Tol_ROD_Res_module_preservation.RData")
#
##ROD_Dermo_mp =modulePreservation(ROD_Dermo_multiExpr, ROD_Dermo_multiColor,
##                                 referenceNetworks=1,
##                                 verbose=3,
##                                 networkType="signed hybrid", # use same signed hybrid as before
##                                 corFnc = "bicor", # use recommended bicor as before 
##                                 nPermutations=100,
##                                 randomSeed = 1, # recommended in langfelder tutorial
##                                 quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
#ROD_Dermo_stats = ROD_Dermo_mp$preservation$Z$ref.Dermo_Tol$inColumnsAlsoPresentIn.ROD_Res
#ROD_Dermo_stats_order <- ROD_Dermo_stats[order(-ROD_Dermo_stats[,2]),c(1:2)]
#ROD_Dermo_stats_preserved <- ROD_Dermo_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
## mod_name moduleSize Zsummary.pres
## 1        pink        552     21.124021 (high median rank)
## 2     magenta        458     14.461520 (medium median rank)
## 3         tan        300     14.281723 (high median rank)
## 4       green        762     14.061856 (lower median rank)
## 5   turquoise       1000     13.976882 (lower median rank)
## 6   darkgreen        218     11.000436 (high median rank)
## 7        blue       1000      9.858004 (not in top 20)
## 8        cyan        261      8.939884 (high median rank)
## 9         red        628      8.528091 (not in top 20)
## 10      black        618      7.999981 (not in top 20)
## 11      plum2        127      6.854030 (good rank)
## 12 lightcyan1        146      6.714021 (medium rank)
## 13     grey60        246      6.669582 (not in top 20)
## 14     violet        174      6.052549 (medium rank)
## 15      white        187      5.852156 (low medium rank)
## 16      brown       1000      5.646930 (not in top twenty)
## 17 darkviolet         70      5.146672 (higest median rank)
#
#ROD_Dermo_stats_MR <- ROD_Dermo_mp$preservation$observed$ref.Dermo_Tol$inColumnsAlsoPresentIn.ROD_Res
#ROD_Dermo_stats_MR <-ROD_Dermo_stats_MR[,c("moduleSize","medianRank.pres")]
#ROD_Dermo_stats_MR_order <- ROD_Dermo_stats_MR[order(ROD_Dermo_stats_MR[,2]),]
#ROD_Dermo_stats_MR_order_less20 <- ROD_Dermo_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)
#
## Save network
##save(ROD_Dermo_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Tol_ROD_Res_module_preservation.RData")
#
## Plot module preservation
## remember - low median rank means high preservation
#ggplot(ROD_Dermo_stats_MR_order_less20  , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
#  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))
#
#ROR_Dermo_stats_preserved_Zsum_medianRank <- ROD_Dermo_stats_MR_order_less20[ROD_Dermo_stats_MR_order_less20$mod_name , ROD_Dermo_stats_preserved,]
## 6 modules
#
## Were any of these preserved modules significant in Pro_RE22 or Dermo_Tol?

#Dermo_Tol_ROD_Res <- ROD_Dermo_stats_preserved[ROD_Dermo_stats_preserved$mod_name %in% Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm,]
## mod_name moduleSize Zsummary.pres
## 1        pink        552     21.124021
## 3         tan        300     14.281723
## 5   turquoise       1000     13.976882
## 11      plum2        127      6.854030
## 12 lightcyan1        146      6.714021
## 17 darkviolet         70      5.146672
#Dermo_Tol_ROD_Res_ROD <- ROD_Dermo_stats_preserved[ROD_Dermo_stats_preserved$mod_name %in% ROD_Res_moduleTraitCor_Pval_df_sig_list_rm,]
##0

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
#ROD_Dermo_Sus_multiExpr = list(Dermo_Sus=list(data=Dermo_Susceptible_dds_vst_matrix_common), ROD_Sus=list(data=ROD_Susceptible_dds_rlog_matrix_common))
#ROD_Dermo_Sus_multiColor = list(Dermo_Sus = Dermo_Sus_moduleColors)
## Load module preservation
##Dermo_Sus_ROD_Sus = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Sus_ROD_Sus_module_preservation.RData")
#
##ROD_Dermo_Sus_mp =modulePreservation(ROD_Dermo_Sus_multiExpr, ROD_Dermo_Sus_multiColor,
##                                     referenceNetworks=1,
##                                     verbose=3,
##                                     networkType="signed hybrid", # use same signed hybrid as before
##                                     corFnc = "bicor", # use recommended bicor as before 
##                                     nPermutations=100,
##                                     randomSeed = 1, # recommended in langfelder tutorial
##                                     quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
#ROD_Dermo_Sus_stats = ROD_Dermo_Sus_mp$preservation$Z$ref.Dermo_Sus$inColumnsAlsoPresentIn.ROD_Sus 
#ROD_Dermo_Sus_stats_order <- ROD_Dermo_Sus_stats[order(-ROD_Dermo_Sus_stats[,2]),c(1:2)] 
#ROD_Dermo_Sus_stats_preserved <- ROD_Dermo_Sus_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
## mod_name moduleSize Zsummary.pres
## 1       black        558     14.052762
## 2   turquoise       1000     10.807381 (only one in the top 20 median rank)
## 3 lightyellow        253      8.239593
## 4      purple        446      6.069997
## 5      salmon        402      5.345025
#
#ROD_Dermo_Sus_stats_MR <- ROD_Dermo_Sus_mp$preservation$observed$ref.Dermo_Sus$inColumnsAlsoPresentIn.ROD_Sus
#ROD_Dermo_Sus_stats_MR <-ROD_Dermo_stats_MR[,c("moduleSize","medianRank.pres")]
#ROD_Dermo_Sus_stats_MR_order <- ROD_Dermo_Sus_stats_MR[order(ROD_Dermo_stats_MR[,2]),]
#ROD_Dermo_Sus_stats_MR_order_less20 <- ROD_Dermo_Sus_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)
#
## Save network
##save(ROD_Dermo_Sus_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Dermo_Sus_ROD_Sus_module_preservation.RData")
#
## Plot module preservation
## remember - low median rank means high preservation
#ggplot(ROD_Dermo_Sus_stats_MR_order_less20  , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
#  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))
#
#ROD_Dermo_Sus_stats_preserved_Zsum_medianRank <- ROD_Dermo_Sus_stats_MR_order_less20[ROD_Dermo_Sus_stats_MR_order_less20$mod_name %in% ROD_Dermo_Sus_stats_preserved,]
## 1 modules, turqoise
#
#Dermo_Sus_ROD_Sus <-ROD_Dermo_Sus_stats_preserved[ROD_Dermo_Sus_stats_preserved$mod_name %in%Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
##0
#Dermo_Sus_ROD_Sus_ROD <-ROD_Dermo_Sus_stats_preserved[ROD_Dermo_Sus_stats_preserved$mod_name %in% ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
## tuqoise

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

### Assess module preservation between Probiotic and Pro_RE22 (just the probiotic network, not the full network)
#Pro_Pro_RE22_multiExpr = list(Pro=list(data=Probiotic_dds_rlog_matrix_common), Pro_RE22=list(data=Pro_RE22_dds_rlog_matrix_common_Pro))
#Pro_Pro_RE22_multiColor = list(Pro = Probiotic_moduleColors)
##Pro_Pro_RE22 = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Probiotic_Pro_module_preservation.RData")
##Pro_Pro_RE22_mp =modulePreservation(Pro_Pro_RE22_multiExpr, Pro_Pro_RE22_multiColor,
##                                     referenceNetworks=1,
##                                     verbose=3,
##                                     networkType="signed hybrid", # use same signed hybrid as before
##                                     corFnc = "bicor", # use recommended bicor as before 
##                                     nPermutations=100,
##                                     randomSeed = 1, # recommended in langfelder tutorial
##                                     quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
#Pro_Pro_RE22_stats = Pro_Pro_RE22_mp$preservation$Z$ref.Pro$inColumnsAlsoPresentIn.Pro_RE22
#Pro_Pro_RE22_stats_order <- Pro_Pro_RE22_stats[order(-Pro_Pro_RE22_stats[,2]),c(1:2)] 
#Pro_Pro_RE22_stats_preserved <- Pro_Pro_RE22_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
## mod_name moduleSize Zsummary.pres
## 1 turquoise       1000      7.814080
## 2     brown       1000      5.279559
#Pro_Pro_RE22_stats_MR <- Pro_Pro_RE22_mp$preservation$observed$ref.Pro$inColumnsAlsoPresentIn.Pro_RE22
#Pro_Pro_RE22_stats_MR <- Pro_Pro_RE22_stats_MR[,c("moduleSize","medianRank.pres")]
#Pro_Pro_RE22_stats_MR_order <- Pro_Pro_RE22_stats_MR[order(Pro_Pro_RE22_stats_MR[,2]),]
#Pro_Pro_RE22_stats_MR_order_less20 <- Pro_Pro_RE22_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)
#
## Save network
##save(Pro_Pro_RE22_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Probiotic_Pro_module_preservation.RData")
#
## Plot module preservation
## remember - low median rank means high preservation
#ggplot(Pro_Pro_RE22_stats_MR_order_less20  , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
#  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))
#
#Pro_Pro_RE22_pro <- Pro_Pro_RE22_stats_preserved[Pro_Pro_RE22_stats_preserved$mod_name %in% Probiotic_moduleTraitCor_Pval_df_sig_list_rm,]
## 0
#Pro_Pro_RE22_Pro_S4 <- Pro_Pro_RE22_stats_preserved[Pro_Pro_RE22_stats_preserved$mod_name %in% Pro_RE22_Pro_moduleTraitCor_Pval_df_sig_list_rm,]
## mod_name moduleSize Zsummary.pres
## 2    brown       1000      5.279559
#Pro_Pro_RE22_Pro_RI <-Pro_Pro_RE22_stats_preserved[Pro_Pro_RE22_stats_preserved$mod_name %in% Pro_RE22_Pro_RI_moduleTraitCor_Pval_df_sig_list_rm,]
## mod_name moduleSize Zsummary.pres
## 1 turquoise       1000      7.814080
## 2     brown       1000      5.279559
#
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

# We work with two sets:
nSets = 2
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("Probiotic", "Pro_RE22")
shortLabels = c("Probiotic", "Pro_RE22")
# Form multi-set expression data: columns starting from 9 contain actual expression data.
C_vir_Pro_multiExpr = vector(mode = "list", length = nSets)

C_vir_Pro_multiExpr[[1]] = list(data = as.data.frame(Probiotic_dds_rlog_matrix_common))
names(C_vir_Pro_multiExpr[[1]]$data) = row.names(Probiotic_dds_rlog_matrix_common)
rownames(C_vir_Pro_multiExpr[[1]]$data) = row.names(Probiotic_dds_rlog_matrix_common)

C_vir_Pro_multiExpr[[2]] = list(data = as.data.frame(Pro_RE22_dds_rlog_matrix_common_Pro))
names(C_vir_Pro_multiExpr[[2]]$data) = row.names(Pro_RE22_dds_rlog_matrix_common_Pro)
rownames(C_vir_Pro_multiExpr[[2]]$data) = row.names(Pro_RE22_dds_rlog_matrix_common_Pro)

# Check that the data has the correct format for many functions operating on multiple sets:
exprSize_Cvir_Pro = checkSets(C_vir_Pro_multiExpr)

# Check that all genes and samples have sufficiently low numbers of missing values.
#gsg = goodSamplesGenesMS(C_vir_Pro_multiExpr, verbose = 3)
gsg$allOK # [1] TRUE

# Assign clinical trait data
# Create data frame with the sample information 
colnames(Probiotic_coldata_collapsed)
colnames(Pro_RE22_coldata_Pro_collapsed)
# all the same

# Recode colnames to be just the sample name and the colname
C_vir_Pro_coldata <- rbind(Probiotic_coldata_collapsed,
                           Pro_RE22_coldata_Pro_collapsed)
levels(C_vir_Pro_coldata$Condition)
#[1] "Untreated_control"    "Bacillus_pumilus_RI0695" "Control_no_treatment"    "Bacillus_pumilus"        "Phaeobacter_inhibens"   

C_vir_Pro_coldata$Condition <- str_replace(C_vir_Pro_coldata$Condition, "Untreated_control", "Control") # ROD Res
C_vir_Pro_coldata$Condition <- str_replace(C_vir_Pro_coldata$Condition, "Bacillus_pumilus_RI0695", "Bacillus_pumilus") # ROD Res
C_vir_Pro_coldata$Condition <- str_replace(C_vir_Pro_coldata$Condition, "Control_no_treatment", "Control") # only recode the control for Pro_RE22, we want to keep 
C_vir_Pro_coldata <- C_vir_Pro_coldata %>% rownames_to_column("sample")
C_vir_Pro_coldata <- C_vir_Pro_coldata[,c("sample","Condition")]

C_vir_Pro_coldata_binarize <- binarizeCategoricalColumns.pairwise(C_vir_Pro_coldata[,2])
row.names(C_vir_Pro_coldata_binarize) <- C_vir_Pro_coldata$sample

# See how big the traits are and what are the trait and sample names
dim(C_vir_Pro_coldata)
names(C_vir_Pro_coldata)
C_vir_Pro_coldata$sample

# Form a multi-set structure that will hold the clinical traits.
C_vir_Pro_coldata_all = vector(mode="list", length = nSets)
setSamples = rownames(C_vir_Pro_multiExpr[[1]]$data)
traitRows = match(setSamples, row.names(C_vir_Pro_coldata_binarize))
C_vir_Pro_coldata_all[[1]] = list(data = as.data.frame(C_vir_Pro_coldata_binarize[traitRows,]))
rownames(C_vir_Pro_coldata_all[[1]]$data) = rownames(C_vir_Pro_multiExpr[[1]]$data)
colnames(C_vir_Pro_coldata_all[[1]]$data) = c("data.Control.vs.Bacillus_pumilus", "data.Phaeobacter_inhibens.vs.Bacillus_pumilus", "data.Phaeobacter_inhibens.vs.Control")

setSamples = rownames(C_vir_Pro_multiExpr[[2]]$data)
traitRows = match(setSamples, row.names(C_vir_Pro_coldata_binarize))
C_vir_Pro_coldata_all[[2]] = list(data = as.data.frame(C_vir_Pro_coldata_binarize[traitRows,]))

# Define data set dimensions
nGenes = exprSize_Cvir_Pro$nGenes
nSamples = exprSize_Cvir_Pro$nSamples

# Choose soft-thresholing powers
# Choose set based on where all three datasets level off - picking 9 since none fit the topology free network criterion

# Network construction 
#C_vir_Pro_net = blockwiseConsensusModules(
#  C_vir_Pro_multiExpr,
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
names(C_vir_Pro_net)
C_vir_Pro_consMEs = C_vir_Pro_net$multiMEs
C_vir_Pro_moduleLabels = C_vir_Pro_net$colors
# Convert the numeric labels to color labels
C_vir_Pro_moduleColors = labels2colors(C_vir_Pro_moduleLabels)
C_vir_Pro_consTree = C_vir_Pro_net$dendrograms[[1]]

save(C_vir_Pro_net, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/C_vir_Pro_consensus_modules_net.RData")

# plot dendrogram
sizeGrWindow(8,6)
plotDendroAndColors(C_vir_Pro_consTree, C_vir_Pro_moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()

# Relate consensus modules to external traits
# Set up variables to contain the module-trait correlations
C_vir_Pro_moduleTraitCor = list()
C_vir_Pro_moduleTraitPvalue = list()
# Calculate the correlations
C_vir_Pro_moduleTraitCor[[1]] = cor(C_vir_Pro_consMEs[[1]]$data, C_vir_Pro_coldata_all[[1]]$data, use = "p")
C_vir_Pro_moduleTraitPvalue[[1]] = corPvalueFisher(C_vir_Pro_moduleTraitCor[[1]], exprSize_Cvir_Pro$nSamples[1])
C_vir_Pro_moduleTraitCor[[2]] = cor(C_vir_Pro_consMEs[[2]]$data, C_vir_Pro_coldata_all[[2]]$data, use = "p")
C_vir_Pro_moduleTraitPvalue[[2]] = corPvalueFisher(C_vir_Pro_moduleTraitCor[[2]], exprSize_Cvir_Pro$nSamples[2])

## DISPLAY MODULE TRAIT CORRELATIONS FOR CONSENSUS MODULES, I care about columns 1 and 3 
# Need to fix and invesigate this!
# Initialize matrices to hold the consensus correlation and p-value
C_vir_Pro_consensusCor = matrix(NA, nrow(C_vir_Pro_moduleTraitCor[[1]]), ncol(C_vir_Pro_moduleTraitCor[[1]]))
C_vir_Pro_consensusPvalue = matrix(NA, nrow(C_vir_Pro_moduleTraitCor[[1]]), ncol(C_vir_Pro_moduleTraitCor[[1]]))
# Find consensus negative correlations
C_vir_Pro_negative = C_vir_Pro_moduleTraitCor[[1]] < 0 & C_vir_Pro_moduleTraitCor[[2]] < 0
C_vir_Pro_consensusCor[C_vir_Pro_negative] = pmax(C_vir_Pro_moduleTraitCor[[1]][C_vir_Pro_negative], C_vir_Pro_moduleTraitCor[[2]][C_vir_Pro_negative])
C_vir_Pro_consensusPvalue[C_vir_Pro_negative] = pmax(C_vir_Pro_moduleTraitPvalue[[1]][C_vir_Pro_negative], C_vir_Pro_moduleTraitPvalue[[2]][C_vir_Pro_negative])
# Find consensus positive correlations
C_vir_Pro_positive = C_vir_Pro_moduleTraitCor[[1]] > 0 & C_vir_Pro_moduleTraitCor[[2]] > 0
C_vir_Pro_consensusCor[C_vir_Pro_positive] = pmin(C_vir_Pro_moduleTraitCor[[1]][C_vir_Pro_positive], C_vir_Pro_moduleTraitCor[[2]][C_vir_Pro_positive])
C_vir_Pro_consensusPvalue[C_vir_Pro_positive] = pmax(C_vir_Pro_moduleTraitPvalue[[1]][C_vir_Pro_positive], C_vir_Pro_moduleTraitPvalue[[2]][C_vir_Pro_positive])

# Convert numerical lables to colors for labeling of modules in the plot
C_vir_Pro_MEColors = labels2colors(as.numeric(substring(names(C_vir_Pro_consMEs[[1]]$data), 3)));
C_vir_Pro_MEColorNames = paste("ME", C_vir_Pro_MEColors, sep="");

# Plot the consensus modules for a set, but is same across sets
C_vir_Pro_textMatrix = paste(signif(C_vir_Pro_consensusCor, 2), "\n(",
                             signif(C_vir_Pro_consensusPvalue, 1), ")", sep = "")
dim(C_vir_Pro_textMatrix) = dim(C_vir_Pro_moduleTraitCor[[1]])
sizeGrWindow(10,7)
#pdf(file = "Plots/ModuleTraitRelationships-consensus.pdf", wi = 10, he = 7);
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = C_vir_Pro_consensusCor,
               xLabels = names(C_vir_Pro_coldata_all[[1]]$data),
               yLabels = C_vir_Pro_MEColorNames,
               ySymbols =C_vir_Pro_MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = C_vir_Pro_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Consensus module--trait relationships across\n",
                            paste(setLabels, collapse = " and ")))

# Conclusion: no significant consensus modules across probiotic challenges


## IDENTIFY MODULES ENRICHED FOR APOPTOSIS USING anRichment


#### C. GIGAS WGCNA ####

#####  DATA FORMATTING, BATCH EFFECT REMOVAL ####

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
## How many modules identified
#table(Zhang_net$colors)
## Plot dendrogram with colors
## open a graphics window
#sizeGrWindow(12, 9)
## Convert labels to colors for plotting
#Zhang_mergedColors = labels2colors(Zhang_net$colors)
## Plot the dendrogram and the module colors underneath
#plotDendroAndColors(Zhang_net$dendrograms[[1]], Zhang_mergedColors[Zhang_net$blockGenes[[1]]],
#                    "Module colors",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05)
#Zhang_moduleLabels = Zhang_net$colors
#Zhang_moduleColors = labels2colors(Zhang_net$colors)
#Zhang_MEs = Zhang_net$MEs
#Zhang_geneTree = Zhang_net$dendrograms[[1]]
#
## Save Zhang network
#save(Zhang_net, Zhang_mergedColors, Zhang_moduleLabels, Zhang_moduleColors, Zhang_MEs,Zhang_geneTree, 
#     file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Zhang_network.RData")

# Load Zhang network
Zhang_network = load(file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/Zhang_network.RData")


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

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
deLorg_Res_moduleTraitCor_Pval_df_sig_list <- deLorg_Res_moduleTraitCor_Pval_df_sig$mod_names
deLorg_Res_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(deLorg_Res_moduleTraitCor_Pval_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
names(deLorg_Res_moduleTraitCor_Pval_df_sig_list_rm) <- c("orangered4" ,     "skyblue3"        ,"indianred2" ,     "lightgreen"   ,   "floralwhite"     ,"lightskyblue4" ,  "coral4"    ,     "lightslateblue" , "green"          ,
                                                           "turquoise" ,      "antiquewhite2"  , "blue"      ,      "plum1"       ,    "darkolivegreen2", "pink"         ,   "skyblue2" ,      "coral1"        ,  "lightpink2"     ,
                                                           "deeppink1" ,      "salmon4"        , "white"     ,      "greenyellow" ,    "royalblue"      , "brown"        ,   "grey" )
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

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
He_moduleTraitCor_Pval_df_sig_list <- He_moduleTraitCor_Pval_df_sig$mod_names
He_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(He_moduleTraitCor_Pval_df_sig_list, "ME")

# Use function to lookup all apop names for each significant module
names(He_moduleTraitCor_Pval_df_sig_list_rm) <- c("black"       ,"skyblue"  ,"darkolivegreen" , "magenta" , "white"   ,"darkgreen", "darkslateblue","salmon4",         "red" ,           
                                                   "tan"        , "yellow"  , "midnightblue"  ,  "salmon" ,  "lightsteelblue1", "brown4"  ,  "orangered4"  , "brown" ,          "blue",           
                                                   "darkorange2", "darkgrey", "turquoise"     ,  "grey"    )
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

## He vs. Delorgeril_Sus have conservation (since He was also susceptible)
He_deLorg_Sus_multiExpr = list(He=list(data=He_dds_vst_matrix_common),deLorg_Sus=list(data=deLorgeril_Susceptible_dds_vst_matrix_common))
He_deLorg_Sus_multiColor = list(He = He_moduleColors) 
He_deLorg_Sus_mp =modulePreservation(He_deLorg_Sus_multiExpr, He_deLorg_Sus_multiColor,
                                   referenceNetworks=1,
                                   verbose=3,
                                   networkType="signed hybrid", # use same signed hybrid as before
                                   corFnc = "bicor", # use recommended bicor as before 
                                   nPermutations=100,
                                   randomSeed = 1, # recommended in langfelder tutorial
                                   quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
He_deLorg_Sus_stats = He_deLorg_Sus_mp$preservation$Z$ref.He$inColumnsAlsoPresentIn.deLorg_Sus
He_deLorg_Sus_stats_order <- He_deLorg_Sus_stats[order(-He_deLorg_Sus_stats[,2]),c(1:2)]
He_deLorg_Sus_stats_preserved <- He_deLorg_Sus_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
He_deLorg_Sus_stats_MR <- He_deLorg_Sus_mp$preservation$observed$ref.He$inColumnsAlsoPresentIn.deLorg_Sus
He_deLorg_Sus_stats_MR <- He_deLorg_Sus_stats_MR  [,c("moduleSize","medianRank.pres")]
He_deLorg_Sus_stats_MR_order <- He_deLorg_Sus_stats_MR [order(-He_deLorg_Sus_stats_MR [,2]),]
He_deLorg_Sus_stats_MR_order_less20 <- He_deLorg_Sus_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save module preservation results
save(He_deLorg_Sus_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/He_deLorg_Sus_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(He_deLorg_Sus_stats_MR_order_less20, aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

He_deLorg_Sus_stats_preserved_Zsum_medianRank <- He_deLorg_Sus_stats_MR_order_less20[He_deLorg_Sus_stats_MR_order_less20$mod_name %in% deLorg_He_stats_preserved$mod_name,]

# Were any of these preserved modules significant in deLorg_Sus and He
He_deLorg_Sus_preserved_all <- He_deLorg_Sus_stats_preserved[He_deLorg_Sus_stats_preserved$mod_name %in% deLorg_Sus_moduleTraitCor_Pval_df_sig_list_rm,]

He_deLorg_Sus_he_preserved_all <- He_deLorg_Sus_stats_preserved[He_deLorg_Sus_stats_preserved$mod_name %in% He_moduleTraitCor_Pval_df_sig_list_rm,]


# Assess whether deLorg_Sus and deLorgTol have conservation 
deLorg_Sus_Tol_multiExpr = list(deLorg_Sus=list(data=deLorgeril_Susceptible_dds_vst_matrix_common),deLorg_Tol=list(data=deLorgeril_Resistant_dds_vst_matrix_common))
deLorg_Sus_Tol_multiColor = list(deLorg_Sus = deLorg_Sus_moduleColors) 
#deLorg_Sus_Tol_mp =modulePreservation(deLorg_Sus_Tol_multiExpr, deLorg_Sus_Tol_multiColor,
#                                 referenceNetworks=1,
#                                 verbose=3,
#                                 networkType="signed hybrid", # use same signed hybrid as before
#                                 corFnc = "bicor", # use recommended bicor as before 
#                                 nPermutations=100,
#                                 randomSeed = 1, # recommended in langfelder tutorial
#                                 quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
deLorg_Sus_Tol_stats = deLorg_Sus_Tol_mp$preservation$Z$ref.deLorg_Sus$inColumnsAlsoPresentIn.deLorg_Tol
deLorg_Sus_Tol_stats_order <- deLorg_Sus_Tol_stats[order(-deLorg_Sus_Tol_stats[,2]),c(1:2)]
deLorg_Sus_Tol_stats_preserved <- deLorg_Sus_Tol_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
deLorg_Sus_Tol_stats_MR <- deLorg_Sus_Tol_mp$preservation$observed$ref.deLorg_Sus$inColumnsAlsoPresentIn.deLorg_Tol
deLorg_Sus_Tol_stats_MR <- deLorg_Sus_Tol_stats_MR[,c("moduleSize","medianRank.pres")]
deLorg_Sus_Tol_stats_MR_order <- deLorg_Sus_Tol_stats_MR[order(-deLorg_Sus_Tol_stats_MR[,2]),]
deLorg_Sus_Tol_stats_MR_order_less20 <- deLorg_Sus_Tol_stats_MR_order %>% rownames_to_column("mod_name") %>% filter(medianRank.pres < 20)

# Save module preservation results
#save(deLorg_Sus_Tol_mp, file="/Volumes/My Passport for Mac/Chapter1_Apoptosis_Paper_Saved_DESeq_WGCNA_Data/deLorg_Sus_deLorg_Res_module_preservation.RData")

# Plot module preservation
# remember - low median rank means high preservation
ggplot(deLorg_Sus_Tol_stats_MR_order_less20 , aes(x= moduleSize, y = medianRank.pres, color = mod_name)) + 
  geom_point() + xlab("Module Size") + ylab("Median Rank") + geom_text(aes(label=mod_name))

deLorg_Sus_Tol_stats_preserved_Zsum_medianRank <- deLorg_Sus_Tol_stats_MR_order_less20[deLorg_Sus_Tol_stats_MR_order_less20$mod_name %in% deLorg_Sus_Tol_stats_preserved$mod_name,]

# Were any of these preserved modules significant in deLorg_Sus and He
deLorg_Sus_Tol_Suspreserved_all <- deLorg_Sus_Tol_stats_preserved[deLorg_Sus_Tol_stats_preserved$mod_name %in% deLorg_Sus_moduleTraitCor_Pval_df_sig_list_rm,]
# mod_name moduleSize Zsummary.pres
# 3  turquoise       1000     20.453268
# 4       blue       1000     17.552898
# 11    purple        208      6.741922
deLorg_Sus_Tol_Tol_preserved_all <- deLorg_Sus_Tol_stats_preserved[deLorg_Sus_Tol_stats_preserved$mod_name %in% deLorg_Res_moduleTraitCor_Pval_df_sig_list_rm,]
# mod_name moduleSize Zsummary.pres
# 2        green       1000     34.124863
# 3    turquoise       1000     20.453268
# 4         blue       1000     17.552898
# 6         pink        406     13.949355
# 9        brown       1000      9.691178
# 10 greenyellow        164      8.095585

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


#### PERFORM CONSENSUS NETWORK ANALYSIS ACROSS C.GIGAS NETWORKS BACTERIA ####

## Perform Functional Enrichment of Apoptosis genes using AnRichment 

# Tutorials: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/

