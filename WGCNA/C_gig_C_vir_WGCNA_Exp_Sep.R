# Script to run WGCNA on experiments separately for each species 
# Erin Roberts, PhD Candidate University of Rhode Island 
# 4/9/2020

# This script perform WGCNA analysis on C. virginica and C. gigas experiments separately 

### LOAD PACKAGES ####
library(tidyverse)
library(limma)
library(WGCNA)
library(cluster)
library(anRichment)
library(anRichmentMethods)
library(dplyr)
library(plyr)
# source("https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/installAnRichment.R"); installAnRichment(); 

#### Helpful tutorials for meta-analysis
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/ModulePreservation/
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/JMiller/Tutorial%20document.pdf
# WGCNA Main tutorial: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
# anRichment: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/#manualInstall
# Package FAQs with some quidelines for running: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html
# tutorials on module preservation: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/ModulePreservation/Tutorials/
# Code for Differential network analysis: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/DifferentialNetworkAnalysis/

#### WGCNA C_VIRGINICA ####
cor <- WGCNA::cor # make sure these are run every time!
options(stringsAsFactors = FALSE) # make sure these are run every time! 
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
class(Dermo_Susceptible_dds_vs_matrix)

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

# Do column names agree between all?
all(colnames(Dermo_Tolerant_dds_vst_matrix_common ) %in% colnames(Probiotic_dds_rlog_matrix_common)) # TRUE
all(colnames(Dermo_Tolerant_dds_vst_matrix_common ) == colnames(Probiotic_dds_rlog_matrix_common)) # TRUE 

#### SELECT SOFT THRESHOLDING POWER FOR EACH EXPERIMENT ####
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function, following general recommendations to set network type to "signed hybrid" and using the "bicor" correlation
Dermo_Tolerant_dds_vst_matrix_common_sft <- pickSoftThreshold(Dermo_Tolerant_dds_vst_matrix_common, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
Dermo_Susceptible_dds_vst_matrix_common_sft <- pickSoftThreshold(Dermo_Susceptible_dds_vst_matrix_common, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
Probiotic_dds_rlog_matrix_common_sft <- pickSoftThreshold(Probiotic_dds_rlog_matrix_common, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
# Warning message:
# executing %dopar% sequentially: no parallel backend registered 
#From Peter Langfelder: https://bioinformatics.stackexchange.com/questions/10555/r-wgcna-error-code:
#What you see is a warning, not an error. 
#Your calculation will run fine, just slower. Unless you see other errors, you should be able to complete all steps of the analysis.

ROD_Resistant_dds_rlog_matrix_common_sft <- pickSoftThreshold(ROD_Resistant_dds_rlog_matrix_common, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
ROD_Susceptible_dds_rlog_matrix_common_sft <- pickSoftThreshold(ROD_Susceptible_dds_rlog_matrix_common, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
Pro_RE22_dds_rlog_matrix_common_sft <- pickSoftThreshold(Pro_RE22_dds_rlog_matrix_common, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 

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

#### ONE STEP NETWORK CONSTRUCTION, MODULE DETECTION, MODULE DENDROGRAM INSPECTION #### 
Dermo_Tol_net = blockwiseModules(Dermo_Tolerant_dds_vst_matrix_common, power = 3, # picked suitable power in the code above 
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
table(Dermo_Tol_net$colors)
# Plot dendrogram with colors
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
Dermo_Tol_mergedColors = labels2colors(Dermo_Tol_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(Dermo_Tol_net$dendrograms[[1]], Dermo_Tol_mergedColors[Dermo_Tol_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
Dermo_Tol_moduleLabels = Dermo_Tol_net$colors
Dermo_Tol_moduleColors = labels2colors(Dermo_Tol_net$colors)
Dermo_Tol_MEs = Dermo_Tol_net$MEs
Dermo_Tol_geneTree = Dermo_Tol_net$dendrograms[[1]]

# Dermo Susceptible repeat 
Dermo_Sus_net = blockwiseModules(Dermo_Susceptible_dds_vst_matrix_common, power = 2, # picked suitable power in the code above 
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
table(Dermo_Sus_net$colors)
# Plot dendrogram with colors
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
Dermo_Sus_mergedColors = labels2colors(Dermo_Sus_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(Dermo_Sus_net$dendrograms[[1]], Dermo_Sus_mergedColors[Dermo_Sus_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
Dermo_Sus_moduleLabels = Dermo_Sus_net$colors
Dermo_Sus_moduleColors = labels2colors(Dermo_Sus_net$colors)
Dermo_Sus_MEs = Dermo_Sus_net$MEs
Dermo_Sus_geneTree = Dermo_Sus_net$dendrograms[[1]]

# Probiotic
Probiotic_net = blockwiseModules(Probiotic_dds_rlog_matrix_common, power = 9, # picked because less than 20 samples and isn't scale free
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
table(Probiotic_net$colors)
# Plot dendrogram with colors
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
Probiotic_mergedColors = labels2colors(Probiotic_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(Probiotic_net$dendrograms[[1]], Probiotic_mergedColors[Probiotic_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
Probiotic_moduleLabels = Probiotic_net$colors
Probiotic_moduleColors = labels2colors(Probiotic_net$colors)
Probiotic_MEs = Probiotic_net$MEs
Probiotic_geneTree = Probiotic_net$dendrograms[[1]]

# ROD Res
ROD_Res_net = blockwiseModules(ROD_Resistant_dds_rlog_matrix_common, power = 9, # picked because less than 20 samples and isn't scale free
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
table(ROD_Res_net$colors)
# Plot dendrogram with colors
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
ROD_Res_mergedColors = labels2colors(ROD_Res_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(ROD_Res_net$dendrograms[[1]], ROD_Res_mergedColors[ROD_Res_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
ROD_Res_moduleLabels = ROD_Res_net$colors
ROD_Res_moduleColors = labels2colors(ROD_Res_net$colors)
ROD_Res_MEs = ROD_Res_net$MEs
ROD_Res_geneTree = ROD_Res_net$dendrograms[[1]]

# ROD Sus
ROD_Sus_net = blockwiseModules(ROD_Susceptible_dds_rlog_matrix_common, power = 9, # picked because less than 20 samples and isn't scale free
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
table(ROD_Sus_net$colors)
# Plot dendrogram with colors
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
ROD_Sus_mergedColors = labels2colors(ROD_Sus_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(ROD_Sus_net$dendrograms[[1]], ROD_Sus_mergedColors[ROD_Sus_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
ROD_Sus_moduleLabels = ROD_Sus_net$colors
ROD_Sus_moduleColors = labels2colors(ROD_Sus_net$colors)
ROD_Sus_MEs = ROD_Sus_net$MEs
ROD_Sus_geneTree = ROD_Sus_net$dendrograms[[1]]

# Pro_RE22
Pro_RE22_net = blockwiseModules(Pro_RE22_dds_rlog_matrix_common, power = 8, 
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
table(Pro_RE22_net$colors)
# Plot dendrogram with colors
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
Pro_RE22_mergedColors = labels2colors(Pro_RE22_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(Pro_RE22_net$dendrograms[[1]], Pro_RE22_mergedColors[Pro_RE22_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
Pro_RE22_moduleLabels = Pro_RE22_net$colors
Pro_RE22_moduleColors = labels2colors(Pro_RE22_net$colors)
Pro_RE22_MEs = Pro_RE22_net$MEs
Pro_RE22_geneTree = Pro_RE22_net$dendrograms[[1]]

#### BINARIZE ALL CATEGORICAL VARIABLES TO TEST DISEASE CHALLENGE ASSOCIATIONS ####
# out = binarizeCategoricalVariable(x,
# includePairwise = TRUE,
# includeLevelVsAll = FALSE);

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

Dermo_Tol_MEthistle1      <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "thistle1"]
Dermo_Tol_MElightgreen    <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "lightgreen" ]
Dermo_Tol_MEpink          <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "pink"]
Dermo_Tol_MEroyalblue     <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "royalblue"]
Dermo_Tol_MEhoneydew1     <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "honeydew1"]
Dermo_Tol_MEdarkseagreen3 <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "darkseagreen3"]
Dermo_Tol_MEtan           <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "tan"]
Dermo_Tol_MEdarkviolet    <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "darkviolet"]
Dermo_Tol_MEorangered4    <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "orangered4"]
Dermo_Tol_MEyellow        <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "yellow"]
Dermo_Tol_MEantiquewhite2 <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "antiquewhite2"]
Dermo_Tol_MElightcyan1    <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "lightcyan1"]
Dermo_Tol_MEplum2         <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "plum2"]
Dermo_Tol_MEturquoise     <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "turquoise"]
Dermo_Tol_MElightsteelblue<- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "lightsteelblue"]

Dermo_Tol_MEthistle1_df <- as.data.frame(Dermo_Tol_MEthistle1)
colnames(Dermo_Tol_MEthistle1_df)[1] <- "ID"
Dermo_Tol_MEthistle1_annot_apop <- merge(Dermo_Tol_MEthistle1_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEthistle1_annot_apop) # 1

Dermo_Tol_MElightgreen_df <- as.data.frame(Dermo_Tol_MElightgreen)
colnames(Dermo_Tol_MElightgreen_df)[1] <- "ID"
Dermo_Tol_MElightgreen_annot_apop <- merge(Dermo_Tol_MElightgreen_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MElightgreen_annot_apop) # 2 GIMAP, MycA

Dermo_Tol_MEpink_df <- as.data.frame(Dermo_Tol_MEpink)
colnames(Dermo_Tol_MEpink_df)[1] <- "ID"
Dermo_Tol_MEpink_annot_apop <- merge(Dermo_Tol_MEpink_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEpink_annot_apop) # 6 AIF, RhoE, heat shock proteins

# royal blue is negative 
Dermo_Tol_MEroyalblue_df <- as.data.frame(Dermo_Tol_MEroyalblue)
colnames(Dermo_Tol_MEroyalblue_df)[1] <- "ID"
Dermo_Tol_MEroyalblue_annot_apop <- merge(Dermo_Tol_MEroyalblue_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEroyalblue_annot_apop) #1 ATF4

Dermo_Tol_MEhoneydew1_df <- as.data.frame(Dermo_Tol_MEhoneydew1)
colnames(Dermo_Tol_MEhoneydew1_df)[1] <- "ID"
Dermo_Tol_MEhoneydew1_annot_apop <- merge(Dermo_Tol_MEhoneydew1_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEhoneydew1_annot_apop) #6 Aven, BIRIAP, MAPK1, cAMP 1

Dermo_Tol_MEdarkseagreen3_df <- as.data.frame(Dermo_Tol_MEdarkseagreen3)
colnames(Dermo_Tol_MEdarkseagreen3_df)[1] <- "ID"
Dermo_Tol_MEdarkseagreen3_annot_apop <- merge(Dermo_Tol_MEdarkseagreen3_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEdarkseagreen3_annot_apop) #1 IP3R

Dermo_Tol_MEtan_df <- as.data.frame(Dermo_Tol_MEtan)
colnames(Dermo_Tol_MEtan_df)[1] <- "ID"
Dermo_Tol_MEtan_annot_apop <- merge(Dermo_Tol_MEtan_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEtan_annot_apop) #0

Dermo_Tol_MEdarkviolet_df <- as.data.frame(Dermo_Tol_MEdarkviolet)
colnames(Dermo_Tol_MEdarkviolet_df)[1] <- "ID"
Dermo_Tol_MEdarkviolet_annot_apop <- merge(Dermo_Tol_MEdarkviolet_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEdarkviolet_annot_apop) #3  BAG, MAPK14, heatshock70

Dermo_Tol_MEorangered4_df <- as.data.frame(Dermo_Tol_MEorangered4)
colnames(Dermo_Tol_MEorangered4_df)[1] <- "ID"
Dermo_Tol_MEorangered4_annot_apop <- merge(Dermo_Tol_MEorangered4_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEorangered4_annot_apop)  # 2 ceramide synthase, ras-like GTP RHoL

Dermo_Tol_MEyellow_df <- as.data.frame(Dermo_Tol_MEyellow)
colnames(Dermo_Tol_MEyellow_df)[1] <- "ID"
Dermo_Tol_MEyellow_annot_apop <- merge(Dermo_Tol_MEyellow_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEyellow_annot_apop)  # 12 TLR1, TLR2, TNFAIP, STAT5B, HTRA2, BTG1, cytoc, CREB3B

# antique white is negative
Dermo_Tol_MEantiquewhite2_df <- as.data.frame(Dermo_Tol_MEantiquewhite2)
colnames(Dermo_Tol_MEantiquewhite2_df)[1] <- "ID"
Dermo_Tol_MEantiquewhite2_annot_apop <- merge(Dermo_Tol_MEantiquewhite2_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEantiquewhite2_annot_apop)  # 1 TRAF6

Dermo_Tol_MElightcyan1_df <- as.data.frame(Dermo_Tol_MElightcyan1)
colnames(Dermo_Tol_MElightcyan1_df)[1] <- "ID"
Dermo_Tol_MElightcyan1_annot_apop <- merge(Dermo_Tol_MElightcyan1_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MElightcyan1_annot_apop)  # 2 STAT5B

Dermo_Tol_MEplum2_df <- as.data.frame(Dermo_Tol_MEplum2)
colnames(Dermo_Tol_MEplum2_df)[1] <- "ID"
Dermo_Tol_MEplum2_annot_apop <- merge(Dermo_Tol_MEplum2_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEplum2_annot_apop)  # 2 IP3R, dynamin 120

# turqoise is negative
Dermo_Tol_MEturquoise_df <- as.data.frame(Dermo_Tol_MEturquoise)
colnames(Dermo_Tol_MEturquoise_df)[1] <- "ID"
Dermo_Tol_MEturquoise_annot_apop <- merge(Dermo_Tol_MEturquoise_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEturquoise_annot_apop)  # 22 IFI44, GIMAP, caspase 1,7, calpain, PCDC, TLR4

# light steel blue is negative
Dermo_Tol_MElightsteelblue_df <- as.data.frame(Dermo_Tol_MElightsteelblue)
colnames(Dermo_Tol_MElightsteelblue_df)[1] <- "ID"
Dermo_Tol_MElightsteelblue_annot_apop <- merge(Dermo_Tol_MElightsteelblue_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MElightsteelblue_annot_apop)  # 0

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
Dermo_S_colorh = c("thistle1","lightgreen","pink","royalblue","honeydew1","darkseagreen3","tan","darkviolet",
"orangered4","yellow","antiquewhite2","lightcyan1","plum2","turquoise","lightsteelblue")

Dermo_Tol_Module_hub_genes <- chooseTopHubInEachModule(
  Dermo_Tolerant_dds_vst_matrix_common, 
  Dermo_Tol_colorh, 
  power = 3,  # power used for the adjacency network
  type = "signed hybrid", 
  corFnc = "bicor"
  )
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

Dermo_Sus_MEgreenyellow      <- colnames(Dermo_Susceptible_dds_vst_matrix_common)[Dermo_Sus_moduleColors == "greenyellow"]
Dermo_Sus_MElightgreen    <- colnames(Dermo_Susceptible_dds_vst_matrix_common)[Dermo_Sus_moduleColors == "lightgreen" ]


Dermo_Sus_MEgreenyellow_df <- as.data.frame(Dermo_Sus_MEgreenyellow)
colnames(Dermo_Sus_MEgreenyellow_df)[1] <- "ID"
Dermo_Sus_MEgreenyellow_annot_apop <- merge(Dermo_Sus_MEgreenyellow_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Sus_MEgreenyellow_annot_apop)  # 2 calpain 9, ATF-4

#  light gren is negative
Dermo_Sus_MElightgreen_df <- as.data.frame(Dermo_Sus_MElightgreen)
colnames(Dermo_Sus_MElightgreen_df)[1] <- "ID"
Dermo_Sus_MElightgreen_annot_apop <- merge(Dermo_Sus_MElightgreen_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Sus_MElightgreen_annot_apop)  # 3 DIAP2, heat shock protein, calpain 5 

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

Dermo_Sus_Module_hub_genes <- chooseTopHubInEachModule(
  Dermo_Susceptible_dds_vst_matrix_common, 
  Dermo_Sus_colorh, 
  power = 2,  # power used for the adjacency network
  type = "signed hybrid", 
  corFnc = "bicor"
)
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

Probiotic_moduleTraitCor_Pval_df_sig$mod_names
#"MEgrey60"     "MElightgreen"

# grey is negative
Probiotic_MEgrey60     <- colnames(Probiotic_dds_rlog_matrix_common)[Probiotic_moduleColors == "grey60"]
Probiotic_MElightgreen    <- colnames(Probiotic_dds_rlog_matrix_common)[Probiotic_moduleColors == "lightgreen" ]

Probiotic_MEgrey60_df <- as.data.frame(Probiotic_MEgrey60)
colnames(Probiotic_MEgrey60_df)[1] <- "ID"
Probiotic_MEgrey60_annot_apop <- merge(Probiotic_MEgrey60_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Probiotic_MEgrey60_annot_apop)  # 4, MAPK, PCDP

Probiotic_MElightgreen_df <- as.data.frame(Probiotic_MElightgreen)
colnames(Probiotic_MElightgreen_df)[1] <- "ID"
Probiotic_MElightgreen_annot_apop <- merge(Probiotic_MElightgreen_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Probiotic_MElightgreen_annot_apop)  # 1 dynamin 120

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

Probiotic_Module_hub_genes <- chooseTopHubInEachModule(
  Probiotic_dds_rlog_matrix_common, 
  Probiotic_colorh, 
  power = 2,  # power used for the adjacency network
  type = "signed hybrid", 
  corFnc = "bicor"
)
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

ROD_Res_moduleTraitCor_Pval_df_sig$mod_names
ROD_Res_moduleTraitCor_Pval_df_sig_list <- ROD_Res_moduleTraitCor_Pval_df_sig$mod_names
#MEplum1"  "MEorange"

# plum1 is negative
ROD_Res_MEplum1     <- colnames(ROD_Resistant_dds_rlog_matrix_common)[ROD_Res_moduleColors == "plum1"]
ROD_Res_MEorange  <- colnames(ROD_Resistant_dds_rlog_matrix_common)[ROD_Res_moduleColors == "orange" ]

ROD_Res_MEplum1_df <- as.data.frame(ROD_Res_MEplum1)
colnames(ROD_Res_MEplum1_df)[1] <- "ID"
ROD_Res_MEplum1_annot_apop <- merge(ROD_Res_MEplum1_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(ROD_Res_MEplum1_annot_apop)  # 2, TLR4

ROD_Res_MEorange_df <- as.data.frame(ROD_Res_MEorange)
colnames(ROD_Res_MEorange_df)[1] <- "ID"
ROD_Res_MEorange_annot_apop <- merge(ROD_Res_MEorange_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(ROD_Res_MEorange_annot_apop)  # 3 casp2, tlr3, rac gamma

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

ROD_Res_Module_hub_genes <- chooseTopHubInEachModule(
  ROD_Resistant_dds_rlog_matrix_common, 
  ROD_Res_colorh, 
  power = 2,  # power used for the adjacency network
  type = "signed hybrid", 
  corFnc = "bicor"
)
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
# MEturquoise  MEbrown (negative) 

ROD_Sus_MEturquoise  <- colnames(ROD_Susceptible_dds_rlog_matrix_common)[ROD_Sus_moduleColors == "turquoise"]
ROD_Sus_MEbrown   <- colnames(ROD_Susceptible_dds_rlog_matrix_common)[ROD_Sus_moduleColors == "brown" ]

ROD_Sus_MEturquoise_df <- as.data.frame(ROD_Sus_MEturquoise)
colnames(ROD_Sus_MEturquoise_df)[1] <- "ID"
ROD_Sus_MEturquoise_annot_apop <- merge(ROD_Sus_MEturquoise_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(ROD_Sus_MEturquoise_annot_apop)  # 46 BIRIAP, GIMAP, TRAF, TLR, caspase

ROD_Sus_MEbrown_df <- as.data.frame(ROD_Sus_MEbrown)
colnames(ROD_Sus_MEbrown_df)[1] <- "ID"
ROD_Sus_MEbrown_annot_apop <- merge(ROD_Sus_MEbrown_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(ROD_Sus_MEbrown_annot_apop)  # 63, IP3R, MAPK, caspase 2,3,8, 9, TRAF, calpain

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

ROD_Sus_Module_hub_genes <- chooseTopHubInEachModule(
  ROD_Susceptible_dds_rlog_matrix_common, 
  ROD_Sus_colorh, 
  power = 2,  # power used for the adjacency network
  type = "signed hybrid", 
  corFnc = "bicor"
)
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

# Which modules have the highest associations with  RI06p5 (high correlation and low P value)?
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
Pro_RE22_RI_moduleTraitCor_Pval_df_sig <- Pro_RE22_RI_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05 & Condition.Control_no_treatment.vs.Bacillus_pumilus > 0)
Pro_RE22_RI_moduleTraitCor_Pval_df_sig # 12

Pro_RE22_S4_moduleTraitCor_Pval_df[order(Pro_RE22_S4_moduleTraitCor_Pval_df$pvalue),]
class(Pro_RE22_S4_moduleTraitCor_Pval_df$pvalue) # numeric
Pro_RE22_S4_moduleTraitCor_Pval_df_sig <- Pro_RE22_S4_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05 & Condition.Phaeobacter_inhibens.vs.Control_no_treatment > 0)
Pro_RE22_S4_moduleTraitCor_Pval_df_sig # 8

Pro_RE22_moduleTraitCor_Pval_df[order(Pro_RE22_moduleTraitCor_Pval_df$pvalue),]
class(Pro_RE22_moduleTraitCor_Pval_df$pvalue) # numeric
Pro_RE22_moduleTraitCor_Pval_df_sig <- Pro_RE22_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05 & Condition.Vibrio_coralliilyticus_RE22_exposure_6h.vs.Control_no_treatment > 0)
Pro_RE22_moduleTraitCor_Pval_df_sig #17

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES

Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list <- Pro_RE22_RI_moduleTraitCor_Pval_df_sig$mod_names
  #[1] "MEmediumpurple3" "MEblack"         "MEsteelblue"     "MEpurple"        "MEsalmon"        "MEmidnightblue"  "MEplum1"         "MEcoral1"       
  #[9] "MEsalmon4"       "MEblue"          "MEred"           "MEorange"     
Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list <- Pro_RE22_S4_moduleTraitCor_Pval_df_sig$mod_names
  #[1] "MEdarkolivegreen" "MEdarkorange"     "MEcyan"           "MEgreenyellow"    "MEmagenta"        "MEgreen"          "MElightgreen"     "MEsaddlebrown" 
Pro_RE22_moduleTraitCor_Pval_df_sig_list <- Pro_RE22_moduleTraitCor_Pval_df_sig$mod_names
  # [1] "MEdarkseagreen4"  "MEdarkolivegreen" "MEcyan"           "MEmagenta"        "MEgreen"          "MElightgreen"     "MEsaddlebrown"    "MElavenderblush3"
  # [9] "MEnavajowhite2"   "MEpink"           "MElightpink4"     "MEcoral2"         "MEmaroon"         "MEdarkred"        "MEdarkgreen"      "MEyellow"        
  #[17] "MEgrey"   

Pro_RE22_MEmediumpurple3 <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "mediumpurple3"]
Pro_RE22_MEblack <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "black"    ]
Pro_RE22_MEsteelblue <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "steelblue"   ]
Pro_RE22_MEpurple <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "purple"      ]
Pro_RE22_MEsalmon <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "salmon"      ]
Pro_RE22_MEmidnightblue <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "midnightblue" ]
Pro_RE22_MEplum1 <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "plum1"        ]
Pro_RE22_MEcoral1 <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "coral1"       ]
Pro_RE22_MEsalmon4 <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "salmon4"    ]
Pro_RE22_MEblue <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "blue"      ]
Pro_RE22_MEred <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "red"     ]
Pro_RE22_MEorange <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "orange"     ]

Pro_RE22_MEdarkolivegreen <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "darkolivegreen" ]
Pro_RE22_MEdarkorange <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "darkorange"     ]
Pro_RE22_MEcyan <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "cyan"          ]
Pro_RE22_MEgreenyellow  <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "greenyellow"  ]
Pro_RE22_MEmagenta <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "magenta"       ]
Pro_RE22_MEgreen  <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "green"         ]
Pro_RE22_MElightgreen <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "lightgreen"   ]
Pro_RE22_MEsaddlebrown <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "saddlebrown" ]
Pro_RE22_MEdarkseagreen4 <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "darkseagreen4" ]
Pro_RE22_MEdarkolivegreen  <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "darkolivegreen" ]
Pro_RE22_MElightgreen   <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "lightgreen"     ]
Pro_RE22_MElavenderblush3 <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "lavenderblush3"]
Pro_RE22_MEnavajowhite2   <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "navajowhite2"   ]
Pro_RE22_MEpink <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "pink"        ]
Pro_RE22_MElightpink4 <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "lightpink4"   ]
Pro_RE22_MEcoral2 <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "coral2"      ]
Pro_RE22_MEmaroon <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "maroon"      ]
Pro_RE22_MEdarkred <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "darkred"    ]
Pro_RE22_MEdarkgreen <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "darkgreen"  ]
Pro_RE22_MEyellow <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "yellow"        ]
Pro_RE22_MEgrey <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "grey" ]

Pro_RE22_MEmediumpurple3_df <- as.data.frame(Pro_RE22_MEmediumpurple3)
Pro_RE22_MEblack_df <- as.data.frame(Pro_RE22_MEblack)
Pro_RE22_MEsteelblue_df <- as.data.frame(Pro_RE22_MEsteelblue)
Pro_RE22_MEpurple_df <- as.data.frame(Pro_RE22_MEpurple)
Pro_RE22_MEsalmon_df <- as.data.frame(Pro_RE22_MEsalmon)
Pro_RE22_MEmidnightblue_df <- as.data.frame(Pro_RE22_MEmidnightblue)
Pro_RE22_MEplum1_df <- as.data.frame(Pro_RE22_MEplum1)
Pro_RE22_MEcoral1_df <- as.data.frame(Pro_RE22_MEcoral1)
Pro_RE22_MEsalmon4_df <- as.data.frame(Pro_RE22_MEsalmon4)
Pro_RE22_MEblue_df <- as.data.frame(Pro_RE22_MEblue)
Pro_RE22_MEred_df <- as.data.frame(Pro_RE22_MEred)
Pro_RE22_MEorange_df <- as.data.frame(Pro_RE22_MEorange)
Pro_RE22_MEdarkolivegreen_df <- as.data.frame(Pro_RE22_MEdarkolivegreen )
Pro_RE22_MEdarkorange_df <- as.data.frame(Pro_RE22_MEdarkorange)
Pro_RE22_MEcyan_df <- as.data.frame(Pro_RE22_MEcyan)
Pro_RE22_MEgreenyellow_df <- as.data.frame(Pro_RE22_MEgreenyellow)
Pro_RE22_MEmagenta_df <- as.data.frame(Pro_RE22_MEmagenta)
Pro_RE22_MEgreen_df <- as.data.frame(Pro_RE22_MEgreen)
Pro_RE22_MElightgreen_df <- as.data.frame(Pro_RE22_MElightgreen)
Pro_RE22_MEsaddlebrown_df <- as.data.frame(Pro_RE22_MEsaddlebrown )
Pro_RE22_MEdarkseagreen4_df <- as.data.frame(Pro_RE22_MEdarkseagreen4)
Pro_RE22_MEdarkolivegreen_df <- as.data.frame(Pro_RE22_MEdarkolivegreen)
Pro_RE22_MElightgreen_df <- as.data.frame(Pro_RE22_MElightgreen)
Pro_RE22_MElavenderblush3_df <- as.data.frame(Pro_RE22_MElavenderblush3)
Pro_RE22_MEnavajowhite2_df <- as.data.frame(Pro_RE22_MEnavajowhite2 )
Pro_RE22_MEpink_df <- as.data.frame(Pro_RE22_MEpink)
Pro_RE22_MElightpink4_df <- as.data.frame(Pro_RE22_MElightpink4)
Pro_RE22_MEcoral2_df <- as.data.frame(Pro_RE22_MEcoral2)
Pro_RE22_MEmaroon_df <- as.data.frame(Pro_RE22_MEmaroon)
Pro_RE22_MEdarkred_df <- as.data.frame(Pro_RE22_MEdarkred)
Pro_RE22_MEdarkgreen_df <- as.data.frame(Pro_RE22_MEdarkgreen)
Pro_RE22_MEyellow_df <- as.data.frame(Pro_RE22_MEyellow)
Pro_RE22_MEgrey_df <- as.data.frame(Pro_RE22_MEgrey)

colnames(Pro_RE22_MEmediumpurple3_df)[1] <- "ID"
colnames(Pro_RE22_MEblack_df)[1] <- "ID"
colnames(Pro_RE22_MEsteelblue_df)[1] <- "ID"
colnames(Pro_RE22_MEpurple_df)[1] <- "ID"
colnames(Pro_RE22_MEsalmon_df)[1] <- "ID"
colnames(Pro_RE22_MEmidnightblue_df)[1] <- "ID"
colnames(Pro_RE22_MEplum1_df)[1] <- "ID"
colnames(Pro_RE22_MEcoral1_df)[1] <- "ID"
colnames(Pro_RE22_MEsalmon4_df)[1] <- "ID"
colnames(Pro_RE22_MEblue_df)[1] <- "ID"
colnames(Pro_RE22_MEred_df)[1] <- "ID"
colnames(Pro_RE22_MEorange_df)[1] <- "ID"
colnames(Pro_RE22_MEdarkolivegreen_df)[1] <- "ID"
colnames(Pro_RE22_MEdarkorange_df)[1] <- "ID"
colnames(Pro_RE22_MEcyan_df)[1] <- "ID"
colnames(Pro_RE22_MEgreenyellow_df)[1] <- "ID"
colnames(Pro_RE22_MEmagenta_df)[1] <- "ID"
colnames(Pro_RE22_MEgreen_df)[1] <- "ID"
colnames(Pro_RE22_MElightgreen_df)[1] <- "ID"
colnames(Pro_RE22_MEsaddlebrown_df)[1] <- "ID"
colnames(Pro_RE22_MEdarkseagreen4_df)[1] <- "ID"
colnames(Pro_RE22_MEdarkolivegreen_df)[1] <- "ID"
colnames(Pro_RE22_MElightgreen_df)[1] <- "ID"
colnames(Pro_RE22_MElavenderblush3_df)[1] <- "ID"
colnames(Pro_RE22_MEnavajowhite2_df)[1] <- "ID"
colnames(Pro_RE22_MEpink_df)[1] <- "ID"
colnames(Pro_RE22_MElightpink4_df)[1] <- "ID"
colnames(Pro_RE22_MEcoral2_df)[1] <- "ID"
colnames(Pro_RE22_MEmaroon_df)[1] <- "ID"
colnames(Pro_RE22_MEdarkred_df)[1] <- "ID"
colnames(Pro_RE22_MEdarkgreen_df)[1] <- "ID"
colnames(Pro_RE22_MEyellow_df)[1] <- "ID"
colnames(Pro_RE22_MEgrey_df)[1] <- "ID"

Pro_RE22_MEmediumpurple3_apop <- merge(Pro_RE22_MEmediumpurple3_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEblack_apop <- merge(Pro_RE22_MEblack_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEsteelblue_apop <- merge(Pro_RE22_MEsteelblue_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEpurple_apop <- merge(Pro_RE22_MEpurple_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEsalmon_apop <- merge(Pro_RE22_MEsalmon_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEmidnightblue_apop <- merge(Pro_RE22_MEmidnightblue_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEplum1_apop <- merge(Pro_RE22_MEplum1_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEcoral1_apop <- merge(Pro_RE22_MEcoral1_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEsalmon4_apop <- merge(Pro_RE22_MEsalmon4_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEblue_apop <- merge(Pro_RE22_MEblue_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEred_apop <- merge(Pro_RE22_MEred_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEorange_apop <- merge(Pro_RE22_MEorange_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEdarkolivegreen_apop <- merge(Pro_RE22_MEdarkolivegreen_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEdarkorange_apop <- merge(Pro_RE22_MEdarkorange_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEcyan_apop <- merge(Pro_RE22_MEcyan_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEgreenyellow_apop <- merge(Pro_RE22_MEgreenyellow_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEmagenta_apop <- merge(Pro_RE22_MEmagenta_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEgreen_apop <- merge(Pro_RE22_MEgreen_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MElightgreen_apop <- merge(Pro_RE22_MElightgreen_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEsaddlebrown_apop <- merge(Pro_RE22_MEsaddlebrown_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEdarkseagreen4_apop <- merge(Pro_RE22_MEdarkseagreen4_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEdarkolivegreen_apop <- merge(Pro_RE22_MEdarkolivegreen_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MElightgreen_apop <- merge(Pro_RE22_MElightgreen_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MElavenderblush3_apop <- merge(Pro_RE22_MElavenderblush3_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEnavajowhite2_apop <- merge(Pro_RE22_MEnavajowhite2_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEpink_apop <- merge(Pro_RE22_MEpink_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MElightpink4_apop <- merge(Pro_RE22_MElightpink4_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEcoral2_apop <- merge(Pro_RE22_MEcoral2_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEmaroon_apop <- merge(Pro_RE22_MEmaroon_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEdarkred_apop <- merge(Pro_RE22_MEdarkred_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEdarkgreen_apop <- merge(Pro_RE22_MEdarkgreen_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEyellow_apop <- merge(Pro_RE22_MEyellow_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEgrey_apop <- merge(Pro_RE22_MEgrey_df, C_vir_rtracklayer_apop_product_final, by = "ID")

nrow(Pro_RE22_MEmediumpurple3_apop) # 4 Rho IFI21
nrow(Pro_RE22_MEblack_apop) # 12 TLR6, programmed cell death protein, STAT5A
nrow(Pro_RE22_MEsteelblue_apop) # 1 cell death-inducing p53-target protein 1
nrow(Pro_RE22_MEpurple_apop) # 4, growth arrest and DNA damage inducible protein GADD45, TGF beta, heat shock protein 
nrow(Pro_RE22_MEsalmon_apop) # 6, CD151, NFkB inhibition, cathepsin B, LITAF
nrow(Pro_RE22_MEmidnightblue_apop) # 2 FADD, calpain 5
nrow(Pro_RE22_MEplum1_apop) # 1 cathepsin L
nrow(Pro_RE22_MEcoral1_apop) # 0
nrow(Pro_RE22_MEsalmon4_apop) # 1 PCDP 10
nrow(Pro_RE22_MEblue_apop) # 23 TRAF, BIRIAP, GIMAP, calpain
nrow(Pro_RE22_MEred_apop) # 12, TNFL, IP3R, PCDP, IFI44, poly ADP polymerase, mitochondrial heat shock proteins
nrow(Pro_RE22_MEorange_apop) # 2 MAPK, IP3R
nrow(Pro_RE22_MEdarkolivegreen_apop) # 3 tollo
nrow(Pro_RE22_MEdarkorange_apop) # 3 caap1
nrow(Pro_RE22_MEcyan_apop) #7 MyD88, lymphotoxin alpha, caspase 7, tLr1, TRAF3 
nrow(Pro_RE22_MEgreenyellow_apop) # 3 IP3R, GADD45, TLR10
nrow(Pro_RE22_MEmagenta_apop) # bcl1, CREB, NFkB, rhoE, calpain 3
nrow(Pro_RE22_MEgreen_apop) # 15, TRAF BIR IAP, CRADD-like, IL17, NR13 like, CD151, TNFRSF16
nrow(Pro_RE22_MElightgreen_apop) # 6, TRAF3, JNK, caspase 8
nrow(Pro_RE22_MEsaddlebrown_apop) #3  GADD 45, TLR2, CREB
nrow(Pro_RE22_MEdarkseagreen4_apop) # 1 TNFAIP3
nrow(Pro_RE22_MElavenderblush3_apop) #0
nrow(Pro_RE22_MEnavajowhite2_apop) #0
nrow(Pro_RE22_MEpink_apop) #9 NFkB, calpain, BIRIAP, JNK, pyrin, UNC5C, MYD88
nrow(Pro_RE22_MElightpink4_apop) #3 STAT5A, P53 damage regulated protein 1, MAPK1
nrow(Pro_RE22_MEcoral2_apop) #0
nrow(Pro_RE22_MEmaroon_apop) # 1 cell death inducing p53 target
nrow(Pro_RE22_MEdarkred_apop) # 3, TLR6, BAG, IP3R
nrow(Pro_RE22_MEdarkgreen_apop) # 1 TLR3
nrow(Pro_RE22_MEyellow_apop) #13 caspase8, BIRIAP, aurora kinase, PDCP 6, TLR4, TLR2, ligase CHIP, caspase9  (definitely extrinsic pathway usage)
nrow(Pro_RE22_MEgrey_apop) # 51 GIMAP, BIRIAP, casp2,3, MAPK, TRAF, programmed cell death proteins , TLRs 

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
          
Pro_RE22_Module_hub_genes <- chooseTopHubInEachModule(
  Pro_RE22_dds_rlog_matrix_common, 
  Pro_RE22_colorh, 
  power = 2,  # power used for the adjacency network
  type = "signed hybrid", 
  corFnc = "bicor"
)
class(Pro_RE22_Module_hub_genes)
Pro_RE22_Module_hub_genes_df <- as.data.frame(Pro_RE22_Module_hub_genes)
colnames(Pro_RE22_Module_hub_genes_df)[1] <- "ID"
Pro_RE22_Module_hub_genes <- merge(Pro_RE22_Module_hub_genes_df, C_vir_rtracklayer, by = "ID")
nrow(Pro_RE22_Module_hub_genes) 
# includes programmed cell death protein 6 as a hub gene, apoptosis regulatory protein Siva 

#### MEASURE MODULE PRESERVATION BETWEEN C. VIRGINICA DIFFERENT EXPERIMENTS ####

# Column names agree between the datasets, so we can go ahead and measure module preservation between all 
## 1. MEASURE OVERLAP BETWEEN RE22 AND RODs

# Assess whether Pro_RE22 and ROD_Res modules are preserved
Pro_RE22_ROD_multiExpr = list(Pro_RE22=list(data=Pro_RE22_dds_rlog_matrix_common),ROD_Res=list(data=ROD_Resistant_dds_rlog_matrix_common))
Pro_RE22_multiColor = list(Pro_RE22 = Pro_RE22_moduleColors) 
Pro_RE22_mp =modulePreservation(Pro_RE22_ROD_multiExpr, Pro_RE22_multiColor,
                      referenceNetworks=1,
                      verbose=3,
                      networkType="signed hybrid", # use same signed hybrid as before
                      corFnc = "bicor", # use recommended bicor as before 
                      nPermutations=200,
                      randomSeed = 1, # recommended in langfelder tutorial
                      quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Pro_RE22_stats = Pro_RE22_mp$preservation$Z$ref.Pro_RE22$inColumnsAlsoPresentIn.ROD_Res
Pro_RE22_stats_order <- Pro_RE22_stats[order(-Pro_RE22_stats[,2]),c(1:2)]
Pro_RE22_stats_preserved <- Pro_RE22_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
## Moderate preservation of between Pro_RE22 and ROD Res
# brown                 1000     7.8915583
# red                   1000     7.1531975
# steelblue              175     6.5939398
# black                 1000     6.3774745
# midnightblue           527     5.7505155
# purple                 731     5.7203439
# blue                  1000     5.1007921
# pink                   787     5.0031835

# Were any of these preserved modules significant in Pro_RE22 or ROD_Res?
Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list, "ME")
Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list, "ME")
Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Pro_RE22_moduleTraitCor_Pval_df_sig_list, "ME")
ROD_Res_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(ROD_Res_moduleTraitCor_Pval_df_sig_list, "ME")

Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm <- as.data.frame(Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm)
Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm <- as.data.frame(Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm)
Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm <- as.data.frame(Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm)
ROD_Res_moduleTraitCor_Pval_df_sig_list_rm <- as.data.frame(  ROD_Res_moduleTraitCor_Pval_df_sig_list_rm)
  
colnames(Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm)[1] <- "mod_name"
colnames(Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm)[1] <- "mod_name"
colnames(Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm)[1] <- "mod_name"
colnames(ROD_Res_moduleTraitCor_Pval_df_sig_list_rm)  [1] <- "mod_name"

Pro_RE22_RI_SIG_ROD_Res_preserved <-merge(Pro_RE22_stats_preserved, Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm)
  # 6 : black,  blue, midnightblue ,purple,red ,steelblue 
Pro_RE22_S4_SIG_ROD_Res_preserved <-merge(Pro_RE22_stats_preserved, Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm)
  # 0
Pro_RE22_RE22_SIG_ROD_Res_preserved <-merge(Pro_RE22_stats_preserved, Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm)
  # pink
Pro_RE22_ROD_RES_SIG_preserved <-merge(Pro_RE22_stats_preserved, ROD_Res_moduleTraitCor_Pval_df_sig_list_rm)
  # 0 

# Assess whether Pro_RE22 and ROD_Sus modules are preserved
Pro_RE22_Sus_multiExpr = list(Pro_RE22=list(data=Pro_RE22_dds_rlog_matrix_common),ROD_Sus=list(data=ROD_Susceptible_dds_rlog_matrix_common))
Pro_RE22_Sus_mp =modulePreservation(Pro_RE22_Sus_multiExpr, Pro_RE22_multiColor,
                                referenceNetworks=1,
                                verbose=3,
                                networkType="signed hybrid", # use same signed hybrid as before
                                corFnc = "bicor", # use recommended bicor as before 
                                nPermutations=100,
                                randomSeed = 1, # recommended in langfelder tutorial
                                quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Pro_RE22_Sus_stats = Pro_RE22_Sus_mp$preservation$Z$ref.Pro_RE22$inColumnsAlsoPresentIn.ROD_Sus 
Pro_RE22_Sus_stats_order <- Pro_RE22_Sus_stats[order(-Pro_RE22_Sus_stats[,2]),c(1:2)]
Pro_RE22_Sus_stats_preserved <- Pro_RE22_Sus_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)

# Were any of these preserved modules significant in Pro_RE22 or ROD_Sus?
ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(ROD_Sus_moduleTraitCor_Pval_df_sig_list, "ME")
ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm <- as.data.frame(  ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm)
colnames(ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm)  [1] <- "mod_name"

Pro_RE22_RI_SIG_ROD_Sus_preserved <-merge(Pro_RE22_Sus_stats_preserved, Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm)
# blue, red
Pro_RE22_S4_SIG_ROD_Sus_preserved <-merge(Pro_RE22_Sus_stats_preserved, Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm)
# 0
Pro_RE22_RE22_SIG_ROD_Sus_preserved <-merge(Pro_RE22_Sus_stats_preserved, Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm)
# 0
Pro_RE22_ROD_SUS_SIG_preserved <-merge(Pro_RE22_Sus_stats_preserved, ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm)
# brown

# Assess whether Pro_RE22 and Dermo_Tol modules are preserved
Pro_RE22_Dermo_multiExpr = list(Pro_RE22=list(data=Pro_RE22_dds_rlog_matrix_common),Dermo_Tol=list(data=Dermo_Tolerant_dds_vst_matrix_common))
Pro_RE22_multiColor = list(Pro_RE22 = Pro_RE22_moduleColors) 
Pro_RE22_Dermo_mp =modulePreservation(Pro_RE22_Dermo_multiExpr, Pro_RE22_multiColor,
                                referenceNetworks=1,
                                verbose=3,
                                networkType="signed hybrid", # use same signed hybrid as before
                                corFnc = "bicor", # use recommended bicor as before 
                                nPermutations=100,
                                randomSeed = 1, # recommended in langfelder tutorial
                                quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Pro_RE22_Dermo_stats = Pro_RE22_Dermo_mp$preservation$Z$ref.Pro_RE22$inColumnsAlsoPresentIn.Dermo_Tol
Pro_RE22_Dermo_stats_order <- Pro_RE22_Dermo_stats[order(-Pro_RE22_Dermo_stats[,2]),c(1:2)]
Pro_RE22_Dermo_stats_preserved <- Pro_RE22_Dermo_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)

# Were any of these preserved modules significant in Pro_RE22 or Dermo_Tol?
Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Dermo_Tol_moduleTraitCor_Pval_df_sig_list, "ME")
Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm <- as.data.frame(  Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm)
colnames(Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm)  [1] <- "mod_name"

Pro_RE22_RI_SIG_Dermo_Tol_preserved <-merge(Pro_RE22_Dermo_stats_preserved, Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm)
#1   black       1000      5.354414
#2         blue       1000      6.243046
#3 midnightblue        527      5.580933
#4       purple        731      9.399914
#5          red       1000      8.038084
#6       salmon        615     10.598431
#7    steelblue        175      9.376713
Pro_RE22_S4_SIG_Dermo_Tol_preserved <-merge(Pro_RE22_Dermo_stats_preserved, Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm)
# green
Pro_RE22_RE22_SIG_RDermo_Tol_preserved <-merge(Pro_RE22_Dermo_stats_preserved, Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm)
# green yellow
Pro_RE22_Dermo_Tol_SIG_preserved <-merge(Pro_RE22_Dermo_stats_preserved, Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm)
# turqoise yellow

# Assess whether Pro_RE22 and Dermo_Sus modules are preserved
Pro_RE22_Dermo_Sus_multiExpr = list(Pro_RE22=list(data=Pro_RE22_dds_rlog_matrix_common),Dermo_Sus=list(data=Dermo_Susceptible_dds_vst_matrix_common))
Pro_RE22_multiColor = list(Pro_RE22 = Pro_RE22_moduleColors) 
Pro_RE22_Dermo_Sus_mp =modulePreservation(Pro_RE22_Dermo_Sus_multiExpr, Pro_RE22_multiColor,
                                      referenceNetworks=1,
                                      verbose=3,
                                      networkType="signed hybrid", # use same signed hybrid as before
                                      corFnc = "bicor", # use recommended bicor as before 
                                      nPermutations=100,
                                      randomSeed = 1, # recommended in langfelder tutorial
                                      quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Pro_RE22_Dermo_Sus_stats = Pro_RE22_Dermo_Sus_mp$preservation$Z$ref.Pro_RE22$inColumnsAlsoPresentIn.Dermo_Sus
Pro_RE22_Dermo_Sus_stats_order <- Pro_RE22_Dermo_Sus_stats[order(-Pro_RE22_Dermo_Sus_stats[,2]),c(1:2)]
Pro_RE22_Dermo_Sus_stats_preserved <- Pro_RE22_Dermo_Sus_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)

# Were any of these preserved modules significant in Pro_RE22 or Dermo_Sus?
Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Dermo_Sus_moduleTraitCor_Pval_df_sig_list, "ME")
Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm <- as.data.frame(  Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm)
colnames(Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm)  [1] <- "mod_name"

Pro_RE22_RI_SIG_Dermo_Sus_preserved <-merge(Pro_RE22_Dermo_Sus_stats_preserved, Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm)
# mod_name moduleSize Zsummary.pres
# 1      blue       1000      6.176316
# 2    purple        731      6.821283
# 3       red       1000      6.060967
# 4    salmon        615      9.461935
# 5 steelblue        175      7.608015
Pro_RE22_S4_SIG_Dermo_Sus_preserved <-merge(Pro_RE22_Dermo_Sus_stats_preserved, Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm)
# 0
Pro_RE22_RE22_SIG_Dermo_Sus_preserved <-merge(Pro_RE22_Dermo_Sus_stats_preserved, Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm)
# 0 
Pro_RE22_Dermo_Sus_SIG_preserved <-merge(Pro_RE22_Dermo_Sus_stats_preserved, Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm)
# 0 

# Assess whether Dermo_Tol and ROD_Res modules are preserved
ROD_Dermo_multiExpr = list(Dermo_Tol=list(data=Dermo_Tolerant_dds_vst_matrix_common), ROD_Res=list(data=ROD_Resistant_dds_rlog_matrix_common))
ROD_Dermo_multiColor = list(Dermo_Tol = Dermo_Tol_moduleColors)
ROD_Dermo_mp =modulePreservation(ROD_Dermo_multiExpr, ROD_Dermo_multiColor,
                                 referenceNetworks=1,
                                 verbose=3,
                                 networkType="signed hybrid", # use same signed hybrid as before
                                 corFnc = "bicor", # use recommended bicor as before 
                                 nPermutations=100,
                                 randomSeed = 1, # recommended in langfelder tutorial
                                 quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
ROD_Dermo_stats = ROD_Dermo_mp$preservation$Z$ref.Dermo_Tol$inColumnsAlsoPresentIn.ROD_Res
ROD_Dermo_stats_order <- ROD_Dermo_stats[order(-ROD_Dermo_stats[,2]),c(1:2)]
ROD_Dermo_stats_preserved <- ROD_Dermo_stats %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)

# Were any of these preserved modules significant in Pro_RE22 or Dermo_Tol?
Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Dermo_Tol_moduleTraitCor_Pval_df_sig_list, "ME")
Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm <- as.data.frame(  Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm)
colnames(Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm)  [1] <- "mod_name"

Dermo_Tol_ROD_Res <-merge(ROD_Dermo_stats_preserved, Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm)
# darkviolet, lightcyan1, pink, plum2,tan,turqoise
Dermo_Tol_ROD_Res_ROD <-merge(ROD_Dermo_stats_preserved, ROD_Res_moduleTraitCor_Pval_df_sig_list_rm)
#0

# Assess whether Dermo_Sus and ROD_Sus modules are preserved
ROD_Dermo_Sus_multiExpr = list(Dermo_Sus=list(data=Dermo_Susceptible_dds_vst_matrix_common), ROD_Sus=list(data=ROD_Susceptible_dds_rlog_matrix_common))
ROD_Dermo_Sus_multiColor = list(Dermo_Sus = Dermo_Sus_moduleColors)
ROD_Dermo_Sus_mp =modulePreservation(ROD_Dermo_Sus_multiExpr, ROD_Dermo_Sus_multiColor,
                                     referenceNetworks=1,
                                     verbose=3,
                                     networkType="signed hybrid", # use same signed hybrid as before
                                     corFnc = "bicor", # use recommended bicor as before 
                                     nPermutations=100,
                                     randomSeed = 1, # recommended in langfelder tutorial
                                     quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
ROD_Dermo_Sus_stats = ROD_Dermo_Sus_mp$preservation$Z$ref.Dermo_Sus$inColumnsAlsoPresentIn.ROD_Sus 
ROD_Dermo_Sus_stats_order <- ROD_Dermo_Sus_stats[order(-ROD_Dermo_Sus_stats[,2]),c(1:2)] 
ROD_Dermo_Sus_stats_preserved <- ROD_Dermo_Sus_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)

Dermo_Sus_ROD_Sus <-merge(ROD_Dermo_Sus_stats_preserved , Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm)
# 0
Dermo_Sus_ROD_Sus_ROD <-merge(ROD_Dermo_Sus_stats_preserved, ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm)
# turqoise

#### SUMMARY TABLE OF APOPTOSIS GENES ACROSS SIGNIFICANT MODULES IN EACH EXPERIMENT ####

#### HEATMAP OF MODULES OF INTEREST ACROSS EXPERIMENTS ####




#### PERFORM CONSENSUS NETWORK ANALYSIS ACROSS C.VIR NETWORKS ####
# CAN DO THIS LATER TO VALIDATE MODULE PRESERVATION IF NECESSARY 



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

# Call the network topology analysis function, following general recommendations to set network type to "signed hybrid" and using the "bicor" correlation
Zhang_dds_rlog_matrix_common_sft <- pickSoftThreshold(Zhang_dds_rlog_matrix_common, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
Rubio_dds_rlog_matrix_common_sft <- pickSoftThreshold(Rubio_dds_rlog_matrix_common, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
deLorgeril_Resistant_dds_vst_matrix_common_sft <- pickSoftThreshold(deLorgeril_Resistant_dds_vst_matrix_common, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
deLorgeril_Susceptible_dds_vst_matrix_common_sft <- pickSoftThreshold(deLorgeril_Susceptible_dds_vst_matrix_common, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
He_dds_vst_matrix_common_sft <- pickSoftThreshold(He_dds_vst_matrix_common, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")

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
Zhang_net = blockwiseModules(Zhang_dds_rlog_matrix_common, power = 4, # picked suitable power in the code above 
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

# Rubio
Rubio_net = blockwiseModules(Rubio_dds_rlog_matrix_common, power = 3, # picked suitable power in the code above 
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
table(Rubio_net$colors)
# Plot dendrogram with colors
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
Rubio_mergedColors = labels2colors(Rubio_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(Rubio_net$dendrograms[[1]], Rubio_mergedColors[Rubio_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
Rubio_moduleLabels = Rubio_net$colors
Rubio_moduleColors = labels2colors(Rubio_net$colors)
Rubio_MEs = Rubio_net$MEs
Rubio_geneTree = Rubio_net$dendrograms[[1]]

# deLorg Res
deLorg_Res_net = blockwiseModules(deLorgeril_Resistant_dds_vst_matrix_common, power = 4, # picked suitable power in the code above 
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
table(deLorg_Res_net$colors)
# Plot dendrogram with colors
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
deLorg_Res_mergedColors = labels2colors(deLorg_Res_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(deLorg_Res_net$dendrograms[[1]], deLorg_Res_mergedColors[deLorg_Res_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
deLorg_Res_moduleLabels = deLorg_Res_net$colors
deLorg_Res_moduleColors = labels2colors(deLorg_Res_net$colors)
deLorg_Res_MEs = deLorg_Res_net$MEs
deLorg_Res_geneTree = deLorg_Res_net$dendrograms[[1]]

# deLorg Sus
deLorg_Sus_net = blockwiseModules(deLorgeril_Susceptible_dds_vst_matrix_common, power = 8, # picked suitable power in the code above 
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
table(deLorg_Sus_net$colors)
# Plot dendrogram with colors
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
deLorg_Sus_mergedColors = labels2colors(deLorg_Sus_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(deLorg_Sus_net$dendrograms[[1]], deLorg_Sus_mergedColors[deLorg_Sus_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
deLorg_Sus_moduleLabels = deLorg_Sus_net$colors
deLorg_Sus_moduleColors = labels2colors(deLorg_Sus_net$colors)
deLorg_Sus_MEs = deLorg_Sus_net$MEs
deLorg_Sus_geneTree = deLorg_Sus_net$dendrograms[[1]]

# He 
He_dds_vst_matrix_common
He_net = blockwiseModules(He_dds_vst_matrix_common, power = 6, # picked suitable power in the code above 
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
table(He_net$colors)
# Plot dendrogram with colors
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
He_mergedColors = labels2colors(He_net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(He_net$dendrograms[[1]], He_mergedColors[He_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
He_moduleLabels = He_net$colors
He_moduleColors = labels2colors(He_net$colors)
He_MEs = He_net$MEs
He_geneTree = He_net$dendrograms[[1]]

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
He_dds_vst_matrix_common

#### QUANTIFYING MODULE TRAIT ASSOCIATIONS ####

## Zhang
Zhang_dds_rlog_matrix_common
Zhang_coldata

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
Zhang_moduleTraitCor_df <- Zhang_moduleTraitCor_df[,c("mod_names","")]
Zhang_moduleTraitPvalue_df <- as.data.frame(Dermo_Tol_moduleTraitPvalue)
Zhang_moduleTraitPvalue_df$mod_names <- row.names(Dermo_Tol_moduleTraitPvalue_df)
Zhang_moduleTraitPvalue_df <- Zhang_moduleTraitPvalue_df[,c("mod_names","")]
colnames(Zhang_moduleTraitPvalue_df)[2] <- "pvalue"

Zhang_moduleTraitCor_Pval_df <- join(Zhang_moduleTraitCor_df, Zhang_moduleTraitPvalue_df, by = "mod_names")

# Significantly correlated modules
Zhang_moduleTraitCor_Pval_df[order(Zhang_moduleTraitCor_Pval_df$pvalue),]
class(Zhang_moduleTraitCor_Pval_df$pvalue) # numeric
Zhang_moduleTraitCor_Pval_df_sig <- Zhang_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05)
Zhang_moduleTraitCor_Pval_df_sig

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES
Zhang_moduleTraitCor_Pval_df_sig_list <- Zhang_moduleTraitCor_Pval_df_sig$mod_names

Dermo_Tol_MEthistle1      <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "thistle1"]

Dermo_Tol_MEthistle1_df <- as.data.frame(Dermo_Tol_MEthistle1)
colnames(Dermo_Tol_MEthistle1_df)[1] <- "ID"
Dermo_Tol_MEthistle1_annot_apop <- merge(Dermo_Tol_MEthistle1_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEthistle1_annot_apop) # 1


## Gene relationship to trait and important modules: Gene Significance and Module Membership

# Define variable injected 
Zhang_injection = as.data.frame(Dermo_Tolerant_coldata_collapsed_binarize$Condition.Injected.vs.Control);
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
Dermo_S_colorh = c("thistle1","lightgreen","pink","royalblue","honeydew1","darkseagreen3","tan","darkviolet",
                   "orangered4","yellow","antiquewhite2","lightcyan1","plum2","turquoise","lightsteelblue")

Dermo_Tol_Module_hub_genes <- chooseTopHubInEachModule(
  Dermo_Tolerant_dds_vst_matrix_common, 
  Dermo_Tol_colorh, 
  power = 3,  # power used for the adjacency network
  type = "signed hybrid", 
  corFnc = "bicor"
)
class(Dermo_Tol_Module_hub_genes)
Dermo_Tol_Module_hub_genes_df <- as.data.frame(Dermo_Tol_Module_hub_genes)
colnames(Dermo_Tol_Module_hub_genes_df)[1] <- "ID"
Dermo_Tol_Module_hub_genes_apop <- merge(Dermo_Tol_Module_hub_genes_df, C_vir_rtracklayer, by = "ID")
nrow(Dermo_Tol_Module_hub_genes_apop) # 15



## Rubio
Rubio_dds_rlog_matrix_common)

## deLorg Res
deLorgeril_Resistant_dds_vst_matrix_common)

## deLorg Sus
deLorgeril_Susceptible_dds_vst_matrix_common

## He
He_dds_vst_matrix_common)



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

Dermo_Tol_MEthistle1      <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "thistle1"]
Dermo_Tol_MElightgreen    <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "lightgreen" ]
Dermo_Tol_MEpink          <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "pink"]
Dermo_Tol_MEroyalblue     <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "royalblue"]
Dermo_Tol_MEhoneydew1     <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "honeydew1"]
Dermo_Tol_MEdarkseagreen3 <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "darkseagreen3"]
Dermo_Tol_MEtan           <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "tan"]
Dermo_Tol_MEdarkviolet    <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "darkviolet"]
Dermo_Tol_MEorangered4    <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "orangered4"]
Dermo_Tol_MEyellow        <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "yellow"]
Dermo_Tol_MEantiquewhite2 <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "antiquewhite2"]
Dermo_Tol_MElightcyan1    <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "lightcyan1"]
Dermo_Tol_MEplum2         <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "plum2"]
Dermo_Tol_MEturquoise     <- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "turquoise"]
Dermo_Tol_MElightsteelblue<- colnames(Dermo_Tolerant_dds_vst_matrix_common)[Dermo_Tol_moduleColors == "lightsteelblue"]

Dermo_Tol_MEthistle1_df <- as.data.frame(Dermo_Tol_MEthistle1)
colnames(Dermo_Tol_MEthistle1_df)[1] <- "ID"
Dermo_Tol_MEthistle1_annot_apop <- merge(Dermo_Tol_MEthistle1_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEthistle1_annot_apop) # 1

Dermo_Tol_MElightgreen_df <- as.data.frame(Dermo_Tol_MElightgreen)
colnames(Dermo_Tol_MElightgreen_df)[1] <- "ID"
Dermo_Tol_MElightgreen_annot_apop <- merge(Dermo_Tol_MElightgreen_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MElightgreen_annot_apop) # 2 GIMAP, MycA

Dermo_Tol_MEpink_df <- as.data.frame(Dermo_Tol_MEpink)
colnames(Dermo_Tol_MEpink_df)[1] <- "ID"
Dermo_Tol_MEpink_annot_apop <- merge(Dermo_Tol_MEpink_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEpink_annot_apop) # 6 AIF, RhoE, heat shock proteins

# royal blue is negative 
Dermo_Tol_MEroyalblue_df <- as.data.frame(Dermo_Tol_MEroyalblue)
colnames(Dermo_Tol_MEroyalblue_df)[1] <- "ID"
Dermo_Tol_MEroyalblue_annot_apop <- merge(Dermo_Tol_MEroyalblue_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEroyalblue_annot_apop) #1 ATF4

Dermo_Tol_MEhoneydew1_df <- as.data.frame(Dermo_Tol_MEhoneydew1)
colnames(Dermo_Tol_MEhoneydew1_df)[1] <- "ID"
Dermo_Tol_MEhoneydew1_annot_apop <- merge(Dermo_Tol_MEhoneydew1_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEhoneydew1_annot_apop) #6 Aven, BIRIAP, MAPK1, cAMP 1

Dermo_Tol_MEdarkseagreen3_df <- as.data.frame(Dermo_Tol_MEdarkseagreen3)
colnames(Dermo_Tol_MEdarkseagreen3_df)[1] <- "ID"
Dermo_Tol_MEdarkseagreen3_annot_apop <- merge(Dermo_Tol_MEdarkseagreen3_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEdarkseagreen3_annot_apop) #1 IP3R

Dermo_Tol_MEtan_df <- as.data.frame(Dermo_Tol_MEtan)
colnames(Dermo_Tol_MEtan_df)[1] <- "ID"
Dermo_Tol_MEtan_annot_apop <- merge(Dermo_Tol_MEtan_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEtan_annot_apop) #0

Dermo_Tol_MEdarkviolet_df <- as.data.frame(Dermo_Tol_MEdarkviolet)
colnames(Dermo_Tol_MEdarkviolet_df)[1] <- "ID"
Dermo_Tol_MEdarkviolet_annot_apop <- merge(Dermo_Tol_MEdarkviolet_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEdarkviolet_annot_apop) #3  BAG, MAPK14, heatshock70

Dermo_Tol_MEorangered4_df <- as.data.frame(Dermo_Tol_MEorangered4)
colnames(Dermo_Tol_MEorangered4_df)[1] <- "ID"
Dermo_Tol_MEorangered4_annot_apop <- merge(Dermo_Tol_MEorangered4_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEorangered4_annot_apop)  # 2 ceramide synthase, ras-like GTP RHoL

Dermo_Tol_MEyellow_df <- as.data.frame(Dermo_Tol_MEyellow)
colnames(Dermo_Tol_MEyellow_df)[1] <- "ID"
Dermo_Tol_MEyellow_annot_apop <- merge(Dermo_Tol_MEyellow_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEyellow_annot_apop)  # 12 TLR1, TLR2, TNFAIP, STAT5B, HTRA2, BTG1, cytoc, CREB3B

# antique white is negative
Dermo_Tol_MEantiquewhite2_df <- as.data.frame(Dermo_Tol_MEantiquewhite2)
colnames(Dermo_Tol_MEantiquewhite2_df)[1] <- "ID"
Dermo_Tol_MEantiquewhite2_annot_apop <- merge(Dermo_Tol_MEantiquewhite2_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEantiquewhite2_annot_apop)  # 1 TRAF6

Dermo_Tol_MElightcyan1_df <- as.data.frame(Dermo_Tol_MElightcyan1)
colnames(Dermo_Tol_MElightcyan1_df)[1] <- "ID"
Dermo_Tol_MElightcyan1_annot_apop <- merge(Dermo_Tol_MElightcyan1_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MElightcyan1_annot_apop)  # 2 STAT5B

Dermo_Tol_MEplum2_df <- as.data.frame(Dermo_Tol_MEplum2)
colnames(Dermo_Tol_MEplum2_df)[1] <- "ID"
Dermo_Tol_MEplum2_annot_apop <- merge(Dermo_Tol_MEplum2_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEplum2_annot_apop)  # 2 IP3R, dynamin 120

# turqoise is negative
Dermo_Tol_MEturquoise_df <- as.data.frame(Dermo_Tol_MEturquoise)
colnames(Dermo_Tol_MEturquoise_df)[1] <- "ID"
Dermo_Tol_MEturquoise_annot_apop <- merge(Dermo_Tol_MEturquoise_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MEturquoise_annot_apop)  # 22 IFI44, GIMAP, caspase 1,7, calpain, PCDC, TLR4

# light steel blue is negative
Dermo_Tol_MElightsteelblue_df <- as.data.frame(Dermo_Tol_MElightsteelblue)
colnames(Dermo_Tol_MElightsteelblue_df)[1] <- "ID"
Dermo_Tol_MElightsteelblue_annot_apop <- merge(Dermo_Tol_MElightsteelblue_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Tol_MElightsteelblue_annot_apop)  # 0

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
Dermo_S_colorh = c("thistle1","lightgreen","pink","royalblue","honeydew1","darkseagreen3","tan","darkviolet",
                   "orangered4","yellow","antiquewhite2","lightcyan1","plum2","turquoise","lightsteelblue")

Dermo_Tol_Module_hub_genes <- chooseTopHubInEachModule(
  Dermo_Tolerant_dds_vst_matrix_common, 
  Dermo_Tol_colorh, 
  power = 3,  # power used for the adjacency network
  type = "signed hybrid", 
  corFnc = "bicor"
)
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

Dermo_Sus_MEgreenyellow      <- colnames(Dermo_Susceptible_dds_vst_matrix_common)[Dermo_Sus_moduleColors == "greenyellow"]
Dermo_Sus_MElightgreen    <- colnames(Dermo_Susceptible_dds_vst_matrix_common)[Dermo_Sus_moduleColors == "lightgreen" ]


Dermo_Sus_MEgreenyellow_df <- as.data.frame(Dermo_Sus_MEgreenyellow)
colnames(Dermo_Sus_MEgreenyellow_df)[1] <- "ID"
Dermo_Sus_MEgreenyellow_annot_apop <- merge(Dermo_Sus_MEgreenyellow_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Sus_MEgreenyellow_annot_apop)  # 2 calpain 9, ATF-4

#  light gren is negative
Dermo_Sus_MElightgreen_df <- as.data.frame(Dermo_Sus_MElightgreen)
colnames(Dermo_Sus_MElightgreen_df)[1] <- "ID"
Dermo_Sus_MElightgreen_annot_apop <- merge(Dermo_Sus_MElightgreen_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Dermo_Sus_MElightgreen_annot_apop)  # 3 DIAP2, heat shock protein, calpain 5 

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

Dermo_Sus_Module_hub_genes <- chooseTopHubInEachModule(
  Dermo_Susceptible_dds_vst_matrix_common, 
  Dermo_Sus_colorh, 
  power = 2,  # power used for the adjacency network
  type = "signed hybrid", 
  corFnc = "bicor"
)
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

Probiotic_moduleTraitCor_Pval_df_sig$mod_names
#"MEgrey60"     "MElightgreen"

# grey is negative
Probiotic_MEgrey60     <- colnames(Probiotic_dds_rlog_matrix_common)[Probiotic_moduleColors == "grey60"]
Probiotic_MElightgreen    <- colnames(Probiotic_dds_rlog_matrix_common)[Probiotic_moduleColors == "lightgreen" ]

Probiotic_MEgrey60_df <- as.data.frame(Probiotic_MEgrey60)
colnames(Probiotic_MEgrey60_df)[1] <- "ID"
Probiotic_MEgrey60_annot_apop <- merge(Probiotic_MEgrey60_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Probiotic_MEgrey60_annot_apop)  # 4, MAPK, PCDP

Probiotic_MElightgreen_df <- as.data.frame(Probiotic_MElightgreen)
colnames(Probiotic_MElightgreen_df)[1] <- "ID"
Probiotic_MElightgreen_annot_apop <- merge(Probiotic_MElightgreen_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(Probiotic_MElightgreen_annot_apop)  # 1 dynamin 120

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

Probiotic_Module_hub_genes <- chooseTopHubInEachModule(
  Probiotic_dds_rlog_matrix_common, 
  Probiotic_colorh, 
  power = 2,  # power used for the adjacency network
  type = "signed hybrid", 
  corFnc = "bicor"
)
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

ROD_Res_moduleTraitCor_Pval_df_sig$mod_names
ROD_Res_moduleTraitCor_Pval_df_sig_list <- ROD_Res_moduleTraitCor_Pval_df_sig$mod_names
#MEplum1"  "MEorange"

# plum1 is negative
ROD_Res_MEplum1     <- colnames(ROD_Resistant_dds_rlog_matrix_common)[ROD_Res_moduleColors == "plum1"]
ROD_Res_MEorange  <- colnames(ROD_Resistant_dds_rlog_matrix_common)[ROD_Res_moduleColors == "orange" ]

ROD_Res_MEplum1_df <- as.data.frame(ROD_Res_MEplum1)
colnames(ROD_Res_MEplum1_df)[1] <- "ID"
ROD_Res_MEplum1_annot_apop <- merge(ROD_Res_MEplum1_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(ROD_Res_MEplum1_annot_apop)  # 2, TLR4

ROD_Res_MEorange_df <- as.data.frame(ROD_Res_MEorange)
colnames(ROD_Res_MEorange_df)[1] <- "ID"
ROD_Res_MEorange_annot_apop <- merge(ROD_Res_MEorange_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(ROD_Res_MEorange_annot_apop)  # 3 casp2, tlr3, rac gamma

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

ROD_Res_Module_hub_genes <- chooseTopHubInEachModule(
  ROD_Resistant_dds_rlog_matrix_common, 
  ROD_Res_colorh, 
  power = 2,  # power used for the adjacency network
  type = "signed hybrid", 
  corFnc = "bicor"
)
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
# MEturquoise  MEbrown (negative) 

ROD_Sus_MEturquoise  <- colnames(ROD_Susceptible_dds_rlog_matrix_common)[ROD_Sus_moduleColors == "turquoise"]
ROD_Sus_MEbrown   <- colnames(ROD_Susceptible_dds_rlog_matrix_common)[ROD_Sus_moduleColors == "brown" ]

ROD_Sus_MEturquoise_df <- as.data.frame(ROD_Sus_MEturquoise)
colnames(ROD_Sus_MEturquoise_df)[1] <- "ID"
ROD_Sus_MEturquoise_annot_apop <- merge(ROD_Sus_MEturquoise_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(ROD_Sus_MEturquoise_annot_apop)  # 46 BIRIAP, GIMAP, TRAF, TLR, caspase

ROD_Sus_MEbrown_df <- as.data.frame(ROD_Sus_MEbrown)
colnames(ROD_Sus_MEbrown_df)[1] <- "ID"
ROD_Sus_MEbrown_annot_apop <- merge(ROD_Sus_MEbrown_df, C_vir_rtracklayer_apop_product_final, by = "ID")
nrow(ROD_Sus_MEbrown_annot_apop)  # 63, IP3R, MAPK, caspase 2,3,8, 9, TRAF, calpain

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

ROD_Sus_Module_hub_genes <- chooseTopHubInEachModule(
  ROD_Susceptible_dds_rlog_matrix_common, 
  ROD_Sus_colorh, 
  power = 2,  # power used for the adjacency network
  type = "signed hybrid", 
  corFnc = "bicor"
)
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

# Which modules have the highest associations with  RI06p5 (high correlation and low P value)?
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
Pro_RE22_RI_moduleTraitCor_Pval_df_sig <- Pro_RE22_RI_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05 & Condition.Control_no_treatment.vs.Bacillus_pumilus > 0)
Pro_RE22_RI_moduleTraitCor_Pval_df_sig # 12

Pro_RE22_S4_moduleTraitCor_Pval_df[order(Pro_RE22_S4_moduleTraitCor_Pval_df$pvalue),]
class(Pro_RE22_S4_moduleTraitCor_Pval_df$pvalue) # numeric
Pro_RE22_S4_moduleTraitCor_Pval_df_sig <- Pro_RE22_S4_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05 & Condition.Phaeobacter_inhibens.vs.Control_no_treatment > 0)
Pro_RE22_S4_moduleTraitCor_Pval_df_sig # 8

Pro_RE22_moduleTraitCor_Pval_df[order(Pro_RE22_moduleTraitCor_Pval_df$pvalue),]
class(Pro_RE22_moduleTraitCor_Pval_df$pvalue) # numeric
Pro_RE22_moduleTraitCor_Pval_df_sig <- Pro_RE22_moduleTraitCor_Pval_df %>% filter(pvalue <= 0.05 & Condition.Vibrio_coralliilyticus_RE22_exposure_6h.vs.Control_no_treatment > 0)
Pro_RE22_moduleTraitCor_Pval_df_sig #17

### ANNOTATE APOPTOSIS GENES IN SIGNIFICANT MODULES

Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list <- Pro_RE22_RI_moduleTraitCor_Pval_df_sig$mod_names
#[1] "MEmediumpurple3" "MEblack"         "MEsteelblue"     "MEpurple"        "MEsalmon"        "MEmidnightblue"  "MEplum1"         "MEcoral1"       
#[9] "MEsalmon4"       "MEblue"          "MEred"           "MEorange"     
Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list <- Pro_RE22_S4_moduleTraitCor_Pval_df_sig$mod_names
#[1] "MEdarkolivegreen" "MEdarkorange"     "MEcyan"           "MEgreenyellow"    "MEmagenta"        "MEgreen"          "MElightgreen"     "MEsaddlebrown" 
Pro_RE22_moduleTraitCor_Pval_df_sig_list <- Pro_RE22_moduleTraitCor_Pval_df_sig$mod_names
# [1] "MEdarkseagreen4"  "MEdarkolivegreen" "MEcyan"           "MEmagenta"        "MEgreen"          "MElightgreen"     "MEsaddlebrown"    "MElavenderblush3"
# [9] "MEnavajowhite2"   "MEpink"           "MElightpink4"     "MEcoral2"         "MEmaroon"         "MEdarkred"        "MEdarkgreen"      "MEyellow"        
#[17] "MEgrey"   

Pro_RE22_MEmediumpurple3 <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "mediumpurple3"]
Pro_RE22_MEblack <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "black"    ]
Pro_RE22_MEsteelblue <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "steelblue"   ]
Pro_RE22_MEpurple <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "purple"      ]
Pro_RE22_MEsalmon <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "salmon"      ]
Pro_RE22_MEmidnightblue <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "midnightblue" ]
Pro_RE22_MEplum1 <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "plum1"        ]
Pro_RE22_MEcoral1 <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "coral1"       ]
Pro_RE22_MEsalmon4 <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "salmon4"    ]
Pro_RE22_MEblue <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "blue"      ]
Pro_RE22_MEred <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "red"     ]
Pro_RE22_MEorange <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "orange"     ]

Pro_RE22_MEdarkolivegreen <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "darkolivegreen" ]
Pro_RE22_MEdarkorange <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "darkorange"     ]
Pro_RE22_MEcyan <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "cyan"          ]
Pro_RE22_MEgreenyellow  <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "greenyellow"  ]
Pro_RE22_MEmagenta <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "magenta"       ]
Pro_RE22_MEgreen  <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "green"         ]
Pro_RE22_MElightgreen <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "lightgreen"   ]
Pro_RE22_MEsaddlebrown <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "saddlebrown" ]
Pro_RE22_MEdarkseagreen4 <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "darkseagreen4" ]
Pro_RE22_MEdarkolivegreen  <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "darkolivegreen" ]
Pro_RE22_MElightgreen   <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "lightgreen"     ]
Pro_RE22_MElavenderblush3 <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "lavenderblush3"]
Pro_RE22_MEnavajowhite2   <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "navajowhite2"   ]
Pro_RE22_MEpink <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "pink"        ]
Pro_RE22_MElightpink4 <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "lightpink4"   ]
Pro_RE22_MEcoral2 <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "coral2"      ]
Pro_RE22_MEmaroon <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "maroon"      ]
Pro_RE22_MEdarkred <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "darkred"    ]
Pro_RE22_MEdarkgreen <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "darkgreen"  ]
Pro_RE22_MEyellow <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "yellow"        ]
Pro_RE22_MEgrey <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "grey" ]

Pro_RE22_MEmediumpurple3_df <- as.data.frame(Pro_RE22_MEmediumpurple3)
Pro_RE22_MEblack_df <- as.data.frame(Pro_RE22_MEblack)
Pro_RE22_MEsteelblue_df <- as.data.frame(Pro_RE22_MEsteelblue)
Pro_RE22_MEpurple_df <- as.data.frame(Pro_RE22_MEpurple)
Pro_RE22_MEsalmon_df <- as.data.frame(Pro_RE22_MEsalmon)
Pro_RE22_MEmidnightblue_df <- as.data.frame(Pro_RE22_MEmidnightblue)
Pro_RE22_MEplum1_df <- as.data.frame(Pro_RE22_MEplum1)
Pro_RE22_MEcoral1_df <- as.data.frame(Pro_RE22_MEcoral1)
Pro_RE22_MEsalmon4_df <- as.data.frame(Pro_RE22_MEsalmon4)
Pro_RE22_MEblue_df <- as.data.frame(Pro_RE22_MEblue)
Pro_RE22_MEred_df <- as.data.frame(Pro_RE22_MEred)
Pro_RE22_MEorange_df <- as.data.frame(Pro_RE22_MEorange)
Pro_RE22_MEdarkolivegreen_df <- as.data.frame(Pro_RE22_MEdarkolivegreen )
Pro_RE22_MEdarkorange_df <- as.data.frame(Pro_RE22_MEdarkorange)
Pro_RE22_MEcyan_df <- as.data.frame(Pro_RE22_MEcyan)
Pro_RE22_MEgreenyellow_df <- as.data.frame(Pro_RE22_MEgreenyellow)
Pro_RE22_MEmagenta_df <- as.data.frame(Pro_RE22_MEmagenta)
Pro_RE22_MEgreen_df <- as.data.frame(Pro_RE22_MEgreen)
Pro_RE22_MElightgreen_df <- as.data.frame(Pro_RE22_MElightgreen)
Pro_RE22_MEsaddlebrown_df <- as.data.frame(Pro_RE22_MEsaddlebrown )
Pro_RE22_MEdarkseagreen4_df <- as.data.frame(Pro_RE22_MEdarkseagreen4)
Pro_RE22_MEdarkolivegreen_df <- as.data.frame(Pro_RE22_MEdarkolivegreen)
Pro_RE22_MElightgreen_df <- as.data.frame(Pro_RE22_MElightgreen)
Pro_RE22_MElavenderblush3_df <- as.data.frame(Pro_RE22_MElavenderblush3)
Pro_RE22_MEnavajowhite2_df <- as.data.frame(Pro_RE22_MEnavajowhite2 )
Pro_RE22_MEpink_df <- as.data.frame(Pro_RE22_MEpink)
Pro_RE22_MElightpink4_df <- as.data.frame(Pro_RE22_MElightpink4)
Pro_RE22_MEcoral2_df <- as.data.frame(Pro_RE22_MEcoral2)
Pro_RE22_MEmaroon_df <- as.data.frame(Pro_RE22_MEmaroon)
Pro_RE22_MEdarkred_df <- as.data.frame(Pro_RE22_MEdarkred)
Pro_RE22_MEdarkgreen_df <- as.data.frame(Pro_RE22_MEdarkgreen)
Pro_RE22_MEyellow_df <- as.data.frame(Pro_RE22_MEyellow)
Pro_RE22_MEgrey_df <- as.data.frame(Pro_RE22_MEgrey)

colnames(Pro_RE22_MEmediumpurple3_df)[1] <- "ID"
colnames(Pro_RE22_MEblack_df)[1] <- "ID"
colnames(Pro_RE22_MEsteelblue_df)[1] <- "ID"
colnames(Pro_RE22_MEpurple_df)[1] <- "ID"
colnames(Pro_RE22_MEsalmon_df)[1] <- "ID"
colnames(Pro_RE22_MEmidnightblue_df)[1] <- "ID"
colnames(Pro_RE22_MEplum1_df)[1] <- "ID"
colnames(Pro_RE22_MEcoral1_df)[1] <- "ID"
colnames(Pro_RE22_MEsalmon4_df)[1] <- "ID"
colnames(Pro_RE22_MEblue_df)[1] <- "ID"
colnames(Pro_RE22_MEred_df)[1] <- "ID"
colnames(Pro_RE22_MEorange_df)[1] <- "ID"
colnames(Pro_RE22_MEdarkolivegreen_df)[1] <- "ID"
colnames(Pro_RE22_MEdarkorange_df)[1] <- "ID"
colnames(Pro_RE22_MEcyan_df)[1] <- "ID"
colnames(Pro_RE22_MEgreenyellow_df)[1] <- "ID"
colnames(Pro_RE22_MEmagenta_df)[1] <- "ID"
colnames(Pro_RE22_MEgreen_df)[1] <- "ID"
colnames(Pro_RE22_MElightgreen_df)[1] <- "ID"
colnames(Pro_RE22_MEsaddlebrown_df)[1] <- "ID"
colnames(Pro_RE22_MEdarkseagreen4_df)[1] <- "ID"
colnames(Pro_RE22_MEdarkolivegreen_df)[1] <- "ID"
colnames(Pro_RE22_MElightgreen_df)[1] <- "ID"
colnames(Pro_RE22_MElavenderblush3_df)[1] <- "ID"
colnames(Pro_RE22_MEnavajowhite2_df)[1] <- "ID"
colnames(Pro_RE22_MEpink_df)[1] <- "ID"
colnames(Pro_RE22_MElightpink4_df)[1] <- "ID"
colnames(Pro_RE22_MEcoral2_df)[1] <- "ID"
colnames(Pro_RE22_MEmaroon_df)[1] <- "ID"
colnames(Pro_RE22_MEdarkred_df)[1] <- "ID"
colnames(Pro_RE22_MEdarkgreen_df)[1] <- "ID"
colnames(Pro_RE22_MEyellow_df)[1] <- "ID"
colnames(Pro_RE22_MEgrey_df)[1] <- "ID"

Pro_RE22_MEmediumpurple3_apop <- merge(Pro_RE22_MEmediumpurple3_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEblack_apop <- merge(Pro_RE22_MEblack_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEsteelblue_apop <- merge(Pro_RE22_MEsteelblue_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEpurple_apop <- merge(Pro_RE22_MEpurple_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEsalmon_apop <- merge(Pro_RE22_MEsalmon_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEmidnightblue_apop <- merge(Pro_RE22_MEmidnightblue_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEplum1_apop <- merge(Pro_RE22_MEplum1_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEcoral1_apop <- merge(Pro_RE22_MEcoral1_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEsalmon4_apop <- merge(Pro_RE22_MEsalmon4_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEblue_apop <- merge(Pro_RE22_MEblue_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEred_apop <- merge(Pro_RE22_MEred_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEorange_apop <- merge(Pro_RE22_MEorange_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEdarkolivegreen_apop <- merge(Pro_RE22_MEdarkolivegreen_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEdarkorange_apop <- merge(Pro_RE22_MEdarkorange_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEcyan_apop <- merge(Pro_RE22_MEcyan_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEgreenyellow_apop <- merge(Pro_RE22_MEgreenyellow_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEmagenta_apop <- merge(Pro_RE22_MEmagenta_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEgreen_apop <- merge(Pro_RE22_MEgreen_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MElightgreen_apop <- merge(Pro_RE22_MElightgreen_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEsaddlebrown_apop <- merge(Pro_RE22_MEsaddlebrown_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEdarkseagreen4_apop <- merge(Pro_RE22_MEdarkseagreen4_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEdarkolivegreen_apop <- merge(Pro_RE22_MEdarkolivegreen_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MElightgreen_apop <- merge(Pro_RE22_MElightgreen_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MElavenderblush3_apop <- merge(Pro_RE22_MElavenderblush3_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEnavajowhite2_apop <- merge(Pro_RE22_MEnavajowhite2_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEpink_apop <- merge(Pro_RE22_MEpink_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MElightpink4_apop <- merge(Pro_RE22_MElightpink4_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEcoral2_apop <- merge(Pro_RE22_MEcoral2_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEmaroon_apop <- merge(Pro_RE22_MEmaroon_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEdarkred_apop <- merge(Pro_RE22_MEdarkred_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEdarkgreen_apop <- merge(Pro_RE22_MEdarkgreen_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEyellow_apop <- merge(Pro_RE22_MEyellow_df, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_MEgrey_apop <- merge(Pro_RE22_MEgrey_df, C_vir_rtracklayer_apop_product_final, by = "ID")

nrow(Pro_RE22_MEmediumpurple3_apop) # 4 Rho IFI21
nrow(Pro_RE22_MEblack_apop) # 12 TLR6, programmed cell death protein, STAT5A
nrow(Pro_RE22_MEsteelblue_apop) # 1 cell death-inducing p53-target protein 1
nrow(Pro_RE22_MEpurple_apop) # 4, growth arrest and DNA damage inducible protein GADD45, TGF beta, heat shock protein 
nrow(Pro_RE22_MEsalmon_apop) # 6, CD151, NFkB inhibition, cathepsin B, LITAF
nrow(Pro_RE22_MEmidnightblue_apop) # 2 FADD, calpain 5
nrow(Pro_RE22_MEplum1_apop) # 1 cathepsin L
nrow(Pro_RE22_MEcoral1_apop) # 0
nrow(Pro_RE22_MEsalmon4_apop) # 1 PCDP 10
nrow(Pro_RE22_MEblue_apop) # 23 TRAF, BIRIAP, GIMAP, calpain
nrow(Pro_RE22_MEred_apop) # 12, TNFL, IP3R, PCDP, IFI44, poly ADP polymerase, mitochondrial heat shock proteins
nrow(Pro_RE22_MEorange_apop) # 2 MAPK, IP3R
nrow(Pro_RE22_MEdarkolivegreen_apop) # 3 tollo
nrow(Pro_RE22_MEdarkorange_apop) # 3 caap1
nrow(Pro_RE22_MEcyan_apop) #7 MyD88, lymphotoxin alpha, caspase 7, tLr1, TRAF3 
nrow(Pro_RE22_MEgreenyellow_apop) # 3 IP3R, GADD45, TLR10
nrow(Pro_RE22_MEmagenta_apop) # bcl1, CREB, NFkB, rhoE, calpain 3
nrow(Pro_RE22_MEgreen_apop) # 15, TRAF BIR IAP, CRADD-like, IL17, NR13 like, CD151, TNFRSF16
nrow(Pro_RE22_MElightgreen_apop) # 6, TRAF3, JNK, caspase 8
nrow(Pro_RE22_MEsaddlebrown_apop) #3  GADD 45, TLR2, CREB
nrow(Pro_RE22_MEdarkseagreen4_apop) # 1 TNFAIP3
nrow(Pro_RE22_MElavenderblush3_apop) #0
nrow(Pro_RE22_MEnavajowhite2_apop) #0
nrow(Pro_RE22_MEpink_apop) #9 NFkB, calpain, BIRIAP, JNK, pyrin, UNC5C, MYD88
nrow(Pro_RE22_MElightpink4_apop) #3 STAT5A, P53 damage regulated protein 1, MAPK1
nrow(Pro_RE22_MEcoral2_apop) #0
nrow(Pro_RE22_MEmaroon_apop) # 1 cell death inducing p53 target
nrow(Pro_RE22_MEdarkred_apop) # 3, TLR6, BAG, IP3R
nrow(Pro_RE22_MEdarkgreen_apop) # 1 TLR3
nrow(Pro_RE22_MEyellow_apop) #13 caspase8, BIRIAP, aurora kinase, PDCP 6, TLR4, TLR2, ligase CHIP, caspase9  (definitely extrinsic pathway usage)
nrow(Pro_RE22_MEgrey_apop) # 51 GIMAP, BIRIAP, casp2,3, MAPK, TRAF, programmed cell death proteins , TLRs 

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

Pro_RE22_Module_hub_genes <- chooseTopHubInEachModule(
  Pro_RE22_dds_rlog_matrix_common, 
  Pro_RE22_colorh, 
  power = 2,  # power used for the adjacency network
  type = "signed hybrid", 
  corFnc = "bicor"
)
class(Pro_RE22_Module_hub_genes)
Pro_RE22_Module_hub_genes_df <- as.data.frame(Pro_RE22_Module_hub_genes)
colnames(Pro_RE22_Module_hub_genes_df)[1] <- "ID"
Pro_RE22_Module_hub_genes <- merge(Pro_RE22_Module_hub_genes_df, C_vir_rtracklayer, by = "ID")
nrow(Pro_RE22_Module_hub_genes) 
# includes programmed cell death protein 6 as a hub gene, apoptosis regulatory protein Siva 

#### MEASURE MODULE PRESERVATION BETWEEN C. VIRGINICA DIFFERENT EXPERIMENTS ####

# Column names agree between the datasets, so we can go ahead and measure module preservation between all 
## 1. MEASURE OVERLAP BETWEEN RE22 AND RODs

# Assess whether Pro_RE22 and ROD_Res modules are preserved
Pro_RE22_ROD_multiExpr = list(Pro_RE22=list(data=Pro_RE22_dds_rlog_matrix_common),ROD_Res=list(data=ROD_Resistant_dds_rlog_matrix_common))
Pro_RE22_multiColor = list(Pro_RE22 = Pro_RE22_moduleColors) 
Pro_RE22_mp =modulePreservation(Pro_RE22_ROD_multiExpr, Pro_RE22_multiColor,
                                referenceNetworks=1,
                                verbose=3,
                                networkType="signed hybrid", # use same signed hybrid as before
                                corFnc = "bicor", # use recommended bicor as before 
                                nPermutations=200,
                                randomSeed = 1, # recommended in langfelder tutorial
                                quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Pro_RE22_stats = Pro_RE22_mp$preservation$Z$ref.Pro_RE22$inColumnsAlsoPresentIn.ROD_Res
Pro_RE22_stats_order <- Pro_RE22_stats[order(-Pro_RE22_stats[,2]),c(1:2)]
Pro_RE22_stats_preserved <- Pro_RE22_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)
## Moderate preservation of between Pro_RE22 and ROD Res
# brown                 1000     7.8915583
# red                   1000     7.1531975
# steelblue              175     6.5939398
# black                 1000     6.3774745
# midnightblue           527     5.7505155
# purple                 731     5.7203439
# blue                  1000     5.1007921
# pink                   787     5.0031835

# Were any of these preserved modules significant in Pro_RE22 or ROD_Res?
Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list, "ME")
Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list, "ME")
Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Pro_RE22_moduleTraitCor_Pval_df_sig_list, "ME")
ROD_Res_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(ROD_Res_moduleTraitCor_Pval_df_sig_list, "ME")

Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm <- as.data.frame(Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm)
Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm <- as.data.frame(Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm)
Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm <- as.data.frame(Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm)
ROD_Res_moduleTraitCor_Pval_df_sig_list_rm <- as.data.frame(  ROD_Res_moduleTraitCor_Pval_df_sig_list_rm)

colnames(Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm)[1] <- "mod_name"
colnames(Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm)[1] <- "mod_name"
colnames(Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm)[1] <- "mod_name"
colnames(ROD_Res_moduleTraitCor_Pval_df_sig_list_rm)  [1] <- "mod_name"

Pro_RE22_RI_SIG_ROD_Res_preserved <-merge(Pro_RE22_stats_preserved, Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm)
# 6 : black,  blue, midnightblue ,purple,red ,steelblue 
Pro_RE22_S4_SIG_ROD_Res_preserved <-merge(Pro_RE22_stats_preserved, Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm)
# 0
Pro_RE22_RE22_SIG_ROD_Res_preserved <-merge(Pro_RE22_stats_preserved, Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm)
# pink
Pro_RE22_ROD_RES_SIG_preserved <-merge(Pro_RE22_stats_preserved, ROD_Res_moduleTraitCor_Pval_df_sig_list_rm)
# 0 

# Assess whether Pro_RE22 and ROD_Sus modules are preserved
Pro_RE22_Sus_multiExpr = list(Pro_RE22=list(data=Pro_RE22_dds_rlog_matrix_common),ROD_Sus=list(data=ROD_Susceptible_dds_rlog_matrix_common))
Pro_RE22_Sus_mp =modulePreservation(Pro_RE22_Sus_multiExpr, Pro_RE22_multiColor,
                                    referenceNetworks=1,
                                    verbose=3,
                                    networkType="signed hybrid", # use same signed hybrid as before
                                    corFnc = "bicor", # use recommended bicor as before 
                                    nPermutations=100,
                                    randomSeed = 1, # recommended in langfelder tutorial
                                    quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Pro_RE22_Sus_stats = Pro_RE22_Sus_mp$preservation$Z$ref.Pro_RE22$inColumnsAlsoPresentIn.ROD_Sus 
Pro_RE22_Sus_stats_order <- Pro_RE22_Sus_stats[order(-Pro_RE22_Sus_stats[,2]),c(1:2)]
Pro_RE22_Sus_stats_preserved <- Pro_RE22_Sus_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)

# Were any of these preserved modules significant in Pro_RE22 or ROD_Sus?
ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(ROD_Sus_moduleTraitCor_Pval_df_sig_list, "ME")
ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm <- as.data.frame(  ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm)
colnames(ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm)  [1] <- "mod_name"

Pro_RE22_RI_SIG_ROD_Sus_preserved <-merge(Pro_RE22_Sus_stats_preserved, Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm)
# blue, red
Pro_RE22_S4_SIG_ROD_Sus_preserved <-merge(Pro_RE22_Sus_stats_preserved, Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm)
# 0
Pro_RE22_RE22_SIG_ROD_Sus_preserved <-merge(Pro_RE22_Sus_stats_preserved, Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm)
# 0
Pro_RE22_ROD_SUS_SIG_preserved <-merge(Pro_RE22_Sus_stats_preserved, ROD_Sus_moduleTraitCor_Pval_df_sig_list_rm)
# brown

# Assess whether Pro_RE22 and Dermo_Tol modules are preserved
Pro_RE22_Dermo_multiExpr = list(Pro_RE22=list(data=Pro_RE22_dds_rlog_matrix_common),Dermo_Tol=list(data=Dermo_Tolerant_dds_vst_matrix_common))
Pro_RE22_multiColor = list(Pro_RE22 = Pro_RE22_moduleColors) 
Pro_RE22_Dermo_mp =modulePreservation(Pro_RE22_Dermo_multiExpr, Pro_RE22_multiColor,
                                      referenceNetworks=1,
                                      verbose=3,
                                      networkType="signed hybrid", # use same signed hybrid as before
                                      corFnc = "bicor", # use recommended bicor as before 
                                      nPermutations=100,
                                      randomSeed = 1, # recommended in langfelder tutorial
                                      quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Pro_RE22_Dermo_stats = Pro_RE22_Dermo_mp$preservation$Z$ref.Pro_RE22$inColumnsAlsoPresentIn.Dermo_Tol
Pro_RE22_Dermo_stats_order <- Pro_RE22_Dermo_stats[order(-Pro_RE22_Dermo_stats[,2]),c(1:2)]
Pro_RE22_Dermo_stats_preserved <- Pro_RE22_Dermo_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)

# Were any of these preserved modules significant in Pro_RE22 or Dermo_Tol?
Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Dermo_Tol_moduleTraitCor_Pval_df_sig_list, "ME")
Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm <- as.data.frame(  Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm)
colnames(Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm)  [1] <- "mod_name"

Pro_RE22_RI_SIG_Dermo_Tol_preserved <-merge(Pro_RE22_Dermo_stats_preserved, Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm)
#1   black       1000      5.354414
#2         blue       1000      6.243046
#3 midnightblue        527      5.580933
#4       purple        731      9.399914
#5          red       1000      8.038084
#6       salmon        615     10.598431
#7    steelblue        175      9.376713
Pro_RE22_S4_SIG_Dermo_Tol_preserved <-merge(Pro_RE22_Dermo_stats_preserved, Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm)
# green
Pro_RE22_RE22_SIG_RDermo_Tol_preserved <-merge(Pro_RE22_Dermo_stats_preserved, Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm)
# green yellow
Pro_RE22_Dermo_Tol_SIG_preserved <-merge(Pro_RE22_Dermo_stats_preserved, Dermo_Tol_moduleTraitCor_Pval_df_sig_list_rm)
# turqoise yellow

# Assess whether Pro_RE22 and Dermo_Sus modules are preserved
Pro_RE22_Dermo_Sus_multiExpr = list(Pro_RE22=list(data=Pro_RE22_dds_rlog_matrix_common),Dermo_Sus=list(data=Dermo_Susceptible_dds_vst_matrix_common))
Pro_RE22_multiColor = list(Pro_RE22 = Pro_RE22_moduleColors) 
Pro_RE22_Dermo_Sus_mp =modulePreservation(Pro_RE22_Dermo_Sus_multiExpr, Pro_RE22_multiColor,
                                          referenceNetworks=1,
                                          verbose=3,
                                          networkType="signed hybrid", # use same signed hybrid as before
                                          corFnc = "bicor", # use recommended bicor as before 
                                          nPermutations=100,
                                          randomSeed = 1, # recommended in langfelder tutorial
                                          quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
Pro_RE22_Dermo_Sus_stats = Pro_RE22_Dermo_Sus_mp$preservation$Z$ref.Pro_RE22$inColumnsAlsoPresentIn.Dermo_Sus
Pro_RE22_Dermo_Sus_stats_order <- Pro_RE22_Dermo_Sus_stats[order(-Pro_RE22_Dermo_Sus_stats[,2]),c(1:2)]
Pro_RE22_Dermo_Sus_stats_preserved <- Pro_RE22_Dermo_Sus_stats_order %>% rownames_to_column("mod_name") %>% filter(Zsummary.pres > 5.0)

# Were any of these preserved modules significant in Pro_RE22 or Dermo_Sus?
Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm <- str_remove(Dermo_Sus_moduleTraitCor_Pval_df_sig_list, "ME")
Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm <- as.data.frame(  Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm)
colnames(Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm)  [1] <- "mod_name"

Pro_RE22_RI_SIG_Dermo_Sus_preserved <-merge(Pro_RE22_Dermo_Sus_stats_preserved, Pro_RE22_RI_moduleTraitCor_Pval_df_sig_list_rm)
# mod_name moduleSize Zsummary.pres
# 1      blue       1000      6.176316
# 2    purple        731      6.821283
# 3       red       1000      6.060967
# 4    salmon        615      9.461935
# 5 steelblue        175      7.608015
Pro_RE22_S4_SIG_Dermo_Sus_preserved <-merge(Pro_RE22_Dermo_Sus_stats_preserved, Pro_RE22_S4_moduleTraitCor_Pval_df_sig_list_rm)
# 0
Pro_RE22_RE22_SIG_Dermo_Sus_preserved <-merge(Pro_RE22_Dermo_Sus_stats_preserved, Pro_RE22_moduleTraitCor_Pval_df_sig_list_rm)
# 0 
Pro_RE22_Dermo_Sus_SIG_preserved <-merge(Pro_RE22_Dermo_Sus_stats_preserved, Dermo_Sus_moduleTraitCor_Pval_df_sig_list_rm)
# 0 

# Assess whether Dermo_Tol and ROD_Res modules are preserved
ROD_Dermo_multiExpr = list(Dermo_Tol=list(data=Dermo_Tolerant_dds_vst_matrix_common), ROD_Res=list(data=ROD_Resistant_dds_rlog_matrix_common))
ROD_Dermo_multiColor = list(Dermo_Tol = Dermo_Tol_moduleColors)
ROD_Dermo_mp =modulePreservation(ROD_Dermo_multiExpr, ROD_Dermo_multiColor,
                                 referenceNetworks=1,
                                 verbose=3,
                                 networkType="signed hybrid", # use same signed hybrid as before
                                 corFnc = "bicor", # use recommended bicor as before 
                                 nPermutations=100,
                                 randomSeed = 1, # recommended in langfelder tutorial
                                 quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
ROD_Dermo_stats = ROD_Dermo_mp$preservation$Z$ref.$inColumnsAlsoPresentIn.A2 
ROD_Dermo_stats_order <- ROD_Dermo_stats[order(-ROD_Dermo_stats[,2]),c(1:2)]

ROD_Dermo_Sus_multiExpr = list(Dermo_Sus=list(data=Dermo_Susceptible_dds_vst_matrix_common), ROD_Sus=list(data=ROD_Susceptible_dds_rlog_matrix_common))
ROD_Dermo_Sus_multiColor = list(Dermo_Sus = Dermo_Sus_moduleColors)
ROD_Dermo_Sus_mp =modulePreservation(ROD_Dermo_Sus_multiExpr, ROD_Dermo_Sus_multiColor,
                                     referenceNetworks=1,
                                     verbose=3,
                                     networkType="signed hybrid", # use same signed hybrid as before
                                     corFnc = "bicor", # use recommended bicor as before 
                                     nPermutations=100,
                                     randomSeed = 1, # recommended in langfelder tutorial
                                     quickCor = 0) # number between 0 and 1 specifying the handling of missing data in calculation of correlation. Zero means exact but potentially slower calculations
ROD_Dermo_Sus_stats = ROD_Dermo_Sus_mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2 
ROD_Dermo_Sus_stats[order(-ROD_Dermo_Sus_stats[,2]),c(1:2)]



## Perform Functional Enrichment of Apoptosis genes using AnRichment 

# Tutorials: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/

