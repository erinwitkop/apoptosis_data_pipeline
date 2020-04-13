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
Dermo_Tol_moduleTraitCor_Pval_df_sig$mod_names

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

Dermo_Sus_moduleTraitCor_Pval_df_sig$mod_names
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

ROD_Sus_moduleTraitCor_Pval_df_sig$mod_names
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

### 1. Are the same module identifiers matched to the same genes across experiments?

Pro_RE22_cyan <- colnames(Pro_RE22_dds_rlog_matrix_common)[Pro_RE22_moduleColors == "cyan"       ]
ROD_Res_cyan <- colnames(ROD_Resistant_dds_rlog_matrix_common)[ROD_Res_moduleColors == "cyan"       ]
ROD_Sus_cyan <- colnames(ROD_Susceptible_dds_rlog_matrix_common)[ROD_Sus_moduleColors == "cyan"       ]

head(Pro_RE22_cyan)

all(head(Pro_RE22_cyan) %in% ROD_Res_cyan) # FALSE, seems like the same genes were not assigned to same modules 

## Which modules match one another?
# Load matchModules function 
matchModules <- function (gn1, mod1, gn2, mod2, omit="grey", allColors = standardColors()){
  ## This function converts the module colors from network 2 to the module colors
  ##  of best match in network 1, given the gene names (gn1,gn2) and corresponding
  ##  module assignments (mod1,mod2).  Non-overlapping modules will have unique labels.
  ## omit      = colors that should not be changed
  ## allColors = color set to choose from (default is all standardColors)
  
  # Write out data from network 2 then read it back in to get module overlaps
  out2 = cbind(gn2, mod2);
  colnames(out2) = c("Gene","Var1")
  write.csv(out2,"eraseMe.csv",row.names=FALSE)
  overlaps = userListEnrichment(gn1, mod1, "eraseMe.csv", "X", "eraseMe.csv")
  c1 = as.character(overlaps$sigOverlaps[,1])
  c2 = as.character(overlaps$sigOverlaps[,2]);
  c2 = substr(c2,3,nchar(c2))
  kp = (!is.element(c1,omit))&(!is.element(c2,omit))
  c1 = c1[kp];  c2 = c2[kp]
  
  # Change the labels in network 2 to the best overlapping module in network 1
  cOld <- cNew <- changed <- sort(unique(mod2))
  cUnused = setdiff(allColors,union(mod1,mod2))
  while(length(c1)>0){
    cNew[cOld==c2[1]]=c1[1]
    changed[cOld==c2[1]]="YES"
    kp = (c1!=c1[1])&(c2!=c2[1])
    c1 = c1[kp];  c2 = c2[kp]
  }
  changed[is.element(cOld,omit)] = "YES"
  cNew[changed!="YES"] = cUnused[1:sum(changed!="YES")]
  modOut = mod2
  for (i in 1:length(cNew)) modOut[mod2==cOld[i]] = cNew[i]
  write(paste("Old - New labels:",cOld,"-",cNew),"")
  return(modOut)
}

Pro_RE22_gene <- rownames(Pro_RE22_dds_rlog_matrix_common)
ROD_Res_gene <- rownames(ROD_Resistant_dds_rlog_matrix_common)
ROD_Sus_gene <- rownames(ROD_Susceptible_dds_rlog_matrix_common)


modulesB2_new = matchModules(Pro_RE22_gene, Pro_RE22_moduleColors, ROD_Res_gene, ROD_Res_moduleColors)
# Old - New labels: black – purple # Your screen output will look like this
# Old - New labels: blue - brown
# Old - New labels: brown – greenyellow # etc.


enrichmentsB2A1 = userListEnrichment(GeneAB,modulesB2,"kMEtable1.csv","A1","enrichmentB2_A1.csv") enrichmentsB2A1$sigOverlaps # To show the significant overlaps on the screen
# InputCategories UserDefinedCategories CorrectedPvalues
# pOut "blue"
# pOut "turquoise"
# pOut "green"
"A1_brown" "A1_yellow" "A1_brown"
"2.41636641189253e-59" "1.51940862185272e-43" "3.67833377321179e-41" # etc.



## 1. MEASURE OVERLAP BETWEEN RE22 AND RODs

# Assess whether RE22 moduldes are preserved
RE22_ROD_multiExpr = list(Pro_RE22=list(data=t(Pro_RE22_dds_rlog_matrix_common)),ROD_Res=list(data=t(ROD_Resistant_coldata_collapsed_binarize)),
                          ROD_Res=list(data=t(ROD_Susceptible_coldata_collapsed_binarize)))
multiColor = list(Pro_RE22 = Pro_RE22_modules) 
mp=modulePreservation(multiExpr,multiColor,
                      referenceNetworks=1,
                      verbose=3,networkType="signed",
                      nPermutations=30
                      ,maxGoldModuleSize=100,
                      maxModuleSize=20000) 
stats = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2 
stats[order(-stats[,2]),c(1:2)]



## IDENTIFY MODULES ENRICHED FOR APOPTOSIS USING anRichment







# Run WGCNA on the datasets with large number of genes 
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-blockwise.pdf
# Choose a set of soft-thresholding powers, # Call the network topology analysis function

# Variation Two
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(Dermo_dds_family_WGCNA_vsd_matrix_transpose,corFnc = "bicor", corOptions=list(maxPOutliers=0.1),networkType = "signed", powerVector = powers, verbose = 5 )
# this setting would suggest power of 14 for the data 
# the signed function is the recommended one
# the authors suggest the biweight mid-correlation as a robust alternative, the "bicor" corFnc
#Bicor note: can produce unwanted results when the data have a bi-modal distribution (e.g., when a gene expression depends heavily on a binary variable 
# such as disease status or genotype) or when one of the variables entering the correlation is itself binary (or ordinal).
# For this reason, we strongly recommend using the argument maxPOutliers = 0.05 or 0.10 whenever the biweight midcorrelation is used. 
# This argument essentially forces bicor to never regard more than the specified proportion of samples as outliers.

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=0.90,col="red")

### Selecting 12 as the softthresholding power

# Perform blockwise clustering
bwnet = blockwiseModules(Dermo_dds_family_WGCNA, maxBlockSize = 5000,
                         power = 12, TOMType = "signed", minModuleSize = 30, # medium sensitivity 
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         corType = "bicor", # suggested by authors
                         maxPOutliers = 0.05,  # suggested in case of disease status examples
                         robustY = FALSE, # when dealing with binary variable, turns off robust treatment
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "Dermo_family_TOM",
                         verbose = 3)


# View modules,  bwnet$colors contains the module assignment, and bwnet$MEs contains the module eigengenes of the modules.
table(bwnet$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(bwnet$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(bwnet$dendrograms[[1]], mergedColors[bwnet$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# save module assignment and module eigengene information necessary for subsequent analysis
moduleLabels = bwnet$colors
moduleColors = labels2colors(bwnet$colors)
MEs = bwnet$MEs;
geneTree = bwnet$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "Dermo_networkConstruction.RData")

## Quantifying module-trait associations

# Load the expression and trait data saved in the first part
lnames = load(file = "Dermo_WGCNA.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Dermo_networkConstruction.RData")
lnames

# Correlate eigengenes with external traits for most significant associations
# Define numbers of genes and samples
nGenes = ncol(Dermo_dds_family_WGCNA)
nSamples = nrow(Dermo_dds_family_WGCNA)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(Dermo_dds_family_WGCNA, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, Dermo_coldata_binary, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Color code each association by the correlation value
sizeGrWindow(15,10)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(Dermo_coldata_binary),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               cex.lab = 0.7,
               zlim = c(-1,1), 
               yColorWidth = 0.2, 
               main = paste("Module-trait relationships"))

## Gene significance and module membership FAMILY 
# Gene significance is the correlation between gene and trait of interest 
# Als define a quantitative measure of module membership 
# MM as the correlation of the module eigengene and the gene expression profile.
# This allows us to quantify the similarity of all genes on the array 
# to every module.

# Define variable FamCode containing the Family column of Dermo_coldata_binary
FamCode = as.data.frame(Dermo_coldata_binary$FamCode)
names(FamCode) = "FamCode"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(Dermo_dds_family_WGCNA, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(Dermo_dds_family_WGCNA, FamCode, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(FamCode), sep="")
names(GSPvalue) = paste("p.GS.", names(FamCode), sep="")

# Identifying genes with high GS and MM
# MElightcyan is the module with the greatest positive correlation with Family
# MEturqoise is the module with greatest negative correlation
module="lightcyan"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")
# significant positive association between FamCode and this module 

# Summary output of network analysis results 
# Merge statistical info with gene annotation and write to file
colnames(Dermo_dds_family_WGCNA)
FamCode_lightcyan <- colnames(Dermo_dds_family_WGCNA)[moduleColors=="lightcyan"]
class(FamCode_lightcyan)
FamCode_lightcyan  <- as.data.frame(FamCode_lightcyan)
head(FamCode_lightcyan)
colnames(FamCode_lightcyan)[1] <- "transcript_id"
head(FamCode_lightcyan)

# Annotate 
FamCode_lightcyan_annot <- FamCode_lightcyan %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product", "gene","unique_go")], by = "transcript_id")
nrow(FamCode_lightcyan_annot) #77 
View(FamCode_lightcyan_annot)
FamCode_lightcyan_annot_apop <- FamCode_lightcyan_annot[grepl(paste(Apoptosis_names,collapse="|"), 
                                                              FamCode_lightcyan_annot$product, ignore.case = TRUE),]

## Gene significance and module membership FAMILY MODULE 2 
# Identifying genes with high GS and MM
# MEpaleturqoise is the module with the greatest positive correlation with Family

module="paleturquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")
# Not significant 

## Gene significance and module membership TREATMENT
# Gene significance is the correlation between gene and trait of interest 
# Als define a quantitative measure of module membership 
# MM as the correlation of the module eigengene and the gene expression profile.
# This allows us to quantify the similarity of all genes on the array 
# to every module.

# Define variable FamCode containing the Family column of Dermo_coldata_binary
Treat = as.data.frame(Dermo_coldata_binary$Treat)
names(Treat) = "Treat"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(Dermo_dds_family_WGCNA, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(Dermo_dds_family_WGCNA, Treat, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(Treat), sep="")
names(GSPvalue) = paste("p.GS.", names(Treat), sep="")

#### Identifying genes with high GS and MM ###
# MElightcyan is the module with the greatest positive correlation with Family
# MEturqoise is the module with greatest negative correlation

# Interpreting the heatmap https://support.bioconductor.org/p/111449/
# The module-trait heatmap usually represents the correlations of the module eigengenes with traits. 
# When that correlation is high, it means the eigengene increases with increasing trait.
# In a signed network (where all genes in a module are positively correlated with the eigengene)
# it will mean that (again if the eigengene-trait correlation is high) pretty much all genes 
# should also follow the same pattern of increasing expression with increasing trait values. 
# In an unsigned network you may also have genes that have the opposite behaviour since in an 
# unsigned network a module can contain also genes strongly negatively correlated with the eigengene.
# The eigengene-trait correlation measures the strength and direction of association between the
# module (more precisely, the representative profile) and the trait.
# If this is positive (negative), it means the trait increases (decreases) with increasing eigengene "expression".
# If this correlation is strong and the network is signed (or signed hybrid) 
# it means that most of the genes in the module will also exhibit a correlation with the trait of the same sign as the eigengene.
# In an unsigned network, the gene-trait correlations can have the same or opposite sign.

module="orange"
column = match(module, modNames)
moduleGenes = moduleColors==module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")

# association not significant 

### Create data frame holding information about all probes ###
probes=colnames(Dermo_dds_family_WGCNA)
probes2annot = match(probes, C_vir_rtracklayer_transcripts_GO$transcript_id)
geneTraitSignificance = as.data.frame(cor(Dermo_dds_family_WGCNA, Treat, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(Treat), sep="")
names(GSPvalue) = paste("p.GS.", names(Treat), sep="")

### Combine with module color, gene significance for weight, and module membership and p-values in all modules ###
# Create the starting data frame
geneInfo0 = data.frame(transcript_id = probes, 
                       product = C_vir_rtracklayer_transcripts_GO$product[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for FamCode
modOrder = order(-abs(cor(MEs, FamCode, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Treat));
geneInfo = geneInfo0[geneOrder, ]

# Write to a file
write.csv(geneInfo, file = "Dermo_geneInfo_WGCNA_matrix.csv") # this value has the significance by Treatment

# Write table with gene significance for Family 
probes=colnames(Dermo_dds_family_WGCNA)
probes2annot = match(probes, C_vir_rtracklayer_transcripts_GO$transcript_id)
geneFamCodeSignificance = as.data.frame(cor(Dermo_dds_family_WGCNA, FamCode, use = "p"))
GSPvalue_FamCode = as.data.frame(corPvalueStudent(as.matrix(geneFamCodeSignificance), nSamples));
names(geneFamCodeSignificance) = paste("GS.", names(FamCode), sep="")
names(GSPvalue_FamCode) = paste("p.GS.", names(FamCode), sep="")

# Create the starting data frame
geneInfo1 = data.frame(transcript_id = probes, 
                       product = C_vir_rtracklayer_transcripts_GO$product[probes2annot],
                       moduleColor = moduleColors,
                       geneFamCodeSignificance,
                       GSPvalue_FamCode)
# Order modules by their significance for FamCode
modOrder = order(-abs(cor(MEs, FamCode, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo1)
  geneInfo1 = data.frame(geneInfo1, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo1) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder1 = order(geneInfo1$moduleColor, -abs(geneInfo1$GS.FamCode));
geneInfo1 = geneInfo1[geneOrder1, ]

# Write to a file
write.csv(geneInfo1, file = "Dermo_geneInfo_WGCNA_matrix_FamCode.csv") # this value has the significance by FamCode

## Create combined gene significance for Treatment and FamCode for each transcript for plotting
gene_info_combined_significance  <- geneInfo0[,c(1:5)] %>% left_join(geneInfo1[,c("transcript_id", "moduleColor","product","GS.FamCode",
                                                                                  "p.GS.FamCode")], by = c("transcript_id","moduleColor",
                                                                                                           "product"))
head(gene_info_combined_significance)

# Subset significant genes
gene_info_combined_significance_Treat <- gene_info_combined_significance %>% filter(p.GS.Treat <= 0.05)
gene_info_combined_significance_FamCode <- gene_info_combined_significance %>% filter(p.GS.FamCode <= 0.05)

nrow(gene_info_combined_significance_Treat)
gene_info_combined_significance_Treat_apop <-gene_info_combined_significance_Treat[grepl(paste(Apoptosis_names,collapse="|"), 
                                                                                         gene_info_combined_significance_Treat$product, ignore.case = TRUE),]
gene_info_combined_significance_FamCode_apop <-gene_info_combined_significance_FamCode[grepl(paste(Apoptosis_names,collapse="|"), 
                                                                                             gene_info_combined_significance_FamCode$product, ignore.case = TRUE),]
nrow(gene_info_combined_significance_Treat_apop )
nrow(gene_info_combined_significance_FamCode_apop)

# plot gene significance of significant genes
ggplot(gene_info_combined_significance_Treat_apop, aes(x=product, y=GS.Treat)) + geom_col(position="dodge") + coord_flip()

gene_info_combined_significance_Treat_apop_plot <- ggplot(gene_info_combined_significance_Treat_apop,
                                                          aes(x=product, y=GS.Treat, fill=GS.Treat)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Gene Significance for Treatment") +
  ylab("WGCNA Gene Significance") + theme(axis.text.x = element_text(size=1) ) 

gene_info_combined_significance_FamCode_apop_plot <- ggplot(gene_info_combined_significance_FamCode_apop,
                                                            aes(x=product, y=GS.FamCode, fill=GS.FamCode)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Gene Significance for Treatment") +
  ylab("WGCNA Gene Significance") + theme(axis.text.x = element_text(size=1) ) 

## Perform Functional Enrichment of Apoptosis genes using AnRichment 

# Tutorials: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/

