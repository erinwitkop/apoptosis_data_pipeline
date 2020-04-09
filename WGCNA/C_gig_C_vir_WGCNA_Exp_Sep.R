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
# source("https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/installAnRichment.R"); installAnRichment(); 


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

# DATA FORMATTING, batch effect removal only necessary for non-biological technical aspects 
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

## SELECT SOFT THRESHOLDING POWER FOR EACH EXPERIMENT ##
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function, following general recommendations to set network type to "signed hybrid" and using the "bicor" correlation
Dermo_Tolerant_dds_vst_matrix_sft <- pickSoftThreshold(Dermo_Tolerant_dds_vst_matrix, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
Dermo_Susceptible_dds_vst_matrix_sft <- pickSoftThreshold(Dermo_Susceptible_dds_vst_matrix, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor")
Probiotic_dds_rlog_matrix_sft <- pickSoftThreshold(Probiotic_dds_rlog_matrix, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
# Warning message:
# executing %dopar% sequentially: no parallel backend registered 
#From Peter Langfelder: https://bioinformatics.stackexchange.com/questions/10555/r-wgcna-error-code:
#What you see is a warning, not an error. 
#Your calculation will run fine, just slower. Unless you see other errors, you should be able to complete all steps of the analysis.

ROD_Resistant_dds_rlog_matrix_sft <- pickSoftThreshold(ROD_Resistant_dds_rlog_matrix, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
ROD_Susceptible_dds_rlog_matrix_sft <- pickSoftThreshold(ROD_Susceptible_dds_rlog_matrix, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 
Pro_RE22_dds_rlog_matrix_sft <- pickSoftThreshold(Pro_RE22_dds_rlog_matrix, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "bicor") 

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

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
# Softthreshold of 2 

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

# Pro_RE22
# Scale-free topology fit index as a function of the soft-thresholding power
plot(Pro_RE22_dds_rlog_matrix_sft$fitIndices[,1], -sign(Pro_RE22_dds_rlog_matrix_sft$fitIndices[,3])*Pro_RE22_dds_rlog_matrix_sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(Pro_RE22_dds_rlog_matrix_sft$fitIndices[,1], -sign(Pro_RE22_dds_rlog_matrix_sft$fitIndices[,3])*Pro_RE22_dds_rlog_matrix_sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(Pro_RE22_dds_rlog_matrix_sft$fitIndices[,1],Pro_RE22_dds_rlog_matrix_sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(Pro_RE22_dds_rlog_matrix_sft$fitIndices[,1], Pro_RE22_dds_rlog_matrix_sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Set soft thresholding to 4

### ONE STEP NETWORK CONSTRUCTION, MODULE DETECTION, MODULE DENDROGRAM INSPECTION
Dermo_Tol_net = blockwiseModules(Dermo_Tolerant_dds_vst_matrix, power = 3, # picked suitable power in the code above 
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
Dermo_Sus_net = blockwiseModules(Dermo_Susceptible_dds_vst_matrix, power = 2, # picked suitable power in the code above 
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
Probiotic_net = blockwiseModules(Probiotic_dds_rlog_matrix, power = 9, # picked because less than 20 samples and isn't scale free
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
plotDendroAndColors(Probiotic_net$dendrograms[[1]], Dermo_Sus_mergedColors[Probiotic_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
Probiotic_moduleLabels = Probiotic_net$colors
Probiotic_moduleColors = labels2colors(Probiotic_net$colors)
Probiotic_MEs = Probiotic_net$MEs
Probiotic_geneTree = Probiotic_net$dendrograms[[1]]

# ROD Res
ROD_Res_net = blockwiseModules(ROD_Resistant_dds_rlog_matrix, power = 9, # picked because less than 20 samples and isn't scale free
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
ROD_Sus_net = blockwiseModules(ROD_Sus_dds_rlog_matrix, power = 9, # picked because less than 20 samples and isn't scale free
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
Pro_RE22_net = blockwiseModules(ROD_Sus_dds_rlog_matrix, power = 9, # picked because less than 20 samples and isn't scale free
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

## BINARIZE ALL CATEGORICAL VARIABLES TO TEST DISEASE CHALLENGE ASSOCIATIONS


## TEST ASSOCIATIONS OF MODULES WITH DISEASE CHALLENGE AND ISOLATE SIGNIFICANT MODULES 


## IDENTIFY MODULES ENRICHED FOR APOPTOSIS USING anRichment

## MEASURE MODULE PRESERVATION BETWEEN DIFFERENT EXPERIMENTS 

### Run WGCNA first to determine differences between treatments (next will do the consensus set)

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

