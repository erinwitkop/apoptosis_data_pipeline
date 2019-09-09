# C_vir_Dermo_Analysis_2015_Script

# Load packages

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("apeglm")

library(DESeq2)  
library(ggplot2)
library(magrittr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(questionr)
library(apeglm)
library(genefilter)
library(fission)
library(tidyr)
library(stringr)
library(rtracklayer)
library(UpSetR)
library(reshape2)
library(plyr)
library(Repitools)
library(purrr)
library(tibble)
  
### Versions
# R version 3.6.1
# Deseq2 DESeq2_1.24.0.tar.gz

## Helpful, recently updataed resources
# Updated 2018 Vignette of DESeq2 : https://bioconductor.org/packages/3.7/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
# Updated Workflow for DESeq2: https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#the-deseqdataset-object-sample-information-and-the-design-formula
# Deseq2 Manual updated July 18th, 2019: https://www.bioconductor.org/packages/devel/bioc/manuals/DESeq2/man/DESeq2.pdf
  # the results section of this is particularly useful

##### Data info #####

# The DA famiy is currently the reference while the 6h timepoint is the reference timepoint
# DA is susceptible and LB is tolerant
# The rownames for the coldata are the library prep IDs (not included in an actual row of the table)
# Transcriptomes provided by Dina Proestou (USDA ARS) from 2015 challenge
  
######## Preparing Data ####

#Load in gene expression data as count matrix
Dermo_counts <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/RAW DATA AND INFO/2015 Dermo Challenge /DATA/6h_36h_7d_DA_LB_countMatrix.csv", 
                                   row.names="X")
head(Dermo_counts)

#Load in sample metadata
Dermo_coldata <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/RAW DATA AND INFO/2015 Dermo Challenge /DATA/6h_36h_7d_DA_LB_metadata.csv",row.names=1 )
head(Dermo_coldata)  
nrow(Dermo_coldata) 

# Make sure the columns of the count matrix and rows of the column data (sample metadata)
# are in the same order. Both of the following should return true

all(rownames(Dermo_coldata) %in% colnames(Dermo_counts))  #Should return TRUE
# returns TRUE
all(rownames(Dermo_coldata) == colnames(Dermo_counts))    # should return TRUE
#returns TRUE

# add column that is specifically rownames
Dermo_coldata$rownames <- rownames(Dermo_coldata)
head(Dermo_coldata)


##### Splitting up the data for each different comparison ####

# 1. Comparison of LB and DA at each time point (three data frames, one for each time point)
Dermo_coldata_6h <- Dermo_coldata %>% filter(Timepoint == "6h")
Dermo_coldata_36h <- Dermo_coldata %>% filter(Timepoint == "36h")
Dermo_coldata_7d <- Dermo_coldata %>% filter(Timepoint == "7d")

#  for each of the groups above and use that list to subset the Counts matrix
SampleID_6h <- Dermo_coldata_6h$rownames
SampleID_36h <- Dermo_coldata_36h$rownames
SampleID_7d <- Dermo_coldata_7d$rownames
  
Dermo_counts_6h <- Dermo_counts[,names(Dermo_counts) %in% SampleID_6h]
head(Dermo_counts_6h)
Dermo_counts_36h <- Dermo_counts[,names(Dermo_counts) %in% SampleID_36h]
Dermo_counts_7d <- Dermo_counts[,names(Dermo_counts) %in% SampleID_7d]

# 2. Comparison of LB to itself across timepoints, and DA to itself across time points

Dermo_coldata_LB_all <- Dermo_coldata %>% filter(FamCode == "LB")
Dermo_coldata_DA_all <- Dermo_coldata %>% filter(FamCode == "DA")

# extract list of samples
SampleID_LB <- Dermo_coldata_LB_all$rownames
SampleID_DA <- Dermo_coldata_DA_all$rownames

# subset counts matrix
Dermo_counts_LB_all <- Dermo_counts[,names(Dermo_counts) %in% SampleID_LB]
Dermo_counts_DA_all <- Dermo_counts[,names(Dermo_counts) %in% SampleID_DA]

# 3. To compare the 36hr to 7d, I need to get rid of the 6h group

Dermo_coldata_LB_middle_late <- Dermo_coldata_LB_all %>% filter(Timepoint !="6h")
Dermo_coldata_DA_middle_late <- Dermo_coldata_DA_all %>% filter(Timepoint !="6h")

# extract list of samples
SampleID_LB_middle_late <- Dermo_coldata_LB_middle_late$rownames
SampleID_DA_middle_late <- Dermo_coldata_DA_middle_late$rownames

# subset counts matrix
Dermo_counts_LB_middle_late <- Dermo_counts[,names(Dermo_counts) %in% SampleID_LB_middle_late]
Dermo_counts_DA_middle_late <- Dermo_counts[,names(Dermo_counts) %in% SampleID_DA_middle_late]



####### Check levels #######
# It is prefered in R that the first level of a factor be the reference level for comparison
# (e.g. control, or untreated samples), so we can relevel the dex factor like so

# Check factor levels, set it so that comparison group is the first

# LB all timepoints
levels(Dermo_coldata_LB_all$Timepoint) # #"6h"  "36h" "7d"
Dermo_coldata_LB_all$Timepoint <- factor(Dermo_coldata_LB_all$Timepoint, levels=c("6h","36h","7d"))

# DA all timepoints
levels(Dermo_coldata_DA_all$Timepoint) # #"6h"  "36h" "7d" 
Dermo_coldata_DA_all$Timepoint <- factor(Dermo_coldata_DA_all$Timepoint, levels=c("6h","36h","7d"))

# DA and LB only 36 and 7d
levels(Dermo_coldata_LB_middle_late$Timepoint) # "6h"  "36h" "7d" , need to droplevels
Dermo_coldata_LB_middle_late$Timepoint <- droplevels(Dermo_coldata_LB_middle_late$Timepoint) # "36h" "7d"
levels(Dermo_coldata_LB_middle_late$Timepoint)
levels(Dermo_coldata_DA_middle_late$Timepoint) # "6h"  "36h" "7d" , need to droplevels
Dermo_coldata_DA_middle_late$Timepoint <- droplevels(Dermo_coldata_DA_middle_late$Timepoint) # "36h" "7d"
levels(Dermo_coldata_DA_middle_late$Timepoint)# "36h" "7d"

# 6h both families
levels(Dermo_coldata_6h$FamCode)
Dermo_coldata_6h$FamCode <- factor(Dermo_coldata_6h$FamCode, levels=c("DA","LB"))
# "DA" "LB"

# 36h both families
levels(Dermo_coldata_36h$FamCode)
Dermo_coldata_36h$FamCode <- factor(Dermo_coldata_36h$FamCode, levels=c("DA","LB"))
# "DA" "LB"

# 7d both families
levels(Dermo_coldata_7d$FamCode)
Dermo_coldata_7d$FamCode <- factor(Dermo_coldata_7d$FamCode, levels=c("DA","LB"))
# "DA" "LB"

#### Build DESeqDataSetFromMatrix ####

#This object specifies the count data and metadata you will work with. The design piece is critical.
#For this, I will include the library prep date because Mary said there are significant batch effects present.
#Unlike their analysis, I will set my second design variable as the Family code as it first
# differences between families that I care about. 

# dds object with all data combined for comparing all sample clustering in exploratory analysis

Dermo_dds_family <- DESeqDataSetFromMatrix(countData = Dermo_counts,
                                           colData = Dermo_coldata,
                                           design = ~Library_Prep_Date + Treat +FamCode)

# This timecourse comparison models the family difference at time zero, the difference through time, and family specific differences through time
Dermo_dds_timecourse_family <- DESeqDataSetFromMatrix(countData = Dermo_counts,
                                                      colData = Dermo_coldata,
                                                      design= ~Library_Prep_Date+FamCode + Timepoint + FamCode:Timepoint + Treat)

# This dds object will allow me to compare both families at each time point separately
  # 6h data all had the same library prep date, so removed this from the formula
Dermo_dds_6h <- DESeqDataSetFromMatrix(countData = Dermo_counts_6h,
                                    colData = Dermo_coldata_6h,
                                    design = ~ FamCode + Treat)

Dermo_dds_36h <- DESeqDataSetFromMatrix(countData = Dermo_counts_36h,
                                       colData = Dermo_coldata_36h,
                                       design = ~Library_Prep_Date + FamCode + Treat)

Dermo_dds_7d <- DESeqDataSetFromMatrix(countData = Dermo_counts_7d,
                                       colData = Dermo_coldata_7d,
                                       design = ~Library_Prep_Date + FamCode + Treat)

# These designs says "test for the effects of family at each time point, controlling for the effect of Prep Date".
# The first factor gets treated as a co-variate and DESeq2 will try to eliminate the fold
# change that results from it. (see this for help: https://www.biostars.org/p/278684/)

# The dds objects below looks at the effect of time on the reference level of Family.
Dermo_dds_LB_time <- DESeqDataSetFromMatrix(countData = Dermo_counts_LB_all,
                                           colData = Dermo_coldata_LB_all,
                                           design = ~Library_Prep_Date + Timepoint + Treat)

Dermo_dds_DA_time <- DESeqDataSetFromMatrix(countData = Dermo_counts_DA_all,
                                            colData = Dermo_coldata_DA_all,
                                            design = ~Library_Prep_Date + Timepoint + Treat)

Dermo_dds_LB_middle_late <- DESeqDataSetFromMatrix(countData = Dermo_counts_LB_middle_late,
                                                   colData = Dermo_coldata_LB_middle_late,
                                                   design = ~Library_Prep_Date + Timepoint + Treat)
  
Dermo_dds_DA_middle_late <- DESeqDataSetFromMatrix(countData = Dermo_counts_DA_middle_late,
                                                   colData = Dermo_coldata_DA_middle_late,
                                                   design = ~Library_Prep_Date + Timepoint + Treat)

# About looking at contrasts: "Alternatively you can group the strain (with its different levels) and time ((with its different levels)
# into one factor, lets call it ALL. and by using contrast () you can look for the difference in log2 fold change between any combination of levels."

##### Prefiltering the data #####

# Data prefiltering helps decrease the size of the data set and get rid of
# rows with no data or very minimal data (<10). Apply a minimal filtering here as more stringent filtering will be applied later

Dermo_dds_family <- Dermo_dds_family[ rowSums(counts(Dermo_dds_family)) > 10, ]
Dermo_dds_timecourse_family <- Dermo_dds_timecourse_family[ rowSums(counts(Dermo_dds_timecourse_family)) > 10, ]
Dermo_dds_6h <- Dermo_dds_6h[ rowSums(counts(Dermo_dds_6h)) > 10, ]
Dermo_dds_36h <- Dermo_dds_36h[ rowSums(counts(Dermo_dds_36h)) > 10, ] 
Dermo_dds_7d <- Dermo_dds_7d[ rowSums(counts(Dermo_dds_7d)) > 10, ]
Dermo_dds_DA_time <- Dermo_dds_DA_time[ rowSums(counts(Dermo_dds_DA_time)) > 10, ]
Dermo_dds_LB_time <- Dermo_dds_LB_time[ rowSums(counts(Dermo_dds_LB_time)) > 10, ]
Dermo_dds_LB_middle_late <- Dermo_dds_LB_middle_late[rowSums(counts(Dermo_dds_LB_middle_late)) > 10,]
Dermo_dds_DA_middle_late <- Dermo_dds_DA_middle_late[rowSums(counts(Dermo_dds_DA_middle_late)) >10,]

#### VST Transformation of count data for visualization #### 

# VST Transforming the count data for visualization
  # this is only going to be performed for the full data set and not each of the individual data sets
# rlog transformation is recommended for small datasets (n<30) while VST is recommended for larger data sets
nrow(Dermo_coldata)
# there are 53 samples, VST is recommended

Dermo_dds_family_vsd <- vst(Dermo_dds_family, blind=FALSE)
Dermo_dds_6h_vsd <- vst(Dermo_dds_6h, blind=FALSE )
Dermo_dds_36h_vsd <- vst(Dermo_dds_36h, blind=FALSE)
Dermo_dds_7d_vsd <- vst(Dermo_dds_7d, blind=FALSE)
Dermo_dds_DA_time_vsd <- vst(Dermo_dds_DA_time, blind=FALSE)
Dermo_dds_LB_time_vsd <- vst(Dermo_dds_LB_time, blind=FALSE)
Dermo_dds_LB_middle_late_vsd <- vst(Dermo_dds_LB_middle_late, blind=FALSE)
Dermo_dds_DA_middle_late_vsd <- vst(Dermo_dds_DA_middle_late , blind=FALSE)


# Calculate sample distances to assess overall sample similarity
  # use function dist to calculate the Euclidean distance between samples. 
  # To ensure we have a roughly equal contribution from all genes, we use it on the VST data. 
  # We need to transpose the matrix of values using t, because the dist function expects the different samples to be rows of its argument, and different dimensions (here, genes) to be columns.

Dermo_dds_family_vsd_dist <- dist(t(assay(Dermo_dds_family_vsd)))
Dermo_dds_family_vsd_dist

#####  plot sample distances using a heatmap with pheatmap.#####
  # manually provide sampledists to the clustering_distance argument
  # make the column names of the matrix the family and the individual
Dermo_dds_family_vsd_dist_matrix <- as.matrix(Dermo_dds_family_vsd_dist)
rownames(Dermo_dds_family_vsd_dist_matrix) <- paste(Dermo_dds_family_vsd$FamCode, Dermo_dds_family_vsd$Ind, sep="-")
#colnames(Dermo_dds_family_vsd_dist_matrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(Dermo_dds_family_vsd_dist_matrix,
         clustering_distance_rows = Dermo_dds_family_vsd_dist,
         clustering_distance_cols = Dermo_dds_family_vsd_dist,
         col=colors)

# For the most part, the LB and the DA cluster together

##### PCA plot visualization of individuals in the family #####
plotPCA(Dermo_dds_family_vsd, intgroup=c("FamCode","Timepoint"))
  #the DA family samples are closely clustered, while the LB family is less close together,
  # but still falling out in the same location 

##### Differential Expression Analysis ####

## RUN DESEQ PIPELINE WITH DESEQ COMMAND

# run pipeline with single command because the formula has already been specified
# Steps: estimation of size factors (controlling for differences in the sequencing depth of the samples), 
# the estimation of dispersion values for each gene, 
# and fitting a generalized linear model.

Dermo_dds_family_deseq <- DESeq(Dermo_dds_family) 
Dermo_dds_6h_deseq <- DESeq(Dermo_dds_6h)
Dermo_dds_36h_deseq <- DESeq(Dermo_dds_36h)
Dermo_dds_7d_deseq <- DESeq(Dermo_dds_7d)
Dermo_dds_DA_time_deseq <- DESeq(Dermo_dds_DA_time)
Dermo_dds_LB_time_deseq <- DESeq(Dermo_dds_LB_time)
Dermo_dds_LB_middle_late_deseq <- DESeq(Dermo_dds_LB_middle_late)
Dermo_dds_DA_middle_late_deseq <- DESeq(Dermo_dds_DA_middle_late)

## Check the resultsNames object of each to look at the available coef

resultsNames(Dermo_dds_family_deseq) # "Intercept" "Library_Prep_Date_17_Dec_vs_15_Dec" "FamCode_LB_vs_DA"
resultsNames(Dermo_dds_6h_deseq) # "Intercept" "FamCode_LB_vs_DA"
resultsNames(Dermo_dds_36h_deseq) # "Intercept"  "Library_Prep_Date_17_Dec_vs_15_Dec" "FamCode_LB_vs_DA"
resultsNames(Dermo_dds_7d_deseq) # "Intercept" "Library_Prep_Date_17_Dec_vs_15_Dec" "FamCode_LB_vs_DA"  
resultsNames(Dermo_dds_DA_time_deseq) # Intercept"  "Library_Prep_Date_17_Dec_vs_15_Dec" "Timepoint_36h_vs_6h"  "Timepoint_7d_vs_6h"
resultsNames(Dermo_dds_LB_time_deseq) # "Intercept" "Library_Prep_Date_17_Dec_vs_15_Dec" "Timepoint_36h_vs_6h"  "Timepoint_7d_vs_6h"
resultsNames(Dermo_dds_LB_middle_late_deseq) # "Intercept"  "Library_Prep_Date_17_Dec_vs_15_Dec" "Timepoint_7d_vs_36h" 
resultsNames(Dermo_dds_DA_middle_late_deseq) # "Intercept" "Library_Prep_Date_17_Dec_vs_15_Dec" "Timepoint_7d_vs_36h"

## BUILD THE RESULTS OBJECT

  # Examining the results object, change alpha to p <0.05, looking at object metadata
  # use mcols to look at metadata for each table

# Full Family DESeq
Dermo_dds_family_res <- results(Dermo_dds_family_deseq, alpha=0.05, name="FamCode_LB_vs_DA")
Dermo_dds_family_res # comparison is log2 fold change (MLE): FamCode LB vs DA

# 6h both family DESeq
Dermo_dds_6h_res <- results(Dermo_dds_6h_deseq, alpha=0.05, name="FamCode_LB_vs_DA")
mcols(Dermo_dds_6h_res, use.names = TRUE) # mcols pulls up the dataframe object metadat about the meaning of each column

# 36h both family DESeq
Dermo_dds_36h_res <- results(Dermo_dds_36h_deseq, alpha=0.05, name="FamCode_LB_vs_DA")
mcols(Dermo_dds_36h_res, use.names = TRUE)

# 7d both family DESeq
Dermo_dds_7d_res <- results(Dermo_dds_7d_deseq, alpha=0.05, name="FamCode_LB_vs_DA")
mcols(Dermo_dds_7d_res, use.names = TRUE)

# All times DA family
resultsNames(Dermo_dds_DA_time_deseq)
Dermo_dds_DA_time_res_early_middle <- results(Dermo_dds_DA_time_deseq, name="Timepoint_36h_vs_6h", alpha=0.05)
mcols(Dermo_dds_DA_time_res_early_middle)
Dermo_dds_DA_time_res_early_late <- results(Dermo_dds_DA_time_deseq, name= "Timepoint_7d_vs_6h", alpha=0.05)
mcols(Dermo_dds_DA_time_res_early_late, use.names = TRUE)

# All times LB family
Dermo_dds_LB_time_res_early_middle <- results(Dermo_dds_LB_time_deseq, name="Timepoint_36h_vs_6h", alpha=0.05)
Dermo_dds_LB_time_res_early_late <- results(Dermo_dds_LB_time_deseq, name= "Timepoint_7d_vs_6h", alpha=0.05)
mcols(Dermo_dds_LB_time_res_early_middle, use.names = TRUE)
mcols(Dermo_dds_LB_time_res_early_late, use.names = TRUE)

# Middle late DA LB
Dermo_dds_DA_time_res_middle_late <- results(Dermo_dds_DA_middle_late_deseq, name="Timepoint_7d_vs_36h", alpha=0.05)
Dermo_dds_LB_time_res_middle_late <- results(Dermo_dds_LB_middle_late_deseq, name="Timepoint_7d_vs_36h", alpha=0.05)
mcols(Dermo_dds_DA_time_res_middle_late, use.names=TRUE)
mcols(Dermo_dds_LB_time_res_middle_late, use.names = TRUE)

#### Perform LFC Shrinkage with ApeGLM ####

## NOTES 
# Before plotting we need to apply an apeglm LFC shrinkage, set the coef as the specific comparison in the ResultsNames function of the deseq object
# Issue that the specific comparisons I set in my results formulas are not available in the ResultsNames coef list. Need to either add these, or go with the 36h vs. 6hr comparison

# NOTES from Michael love on Lfcshrinkage (https://support.bioconductor.org/p/77461/): 
# https://support.bioconductor.org/p/110307/ # very helpful distinction between lfcestimate and lfc shrinkage
# The difference between results() and lfcShrink() is that the former does not provide fold change shrinkage. 
# The latter function calls results() internally to create the p-value and adjusted p-value columns, 
# which provide inference on the maximum likelihood LFC. The shrunken fold changes are useful for ranking genes by 
# effect size and for visualization.
# The shrinkage is generally useful, which is why it is enabled by default. Full methods are described in the DESeq2 paper (see DESeq2 citation),
# but in short, it looks at the largest fold changes that are not due to low counts and uses these to inform a prior distribution. 
# So the large fold changes from genes with lots of statistical information are not shrunk, while the imprecise fold changes are shrunk. 
# This allows you to compare all estimated LFC across experiments, for example, which is not really feasible without the use of a prior.
# THE lfcshrinkage is not Affecting the p values at all, but its just shrinking the log2 fold change and calculating a new standard error for it 
# https://support.bioconductor.org/p/95695/

#notes on setting up coefficients for apeglm, https://support.bioconductor.org/p/115435/ , https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#extended-section-on-shrinkage-estimators

# Although apeglm cannot be used with contrast, we note that many designs can be easily rearranged such that what was a contrast becomes its own coefficient.
# In this case, the dispersion does not have to be estimated again, as the designs are equivalent, up to the meaning of the coefficients. 
# Instead, one need only run nbinomWaldTest to re-estimate MLE coefficients – these are necessary for apeglm – and then run lfcShrink specifying 
# the coefficient of interest in resultsNames(dds)
# The user would for example, either change the levels of dds$condition or replace the design using design(dds)<-, then run nbinomWaldTest followed by lfcShrink

# For each LFCshrink I can pass to it my res object for each so that I can keep my alpha setting at 0.05. Doing this procedure will 
# keep the p-values and padj from the results() call, and simply update the LFCs so they are posterior estimates.

## DECISION: GO WITH AVAILABLE CONTRASTS FROM THE RESULTSNAMES OBJECT FOR EACH, AND USE RES OBJECT TO KEEP ALPHA ADJUSTMENT

Dermo_dds_family_res_LFC <- lfcShrink(Dermo_dds_family_deseq, coef="FamCode_LB_vs_DA", type="apeglm", res=Dermo_dds_family_res)

Dermo_dds_6h_res_LFC <- lfcShrink(Dermo_dds_6h_deseq, coef="FamCode_LB_vs_DA", type="apeglm", res=Dermo_dds_6h_res)
Dermo_dds_36h_res_LFC <- lfcShrink(Dermo_dds_36h_deseq, coef="FamCode_LB_vs_DA", type="apeglm", res=Dermo_dds_36h_res)
Dermo_dds_7d_res_LFC <- lfcShrink(Dermo_dds_7d_deseq, coef="FamCode_LB_vs_DA", type="apeglm", res=Dermo_dds_7d_res)

Dermo_dds_DA_time_res_early_middle_LFC <- lfcShrink(Dermo_dds_DA_time_deseq, coef="Timepoint_36h_vs_6h", type="apeglm", res=Dermo_dds_DA_time_res_early_middle)
Dermo_dds_DA_time_res_early_late_LFC <- lfcShrink(Dermo_dds_DA_time_deseq, coef="Timepoint_7d_vs_6h", type="apeglm", res=Dermo_dds_DA_time_res_early_late)

Dermo_dds_LB_time_res_early_middle_LFC <- lfcShrink(Dermo_dds_LB_time_deseq,coef="Timepoint_36h_vs_6h", type="apeglm", res=Dermo_dds_LB_time_res_early_middle)
Dermo_dds_LB_time_res_early_late_LFC <- lfcShrink(Dermo_dds_LB_time_deseq, coef="Timepoint_7d_vs_6h", type="apeglm", res=Dermo_dds_LB_time_res_early_late)

Dermo_dds_DA_time_res_middle_late_LFC <- lfcShrink(Dermo_dds_DA_middle_late_deseq, coef="Timepoint_7d_vs_36h", type="apeglm", res=Dermo_dds_DA_time_res_middle_late)
Dermo_dds_LB_time_res_middle_late_LFC <- lfcShrink(Dermo_dds_LB_middle_late_deseq, coef="Timepoint_7d_vs_36h", type="apeglm", res=Dermo_dds_LB_time_res_middle_late)

## REVIEW RESULTS OBJECT SUMMARY ####
# SHOWS NUMBER OF SIGNIFICANT GENES

summary(Dermo_dds_timecourse_family_res) # NO LFC performed for this one yet as it may not end up being used
  # out of 53327 with nonzero total read count
  # adjusted p-value < 0.05
  # LFC > 0 (up)       : 18, 0.034%
  # LFC < 0 (down)     : 15, 0.028%
  # outliers [1]       : 491, 0.92%
  # low counts [2]     : 20306, 38%
  # (mean count < 8)

summary(Dermo_dds_family_res_LFC)
  #out of 53188 with nonzero total read count
  #adjusted p-value < 0.05
  #LFC > 0 (up)       : 5024, 9.4%
  #LFC < 0 (down)     : 5386, 10%
  #outliers [1]       : 0, 0%
  #low counts [2]     : 191, 0.36%
  #(mean count < 0)

summary(Dermo_dds_6h_res_LFC) # few genes differ at 6hr 
  #out of 46991 with nonzero total read count
  #adjusted p-value < 0.05
  #LFC > 0 (up)       : 654, 1.4%
  #LFC < 0 (down)     : 752, 1.6%
  #outliers [1]       : 2111, 4.5%
  #low counts [2]     : 4273, 9.1%
  #(mean count < 2)
  
summary(Dermo_dds_36h_res_LFC)
  #out of 49408 with nonzero total read count
  #adjusted p-value < 0.05
  #LFC > 0 (up)       : 2991, 6.1%
  #LFC < 0 (down)     : 3111, 6.3%
  #outliers [1]       : 1224, 2.5%
  #low counts [2]     : 2678, 5.4%
  #(mean count < 1)
  
summary(Dermo_dds_7d_res_LFC)
  #out of 49681 with nonzero total read count
  #adjusted p-value < 0.05
  #LFC > 0 (up)       : 3508, 7.1%
  #LFC < 0 (down)     : 3323, 6.7%
  #outliers [1]       : 957, 1.9%
  #low counts [2]     : 1808, 3.6%
  #(mean count < 1)
  
summary(Dermo_dds_DA_time_res_early_middle_LFC)
#out of 49827 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 2998, 6%
#LFC < 0 (down)     : 4182, 8.4%
#outliers [1]       : 1157, 2.3%
#low counts [2]     : 7107, 14%
#(mean count < 2)
  
summary(Dermo_dds_DA_time_res_early_late_LFC)
  #out of 49870 with nonzero total read count
  #adjusted p-value < 0.1
  #LFC > 0 (up)       : 5280, 11%
  #LFC < 0 (down)     : 4903, 9.8%
  #outliers [1]       : 868, 1.7%
  #low counts [2]     : 7257, 15%
  #(mean count < 2)

summary(Dermo_dds_DA_time_res_middle_late_LFC) # Comparison of middle to late has very few significant genes!! Suggesting that 36hr
# and 7 days are not that interesting compared to one another...this is really cool..looking at gene expression at 36hr may be as good as looking at it at 7days!
    #out of 48162 with nonzero total read count
    #adjusted p-value < 0.05
    #LFC > 0 (up)       : 354, 0.74%
    #LFC < 0 (down)     : 248, 0.51%
    #outliers [1]       : 994, 2.1%
    #low counts [2]     : 15891, 33%
    #(mean count < 9)

summary(Dermo_dds_LB_time_res_early_middle_LFC)
  #out of 49870 with nonzero total read count
  #adjusted p-value < 0.05
  #LFC > 0 (up)       : 4356, 8.7%
  #LFC < 0 (down)     : 4298, 8.6%
  #outliers [1]       : 868, 1.7%
  #low counts [2]     : 6353, 13%
  #(mean count < 2)
  
summary(Dermo_dds_LB_time_res_early_late_LFC)
  #out of 49870 with nonzero total read count
  #adjusted p-value < 0.05
  #LFC > 0 (up)       : 4137, 8.3%
  #LFC < 0 (down)     : 4103, 8.2%
  #outliers [1]       : 868, 1.7%
  #low counts [2]     : 9081, 18%
  #(mean count < 3)
 
summary(Dermo_dds_LB_time_res_middle_late_LFC)  # Comparison of middle to late has very few significant genes!! Suggesting that 36hr
  #out of 48660 with nonzero total read count
  #adjusted p-value < 0.05
  #LFC > 0 (up)       : 17, 0.035%
  #LFC < 0 (down)     : 0, 0%
  #outliers [1]       : 1114, 2.3%
  #low counts [2]     : 8714, 18%
  #(mean count < 3)


#### Exploratory Plotting of Results ####

## MA Plotting

plotMA(Dermo_dds_LB_time_res_early_late_LFC , ylim = c(-5, 5))
plotMA(Dermo_dds_LB_time_res_early_middle_LFC , ylim = c(-5, 5))
plotMA(Dermo_dds_DA_time_res_early_late_LFC , ylim = c(-5, 5))
plotMA(Dermo_dds_DA_time_res_early_middle_LFC, ylim = c(-5, 5))
plotMA(Dermo_dds_DA_time_res_middle_late_LFC, ylim = c(-5, 5))
plotMA(Dermo_dds_LB_time_res_middle_late_LFC, ylim = c(-5, 5))
plotMA(Dermo_dds_7d_res_LFC, ylim = c(-5, 5))
plotMA(Dermo_dds_36h_res_LFC, ylim = c(-5, 5))
plotMA(Dermo_dds_family_res_LFC, ylim = c(-5, 5)) # the plot looks strange!
plotMA(Dermo_dds_6h_res_LFC, ylim = c(-5, 5)) # this plot looks normal

## Histogram of P values ##
  # exclude genes with very small counts to avoid spikes and plot using the LFCshrinkage

hist(Dermo_dds_LB_time_res_early_late_LFC$padj[Dermo_dds_LB_time_res_early_late_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
hist(Dermo_dds_LB_time_res_early_middle_LFC$padj[Dermo_dds_LB_time_res_early_middle_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
hist(Dermo_dds_DA_time_res_early_late_LFC$padj[Dermo_dds_DA_time_res_early_late_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
hist(Dermo_dds_DA_time_res_early_middle_LFC$padj[Dermo_dds_DA_time_res_early_middle_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
hist(Dermo_dds_DA_time_res_middle_late_LFC$padj[Dermo_dds_DA_time_res_middle_late_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
hist(Dermo_dds_LB_time_res_middle_late_LFC$padj[Dermo_dds_LB_time_res_middle_late_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
hist(Dermo_dds_7d_res_LFC$padj[Dermo_dds_7d_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
hist(Dermo_dds_36h_res_LFC$padj[Dermo_dds_36h_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
hist(Dermo_dds_family_res_LFC$padj[Dermo_dds_family_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
hist(Dermo_dds_6h_res_LFC$padj[Dermo_dds_6h_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white") # few are significantly different at 6h

#### Subsetting Significant Genes by padj < 0.05 and L2FC of greater than 1.0 or less than -1.0 #####

# again, only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
Dermo_dds_family_res_LFC_sig <- subset(Dermo_dds_family_res_LFC, padj < 0.05)
Dermo_dds_family_res_LFC_sig$transcript_id <- row.names(Dermo_dds_family_res_LFC_sig)
Dermo_dds_family_res_LFC_sig <- as.data.frame(Dermo_dds_family_res_LFC_sig)
Dermo_dds_family_res_LFC_sig <- Dermo_dds_family_res_LFC_sig %>% filter(abs(log2FoldChange) >= 1.0)
nrow(Dermo_dds_family_res_LFC_sig) #3636
nrow(Dermo_dds_family_res_LFC) # 53379

Dermo_dds_6h_res_LFC_sig <- subset(Dermo_dds_6h_res_LFC, padj < 0.05)
Dermo_dds_6h_res_LFC_sig$transcript_id <- row.names(Dermo_dds_6h_res_LFC_sig)
Dermo_dds_6h_res_LFC_sig <- as.data.frame(Dermo_dds_6h_res_LFC_sig)
Dermo_dds_6h_res_LFC_sig <- Dermo_dds_6h_res_LFC_sig %>% filter(abs(log2FoldChange) >=1.0)
nrow(Dermo_dds_6h_res_LFC) # 46991
nrow(Dermo_dds_6h_res_LFC_sig) # 954

Dermo_dds_36h_res_LFC_sig <- subset(Dermo_dds_36h_res_LFC,padj < 0.05)
Dermo_dds_36h_res_LFC_sig$transcript_id <- row.names(Dermo_dds_36h_res_LFC_sig)
Dermo_dds_36h_res_LFC_sig <- as.data.frame(Dermo_dds_36h_res_LFC_sig)
Dermo_dds_36h_res_LFC_sig <- Dermo_dds_36h_res_LFC_sig %>% filter(abs(log2FoldChange) >=1.0)
nrow(Dermo_dds_36h_res_LFC) # 49408
nrow(Dermo_dds_36h_res_LFC_sig) # 4258

Dermo_dds_7d_res_LFC_sig <- subset(Dermo_dds_7d_res_LFC,padj < 0.05)
Dermo_dds_7d_res_LFC_sig$transcript_id <- row.names(Dermo_dds_7d_res_LFC_sig)
Dermo_dds_7d_res_LFC_sig <- as.data.frame(Dermo_dds_7d_res_LFC_sig)
Dermo_dds_7d_res_LFC_sig<- Dermo_dds_7d_res_LFC_sig  %>% filter(abs(log2FoldChange) >=1.0)
nrow(Dermo_dds_7d_res_LFC) # 49681
nrow(Dermo_dds_7d_res_LFC_sig) # 4299

Dermo_dds_LB_time_res_early_late_LFC_sig   <- subset(Dermo_dds_LB_time_res_early_late_LFC, padj < 0.05)
Dermo_dds_LB_time_res_early_late_LFC_sig$transcript_id <- row.names(Dermo_dds_LB_time_res_early_late_LFC_sig)
Dermo_dds_LB_time_res_early_late_LFC_sig <- as.data.frame(Dermo_dds_LB_time_res_early_late_LFC_sig)
Dermo_dds_LB_time_res_early_late_LFC_sig <- Dermo_dds_LB_time_res_early_late_LFC_sig %>% filter(abs(log2FoldChange) >=1.0)
nrow(Dermo_dds_LB_time_res_early_late_LFC) # 49870
nrow(Dermo_dds_LB_time_res_early_late_LFC_sig) # 5513
 
Dermo_dds_LB_time_res_early_middle_LFC_sig <- subset(Dermo_dds_LB_time_res_early_middle_LFC, padj < 0.05)
Dermo_dds_LB_time_res_early_middle_LFC_sig$transcript_id <- row.names(Dermo_dds_LB_time_res_early_middle_LFC_sig)
Dermo_dds_LB_time_res_early_middle_LFC_sig <- as.data.frame(Dermo_dds_LB_time_res_early_middle_LFC_sig)
Dermo_dds_LB_time_res_early_middle_LFC_sig <- Dermo_dds_LB_time_res_early_middle_LFC_sig %>% filter(abs(log2FoldChange) >=1.0)
nrow(Dermo_dds_LB_time_res_early_middle_LFC) # 49870
nrow(Dermo_dds_LB_time_res_early_middle_LFC_sig) # 5976

Dermo_dds_LB_time_res_middle_late_LFC_sig <- subset(Dermo_dds_LB_time_res_middle_late_LFC, padj < 0.05)
Dermo_dds_LB_time_res_middle_late_LFC_sig$transcript_id <- row.names(Dermo_dds_LB_time_res_middle_late_LFC_sig)
Dermo_dds_LB_time_res_middle_late_LFC_sig <- as.data.frame(Dermo_dds_LB_time_res_middle_late_LFC_sig)
Dermo_dds_LB_time_res_middle_late_LFC_sig <- Dermo_dds_LB_time_res_middle_late_LFC_sig %>% filter(abs(log2FoldChange) >=1.0)
nrow(Dermo_dds_LB_time_res_middle_late_LFC) # 48660
nrow(Dermo_dds_LB_time_res_middle_late_LFC_sig) # 7

Dermo_dds_DA_time_res_early_late_LFC_sig <- subset(Dermo_dds_DA_time_res_early_late_LFC, padj < 0.05)
Dermo_dds_DA_time_res_early_late_LFC_sig$transcript_id <- row.names(Dermo_dds_DA_time_res_early_late_LFC_sig)
Dermo_dds_DA_time_res_early_late_LFC_sig <- as.data.frame(Dermo_dds_DA_time_res_early_late_LFC_sig)
Dermo_dds_DA_time_res_early_late_LFC_sig <- Dermo_dds_DA_time_res_early_late_LFC_sig  %>% filter(abs(log2FoldChange) >=1.0)
nrow(Dermo_dds_DA_time_res_early_late_LFC) # 49827
nrow(Dermo_dds_DA_time_res_early_late_LFC_sig) # 4828

Dermo_dds_DA_time_res_early_middle_LFC_sig <- subset(Dermo_dds_DA_time_res_early_middle_LFC, padj < 0.05)
Dermo_dds_DA_time_res_early_middle_LFC_sig$transcript_id <- row.names(Dermo_dds_DA_time_res_early_middle_LFC_sig)
Dermo_dds_DA_time_res_early_middle_LFC_sig <- as.data.frame(Dermo_dds_DA_time_res_early_middle_LFC_sig)
Dermo_dds_DA_time_res_early_middle_LFC_sig <- Dermo_dds_DA_time_res_early_middle_LFC_sig %>% filter(abs(log2FoldChange) >=1.0)
nrow(Dermo_dds_DA_time_res_early_middle_LFC) # 49827
nrow(Dermo_dds_DA_time_res_early_middle_LFC_sig) # 4378
 
Dermo_dds_DA_time_res_middle_late_LFC_sig <- subset(Dermo_dds_DA_time_res_middle_late_LFC, padj < 0.05)
Dermo_dds_DA_time_res_middle_late_LFC_sig$transcript_id <- row.names(Dermo_dds_DA_time_res_middle_late_LFC_sig)
Dermo_dds_DA_time_res_middle_late_LFC_sig <- as.data.frame(Dermo_dds_DA_time_res_middle_late_LFC_sig)
Dermo_dds_DA_time_res_middle_late_LFC_sig <- Dermo_dds_DA_time_res_middle_late_LFC_sig %>% filter(abs(log2FoldChange) >= 1.0)
nrow(Dermo_dds_DA_time_res_middle_late_LFC) # 48162
nrow(Dermo_dds_DA_time_res_middle_late_LFC_sig) # 168

#### Which transcripts are different between timepoint comparisons? #####

## Comparison of transcripts within families
# Transcripts in DA in the 6 vs 7d comparison not in the 6 vs 36r comparison
DA_early_late_not_early_middle <- Dermo_dds_DA_time_res_early_late_LFC_sig$transcript_id[!(Dermo_dds_DA_time_res_early_late_LFC_sig$transcript_id %in% Dermo_dds_DA_time_res_early_middle_LFC_sig$transcript_id)] 
length(DA_early_late_not_early_middle) #1561 , DA_early_late has 4828

# Transcripts in DA 6 vs 36hr not in 6 vs 7d
DA_early_middle_not_early_late <- Dermo_dds_DA_time_res_early_middle_LFC_sig$transcript_id[!(Dermo_dds_DA_time_res_early_middle_LFC_sig$transcript_id %in% Dermo_dds_DA_time_res_early_late_LFC_sig$transcript_id)] 
length(DA_early_middle_not_early_late) #1111, DA_early_middle has 4375

# Transcripts in LB in the 6 vs 7d comparison not in the 6 vs 36r comparison
LB_early_late_not_early_middle <- Dermo_dds_LB_time_res_early_late_LFC_sig$transcript_id[!(Dermo_dds_LB_time_res_early_late_LFC_sig$transcript_id %in% Dermo_dds_LB_time_res_early_middle_LFC_sig$transcript_id)] 
length(LB_early_late_not_early_middle) #1139 , LB_early_late has 5305

# Transcripts in LB 6 vs 36hr not in 6 vs 7d
LB_early_middle_not_early_late <- Dermo_dds_LB_time_res_early_middle_LFC_sig$transcript_id[!(Dermo_dds_LB_time_res_early_middle_LFC_sig$transcript_id %in% Dermo_dds_LB_time_res_early_late_LFC_sig$transcript_id)] 
length(LB_early_middle_not_early_late) #1602, LB_early_middle has 5644

## Comparison of transcripts BETWEEN families at different timepoints
sixhr_vs_36hr <- Dermo_dds_6h_res_LFC_sig$transcript_id[!(Dermo_dds_6h_res_LFC_sig$transcript_id %in% Dermo_dds_36h_res_LFC_sig$transcript_id)] 
length(sixhr_vs_36hr) #445

sixhr_vs_7d <- Dermo_dds_6h_res_LFC_sig$transcript_id[!(Dermo_dds_6h_res_LFC_sig$transcript_id %in% Dermo_dds_7d_res_LFC_sig$transcript_id)] 
length(sixhr_vs_7d ) #378

thirtysix_vs_7d <- Dermo_dds_36h_res_LFC_sig$transcript_id[!(Dermo_dds_36h_res_LFC_sig$transcript_id %in% Dermo_dds_7d_res_LFC_sig$transcript_id)] 
length(thirtysix_vs_7d) # 2133


#### Graphing Significant Genes Between Timepoints ####

# Make a table with the numbers for graphing 

Coefficient <- c("FamCode_LB_vs_DA",
                  "FamCode_LB_vs_DA",
                  "FamCode_LB_vs_DA",
                  "FamCode_LB_vs_DA",
                  "Timepoint_36h_vs_6h",
                 "Timepoint_7d_vs_6h",
                 "Timepoint_36h_vs_6h",
                 "Timepoint_7d_vs_6h",
                 "Timepoint_7d_vs_36h",
                 "Timepoint_7d_vs_36h",
                 "FamCode_LB_vs_DA",
                 "FamCode_LB_vs_DA",
                 "FamCode_LB_vs_DA",
                 "FamCode_LB_vs_DA",
                 "Timepoint_36h_vs_6h",
                 "Timepoint_7d_vs_6h",
                 "Timepoint_36h_vs_6h",
                 "Timepoint_7d_vs_6h",
                 "Timepoint_7d_vs_36h",
                 "Timepoint_7d_vs_36h")
Variable_Tested <- c("FamCode",
                     "FamCode",
                     "FamCode",
                     "FamCode",
                     "Timepoint",                     "Timepoint",
                     "Timepoint",                     "Timepoint",
                     "Timepoint",                     "Timepoint",
                     "FamCode",
                     "FamCode",
                     "FamCode",
                     "FamCode",
                     "Timepoint",                     "Timepoint",
                     "Timepoint",                     "Timepoint",
                     "Timepoint",                     "Timepoint")
Families_Included <- c("All",
                       "All",
                       "All",
                       "All",
                       "DA",
                       "DA",
                       "LB",
                       "LB",
                       "DA",
                       "LB",
                       "All",
                       "All",
                       "All",
                       "All",
                       "DA",
                       "DA",
                       "LB",
                       "LB",
                       "DA",
                       "LB")
Timepoint_Compared <- c("All",
                        "6hr",
                        "36hr",
                        "7d",
                        "6hr_36hr",
                        "6hr_7d",
                        "6hr_36hr",
                        "6hr_7d",
                        "36hr_7d",
                        "36hr_7d",
                        "All",
                        "6hr",
                        "36hr",
                        "7d",
                        "6hr_36hr",
                        "6hr_7d",
                        "6hr_36hr",
                        "6hr_7d",
                        "36hr_7d",
                        "36hr_7d")
DataSet <- c("Dermo_dds_family_res_LFC",
             "Dermo_dds_6h_res_LFC",
             "Dermo_dds_36h_res_LFC",
             "Dermo_dds_7d_res_LFC",
             "Dermo_dds_DA_time_res_early_middle_LFC",
             "Dermo_dds_DA_time_res_early_late_LFC ",
             "Dermo_dds_LB_time_res_early_middle_LFC",
             "Dermo_dds_LB_time_res_early_late_LFC",
             "Dermo_dds_DA_time_res_middle_late_LFC",
             "Dermo_dds_LB_time_res_middle_late_LFC",
             "Dermo_dds_family_res_LFC",
             "Dermo_dds_6h_res_LFC",
             "Dermo_dds_36h_res_LFC",
             "Dermo_dds_7d_res_LFC",
             "Dermo_dds_DA_time_res_early_middle_LFC",
             "Dermo_dds_DA_time_res_early_late_LFC ",
             "Dermo_dds_LB_time_res_early_middle_LFC",
             "Dermo_dds_LB_time_res_early_late_LFC",
             "Dermo_dds_DA_time_res_middle_late_LFC",
             "Dermo_dds_LB_time_res_middle_late_LFC")

Value_type <- c("Sig_Genes",
           "Sig_Genes",
           "Sig_Genes",
           "Sig_Genes",
           "Sig_Genes",
           "Sig_Genes",
           "Sig_Genes",
           "Sig_Genes",
           "Sig_Genes",
           "Sig_Genes",
           "Total_Genes",
           "Total_Genes",
           "Total_Genes",
           "Total_Genes",
           "Total_Genes",
           "Total_Genes",
           "Total_Genes",
           "Total_Genes",
           "Total_Genes",
           "Total_Genes") 


Value <- c("3843",
           "906",
           "4261",
           "4419",
           "4375",
           "4828",
           "5644",
           "5305",
           "146",
           "6",
           "53379",
           "46991",
           "49408",
           "49618",
           "49827",
           "49827",
           "49870",
           "49870",
           "48162",
           "48660")


Comparison_Sig_Genes <- data.frame(Coefficient, Variable_Tested, Families_Included, DataSet, Value_type, Value, Timepoint_Compared)
Comparison_Sig_Genes$Value <- as.numeric(as.character(Comparison_Sig_Genes$Value))

ggplot(Comparison_Sig_Genes, aes(x=DataSet, y=Value, fill = Value_type)) + geom_col(position="dodge") +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + ylab("Number of Transcripts")

Comparison_Sig_Genes_FamCode <- Comparison_Sig_Genes %>% filter(Variable_Tested == "FamCode")
Comparison_DA <- Comparison_Sig_Genes %>% filter(Families_Included == "DA")
Comparison_LB <- Comparison_Sig_Genes %>% filter(Families_Included == "LB")
class(Comparison_Sig_Genes$Value)

# Set factor levels
levels(Comparison_Sig_Genes_FamCode$Timepoint_Compared)
Comparison_Sig_Genes_FamCode$Timepoint_Compared <- factor(Comparison_Sig_Genes_FamCode$Timepoint_Compared, levels=c("6hr","36hr","7d", 
                                                                                                                     "6hr_36hr",  "36hr_7d", "6hr_7d", "All"))
levels(Comparison_DA$Timepoint_Compared)
Comparison_DA$Timepoint_Compared <- droplevels(Comparison_DA$Timepoint_Compared)
Comparison_DA$Timepoint_Compared <- factor(Comparison_DA$Timepoint_Compared, levels=c("6hr_36hr","36hr_7d","6hr_7d"))
levels(Comparison_LB$Timepoint_Compared)
Comparison_LB$Timepoint_Compared <- droplevels(Comparison_LB$Timepoint_Compared)
Comparison_LB$Timepoint_Compared <- factor(Comparison_LB$Timepoint_Compared, levels=c("6hr_36hr","36hr_7d","6hr_7d"))

## Comparison of Gene lists ##
# Differentially Expressed Transcripts Between DA and LB
ggplot(Comparison_Sig_Genes_FamCode, aes(x=Timepoint_Compared, y=Value, fill = Value_type)) + geom_col(position="dodge") +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + ylab("Number of Transcripts") + 
  ggtitle("Differentially Expressed Transcripts Between DA and LB")  

# Differentially Expressed Transcripts Over Time in DA
ggplot(Comparison_DA, aes(x=Timepoint_Compared, y=Value, fill = Value_type)) + geom_col(position="dodge") +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + ylab("Number of Transcripts") +  xlab("Timepoint Compared") +
  ggtitle("Differentially Expressed Transcripts Over Time in DA")

# Differentially Expressed Transcripts Over Time in LB
ggplot(Comparison_LB, aes(x=Timepoint_Compared, y=Value, fill = Value_type)) + geom_col(position="dodge") +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + ylab("Number of Transcripts") + xlab("Timepoint Compared") +
  ggtitle("Differentially Expressed Transcripts Over Time in LB")

#### Import Genome Annotation and add GO terms ####

# Import gff file with rtracklayer
C_vir_rtracklayer <- import("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/C_VIRGINICA_PIPELINE/ref_C_virginica-3.0_top_level.gff3")
C_vir_rtracklayer <- as.data.frame(C_vir_rtracklayer)

# Isolate transcript lines in GFF annotation so that I can use Batch Entrez lookup for their parent proteins
C_vir_rtracklayer_transcripts <- filter(C_vir_rtracklayer, grepl("XM", transcript_id))
C_vir_rtracklayer_transcripts_unique <- subset(C_vir_rtracklayer_transcripts, !duplicated(transcript_id))
write.table(file="./Dermo_2015_Analysis/OUTPUT/C_vir_unique_transcripts.txt", C_vir_rtracklayer_transcripts_unique$transcript_id) 
  # remove first column and row with the rownames: cut -d' ' -f2 C_vir_unique_transcripts.txt | tail -n +2 > C_vir_unique_transcripts_cut.txt
  # remove quotes: sed 's/\"//g' C_vir_unique_transcripts_cut.txt > C_vir_unique_transcripts_cut_no_quotes.txt
  # Split in terminal to multiple files: split -l 10000 C_vir_unique_transcripts_cut_no_quotes.txt
  # mv output to C_vir_unique_transcripts*.txt
  # In batch entrez select the nucleotide format ,and then don't highlight anything and click on "send to" and select, "complete record" then "file" then "GFF3"
  # extract XM and XP information
  # Combined all files using 'for i in x*_batch.txt ; do cat $i >> combined_batch_lookup.txt; done'
  # Extracting Gnomon CDS entry for each transcript by grepping for "cds-" which comes before the XP name for every line: grep "cds-" combined_batch_lookup.txt > combined_batch_lookup_cds.txt
  # changing file suffix to gff3 so it can be imported 

#upload Batch entrez NCBI gff2 format 
C_vir_XM_with_XP <- import("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/C_VIRGINICA_PIPELINE/Dermo_2015_Analysis/OUTPUT/combined_batch_lookup_cds.gff")
C_vir_XM_with_XP <- as.data.frame(C_vir_XM_with_XP)
C_vir_XM_with_XP_index <- C_vir_XM_with_XP[,c("seqnames","protein_id")]
colnames(C_vir_XM_with_XP_index)[1] <- "transcript_id"

# Load in Interproscan GO annotation from LSU Kevin using Rtracklayer
# cat all the edited header removed files from Interproscan Kevin LSU files Import gff file with rtracklayer
# for i in edited*.gff3; do cat $i >> combined_CV_prot_id_interproscan.gff3 ; done

# In terminal, combine all the gff files and remove the comment lines from the script 
Cvir_Interproscan1 <- import("/Users/erinroberts/Documents/PhD_Research/GOMEZCHIARRI2_FILES_DAILYNOTES/GENOME_TRANSCRIPTOME_SEQUENCES/IP5_out/protein_id/combined_CV_prot_id_interproscan.gff3")
class(Cvir_Interproscan1)
# convert to dataframe object using annoGR2DF (repitools)
Cvir_Interproscan_DF <- annoGR2DF(Cvir_Interproscan1)
class(Cvir_Interproscan_DF) #data.frame

# Interproscan file has multiple lines per entry, need to subset for the lines that do have a GO entry because not all of them do, then subset for unique protein lines
class(Cvir_Interproscan_DF$Ontology_term)
Cvir_Interproscan_DF$Ontology_term <- as.character(Cvir_Interproscan_DF$Ontology_term)
C_vir_Interproscan_GO <- Cvir_Interproscan_DF %>% filter(Ontology_term !="character(0)")
head(C_vir_Interproscan_GO)

# keep lines with unique protein and GO terms 
C_vir_Interproscan_GO_unique <- C_vir_Interproscan_GO[!duplicated(C_vir_Interproscan_GO[,c("chr","Ontology_term")]),]

# Format GO term column correctly 
C_vir_Interproscan_GO_unique$Ontology_term <- gsub("[^[:alnum:][:blank:]+:/\\-]", "", C_vir_Interproscan_GO_unique$Ontology_term )
C_vir_Interproscan_GO_unique$Ontology_term <- gsub(C_vir_Interproscan_GO_unique$Ontology_term, pattern="c\\", replacement="", fixed=TRUE)
C_vir_Interproscan_GO_unique$Ontology_term <- gsub(C_vir_Interproscan_GO_unique$Ontology_term, pattern="\\", replacement="", fixed=TRUE)
C_vir_Interproscan_GO_unique$Ontology_term <- gsub(C_vir_Interproscan_GO_unique$Ontology_term, pattern=" ",replacement =",", fixed=TRUE)
class(C_vir_Interproscan_GO_unique)

# merge GO terms with the same protein name so there aren't multiple lines for a transcript
C_vir_Interproscan_GO_unique <- ddply(C_vir_Interproscan_GO_unique, "chr", summarize, Combined_ontology = toString(Ontology_term))

# Remove duplicate GO strings in the same column and create new column
C_vir_Interproscan_GO_unique$Combined_ontology <- gsub(" ", "", C_vir_Interproscan_GO_unique$Combined_ontology)
C_vir_Interproscan_GO_unique <- mutate(C_vir_Interproscan_GO_unique, unique_go = map(str_split(Combined_ontology, ","), unique))

# make new column that removes the "_ORF" from the end of the seqnames columns
C_vir_Interproscan_GO_unique$protein_id <- str_remove(C_vir_Interproscan_GO_unique$chr, "_ORF")

# merge XM and XP list with Interproscan list
C_vir_Interproscan_GO_unique_XM_merged <- C_vir_Interproscan_GO_unique %>% left_join(select(C_vir_XM_with_XP_index, "transcript_id","protein_id"), by = "protein_id")

# Merge Interproscan list onto transcript annotation using left_join of the TRANSCRIPT ID so that it joins correctly at the level of transcripts
# non unique file
C_vir_rtracklayer_GO <-  C_vir_rtracklayer %>% left_join(select(C_vir_Interproscan_GO_unique_XM_merged,"unique_go","transcript_id"), by = "transcript_id")

#unique file with one line per transcript
C_vir_rtracklayer_transcripts_GO <- C_vir_rtracklayer_transcripts_unique %>% left_join(select(C_vir_Interproscan_GO_unique_XM_merged,"unique_go","transcript_id"), by = "transcript_id")

# are they all NULL?
C_vir_rtracklayer_GO %>% filter(unique_go !="NULL") # nope they are not!!

# Export for use in cluster later
class(C_vir_rtracklayer_transcripts_GO)
# coerce data frame to be all character
# First coerce the data.frame to all-character
C_vir_rtracklayer_transcripts_GO_2 = data.frame(lapply(C_vir_rtracklayer_transcripts_GO, as.character), stringsAsFactors=FALSE)
head(C_vir_rtracklayer_transcripts_GO_2)
# now write to csv
write.csv(C_vir_rtracklayer_transcripts_GO_2, file="C_vir_rtracklayer_transcripts_GO.csv")

#### Annotate significant genes in each results data frame ####

Dermo_dds_family_res_LFC_sig_annot <- Dermo_dds_family_res_LFC_sig %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product", "gene")], by = "transcript_id")
head(arrange(Dermo_dds_family_res_LFC_sig_annot,-log2FoldChange),n=50) # arrange in descending order

Dermo_dds_6h_res_LFC_sig_annot <- Dermo_dds_6h_res_LFC_sig %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product", "gene")], by = "transcript_id")
head(arrange(Dermo_dds_6h_res_LFC_sig_annot,-log2FoldChange),n=50) #most differentially expressed compared between families at 6hr

Dermo_dds_36h_res_LFC_sig_annot <- Dermo_dds_36h_res_LFC_sig %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product", "gene")], by = "transcript_id")
head(arrange(Dermo_dds_36h_res_LFC_sig_annot ,-log2FoldChange),n=50) #most differentially expressed compared between families at 36hr

Dermo_dds_7d_res_LFC_sig_annot <- Dermo_dds_7d_res_LFC_sig %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product", "gene")], by = "transcript_id")
head(arrange(Dermo_dds_7d_res_LFC_sig_annot, -log2FoldChange),n=50) #most differentially expressed compared between families at 36hr

Dermo_dds_LB_time_res_early_late_LFC_sig_annot <- Dermo_dds_LB_time_res_early_late_LFC_sig %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product", "gene")], by = "transcript_id")
head(arrange(Dermo_dds_LB_time_res_early_late_LFC_sig_annot, -log2FoldChange),n=50)

Dermo_dds_LB_time_res_early_middle_LFC_sig_annot <- Dermo_dds_LB_time_res_early_middle_LFC_sig %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product", "gene")], by = "transcript_id")
head(arrange(Dermo_dds_LB_time_res_early_middle_LFC_sig_annot,-log2FoldChange),n=50)

Dermo_dds_LB_time_res_middle_late_LFC_sig_annot <- Dermo_dds_LB_time_res_middle_late_LFC_sig %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product", "gene")], by = "transcript_id")
head(arrange(Dermo_dds_LB_time_res_middle_late_LFC_sig_annot,-log2FoldChange),n=50)

Dermo_dds_DA_time_res_early_late_LFC_sig_annot <- Dermo_dds_DA_time_res_early_late_LFC_sig %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product", "gene")], by = "transcript_id")
head(arrange(Dermo_dds_DA_time_res_early_late_LFC_sig_annot, -log2FoldChange),n=50)

Dermo_dds_DA_time_res_early_middle_LFC_sig_annot <- Dermo_dds_DA_time_res_early_middle_LFC_sig %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product", "gene")], by = "transcript_id")
head(arrange(Dermo_dds_DA_time_res_early_middle_LFC_sig_annot, -log2FoldChange),n=50)

Dermo_dds_DA_time_res_middle_late_LFC_sig_annot <- Dermo_dds_DA_time_res_middle_late_LFC_sig %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product", "gene")], by = "transcript_id")
head(arrange(Dermo_dds_DA_time_res_middle_late_LFC_sig_annot,-log2FoldChange),n=50)

#### Upset plots of overall gene expression changes between samples ####
# helpful tutorial for doing this: http://genomespot.blogspot.com/2017/09/upset-plots-as-replacement-to-venn.html
# http://crazyhottommy.blogspot.com/2016/01/upset-plot-for-overlapping-chip-seq.html
# UpsetR vignette: https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
# first extract gene list for each set
Dermo_dds_family_res_LFC_sig_id <- Dermo_dds_family_res_LFC_sig$transcript_id
Dermo_dds_6h_res_LFC_sig_id  <-Dermo_dds_6h_res_LFC_sig$transcript_id
Dermo_dds_36h_res_LFC_sig_id <- Dermo_dds_36h_res_LFC_sig$transcript_id
Dermo_dds_7d_res_LFC_sig_id <- Dermo_dds_7d_res_LFC_sig$transcript_id
Dermo_dds_DA_time_res_early_middle_LFC_sig_id <- Dermo_dds_DA_time_res_early_middle_LFC_sig$transcript_id
Dermo_dds_DA_time_res_early_late_LFC_sig_id <- Dermo_dds_DA_time_res_early_late_LFC_sig$transcript_id 
Dermo_dds_LB_time_res_early_middle_LFC_sig_id <- Dermo_dds_LB_time_res_early_middle_LFC_sig$transcript_id
Dermo_dds_LB_time_res_early_late_LFC_sig_id <- Dermo_dds_LB_time_res_early_late_LFC_sig$transcript_id
Dermo_dds_DA_time_res_middle_late_LFC_sig_id <- Dermo_dds_DA_time_res_middle_late_LFC_sig$transcript_id
Dermo_dds_LB_time_res_middle_late_LFC_sig_id <- Dermo_dds_LB_time_res_middle_late_LFC_sig$transcript_id

# Make each into a dataframe with two columns, gene name and sample.
Dermo_dds_family_res_LFC_sig_id <- as.data.frame(Dermo_dds_family_res_LFC_sig_id)
Dermo_dds_family_res_LFC_sig_id$Comparison <- "Family_Comparison_all_timepoints"
colnames(Dermo_dds_family_res_LFC_sig_id)[1] <- "Transcript_id"

Dermo_dds_6h_res_LFC_sig_id <- as.data.frame(Dermo_dds_6h_res_LFC_sig_id)
Dermo_dds_6h_res_LFC_sig_id$Comparison <- "Dermo_6hr_LB_vs_DA"
colnames(Dermo_dds_6h_res_LFC_sig_id)[1] <- "Transcript_id"
  
Dermo_dds_36h_res_LFC_sig_id <- as.data.frame(Dermo_dds_36h_res_LFC_sig_id)
Dermo_dds_36h_res_LFC_sig_id$Comparison <- "Dermo_36hr_LB_vs_DA"
colnames(Dermo_dds_36h_res_LFC_sig_id)[1] <- "Transcript_id"

Dermo_dds_7d_res_LFC_sig_id <- as.data.frame(Dermo_dds_7d_res_LFC_sig_id)
Dermo_dds_7d_res_LFC_sig_id$Comparison <- "Dermo_7d_LB_vs_DA"
colnames(Dermo_dds_7d_res_LFC_sig_id)[1] <- "Transcript_id"

Dermo_dds_DA_time_res_early_middle_LFC_sig_id <- as.data.frame(Dermo_dds_DA_time_res_early_middle_LFC_sig_id)
Dermo_dds_DA_time_res_early_middle_LFC_sig_id$Comparison <- "Dermo_DA_6hr_vs_36hr"
colnames(Dermo_dds_DA_time_res_early_middle_LFC_sig_id)[1] <- "Transcript_id"

Dermo_dds_DA_time_res_early_late_LFC_sig_id <- as.data.frame(Dermo_dds_DA_time_res_early_late_LFC_sig_id)
Dermo_dds_DA_time_res_early_late_LFC_sig_id$Comparison <- "Dermo_DA_6hr_vs_7d"
colnames(Dermo_dds_DA_time_res_early_late_LFC_sig_id)[1] <- "Transcript_id"

Dermo_dds_LB_time_res_early_middle_LFC_sig_id <- as.data.frame(Dermo_dds_LB_time_res_early_middle_LFC_sig_id)
Dermo_dds_LB_time_res_early_middle_LFC_sig_id$Comparison <- "Dermo_LB_6hr_vs_36hr"
colnames(Dermo_dds_LB_time_res_early_middle_LFC_sig_id)[1] <- "Transcript_id"

Dermo_dds_LB_time_res_early_late_LFC_sig_id <- as.data.frame(Dermo_dds_LB_time_res_early_late_LFC_sig_id)
Dermo_dds_LB_time_res_early_late_LFC_sig_id$Comparison <- "Dermo_LB_6hr_vs_7d"
colnames(Dermo_dds_LB_time_res_early_late_LFC_sig_id)[1] <- "Transcript_id"

Dermo_dds_DA_time_res_middle_late_LFC_sig_id <- as.data.frame(Dermo_dds_DA_time_res_middle_late_LFC_sig_id)
Dermo_dds_DA_time_res_middle_late_LFC_sig_id$Comparison <- "Dermo_DA_36hr_vs_7d"
colnames(Dermo_dds_DA_time_res_middle_late_LFC_sig_id)[1] <- "Transcript_id"

Dermo_dds_LB_time_res_middle_late_LFC_sig_id <- as.data.frame(Dermo_dds_LB_time_res_middle_late_LFC_sig_id)
Dermo_dds_LB_time_res_middle_late_LFC_sig_id$Comparison <- "Dermo_LB_36hr_vs_7d"
colnames(Dermo_dds_LB_time_res_middle_late_LFC_sig_id)[1] <- "Transcript_id"

# Combine the data frames using rbind()
upset_all_sig <- rbind(Dermo_dds_family_res_LFC_sig_id, Dermo_dds_6h_res_LFC_sig_id,  Dermo_dds_36h_res_LFC_sig_id,
                       Dermo_dds_7d_res_LFC_sig_id, Dermo_dds_DA_time_res_early_middle_LFC_sig_id, Dermo_dds_DA_time_res_early_late_LFC_sig_id,
                         Dermo_dds_LB_time_res_early_middle_LFC_sig_id, Dermo_dds_LB_time_res_early_late_LFC_sig_id, Dermo_dds_DA_time_res_middle_late_LFC_sig_id, 
                       Dermo_dds_LB_time_res_middle_late_LFC_sig_id )

# Convert into wide format using reshape
upset_all_sig_wide <- upset_all_sig %>% mutate(value=1) %>% spread(Comparison, value, fill =0 )
head(upset_all_sig_wide)

# Make upset plot
upset(upset_all_sig_wide,  mainbar.y.label = "Transcript id Intersections", 
      sets.x.label = "# Significantly Differentially Expressed Transcripts", text.scale = c(1.3, 1.3, 1.3, 1.3, 2, 1.0),
      sets= c("Family_Comparison_all_timepoints","6hr_LB_vs_DA", "36hr_LB_vs_DA","7d_LB_vs_DA","DA_6hr_vs_36hr","DA_6hr_vs_7d",
      "DA_36hr_vs_7d","LB_6hr_vs_36hr","LB_6hr_vs_7d","LB_36hr_vs_7d"), order.by="freq")

# Make upset plot prettier with complex heatmap

##### Isolate Annotated Genes from Different Comparisons #####
# Find matches between groups below in their genes using inner join of transcript id 
# how can I export interactions from my upset plots?
# inner join won't work because though it will show me common things, it won't show me what in unique as compared to all other comparisons

# Code written at https://github.com/hms-dbmi/UpSetR/issues/85 to help get the list of all intersections 
# to format for code I need a binary table from my data first

# example data set
movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
                   header = T, sep = ";")

upset_all_sig_binary <- upset_all_sig %>% mutate(value=1) %>% spread(Comparison, value, fill=0)
class(upset_all_sig_binary )
row.names(upset_all_sig_binary) <- upset_all_sig_binary$Transcript_id

# Get intersections from Binary table using the following function formatted for numerical input 
# get_intersect_members() takes as arguments a dataframe that has been formatted as a binary 
#     table, such as movies from the UpSetR vignette; as well as a series of strings with the names of 
#     columns you wish to test membership for.
get_intersect_members <- function (x, ...){
  require(dplyr)
  require(tibble)
  # the following makes sure that we don't have any weird values in the dataframe
  x <- x[,sapply(x, is.numeric)][,0<=colMeans(x[,sapply(x, is.numeric)],na.rm=T) & colMeans(x[,sapply(x, is.numeric)],na.rm=T)<=1]
  n <- names(x)
  #convert rownames to a column to prevent mulching by tidyr
  x %>% rownames_to_column() -> x
  l <- c(...)
  a <- intersect(names(x), l)
  ar <- vector('list',length(n)+1)
  ar[[1]] <- x
  i=2
  for (item in n) {
    if (item %in% a){
      if (class(x[[item]])=='numeric'){   #Now uses numeric instead of integer
        ar[[i]] <- paste(item, '>= 1')
        i <- i + 1
      }
    } else {
      if (class(x[[item]])=='numeric'){
        ar[[i]] <- paste(item, '== 0')
        i <- i + 1
      }
    }
  }
  do.call(filter_, ar) %>% column_to_rownames() -> x
  return(x)
}


# Comparison 1 6hr_LB_vs_DA, 7d_LB_vs_DA, 36hr_LB_vs_DA 
comp1 <- get_intersect_members(upset_all_sig_binary, "Dermo_6hr_LB_vs_DA", "Dermo_36hr_LB_vs_DA","Dermo_7d_LB_vs_DA")
comp1_id <- row.names(comp1)  
comp1_id <- as.data.frame(comp1_id)
colnames(comp1_id)[1] <- "transcript_id"
comp1_annot <- left_join(comp1_id, C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product","gene")], by = "transcript_id")

# Comparison 2 6hr_LB_vs_DA, 7d_LB_vs_DA
comp2 <- get_intersect_members(upset_all_sig_binary, "Dermo_6hr_LB_vs_DA","Dermo_7d_LB_vs_DA")
comp2_id <- row.names(comp2)  
comp2_id <- as.data.frame(comp2_id)
colnames(comp2_id)[1] <- "transcript_id"
comp2_annot <- left_join(comp2_id, C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product","gene")], by = "transcript_id")

# Comparison 3 LB_6hr_vs_7d, LB_6hr_vs_36hr

comp3 <- get_intersect_members(upset_all_sig_binary, "Dermo_LB_6hr_vs_7d","Dermo_LB_6hr_vs_36hr")
comp3_id <- row.names(comp3)  
comp3_id <- as.data.frame(comp3_id)
colnames(comp3_id)[1] <- "transcript_id"
comp3_annot <- left_join(comp3_id, C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product","gene")], by = "transcript_id")
nrow(comp3_annot) # 1298, these are the genes shared between both timepoint comparisons
# Heatmaps for specific comparisons of vsd data 

comp3_annot_APOP <- comp3_annot[grepl(paste(Apoptosis_names,collapse="|"), 
                                      comp3_annot$product, ignore.case = TRUE),]
nrow(comp3_annot_APOP ) # 35

# These are the significantly differentially expressed genes unique to the LB 6hr vs 36hr nd 6hr vs 7d comparison. 

##### Gene Clustering Analysis Heatmaps ####

# plot 40 genes with the highest variance across samples for each comparison, we will work with the vsd data for this
# then we will generate a heatmap with the 40 most variable VST-transformed genes 

# Subset the top 40 variable genes for plotting 

# example codes from RNAseq workflow: https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#other-comparisons
  # topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 40)
  # mat  <- assay(vsd)[ topVarGenes, ]
  # mat  <- mat - rowMeans(mat)
  # anno <- as.data.frame(colData(vsd)[, c("cell","dex")])
  # pheatmap(mat, annotation_col = anno)

# Family all timepoints
topVarGenes_Dermo_dds_family_vsd <-  head(order(rowVars(assay(Dermo_dds_family_vsd )), decreasing = TRUE), 40)
family_res_mat <- assay(Dermo_dds_family_vsd)[topVarGenes_Dermo_dds_family_vsd,]
family_res_mat <- family_res_mat - rowMeans(family_res_mat)
family_res_anno <- as.data.frame(colData(Dermo_dds_family_vsd)[, c("FamCode", "Timepoint")])
family_heatmap <- pheatmap(family_res_mat, annotation_col = family_res_anno)
head(family_res_mat)
# reorder annotation table to match ordering in heatmap 
family_heatmap_reorder <-rownames(family_res_mat[family_heatmap$tree_row[["order"]],])
# annotate the row.names
family_res_mat_prot <- as.data.frame(family_heatmap_reorder)
colnames(family_res_mat_prot)[1] <- "transcript_id"
family_res_mat_prot_annot <- family_res_mat_prot %>% left_join(select(C_vir_rtracklayer_transcripts_GO, transcript_id, product, gene), by = "transcript_id")
#isolate interesting clusters

# LB vs. DA cluster
LB_DA_cluster <- c("XM_022479343.1","XM_022473887.1","XM_022477354.1","XM_022477356.1")
LB_DA_cluster <- as.data.frame(LB_DA_cluster)
LB_DA_cluster <- rename(LB_DA_cluster, "transcript_id"=LB_DA_cluster)
LB_DA_cluster_subset <- subset(family_res_mat_prot_annot, transcript_id %in% LB_DA_cluster$transcript_id)
LB_DA_cluster_subset
# grab 10 clusters assigned by pheatmap using cutree
family_all_clusters <- cbind(family_res_mat, cluster=cutree(family_heatmap$tree_row, k=10))

# 6hr
topVarGenes_Dermo_dds_6h_vsd <- head(order(rowVars(assay(Dermo_dds_6h_vsd)), decreasing = TRUE), 40)
Res_mat_6hr <- assay(Dermo_dds_6h_vsd)[topVarGenes_Dermo_dds_6h_vsd ,]
Res_mat_6hr <- Res_mat_6hr - rowMeans(Res_mat_6hr)
Res_mat_6hr_anno <- as.data.frame(colData(Dermo_dds_6h_vsd)[, c("FamCode","Timepoint")])
six_hr_heatmap <- pheatmap(Res_mat_6hr, annotation_col = Res_mat_6hr_anno)
# reorder annotation table to match ordering in heatmap 
six_hr_heatmap_reorder <-rownames(Res_mat_6hr[six_hr_heatmap$tree_row[["order"]],])
# annotate the row.names
Res_mat_6hr_prot <- as.data.frame(six_hr_heatmap_reorder)
colnames(Res_mat_6hr_prot)[1] <- "transcript_id"
Res_mat_6hr_prot_annot <- Res_mat_6hr_prot %>% left_join(select(C_vir_rtracklayer_transcripts_GO, transcript_id, product, gene), by = "transcript_id")
#isolate interesting clusters
six_hr_comparison_cluster <- c("XM_022455505.1", "XM_022484575.1", "XM_022461506.1", "XM_022430618.1", "XM_022490512.1",
                              "XM_022461508.1", "XM_022464459.1", "XM_022483469.1", "XM_022483473.1", "XM_022442223.1", "XM_022457463.1", "XM_022442224.1")
six_hr_comparison_cluster <- as.data.frame(six_hr_comparison_cluster)
colnames(six_hr_comparison_cluster)[1] <- "transcript_id"
six_hr_comparison_cluster_subset <- subset(Res_mat_6hr_prot_annot, transcript_id %in% six_hr_comparison_cluster$transcript_id)

#36hr
topVarGenes_Dermo_dds_36h_vsd <- head(order(rowVars(assay(Dermo_dds_36h_vsd)), decreasing = TRUE), 40)
Res_mat_36hr <- assay(Dermo_dds_36h_vsd)[topVarGenes_Dermo_dds_36h_vsd ,]
Res_mat_36hr <- Res_mat_36hr - rowMeans(Res_mat_36hr)
Res_mat_36hr_anno <- as.data.frame(colData(Dermo_dds_36h_vsd)[, c("FamCode","Timepoint")])
thirty_six_hr_heatmap <- pheatmap(Res_mat_36hr, annotation_col = Res_mat_36hr_anno)
# reorder annotation table to match ordering in heatmap 
thirty_six_hr_heatmap_reorder <-rownames(Res_mat_36hr[thirty_six_hr_heatmap$tree_row[["order"]],])
# annotate the row.names
Res_mat_36hr_prot <- as.data.frame(thirty_six_hr_heatmap_reorder)
colnames(Res_mat_36hr_prot)[1] <- "transcript_id"
Res_mat_36hr_prot_annot <- Res_mat_36hr_prot %>% left_join(select(C_vir_rtracklayer_transcripts_GO, transcript_id, product, gene), by = "transcript_id")
#isolate interesting clusters
thirty_six_cluster <- c("XM_022474564.1", "XM_022473887.1", "XM_022477354.1", "XM_022477356.1")
thirty_six_cluster <- as.data.frame(thirty_six_cluster)
colnames(thirty_six_cluster)[1] <- "transcript_id"
thirty_six_cluster_subset <- subset(Res_mat_36hr_prot_annot , transcript_id %in% thirty_six_cluster$transcript_id)

#7d # by 7 days they arent clustering out like we would like 
topVarGenes_Dermo_dds_7d_vsd <- head(order(rowVars(assay(Dermo_dds_7d_vsd)), decreasing = TRUE), 40)
Res_mat_7d <- assay(Dermo_dds_7d_vsd)[topVarGenes_Dermo_dds_7d_vsd ,]
Res_mat_7d <- Res_mat_7d - rowMeans(Res_mat_7d)
Res_mat_7d_anno <- as.data.frame(colData(Dermo_dds_7d_vsd)[, c("FamCode","Timepoint")])
seven_d_heatmap <- pheatmap(Res_mat_7d, annotation_col = Res_mat_7d_anno)
# reorder annotation table to match ordering in heatmap 
seven_d_heatmap_reorder <-rownames(Res_mat_7d[seven_d_heatmap$tree_row[["order"]],])

#LB_time can do later if needed
#topVarGenes_Dermo_dds_LB_time_res_early_late_LFC_sig_vsd   <- head(order(rowVars(assay(Dermo_dds_LB_time_res_early_late_LFC_sig_vsd)), decreasing = TRUE), 40)
#topVarGenes_Dermo_dds_LB_time_res_early_middle_LFC_sig_vsd <- head(order(rowVars(assay(Dermo_dds_LB_time_res_early_middle_LFC_sig_vsd)), decreasing = TRUE), 40)
#topVarGenes_Dermo_dds_LB_time_res_middle_late_LFC_sig_vsd  <- head(order(rowVars(assay(Dermo_dds_LB_time_res_middle_late_LFC_sig_vsd)), decreasing = TRUE), 40)
#topVarGenes_Dermo_dds_DA_time_res_early_late_LFC_sig_vsd   <- head(order(rowVars(assay(Dermo_dds_DA_time_res_early_late_LFC_sig_vsd )), decreasing = TRUE), 40)
#topVarGenes_Dermo_dds_DA_time_res_early_middle_LFC_sig_vsd <- head(order(rowVars(assay(Dermo_dds_DA_time_res_early_middle_LFC_sig_vsd)), decreasing = TRUE), 40)
#topVarGenes_Dermo_dds_DA_time_res_middle_late_LFC_sig_vsd  <- head(order(rowVars(assay(Dermo_dds_DA_time_res_middle_late_LFC_sig_vsd)), decreasing = TRUE), 40)

####  Apoptosis Gene List Upload ####

Apoptosis_names <- c('bcl-2-related protein A1',
                     'apoptosis-inducing factor 1',
                     'RAC-alpha',
                     'RAC-gamma serine/threonine-protein kinase-',
                     'APAF-1 interacting protein',
                     'tumor necrosis factor ligand superfamily member 10',
                     'tumor necrosis factor superfamily member 12',
                     'Aven',
                     'BCL2 associated agonist of cell death',
                     'BAG family molecular chaperone regulator 1',
                     'BAG family molecular chaperone regulator 2',
                     'BAG family molecular chaperone regulator 3',
                     'bcl-2 homologous antagonist/killer',
                     'apoptosis regulator BAX',
                     'bcl-2-like protein 2',
                     'bcl-2-like protein 1',
                     'bcl-2',
                     'bcl2 modifying factor',
                     'Bax inhibitor 1',
                     'BH3 interacting domain death agonist',
                     'Bcl-2 interacting killer',
                     'bcl-2 interacting protein BIM',
                     'Bik-like killer protein',
                     'Bcl-2 related ovarian killer',
                     'CASP8 and FADD like apoptosis regulator',
                     'transcription factor AP-1',
                     'caspase activity and apoptosis inhibitor 1',
                     'DNA fragmentation factor subunit beta',
                     'adenylate cyclase',
                     'caspase-1',
                     'caspase-2',
                     'caspase-3, caspase-7',
                     'caspase-8',
                     'caspase-9',
                     'caspase-10',
                     'caspase-11',
                     'caspase-6',
                     'caspase-4',
                     'caspase-5',
                     'cell division cycle and apoptosis regulator protein 1',
                     'CD151 antigen',
                     'protein BTG1',
                     'interferon alpha',
                     'adenylate cyclase 10',
                     'caspase activity and apoptosis inhibitor 1',
                     'baculoviral IAP repeat-containing protein 2',
                     'baculoviral IAP repeat-containing protein 3',
                     'cAMP-responsive element modulator',
                     'cytochrome c',
                     'death-associated inhibitor of apoptosis 2',
                     'tumor necrosis factor receptor superfamily member 25',
                     'tumor necrosis factor receptor superfamily member 10A',
                     'tumor necrosis factor receptor superfamily member 10B',
                     'endonuclease G',
                     'FAS-associated death domain protein',
                     'fas apoptotic inhibitory molecule 1',
                     'tumor necrosis factor receptor superfamily member 6',
                     'fas cell surface death receptor',
                     'gasdermin E',
                     'GTPase IMAP family member 4',
                     'GTPase IMAP family member 7',
                     'GTPase IMAP family 8',
                     'harakiri',
                     'baculoviral IAP repeat-containing protein 5',
                     'baculoviral IAP repeat-containing protein 6',
                     'baculoviral IAP repeat-containing protein 1',
                     'baculoviral IAP repeat-containing protein 8',
                     'baculoviral IAP repeat-containing protein 7',
                     'DNA fragmentation factor subunit alpha',
                     'interferon-induced protein 44',
                     'NF-kappa-B inhibitor alpha',
                     'NF-kappa-B inhibitor epsilon',
                     'inositol 1,4,5-trisphosphate receptor',
                     'stress-activated protein kinase JNK',
                     'lipopolysaccharide-induced tumor necrosis factor-alpha',
                     'induced myeloid leukemia cell differentiation protein Mcl-1-like',
                     'mitogen-activated protein kinase kinase kinase 1',
                     'mitogen-activated protein kinase 1',
                     'mitogen-activated protein kinase kinase kinase 7',
                     'MYC proto-oncogene',
                     'myeloid differentiation primary response protein MyD88',
                     'phorbol-12-myristate-13-acetate-induced protein 1',
                     'nuclear factor NF-kappa-B p105 subunit',
                     'nuclear factor NF-kappa-B p100 subunit',
                     'transcription factor p65',
                     'RELB proto-oncogene, NF-kB subunit',
                     'reticuloendotheliosis oncogene',
                     'anti-apoptotic protein NR13',
                     'nuclear mitotic apparatus protein 1',
                     'dynamin-like 120 kDa protein',
                     'Early 35 kDa protein',
                     'mitogen-activated protein kinase 11`',
                     'mitogen-activated protein kinase 12',
                     'mitogen-activated protein kinase 13',
                     'mitogen-activated protein kinase 14A',
                     'cellular tumor antigen p53',
                     'programmed cell death protein',
                     'programmed cell death protein 1',
                     'p53 and DNA damage-regulated protein 1',
                     'phosphatidylinositol 3-kinase',
                     'cAMP-dependent protein kinase',
                     'protein kinase C',
                     'BCL2 binding component 3',
                     'cdc42',
                     'ras-like GTP-binding protein rhoA',
                     'ras-like GTP-binding protein RHO',
                     'rho-related GTP-binding protein RhoE-like',
                     'ras-related C3 botulinum toxin substrate 1',
                     'rho-related protein racA',
                     'ras-like GTP-binding protein Rho1',
                     'mitochondrial Rho GTPase 1',
                     'receptor-interacting serine/threonine-protein kinase 1',
                     'receptor-interacting serine/threonine-protein kinase 4',
                     'diablo homolog, mitochondrial',
                     'toll-like receptor',
                     'tumor necrosis factor',
                     'lymphotoxin-alpha',
                     'tumor necrosis factor receptor superfamily member 1A',
                     'CD40 ligand',
                     'tumor necrosis factor receptor superfamily member',
                     'TNFRSF1A associated via death domain',
                     'TNF receptor-associated factor 2',
                     'TNF receptor-associated factor',
                     'E3 ubiquitin-protein ligase XIAP',
                     'netrin receptor DCC',
                     'netrin receptor UNC5A',
                     'netrin receptor UNC5B',
                     'netrin receptor UNC5C',
                     'netrin receptor UNC5D',
                     'neurotrophic receptor tyrosine kinase 1',
                     'sonic hedgehog receptor',
                     'peptidyl-prolyl cis-trans isomerase',
                     'receptor-interacting serine/threonine-protein kinase',
                     'mixed lineage kinase domain',
                     'HSP90',
                     'E3 ubiquitin-protein ligase CHIP',
                     'tumor necrosis factor alpha-induced protein 3',
                     'protein phosphatase 1B',
                     'aurora kinase A',
                     'glutathione peroxidase 4',
                     'gasdermin',
                     'poly [ADP-ribose] polymerase 1',
                     'macrophage migration inhibitory factor',
                     'hexokinase-1',
                     'Raf-1 protooncogene serine/threonine kinase',
                     'elastase, neutrophil expressed',
                     'cathepsin',
                     'PRKC apoptosis WT1 regulator protein',
                     'apoptosis-stimulating of p53 protein 1',
                     'apoptosis-stimulating of p53 protein 2',
                     'apoptosis inhibitory protein 5',
                     'apoptotic chromatin condensation inducer in the nucleus',
                     'interleukin 17-like protein',
                     'heat shock',
                     'calreticulin',
                     'serine protease inhibitor') # made some additions of interleukins 

#### Extract list of significant Apoptosis Genes ####
Dermo_dds_family_res_LFC_sig_annot_APOP <- Dermo_dds_family_res_LFC_sig_annot[grepl(paste(Apoptosis_names,collapse="|"), 
      Dermo_dds_family_res_LFC_sig_annot$product, ignore.case = TRUE),]
arrange(Dermo_dds_family_res_LFC_sig_annot_APOP, -log2FoldChange) 
nrow(Dermo_dds_family_res_LFC_sig_annot_APOP) # 78

Dermo_dds_6h_res_LFC_sig_annot_APOP <- Dermo_dds_6h_res_LFC_sig_annot[grepl(paste(Apoptosis_names,collapse="|"), 
Dermo_dds_6h_res_LFC_sig_annot$product, ignore.case = TRUE),]
arrange(Dermo_dds_6h_res_LFC_sig_annot_APOP , -log2FoldChange) 
nrow(Dermo_dds_6h_res_LFC_sig_annot_APOP) # 25

#The 36hr and response are perhaps not very different from one another, 146 rows of genes
Dermo_dds_36h_res_LFC_sig_annot_APOP <- Dermo_dds_36h_res_LFC_sig_annot[grepl(paste(Apoptosis_names,collapse="|"), 
Dermo_dds_36h_res_LFC_sig_annot$product, ignore.case = TRUE),]
arrange(Dermo_dds_36h_res_LFC_sig_annot_APOP , -log2FoldChange) 
nrow(Dermo_dds_36h_res_LFC_sig_annot_APOP) # 102

#94
Dermo_dds_7d_res_LFC_sig_annot_APOP <- Dermo_dds_7d_res_LFC_sig_annot[grepl(paste(Apoptosis_names,collapse="|"), 
Dermo_dds_7d_res_LFC_sig_annot$product, ignore.case = TRUE),]
arrange(Dermo_dds_7d_res_LFC_sig_annot_APOP  , -log2FoldChange) 
nrow(Dermo_dds_7d_res_LFC_sig_annot_APOP) #94

#DA comparisons (Capturing more the changes through time), #103 rows
Dermo_dds_DA_time_res_early_late_LFC_sig_annot_APOP <- Dermo_dds_DA_time_res_early_late_LFC_sig_annot[grepl(paste(Apoptosis_names,collapse="|"), 
Dermo_dds_DA_time_res_early_late_LFC_sig_annot$product, ignore.case = TRUE),]
arrange(Dermo_dds_DA_time_res_early_late_LFC_sig_annot_APOP  , -log2FoldChange) 
nrow(Dermo_dds_DA_time_res_early_late_LFC_sig_annot_APOP) #97

#106
Dermo_dds_DA_time_res_early_middle_LFC_sig_annot_APOP <- Dermo_dds_DA_time_res_early_middle_LFC_sig_annot[grepl(paste(Apoptosis_names,collapse="|"), 
Dermo_dds_DA_time_res_early_middle_LFC_sig_annot$product, ignore.case = TRUE),]
arrange(Dermo_dds_DA_time_res_early_middle_LFC_sig_annot_APOP  , -log2FoldChange) 
nrow(Dermo_dds_DA_time_res_early_middle_LFC_sig_annot_APOP) #87

#4, reiterating that the major changes happening are not between 36 and 7d, its early 
Dermo_dds_DA_time_res_middle_late_LFC_sig_annot_APOP <- Dermo_dds_DA_time_res_middle_late_LFC_sig_annot[grepl(paste(Apoptosis_names,collapse="|"), 
Dermo_dds_DA_time_res_middle_late_LFC_sig_annot$product, ignore.case = TRUE),]
arrange(Dermo_dds_DA_time_res_middle_late_LFC_sig_annot_APOP  , -log2FoldChange) 
nrow(Dermo_dds_DA_time_res_middle_late_LFC_sig_annot_APOP ) # 4

##LB comparisons
# 126 genes
Dermo_dds_LB_time_res_early_late_LFC_sig_annot_APOP <- Dermo_dds_LB_time_res_early_late_LFC_sig_annot[grepl(paste(Apoptosis_names,collapse="|"), 
Dermo_dds_LB_time_res_early_late_LFC_sig_annot$product, ignore.case = TRUE),]
arrange(Dermo_dds_LB_time_res_early_late_LFC_sig_annot_APOP  , -log2FoldChange) 
nrow(Dermo_dds_LB_time_res_early_late_LFC_sig_annot_APOP) 

# 136
Dermo_dds_LB_time_res_early_middle_LFC_sig_annot_APOP <- Dermo_dds_LB_time_res_early_middle_LFC_sig_annot[grepl(paste(Apoptosis_names,collapse="|"), 
Dermo_dds_LB_time_res_early_middle_LFC_sig_annot$product, ignore.case = TRUE),]
arrange(Dermo_dds_LB_time_res_early_middle_LFC_sig_annot_APOP  , -log2FoldChange) 
nrow(Dermo_dds_LB_time_res_early_middle_LFC_sig_annot_APOP)

# 0 genes in the middle to late comparison 
Dermo_dds_LB_time_res_middle_late_LFC_sig_annot_APOP <- Dermo_dds_LB_time_res_middle_late_LFC_sig_annot[grepl(paste(Apoptosis_names,collapse="|"), 
Dermo_dds_LB_time_res_middle_late_LFC_sig_annot$product, ignore.case = TRUE),]
arrange(Dermo_dds_LB_time_res_middle_late_LFC_sig_annot_APOP  , -log2FoldChange) 

### EARLY responses are the most critical!! Major changes in either family are not happening between 36hr and 7days.
  # We are still capturing changes by looking at the two days, but the 36hr response is much more interesting

#### Apoptosis Data Visualizations####

### Upset Plot of Apoptosis Gene Data 

Dermo_dds_family_res_LFC_sig_annot_APOP_id <- Dermo_dds_family_res_LFC_sig_annot_APOP$transcript_id
Dermo_dds_family_res_LFC_sig_annot_APOP_id <- as.data.frame(Dermo_dds_family_res_LFC_sig_annot_APOP_id)
Dermo_dds_family_res_LFC_sig_annot_APOP_id$Comparison <- "All_family_Res"
colnames(Dermo_dds_family_res_LFC_sig_annot_APOP_id)[1] <- "Transcript_id"

Dermo_dds_6h_res_LFC_sig_annot_APOP_id <- Dermo_dds_6h_res_LFC_sig_annot_APOP$transcript_id
Dermo_dds_6h_res_LFC_sig_annot_APOP_id <- as.data.frame(Dermo_dds_6h_res_LFC_sig_annot_APOP_id)
Dermo_dds_6h_res_LFC_sig_annot_APOP_id$Comparison <- "6hr_LB_vs_DA"
colnames(Dermo_dds_6h_res_LFC_sig_annot_APOP_id)[1] <- "Transcript_id"

Dermo_dds_36h_res_LFC_sig_annot_APOP_id <-  Dermo_dds_36h_res_LFC_sig_annot_APOP$transcript_id
Dermo_dds_36h_res_LFC_sig_annot_APOP_id <- as.data.frame(Dermo_dds_36h_res_LFC_sig_annot_APOP_id)
Dermo_dds_36h_res_LFC_sig_annot_APOP_id$Comparison <- "36hr_LB_vs_DA"
colnames(Dermo_dds_36h_res_LFC_sig_annot_APOP_id)[1] <- "Transcript_id"

Dermo_dds_7d_res_LFC_sig_annot_APOP_id <- Dermo_dds_7d_res_LFC_sig_annot_APOP$transcript_id
Dermo_dds_7d_res_LFC_sig_annot_APOP_id <- as.data.frame(Dermo_dds_7d_res_LFC_sig_annot_APOP_id)
Dermo_dds_7d_res_LFC_sig_annot_APOP_id$Comparison <- "7d_LB_vs_DA"
colnames(Dermo_dds_7d_res_LFC_sig_annot_APOP_id)[1] <-  "Transcript_id"

Dermo_dds_DA_time_res_early_late_LFC_sig_annot_APOP_id <- Dermo_dds_DA_time_res_early_late_LFC_sig_annot_APOP$transcript_id
Dermo_dds_DA_time_res_early_late_LFC_sig_annot_APOP_id <- as.data.frame(Dermo_dds_DA_time_res_early_late_LFC_sig_annot_APOP_id)
Dermo_dds_DA_time_res_early_late_LFC_sig_annot_APOP_id$Comparison <- "DA_6hr_vs_7d"
colnames(Dermo_dds_DA_time_res_early_late_LFC_sig_annot_APOP_id)[1] <- "Transcript_id"

Dermo_dds_DA_time_res_early_middle_LFC_sig_annot_APOP_id <- Dermo_dds_DA_time_res_early_middle_LFC_sig_annot_APOP$transcript_id
Dermo_dds_DA_time_res_early_middle_LFC_sig_annot_APOP_id <- as.data.frame(Dermo_dds_DA_time_res_early_middle_LFC_sig_annot_APOP_id)
Dermo_dds_DA_time_res_early_middle_LFC_sig_annot_APOP_id$Comparison <- "DA_6hr_vs_36hr"
colnames(Dermo_dds_DA_time_res_early_middle_LFC_sig_annot_APOP_id)[1] <- "Transcript_id"

Dermo_dds_DA_time_res_middle_late_LFC_sig_annot_APOP_id <- Dermo_dds_DA_time_res_middle_late_LFC_sig_annot_APOP$transcript_id
Dermo_dds_DA_time_res_middle_late_LFC_sig_annot_APOP_id <- as.data.frame(Dermo_dds_DA_time_res_middle_late_LFC_sig_annot_APOP_id)
Dermo_dds_DA_time_res_middle_late_LFC_sig_annot_APOP_id$Comparison <- "DA_36hr_vs_7d"
colnames(Dermo_dds_DA_time_res_middle_late_LFC_sig_annot_APOP_id)[1] <- "Transcript_id"

Dermo_dds_LB_time_res_early_late_LFC_sig_annot_APOP_id <- Dermo_dds_LB_time_res_early_late_LFC_sig_annot_APOP$transcript_id
Dermo_dds_LB_time_res_early_late_LFC_sig_annot_APOP_id <- as.data.frame(Dermo_dds_LB_time_res_early_late_LFC_sig_annot_APOP_id)
Dermo_dds_LB_time_res_early_late_LFC_sig_annot_APOP_id$Comparison <-"LB_6hr_vs_7d"
colnames(Dermo_dds_LB_time_res_early_late_LFC_sig_annot_APOP_id)[1] <-"Transcript_id"

Dermo_dds_LB_time_res_early_middle_LFC_sig_annot_APOP_id <- Dermo_dds_LB_time_res_early_middle_LFC_sig_annot_APOP$transcript_id
Dermo_dds_LB_time_res_early_middle_LFC_sig_annot_APOP_id <- as.data.frame(Dermo_dds_LB_time_res_early_middle_LFC_sig_annot_APOP_id)
Dermo_dds_LB_time_res_early_middle_LFC_sig_annot_APOP_id$Comparison <- "LB_6hr_vs_36hr" 
colnames(Dermo_dds_LB_time_res_early_middle_LFC_sig_annot_APOP_id)[1] <- "Transcript_id"

# No data here
Dermo_dds_LB_time_res_middle_late_LFC_sig_annot_APOP_id <- Dermo_dds_LB_time_res_middle_late_LFC_sig_annot_APOP$transcript_id
Dermo_dds_LB_time_res_middle_late_LFC_sig_annot_APOP_id <- as.data.frame(Dermo_dds_LB_time_res_middle_late_LFC_sig_annot_APOP_id)
Dermo_dds_LB_time_res_middle_late_LFC_sig_annot_APOP_id$Comparison <- "LB_36hr_vs_7d"

# Combine the data frames using rbind()
upset_all_sig_APOP <- rbind(Dermo_dds_family_res_LFC_sig_annot_APOP_id, Dermo_dds_6h_res_LFC_sig_annot_APOP_id,  Dermo_dds_36h_res_LFC_sig_annot_APOP_id,
                       Dermo_dds_7d_res_LFC_sig_annot_APOP_id, Dermo_dds_DA_time_res_early_middle_LFC_sig_annot_APOP_id, Dermo_dds_DA_time_res_early_late_LFC_sig_annot_APOP_id,
                       Dermo_dds_LB_time_res_early_middle_LFC_sig_annot_APOP_id, Dermo_dds_LB_time_res_early_late_LFC_sig_annot_APOP_id, Dermo_dds_DA_time_res_middle_late_LFC_sig_annot_APOP_id, 
                       Dermo_dds_LB_time_res_middle_late_LFC_sig_annot_APOP_id )

# Convert into wide format using reshape
upset_all_sig_wide_APOP <- upset_all_sig_APOP %>% mutate(value=1) %>% spread(Comparison, value, fill =0 )
head(upset_all_sig_wide_APOP)

# Make upset plot
upset(upset_all_sig_wide_APOP,  mainbar.y.label = "Apoptosis Pathway Transcript id Intersections", 
      sets.x.label = "# Significantly Differentially Expressed Transcripts", text.scale = c(1.3, 1.3, 1, 1, 2, 1.3),
      sets= c("All_family_Res","6hr_LB_vs_DA", "36hr_LB_vs_DA","7d_LB_vs_DA","DA_6hr_vs_36hr","DA_6hr_vs_7d",
              "DA_36hr_vs_7d","LB_6hr_vs_36hr","LB_6hr_vs_7d"), order.by="degree")

upset(upset_all_sig_wide_APOP,  mainbar.y.label = "Apoptosis Pathway Transcript id Intersections", 
      sets.x.label = "# Significantly Differentially Expressed Transcripts", text.scale = c(1.3, 1.3, 1.3, .9, 2, 1.3),
      sets= c("All_family_Res","6hr_LB_vs_DA", "36hr_LB_vs_DA","7d_LB_vs_DA","DA_6hr_vs_36hr","DA_6hr_vs_7d",
              "DA_36hr_vs_7d","LB_6hr_vs_36hr","LB_6hr_vs_7d"), order.by="freq")


# Heatmap of vst count data of apoptosis genes

# Use sig_annot_APOP data to subset the vsd dataframe for each comparison, then plot heatmap of it for each column 

#### TIMECOURSE ANALYSIS (all timepoints, LB vs. DA) ####

# The following chunk of code performs a likelihood ratio test,
# where we remove the strain-specific differences over time. 
# Genes with small p values from this test are those which at one 
# or more time points after time 0 showed a strain-specific effect. 
# Note therefore that this will not give small p values to genes that moved up or down over time in the same way in both strains

# An LRT tests all the changes across time while the Wald test is for a specific time point. Run the LRT first
Dermo_dds_timecourse_family_deseq <-  DESeq(Dermo_dds_timecourse_family, test="LRT", reduced=~FamCode+Timepoint)
Dermo_dds_timecourse_family_res <- results(Dermo_dds_timecourse_family_deseq, alpha=0.05)
Dermo_dds_timecourse_family_res$symbol <- mcols(Dermo_dds_timecourse_family_deseq)$symbol
head(Dermo_dds_timecourse_family_res[order(Dermo_dds_timecourse_family_res$padj),], 4)
resultsNames(Dermo_dds_timecourse_family_res)

# Wald tests for the log2 fold changes at individual time points can be investigated using the test argument to results
resultsNames(Dermo_dds_timecourse_family_deseq) # "Intercept", "FamCode_LB_vs_DA", "Timepoint_6h_vs_36h", "Timepoint_7d_vs_36h", "FamCodeLB.Timepoint6h", "FamCodeLB.Timepoint7d" "Treat_Injected_vs_Control"

# Perform LFC shrinkage with apeglm for contrasts
resTC_LB_vs_DA <- results(Dermo_dds_timecourse_family_deseq, name="FamCode_LB_vs_DA", test="Wald")
resTC_LB_vs_DA_LFC <- lfcShrink(Dermo_dds_timecourse_family_deseq, coef="FamCode_LB_vs_DA", res=resTC_LB_vs_DA, type="apeglm")

# Annotate full list
resTC_LB_vs_DA_LFC$transcript_id <- row.names(resTC_LB_vs_DA_LFC)
resTC_LB_vs_DA_LFC <- as.data.frame(resTC_LB_vs_DA_LFC)
resTC_LB_vs_DA_LFC_annot <- left_join(resTC_LB_vs_DA_LFC, C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product", "unique_go")], by ="transcript_id")
head(resTC_LB_vs_DA_LFC_annot)

# Extract significant results

resTC_LB_vs_DA_LFC_sig <- subset(resTC_LB_vs_DA_LFC, padj < 0.05)
resTC_LB_vs_DA_LFC_sig$transcript_id <- row.names(resTC_LB_vs_DA_LFC_sig)
resTC_LB_vs_DA_LFC_sig <- as.data.frame(resTC_LB_vs_DA_LFC_sig)
resTC_LB_vs_DA_LFC_sig <- resTC_LB_vs_DA_LFC_sig %>% filter(abs(log2FoldChange) >= 1.0)
nrow(resTC_LB_vs_DA_LFC_sig) #4059
nrow(resTC_LB_vs_DA_LFC) # 53379
     
# Annotate all significant genes and pull out apoptosis genes

resTC_LB_vs_DA_LFC_sig_annot <- left_join(resTC_LB_vs_DA_LFC_sig, C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product", "unique_go")], by ="transcript_id")
head(arrange(resTC_LB_vs_DA_LFC_sig_annot ,-log2FoldChange),n=50) # arrange in descending order

resTC_LB_vs_DA_LFC_sig_annot_apop <- resTC_LB_vs_DA_LFC_sig_annot[grepl(paste(Apoptosis_names,collapse="|"), 
                                                                        resTC_LB_vs_DA_LFC_sig_annot$product, ignore.case = TRUE),]
arrange(resTC_LB_vs_DA_LFC_sig_annot_apop , -log2FoldChange) 
nrow(resTC_LB_vs_DA_LFC_sig_annot_apop ) # 100

# Terms to remove

Terms_to_remove <- c("peptidyl-prolyl cis-trans isomerase G-like",
                     "peptidyl-prolyl cis-trans isomerase B-like",                                          
                     "endothelial zinc finger protein induced by tumor necrosis factor alpha",
                     "cAMP-dependent protein kinase type II regulatory subunit",
                     "scavenger receptor class F member")                    

resTC_LB_vs_DA_LFC_sig_annot_apop <- resTC_LB_vs_DA_LFC_sig_annot_apop[!grepl(paste(Terms_to_remove,collapse="|"), 
                                                                        resTC_LB_vs_DA_LFC_sig_annot_apop$product, ignore.case = TRUE),]

# Perform vsd calculation for gene clustering analysis 

Dermo_dds_timecourse_family_deseq_vsd <- vst(Dermo_dds_timecourse_family_deseq, blind=FALSE)

# Heatmap of top 40 most variable genes 
topVarGenes_Dermo_dds_timecourse_family_deseq_vsd <-  head(order(rowVars(assay(Dermo_dds_timecourse_family_deseq_vsd )), decreasing = TRUE), 40)
TC_res_mat <- assay(Dermo_dds_timecourse_family_deseq_vsd)[topVarGenes_Dermo_dds_timecourse_family_deseq_vsd,]
TC_res_mat <- TC_res_mat - rowMeans(TC_res_mat)
TC_res_anno <- as.data.frame(colData(Dermo_dds_timecourse_family_deseq_vsd)[, c("FamCode", "Timepoint")])
TC_heatmap <- pheatmap(TC_res_mat, annotation_col = TC_res_anno)
head(TC_res_mat)
# reorder annotation table to match ordering in heatmap 
TC_heatmap_reorder <-rownames(TC_res_mat[TC_heatmap$tree_row[["order"]],])
# annotate the row.names
TC_res_mat_prot <- as.data.frame(TC_heatmap_reorder)
colnames(TC_res_mat_prot)[1] <- "transcript_id"
TC_res_mat_prot_annot <- TC_res_mat_prot %>% left_join(select(C_vir_rtracklayer_transcripts_GO, transcript_id, product, gene), by = "transcript_id")

# Heatmap of apoptosis genes 
  # Subset vsd data frame based on transcript_id of apop genes because vsd rownames are the transcript_id names
apop_TC_res_mat <- assay(Dermo_dds_timecourse_family_deseq_vsd)[resTC_LB_vs_DA_LFC_sig_annot_apop$transcript_id,]
apop_TC_res_anno <- as.data.frame(colData(Dermo_dds_timecourse_family_deseq_vsd)[, c("FamCode", "Timepoint", "Treat")])
apop_TC_heatmap <- pheatmap(apop_TC_res_mat, annotation_col = apop_TC_res_anno)
# reorder annotation table to match ordering in heatmap 
apop_TC_heatmap_reorder <-rownames(apop_TC_res_mat[apop_TC_heatmap$tree_row[["order"]],])
# annotate the row.names
apop_TC_res_mat_prot <- as.data.frame(apop_TC_heatmap_reorder)
colnames(apop_TC_res_mat_prot)[1] <- "transcript_id"
apop_TC_res_mat_prot_annot <- apop_TC_res_mat_prot %>% left_join(select(C_vir_rtracklayer_transcripts_GO, transcript_id, product, gene), by = "transcript_id")
apop_TC_res_mat_prot_annot_list <- apop_TC_res_mat_prot_annot$product


# Log 2 fold change plot

TC_log2foldchange <- ggplot(resTC_LB_vs_DA_LFC_sig_annot_apop, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Time Course Analysis Regulated Cell Death Transcript \nLog2 Fold Change LB vs DA") +
  ylab("Log2 Fold Change between LB and DA")

## Apoptosis gene pathway analysis 

apop_TC_res_mat_prot_annot_list[grepl("TNF receptor", apop_TC_res_mat_prot_annot_list)] 
apop_TC_res_mat_prot_annot_list[grepl("toll", apop_TC_res_mat_prot_annot_list)] 
apop_TC_res_mat_prot_annot_list[grepl("caspase", apop_TC_res_mat_prot_annot_list)] 
apop_TC_res_mat_prot_annot_list[grepl("netrin", apop_TC_res_mat_prot_annot_list)] 
apop_TC_res_mat_prot_annot_list[grepl("GTPase", apop_TC_res_mat_prot_annot_list)] 
apop_TC_res_mat_prot_annot_list[grepl("baculoviral", apop_TC_res_mat_prot_annot_list)] 
apop_TC_res_mat_prot_annot_list[grepl("interferon", apop_TC_res_mat_prot_annot_list)] 
apop_TC_res_mat_prot_annot_list[grepl("macrophage migration", apop_TC_res_mat_prot_annot_list)] 
apop_TC_res_mat_prot_annot_list[grepl("cytochrome", apop_TC_res_mat_prot_annot_list)] 
apop_TC_res_mat_prot_annot_list[grepl("alpha-induced", apop_TC_res_mat_prot_annot_list)] 
apop_TC_res_mat_prot_annot_list[grepl("kinase", apop_TC_res_mat_prot_annot_list)] 
apop_TC_res_mat_prot_annot_list[grepl("cathepsin", apop_TC_res_mat_prot_annot_list)] 

# split heatmaps into separate based on pathways

TC_intrinsic <- c("IAP","mitogen-activated protein kinase","interferon-induced protein 44-like","interferon alpha-inducible protein 27-like",
                  "CD151","GTPase IMAP","nuclear apoptosis-inducing factor","cytochrome c1-2")
  
res_TC_intrinsic <- resTC_LB_vs_DA_LFC_sig_annot[grepl(paste(TC_intrinsic,collapse="|"), 
                                                       resTC_LB_vs_DA_LFC_sig_annot$product, ignore.case = TRUE),]

TC_extrinsic <- c("toll-like receptor","caspase-8-like")
res_TC_extrinsic <- resTC_LB_vs_DA_LFC_sig_annot[grepl(paste(TC_extrinsic,collapse="|"), 
                                                       resTC_LB_vs_DA_LFC_sig_annot$product, ignore.case = TRUE),]
TC_inflammation <- c("TNF receptor-associated factor","NF-kappa-B inhibitor","caspase-1-like")
res_TC_inflammation <- resTC_LB_vs_DA_LFC_sig_annot[grepl(paste(TC_inflammation,collapse="|"), 
                                                          resTC_LB_vs_DA_LFC_sig_annot$product, ignore.case = TRUE),]

TC_alternative <- c("caspase-1-like","macrophage migration inhibitory","cathepsin","complement","netrin","caspase-14-like")
res_TC_alternative <- resTC_LB_vs_DA_LFC_sig_annot[grepl(paste(TC_alternative,collapse="|"), 
                                                         resTC_LB_vs_DA_LFC_sig_annot$product, ignore.case = TRUE),]

# individual heatmaps

Tc_intrinsic_res_mat <- assay(Dermo_dds_timecourse_family_deseq_vsd)[res_TC_intrinsic$transcript_id,]
TC_extrinsic_res_mat <- assay(Dermo_dds_timecourse_family_deseq_vsd)[res_TC_extrinsic$transcript_id,]
TC_inflammation_res_mat <- assay(Dermo_dds_timecourse_family_deseq_vsd)[res_TC_inflammation$transcript_id,]
TC_alternative_res_mat <- assay(Dermo_dds_timecourse_family_deseq_vsd)[res_TC_alternative$transcript_id,]

TC_intrinsic_heatmap <- pheatmap(Tc_intrinsic_res_mat, annotation_col = apop_TC_res_anno)
TC_extrinsic_heatmap <- pheatmap(TC_extrinsic_res_mat, annotation_col = apop_TC_res_anno)
TC_inflammation_heatmap <- pheatmap(TC_inflammation_res_mat, annotation_col = apop_TC_res_anno)
TC_alternative_heatmap <- pheatmap(TC_alternative_res_mat, annotation_col = apop_TC_res_anno)

# reorder annotation table to match ordering in heatmap 
TC_intrinsic_reorder <- rownames(Tc_intrinsic_res_mat[TC_intrinsic_heatmap$tree_row[["order"]],])
TC_extrinsic_reorder <- rownames(TC_extrinsic_res_mat[TC_extrinsic_heatmap$tree_row[["order"]],])
TC_inflammation_reorder <- rownames(TC_inflammation_res_mat[TC_inflammation_heatmap$tree_row[["order"]],])
TC_alternative_reorder <- rownames(TC_alternative_res_mat[TC_alternative_heatmap$tree_row[["order"]],])

# annotate the row.names
TC_intrinsic_res_mat_prot <- as.data.frame(TC_intrinsic_reorder)
TC_extrinsic_res_mat_prot <- as.data.frame(TC_extrinsic_reorder )
TC_inflammation_res_mat_prot <- as.data.frame(TC_inflammation_reorder)
TC_alternative_res_mat_prot <- as.data.frame(TC_alternative_reorder)

colnames(TC_intrinsic_res_mat_prot)[1] <- "transcript_id"
colnames(TC_extrinsic_res_mat_prot)[1] <- "transcript_id"
colnames(TC_inflammation_res_mat_prot)[1] <- "transcript_id"
colnames(TC_alternative_res_mat_prot)[1] <- "transcript_id"

TC_intrinsic_res_mat_prot_annot <- TC_intrinsic_res_mat_prot %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id","product","gene")], by = "transcript_id")
TC_extrinsic_res_mat_prot_annot <- TC_extrinsic_res_mat_prot %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id","product","gene")], by = "transcript_id")
TC_inflammation_res_mat_prot_annot <- TC_inflammation_res_mat_prot %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id","product","gene")], by = "transcript_id")
TC_alternative_res_mat_prot_annot <- TC_alternative_res_mat_prot %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id","product","gene")], by = "transcript_id")

TC_intrinsic_res_mat_prot_annot_list <- TC_intrinsic_res_mat_prot_annot$product
TC_extrinsic_res_mat_prot_annot_list <- TC_extrinsic_res_mat_prot_annot$product
TC_inflammation_res_mat_prot_annot_list <- TC_inflammation_res_mat_prot_annot$product
TC_alternative_res_mat_prot_annot_list <- TC_alternative_res_mat_prot_annot$product

#### 7d INJECTED ANALYSIS LB vs. DA ######

# PREPARE AND SUBSET THE DATA
# Compare 7d injected to one another 

Dermo_coldata_7d_injected <- Dermo_coldata %>% filter(Treat =="Injected" & Timepoint == "7d")
rownames(Dermo_coldata_7d_injected) <- Dermo_coldata_7d_injected$rownames

# extract list of samples
SampleID_Dermo_coldata_7d_injected <- Dermo_coldata_7d_injected$rownames

# subset counts matrix
Dermo_counts_7d_injected <- Dermo_counts[,names(Dermo_counts) %in% SampleID_Dermo_coldata_7d_injected]

# Check all sample IDs in ColData are also in CountData and match their orders
all(rownames(Dermo_coldata_7d_injected) %in% colnames(Dermo_counts_7d_injected))  #Should return TRUE
# returns TRUE
all(rownames(Dermo_coldata_7d_injected) == colnames(Dermo_counts_7d_injected))   # should return TRUE
#returns TRUE

# Plot Dermo levels in the data set 
ggplot(Dermo_coldata_7d_injected, aes(x=as.factor(Ind),y=LogConc, color=FamCode)) + geom_point() + ylim(c(0,10))
# all are between 4 and 5 

# CHECK LEVELS
# 7 day just injected, 
levels(Dermo_coldata_7d_injected$FamCode) # "DA" "LB"
levels(Dermo_coldata_7d_injected$Treat) # "Control"  "Injected"
Dermo_coldata_7d_injected$Treat <- droplevels(Dermo_coldata_7d_injected$Treat) # drop control level
levels(Dermo_coldata_7d_injected$Treat)

# MAKE DESEQ DATA SET FROM MATRIX
# 7 day injected
Dermo_dds_7d_injected <- DESeqDataSetFromMatrix(countData= Dermo_counts_7d_injected,
                                                colData = Dermo_coldata_7d_injected,
                                                design = ~Library_Prep_Date + FamCode)
# PREFILTER THE DATA
Dermo_dds_7d_injected <- Dermo_dds_7d_injected[rowSums(counts(Dermo_dds_7d_injected)) >10,]

# VST TRANSFORMATION OF COUNTS FOR VISUALIZATION 

Dermo_dds_7d_injected_vsd <- vst(Dermo_dds_7d_injected, blind=FALSE)

# Calculate sample distances to assess overall sample similarity
# use function dist to calculate the Euclidean distance between samples. 
# To ensure we have a roughly equal contribution from all genes, we use it on the VST data. 
# We need to transpose the matrix of values using t, because the dist function expects the different samples to be rows of its argument, and different dimensions (here, genes) to be columns.

Dermo_dds_7d_injected_vsd_dist <- dist(t(assay(Dermo_dds_7d_injected_vsd)))
Dermo_dds_7d_injected_vsd_dist

# PCA plot visualization of individuals in the family #
plotPCA(Dermo_dds_7d_injected_vsd, intgroup=c("FamCode")) # still separation by family 
#DA and LB are separate but they don't seem to be very highly clustered together

# Calculate Differential expression
Dermo_dds_7d_injected_deseq <- DESeq(Dermo_dds_7d_injected) 

# Inspect Results names object
resultsNames(Dermo_dds_7d_injected_deseq) # [1] "Intercept"  "Library_Prep_Date_17_Dec_vs_15_Dec" "FamCode_LB_vs_DA"

# Build Results object
Dermo_dds_7d_injected_res <- results(Dermo_dds_7d_injected_deseq, alpha=0.05, name="FamCode_LB_vs_DA")
Dermo_dds_7d_injected_res # comparison is log2 fold change (MLE): FamCode LB vs DA

# LFC Shrinkage with Apeglm
Dermo_dds_7d_injected_res_LFC <- lfcShrink(Dermo_dds_7d_injected_deseq, coef="FamCode_LB_vs_DA", type="apeglm", res=Dermo_dds_7d_injected_res)

# Inspect results summary
summary(Dermo_dds_7d_injected_res_LFC)
#out of 44881 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 1188, 2.6%
#LFC < 0 (down)     : 1527, 3.4%
#outliers [1]       : 1691, 3.8%
#low counts [2]     : 4207, 9.4%
#(mean count < 2)

# Susbet results with padj < 0.05 and Logfold change of >1.0 or <-1.0
# only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
Dermo_dds_7d_injected_res_LFC_sig <- subset(Dermo_dds_7d_injected_res_LFC, padj < 0.05)
Dermo_dds_7d_injected_res_LFC_sig$transcript_id <- row.names(Dermo_dds_7d_injected_res_LFC_sig)
Dermo_dds_7d_injected_res_LFC_sig<- as.data.frame(Dermo_dds_7d_injected_res_LFC_sig)
Dermo_dds_7d_injected_res_LFC_sig <- Dermo_dds_7d_injected_res_LFC_sig %>% filter(abs(log2FoldChange) >= 1.0)
nrow(Dermo_dds_7d_injected_res_LFC_sig) #1836
nrow(Dermo_dds_family_res_LFC) # 53379

# Annotate transcripts
Dermo_dds_7d_injected_res_LFC_sig_annot <-Dermo_dds_7d_injected_res_LFC_sig %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product", "gene","unique_go")], by = "transcript_id")
head(arrange(Dermo_dds_7d_injected_res_LFC_sig_annot,-log2FoldChange),n=50) # arrange in descending order

# Isolate Apoptosis specific transcripts

Dermo_dds_7d_injected_res_LFC_sig_annot_APOP <- Dermo_dds_7d_injected_res_LFC_sig_annot[grepl(paste(Apoptosis_names,collapse="|"), 
                                                                                              Dermo_dds_7d_injected_res_LFC_sig_annot$product, ignore.case = TRUE),]
arrange(Dermo_dds_7d_injected_res_LFC_sig_annot_APOP , -log2FoldChange) 
nrow(Dermo_dds_7d_injected_res_LFC_sig_annot_APOP) # 65

# terms to remove 
Terms_to_remove <- c("peptidyl-prolyl cis-trans isomerase G-like",
                     "peptidyl-prolyl cis-trans isomerase B-like",
                     "peptidyl-prolyl cis-trans isomerase FKBP8",
                     "endothelial zinc finger protein induced by tumor necrosis factor alpha",
                     "cAMP-dependent protein kinase type II regulatory subunit",
                     "scavenger receptor class F member"
                     )       

Dermo_dds_7d_injected_res_LFC_sig_annot_APOP <- Dermo_dds_7d_injected_res_LFC_sig_annot_APOP[!grepl(paste(Terms_to_remove,collapse="|"), 
                                Dermo_dds_7d_injected_res_LFC_sig_annot_APOP$product, ignore.case = TRUE),]

# Heatmap of top 40 most variable genes 
topVarGenes_Dermo_dds_7d_injected_vsd <-  head(order(rowVars(assay(Dermo_dds_7d_injected_vsd )), decreasing = TRUE), 100)
Injected_7d_res_mat <- assay(Dermo_dds_7d_injected_vsd)[topVarGenes_Dermo_dds_7d_injected_vsd ,]
Injected_7d_res_mat <- Injected_7d_res_mat- rowMeans(Injected_7d_res_mat)
Injected_7d_res_anno <- as.data.frame(colData(Dermo_dds_7d_injected_vsd)[, c("FamCode","Timepoint")]) # needed to take a second column for some reason to get it to work 
Injected_7d_heatmap <- pheatmap(Injected_7d_res_mat, annotation_col = Injected_7d_res_anno)
head(Injected_7d_res_mat)
# reorder annotation table to match ordering in heatmap 
Injected_7d_reorder <-rownames(Injected_7d_res_mat[Injected_7d_heatmap$tree_row[["order"]],])
# annotate the row.names
Injected_7d_res_mat_prot <- as.data.frame(Injected_7d_reorder)
colnames(Injected_7d_res_mat_prot)[1] <- "transcript_id"
Injected_7d_res_mat_prot_annot <- Injected_7d_res_mat_prot %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product", "gene")], by = "transcript_id")
Injected_7d_res_mat_prot_annot$product

# Heatmap of apoptosis genes 
# Subset vsd data frame based on transcript_id of apop genes because vsd rownames are the transcript_id names
Injected_7d_Res_mat_apop <- assay(Dermo_dds_7d_injected_vsd)[Dermo_dds_7d_injected_res_LFC_sig_annot_APOP$transcript_id,]
Injected_7d_res_anno <- as.data.frame(colData(Dermo_dds_7d_injected_vsd)[, c("FamCode", "Timepoint")])
Injected_7d_heatmap <- pheatmap(Injected_7d_Res_mat_apop, annotation_col = Injected_7d_res_anno)
# reorder annotation table to match ordering in heatmap 
Injected_7d_heatmap_reorder <-rownames(Injected_7d_Res_mat_apop[Injected_7d_heatmap$tree_row[["order"]],])
# annotate the row.names
Injected_7d_Res_mat_apop_prot <- as.data.frame(Injected_7d_heatmap_reorder)
colnames(Injected_7d_Res_mat_apop_prot)[1] <- "transcript_id"
Injected_7d_Res_mat_apop_prot_annot <- Injected_7d_Res_mat_apop_prot %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product", "gene")], by = "transcript_id")
Injected_7d_Res_mat_apop_prot_annot_list <- Injected_7d_Res_mat_apop_prot_annot$product
Injected_7d_Res_mat_apop_prot_annot_list

# Split up genes by pathway 

Injected_7d_intrinsic <- c("IAP","interferon-induced protein 44-like","interferon alpha-inducible protein 27-like",
                  "CD151","GTPase IMAP","nuclear apoptosis-inducing factor","cytochrome c1-2", "E3 ubiquitin-protein ligase XIAP")

res_Injected_7d_intrinsic <- Dermo_dds_7d_injected_res_LFC_sig_annot_APOP[grepl(paste(Injected_7d_intrinsic,collapse="|"), 
                                   Dermo_dds_7d_injected_res_LFC_sig_annot_APOP$product, ignore.case = TRUE),]

Injected_7d_extrinsic <- c("toll-like receptor","caspase-8-like","TNF receptor-associated factor 2","ras-like GTP-binding protein RHO", "mitogen-activated protein kinase")
res_Injected_7d_extrinsic <- Dermo_dds_7d_injected_res_LFC_sig_annot_APOP[grepl(paste(Injected_7d_extrinsic,collapse="|"), 
                                                                            Dermo_dds_7d_injected_res_LFC_sig_annot_APOP$product, ignore.case = TRUE),]
Injected_7d_inflammation <- c("TNF receptor-associated factor 3","receptor-interacting serine/threonine-protein kinase 4", "TNF receptor-associated factor 6")
res_Injected_7d_inflammation <- Dermo_dds_7d_injected_res_LFC_sig_annot_APOP[grepl(paste(Injected_7d_inflammation,collapse="|"), 
                                                                      Dermo_dds_7d_injected_res_LFC_sig_annot_APOP$product, ignore.case = TRUE),]

Injected_7d_alternative <- c("caspase-1-like","macrophage migration inhibitory","cathepsin","complement","netrin","E3 ubiquitin-protein ligase CHIP")
res_Injected_7d_alternative <- Dermo_dds_7d_injected_res_LFC_sig_annot_APOP[grepl(paste(Injected_7d_alternative,collapse="|"), 
                                                                         Dermo_dds_7d_injected_res_LFC_sig_annot_APOP$product, ignore.case = TRUE),]
# To look for serine protease inhibitors
Day7_SPI <- Dermo_dds_7d_injected_res_LFC_sig_annot[grepl("serine protease inhibitor", Dermo_dds_7d_injected_res_LFC_sig_annot$product, ignore.case = TRUE),]
  # 

# Look for superoxide dismutase
SOD <- Dermo_dds_7d_injected_res_LFC_sig_annot[grepl("superoxide dismutase", Dermo_dds_7d_injected_res_LFC_sig_annot$product, ignore.case = TRUE),]
  # no sig genes
# individual heatmaps

Injected_7d_intrinsic_res_mat <- assay(Dermo_dds_7d_injected_vsd)[res_Injected_7d_intrinsic$transcript_id,]
Injected_7d_extrinsic_res_mat <- assay(Dermo_dds_7d_injected_vsd)[res_Injected_7d_extrinsic$transcript_id,]
Injected_7d_inflammation_res_mat <- assay(Dermo_dds_7d_injected_vsd)[res_Injected_7d_inflammation$transcript_id,]
Injected_7d_alternative_res_mat <- assay(Dermo_dds_7d_injected_vsd)[res_Injected_7d_alternative$transcript_id,]

Injected_7d_intrinsic_heatmap <- pheatmap(Injected_7d_intrinsic_res_mat, annotation_col = Injected_7d_res_anno)
Injected_7d_extrinsic_heatmap <- pheatmap(Injected_7d_extrinsic_res_mat, annotation_col = Injected_7d_res_anno)
Injected_7d_inflammation_heatmap <- pheatmap(Injected_7d_inflammation_res_mat, annotation_col = Injected_7d_res_anno)
Injected_7d_alternative_heatmap <- pheatmap(Injected_7d_alternative_res_mat, annotation_col = Injected_7d_res_anno)

# reorder annotation table to match ordering in heatmap 
Injected_7d_intrinsic_reorder <- rownames(Injected_7d_intrinsic_res_mat[Injected_7d_intrinsic_heatmap$tree_row[["order"]],])
Injected_7d_extrinsic_reorder <- rownames(Injected_7d_extrinsic_res_mat[Injected_7d_extrinsic_heatmap$tree_row[["order"]],])
Injected_7d_inflammation_reorder <- rownames(Injected_7d_inflammation_res_mat[Injected_7d_inflammation_heatmap$tree_row[["order"]],])
Injected_7d_alternative_reorder <- rownames(Injected_7d_alternative_res_mat[Injected_7d_alternative_heatmap$tree_row[["order"]],])

# annotate the row.names
Injected_7d_intrinsic_res_mat_prot <- as.data.frame(Injected_7d_intrinsic_reorder)
Injected_7d_extrinsic_res_mat_prot <- as.data.frame(Injected_7d_extrinsic_reorder )
Injected_7d_inflammation_res_mat_prot <- as.data.frame(Injected_7d_inflammation_reorder)
Injected_7d_alternative_res_mat_prot <- as.data.frame(Injected_7d_alternative_reorder)

colnames(Injected_7d_intrinsic_res_mat_prot)[1] <- "transcript_id"
colnames(Injected_7d_extrinsic_res_mat_prot)[1] <- "transcript_id"
colnames(Injected_7d_inflammation_res_mat_prot)[1] <- "transcript_id"
colnames(Injected_7d_alternative_res_mat_prot)[1] <- "transcript_id"

Injected_7d_intrinsic_res_mat_prot_annot <-    Injected_7d_intrinsic_res_mat_prot %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id","product","gene")], by = "transcript_id")
Injected_7d_extrinsic_res_mat_prot_annot <-    Injected_7d_extrinsic_res_mat_prot %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id","product","gene")], by = "transcript_id")
Injected_7d_inflammation_res_mat_prot_annot <- Injected_7d_inflammation_res_mat_prot %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id","product","gene")], by = "transcript_id")
Injected_7d_alternative_res_mat_prot_annot <-  Injected_7d_alternative_res_mat_prot %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id","product","gene")], by = "transcript_id")

Injected_7d_intrinsic_res_mat_prot_annot_list <-    Injected_7d_intrinsic_res_mat_prot_annot$product
Injected_7d_extrinsic_res_mat_prot_annot_list <-    Injected_7d_extrinsic_res_mat_prot_annot$product
Injected_7d_inflammation_res_mat_prot_annot_list <- Injected_7d_inflammation_res_mat_prot_annot$product
Injected_7d_alternative_res_mat_prot_annot_list <-  Injected_7d_alternative_res_mat_prot_annot$product

# Plots of apop log2fold change

Dermo_dds_7d_injected_res_LFC_sig_annot_APOP

head(Dermo_dds_7d_injected_res_LFC_sig_annot_APOP)

# Expand ggplot color plot 

Day7_injected_log2foldchange <- ggplot(Dermo_dds_7d_injected_res_LFC_sig_annot_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Day 7 Injected Regulated Cell Death Transcript \nLog2 Fold Change LB vs DA") +
  ylab("Log2 Fold Change between LB and DA")

#### 6h INJECTED ANALYSIS LB vs. DA ######

# PREPARE AND SUBSET THE DATA
# Compare 6h injected to one another 

Dermo_coldata_6hr_injected <- Dermo_coldata %>% filter(Treat =="Injected" & Timepoint == "6h")
row.names(Dermo_coldata_6hr_injected) <- Dermo_coldata_6hr_injected$rownames

# extract list of samples
SampleID_Dermo_coldata_6hr_injected <- row.names(Dermo_coldata_6hr_injected)

# subset counts matrix
Dermo_counts_6hr_injected <- Dermo_counts[,names(Dermo_counts) %in% SampleID_Dermo_coldata_6hr_injected]

# Check all sample IDs in ColData are also in CountData and match their orders
all(rownames(Dermo_coldata_6hr_injected) %in% colnames(Dermo_counts_6hr_injected))  #Should return TRUE
# returns TRUE
all(rownames(Dermo_coldata_6hr_injected) == colnames(Dermo_counts_6hr_injected))   # should return TRUE
#returns TRUE

# CHECK LEVELS
# 6hr injected, 
levels(Dermo_coldata_6hr_injected$FamCode) # "DA" "LB"
levels(Dermo_coldata_6hr_injected$Treat) # "Control"  "Injected"
Dermo_coldata_6hr_injected$Treat <- droplevels(Dermo_coldata_6hr_injected$Treat) # drop control level
levels(Dermo_coldata_6hr_injected$Treat)
levels(Dermo_coldata_6hr_injected$Timepoint)  # "6h"  "36h" "7d" 
Dermo_coldata_6hr_injected$Timepoint <- droplevels(Dermo_coldata_6hr_injected$Timepoint) # drop control level
levels(Dermo_coldata_6hr_injected$Timepoint) 

# MAKE DESEQ DATA SET FROM MATRIX
# 6hr injected
Dermo_dds_6hr_injected <- DESeqDataSetFromMatrix(countData= Dermo_counts_6hr_injected,
                                                colData = Dermo_coldata_6hr_injected,
                                                design =  ~ FamCode)
# PREFILTER THE DATA
Dermo_dds_6hr_injected <- Dermo_dds_6hr_injected[rowSums(counts(Dermo_dds_6hr_injected)) >10,]

# VST TRANSFORMATION OF COUNTS FOR VISUALIZATION 

Dermo_dds_6hr_injected_vsd <- vst(Dermo_dds_6hr_injected, blind=FALSE)

# Calculate sample distances to assess overall sample similarity
# use function dist to calculate the Euclidean distance between samples. 
# To ensure we have a roughly equal contribution from all genes, we use it on the VST data. 
# We need to transpose the matrix of values using t, because the dist function expects the different samples to be rows of its argument, and different dimensions (here, genes) to be columns.

Dermo_dds_6hr_injected_vsd_dist <- dist(t(assay(Dermo_dds_6hr_injected_vsd)))
Dermo_dds_6hr_injected_vsd_dist

# PCA plot visualization of individuals in the family #
plotPCA(Dermo_dds_6hr_injected_vsd, intgroup=c("FamCode")) # still separation by family 
#DA and LB are separate but they don't seem to be very highly clustered together

# Calculate Differential expression
Dermo_dds_6hr_injected_deseq <- DESeq(Dermo_dds_6hr_injected) 

# Inspect Results names object
resultsNames(Dermo_dds_6hr_injected_deseq) # [1] "Intercept" "FamCode_LB_vs_DA"

# Build Results object
Dermo_dds_6hr_injected_res <- results(Dermo_dds_6hr_injected_deseq, alpha=0.05, name="FamCode_LB_vs_DA")
Dermo_dds_6hr_injected_res # comparison is log2 fold change (MLE): FamCode LB vs DA

# LFC Shrinkage with Apeglm
Dermo_dds_6hr_injected_res_LFC <- lfcShrink(Dermo_dds_6hr_injected_deseq, coef="FamCode_LB_vs_DA", type="apeglm", res=Dermo_dds_6hr_injected_res)

# Inspect results summary
summary(Dermo_dds_6hr_injected_res_LFC)
#out of 41971 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 392, 0.93%
#LFC < 0 (down)     : 388, 0.92%
#outliers [1]       : 737, 1.8%
#low counts [2]     : 6510, 16%
#(mean count < 5)

# Susbet results with padj < 0.05 and Logfold change of >1.0 or <-1.0
# only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
Dermo_dds_6hr_injected_res_LFC_sig <- subset(Dermo_dds_6hr_injected_res_LFC, padj < 0.05)
Dermo_dds_6hr_injected_res_LFC_sig$transcript_id <- row.names(Dermo_dds_6hr_injected_res_LFC_sig)
Dermo_dds_6hr_injected_res_LFC_sig<- as.data.frame(Dermo_dds_6hr_injected_res_LFC_sig)
Dermo_dds_6hr_injected_res_LFC_sig <- Dermo_dds_6hr_injected_res_LFC_sig %>% filter(abs(log2FoldChange) >= 1.0)
nrow(Dermo_dds_6hr_injected_res_LFC_sig) #449
nrow(Dermo_dds_6hr_injected_res_LFC) # 41971

# Annotate transcripts
Dermo_dds_6hr_injected_res_LFC_sig_annot <-Dermo_dds_6hr_injected_res_LFC_sig %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product", "gene","unique_go")], by = "transcript_id")
head(arrange(Dermo_dds_6hr_injected_res_LFC_sig_annot,-log2FoldChange),n=50) # arrange in descending order

# Isolate Apoptosis specific transcripts

Dermo_dds_6hr_injected_res_LFC_sig_annot_APOP <- Dermo_dds_6hr_injected_res_LFC_sig_annot[grepl(paste(Apoptosis_names,collapse="|"), 
                                                                                              Dermo_dds_6hr_injected_res_LFC_sig_annot$product, ignore.case = TRUE),]
arrange(Dermo_dds_6hr_injected_res_LFC_sig_annot_APOP , -log2FoldChange) 
nrow(Dermo_dds_6hr_injected_res_LFC_sig_annot_APOP) # 15

# terms to remove 
Terms_to_remove <- c("peptidyl-prolyl cis-trans isomerase G-like",
                     "peptidyl-prolyl cis-trans isomerase B-like",
                     "peptidyl-prolyl cis-trans isomerase FKBP8",
                     "endothelial zinc finger protein induced by tumor necrosis factor alpha",
                     "cAMP-dependent protein kinase type II regulatory subunit",
                     "scavenger receptor class F member"
)       

Dermo_dds_6hr_injected_res_LFC_sig_annot_APOP <- Dermo_dds_6hr_injected_res_LFC_sig_annot_APOP[!grepl(paste(Terms_to_remove,collapse="|"), 
                                                                                                    Dermo_dds_6hr_injected_res_LFC_sig_annot_APOP$product, ignore.case = TRUE),]

# Heatmap of top 40 most variable genes 
topVarGenes_Dermo_dds_6hr_injected_vsd <-  head(order(rowVars(assay(Dermo_dds_6hr_injected_vsd )), decreasing = TRUE), 100)
Injected_6hr_res_mat <- assay(Dermo_dds_6hr_injected_vsd)[topVarGenes_Dermo_dds_6hr_injected_vsd ,]
Injected_6hr_res_mat <- Injected_6hr_res_mat- rowMeans(Injected_6hr_res_mat)
Injected_6hr_res_anno <- as.data.frame(colData(Dermo_dds_6hr_injected_vsd)[, c("FamCode","Timepoint")]) # needed to take a second column for some reason to get it to work 
Injected_6hr_heatmap <- pheatmap(Injected_6hr_res_mat, annotation_col = Injected_6hr_res_anno)
head(Injected_6hr_res_mat)
# reorder annotation table to match ordering in heatmap 
Injected_6hr_reorder <-rownames(Injected_6hr_res_mat[Injected_6hr_heatmap$tree_row[["order"]],])
# annotate the row.names
Injected_6hr_res_mat_prot <- as.data.frame(Injected_6hr_reorder)
colnames(Injected_6hr_res_mat_prot)[1] <- "transcript_id"
Injected_6hr_res_mat_prot_annot <- Injected_6hr_res_mat_prot %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product", "gene")], by = "transcript_id")
Injected_6hr_res_mat_prot_annot$product

# Heatmap of apoptosis genes 
# Subset vsd data frame based on transcript_id of apop genes because vsd rownames are the transcript_id names
Injected_6hr_Res_mat_apop <- assay(Dermo_dds_6hr_injected_vsd)[Dermo_dds_6hr_injected_res_LFC_sig_annot_APOP$transcript_id,]
Injected_6hr_res_anno <- as.data.frame(colData(Dermo_dds_6hr_injected_vsd)[, c("FamCode", "Timepoint")])
Injected_6hr_heatmap <- pheatmap(Injected_6hr_Res_mat_apop, annotation_col = Injected_6hr_res_anno)
# reorder annotation table to match ordering in heatmap 
Injected_6hr_heatmap_reorder <-rownames(Injected_6hr_Res_mat_apop[Injected_6hr_heatmap$tree_row[["order"]],])
# annotate the row.names
Injected_6hr_Res_mat_apop_prot <- as.data.frame(Injected_6hr_heatmap_reorder)
colnames(Injected_6hr_Res_mat_apop_prot)[1] <- "transcript_id"
Injected_6hr_Res_mat_apop_prot_annot <- Injected_6hr_Res_mat_apop_prot %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product", "gene")], by = "transcript_id")
Injected_6hr_Res_mat_apop_prot_annot_list <- Injected_6hr_Res_mat_apop_prot_annot$product
Injected_6hr_Res_mat_apop_prot_annot_list

# Plots of apop log2fold change

Sixhr_injected_log2foldchange <- ggplot(Dermo_dds_6hr_injected_res_LFC_sig_annot_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("6 hr Injected Regulated Cell Death Transcript \nLog2 Fold Change LB vs DA") +
  ylab("Log2 Fold Change between LB and DA")


##### DAY 7 DA Control vs. Injected #######

# Comparing the Day 7 control to day 7 injected will give a better idea of the response of susceptible oysters alone,
# this may help me better plan my assays, as now the "resitant" family is out

# PREPARE AND SUBSET THE DATA
# Compare 7d injected to one another 

Dermo_coldata_7d_suscept <- Dermo_coldata %>% filter(FamCode =="DA" & Timepoint == "7d")
rownames(Dermo_coldata_7d_suscept) <- Dermo_coldata_7d_suscept$rownames

# extract list of samples
SampleID_Dermo_coldata_7d_suscept <- Dermo_coldata_7d_suscept$rownames

# subset counts matrix
Dermo_counts_7d_suscept <- Dermo_counts[,names(Dermo_counts) %in% SampleID_Dermo_coldata_7d_suscept]

# Check all sample IDs in ColData are also in CountData and match their orders
all(rownames(Dermo_coldata_7d_suscept) %in% colnames(Dermo_counts_7d_suscept))  #Should return TRUE
# returns TRUE
all(rownames(Dermo_coldata_7d_suscept) == colnames(Dermo_counts_7d_suscept))   # should return TRUE
#returns TRUE

# Plot Dermo levels in the data set 
ggplot(Dermo_coldata_7d_suscept, aes(x=as.factor(Ind),y=LogConc, color=Treat)) + geom_point() + ylim(c(0,10))
# all are between 4 and 5 

# CHECK LEVELS
# 7 day just susceptible, 
levels(Dermo_coldata_7d_suscept$FamCode) # "DA" "LB"
levels(Dermo_coldata_7d_suscept$Treat) # "Control"  "Injected"
Dermo_coldata_7d_suscept$FamCode <- droplevels(Dermo_coldata_7d_suscept$FamCode) # drop control level
levels(Dermo_coldata_7d_suscept$FamCode) # DA

# MAKE DESEQ DATA SET FROM MATRIX
# 7 day injected
Dermo_dds_7d_suscept <- DESeqDataSetFromMatrix(countData= Dermo_counts_7d_suscept,
                                                colData = Dermo_coldata_7d_suscept,
                                                design = ~Library_Prep_Date + Treat)
# PREFILTER THE DATA
Dermo_dds_7d_suscept <- Dermo_dds_7d_suscept[rowSums(counts(Dermo_dds_7d_suscept)) >10,]

# VST TRANSFORMATION OF COUNTS FOR VISUALIZATION 

Dermo_dds_7d_suscept_vsd <- vst(Dermo_dds_7d_suscept, blind=FALSE)

# Calculate sample distances to assess overall sample similarity
# use function dist to calculate the Euclidean distance between samples. 
# To ensure we have a roughly equal contribution from all genes, we use it on the VST data. 
# We need to transpose the matrix of values using t, because the dist function expects the different samples to be rows of its argument, and different dimensions (here, genes) to be columns.

Dermo_dds_7d_suscept_vsd_dist <- dist(t(assay(Dermo_dds_7d_suscept_vsd)))
Dermo_dds_7d_suscept_vsd_dist

# PCA plot visualization of individuals in DA
plotPCA(Dermo_dds_7d_suscept_vsd, intgroup=c("FamCode", "Treat")) # some clustering by treatment, explains about 40% of the  variance
#DA and LB are separate but they don't seem to be very highly clustered together

# Calculate Differential expression
Dermo_dds_7d_suscept_deseq <- DESeq(Dermo_dds_7d_suscept) 

# Inspect Results names object
resultsNames(Dermo_dds_7d_suscept_deseq) # [1] "Intercept"  "Library_Prep_Date_17_Dec_vs_15_Dec" "Treat_Injected_vs_Control" 

# Build Results object
Dermo_dds_7d_suscept_res <- results(Dermo_dds_7d_suscept_deseq, alpha=0.05, name="Treat_Injected_vs_Control")
Dermo_dds_7d_suscept_res # comparison is log2 fold change (MLE): Treat Injected vs Control

# LFC Shrinkage with Apeglm
Dermo_dds_7d_suscept_res_LFC <- lfcShrink(Dermo_dds_7d_suscept_deseq, coef="Treat_Injected_vs_Control", type="apeglm", res=Dermo_dds_7d_suscept_res)

# Inspect results summary
#summary(Dermo_dds_7d_suscept_res_LFC)
#out of 44740 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 59, 0.13%
#LFC < 0 (down)     : 110, 0.25%
#outliers [1]       : 1614, 3.6%
#low counts [2]     : 13736, 31%
#(mean count < 10)

# Susbet results with padj < 0.05 and Logfold change of >1.0 or <-1.0
# only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
Dermo_dds_7d_suscept_res_LFC_sig <- subset(Dermo_dds_7d_suscept_res_LFC, padj < 0.05)
Dermo_dds_7d_suscept_res_LFC_sig$transcript_id <- row.names(Dermo_dds_7d_suscept_res_LFC_sig)
Dermo_dds_7d_suscept_res_LFC_sig<- as.data.frame(Dermo_dds_7d_suscept_res_LFC_sig)
Dermo_dds_7d_suscept_res_LFC_sig <- Dermo_dds_7d_suscept_res_LFC_sig%>% filter(abs(log2FoldChange) >= 1.0)
nrow(Dermo_dds_7d_suscept_res_LFC_sig) #61
nrow(Dermo_dds_7d_suscept_res_LFC) # 44740

# Annotate transcripts
Dermo_dds_7d_suscept_res_LFC_sig_annot <-Dermo_dds_7d_suscept_res_LFC_sig %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product", "gene","unique_go")], by = "transcript_id")
head(arrange(Dermo_dds_7d_suscept_res_LFC_sig_annot ,-log2FoldChange),n=50) # arrange in descending order

# Isolate Apoptosis specific transcripts

Dermo_dds_7d_suscept_res_LFC_sig_annot_APOP <- Dermo_dds_7d_suscept_res_LFC_sig_annot[grepl(paste(Apoptosis_names,collapse="|"), 
                                                                                            Dermo_dds_7d_suscept_res_LFC_sig_annot$product, ignore.case = TRUE),]
arrange(Dermo_dds_7d_suscept_res_LFC_sig_annot_APOP , -log2FoldChange) 
nrow(Dermo_dds_7d_injected_res_LFC_sig_annot_APOP) # 1 

# terms to remove 
Terms_to_remove <- c("peptidyl-prolyl cis-trans isomerase G-like",
                     "peptidyl-prolyl cis-trans isomerase B-like",
                     "peptidyl-prolyl cis-trans isomerase FKBP8",
                     "endothelial zinc finger protein induced by tumor necrosis factor alpha",
                     "cAMP-dependent protein kinase type II regulatory subunit",
                     "scavenger receptor class F member"
)       

Dermo_dds_7d_suscept_res_LFC_sig_annot_APOP <- Dermo_dds_7d_suscept_res_LFC_sig_annot_APOP[!grepl(paste(Terms_to_remove,collapse="|"), 
                                                                                                  Dermo_dds_7d_suscept_res_LFC_sig_annot_APOP$product, ignore.case = TRUE),]
## no significant apoptosis genes 

### LB control vs. challenge all timepoints, DA control vs. challenge all timepoints ####

Dermo_coldata_LB_all <- Dermo_coldata %>% filter(FamCode == "LB")
row.names(Dermo_coldata_LB_all) <- Dermo_coldata_LB_all$rownames
Dermo_coldata_DA_all <- Dermo_coldata %>% filter(FamCode == "DA")
row.names(Dermo_coldata_DA_all) <- Dermo_coldata_DA_all$rownames

# extract list of samples
SampleID_LB <- Dermo_coldata_LB_all$rownames
SampleID_DA <- Dermo_coldata_DA_all$rownames

# subset counts matrix
Dermo_counts_LB_all <- Dermo_counts[,names(Dermo_counts) %in% SampleID_LB]
Dermo_counts_DA_all <- Dermo_counts[,names(Dermo_counts) %in% SampleID_DA]

# Check all sample IDs in ColData are also in CountData and match their orders
all(rownames(Dermo_coldata_LB_all) %in% colnames(Dermo_counts_LB_all))  #Should return TRUE
all(rownames(Dermo_coldata_LB_all) == colnames(Dermo_counts_LB_all))
# returns TRUE
all(rownames(Dermo_coldata_DA_all) %in% colnames(Dermo_counts_DA_all))    # should return TRUE
all(rownames(Dermo_coldata_DA_all) == colnames(Dermo_counts_DA_all))
#returns TRUE

# Check factor levels, set it so that comparison group is the first

# LB all timepoints
levels(Dermo_coldata_LB_all$Timepoint) #"6h"  "36h" "7d"
levels(Dermo_coldata_LB_all$FamCode) # "DA" "LB"
Dermo_coldata_LB_all$FamCode <- droplevels(Dermo_coldata_LB_all$FamCode)
levels(Dermo_coldata_LB_all$FamCode) # LB

# DA all timepoints
levels(Dermo_coldata_DA_all$Timepoint) # #"6h"  "36h" "7d" 
levels(Dermo_coldata_DA_all$FamCode)
Dermo_coldata_DA_all$FamCode <- droplevels(Dermo_coldata_DA_all$FamCode)
levels(Dermo_coldata_DA_all$FamCode) # DA

# Make DESeq data set from matrix
# The dds objects below looks at the effect of treatment within family
Dermo_dds_LB_treat <- DESeqDataSetFromMatrix(countData = Dermo_counts_LB_all,
                                            colData = Dermo_coldata_LB_all,
                                            design = ~Library_Prep_Date + Treat)

Dermo_dds_DA_treat <- DESeqDataSetFromMatrix(countData = Dermo_counts_DA_all,
                                            colData = Dermo_coldata_DA_all,
                                            design = ~Library_Prep_Date  + Treat)
# PREFILTER THE DATA
Dermo_dds_LB_treat <- Dermo_dds_LB_treat[rowSums(counts(Dermo_dds_LB_treat)) >10,]
Dermo_dds_DA_treat <- Dermo_dds_DA_treat[rowSums(counts(Dermo_dds_DA_treat)) >10,]

# VST TRANSFORMATION OF COUNTS FOR VISUALIZATION 
Dermo_dds_LB_treat_vsd <- vst(Dermo_dds_LB_treat, blind=FALSE)
Dermo_dds_DA_treat_vsd <- vst(Dermo_dds_DA_treat, blind=FALSE)

# Calculate sample distances to assess overall sample similarity
# use function dist to calculate the Euclidean distance between samples. 
# To ensure we have a roughly equal contribution from all genes, we use it on the VST data. 
# We need to transpose the matrix of values using t, because the dist function expects the different samples to be rows of its argument, and different dimensions (here, genes) to be columns.

Dermo_dds_LB_treat_vsd_dist <- dist(t(assay(Dermo_dds_LB_treat_vsd)))
Dermo_dds_DA_treat_vsd_dist <- dist(t(assay(Dermo_dds_DA_treat_vsd)))

# PCA plot visualization of individuals in DA
plotPCA(Dermo_dds_LB_treat_vsd, intgroup=c("FamCode", "Treat")) # some clustering by treatment, explains about 40% of the  variance
plotPCA(Dermo_dds_DA_treat_vsd, intgroup=c("FamCode", "Treat"))

# Calculate Differential expression
Dermo_dds_LB_treat_deseq <- DESeq(Dermo_dds_LB_treat) 
Dermo_dds_DA_treat_deseq <- DESeq(Dermo_dds_DA_treat)

# Inspect Results names object
resultsNames(Dermo_dds_LB_treat_deseq) # [1] "Intercept"  "Library_Prep_Date_17_Dec_vs_15_Dec" "Treat_Injected_vs_Control" 
resultsNames(Dermo_dds_DA_treat_deseq) # [1] "Intercept"  "Library_Prep_Date_17_Dec_vs_15_Dec" "Treat_Injected_vs_Control" 

# Build Results object
Dermo_dds_LB_treat_deseq_res <- results(Dermo_dds_LB_treat_deseq, alpha=0.05, name="Treat_Injected_vs_Control")
Dermo_dds_LB_treat_deseq_res # comparison is log2 fold change (MLE): Treat Injected vs Control

Dermo_dds_DA_treat_deseq_res <- results(Dermo_dds_DA_treat_deseq, alpha=0.05, name="Treat_Injected_vs_Control")
Dermo_dds_DA_treat_deseq_res # comparison is log2 fold change (MLE): Treat Injected vs Control

# LFC Shrinkage with Apeglm
Dermo_dds_LB_treat_deseq_res_LFC <- lfcShrink(Dermo_dds_LB_treat_deseq, coef="Treat_Injected_vs_Control", type="apeglm", res=Dermo_dds_LB_treat_deseq_res)
Dermo_dds_DA_treat_deseq_res_LFC <- lfcShrink(Dermo_dds_DA_treat_deseq, coef="Treat_Injected_vs_Control", type="apeglm", res=Dermo_dds_DA_treat_deseq_res)

# Inspect results summary
summary(Dermo_dds_LB_treat_deseq_res_LFC)
#out of 49668 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 154, 0.31%
#LFC < 0 (down)     : 150, 0.3%
#outliers [1]       : 531, 1.1%
#low counts [2]     : 12335, 25%
#(mean count < 4)

summary(Dermo_dds_DA_treat_deseq_res_LFC)
#out of 49533 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 8, 0.016%
#LFC < 0 (down)     : 25, 0.05%
#outliers [1]       : 500, 1%
#low counts [2]     : 9584, 19%
#(mean count < 3)

# Susbet results with padj < 0.05 and Logfold change of >1.0 or <-1.0
# only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
Dermo_dds_LB_treat_deseq_res_LFC_sig <- subset(Dermo_dds_LB_treat_deseq_res_LFC, padj < 0.05)
Dermo_dds_LB_treat_deseq_res_LFC_sig$transcript_id <- row.names(Dermo_dds_LB_treat_deseq_res_LFC_sig)
Dermo_dds_LB_treat_deseq_res_LFC_sig<- as.data.frame(Dermo_dds_LB_treat_deseq_res_LFC_sig)
Dermo_dds_LB_treat_deseq_res_LFC_sig <- Dermo_dds_LB_treat_deseq_res_LFC_sig%>% filter(abs(log2FoldChange) >= 1.0)
nrow(Dermo_dds_LB_treat_deseq_res_LFC_sig) #115
nrow(Dermo_dds_LB_treat_deseq_res_LFC) # 49870

Dermo_dds_DA_treat_deseq_res_LFC_sig <- subset(Dermo_dds_DA_treat_deseq_res_LFC, padj < 0.05)
Dermo_dds_DA_treat_deseq_res_LFC_sig$transcript_id <- row.names(Dermo_dds_DA_treat_deseq_res_LFC_sig)
Dermo_dds_DA_treat_deseq_res_LFC_sig<- as.data.frame(Dermo_dds_DA_treat_deseq_res_LFC_sig)
Dermo_dds_DA_treat_deseq_res_LFC_sig <- Dermo_dds_DA_treat_deseq_res_LFC_sig%>% filter(abs(log2FoldChange) >= 1.0)
nrow(Dermo_dds_DA_treat_deseq_res_LFC_sig) #11
nrow(Dermo_dds_DA_treat_deseq_res_LFC) # 49827

# Annotate transcripts
Dermo_dds_LB_treat_deseq_res_LFC_sig_annot <- Dermo_dds_LB_treat_deseq_res_LFC_sig %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product", "gene","unique_go")], by = "transcript_id")
head(arrange(Dermo_dds_LB_treat_deseq_res_LFC_sig_annot ,-log2FoldChange),n=50) # arrange in descending order by absolute value 
Dermo_dds_LB_treat_deseq_res_LFC_sig_annot[sort(abs(Dermo_dds_LB_treat_deseq_res_LFC_sig_annot$log2FoldChange),decreasing=T,index.return=T)[[2]],]

Dermo_dds_DA_treat_deseq_res_LFC_sig_annot <-Dermo_dds_DA_treat_deseq_res_LFC_sig %>% left_join(C_vir_rtracklayer_transcripts_GO[,c("transcript_id", "product", "gene","unique_go")], by = "transcript_id")
Dermo_dds_DA_treat_deseq_res_LFC_sig_annot[sort(abs(Dermo_dds_DA_treat_deseq_res_LFC_sig_annot$log2FoldChange),decreasing=T,index.return=T)[[2]],] # arrange in descending order by absolute value

# Isolate Apoptosis specific transcripts
Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_APOP <- Dermo_dds_LB_treat_deseq_res_LFC_sig_annot[grepl(paste(Apoptosis_names,collapse="|"), 
                                                                                                    Dermo_dds_LB_treat_deseq_res_LFC_sig_annot$product, ignore.case = TRUE),]
arrange(Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_APOP , -log2FoldChange) 
nrow(Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_APOP)  #6

Dermo_dds_DA_treat_deseq_res_LFC_sig_annot_APOP <- Dermo_dds_DA_treat_deseq_res_LFC_sig_annot[grepl(paste(Apoptosis_names,collapse="|"), 
                                                                                                    Dermo_dds_DA_treat_deseq_res_LFC_sig_annot$product, ignore.case = TRUE),]
arrange(Dermo_dds_DA_treat_deseq_res_LFC_sig_annot_APOP , -log2FoldChange) #0 

# terms to remove 
Terms_to_remove <- c("peptidyl-prolyl cis-trans isomerase G-like",
                     "peptidyl-prolyl cis-trans isomerase B-like",
                     "peptidyl-prolyl cis-trans isomerase FKBP8",
                     "endothelial zinc finger protein induced by tumor necrosis factor alpha",
                     "cAMP-dependent protein kinase type II regulatory subunit",
                     "scavenger receptor class F member"
)       

Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_APOP <- Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_APOP[!grepl(paste(Terms_to_remove,collapse="|"), 
                                                                                                          Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_APOP$product, ignore.case = TRUE),]
Dermo_dds_DA_treat_deseq_res_LFC_sig_annot_APOP <- Dermo_dds_DA_treat_deseq_res_LFC_sig_annot_APOP[!grepl(paste(Terms_to_remove,collapse="|"), 
                                                                                                          Dermo_dds_DA_treat_deseq_res_LFC_sig_annot_APOP$product, ignore.case = TRUE),]
# Plot LogFC for significant genes
# create new column with product and XM id combined to make unique for factor levels using unite
Dermo_dds_LB_treat_deseq_res_LFC_sig_annot <- unite(Dermo_dds_LB_treat_deseq_res_LFC_sig_annot, product_unite, c(7,6), sep="_", remove=FALSE)

# create new column with product and XM id combined to make unique for factor levels using unite
Dermo_dds_DA_treat_deseq_res_LFC_sig_annot <- unite(Dermo_dds_DA_treat_deseq_res_LFC_sig_annot, product_unite, c(7,6), sep="_", remove=FALSE)

# arrange in descending order with absolute value
Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_sorted <- Dermo_dds_LB_treat_deseq_res_LFC_sig_annot[sort(abs(Dermo_dds_LB_treat_deseq_res_LFC_sig_annot$log2FoldChange),decreasing=T,index.return=T)[[2]],]
Dermo_dds_DA_treat_deseq_res_LFC_sig_annot_sorted <- Dermo_dds_DA_treat_deseq_res_LFC_sig_annot[sort(abs(Dermo_dds_DA_treat_deseq_res_LFC_sig_annot$log2FoldChange),decreasing=T,index.return=T)[[2]],]

# set factor levels explicitly
Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_sorted$product_unite <- factor(Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_sorted$product_unite, levels = Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_sorted$product_unite[order(Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_sorted$log2FoldChange)])
# not all would fit on one graph 
Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_sorted_up <- Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_sorted %>% filter(log2FoldChange >= 0.0) 
Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_sorted_up$direction <- "up"
Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_sorted_down <- Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_sorted %>% filter(log2FoldChange < 0.0) 
Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_sorted_down$direction <- "down"
Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_sorted_combined <- rbind(Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_sorted_up, Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_sorted_down)

Dermo_dds_DA_treat_deseq_res_LFC_sig_annot_sorted$product_unite <- factor(Dermo_dds_DA_treat_deseq_res_LFC_sig_annot_sorted$product_unite, levels = Dermo_dds_DA_treat_deseq_res_LFC_sig_annot_sorted$product_unite[order(Dermo_dds_DA_treat_deseq_res_LFC_sig_annot_sorted$log2FoldChange)])

LB_treat_log2foldchange <- ggplot(Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_sorted_combined, aes(x=product_unite, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("LB control vs. challenge LogFC") +
  ylab("Log2 Fold Change between control and challenge") + theme(axis.text.x = element_text(size=1) ) + facet_grid(direction~., scales="free")

LB_treat_log2foldchange_up <- ggplot(Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_sorted_up, aes(x=product_unite, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("LB control vs. challenge LogFC") +
  ylab("Log2 Fold Change between control and challenge") + theme(axis.text.x = element_text(size=1) ) 

LB_treat_log2foldchange_down <- ggplot(Dermo_dds_LB_treat_deseq_res_LFC_sig_annot_sorted_down, aes(x=product_unite, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("LB control vs. challenge LogFC") +
  ylab("Log2 Fold Change between control and challenge") + theme(axis.text.x = element_text(size=1) ) 

DA_treat_log2foldchange <- ggplot(Dermo_dds_DA_treat_deseq_res_LFC_sig_annot_sorted, aes(x=product_unite, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("DA control vs. challenge LogFC") +
  ylab("Log2 Fold Change between control and challenge") + theme_bw()

##### Test of significant enrichment of apopsois GO terms in the datasets to confirm #####
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("topGO")
# install.packages("genefilter")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  + install.packages("BiocManager")
# BiocManager::install("Rgraphviz")
library(topGO)
library(genefilter)
library(Rgraphviz)
library(GO.db)
library(qdapTools)


# GO ENRICHMENT IN THE TIMECOURSE ANALYSIS DATA SET 

# Prepare two input data tables, input is a table that contains the genes, the GO mappings, and the p-values
  # one table will be the full list of GO annotated genes
  # second table will be the significant genes that I have decided above 

# Use timecourse all significant annotated genes and subset the columns for GO, padj, and XM_name
resTC_LB_vs_DA_LFC_sig_annot_GO_analysis <- resTC_LB_vs_DA_LFC_sig_annot[,c("transcript_id","unique_go","padj")]
resTC_LB_vs_DA_LFC_annot_GO_analysis <- resTC_LB_vs_DA_LFC_annot[,c("transcript_id","unique_go","padj")]
#only 62 lines with significant GO terms, not quite enough to make any interesting conclusions

# Extract only the lines that have a GO annotation
resTC_LB_vs_DA_LFC_sig_annot_GO_analysis <- resTC_LB_vs_DA_LFC_sig_annot_GO_analysis %>% filter(unique_go != "NULL")
nrow(resTC_LB_vs_DA_LFC_sig_annot_GO_analysis ) #62
resTC_LB_vs_DA_LFC_annot_GO_analysis <- resTC_LB_vs_DA_LFC_annot_GO_analysis %>% filter(unique_go != "NULL")
nrow(resTC_LB_vs_DA_LFC_annot_GO_analysis) #1691

# use resTC_LB_vs_DA_LFC_sig_annot_GO_analysis GO terms column and p value column to load into REVIGO
resTC_LB_vs_DA_LFC_sig_annot_GO_analysis[,c("unique_go","padj")]
write.csv()

#Create topGO objects
#This object will contain all information necessary for the GO analysis, 
# namely the list of genes, the list of interesting genes, the gene
#scores (if available) and the part of the GO ontology (the GO graph) 
#which needs to be used in the analysis. The topGOdata object will be the
#input of the testing procedures, the evaluation and visualisation functions.
#Can input custom annotations already found, need a gene to GO mapping approach
#parse the annotations file with function readMappings()
#file format required by readMappings = gene_ID<TAB>GO_ID1, GO_ID2, GO_ID3,

#format GO mapping files to readMapping format
TC_GO_readmapping_sig <- resTC_LB_vs_DA_LFC_sig_annot_GO_analysis[,c("transcript_id", "unique_go")]
class(TC_GO_readmapping_sig$transcript_id)
class(TC_GO_readmapping_sig$unique_go)
class(TC_GO_readmapping_sig)

TC_GO_readmapping <- resTC_LB_vs_DA_LFC_annot_GO_analysis[,c("transcript_id", "unique_go")]
class(TC_GO_readmapping$transcript_id)
class(TC_GO_readmapping$unique_go)
class(TC_GO_readmapping)

# First coerce the data.frame to all-character
TC_GO_readmapping = TC_GO_readmapping %>% mutate(unique_go = sapply(unique_go, toString))
class(TC_GO_readmapping$transcript_id)
class(TC_GO_readmapping$unique_go)
class(TC_GO_readmapping)
write.table(TC_GO_readmapping, file="TC_GO_readMapping.txt", row.names=FALSE, sep= "\t")

TC_GO_readmapping_sig = TC_GO_readmapping_sig %>% mutate(unique_go = sapply(unique_go, toString))
class(TC_GO_readmapping_sig$transcript_id)
class(TC_GO_readmapping_sig$unique_go)
class(TC_GO_readmapping_sig)
write.table(TC_GO_readmapping_sig, file="TC_GO_readMapping_sig.txt", row.names=FALSE, sep= "\t")

#must be tab delimited
#remove first line in bash:
#tail -n +2 TC_GO_readMapping.txt > TC_GO_readMapping_header_gone.txt
#tail -n +2 TC_GO_readMapping_sig.txt > TC_GO_readMapping_sig_header_gone.txt 
#sed 's/\"//g' TC_GO_readMapping_header_gone.txt > TC_GO_readMapping_header_gone_nocommas.txt
#sed 's/\"//g' TC_GO_readMapping_sig_header_gone.txt > TC_GO_readMapping_sig_header_gone_nocommas.txt

#Use readMappings() to great gene to GO mapping, you have to import it as a file. It won't let you not import it as an outside file
# Use the full list and not the sig list
TC_GO_geneID2GO <- readMappings(file="TC_GO_readMapping_header_gone_nocommas.txt")
str(head(TC_GO_geneID2GO))

# define gene universe and gene list (gene list)
TC_GO_Gene_Universe <- names(TC_GO_geneID2GO) # gene universe is the full file
TC_GO_genes_of_interest <- TC_GO_readmapping_sig$transcript_id  
TC_GO_genelist <- factor(as.integer(TC_GO_Gene_Universe %in% TC_GO_genes_of_interest))
names(TC_GO_genelist) <- TC_GO_Gene_Universe
nrow(TC_GO_genelist)

TC_GOdata <- new("topGOdata", description = "Timecourse Dermo Gene Enrichment", ontology = "BP",
                   allGenes = TC_GO_genelist,  annot = annFUN.gene2GO, gene2GO = TC_GO_geneID2GO)
      # Node size trims the GO hierarchy, we don't need this 
TC_GOdata

## Enrichment Analysis
# Run Fisher's exact test using Weight01 algorithm which takes into account the GO hierarchy, "classic" does not
resultWeight <- runTest(TC_GOdata, algorithm="weight01", statistic="fisher")
resultFisher <- runTest(TC_GOdata, algorithm = "classic", statistic = "fisher")
score(resultFisher)
# Print out resultsing significant node 
allRes <- GenTable(TC_GOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher")
allRes

#GO.ID                                        Term Annotated Significant Expected classicFisher
#1  GO:0010647 positive regulation of cell communicatio...         1           1     0.03         0.026
#2  GO:0023056            positive regulation of signaling         1           1     0.03         0.026
#3  GO:0046928    regulation of neurotransmitter secretion         1           1     0.03         0.026
#4  GO:0048167           regulation of synaptic plasticity         1           1     0.03         0.026
#5  GO:0050806 positive regulation of synaptic transmis...         1           1     0.03         0.026
#6  GO:0051588    regulation of neurotransmitter transport         1           1     0.03         0.026
#7  GO:0060291             long-term synaptic potentiation         1           1     0.03         0.026
#8  GO:0060341         regulation of cellular localization         1           1     0.03         0.026
#9  GO:0016567                      protein ubiquitination        11           2     0.28         0.030
#10 GO:0032446 protein modification by small protein co...        11           2     0.28         0.030

# REVIGO input will be the result of the score(resultFisher) named vector after converting to data frame
  # p-values here are a value of the enrichment

resultsFisher_score <- score(resultFisher)
TC_REVIGO_input <- as.data.frame(resultsFisher_score)
head(TC_REVIGO_input)
TC_REVIGO_input$GO_term <- row.names(TC_REVIGO_input)
head(TC_REVIGO_input)
TC_REVIGO_input <- TC_REVIGO_input[,c(2,1)]

# filter out everything with p value not below 0.1
TC_REVIGO_input_subset <- TC_REVIGO_input %>% filter(resultsFisher_score <= 0.1)
#only 27 GO terms fit this criteria 

write.table(file="./Dermo_2015_Analysis/OUTPUT/TC_REVIGO_input_subset.txt", TC_REVIGO_input_subset)
# in terminal remove quotes are the GO term: sed 's/\"//g' TC_REVIGO_input_subset.txt > TC_REVIGO_input_subset_noquotes.txt
# cut -d' ' -f2,3 TC_REVIGO_input_subset_noquotes.txt
# select in REVIGO that the numbers associated with the GO term are p-values 

# REVIGO plotting script
library( ggplot2 )
library( scales )

# Here is your data from REVIGO. Scroll down for plot configuration options.
revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0008152","metabolic process",75.387, 2.930, 5.930, 6.986,-1.1957,0.980,0.000),
                     c("GO:0010647","positive regulation of cell communication", 0.359, 5.010,-3.587, 4.663,-1.5904,0.555,0.000),
                     c("GO:0051259","protein oligomerization", 0.188,-0.246, 5.372, 4.382,-1.1239,0.843,0.027),
                     c("GO:0016567","protein ubiquitination", 0.523,-4.761,-4.523, 4.827,-1.5246,0.726,0.030),
                     c("GO:0043170","macromolecule metabolic process",39.491,-5.313, 0.579, 6.705,-1.3178,0.879,0.073),
                     c("GO:0006807","nitrogen compound metabolic process",38.744,-4.571, 3.506, 6.696,-1.1008,0.920,0.074),
                     c("GO:0044238","primary metabolic process",53.743,-2.464, 3.058, 6.839,-1.1115,0.929,0.090),
                     c("GO:0006310","DNA recombination", 1.641,-2.302,-7.167, 5.323,-1.3258,0.764,0.150),
                     c("GO:0019538","protein metabolic process",18.489,-4.191,-2.479, 6.375,-1.1537,0.814,0.194),
                     c("GO:0060341","regulation of cellular localization", 0.194, 5.627, 0.013, 4.396,-1.5904,0.688,0.223),
                     c("GO:0023051","regulation of signaling", 0.934, 6.431,-2.774, 5.079,-1.2947,0.681,0.257),
                     c("GO:0055085","transmembrane transport", 8.916, 5.459, 3.623, 6.058,-1.0856,0.866,0.318),
                     c("GO:0006508","proteolysis", 5.223,-4.446,-3.656, 5.826,-1.0456,0.781,0.335),
                     c("GO:0006259","DNA metabolic process", 5.607,-2.912,-6.126, 5.857,-1.5069,0.770,0.346),
                     c("GO:0070647","protein modification by small protein conjugation or removal", 0.821,-4.258,-4.759, 5.023,-1.2188,0.752,0.470),
                     c("GO:0010646","regulation of cell communication", 0.929, 4.502,-2.760, 5.076,-1.2947,0.668,0.473),
                     c("GO:0015074","DNA integration", 0.682,-1.862,-7.628, 4.942,-1.1661,0.775,0.598),
                     c("GO:0051588","regulation of neurotransmitter transport", 0.014, 5.643, 0.660, 3.240,-1.5904,0.704,0.645),
                     c("GO:0048167","regulation of synaptic plasticity", 0.027, 5.547,-4.148, 3.543,-1.5904,0.515,0.670),
                     c("GO:0051260","protein homooligomerization", 0.106,-0.553, 5.354, 4.132,-1.1239,0.844,0.701),
                     c("GO:0023056","positive regulation of signaling", 0.356, 5.892,-3.537, 4.660,-1.5904,0.569,0.706),
                     c("GO:0050804","modulation of synaptic transmission", 0.057, 5.478,-3.772, 3.866,-1.2947,0.525,0.800),
                     c("GO:0060291","long-term synaptic potentiation", 0.010, 5.317,-4.500, 3.124,-1.5904,0.533,0.858),
                     c("GO:0046928","regulation of neurotransmitter secretion", 0.012, 5.537,-1.951, 3.184,-1.5904,0.456,0.865),
                     c("GO:0050806","positive regulation of synaptic transmission", 0.020, 5.190,-4.213, 3.415,-1.5904,0.521,0.893));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
head(one.data)

# Names of the axes, sizes of the numbers and letters, names of the columns, etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ]; 
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);
p1 <- p1 + ggtitle("Significantly Encriched terms over time LB vs DA")

# Output the plot to screen

p1

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("C:/Users/path_to_your_file/revigo-plot.pdf");

# CONCLUSION: GO enrichment is not extremely informative because so few genes have GO terms

#### Import Perkinsus marinus genome ######

# Import Perkinsus gff file from NCBI 
PMG_annotation <- import("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/C_VIRGINICA_PIPELINE/GCF_000006405.1_JCVI_PMG_1.0_genomic.gff")
PMG_annotation <- as.data.frame(PMG_annotation)

# Isolate the lines with Transcripts
PMG_annotation_transcripts <- PMG_annotation %>% filter(transcript_id != "NA")
head(PMG_annotation_transcripts)

# Isolate unique transcript lines
PMG_annotation_transcripts_unique <- PMG_annotation_transcripts[!duplicated(PMG_annotation_transcripts$transcript_id),]
head(PMG_annotation_transcripts_unique)

##### Perkinsus Differential Expression Day 7 Injected LB vs. DA ######

# Load counts Matrix file
Perkinsus_counts_data <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/RAW DATA AND INFO/2015 Dermo Challenge /2015_transcriptomes_CLC_PMG_mapping/Perkinsus_mapping_TE_Counts_Matrix.csv", header=TRUE)
row.names(Perkinsus_counts_data) <- Perkinsus_counts_data$X # loading in from csv made the rownames a column
Perkinsus_counts_data <- Perkinsus_counts_data[,-1]
head(Perkinsus_counts_data)

# Load Edited ColData with "PMG" tag added to samples 
PMG_coldata <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/RAW DATA AND INFO/2015 Dermo Challenge /6h_36h_7d_DA_LB_PMG_metadata.csv", header=TRUE) 
head(PMG_coldata)
row.names(PMG_coldata) <- PMG_coldata$X

# PREPARE AND SUBSET THE DATA
PMG_coldata_7d_injected <- PMG_coldata %>% filter(Treat =="Injected_PMG" & Timepoint == "7d")
row.names(PMG_coldata_7d_injected) <- PMG_coldata_7d_injected$X
head(PMG_coldata_7d_injected)

# extract list of samples
SampleID_PMG_coldata_7d_injected <- row.names(PMG_coldata_7d_injected)

# subset counts matrix
Perkinsus_counts_7d_injected <- Perkinsus_counts_data[,names(Perkinsus_counts_data) %in% SampleID_PMG_coldata_7d_injected]

# Check all sample IDs in ColData are also in CountData and match their orders
all(rownames(PMG_coldata_7d_injected) %in% colnames(Perkinsus_counts_7d_injected))  #Should return TRUE
# returns TRUE
all(rownames(PMG_coldata_7d_injected) == colnames(Perkinsus_counts_7d_injected))   # should return TRUE
#returns TRUE! 

# CHECK LEVELS
# 7 day just injected, 
levels(PMG_coldata_7d_injected$FamCode) # "DA" "LB"
levels(PMG_coldata_7d_injected$Treat) # "Control"  "Injected"
PMG_coldata_7d_injected$Treat <- droplevels(PMG_coldata_7d_injected$Treat) # drop control level
levels(PMG_coldata_7d_injected$Treat) # Injected

# MAKE DESEQ DATA SET FROM MATRIX
# Perkinsus 7 day injected
Perkinsus_dds_7d_injected <- DESeqDataSetFromMatrix(countData= Perkinsus_counts_7d_injected,
                                                colData = PMG_coldata_7d_injected,
                                                design = ~Library_Prep_Date + FamCode)
# PREFILTER THE DATA
Perkinsus_dds_7d_injected <- Perkinsus_dds_7d_injected[rowSums(counts(Perkinsus_dds_7d_injected)) >10,] 
nrow(Perkinsus_dds_7d_injected) # changing to less than 10 doesn't make a difference

# VST TRANSFORMATION OF COUNTS FOR VISUALIZATION 
# Need to use the varianceStabilizingTransformation directly because nsub rows because the default is 1000 and there are only 738
Perkinsus_dds_7d_injected_vsd <- varianceStabilizingTransformation(Perkinsus_dds_7d_injected)

# Calculate sample distances to assess overall sample similarity
# use function dist to calculate the Euclidean distance between samples. 
# To ensure we have a roughly equal contribution from all genes, we use it on the VST data. 
# We need to transpose the matrix of values using t, because the dist function expects the different samples to be rows of its argument, and different dimensions (here, genes) to be columns.

Perkinsus_dds_7d_injected_vsd_dist <- dist(t(assay(Perkinsus_dds_7d_injected_vsd)))
Perkinsus_dds_7d_injected_vsd_dist

# PCA plot visualization of individuals in the family #
plotPCA(Perkinsus_dds_7d_injected_vsd, intgroup=c("FamCode")) # There is separation by family! 
#DA and LB are separate but they don't seem to be very highly clustered together

# Calculate Differential expression
Perkinsus_dds_7d_injected_deseq <- DESeq(Perkinsus_dds_7d_injected ) 

# Inspect Results names object
resultsNames(Perkinsus_dds_7d_injected_deseq) # [1] "Intercept"  "Library_Prep_Date_17_Dec_vs_15_Dec" "FamCode_LB_vs_DA"

# Build Results object
Perkinsus_dds_7d_injected_res <- results(Perkinsus_dds_7d_injected_deseq, alpha=0.05, name="FamCode_LB_vs_DA")
Perkinsus_dds_7d_injected_res # comparison is log2 fold change (MLE): FamCode LB vs DA

# LFC Shrinkage with Apeglm
Perkinsus_dds_7d_injected_res_LFC <- lfcShrink(Perkinsus_dds_7d_injected_deseq, coef="FamCode_LB_vs_DA", type="apeglm", res=Perkinsus_dds_7d_injected_res)

# Inspect results summary
summary(Perkinsus_dds_7d_injected_res_LFC)
#out of 738 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 9, 1.2%
#LFC < 0 (down)     : 13, 1.8%
#outliers [1]       : 6, 0.81%
#low counts [2]     : 401, 54%
#(mean count < 3)

# Susbet results with padj < 0.05 and Logfold change of >1.0 or <-1.0
# only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
Perkinsus_dds_7d_injected_res_LFC_sig <- subset(Perkinsus_dds_7d_injected_res_LFC, padj < 0.05)
Perkinsus_dds_7d_injected_res_LFC_sig$transcript_id <- row.names(Perkinsus_dds_7d_injected_res_LFC_sig)
Perkinsus_dds_7d_injected_res_LFC_sig<- as.data.frame(Perkinsus_dds_7d_injected_res_LFC_sig)
Perkinsus_dds_7d_injected_res_LFC_sig <- Perkinsus_dds_7d_injected_res_LFC_sig %>% filter(abs(log2FoldChange) >= 1.0)
nrow(Perkinsus_dds_7d_injected_res_LFC_sig) #16
nrow(Perkinsus_dds_7d_injected_res_LFC) # 738

# Annotate transcripts
Perkinsus_dds_7d_injected_res_LFC_sig_annot <- Perkinsus_dds_7d_injected_res_LFC_sig %>% left_join(PMG_annotation_transcripts_unique[,c("transcript_id", "product")], by = "transcript_id")
head(arrange(Perkinsus_dds_7d_injected_res_LFC_sig_annot,-log2FoldChange),n=50) # arrange in descending order

# Log2 Fold Change 
Day7_PMG_injected_log2foldchange <- ggplot(Perkinsus_dds_7d_injected_res_LFC_sig_annot, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Day 7 Injected Significant Perkinsus expression \nLog2 Fold Change LB vs DA") +
  ylab("Log2 Fold Change between LB and DA")



##### Perkinsus Differential Expression Analysis All Timepoints Injected LB vs. DA ######

# PREPARE AND SUBSET THE DATA
PMG_coldata_All_injected <- PMG_coldata %>% filter(Treat =="Injected_PMG")
row.names(PMG_coldata_All_injected) <- PMG_coldata_All_injected$X
head(PMG_coldata_All_injected)

# extract list of samples
SampleID_PMG_coldata_All_injected <- row.names(PMG_coldata_All_injected)

# subset counts matrix
Perkinsus_counts_All_injected <- Perkinsus_counts_data[,names(Perkinsus_counts_data) %in% SampleID_PMG_coldata_All_injected]

# Check all sample IDs in ColData are also in CountData and match their orders
all(rownames(PMG_coldata_All_injected) %in% colnames(Perkinsus_counts_All_injected))  #Should return TRUE
# returns TRUE
all(rownames(PMG_coldata_All_injected) == colnames(Perkinsus_counts_All_injected))   # should return TRUE
#returns TRUE! 

# CHECK LEVELS
# All injected 
levels(PMG_coldata_All_injected$FamCode) # "DA" "LB"
levels(PMG_coldata_All_injected$Treat) # "Control"  "Injected"
PMG_coldata_All_injected$Treat <- droplevels(PMG_coldata_All_injected$Treat) # drop control level
levels(PMG_coldata_All_injected$Treat) # Injected

# MAKE DESEQ DATA SET FROM MATRIX
# Perkinsus all injected
Perkinsus_dds_All_injected <- DESeqDataSetFromMatrix(countData= Perkinsus_counts_All_injected,
                                                    colData = PMG_coldata_All_injected,
                                                    design = ~Library_Prep_Date + FamCode)
# PREFILTER THE DATA
Perkinsus_dds_All_injected <- Perkinsus_dds_All_injected[rowSums(counts(Perkinsus_dds_All_injected)) >10,] 
nrow(Perkinsus_dds_All_injected) # 5169

# VST TRANSFORMATION OF COUNTS FOR VISUALIZATION 
# Need to use the varianceStabilizingTransformation directly because nsub rows because the default is 1000 and there are only 738
Perkinsus_dds_All_injected_vsd <- varianceStabilizingTransformation(Perkinsus_dds_All_injected)

# Calculate sample distances to assess overall sample similarity
# use function dist to calculate the Euclidean distance between samples. 
# To ensure we have a roughly equal contribution from all genes, we use it on the VST data. 
# We need to transpose the matrix of values using t, because the dist function expects the different samples to be rows of its argument, and different dimensions (here, genes) to be columns.

Perkinsus_dds_All_injected_vsd_dist <- dist(t(assay(Perkinsus_dds_All_injected_vsd)))
Perkinsus_dds_All_injected_vsd_dist

# PCA plot visualization of individuals in the family #
plotPCA(Perkinsus_dds_All_injected_vsd, intgroup=c("FamCode")) # PC1 explains almost all variance, may be just correlated with load

# Calculate Differential expression
Perkinsus_dds_All_injected_deseq <- DESeq(Perkinsus_dds_All_injected ) 

# Inspect Results names object
resultsNames(Perkinsus_dds_All_injected_deseq) # [1] "Intercept"  "Library_Prep_Date_17_Dec_vs_15_Dec" "FamCode_LB_vs_DA"

# Build Results object
Perkinsus_dds_All_injected_res <- results(Perkinsus_dds_All_injected_deseq, alpha=0.05, name="FamCode_LB_vs_DA")
Perkinsus_dds_All_injected_res # comparison is log2 fold change (MLE): FamCode LB vs DA

# LFC Shrinkage with Apeglm
Perkinsus_dds_All_injected_res_LFC <- lfcShrink(Perkinsus_dds_All_injected_deseq, coef="FamCode_LB_vs_DA", type="apeglm", res=Perkinsus_dds_All_injected_res)

# Inspect results summary
summary(Perkinsus_dds_All_injected_res_LFC)
#out of 5165 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 15, 0.29%
#LFC < 0 (down)     : 16, 0.31%
#outliers [1]       : 3, 0.058%
#low counts [2]     : 3709, 72%
#(mean count < 1)

# Susbet results with padj < 0.05 and Logfold change of >1.0 or <-1.0
# only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
Perkinsus_dds_All_injected_res_LFC_sig <- subset(Perkinsus_dds_All_injected_res_LFC, padj < 0.05)
Perkinsus_dds_All_injected_res_LFC_sig_not_less_1 <- Perkinsus_dds_All_injected_res_LFC_sig
Perkinsus_dds_All_injected_res_LFC_sig_not_less_1$transcript_id <- row.names(Perkinsus_dds_All_injected_res_LFC_sig_not_less_1)
Perkinsus_dds_All_injected_res_LFC_sig_not_less_1<- as.data.frame(Perkinsus_dds_All_injected_res_LFC_sig_not_less_1)
nrow(Perkinsus_dds_All_injected_res_LFC_sig_not_less_1) # 31

Perkinsus_dds_All_injected_res_LFC_sig$transcript_id <- row.names(Perkinsus_dds_All_injected_res_LFC_sig)
Perkinsus_dds_All_injected_res_LFC_sig<- as.data.frame(Perkinsus_dds_All_injected_res_LFC_sig)
Perkinsus_dds_All_injected_res_LFC_sig <- Perkinsus_dds_All_injected_res_LFC_sig %>% filter(abs(log2FoldChange) >= 1.0)
nrow(Perkinsus_dds_All_injected_res_LFC_sig) #12
nrow(Perkinsus_dds_All_injected_res_LFC) # 5169

# Annotate transcripts
Perkinsus_dds_All_injected_res_LFC_sig_annot <- Perkinsus_dds_All_injected_res_LFC_sig %>% left_join(PMG_annotation_transcripts_unique[,c("transcript_id", "product")], by = "transcript_id")
head(arrange(Perkinsus_dds_All_injected_res_LFC_sig_annot,-log2FoldChange),n=50) # arrange in descending order

Perkinsus_dds_All_injected_res_LFC_sig_not_less_1_annot <- Perkinsus_dds_All_injected_res_LFC_sig_not_less_1 %>% left_join(PMG_annotation_transcripts_unique[,c("transcript_id", "product")], by = "transcript_id")
head(arrange(Perkinsus_dds_All_injected_res_LFC_sig_not_less_1_annot,-log2FoldChange),n=50) # arrange in descending order

# Log2 Fold Change 
All_timepoints_PMG_injected_log2foldchange <- ggplot(Perkinsus_dds_All_injected_res_LFC_sig_annot, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("All Injected Significant Perkinsus expression \nLog2 Fold Change LB vs DA") +
  ylab("Log2 Fold Change between LB and DA")
  # similar results as the other one
  # issue may just be low read counts for Perkinsus 

# Compare this to the Day 7 injected DA vs LB
Perkinsus_dds_All_injected_res_LFC_sig_annot[!(Perkinsus_dds_All_injected_res_LFC_sig_annot$product %in% Perkinsus_dds_7d_injected_res_LFC_sig_annot$product),]
  # transcripts not in this all timepoint but in Day 7
#baseMean log2FoldChange     lfcSE       pvalue         padj  transcript_id                                        product
#3  143.665815       1.518006 0.6761758 6.065041e-04 3.240492e-02 XM_002774352.1                   tubulin beta chain, putative
#6    5.639496       1.940772 0.6298091 6.470120e-05 7.108106e-03 XM_002779366.1                             helicase, putative
#8   46.831470      -1.285116 0.3761378 1.862502e-05 3.876664e-03 XM_002776817.1                 conserved hypothetical protein
#9   16.515663       6.262025 1.1496807 2.521513e-11 1.836923e-08 XM_002779492.1                 proteosome component, putative
#10   4.841730       2.116334 0.6783346 6.146912e-05 7.108106e-03 XM_002772938.1  40S ribosomal protein S4, X isoform, putative
#11  40.249652      -2.066824 0.7324291 4.701465e-05 6.850035e-03 XM_002773276.1 thioredoxin reductase 1, cytoplasmic, putative
Perkinsus_dds_7d_injected_res_LFC_sig_annot[!(Perkinsus_dds_7d_injected_res_LFC_sig_annot$product %in% Perkinsus_dds_All_injected_res_LFC_sig_annot$product),]
  # transcript in Day 7 injected analysis not in the all timepoints comparison 
  #product
  #1             eukaryotic translation initiation factor, putative
  #2               ubiquitin-conjugating enzyme E2-17 kDa, putative
  #7                                         beta-tubulin, putative
  #8                      heat shock protein 70 precursor, putative
  #10                            isocitrate dehydrogenase, putative
  #11  translation elongation factor EF-1, subunit alpha,, putative
  #12 ribonucleoside-diphosphate reductase, large subunit, putative
  #15                               ribosomal protein S20, putative
  #16                delta-aminolevulinic acid synthetase, putative


##### Perkinsus Differential Expression control vs. injected (no samples removed) ####

# Check all sample IDs in ColData are also in CountData and match their orders
all(rownames(PMG_coldata) %in% colnames(Perkinsus_counts_data))  #Should return TRUE
# returns TRUE
all(rownames(PMG_coldata) == colnames(Perkinsus_counts_data))   # should return TRUE
#returns TRUE! 

# CHECK LEVELS
# All injected 
levels(PMG_coldata$FamCode) # "DA" "LB"
levels(PMG_coldata$Treat) # "Control"  "Injected"

# MAKE DESEQ DATA SET FROM MATRIX
# Perkinsus all injected
Perkinsus_dds <- DESeqDataSetFromMatrix(countData= Perkinsus_counts_data,
                                        colData = PMG_coldata,
                                        design = ~Library_Prep_Date + Treat + FamCode)
# PREFILTER THE DATA
Perkinsus_dds <- Perkinsus_dds[rowSums(counts(Perkinsus_dds)) >10,] 
nrow(Perkinsus_dds) # 5236

# VST TRANSFORMATION OF COUNTS FOR VISUALIZATION 
# Need to use the varianceStabilizingTransformation directly because nsub rows because the default is 1000 and there are only 738
Perkinsus_dds_vsd <- varianceStabilizingTransformation(Perkinsus_dds)

# Calculate sample distances to assess overall sample similarity
# use function dist to calculate the Euclidean distance between samples. 
# To ensure we have a roughly equal contribution from all genes, we use it on the VST data. 
# We need to transpose the matrix of values using t, because the dist function expects the different samples to be rows of its argument, and different dimensions (here, genes) to be columns.

Perkinsus_dds_vsd_dist <- dist(t(assay(Perkinsus_dds_vsd)))
Perkinsus_dds_vsd_dist

# PCA plot visualization of individuals in the family #
plotPCA(Perkinsus_dds_vsd, intgroup=c("LogConc")) # PC1 explains almost all variance, may be just correlated with load
plotPCA(Perkinsus_dds_vsd, intgroup=c("Treat","FamCode","Timepoint"))

# Is there a correlation between LogConc and Number of reads? (this is an associated question) 

# Calculate Differential expression
Perkinsus_dds_deseq <- DESeq(Perkinsus_dds) 

# Inspect Results names object
resultsNames(Perkinsus_dds_deseq) # [1] "Intercept" "Library_Prep_Date_17_Dec_vs_15_Dec" "Treat_Injected_PMG_vs_Control_PMG" "FamCode_LB_vs_DA"  

# Build Results object
Perkinsus_dds_res <- results(Perkinsus_dds_deseq, alpha=0.05, name="Treat_Injected_PMG_vs_Control_PMG")
Perkinsus_dds_res # comparison is log2 fold change (MLE): FamCode LB vs DA

# LFC Shrinkage with Apeglm
Perkinsus_dds_res_LFC <- lfcShrink(Perkinsus_dds_deseq, coef="Treat_Injected_PMG_vs_Control_PMG", type="apeglm", res=Perkinsus_dds_res)

# Inspect results summary
summary(Perkinsus_dds_res_LFC)
#out of 5257 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 830, 16%
#LFC < 0 (down)     : 6, 0.11%
#outliers [1]       : 3, 0.057%
#low counts [2]     : 3366, 64%
#(mean count < 1)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

# Susbet results with padj < 0.05 and Logfold change of >1.0 or <-1.0
# only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
Perkinsus_dds_res_LFC_sig <- subset(Perkinsus_dds_res_LFC, padj < 0.05)
Perkinsus_dds_res_LFC_sig$transcript_id <- row.names(Perkinsus_dds_res_LFC_sig)
Perkinsus_dds_res_LFC_sig<- as.data.frame(Perkinsus_dds_res_LFC_sig)
Perkinsus_dds_res_LFC_sig <- Perkinsus_dds_res_LFC_sig %>% filter(abs(log2FoldChange) >= 1.0)
nrow(Perkinsus_dds_res_LFC_sig) #8
nrow(Perkinsus_dds_res_LFC) # 5260

# Annotate transcripts
Perkinsus_dds_res_LFC_sig_annot <- Perkinsus_dds_res_LFC_sig %>% left_join(PMG_annotation_transcripts_unique[,c("transcript_id", "product")], by = "transcript_id")
head(arrange(Perkinsus_dds_res_LFC_sig_annot,-log2FoldChange),n=50) # arrange in descending order

# Log2 Fold Change 
All_PMG_log2foldchange <- ggplot(Perkinsus_dds_res_LFC_sig_annot, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("All Injected Significant Perkinsus expression \nLog2 Fold Change LB vs DA") +
  ylab("Log2 Fold Change between Control and Injected PMG")


##### Perkinsus Differential Expression control vs. injected (Zero load removed) ####
# Goal of this comparison is to see if there are differences in expression between low load and high,
  # which mostly corresponds to the family and treatment groups 

# Rank samples by LOAD
PMG_coldata_infected <- PMG_coldata %>% filter(LogConc > 0.00)
PMG_coldata_infected[order(PMG_coldata_infected$LogConc),] 
row.names(PMG_coldata_infected) <- PMG_coldata_infected$X

# extract list of samples
SampleID_PMG_coldata_infected <- row.names(PMG_coldata_infected)

# Extract those samples from the counts matrix 
Perkinsus_counts_infected <- Perkinsus_counts_data[,names(Perkinsus_counts_data) %in% SampleID_PMG_coldata_infected]

# Check all sample IDs in ColData are also in CountData and match their orders
all(rownames(PMG_coldata_infected) %in% colnames(Perkinsus_counts_infected))  #Should return TRUE
# returns TRUE
all(rownames(PMG_coldata_infected) == colnames(Perkinsus_counts_infected))   # should return TRUE
#returns TRUE! 

# CHECK LEVELS
# All injected 
levels(PMG_coldata_infected$FamCode) # "DA" "LB"
levels(PMG_coldata_infected$Treat) # "Control"  "Injected"

# MAKE DESEQ DATA SET FROM MATRIX
# Perkinsus all injected
Perkinsus_dds_infected <- DESeqDataSetFromMatrix(countData= Perkinsus_counts_infected,
                                                     colData = PMG_coldata_infected,
                                                     design = ~Library_Prep_Date + Treat + FamCode)
# PREFILTER THE DATA
Perkinsus_dds_infected <- Perkinsus_dds_infected[rowSums(counts(Perkinsus_dds_infected)) >10,] 
nrow(Perkinsus_dds_infected) # 5236

# VST TRANSFORMATION OF COUNTS FOR VISUALIZATION 
# Need to use the varianceStabilizingTransformation directly because nsub rows because the default is 1000 and there are only 738
Perkinsus_dds_infected_vsd <- varianceStabilizingTransformation(Perkinsus_dds_infected)

# Calculate sample distances to assess overall sample similarity
# use function dist to calculate the Euclidean distance between samples. 
# To ensure we have a roughly equal contribution from all genes, we use it on the VST data. 
# We need to transpose the matrix of values using t, because the dist function expects the different samples to be rows of its argument, and different dimensions (here, genes) to be columns.

Perkinsus_dds_infected_vsd_dist <- dist(t(assay(Perkinsus_dds_infected_vsd)))
Perkinsus_dds_infected_vsd_dist

# PCA plot visualization of individuals in the family #
plotPCA(Perkinsus_dds_infected_vsd, intgroup=c("LogConc")) # PC1 explains almost all variance, may be just correlated with load
plotPCA(Perkinsus_dds_infected_vsd, intgroup=c("Treat","FamCode","Timepoint"))

# Is there a correlation between LogConc and Number of reads? (this is an associated question) 

# Calculate Differential expression
Perkinsus_dds_infected_deseq <- DESeq(Perkinsus_dds_infected) 

# Inspect Results names object
resultsNames(Perkinsus_dds_infected_deseq) # [1] "Intercept" "Library_Prep_Date_17_Dec_vs_15_Dec" "Treat_Injected_PMG_vs_Control_PMG" "FamCode_LB_vs_DA"  

# Build Results object
Perkinsus_dds_infected_res <- results(Perkinsus_dds_infected_deseq, alpha=0.05, name="Treat_Injected_PMG_vs_Control_PMG")
Perkinsus_dds_infected_res # comparison is log2 fold change (MLE): Treat Injected PMG vs Control PMG 

# LFC Shrinkage with Apeglm
Perkinsus_dds_infected_res_LFC <- lfcShrink(Perkinsus_dds_infected_deseq, coef="Treat_Injected_PMG_vs_Control_PMG", type="apeglm", res=Perkinsus_dds_infected_res)

# Inspect results summary
#summary(Perkinsus_dds_infected_res_LFC)
#out of 5233 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 933, 18%
#LFC < 0 (down)     : 6, 0.11%
#outliers [1]       : 3, 0.057%
#low counts [2]     : 3554, 68%
(mean count < 1)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

# Susbet results with padj < 0.05 and Logfold change of >1.0 or <-1.0
# only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
Perkinsus_dds_infected_res_LFC_sig <- subset(Perkinsus_dds_infected_res_LFC, padj < 0.05)
Perkinsus_dds_infected_res_LFC_sig_less_1_kept <- Perkinsus_dds_infected_res_LFC_sig
Perkinsus_dds_infected_res_LFC_sig_less_1_kept$transcript_id <- row.names(Perkinsus_dds_infected_res_LFC_sig_less_1_kept )
Perkinsus_dds_infected_res_LFC_sig_less_1_kept<- as.data.frame(Perkinsus_dds_infected_res_LFC_sig_less_1_kept )
Perkinsus_dds_infected_res_LFC_sig_less_1_kept <- Perkinsus_dds_infected_res_LFC_sig_less_1_kept  %>% filter(abs(log2FoldChange) > 0.1)
nrow(Perkinsus_dds_infected_res_LFC_sig_less_1_kept) #95

Perkinsus_dds_infected_res_LFC_sig$transcript_id <- row.names(Perkinsus_dds_infected_res_LFC_sig )
Perkinsus_dds_infected_res_LFC_sig<- as.data.frame(Perkinsus_dds_infected_res_LFC_sig )
Perkinsus_dds_infected_res_LFC_sig <- Perkinsus_dds_infected_res_LFC_sig  %>% filter(abs(log2FoldChange) >= 1.0)
nrow(Perkinsus_dds_infected_res_LFC_sig ) #93
nrow(Perkinsus_dds_infected_res_LFC) # 5236

# Annotate transcripts
Perkinsus_dds_infected_res_LFC_sig_annot <- Perkinsus_dds_infected_res_LFC_sig %>% left_join(PMG_annotation_transcripts_unique[,c("transcript_id", "product")], by = "transcript_id")
head(arrange(Perkinsus_dds_infected_res_LFC_sig_annot,-log2FoldChange),n=50) # arrange in descending order

# Log2 Fold Change 
PMG_infected_log2foldchange <- ggplot(Perkinsus_dds_infected_res_LFC_sig_annot, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Significant P. marinus \nLog2 Fold Change Control vs. Injected, zero load removed") +
  ylab("Log2 Fold Change between Control and Injected PMG")

##### Perkinsus Differential Expression Day 7 control vs. injected (Zero load removed) ####
# Goal of this comparison is to see if there are differences in expression between low load and high,
# which mostly corresponds to the family and treatment groups 

# Rank samples by LOAD
PMG_coldata_infected <- PMG_coldata %>% filter(LogConc > 0.00)
PMG_coldata_infected[order(PMG_coldata_infected$LogConc),] 
row.names(PMG_coldata_infected) <- PMG_coldata_infected$X

PMG_coldata_infected_Day7 <- PMG_coldata_infected %>% filter(Timepoint == "7d")
row.names(PMG_coldata_infected_Day7) <- PMG_coldata_infected_Day7$X

# extract list of samples
SampleID_PMG_coldata_infected_Day7 <- row.names(PMG_coldata_infected_Day7)

# Extract those samples from the counts matrix 
Perkinsus_counts_infected_Day7 <- Perkinsus_counts_data[,names(Perkinsus_counts_data) %in% SampleID_PMG_coldata_infected_Day7]

# Check all sample IDs in ColData are also in CountData and match their orders
all(rownames(PMG_coldata_infected_Day7 ) %in% colnames(Perkinsus_counts_infected_Day7))  #Should return TRUE
# returns TRUE
all(rownames(PMG_coldata_infected_Day7 ) == colnames(Perkinsus_counts_infected_Day7))   # should return TRUE
#returns TRUE! 

# CHECK LEVELS
# All injected 
levels(PMG_coldata_infected_Day7$FamCode) # "DA" "LB"
levels(PMG_coldata_infected_Day7$Treat) # "Control"  "Injected"
levels(PMG_coldata_infected_Day7$Timepoint)
PMG_coldata_infected_Day7$Timepoint <- droplevels(PMG_coldata_infected_Day7$Timepoint)
levels(PMG_coldata_infected_Day7$Timepoint) #7d

# MAKE DESEQ DATA SET FROM MATRIX
# Perkinsus all injected
Perkinsus_dds_infected_Day7 <- DESeqDataSetFromMatrix(countData= Perkinsus_counts_infected_Day7,
                                                 colData = PMG_coldata_infected_Day7,
                                                 design = ~Library_Prep_Date + Treat + FamCode)
# PREFILTER THE DATA
Perkinsus_dds_infected_Day7 <- Perkinsus_dds_infected_Day7[rowSums(counts(Perkinsus_dds_infected_Day7)) >10,] 
nrow(Perkinsus_dds_infected_Day7) # 795

# VST TRANSFORMATION OF COUNTS FOR VISUALIZATION 
# Need to use the varianceStabilizingTransformation directly because nsub rows because the default is 1000 and there are only 738
Perkinsus_dds_infected_Day7_vsd <- varianceStabilizingTransformation(Perkinsus_dds_infected_Day7)

# Calculate sample distances to assess overall sample similarity
# use function dist to calculate the Euclidean distance between samples. 
# To ensure we have a roughly equal contribution from all genes, we use it on the VST data. 
# We need to transpose the matrix of values using t, because the dist function expects the different samples to be rows of its argument, and different dimensions (here, genes) to be columns.

Perkinsus_dds_infected_Day7_vsd_dist <- dist(t(assay(Perkinsus_dds_infected_Day7_vsd)))
Perkinsus_dds_infected_Day7_vsd_dist

# PCA plot visualization of individuals in the family #
plotPCA(Perkinsus_dds_infected_Day7_vsd, intgroup=c("LogConc")) # PC1 explains almost all variance, may be just correlated with load
plotPCA(Perkinsus_dds_infected_Day7_vsd, intgroup=c("FamCode"))

# Is there a correlation between LogConc and Number of reads? (this is an associated question) 

# Calculate Differential expression
Perkinsus_dds_infected_Day7_deseq <- DESeq(Perkinsus_dds_infected_Day7) 

# Inspect Results names object
resultsNames(Perkinsus_dds_infected_Day7_deseq ) # [1] "Intercept" "Library_Prep_Date_17_Dec_vs_15_Dec" "Treat_Injected_PMG_vs_Control_PMG" "FamCode_LB_vs_DA"  

# Build Results object
Perkinsus_dds_infected_Day7_res <- results(Perkinsus_dds_infected_Day7_deseq, alpha=0.05, name="Treat_Injected_PMG_vs_Control_PMG")
Perkinsus_dds_infected_Day7_res # comparison is log2 fold change (MLE): Treat Injected PMG vs Control PMG 

# LFC Shrinkage with Apeglm
Perkinsus_dds_infected_Day7_res_LFC <- lfcShrink(Perkinsus_dds_infected_Day7_deseq, coef="Treat_Injected_PMG_vs_Control_PMG", type="apeglm", res=Perkinsus_dds_infected_Day7_res)

# Inspect results summary
summary(Perkinsus_dds_infected_Day7_res_LFC)
#out of 795 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 12, 1.5%
#LFC < 0 (down)     : 1, 0.13%
#outliers [1]       : 8, 1%
#low counts [2]     : 216, 27%
#(mean count < 1)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

# Susbet results with padj < 0.05 and Logfold change of >1.0 or <-1.0
# only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
Perkinsus_dds_infected_Day7_res_LFC_sig <- subset(Perkinsus_dds_infected_Day7_res_LFC, padj < 0.05)
Perkinsus_dds_infected_Day7_res_LFC_sig_less_1_kept <- Perkinsus_dds_infected_Day7_res_LFC_sig
Perkinsus_dds_infected_Day7_res_LFC_sig_less_1_kept$transcript_id <- row.names(Perkinsus_dds_infected_Day7_res_LFC_sig_less_1_kept )
Perkinsus_dds_infected_Day7_res_LFC_sig_less_1_kept<- as.data.frame(Perkinsus_dds_infected_Day7_res_LFC_sig_less_1_kept )
Perkinsus_dds_infected_Day7_res_LFC_sig_less_1_kept <- Perkinsus_dds_infected_Day7_res_LFC_sig_less_1_kept %>% filter(abs(log2FoldChange) > 0.1)
nrow(Perkinsus_dds_infected_Day7_res_LFC_sig_less_1_kept) #3

Perkinsus_dds_infected_Day7_res_LFC_sig$transcript_id <- row.names(Perkinsus_dds_infected_Day7_res_LFC_sig)
Perkinsus_dds_infected_Day7_res_LFC_sig<- as.data.frame(Perkinsus_dds_infected_Day7_res_LFC_sig )
Perkinsus_dds_infected_Day7_res_LFC_sig <- Perkinsus_dds_infected_Day7_res_LFC_sig %>% filter(abs(log2FoldChange) >= 1.0)
nrow(Perkinsus_dds_infected_Day7_res_LFC_sig ) #3
nrow(Perkinsus_dds_infected_Day7_res_LFC) # 795

# Annotate transcripts
Perkinsus_dds_infected_Day7_res_LFC_sig_annot <- Perkinsus_dds_infected_Day7_res_LFC_sig %>% left_join(PMG_annotation_transcripts_unique[,c("transcript_id", "product")], by = "transcript_id")
head(arrange(Perkinsus_dds_infected_Day7_res_LFC_sig_annot,-log2FoldChange),n=50) # arrange in descending order

# Log2 Fold Change 
PMG_infected_Day7_log2foldchange <- ggplot(Perkinsus_dds_infected_Day7_res_LFC_sig_annot, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Significant P. marinus \nLog2 Fold Change Control vs. Injected Day7, zero load removed") +
  ylab("Log2 Fold Change between Control and Injected PMG")

##### Perkinsus Differential Expression LB vs DA (Zero load removed) Control in tbe formula ####
# Goal of this comparison is to see if there are differences in expression between low load and high,
# which mostly corresponds to the family and treatment groups 

# Rank samples by LOAD
PMG_coldata_infected <- PMG_coldata %>% filter(LogConc > 0.00)
PMG_coldata_infected[order(PMG_coldata_infected$LogConc),] 
row.names(PMG_coldata_infected) <- PMG_coldata_infected$X

# extract list of samples
SampleID_PMG_coldata_infected <- row.names(PMG_coldata_infected)

# Extract those samples from the counts matrix 
Perkinsus_counts_infected <- Perkinsus_counts_data[,names(Perkinsus_counts_data) %in% SampleID_PMG_coldata_infected]

# Check all sample IDs in ColData are also in CountData and match their orders
all(rownames(PMG_coldata_infected) %in% colnames(Perkinsus_counts_infected))  #Should return TRUE
# returns TRUE
all(rownames(PMG_coldata_infected) == colnames(Perkinsus_counts_infected))   # should return TRUE
#returns TRUE! 

# CHECK LEVELS
# All injected 
levels(PMG_coldata_infected$FamCode) # "DA" "LB"
levels(PMG_coldata_infected$Treat) # "Control"  "Injected"

# MAKE DESEQ DATA SET FROM MATRIX
# Perkinsus all injected
Perkinsus_dds_infected <- DESeqDataSetFromMatrix(countData= Perkinsus_counts_infected,
                                                 colData = PMG_coldata_infected,
                                                 design = ~Library_Prep_Date + Treat + FamCode)
# PREFILTER THE DATA
Perkinsus_dds_infected <- Perkinsus_dds_infected[rowSums(counts(Perkinsus_dds_infected)) >10,] 
nrow(Perkinsus_dds_infected) # 5236

# VST TRANSFORMATION OF COUNTS FOR VISUALIZATION 
# Need to use the varianceStabilizingTransformation directly because nsub rows because the default is 1000 and there are only 738
Perkinsus_dds_infected_vsd <- varianceStabilizingTransformation(Perkinsus_dds_infected)

# Calculate sample distances to assess overall sample similarity
# use function dist to calculate the Euclidean distance between samples. 
# To ensure we have a roughly equal contribution from all genes, we use it on the VST data. 
# We need to transpose the matrix of values using t, because the dist function expects the different samples to be rows of its argument, and different dimensions (here, genes) to be columns.

Perkinsus_dds_infected_vsd_dist <- dist(t(assay(Perkinsus_dds_infected_vsd)))
Perkinsus_dds_infected_vsd_dist

# PCA plot visualization of individuals in the family #
plotPCA(Perkinsus_dds_infected_vsd, intgroup=c("LogConc")) # PC1 explains almost all variance, may be just correlated with load
plotPCA(Perkinsus_dds_infected_vsd, intgroup=c("Treat","FamCode","Timepoint"))

# Is there a correlation between LogConc and Number of reads? (this is an associated question) 

# Calculate Differential expression
Perkinsus_dds_infected_deseq <- DESeq(Perkinsus_dds_infected) 

# Inspect Results names object
resultsNames(Perkinsus_dds_infected_deseq) # [1] "Intercept" "Library_Prep_Date_17_Dec_vs_15_Dec" "Treat_Injected_PMG_vs_Control_PMG" "FamCode_LB_vs_DA"  

# Build Results object
Perkinsus_dds_infected_LB_vs_DA_res <- results(Perkinsus_dds_infected_deseq, alpha=0.05, name="FamCode_LB_vs_DA")
Perkinsus_dds_infected_LB_vs_DA_res # comparison is log2 fold change (MLE): FamCode LB vs DA  

# LFC Shrinkage with Apeglm
Perkinsus_dds_infected_LB_vs_DA_res_LFC <- lfcShrink(Perkinsus_dds_infected_deseq, coef="FamCode_LB_vs_DA", type="apeglm", res=Perkinsus_dds_infected_LB_vs_DA_res)

# Inspect results summary
#summary(Perkinsus_dds_infected_res_LFC)
#out of 5233 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 933, 18%
#LFC < 0 (down)     : 6, 0.11%
#outliers [1]       : 3, 0.057%
#low counts [2]     : 3554, 68%
#(mean count < 1)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

# Susbet results with padj < 0.05 and Logfold change of >1.0 or <-1.0
# only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig <- subset(Perkinsus_dds_infected_LB_vs_DA_res_LFC, padj < 0.05)
Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig_less_1_kept <- Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig
Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig_less_1_kept$transcript_id <- row.names(Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig_less_1_kept )
Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig_less_1_kept<- as.data.frame(Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig_less_1_kept )
nrow(Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig_less_1_kept) #45

Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig$transcript_id <- row.names(Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig )
Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig<- as.data.frame(Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig)
Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig <- Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig  %>% filter(abs(log2FoldChange) >= 1.0)
nrow(Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig ) #15
nrow(Perkinsus_dds_infected_LB_vs_DA_res_LFC) # 5236

# Annotate transcripts
Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig_less_1_kept_annot <- Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig_less_1_kept %>% left_join(PMG_annotation_transcripts_unique[,c("transcript_id", "product")], by = "transcript_id")
head(arrange(Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig_less_1_kept_annot,-log2FoldChange),n=50) # arrange in descending order

Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig_annot <- Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig %>% left_join(PMG_annotation_transcripts_unique[,c("transcript_id", "product")], by = "transcript_id")
head(arrange(Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig_annot,-log2FoldChange),n=50) # arrange in descending order

# Log2 Fold Change 
PMG_infected_DA_LB_log2foldchange <- ggplot(Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig_less_1_kept_annot, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Significant P. marinus \nLog2 Fold Change LB vs. DA all timepoints, zero load removed") +
  ylab("Log2 Fold Change between LB vs. DA")
PMG_infected_DA_LB_subset_log2foldchange <- ggplot(Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig_annot, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Significant P. marinus \nLog2 Fold Change LB vs. DA all timepoints, zero load removed") +
  ylab("Log2 Fold Change between LB vs. DA")

Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig_annot[!(Perkinsus_dds_infected_LB_vs_DA_res_LFC_sig_annot %in% Perkinsus_dds_All_injected_res_LFC_sig_annot ),]

##### WGCNA Analysis between LB control vs challenge and DA control vs. challenge #####
library(WGCNA)
library(limma)
cor <- WGCNA::cor
options(stringsAsFactors = FALSE)
# WGCNA analysis resources
  # https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
  # https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/JMiller/
  # Input data is variance stabilizing transformation counts data 

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

#Load expression data
Dermo_counts

# Import subset coldata that has categorical variables as binary 
Dermo_coldata_binary <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/RAW DATA AND INFO/2015 Dermo Challenge /DATA/6h_36h_7d_DA_LB_binary.csv", header=TRUE )
head(Dermo_coldata_binary)
rownames(Dermo_coldata_binary) <- Dermo_coldata_binary$rownames
Dermo_coldata_binary <- Dermo_coldata_binary[,-1]
head(Dermo_coldata_binary)
class(Dermo_coldata_binary$FamCode)
class(Dermo_coldata_binary$Treat)
class(Dermo_coldata_binary$Timepoint)
class(Dermo_coldata_binary$LogConc)
Dermo_coldata_binary$FamCode <- as.numeric(Dermo_coldata_binary$FamCode)
Dermo_coldata_binary$Treat <- as.numeric(Dermo_coldata_binary$Treat)
Dermo_coldata_binary$Timepoint <- as.numeric(Dermo_coldata_binary$Timepoint)

# check c
all(rownames(Dermo_coldata_binary) %in% colnames(Dermo_counts))  # Should return TRUE
# returns TRUE
all(rownames(Dermo_coldata_binary) == colnames(Dermo_counts))    # Should return TRUE
# returns FALSE

# Visualize how sample data fits into dendrogram
  # Can't do because data is too large

# Load in Data as one vst transformed count matrix that have been filtered for low counts (see https://support.bioconductor.org/p/115583/)
Dermo_dds_family_WGCNA <- DESeqDataSetFromMatrix(countData = Dermo_counts,
                                           colData = Dermo_coldata, # need to use original, then reorder
                                           design = ~Library_Prep_Date + Treat +FamCode)
Dermo_dds_family_WGCNA <- Dermo_dds_family_WGCNA[ rowSums(counts(Dermo_dds_family_WGCNA)) > 20, ]

# reorder before VST 
# reorder coldata 
Dermo_dds_family_WGCNA_reorder <- Dermo_dds_family_WGCNA[, rownames(Dermo_coldata_binary)]

# recheck 
all(rownames(Dermo_coldata_binary) %in% colnames(Dermo_dds_family_WGCNA_reorder))  #Should return TRUE
# returns TRUE
all(rownames(Dermo_coldata_binary) == colnames(Dermo_dds_family_WGCNA_reorder))    # Should return TRUE
# TRUE

# perform VST analysis 
Dermo_dds_family_WGCNA_vsd <- vst(Dermo_dds_family_WGCNA, blind=FALSE)

# Remove Batch effect from library prep date with limma
colData(Dermo_dds_family_WGCNA)
plotPCA(Dermo_dds_family_WGCNA_vsd, "Library_Prep_Date")
plotPCA(Dermo_dds_family_WGCNA_vsd, "RNA_Extract_Date") # no effect
plotPCA(Dermo_dds_family_WGCNA_vsd, "RNA.Extracted_by")
plotPCA(Dermo_dds_family_WGCNA_vsd, "Library_Prep_by")
plotPCA(Dermo_dds_family_WGCNA_vsd, "FamCode")
plotPCA(Dermo_dds_family_WGCNA_vsd, "Sequence_Pool") # no effect
plotPCA(Dermo_dds_family_WGCNA_vsd, "TruSeq_Kit")
plotPCA(Dermo_dds_family_WGCNA_vsd, "Index") # no effect 

mat <- assay(Dermo_dds_family_WGCNA_vsd)
mat <- limma::removeBatchEffect(mat, Dermo_dds_family_WGCNA_vsd$Library_Prep_Date)
mat <- limma::removeBatchEffect(mat, Dermo_dds_family_WGCNA_vsd$RNA.Extracted_by)
mat <- limma::removeBatchEffect(mat, Dermo_dds_family_WGCNA_vsd$Library_Prep_by)
mat <- limma::removeBatchEffect(mat, Dermo_dds_family_WGCNA_vsd$TruSeq_Kit)

assay(Dermo_dds_family_WGCNA_vsd) <- mat
class(Dermo_dds_family_WGCNA_vsd)
Dermo_dds_family_WGCNA_vsd_matrix <- assay(Dermo_dds_family_WGCNA_vsd)
plotPCA(Dermo_dds_family_WGCNA_vsd, "Library_Prep_Date")
plotPCA(Dermo_dds_family_WGCNA_vsd, "RNA_Extract_Date") # no effect
plotPCA(Dermo_dds_family_WGCNA_vsd, "RNA.Extracted_by")
plotPCA(Dermo_dds_family_WGCNA_vsd, "Library_Prep_by")
plotPCA(Dermo_dds_family_WGCNA_vsd, "FamCode")
plotPCA(Dermo_dds_family_WGCNA_vsd, "Sequence_Pool") # no effect
plotPCA(Dermo_dds_family_WGCNA_vsd, "TruSeq_Kit")
plotPCA(Dermo_dds_family_WGCNA_vsd, "Index") # no effect 

# Expression data the columns should be genes and rows should be sample
head(Dermo_dds_family_WGCNA_vsd_matrix)
class(Dermo_dds_family_WGCNA_vsd_matrix)
Dermo_dds_family_WGCNA_vsd_matrix_transpose <- t(Dermo_dds_family_WGCNA_vsd_matrix)
Dermo_dds_family_WGCNA <- Dermo_dds_family_WGCNA_vsd_matrix_transpose

# export data to be loaded later
write.table(Dermo_dds_family_WGCNA, file="Dermo_dds_family_WGCNA.txt")
save(Dermo_dds_family_WGCNA, Dermo_coldata_binary, file = "Dermo_WGCNA.RData")

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
FamCode_lightcyan_annot_apop

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

# Identifying genes with high GS and MM
# MElightcyan is the module with the greatest positive correlation with Family
# MEturqoise is the module with greatest negative correlation
module="orange"
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
# association not significant 

### Create data frame holding information about all probes ###
probes=colnames(Dermo_dds_family_WGCNA)
probes2annot = match(probes, C_vir_rtracklayer_transcripts_GO$transcript_id)
geneTraitSignificance = as.data.frame(cor(Dermo_dds_family_WGCNA, Treat, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(Treat), sep="")
names(GSPvalue) = paste("p.GS.", names(Treat), sep="")

# Combine with module color, gene significance for weight, and module membership and p-values in all modules
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
gene_info_combined_significance_Treat_apop <- 

# plot gene significance of significant genes
ggplot(gene_info_combined_significance_Treat, aes(x=transcript_id, y=))

#### Exporting WGCNA results to Cytoscape ####

# Cytoscape wants specific modules to be input into the software. Need to select particular modules of interest

# Recalculate topological overlap using TOMsimilarityFromExpr in bluewaves

  ## in terminal
  ## scp Dermo_WGCNA.RData erin_roberts@bluewaves:/data3/marine_diseases_lab/erin/2015_Dermo_transcriptome_analysis
  #  $ cd /data3/marine_diseases_lab/erin/2015_Dermo_transcriptome_analysis
  #  $ interactive
  #  $ module load R/3.3.1-foss-2016b
  # $ module load R-bundle-Bioconductor/3.3-foss-2016b-R-3.3.1
  # installed dependecies and installed WGCNA using 
      # > install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival")) 
      #  > source("http://bioconductor.org/biocLite.R") 
      #  > biocLite(c("GO.db", "preprocessCore", "impute", "WGCNA"))
  # Wrote script to run in R via sbatch script
          # library(WGCNA)
          # lnames = load(file="Dermo_WGCNA.Rdata")
          # TOM = TOMsimilarityFromExpr(Dermo_dds_family_WGCNA, power=12)
          # TOM <- as.matrix(TOM)
          # write.txt(TOM, file ="Dermo_TOM.txt", row.names=TRUE)

  

# Select modules
modules = c("lightcyan", "paleturqoise", "green", "orange","plum1")

# Select module probes
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]
modGenes = C_vir_rtracklayer_transcripts_GO$product[match(modProbes, C_vir_rtracklayer_transcripts_GO$transcript_id)]

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])
