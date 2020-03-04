# Script to Analyze C. virginica and C. gigas orthologs involved in apoptosis
# Orthologs for apoptosis will be input into PCA and WGCNA
# Erin Roberts, PhD Candidate University of Rhode Island 
# 2/26/2020


# Helpful tutorial for OrthoFinder
# OrthoFinder Manual: https://github.com/davidemms/OrthoFinder/blob/master/README.md
# http://arken.nmbu.no/~larssn/teach/bin310/week8.html

# Load packages 
library(tidyverse)
library(UpSetR)
library(ComplexHeatmap)

#### LOAD ORTHOGROUP INFORMATION FROM ORTHOFINDER####
# OrthoFinder 2.3.3 was used to acquire OrthoFinder information 
# Load Orthogroups.tsv containing proteins in each orthogroup,Orthogroups_SingleCopyOrthologues.txt which lists the 1:1 orthologs,
  # and Orthogroups.GeneCount.tsv which lists the number of genes in each orthogroup

Single_copy_orthogroups <- read_table("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/OrthoFinder/Results_Feb25/Orthogroups/Orthogroups_SingleCopyOrthologues.txt",
                                      col_names = "Orthogroup")
Single_copy_orthogroups <- as.data.frame(Single_copy_orthogroups)
class(Single_copy_orthogroups)

# Load and Parse orthogroup tsv
Orthogroups <- read_tsv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter_1_Apoptosis_Annotation_Data_Analyses_2019/DATA/OrthoFinder/Results_Feb25/Orthogroups/Orthogroups.tsv", 
              col_names = c("Orthogroup", "C_gigas","C_virginica" ))
# Drop groups that have only one-species.
Orthogroups <- Orthogroups %>% drop_na()

# Split into two data frames
Orthogroups_Cgig <- Orthogroups[,c(1:2)]
Orthogroups_Cvir <- Orthogroups[,c(1,3)]

G <- strsplit(Orthogroups_Cgig$C_gigas, split = ",")
C_gig_orthogroups <- data.frame(Orthogroup = rep(Orthogroups_Cgig$Orthogroup, sapply(G, length)), C_gigas = unlist(G))

C <- strsplit(Orthogroups_Cvir$C_virginica, split = ",")
C_vir_orthogroups <-data.frame(Orthogroup = rep(Orthogroups_Cvir$Orthogroup, sapply(C, length)), C_virginica = unlist(C))

#### ISOLATE SINGLE COPY ORTHOLOGS, XPs FROM ORTHOGROUPS #### 
# NOTE: decided not to try to match the XPs to the XMs. It is too big of a jump to say that because XPs are orthologous that the XMs are orthologous 

## Isolate single copy orthologs from each table
C_gig_orthogroups_single_copy <- C_gig_orthogroups[C_gig_orthogroups$Orthogroup %in% Single_copy_orthogroups$Orthogroup,]
C_vir_orthogroups_single_copy <- C_vir_orthogroups[C_vir_orthogroups$Orthogroup %in% Single_copy_orthogroups$Orthogroup,]

## Isolate XPs and gene info from C_vir_rtracklayer and C_gig_rtracklayer
colnames(C_gig_orthogroups_single_copy)[2] <- "protein_id"
C_gig_orthogroups_single_copy <- left_join(C_gig_orthogroups_single_copy, C_gig_rtracklayer, by = "protein_id")

colnames(C_vir_orthogroups_single_copy) [2] <- "protein_id"
C_vir_orthogroups_single_copy  <- left_join(C_vir_orthogroups_single_copy ,  C_vir_rtracklayer, by="protein_id")

##### LOAD GENE_COUNT_MATRICES ####

# Load gene_count_matrix.csv
Dermo_gene_counts <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/Dermo_gene_count_matrix.csv", row.names = 1 )
Probiotic_gene_counts <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/Probiotic_gene_count_matrix.csv", row.names = 1 )
Pro_RE22_gene_counts <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/Probiotic_RE22_gene_count_matrix.csv", row.names = 1 )
Rubio_gene_counts <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/Rubio_gene_count_matrix.csv", row.names = 1 )
deLorgeril_gene_counts <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/deLorgeril_gene_count_matrix.csv", row.names = 1 )
He_gene_counts <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/He_gene_count_matrix.csv", row.names = 1 )
ROD_gene_counts <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/ROD_gene_count_matrix.csv", row.names = 1 )
Zhang_gene_counts <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/Zhang_gene_count_matrix.csv", row.names = 1 )

head(Dermo_gene_counts )
head(Probiotic_gene_counts) # header need the "_1" removed
head(Pro_RE22_gene_counts) # header need the "_1" removed
head(Rubio_gene_counts) # header need the "_1" removed
head(deLorgeril_gene_counts) # header need the "_1" removed
head(He_gene_counts )
head(ROD_gene_counts)
head(Zhang_gene_counts )

# # colnames all have "_1", remove this. It was an artifact of the PE sample names
colnames(Probiotic_gene_counts ) <- sub('\\_[^_]+$', '', colnames(Probiotic_gene_counts))
colnames(Rubio_gene_counts ) <- sub('\\_[^_]+$', '', colnames(Rubio_gene_counts))
colnames(deLorgeril_gene_counts) <- sub('\\_[^_]+$', '', colnames(deLorgeril_gene_counts))
colnames(Pro_RE22_gene_counts) <- sub('\\_[^_]+$', '', colnames(Pro_RE22_gene_counts))

# C_vir gene file headers are gene#|LOC...
head(Dermo_gene_counts)
head(Probiotic_gene_counts)
head(ROD_gene_counts)
head(Pro_RE22_gene_counts)

# C_gig gene file headers are gene-LOC|LOC..
head(Rubio_gene_counts)
head(deLorgeril_gene_counts)
head(He_gene_counts)
head(Zhang_gene_counts)

# remove MSTRG novel gene lines (can assess these later if necessary)
Dermo_gene_counts <- Dermo_gene_counts[!grepl("MSTRG", row.names(Dermo_gene_counts)),]
Probiotic_gene_counts  <- Probiotic_gene_counts [!grepl("MSTRG", row.names(Probiotic_gene_counts )),]
Pro_RE22_gene_counts <- Pro_RE22_gene_counts[!grepl("MSTRG", row.names(Pro_RE22_gene_counts)),]
Rubio_gene_counts <- Rubio_gene_counts[!grepl("MSTRG", row.names(Rubio_gene_counts)),]
deLorgeril_gene_counts <- deLorgeril_gene_counts[!grepl("MSTRG", row.names(deLorgeril_gene_counts)),]
He_gene_counts <- He_gene_counts[!grepl("MSTRG", row.names(He_gene_counts)),]
ROD_gene_counts <- ROD_gene_counts[!grepl("MSTRG", row.names(ROD_gene_counts)),]
Zhang_gene_counts <- Zhang_gene_counts[!grepl("MSTRG", row.names(Zhang_gene_counts)),]

# add row.names to new column for merging 
Dermo_gene_counts$gene<- row.names(Dermo_gene_counts)
Probiotic_gene_counts$gene <- row.names(Probiotic_gene_counts)
ROD_gene_counts$gene <- row.names(ROD_gene_counts)
Pro_RE22_gene_counts$gene <- row.names(Pro_RE22_gene_counts)
Rubio_gene_counts$gene <- row.names(Rubio_gene_counts)
deLorgeril_gene_counts$gene <- row.names(deLorgeril_gene_counts)
He_gene_counts$gene <- row.names(He_gene_counts)
Zhang_gene_counts$gene <- row.names(Zhang_gene_counts)

# Remove everything before | 
remove_gene_before = function(x){
  return(gsub(".*\\|","",x)) # need to escape the | symbol first
}

Dermo_gene_counts$gene <- remove_gene_before(Dermo_gene_counts$gene)
Probiotic_gene_counts$gene <- remove_gene_before(Probiotic_gene_counts$gene)
ROD_gene_counts$gene<- remove_gene_before(ROD_gene_counts$gene)
Pro_RE22_gene_counts$gene<- remove_gene_before(Pro_RE22_gene_counts$gene)
Rubio_gene_counts$gene <- remove_gene_before(Rubio_gene_counts$gene)
deLorgeril_gene_counts$gene<- remove_gene_before(deLorgeril_gene_counts$gene)
He_gene_counts$gene <- remove_gene_before(He_gene_counts$gene)
Zhang_gene_counts$gene <- remove_gene_before(Zhang_gene_counts$gene)

#### SUBSET GENE COUNT MATRICES FOR ORTHOLOGS AND COMBINE #####

C_vir_orthogroups_single_copy_unique <- unique(C_vir_orthogroups_single_copy[,c("gene","Orthogroup")])
C_gig_orthogroups_single_copy_unique <- unique(C_gig_orthogroups_single_copy[,c("gene","Orthogroup")])

Dermo_gene_counts_single_copy <- Dermo_gene_counts[Dermo_gene_counts$gene %in% C_vir_orthogroups_single_copy_unique$gene,]
Dermo_gene_counts_single_copy <- left_join(Dermo_gene_counts_single_copy, C_vir_orthogroups_single_copy_unique)
nrow(Dermo_gene_counts_single_copy) # 347

Probiotic_gene_counts_single_copy <- Probiotic_gene_counts[Probiotic_gene_counts$gene %in% C_vir_orthogroups_single_copy_unique$gene,]
Probiotic_gene_counts_single_copy <- left_join(Probiotic_gene_counts_single_copy, C_vir_orthogroups_single_copy_unique)
nrow(Probiotic_gene_counts_single_copy) #772

Pro_RE22_gene_counts_single_copy <- Pro_RE22_gene_counts[Pro_RE22_gene_counts$gene %in% C_vir_orthogroups_single_copy_unique$gene,]
Pro_RE22_gene_counts_single_copy <- left_join(Pro_RE22_gene_counts_single_copy, C_vir_orthogroups_single_copy_unique)
nrow(Pro_RE22_gene_counts_single_copy) #610

ROD_gene_counts_single_copy <- ROD_gene_counts[ROD_gene_counts$gene %in% C_vir_orthogroups_single_copy_unique$gene,]
ROD_gene_counts_single_copy <- left_join(ROD_gene_counts_single_copy, C_vir_orthogroups_single_copy_unique)
nrow(ROD_gene_counts_single_copy) #625

Rubio_gene_counts_single_copy  <- Rubio_gene_counts[Rubio_gene_counts$gene %in% C_gig_orthogroups_single_copy_unique$gene,]
Rubio_gene_counts_single_copy  <- left_join(Rubio_gene_counts_single_copy, C_gig_orthogroups_single_copy_unique)
nrow(Rubio_gene_counts_single_copy ) #707

deLorgeril_gene_counts_single_copy <- deLorgeril_gene_counts[deLorgeril_gene_counts$gene %in% C_gig_orthogroups_single_copy_unique$gene,]
deLorgeril_gene_counts_single_copy <- left_join(deLorgeril_gene_counts_single_copy, C_gig_orthogroups_single_copy_unique)
nrow(deLorgeril_gene_counts_single_copy ) #557

He_gene_counts_single_copy <- He_gene_counts[He_gene_counts$gene %in% C_gig_orthogroups_single_copy_unique$gene,]
He_gene_counts_single_copy <- left_join(He_gene_counts_single_copy, C_gig_orthogroups_single_copy_unique)
nrow(He_gene_counts_single_copy ) #1071

Zhang_gene_counts_single_copy <- Zhang_gene_counts[Zhang_gene_counts$gene %in% C_gig_orthogroups_single_copy_unique$gene,]
Zhang_gene_counts_single_copy <- left_join(Zhang_gene_counts_single_copy, C_gig_orthogroups_single_copy_unique)
nrow(Zhang_gene_counts_single_copy ) #1808

# Join within species based on gene and orthogroup, starting with one that has the least (only can compare ones shared by all in expression)
C_vir_ortholog_gene_counts <- left_join(Dermo_gene_counts_single_copy,ROD_gene_counts_single_copy, by =c("gene","Orthogroup"))
C_vir_ortholog_gene_counts <- left_join(C_vir_ortholog_gene_counts ,Pro_RE22_gene_counts_single_copy, by =c("gene","Orthogroup"))
C_vir_ortholog_gene_counts <- left_join(C_vir_ortholog_gene_counts,Probiotic_gene_counts_single_copy, by =c("gene","Orthogroup"))
nrow(C_vir_ortholog_gene_counts) # 347 

C_gig_ortholog_gene_counts <- left_join(deLorgeril_gene_counts_single_copy,  Rubio_gene_counts_single_copy , by =c("gene","Orthogroup"))
C_gig_ortholog_gene_counts <- left_join(C_gig_ortholog_gene_counts, He_gene_counts_single_copy , by =c("gene","Orthogroup"))
C_gig_ortholog_gene_counts <- left_join(C_gig_ortholog_gene_counts, Zhang_gene_counts_single_copy, by =c("gene","Orthogroup"))
nrow(C_gig_ortholog_gene_counts) # 557

## Merge data frames based on matching orthologroup IDs to get full table of counts 
Full_ortholog_gene_count <- left_join(C_vir_ortholog_gene_counts,C_gig_ortholog_gene_counts,  by ="Orthogroup")
nrow(Full_ortholog_gene_count) # 347

# set colnames as the unique orthogroup values (there are some duplicated gene names)
row.names(Full_ortholog_gene_count) <- Full_ortholog_gene_count$Orthogroup
colnames(Full_ortholog_gene_count)
# remove orthogroup and gene name groups so I only have the counts in the matrix 
Full_ortholog_gene_count_only <- Full_ortholog_gene_count[,c(1:97,100:177,179:237)]
ncol(Full_ortholog_gene_count_only) #234
Full_ortholog_gene_count_only[is.na(Full_ortholog_gene_count_only)] <-0
head(Full_ortholog_gene_count_only)

#### Ortholog Upset plot ####

# Using UpsetR
# Input data can be a list of sets where each set is a vector: https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html#input-data
Dermo_gene_counts_single_copy_vector <- as.vector(row.names(Dermo_gene_counts_single_copy))
Probiotic_gene_counts_single_copy_vector <- as.vector(row.names(Probiotic_gene_counts_single_copy))
ROD_gene_counts_single_copy_vector <- as.vector(row.names(ROD_gene_counts_single_copy))
Pro_RE22_gene_counts_single_copy_vector <- as.vector(row.names(Pro_RE22_gene_counts_single_copy))
Rubio_gene_counts_single_copy_vector <- as.vector(row.names(Rubio_gene_counts_single_copy))
deLorgeril_gene_counts_single_copy_vector <- as.vector(row.names(deLorgeril_gene_counts_single_copy))
He_gene_counts_single_copy_vector <- as.vector(row.names(He_gene_counts_single_copy))
Zhang_gene_counts_single_copy_vector <- as.vector(row.names(Zhang_gene_counts_single_copy))

Ortholog_list <- 
  list(Dermo_gene_counts_single_copy_vector = Dermo_gene_counts_single_copy_vector,
       Probiotic_gene_counts_single_copy_vector = Probiotic_gene_counts_single_copy_vector,
       ROD_gene_counts_single_copy_vector  =ROD_gene_counts_single_copy_vector ,
       Pro_RE22_gene_counts_single_copy_vector  =Pro_RE22_gene_counts_single_copy_vector ,
       Rubio_gene_counts_single_copy_vector =Rubio_gene_counts_single_copy_vector,
       deLorgeril_gene_counts_single_copy_vector =deLorgeril_gene_counts_single_copy_vector,
       He_gene_counts_single_copy_vector  =He_gene_counts_single_copy_vector ,
       Zhang_gene_counts_single_copy_vector =Zhang_gene_counts_single_copy_vector
     )

Ortholog_list <- list_to_matrix(Ortholog_list)


##### DIFFERENTIAL EXPRESSION OF GENE COUNT MATRICES FOR EACH EXPERIMENT ######

#### Zhang ####
# For each experiment, the same models are going to be used as before (ensure that I'm happy with them all before proceeding)
# Zhang 
Zhang_gene_counts
ncol(Zhang_gene_counts)
colnames(Zhang_gene_counts)
Zhang_gene_counts <- Zhang_gene_counts[,-10]
colnames(Zhang_gene_counts)
Zhang_coldata <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/Zhang_coldata2.csv", row.names = 1 )
View(Zhang_coldata)  
nrow(Zhang_coldata) 

# Make sure the columns of the count matrix and rows of the column data (sample metadata) are in the same order. 
all(rownames(Zhang_coldata) %in% colnames(Zhang_gene_counts))  #Should return TRUE
# returns TRUE
all(colnames(Zhang_gene_counts) %in% rownames(Zhang_coldata))  
# returns TRUE
all(rownames(Zhang_coldata) == colnames(Zhang_gene_counts))    # should return TRUE
# returns TRUE

### DATA QC PCA PLOT 
Zhang_gene_counts_matrix <- as.matrix(Zhang_gene_counts)
Zhangrlog_gene_counts <- rlog(Zhang_gene_counts_matrix, blind =TRUE)
# run PCA
pcZhang_gene <- prcomp(t(Zhangrlog_gene_counts))

# colour to look at the clustering 
autoplot(pcZhang_gene,
         data = Zhang_coldata, 
         colour="condition", 
         size=5) #
autoplot(pcZhang,
         data = Zhang_coldata, 
         colour="group_by_sim", 
         size=5) 

## MAKE DESEQ DATA SET FROM MATRIX
## Creating three here so I can compare the results
Zhang_gene_dds <- DESeqDataSetFromMatrix(countData = Zhang_gene_counts,
                                    colData = Zhang_coldata,
                                    design = ~time + group_by_sim) # add time to control for injection and time effect
## Prefiltering the data
Zhang_gene_dds <- Zhang_gene_dds[ rowSums(counts(Zhang_gene_dds)) > 10, ]

## Check levels 
levels(Zhang_coldata$condition) #"Control" 
levels(Zhang_coldata$group_by_sim) #  "control"             "LPS_M_lut"           "V_aes_V_alg1_V_alg2" "V_tub_V_ang"   
levels(Zhang_coldata$time) # 12hr, no injection
Zhang_gene_dds$time <- factor(Zhang_gene_dds$time , levels = c("No_injection","12h"))

## DATA TRANSFORMATION AND VISUALIZATION
Zhang_gene_dds_rlog <- rlog(Zhang_gene_dds, blind = TRUE) # keep blind = true before deseq function has been run

## PCA plot visualization of individuals in the family 
plotPCA(Zhang_gene_dds_rlog, intgroup=c("group_by_sim", "condition"))

### DIFFERENTIAL EXPRESSION ANALYSIS
Zhang_gene_dds_deseq <- DESeq(Zhang_gene_dds) 

## Check the resultsNames object of each to look at the available coefficients for use in lfcShrink command
resultsNames(Zhang_gene_dds_deseq)  # [1] "Intercept" "time_12h_vs_No_injection" "group_by_sim_LPS_M_lut_vs_control"          
# "group_by_sim_V_aes_V_alg1_V_alg2_vs_control" "group_by_sim_V_tub_V_ang_vs_contro  

## BUILD THE RESULTS OBJECT
# Examining the results object, change alpha to p <0.05, looking at object metadata
Zhang_gene_dds_deseq_res_V_alg1 <- results(Zhang_gene_dds_deseq, alpha=0.1, name = "group_by_sim_V_aes_V_alg1_V_alg2_vs_control"  )
Zhang_gene_dds_deseq_res_V_tub <- results(Zhang_gene_dds_deseq, alpha=0.1, name= "group_by_sim_V_tub_V_ang_vs_control" )
Zhang_gene_dds_deseq_res_LPS <- results(Zhang_gene_dds_deseq, alpha=0.1, name= "group_by_sim_LPS_M_lut_vs_control")

### Perform LFC Shrinkage with apeglm
Zhang_gene_dds_deseq_res_V_alg1_LFC <- lfcShrink(Zhang_gene_dds_deseq, coef="group_by_sim_V_aes_V_alg1_V_alg2_vs_control" , type= "apeglm", res=Zhang_gene_dds_deseq_res_V_alg1)
summary(Zhang_gene_dds_deseq_res_V_alg1_LFC) # 1 significant gene

Zhang_gene_dds_deseq_res_V_tub_LFC <- lfcShrink(Zhang_gene_dds_deseq, coef="group_by_sim_V_tub_V_ang_vs_control" , type= "apeglm", res=Zhang_gene_dds_deseq_res_V_tub)
summary(Zhang_gene_dds_deseq_res_V_tub_LFC) # 0

Zhang_gene_dds_deseq_res_LFC_LPS <- lfcShrink(Zhang_gene_dds_deseq, coef="group_by_sim_LPS_M_lut_vs_control", type= "apeglm", res=Zhang_gene_dds_deseq_res_LPS)
summary(Zhang_gene_dds_deseq_res_LFC_LPS) # 0

#### RUBIO ####
Rubio_gene_counts
colnames(Rubio_gene_counts)
Rubio_gene_counts <- Rubio_gene_counts[,-19]
Rubio_coldata <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/Rubio_coldata.csv", row.names = 1 )
View(Rubio_coldata)  
nrow(Rubio_coldata) 

# Make sure the columns of the count matrix and rows of the column data (sample metadata) are in the same order. 
all(rownames(Rubio_coldata) %in% colnames(Rubio_gene_counts ))  #Should return TRUE
# returns TRUE
all(colnames(Rubio_gene_counts) %in% rownames(Rubio_coldata))  
# returns TRUE
all(rownames(Rubio_coldata) == colnames(Rubio_gene_counts))    # should return TRUE
# returns FALSE

# Fix the order
Rubio_gene_counts <-Rubio_gene_counts[,row.names(Rubio_coldata)]
row.names(Rubio_coldata)

all(rownames(Rubio_coldata) %in% colnames(Rubio_gene_counts ))  #Should return TRUE
# returns TRUE
all(colnames(Rubio_gene_counts) %in% rownames(Rubio_coldata))  
# returns TRUE
all(rownames(Rubio_coldata) == colnames(Rubio_gene_counts))    # should return TRUE
# returns FALSE

### DATA QC PCA PLOT 
# rlog transform data is recommended over vst for small data sets 
# PCA plots of data (https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/02_Preprocessing_Data.nb.html#count-distribution-boxplots)
Rubio_gene_counts_matrix <- as.matrix(Rubio_gene_counts)
Rubiorlog_gene_counts <- rlog(Rubio_gene_counts_matrix, blind =TRUE)
#-- note: fitType='parametric', but the dispersion trend was not well captured by the
#function: y = a/x + b, and a local regression fit was automatically substituted.
#specify fitType='local' or 'mean' to avoid this message next time.

# run PCA
pcRubio_gene <- prcomp(t(Rubiorlog_gene_counts))

# Plot PCA
autoplot(pcRubio_gene,
         data = Rubio_coldata, 
         colour="Condition", 
         size=5) # PCA axes explain little of the variation 

## MAKE DESEQ DATA SET FROM MATRIX
## Check levels 
levels(Rubio_coldata$Condition) # "Control_anesthesis" "Control_untreated"  "Vcrass_J2_8"        "Vcrass_J2_9"        "Vtasma_LGP32"       "Vtasma_LMG20012T"  
Rubio_coldata$Condition <- factor(Rubio_coldata$Condition , levels = c("Control_untreated","Control_anesthesis", "Vcrass_J2_8", "Vcrass_J2_9","Vtasma_LGP32", "Vtasma_LMG20012T"  ))
levels(Rubio_coldata$Condition)
levels(Rubio_coldata$Group) # "Control"      "Non_virulent" "Virulent"  

## Creating three here so I can compare the results
Rubio_gene_dds <- DESeqDataSetFromMatrix(countData = Rubio_gene_counts,
                                    colData = Rubio_coldata,
                                    design = ~ Condition) 

## Prefiltering the data
Rubio_gene_dds <- Rubio_gene_dds [ rowSums(counts(Rubio_gene_dds )) > 10, ]

## DATA TRANSFORMATION AND VISUALIZATION
Rubio_gene_dds_rlog <- rlog(Rubio_gene_dds, blind = TRUE) # keep blind = true before deseq function has been run

## PCA plot visualization of individuals in the family 
plotPCA(Rubio_gene_dds_rlog, intgroup=c("Sample", "Condition")) # more of the variation explained

### DIFFERENTIAL EXPRESSION ANALYSIS
Rubio_gene_dds_deseq <- DESeq(Rubio_gene_dds) 
resultsNames(Rubio_gene_dds_deseq)  
#[1] "Intercept"                                         "Condition_Control_anesthesis_vs_Control_untreated"
#[3] "Condition_Vcrass_J2_8_vs_Control_untreated"        "Condition_Vcrass_J2_9_vs_Control_untreated"       
#[5] "Condition_Vtasma_LGP32_vs_Control_untreated"       "Condition_Vtasma_LMG20012T_vs_Control_untreated"

## BUILD THE RESULTS OBJECT
Rubio_gene_dds_deseq_J2_8_res <- results(Rubio_gene_dds_deseq, alpha=0.05, name="Condition_Vcrass_J2_8_vs_Control_untreated" )
Rubio_gene_dds_deseq_J2_9_res <- results(Rubio_gene_dds_deseq, alpha=0.05, name="Condition_Vcrass_J2_9_vs_Control_untreated"  )
Rubio_gene_dds_deseq_LGP32_res <- results(Rubio_gene_dds_deseq, alpha=0.05, name="Condition_Vtasma_LGP32_vs_Control_untreated"   )
Rubio_gene_dds_deseq_LMG20012T_res <- results(Rubio_gene_dds_deseq, alpha=0.05, name="Condition_Vtasma_LMG20012T_vs_Control_untreated"  )

### Perform LFC Shrinkage with apeglm
Rubio_gene_dds_deseq_J2_8_res_LFC<- lfcShrink(Rubio_gene_dds_deseq, coef="Condition_Vcrass_J2_8_vs_Control_untreated", type="apeglm", res= Rubio_gene_dds_deseq_J2_8_res)
Rubio_gene_dds_deseq_J2_9_res_LFC  <- lfcShrink(Rubio_gene_dds_deseq, coef="Condition_Vcrass_J2_9_vs_Control_untreated" , type= "apeglm", res= Rubio_gene_dds_deseq_J2_9_res   )
Rubio_gene_dds_deseq_LGP32_res_LFC<- lfcShrink(Rubio_gene_dds_deseq, coef="Condition_Vtasma_LGP32_vs_Control_untreated" , type= "apeglm", res= Rubio_gene_dds_deseq_LGP32_res   )
Rubio_gene_dds_deseq_LMG20012T_res_LFC<- lfcShrink(Rubio_gene_dds_deseq, coef="Condition_Vtasma_LMG20012T_vs_Control_untreated", type= "apeglm", res= Rubio_gene_dds_deseq_LMG20012T_res )

summary(Rubio_gene_dds_deseq_J2_8_res_LFC) 
summary(Rubio_gene_dds_deseq_J2_9_res_LFC  )
summary(Rubio_gene_dds_deseq_LGP32_res_LFC)
summary(Rubio_gene_dds_deseq_LMG20012T_res_LFC)

## EXPLORATORY PLOTTING OF RESULTS 
## MA Plotting
plotMA(Rubio_gene_dds_deseq_J2_8_res_LFC, ylim = c(-5, 5))
plotMA(Rubio_gene_dds_deseq_J2_9_res_LFC, ylim = c(-5, 5))
plotMA(Rubio_gene_dds_deseq_LGP32_res_LFC, ylim = c(-5, 5))
plotMA(Rubio_gene_dds_deseq_LMG20012T_res_LFC, ylim = c(-5, 5))
## Histogram of P values 
# exclude genes with very small counts to avoid spikes and plot using the LFCshrinkage
hist(Rubio_gene_dds_deseq_J2_8_res_LFC$padj[Rubio_gene_dds_deseq_J2_8_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
hist(Rubio_gene_dds_deseq_J2_9_res_LFC$padj[Rubio_gene_dds_deseq_J2_9_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
hist(Rubio_gene_dds_deseq_LGP32_res_LFC$padj[Rubio_gene_dds_deseq_LGP32_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
hist(Rubio_gene_dds_deseq_LMG20012T_res_LFC$padj[Rubio_gene_dds_deseq_LMG20012T_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

### Subsetting Significant Genes by padj < 0.05
Rubio_gene_dds_deseq_J2_8_res_LFC_sig <-  subset(Rubio_gene_dds_deseq_J2_8_res_LFC , padj < 0.05)
Rubio_gene_dds_deseq_J2_9_res_LFC_sig <-  subset(Rubio_gene_dds_deseq_J2_9_res_LFC , padj < 0.05)
Rubio_gene_dds_deseq_LGP32_res_LFC_sig <-  subset(Rubio_gene_dds_deseq_LGP32_res_LFC , padj < 0.05)
Rubio_gene_dds_deseq_LMG20012T_res_LFC_sig <-  subset(Rubio_gene_dds_deseq_LMG20012T_res_LFC, padj < 0.05)

Rubio_gene_dds_deseq_J2_8_res_LFC_sig$transcript_id <- row.names(Rubio_gene_dds_deseq_J2_8_res_LFC_sig)
Rubio_gene_dds_deseq_J2_9_res_LFC_sig$transcript_id <- row.names(Rubio_gene_dds_deseq_J2_9_res_LFC_sig)
Rubio_gene_dds_deseq_LGP32_res_LFC_sig$transcript_id <- row.names(Rubio_gene_dds_deseq_LGP32_res_LFC_sig)
Rubio_gene_dds_deseq_LMG20012T_res_LFC_sig$transcript_id <- row.names(Rubio_gene_dds_deseq_LMG20012T_res_LFC_sig)

Rubio_gene_dds_deseq_J2_8_res_LFC_sig <- as.data.frame(Rubio_gene_dds_deseq_J2_8_res_LFC_sig)
Rubio_gene_dds_deseq_J2_9_res_LFC_sig <- as.data.frame(Rubio_gene_dds_deseq_J2_9_res_LFC_sig)
Rubio_gene_dds_deseq_LGP32_res_LFC_sig <- as.data.frame(Rubio_gene_dds_deseq_LGP32_res_LFC_sig)
Rubio_gene_dds_deseq_LMG20012T_res_LFC_sig<- as.data.frame(Rubio_gene_dds_deseq_LMG20012T_res_LFC_sig)

nrow(Rubio_gene_dds_deseq_J2_8_res_LFC_sig)  # 33
nrow(Rubio_gene_dds_deseq_J2_9_res_LFC_sig) # 53
nrow(Rubio_gene_dds_deseq_LGP32_res_LFC_sig) # 38
nrow(Rubio_gene_dds_deseq_LMG20012T_res_LFC_sig)  # 33

### Extract list of significant Apoptosis Genes (not less than or greater than 1 LFC) using merge
Rubio_gene_dds_deseq_J2_8_res_LFC_sig_APOP <-    merge(Rubio_gene_dds_deseq_J2_8_res_LFC_sig , C_gig_rtracklayer_apop_product_final, by = "transcript_id")
Rubio_gene_dds_deseq_J2_9_res_LFC_sig_APOP <-    merge(Rubio_gene_dds_deseq_J2_9_res_LFC_sig , C_gig_rtracklayer_apop_product_final, by = "transcript_id")
Rubio_gene_dds_deseq_LGP32_res_LFC_sig_APOP <-   merge(Rubio_gene_dds_deseq_LGP32_res_LFC_sig , C_gig_rtracklayer_apop_product_final, by = "transcript_id")
Rubio_gene_dds_deseq_LMG20012T_res_LFC_sig_APOP<-merge(Rubio_gene_dds_deseq_LMG20012T_res_LFC_sig , C_gig_rtracklayer_apop_product_final, by = "transcript_id")

nrow(Rubio_gene_dds_deseq_J2_8_res_LFC_sig_APOP ) # 0
nrow(Rubio_gene_dds_deseq_J2_9_res_LFC_sig_APOP ) # 0
nrow(Rubio_gene_dds_deseq_LGP32_res_LFC_sig_APOP ) # 0
nrow(Rubio_gene_dds_deseq_LMG20012T_res_LFC_sig_APOP) # 0

#### HE ####
He_gene_counts
colnames(He_gene_counts)
He_gene_counts <- He_gene_counts[,-33]

He_coldata <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/He_coldata.csv", row.names = 1 )
View(He_coldata)  
nrow(He_coldata) 

# Make sure the columns of the count matrix and rows of the column data (sample metadata) are in the same order. 
all(rownames(He_coldata) %in% colnames(He_gene_counts ))  #Should return TRUE
# returns TRUE
all(colnames(He_gene_counts ) %in% rownames(He_coldata))  
# returns TRUE
all(rownames(He_coldata) == colnames(He_gene_counts ))  # should return TRUE
# returns true

### DATA QC PCA PLOT 
He_gene_counts_matrix <- as.matrix(He_gene_counts)
Hevstgene_counts <- varianceStabilizingTransformation(He_gene_counts_matrix, blind =TRUE)
# run PCA
pcHe_gene <- prcomp(t(Hevstgene_counts))

# Plot PCA
autoplot(pcHe_gene,
         data = He_coldata, 
         colour="Condition", 
         size=5) #  clustering by condition
autoplot(pcHe_gene,
         data = He_coldata, 
         colour="Time", 
         size=5) # less clustering by time (the transcripts more clustered by time)

## MAKE DESEQ DATA SET FROM MATRIX
## Check levels 
levels(He_coldata$Condition) # "control" "OsHV1"  
levels(He_coldata$Time) # "120hr" "12h"   "24h"   "48h"   "6h"    "Time0"
He_coldata$Time <- factor(He_coldata$Time, levels=c("Time0","6h", "12h", "24h" , "48h" ,"120hr" ))

## Creating three here so I can compare the results
He_gene_dds <- DESeqDataSetFromMatrix(countData = He_gene_counts,
                                 colData = He_coldata,
                                 design = ~ Condition + Time)  # again accounting for time in my formula, keeping condition as the last
## Prefiltering the data
He_gene_dds <- He_dds [ rowSums(counts(He_gene_dds )) > 10, ]

## DATA TRANSFORMATION AND VISUALIZATION
He_gene_dds_vst <- vst(He_gene_dds, blind = TRUE) # keep blind = true before deseq function has been run

## PCA plot visualization of individuals in the family, use vst because greater than 30 samples 
plotPCA(He_gene_dds_vst, intgroup=c("Time", "Condition")) # more variance explained some clustering of time and condition

### DIFFERENTIAL EXPRESSION ANALYSIS
He_gene_dds_deseq <- DESeq(He_gene_dds) 

## Check the resultsNames object of each to look at the available coefficients for use in lfcShrink command
resultsNames(He_gene_dds_deseq) 
#[1] "Intercept"                  "Condition_OsHV1_vs_control" "Time_6h_vs_Time0"           "Time_12h_vs_Time0"         
#[5] "Time_24h_vs_Time0"          "Time_48h_vs_Time0"          "Time_120hr_vs_Time0" 

## BUILD THE RESULTS OBJECT
He_gene_dds_res <- results(He_gene_dds_deseq, alpha=0.05, name= "Condition_OsHV1_vs_control")
He_gene_dds_res_6hr <- results(He_gene_dds_deseq, alpha=0.05, name= "Time_6h_vs_Time0" )
He_gene_dds_res_12hr <- results(He_gene_dds_deseq, alpha=0.05, name= "Time_12h_vs_Time0")
He_gene_dds_res_24hr <- results(He_gene_dds_deseq, alpha=0.05, name=  "Time_24h_vs_Time0")
He_gene_dds_res_48hr <- results(He_gene_dds_deseq, alpha=0.05, name= "Time_48h_vs_Time0")
He_gene_dds_res_120hr <- results(He_gene_dds_deseq, alpha=0.05, name= "Time_120hr_vs_Time0")

summary(He_gene_dds_res)
summary(He_gene_dds_res_6hr )
summary(He_gene_dds_res_12hr )
summary(He_gene_dds_res_24hr )
summary(He_gene_dds_res_48hr )
summary(He_gene_dds_res_120hr)

### Perform LFC Shrinkage with apeglm
He_gene_dds_res_LFC<-   lfcShrink(He_gene_dds_deseq, coef="Condition_OsHV1_vs_control", type="apeglm", res= He_gene_dds_res)
He_gene_dds_res_6hr  <- lfcShrink(He_gene_dds_deseq, coef="Time_6h_vs_Time0"   , type= "apeglm", res=He_gene_dds_res_6hr)
He_gene_dds_res_12hr <- lfcShrink(He_gene_dds_deseq, coef="Time_12h_vs_Time0"  , type= "apeglm", res=He_gene_dds_res_12hr )
He_gene_dds_res_24hr <- lfcShrink(He_gene_dds_deseq, coef= "Time_24h_vs_Time0" , type= "apeglm", res=He_gene_dds_res_24hr )
He_gene_dds_res_48hr <- lfcShrink(He_gene_dds_deseq, coef="Time_48h_vs_Time0"  , type= "apeglm", res=He_gene_dds_res_48hr )
He_gene_dds_res_120hr <-lfcShrink(He_gene_dds_deseq, coef="Time_120hr_vs_Time0", type= "apeglm", res=He_gene_dds_res_120hr)

## EXPLORATORY PLOTTING OF RESULTS 
## MA Plotting
plotMA(He_gene_dds_res_LFC, ylim = c(-5, 5))
plotMA(He_gene_dds_res_6hr  , ylim=c(-5,5))
plotMA(He_gene_dds_res_12hr , ylim=c(-5,5))
plotMA(He_gene_dds_res_24hr , ylim=c(-5,5))
plotMA(He_gene_dds_res_48hr , ylim=c(-5,5))
plotMA(He_gene_dds_res_120hr, ylim=c(-5,5))

## Histogram of P values 
# exclude genes with very small counts to avoid spikes and plot using the LFCshrinkage
hist(He_gene_dds_res_LFC$padj[He_gene_dds_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

### Subsetting Significant Genes by padj < 0.05
He_gene_dds_res_LFC_sig <-  subset(He_gene_dds_res_LFC , padj < 0.05)
He_gene_dds_res_LFC_sig $transcript_id <- row.names(He_gene_dds_res_LFC_sig )
He_gene_dds_res_LFC_sig  <- as.data.frame(He_gene_dds_res_LFC_sig )
nrow(He_gene_dds_res_LFC_sig ) # 1110

He_gene_dds_res_6hr_sig <-  subset(He_gene_dds_res_6hr, padj < 0.05)
He_gene_dds_res_12hr_sig <- subset(He_gene_dds_res_12hr, padj < 0.05)
He_gene_dds_res_24hr_sig <- subset(He_gene_dds_res_24hr, padj < 0.05)
He_gene_dds_res_48hr_sig <- subset(He_gene_dds_res_48hr, padj < 0.05)
He_gene_dds_res_120hr_sig <-subset(He_gene_dds_res_120hr, padj < 0.05)

He_gene_dds_res_6hr_sig $transcript_id <-  row.names(He_gene_dds_res_6hr_sig )
He_gene_dds_res_12hr_sig $transcript_id <- row.names(He_gene_dds_res_12hr_sig )
He_gene_dds_res_24hr_sig $transcript_id <- row.names(He_gene_dds_res_24hr_sig )
He_gene_dds_res_48hr_sig $transcript_id <- row.names(He_gene_dds_res_48hr_sig )
He_gene_dds_res_120hr_sig$transcript_id <- row.names(He_gene_dds_res_120hr_sig)

He_gene_dds_res_6hr_sig <-  as.data.frame(He_gene_dds_res_6hr_sig)
He_gene_dds_res_12hr_sig <- as.data.frame(He_gene_dds_res_12hr_sig)
He_gene_dds_res_24hr_sig <- as.data.frame(He_gene_dds_res_24hr_sig)
He_gene_dds_res_48hr_sig <- as.data.frame(He_gene_dds_res_48hr_sig)
He_gene_dds_res_120hr_sig <-as.data.frame(He_gene_dds_res_120hr_sig)

nrow(He_gene_dds_res_6hr_sig) # 246
nrow(He_gene_dds_res_12hr_sig)  # 170
nrow(He_gene_dds_res_24hr_sig) # 217
nrow(He_gene_dds_res_48hr_sig) # 209
nrow(He_gene_dds_res_120hr_sig) # 171


### Extract list of significant Apoptosis Genes (not less than or greater than 1 LFC) using merge
He_gene_dds_res_LFC_sig_APOP <- merge(He_gene_dds_res_LFC_sig , C_gig_rtracklayer_apop_product_final, by = "transcript_id")
He_gene_dds_res_LFC_sig_arranged <-arrange(He_gene_dds_res_LFC_sig_APOP, -log2FoldChange)
nrow(He_gene_dds_res_LFC_sig_arranged ) #26
View(He_gene_dds_res_LFC_sig_arranged)

He_gene_dds_res_6hr_sig_APOP <-  merge(He_gene_dds_res_6hr_sig,    C_gig_rtracklayer_apop_product_final, by = "transcript_id")
He_gene_dds_res_12hr_sig_APOP <- merge(He_gene_dds_res_12hr_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")
He_gene_dds_res_24hr_sig_APOP <- merge(He_gene_dds_res_24hr_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")
He_gene_dds_res_48hr_sig_APOP <- merge(He_gene_dds_res_48hr_sig , C_gig_rtracklayer_apop_product_final, by = "transcript_id")
He_gene_dds_res_120hr_sig_APOP <-merge(He_gene_dds_res_120hr_sig ,C_gig_rtracklayer_apop_product_final, by = "transcript_id" )

He_gene_dds_res_6hr_sig_arranged <-  arrange(He_gene_dds_res_6hr_sig_APOP , -log2FoldChange)
He_gene_dds_res_12hr_sig_arranged <- arrange(He_gene_dds_res_12hr_sig_APOP , -log2FoldChange)
He_gene_dds_res_24hr_sig_arranged <- arrange(He_gene_dds_res_24hr_sig_APOP , -log2FoldChange)
He_gene_dds_res_48hr_sig_arranged <- arrange(He_gene_dds_res_48hr_sig_APOP , -log2FoldChange)
He_gene_dds_res_120hr_sig_arranged <-arrange(He_gene_dds_res_120hr_sig_APOP, -log2FoldChange)

nrow(He_gene_dds_res_6hr_sig_arranged ) # 6
nrow(He_gene_dds_res_12hr_sig_arranged) # 5
nrow(He_gene_dds_res_24hr_sig_arranged) # 7
nrow(He_gene_dds_res_48hr_sig_arranged) # 5
nrow(He_gene_dds_res_120hr_sig_arranged) # 3

He_gene_dds_res_6hr_sig_APOP$group_by_sim <- "He_gene_dds_res_6hr_sig_APOP"
He_gene_dds_res_12hr_sig_APOP$group_by_sim <-"He_gene_dds_res_12hr_sig_APOP"
He_gene_dds_res_24hr_sig_APOP$group_by_sim <-"He_gene_dds_res_24hr_sig_APOP"
He_gene_dds_res_48hr_sig_APOP$group_by_sim <-"He_gene_dds_res_48hr_sig_APOP"
He_gene_dds_res_120hr_sig_APOP$group_by_sim<-"He_gene_dds_res_120hr_sig_APOP"

# combine data frames 
He_gene_all_sig_APOP <- rbind(He_gene_dds_res_6hr_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
                         He_gene_dds_res_12hr_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
                         He_gene_dds_res_24hr_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
                         He_gene_dds_res_48hr_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
                         He_gene_dds_res_120hr_sig_APOP[,c("product","group_by_sim","log2FoldChange")])

# Make plot
He_gene_full_LFC_plot <- ggplot(He_gene_all_sig_APOP , aes(x=product,y=log2FoldChange, fill=group_by_sim )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()

#### DELORGERIL ####
deLorgeril_gene_counts
colnames(deLorgeril_gene_counts)
deLorgeril_gene_counts <- deLorgeril_gene_counts[,-43]

deLorgeril_coldata <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/deLorgeril_coldata.csv", row.names = 1 )
View(deLorgeril_coldata)  
nrow(deLorgeril_coldata) 

# Make sure the columns of the count matrix and rows of the column data (sample metadata) are in the same order. 
all(rownames(deLorgeril_coldata) %in% colnames(deLorgeril_gene_counts))  #Should return TRUE
# returns TRUE
all(colnames(deLorgeril_gene_counts) %in% rownames(deLorgeril_coldata))  
# returns TRUE
all(rownames(deLorgeril_coldata) == colnames(deLorgeril_gene_counts))  # FALSE
# returns true

# Change order
deLorgeril_gene_counts <- deLorgeril_gene_counts[,rownames(deLorgeril_coldata)]
all(rownames(deLorgeril_coldata) %in% colnames(deLorgeril_gene_counts))  #Should return TRUE
# returns TRUE
all(colnames(deLorgeril_gene_counts) %in% rownames(deLorgeril_coldata))  
# returns TRUE
all(rownames(deLorgeril_coldata) == colnames(deLorgeril_gene_counts))  # TRUE

# split up counts and coldata into the resistant and sucsceptible families (since comparing families is not what I want)
deLorgeril_Resistant_coldata <- deLorgeril_coldata %>% subset(Condition == "AF21_Resistant" | Condition == "AF21_Resistant_control")
deLorgeril_Resistant_gene_counts <- deLorgeril_gene_counts[,row.names(deLorgeril_Resistant_coldata)]
colnames(deLorgeril_Resistant_gene_counts)

deLorgeril_Susceptible_coldata <- deLorgeril_coldata %>% subset(Condition == "AF11_Susceptible" | Condition == "AF11_Susceptible_control")
deLorgeril_Susceptible_gene_counts <- deLorgeril_gene_counts[,row.names(deLorgeril_Susceptible_coldata)]
colnames(deLorgeril_Susceptible_gene_counts)

### DATA QC PCA PLOT 
deLorgeril_Resistant_gene_counts_matrix <- as.matrix(deLorgeril_Resistant_gene_counts)
deLorgeril_Susceptible_gene_counts_matrix <- as.matrix(deLorgeril_Susceptible_gene_counts)

deLorgeril_Resistant_vstgene_counts <- vst(deLorgeril_Resistant_gene_counts_matrix, blind =TRUE)
deLorgeril_Susceptible_vstgene_counts <- vst(deLorgeril_Susceptible_gene_counts_matrix, blind =TRUE)

# run PCA
pcdeLorgeril_gene <- prcomp(t(deLorgeril_Resistant_vstgene_counts))
pcdeLorgeril_susceptible_gene <- prcomp(t(deLorgeril_Susceptible_vstgene_counts ))

# Plot PCA
autoplot(pcdeLorgeril_gene,
         data = deLorgeril_Resistant_coldata, 
         colour="Condition", 
         size=5)# all but one control cluster together on PC1 
autoplot(pcdeLorgeril_gene,
         data = deLorgeril_Resistant_coldata, 
         colour="Time", 
         size=5) # some clustering by time, 12h and 60h are similar
autoplot(pcdeLorgeril_susceptible_gene,
         data = deLorgeril_Susceptible_coldata, 
         colour="Condition", 
         size=5) # ~25% of the variance explained 
autoplot(pcdeLorgeril_susceptible_gene,
         data = deLorgeril_Susceptible_coldata, 
         colour="Time", 
         size=5) # more clustering by time here
## MAKE DESEQ DATA SET FROM MATRIX
## Check levels 
levels(deLorgeril_Resistant_coldata$Condition)# "AF11_Susceptible"         "AF11_Susceptible_control" "AF21_Resistant"           "AF21_Resistant_control" 
levels(deLorgeril_Resistant_coldata$Condition) #  "AF21_Resistant"         "AF21_Resistant_control"       
levels(deLorgeril_Susceptible_coldata$Condition)  
levels(deLorgeril_Susceptible_coldata$Condition)   #    "AF11_Susceptible"         "AF11_Susceptible_control"    
levels(deLorgeril_Resistant_coldata$Time) # "0h"  "12h" "24h" "48h" "60h" "6h"  "72h"
levels(deLorgeril_Susceptible_coldata$Time) # "0h"  "12h" "24h" "48h" "60h" "6h"  "72h"

## Creating three here so I can compare the results
deLorgeril_Resistant_gene_dds <- DESeqDataSetFromMatrix(countData = deLorgeril_Resistant_gene_counts,
                                                   colData = deLorgeril_Resistant_coldata,
                                                   design = ~ Time) 
deLorgeril_Susceptible_gene_dds <- DESeqDataSetFromMatrix(countData = deLorgeril_Susceptible_gene_counts,
                                                     colData = deLorgeril_Susceptible_coldata,
                                                     design = ~ Time)  

## Prefiltering the data
deLorgeril_Resistant_gene_dds <- deLorgeril_Resistant_dds[ rowSums(counts(deLorgeril_Resistant_gene_dds)) > 10, ]
deLorgeril_Susceptible_gene_dds <- deLorgeril_Susceptible_dds[ rowSums(counts(deLorgeril_Susceptible_gene_dds)) > 10, ]

## DATA TRANSFORMATION AND VISUALIZATION
deLorgeril_Resistant_gene_dds_vst <- vst(deLorgeril_Resistant_gene_dds, blind = TRUE) # keep blind = true before deseq function has been run
deLorgeril_Susceptible_gene_dds_vst <- vst(deLorgeril_Susceptible_gene_dds , blind = TRUE)

## PCA plot visualization of individuals in the family, use vst because greater than 30 samples 
plotPCA(deLorgeril_Resistant_gene_dds_vst , intgroup=c("Time", "Condition")) # less variation explained than before
plotPCA(deLorgeril_Susceptible_gene_dds_vst , intgroup=c("Time", "Condition")) # close clustering of 12, 24h, 60 hr

### DIFFERENTIAL EXPRESSION ANALYSIS
deLorgeril_Resistant_gene_dds_deseq <- DESeq(deLorgeril_Resistant_gene_dds) 
deLorgeril_Susceptible_gene_dds_deseq <- DESeq(deLorgeril_Susceptible_gene_dds) 

## Check the resultsNames object of each to look at the available coefficients for use in lfcShrink command
resultsNames(deLorgeril_Resistant_gene_dds_deseq ) #"Intercept""Time_6h_vs_0h"  "Time_12h_vs_0h" "Time_24h_vs_0h" "Time_48h_vs_0h" "Time_60h_vs_0h" "Time_72h_vs_0h"
resultsNames(deLorgeril_Susceptible_gene_dds_deseq ) #"Intercept"      "Time_6h_vs_0h"  "Time_12h_vs_0h" "Time_24h_vs_0h" "Time_48h_vs_0h" "Time_60h_vs_0h" "Time_72h_vs_0h"

## BUILD THE RESULTS OBJECT
deLorgeril_Resistant_gene_dds_res_12 <- results(deLorgeril_Resistant_gene_dds_deseq , alpha=0.05, name= "Time_12h_vs_0h" )
deLorgeril_Resistant_gene_dds_res_24 <- results(deLorgeril_Resistant_gene_dds_deseq , alpha=0.05, name= "Time_24h_vs_0h" )
deLorgeril_Resistant_gene_dds_res_48 <- results(deLorgeril_Resistant_gene_dds_deseq , alpha=0.05, name= "Time_48h_vs_0h" )
deLorgeril_Resistant_gene_dds_res_60 <- results(deLorgeril_Resistant_gene_dds_deseq , alpha=0.05, name= "Time_60h_vs_0h" )
deLorgeril_Resistant_gene_dds_res_72 <- results(deLorgeril_Resistant_gene_dds_deseq , alpha=0.05, name= "Time_72h_vs_0h" )

mcols(deLorgeril_Susceptible_dds_deseq)
deLorgeril_Susceptible_dds_res_48 <- results(deLorgeril_Susceptible_dds_deseq , alpha=0.05, name= "Time_48h_vs_0h" )
deLorgeril_Susceptible_dds_res_60 <- results(deLorgeril_Susceptible_dds_deseq , alpha=0.05, name= "Time_60h_vs_0h" )
deLorgeril_Susceptible_dds_res_72 <- results(deLorgeril_Susceptible_dds_deseq , alpha=0.05, name= "Time_72h_vs_0h" )

### Perform LFC Shrinkage with apeglm
## NOTES 
# Before plotting we need to apply an LFC shrinkage, set the coef as the specific comparison in the ResultsNames function of the deseq object
# Issue that the specific comparisons I set in my results formulas are not available in the ResultsNames coef list.
# More detailed notes about LFC Shrinkage are in the code for the Zhang Vibrio

## DECISION: USE SAME RES OBJECT TO KEEP ALPHA ADJUSTMENT, and use LFCShrink apeglm
deLorgeril_Resistant_dds_res_48_LFC <- lfcShrink(deLorgeril_Resistant_dds_deseq , coef="Time_48h_vs_0h", type="apeglm",res=deLorgeril_Resistant_dds_res_48)
deLorgeril_Resistant_dds_res_60_LFC <- lfcShrink(deLorgeril_Resistant_dds_deseq , coef="Time_60h_vs_0h", type="apeglm",res=deLorgeril_Resistant_dds_res_60)
deLorgeril_Resistant_dds_res_72_LFC <- lfcShrink(deLorgeril_Resistant_dds_deseq , coef="Time_72h_vs_0h", type="apeglm",res=deLorgeril_Resistant_dds_res_72)
deLorgeril_Susceptible_dds_res_48_LFC <- lfcShrink(deLorgeril_Susceptible_dds_deseq , coef= "Time_48h_vs_0h", type="apeglm", res=deLorgeril_Susceptible_dds_res_48)
deLorgeril_Susceptible_dds_res_60_LFC <- lfcShrink(deLorgeril_Susceptible_dds_deseq , coef= "Time_60h_vs_0h", type="apeglm", res=deLorgeril_Susceptible_dds_res_60)
deLorgeril_Susceptible_dds_res_72_LFC <- lfcShrink(deLorgeril_Susceptible_dds_deseq , coef= "Time_72h_vs_0h", type="apeglm", res=deLorgeril_Susceptible_dds_res_72)

## EXPLORATORY PLOTTING OF RESULTS 
## MA Plotting
plotMA(deLorgeril_Resistant_dds_res_48_LFC , ylim = c(-5, 5))
plotMA(deLorgeril_Resistant_dds_res_60_LFC , ylim = c(-5, 5)) # 60hr has more significant genes
plotMA(deLorgeril_Resistant_dds_res_72_LFC , ylim = c(-5, 5)) # 72 hrs has less
plotMA(deLorgeril_Susceptible_dds_res_48_LFC, ylim = c(-5, 5))
plotMA(deLorgeril_Susceptible_dds_res_60_LFC, ylim = c(-5, 5)) # 60hr has many more significant genes
plotMA(deLorgeril_Susceptible_dds_res_72_LFC, ylim = c(-5, 5)) # 72hr has less

## Histogram of P values 
# exclude genes with very small counts to avoid spikes and plot using the LFCshrinkage
hist(deLorgeril_Resistant_dds_res_48_LFC $padj[deLorgeril_Resistant_dds_res_48_LFC$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white")
hist(deLorgeril_Resistant_dds_res_60_LFC $padj[deLorgeril_Resistant_dds_res_60_LFC$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white")
hist(deLorgeril_Resistant_dds_res_72_LFC $padj[deLorgeril_Resistant_dds_res_72_LFC$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white")
hist(deLorgeril_Susceptible_dds_res_48_LFC$padj[deLorgeril_Susceptible_dds_res_48_LFC$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white")
hist(deLorgeril_Susceptible_dds_res_60_LFC$padj[deLorgeril_Susceptible_dds_res_60_LFC$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white")
hist(deLorgeril_Susceptible_dds_res_72_LFC$padj[deLorgeril_Susceptible_dds_res_72_LFC$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white")

### Subsetting Significant Genes by padj < 0.05
# again, only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
deLorgeril_Resistant_dds_res_48_LFC_sig <- subset(  deLorgeril_Resistant_dds_res_48_LFC, padj < 0.05)
deLorgeril_Resistant_dds_res_60_LFC_sig <- subset(  deLorgeril_Resistant_dds_res_60_LFC, padj < 0.05)
deLorgeril_Resistant_dds_res_72_LFC_sig <- subset(  deLorgeril_Resistant_dds_res_72_LFC, padj < 0.05)
deLorgeril_Susceptible_dds_res_48_LFC_sig <- subset(deLorgeril_Susceptible_dds_res_48_LFC, padj < 0.05)
deLorgeril_Susceptible_dds_res_60_LFC_sig <- subset(deLorgeril_Susceptible_dds_res_60_LFC, padj < 0.05)
deLorgeril_Susceptible_dds_res_72_LFC_sig <- subset(deLorgeril_Susceptible_dds_res_72_LFC, padj < 0.05)

deLorgeril_Resistant_dds_res_48_LFC_sig$transcript_id <- row.names(  deLorgeril_Resistant_dds_res_48_LFC_sig) 
deLorgeril_Resistant_dds_res_60_LFC_sig$transcript_id <- row.names(  deLorgeril_Resistant_dds_res_60_LFC_sig) 
deLorgeril_Resistant_dds_res_72_LFC_sig$transcript_id <- row.names(  deLorgeril_Resistant_dds_res_72_LFC_sig) 
deLorgeril_Susceptible_dds_res_48_LFC_sig$transcript_id <- row.names(deLorgeril_Susceptible_dds_res_48_LFC_sig) 
deLorgeril_Susceptible_dds_res_60_LFC_sig$transcript_id <- row.names(deLorgeril_Susceptible_dds_res_60_LFC_sig) 
deLorgeril_Susceptible_dds_res_72_LFC_sig$transcript_id <- row.names(deLorgeril_Susceptible_dds_res_72_LFC_sig) 

deLorgeril_Resistant_dds_res_48_LFC_sig   <-as.data.frame(deLorgeril_Resistant_dds_res_48_LFC_sig)
deLorgeril_Resistant_dds_res_60_LFC_sig   <-as.data.frame(deLorgeril_Resistant_dds_res_60_LFC_sig)
deLorgeril_Resistant_dds_res_72_LFC_sig   <-as.data.frame(deLorgeril_Resistant_dds_res_72_LFC_sig)
deLorgeril_Susceptible_dds_res_48_LFC_sig<-as.data.frame(deLorgeril_Susceptible_dds_res_48_LFC_sig)
deLorgeril_Susceptible_dds_res_60_LFC_sig<-as.data.frame(deLorgeril_Susceptible_dds_res_60_LFC_sig)
deLorgeril_Susceptible_dds_res_72_LFC_sig<-as.data.frame(deLorgeril_Susceptible_dds_res_72_LFC_sig)

nrow(deLorgeril_Resistant_dds_res_48_LFC_sig) #1593
nrow(deLorgeril_Resistant_dds_res_60_LFC_sig) # 3403
nrow(deLorgeril_Resistant_dds_res_72_LFC_sig) # 2309
nrow(deLorgeril_Susceptible_dds_res_48_LFC_sig) # 1778
nrow(deLorgeril_Susceptible_dds_res_60_LFC_sig) # 10425
nrow(deLorgeril_Susceptible_dds_res_72_LFC_sig) # 2991

### GENE CLUSTERING ANALYSIS HEATMAPS  
# Extract genes with the highest variance across samples for each comparison using either vst or rlog transformed data
# This heatmap rather than plotting absolute expression strength plot the amount by which each gene deviates in a specific sample from the gene’s average across all samples. 
# example codes from RNAseq workflow: https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#other-comparisons
deLorgeril_Resistant_dds_vst_assay  <-  head(order(rowVars(assay(deLorgeril_Resistant_dds_vst  )), decreasing = TRUE), 200)
deLorgeril_res_mat <- assay(deLorgeril_Resistant_dds_vst )[deLorgeril_Resistant_dds_vst_assay ,]
deLorgeril_res_mat <- deLorgeril_res_mat - rowMeans(deLorgeril_res_mat)
deLorgeril_res_anno <- as.data.frame(colData(deLorgeril_Resistant_dds_vst )[, c("Condition","Time")])
deLorgeril_res_heatmap <- pheatmap(deLorgeril_res_mat, annotation_col = deLorgeril_res_anno)
head(deLorgeril_res_mat) # largely clustering by time

deLorgeril_Susceptible_dds_vst_assay  <-  head(order(rowVars(assay(deLorgeril_Susceptible_dds_vst  )), decreasing = TRUE), 200)
deLorgeril_mat <- assay(deLorgeril_Susceptible_dds_vst )[deLorgeril_Resistant_dds_vst_assay ,]
deLorgeril_mat <- deLorgeril_mat - rowMeans(deLorgeril_mat)
deLorgeril_anno <- as.data.frame(colData(deLorgeril_Susceptible_dds_vst )[, c("Condition","Time")])
deLorgeril_heatmap <- pheatmap(deLorgeril_mat, annotation_col = deLorgeril_anno)
head(deLorgeril_mat)  # largely clustering by time

# Gene clustering heatmap with only apoptosis genes #
# vector C_gig_apop transcript IDs
C_gig_rtracklayer_apop_product_final_transcript_id 
# Search original counts table for apoptosis genes and do vst on just these
deLorgeril_Resistant_counts_apop <- deLorgeril_Resistant_counts[row.names(deLorgeril_Resistant_counts) %in% C_gig_rtracklayer_apop_product_final_transcript_id,]
deLorgeril_Susceptible_counts_apop <- deLorgeril_Susceptible_counts[row.names(deLorgeril_Susceptible_counts) %in% C_gig_rtracklayer_apop_product_final_transcript_id,]
nrow(deLorgeril_Resistant_counts_apop ) #659
nrow(deLorgeril_Susceptible_counts_apop ) #659

deLorgeril_Resistant_counts_apop_dds <- DESeqDataSetFromMatrix(countData = deLorgeril_Resistant_counts_apop,
                                                               colData = deLorgeril_Resistant_coldata,
                                                               design = ~ Time) # using same formula as before
deLorgeril_Susceptible_counts_apop_dds <- DESeqDataSetFromMatrix(countData = deLorgeril_Susceptible_counts_apop,
                                                                 colData = deLorgeril_Susceptible_coldata,
                                                                 design = ~ Time) # using same formula as before

# Prefiltering the data and running vst
deLorgeril_Resistant_counts_apop_dds <- deLorgeril_Resistant_counts_apop_dds[ rowSums(counts(deLorgeril_Resistant_counts_apop_dds)) > 10, ]
deLorgeril_Susceptible_counts_apop_dds <- deLorgeril_Susceptible_counts_apop_dds[ rowSums(counts(deLorgeril_Susceptible_counts_apop_dds)) > 10, ]
deLorgeril_Resistant_counts_apop_dds_vst <- varianceStabilizingTransformation(deLorgeril_Resistant_counts_apop_dds, blind=TRUE)
deLorgeril_Susceptible_counts_apop_dds_vst <- varianceStabilizingTransformation(deLorgeril_Susceptible_counts_apop_dds, blind=TRUE)

## PCA plot of rlog transformed counts for apoptosis
plotPCA(deLorgeril_Resistant_counts_apop_dds_vst , intgroup="Time") # 0hr are outliers in apoptosis gene expression
# all timepoints other than time 0 cluster pretty closely together with not much separation
plotPCA(deLorgeril_Susceptible_counts_apop_dds_vst, intgroup="Time") # greater separation of timepoints

# heatmap of apoptosis genes 
deLorgeril_Resistant_counts_apop_dds_vst_assay <-  assay(deLorgeril_Resistant_counts_apop_dds_vst)[,]
deLorgeril_Resistant_apop_assay_mat <- deLorgeril_Resistant_counts_apop_dds_vst_assay  - rowMeans(deLorgeril_Resistant_counts_apop_dds_vst_assay )
deLorgeril_Resistant_apop_assay_anno <- as.data.frame(colData(deLorgeril_Resistant_counts_apop_dds_vst)[, c("Condition","Time")])
deLorgeril_Resistant_apop_assay_heatmap <- pheatmap(deLorgeril_Resistant_apop_assay_mat  , annotation_col = deLorgeril_Resistant_apop_assay_anno)
head(deLorgeril_Resistant_apop_assay_mat ) 

deLorgeril_Susceptible_counts_apop_dds_vst_assay <-  assay(deLorgeril_Susceptible_counts_apop_dds_vst)[,]
deLorgeril_Susceptible_apop_assay_mat <- deLorgeril_Susceptible_counts_apop_dds_vst_assay  - rowMeans(deLorgeril_Susceptible_counts_apop_dds_vst_assay )
deLorgeril_Susceptible_apop_assay_anno <- as.data.frame(colData(deLorgeril_Susceptible_counts_apop_dds_vst)[, c("Condition","Time")])
deLorgeril_Susceptible_apop_assay_heatmap <- pheatmap(deLorgeril_Susceptible_apop_assay_mat  , annotation_col = deLorgeril_Susceptible_apop_assay_anno)
head(deLorgeril_Susceptible_apop_assay_mat ) 

# heatmap of most variable apoptosis genes (this selects genes with the greatest variance in the sample)
topVarGenes_Resistant_counts_apop_assay <-  head(order(rowVars(assay(deLorgeril_Resistant_counts_apop_dds_vst )), decreasing = TRUE), 100) 
top_Var_Resistant_counts_apop_assay_mat<- assay(deLorgeril_Resistant_counts_apop_dds_vst)[topVarGenes_Resistant_counts_apop_assay,]
top_Var_Resistant_counts_apop_assay_mat <- top_Var_Resistant_counts_apop_assay_mat - rowMeans(top_Var_Resistant_counts_apop_assay_mat)
top_Var_Resistant_counts_apop_assay_anno <- as.data.frame(colData(deLorgeril_Resistant_counts_apop_dds_vst)[, c("Condition","Time")])
top_Var_Resistant_counts_apop_assay_heatmap <- pheatmap(top_Var_Resistant_counts_apop_assay_mat  , annotation_col = top_Var_Resistant_counts_apop_assay_anno)
head(top_Var_Resistant_counts_apop_assay_mat )

topVarGenes_Susceptible_counts_apop_assay <-  head(order(rowVars(assay(deLorgeril_Susceptible_counts_apop_dds_vst )), decreasing = TRUE), 100) 
top_Var_Susceptible_counts_apop_assay_mat<- assay(deLorgeril_Susceptible_counts_apop_dds_vst)[topVarGenes_Susceptible_counts_apop_assay,]
top_Var_Susceptible_counts_apop_assay_mat <- top_Var_Susceptible_counts_apop_assay_mat - rowMeans(top_Var_Susceptible_counts_apop_assay_mat)
top_Var_Susceptible_counts_apop_assay_anno <- as.data.frame(colData(deLorgeril_Susceptible_counts_apop_dds_vst)[, c("Condition","Time")])
top_Var_Susceptible_counts_apop_assay_heatmap <- pheatmap(top_Var_Susceptible_counts_apop_assay_mat  , annotation_col = top_Var_Susceptible_counts_apop_assay_anno)
head(top_Var_Susceptible_counts_apop_assay_mat )

# annotate the top 100 most variable genes  
# reorder annotation table to match ordering in heatmap 
top_Var_Susceptible_counts_apop_assay_heatmap_reorder <-rownames(top_Var_Susceptible_counts_apop_assay_mat[top_Var_Susceptible_counts_apop_assay_heatmap$tree_row[["order"]],])
# annotate the row.names
top_Var_Susceptible_counts_apop_assay_prot <- as.data.frame(top_Var_Susceptible_counts_apop_assay_heatmap_reorder)
colnames(top_Var_Susceptible_counts_apop_assay_prot)[1] <- "transcript_id"
top_Var_Susceptible_counts_apop_assay_prot_annot <- left_join(top_Var_Susceptible_counts_apop_assay_prot, select(C_gig_rtracklayer_transcripts, transcript_id, product, gene), by = "transcript_id")

top_Var_Resistant_counts_apop_assay_heatmap_reorder <-rownames(top_Var_Resistant_counts_apop_assay_mat[top_Var_Resistant_counts_apop_assay_heatmap$tree_row[["order"]],])
# annotate the row.names
top_Var_Resistant_counts_apop_assay_prot <- as.data.frame(top_Var_Resistant_counts_apop_assay_heatmap_reorder)
colnames(top_Var_Resistant_counts_apop_assay_prot)[1] <- "transcript_id"
top_Var_Resistant_counts_apop_assay_prot_annot <- left_join(top_Var_Resistant_counts_apop_assay_prot, select(C_gig_rtracklayer_transcripts, transcript_id, product, gene), by = "transcript_id")

### Extract list of significant Apoptosis Genes (not less than or greater than 1 LFC) using merge
deLorgeril_Resistant_dds_res_48_LFC_sig_APOP <- merge(deLorgeril_Resistant_dds_res_48_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id" )
deLorgeril_Resistant_dds_res_60_LFC_sig_APOP <- merge(deLorgeril_Resistant_dds_res_60_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id" )
deLorgeril_Resistant_dds_res_72_LFC_sig_APOP <- merge(deLorgeril_Resistant_dds_res_72_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id" )
deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP <- merge(deLorgeril_Susceptible_dds_res_48_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")
deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP <- merge(deLorgeril_Susceptible_dds_res_60_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")
deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP <- merge(deLorgeril_Susceptible_dds_res_72_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")

deLorgeril_Resistant_dds_res_48_LFC_sig_APOP  <- arrange(deLorgeril_Resistant_dds_res_48_LFC_sig_APOP , -log2FoldChange)
deLorgeril_Resistant_dds_res_60_LFC_sig_APOP  <- arrange(deLorgeril_Resistant_dds_res_60_LFC_sig_APOP , -log2FoldChange)
deLorgeril_Resistant_dds_res_72_LFC_sig_APOP  <- arrange(deLorgeril_Resistant_dds_res_72_LFC_sig_APOP , -log2FoldChange)
deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP <- arrange(deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP, -log2FoldChange)
deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP <- arrange(deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP, -log2FoldChange)
deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP <- arrange(deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP, -log2FoldChange)

nrow(deLorgeril_Resistant_dds_res_48_LFC_sig_APOP) #25
nrow(deLorgeril_Resistant_dds_res_60_LFC_sig_APOP) # 54
nrow(deLorgeril_Resistant_dds_res_72_LFC_sig_APOP) # 33
nrow(deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP) #34
nrow(deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP) # 186
nrow(deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP) # 47

View(deLorgeril_Resistant_dds_res_48_LFC_sig_APOP  $product)
View(deLorgeril_Resistant_dds_res_60_LFC_sig_APOP  $product)
View(deLorgeril_Resistant_dds_res_72_LFC_sig_APOP  $product)
View(deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP$product)
View(deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP$product) 
View(deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP$product)

# Compare apoptosis genes between group_by_sim groups
deLorgeril_Resistant_dds_res_48_LFC_sig_APOP  $group_by_sim <- "Resistant_dds_res_48"
deLorgeril_Resistant_dds_res_60_LFC_sig_APOP  $group_by_sim <- "Resistant_dds_res_60"
deLorgeril_Resistant_dds_res_72_LFC_sig_APOP  $group_by_sim <- "Resistant_dds_res_72"
deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP$group_by_sim <- "Susceptible_dds_res_48"
deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP$group_by_sim <- "Susceptible_dds_res_60"
deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP$group_by_sim <- "Susceptible_dds_res_72"

# combine data frames 
deLorgeril_all_sig_APOP <- rbind( 
  deLorgeril_Resistant_dds_res_48_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
  deLorgeril_Resistant_dds_res_60_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
  deLorgeril_Resistant_dds_res_72_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
  deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
  deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
  deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")])

# Make plot or up and downregulated
deLorgeril_all_sig_APOP_downregulated <- deLorgeril_all_sig_APOP %>% filter(log2FoldChange <= 0)
deLorgeril_all_sig_APOP_upregulated <- deLorgeril_all_sig_APOP %>% filter(log2FoldChange > 0)

deLorgeril_all_sig_APOP_downregulated_plot <- ggplot(deLorgeril_all_sig_APOP_downregulated , aes(x=product,y=log2FoldChange, fill=group )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()
deLorgeril_all_sig_APOP_upregulated_plot <- ggplot(deLorgeril_all_sig_APOP_upregulated, aes(x=product,y=log2FoldChange, fill=group )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()

deLorgeril_Resistant_dds_res_48_LFC_sig_APOP <- ggplot(deLorgeril_Resistant_dds_res_48_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Resistant 48hr vs Control") +
  ylab("Log2 Fold Change")
deLorgeril_Resistant_dds_res_60_LFC_sig_APOP   <- ggplot(deLorgeril_Resistant_dds_res_60_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Resistant 60hr vs Control") +
  ylab("Log2 Fold Change") # mostly upregulation
deLorgeril_Resistant_dds_res_72_LFC_sig_APOP   <- ggplot(deLorgeril_Resistant_dds_res_72_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Resistant 72hr vs Control") +
  ylab("Log2 Fold Change") # downregulation of cathepsins

deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP <- ggplot(deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Susceptible 48hr vs Control") +
  ylab("Log2 Fold Change")
deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP <- ggplot(deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Susceptible 60hr vs Control") +
  ylab("Log2 Fold Change") # TLR downregulation, a lot of upregulation
deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP <- ggplot(deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Resistant 72hr vs Control") +
  ylab("Log2 Fold Change") 


Dermo_gene_counts
Probiotic_gene_counts
ROD_gene_counts
Pro_RE22_gene_counts






#### PCA AND HEATMAPS OF ORTHOLOGOUS GENE COUNTS ####


## Remove Batch effects from experiment for C_vir
plotPCA(C_vir_full_counts_vst, "Experiment") # grouping by experiment, ROD and probiotic cluster more closely
mat_C_vir <- assay(C_vir_full_counts_vst)
mat_C_vir <- limma::removeBatchEffect(mat_C_vir, C_vir_full_counts_vst$Experiment)
#Coefficients not estimable: batch1 batch3 batch5 batch6 
#Warning message:
#  Partial NA coefficients for 67876 probe(s) 
assay(C_vir_full_counts_vst) <- mat_C_vir
plotPCA(C_vir_full_counts_vst, "Experiment") # Probiotic and ROD now cluster together
plotPCA(C_vir_full_counts_vst, "Sample")
plotPCA(C_vir_full_counts_vst, "Time") # no clustering by time
plotPCA(C_vir_full_counts_vst, "Family")

plotPCA(C_gig_full_counts_vst, "Experiment") # grouping by experiment, Rubio delorgeril cluster, HE and Zhang far apart
mat_C_gig <- assay(C_gig_full_counts_vst)
mat_C_gig <- limma::removeBatchEffect(mat_C_gig, C_gig_full_counts_vst$Experiment)
#Coefficients not estimable: batch4 batch5 batch6 
#Warning message:
#  Partial NA coefficients for 86859 probe(s) 
assay(C_gig_full_counts_vst) <- mat_C_gig
plotPCA(C_gig_full_counts_vst, "Experiment") # He, Rubio and Zhang cluster closely 
plotPCA(C_gig_full_counts_vst, "Sample")
plotPCA(C_gig_full_counts_vst, "Time") # no clustering by time
plotPCA(C_gig_full_counts_vst, "Family")

# subset PCA for apoptosis genes
C_vir_full_counts_apop <- C_vir_full_counts[row.names(C_vir_full_counts) %in% C_vir_rtracklayer_apop_product_final_ID,]
nrow(C_vir_full_counts_apop) # 1026
all(row.names(C_vir_full_counts_apop) %in% C_vir_rtracklayer_apop_product_final_ID) # TRUE
all(C_vir_rtracklayer_apop_product_final_ID %in% row.names(C_vir_full_counts_apop)) # TRUE

C_vir_full_counts_apop_list <- row.names(C_vir_full_counts_apop)
plotPCA(C_vir_full_counts_vst[C_vir_full_counts_apop_list,], "Experiment") # similar clustering where ROD and Probiotic cluster

C_gig_full_counts_apop <- C_gig_full_counts[row.names(C_gig_full_counts) %in% C_gig_rtracklayer_apop_product_final_transcript_id,]
C_gig_full_counts_apop_list <- row.names(C_gig_full_counts_apop)
plotPCA(C_gig_full_counts_vst[C_gig_full_counts_apop_list,], "Experiment") # some more overlap between He and deLorgeril

### Heatmaps of apoptosis gene vsts ###
## C_virginica
# heatmap of all apoptosis genes 
C_vir_full_counts_apop_assay <-  assay(C_vir_full_counts_vst)[C_vir_full_counts_apop_list,]
C_vir_full_counts_apop_assay_mat <- C_vir_full_counts_apop_assay - rowMeans(C_vir_full_counts_apop_assay)
C_vir_full_counts_apop_assay_anno <- as.data.frame(colData(C_vir_full_counts_vst)[, c("Condition","Experiment")])
C_vir_full_counts_apop_assay_heatmap <- pheatmap(C_vir_full_counts_apop_assay_mat  , annotation_col = C_vir_full_counts_apop_assay_anno)
# Probiotic and ROD cluster closely, CGX and GX cluster more with Dermo DA (susceptible), while F3L and LB cluster more closely
# the sucsceptible to Dermo are more tolerant to ROD (GX), and vice versa..interesting
# The ROD Tolerant (GX and CGX), Dermo Susceptible DA and Probiotic are all in the same larger cluster
# the ROD susceptible (F3L) and the Dermo tolerant LB are in the same larger cluster

# heatmap of most variable apoptosis genes for both C_vir and C_gug (this selects genes with the greatest variance in the sample)
topVarGenes_C_vir_full_counts_apop_assay <-  head(order(rowVars(C_vir_full_counts_apop_assay), decreasing = TRUE), 100) 
top_Var_C_vir_full_counts_apop_assay_mat <- C_vir_full_counts_apop_assay[topVarGenes_C_vir_full_counts_apop_assay,]
top_Var_C_vir_full_counts_apop_assay_mat <- top_Var_C_vir_full_counts_apop_assay_mat - rowMeans(top_Var_C_vir_full_counts_apop_assay_mat)
top_Var_C_vir_full_counts_apop_assay_anno <- as.data.frame(colData(C_vir_full_counts_vst)[, c("Experiment","Condition")])
top_Var_C_vir_full_counts_apop_assay_heatmap <- pheatmap(top_Var_C_vir_full_counts_apop_assay_mat  , annotation_col = top_Var_C_vir_full_counts_apop_assay_anno)
head(top_Var_C_vir_full_counts_apop_assay_mat ) 
# Same relationship as above, but tighter clustering

# reorder annotation table to match ordering in heatmap 
top_Var_C_vir_full_counts_apop_assay_heatmap_reorder <- rownames(top_Var_C_vir_full_counts_apop_assay_mat[top_Var_C_vir_full_counts_apop_assay_heatmap$tree_row[["order"]],])
class(top_Var_C_vir_full_counts_apop_assay_heatmap_reorder)
# annotate the row.names
top_Var_C_virginica_apop_assay_mat_prot <- as.data.frame(top_Var_C_vir_full_counts_apop_assay_heatmap_reorder )
class(top_Var_C_virginica_apop_assay_mat_prot)
colnames(top_Var_C_virginica_apop_assay_mat_prot)[1] <- "ID"
class(top_Var_C_virginica_apop_assay_mat_prot$ID)
top_Var_C_virginica_apop_assay_mat_prot$ID <- as.character(top_Var_C_virginica_apop_assay_mat_prot$ID)
top_Var_C_virginica_apop_assay_mat_prot_annot <- left_join(top_Var_C_virginica_apop_assay_mat_prot, select(C_vir_rtracklayer, ID, product), by = "ID")

## C_gigas genes 
C_gig_full_counts_apop_assay <-  assay(C_gig_full_counts_vst)[C_gig_full_counts_apop_list,]
C_gig_full_counts_apop_assay_mat <- C_gig_full_counts_apop_assay - rowMeans(C_gig_full_counts_apop_assay)
C_gig_full_counts_apop_assay_anno <- as.data.frame(colData(C_gig_full_counts_vst)[, c("Experiment","Condition")])
C_gig_full_counts_apop_assay_heatmap <- pheatmap(C_gig_full_counts_apop_assay_mat  , annotation_col = C_gig_full_counts_apop_assay_anno)
# clustering of delorgeril and HE: deLorgeril AF11 susceptible, HE OsHV1 challenge, 
# delorgeril OsHV1 resistant are in their own separate cluster
# Zhang and Rubio cluster very closely together
# some delorgeril and rubio clustering

# heatmap of most variable apoptosis genes for both C_vir and C_gug (this selects genes with the greatest variance in the sample)
topVarGenes_C_gig_full_counts_apop_assay <-  head(order(rowVars(C_gig_full_counts_apop_assay), decreasing = TRUE), 100) 
top_Var_C_gig_full_counts_apop_assay_mat <- C_gig_full_counts_apop_assay[topVarGenes_C_gig_full_counts_apop_assay,]
top_Var_C_gig_full_counts_apop_assay_mat <- top_Var_C_gig_full_counts_apop_assay_mat - rowMeans(top_Var_C_gig_full_counts_apop_assay_mat)
top_Var_C_gig_full_counts_apop_assay_anno <- as.data.frame(colData(C_gig_full_counts_vst)[, c("Experiment","Condition")])
top_Var_C_gig_full_counts_apop_assay_heatmap <- pheatmap(top_Var_C_gig_full_counts_apop_assay_mat  , annotation_col = top_Var_C_gig_full_counts_apop_assay_anno)
head(top_Var_C_gig_full_counts_apop_assay_mat ) 
# Same relationship as above, but tighter clustering

# reorder annotation table to match ordering in heatmap 
top_Var_C_gig_full_counts_apop_assay_heatmap_reorder <- rownames(top_Var_C_gig_full_counts_apop_assay_mat[top_Var_C_gig_full_counts_apop_assay_heatmap$tree_row[["order"]],])
class(top_Var_C_gig_full_counts_apop_assay_heatmap_reorder)
# annotate the row.names
top_Var_C_gig_apop_assay_mat_prot <- as.data.frame(top_Var_C_gig_full_counts_apop_assay_heatmap_reorder )
class(top_Var_C_gig_apop_assay_mat_prot)
colnames(top_Var_C_gig_apop_assay_mat_prot)[1] <- "transcript_id"
class(top_Var_C_gig_apop_assay_mat_prot$transcript_id)
top_Var_C_gig_apop_assay_mat_prot$transcript_id <- as.character(top_Var_C_gig_apop_assay_mat_prot$transcript_id)
top_Var_C_gig_assay_mat_prot_annot <- left_join(top_Var_C_gig_apop_assay_mat_prot, select(C_gig_rtracklayer_apop_product_final, transcript_id, product), by = "transcript_id")

#isolate interesting clusters
#six_hr_comparison_cluster <- c("XM_022455505.1", "XM_022484575.1", "XM_022461506.1", "XM_022430618.1", "XM_022490512.1",
#                               "XM_022461508.1", "XM_022464459.1", "XM_022483469.1", "XM_022483473.1", "XM_022442223.1", "XM_022457463.1", "XM_022442224.1")
#six_hr_comparison_cluster <- as.data.frame(six_hr_comparison_cluster)
#colnames(six_hr_comparison_cluster)[1] <- "transcript_id"
#six_hr_comparison_cluster_subset <- subset(Res_mat_6hr_prot_annot, transcript_id %in% six_hr_comparison_cluster$transcript_id)


## Extract apop genes from each and match to names ###
# For C_vir the gene# is the "Parent" column  vector C_vir_apop transcript Parent
C_vir_rtracklayer_apop_product_final_parent <- unique(C_vir_rtracklayer_apop_product_final$Parent)
View(C_vir_rtracklayer_apop_product_final_parent)

# C_vir: Search counts tables for apoptosis genes using the Parent column from genome
Dermo_gene_counts_apop <- Dermo_gene_counts[row.names(Dermo_gene_counts) %in% C_vir_rtracklayer_apop_product_final_parent,]
nrow(Dermo_gene_counts_apop) # 45
View(Dermo_gene_counts_apop)
Probiotic_gene_counts_apop <- Probiotic_gene_counts[row.names(Probiotic_gene_counts) %in% C_vir_rtracklayer_apop_product_final_parent,]
nrow(Probiotic_gene_counts_apop) # 198
head(Probiotic_gene_counts_apop)
ROD_gene_counts_apop <- ROD_gene_counts[row.names(ROD_gene_counts) %in% C_vir_rtracklayer_apop_product_final_parent,]
nrow(ROD_gene_counts_apop) # 98
head(ROD_gene_counts_apop)
# Make the rownames a column for each so I can merge on it 
Dermo_gene_counts_apop$Parent <- row.names(Dermo_gene_counts_apop)
Probiotic_gene_counts_apop$Parent <- row.names(Probiotic_gene_counts_apop)
ROD_gene_counts_apop$Parent <- row.names(ROD_gene_counts_apop)

# Combine dataframes starting with largest first (probiotic)
C_vir_apop_gene_counts <- left_join(Probiotic_gene_counts_apop, ROD_gene_counts_apop, by = "Parent")
C_vir_apop_gene_counts <- left_join(C_vir_apop_gene_counts,Dermo_gene_counts_apop, by = "Parent")
class(C_vir_apop_gene_counts$Parent) # character
class(C_vir_rtracklayer_apop_product_final$Parent) # AsIs
C_vir_rtracklayer_apop_product_final$Parent <- as.character(C_vir_rtracklayer_apop_product_final$Parent) # change the class to character
View(C_vir_apop_gene_counts)
nrow(C_vir_apop_gene_counts) # 198

# Combine product name 
C_vir_apop_gene_counts <- left_join(C_vir_apop_gene_counts,C_vir_rtracklayer_apop_product_final[,c("Parent","product")], by = "Parent")
View(C_vir_apop_gene_counts)
# extra rows added because of the ", transcript id" 
nrow(C_vir_apop_gene_counts) #283
# Annotate C_vir rownames



# For C_gig: the LOC is the gene column C_gig: Search counts tables for apoptosis genes using the Parent column from genome
# vector C_gig_apop transcript IDs
C_gig_rtracklayer_apop_product_final_gene <- unique(C_gig_rtracklayer_apop_product_final$gene)

nrow(Rubio_gene_counts) #6888
nrow(deLorgeril_gene_counts) # 5023
nrow(He_gene_counts) # 7888
nrow(Zhang_gene_counts)


deLorgeril_Resistant_counts_apop <- deLorgeril_Resistant_counts[row.names(deLorgeril_Resistant_counts) %in% C_gig_rtracklayer_apop_product_final_transcript_id,]
deLorgeril_Susceptible_counts_apop <- deLorgeril_Susceptible_counts[row.names(deLorgeril_Susceptible_counts) %in% C_gig_rtracklayer_apop_product_final_transcript_id,]
nrow(deLorgeril_Resistant_counts_apop ) #659
nrow(deLorgeril_Susceptible_counts_apop ) #659



# Merge tables 

# Search original counts for apoptosis genes and do rlog on just these
Dermo_Tolerant_counts_apop <- Dermo_Tolerant_counts[row.names(Dermo_Tolerant_counts) %in% C_vir_rtracklayer_apop_product_final_ID,]
Dermo_Susceptible_counts_apop <- Dermo_Susceptible_counts[row.names(Dermo_Susceptible_counts) %in% C_vir_rtracklayer_apop_product_final_ID,]
nrow(Dermo_Tolerant_counts_apop) #1026
nrow( Dermo_Susceptible_counts_apop) # 1026

top_Var_Dermo_Tolerant_apop_assay_prot <- as.data.frame(top_Var_Dermo_Tolerant_apop_assay_heatmap_reorder)
colnames(top_Var_Dermo_Tolerant_apop_assay_prot)[1] <- "ID"
top_Var_Dermo_Tolerant_apop_assay_prot_annot <- left_join(top_Var_Dermo_Tolerant_apop_assay_prot, select(C_vir_rtracklayer_apop_product_final, ID, product, gene), by = "ID")



top_Var_Dermo_Susceptible_apop_assay_prot_annot <- left_join(top_Var_Dermo_Susceptible_apop_assay_prot, select(C_vir_rtracklayer_apop_product_final, ID, product, gene), by = "ID")

# unique gene name lists for both species
C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name <- C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique[!duplicated(C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique$gene_name),]
C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name <- C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique[!duplicated(C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique$gene_name),]

# rename gene column for joining
colnames(C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name)[1] <- "C_vir_gene_LOC"
colnames(C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name)[1] <- "C_gig_gene_LOC"

#remove "-like" from  gene name for correct name joining using regex
# this isn't exactly working
# https://datascience.stackexchange.com/questions/8922/removing-strings-after-a-certain-character-in-a-given-text 
C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name$gene_name <-  stringr::str_remove(C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name$gene_name, "\\like$")
C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name$gene_name <-  stringr::str_remove(C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name$gene_name, "\\like$")
C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name$gene_name <-  stringr::str_remove(C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name$gene_name, "\\-$")
C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name$gene_name <-  stringr::str_remove(C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name$gene_name, "\\-$")

# https://stackoverflow.com/questions/50861626/removing-dot-from-the-end-of-string
# full join and places where names don't match will get an NA
combined_gene_name_yes_no_table <- full_join(C_vir_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name, C_gig_rtracklayer_apop_product_final_product_joined_split_product_unique_gene_name,
                                             by = "gene_name")
#remove duplicate gene names
combined_gene_name_yes_no_table_unique <- combined_gene_name_yes_no_table[!duplicated(combined_gene_name_yes_no_table$gene_name),]

#join on the pathway descriptions for the molecules 
# load in table with pathway descriptions for each 

# Still a few discrepancies in the gene names, namely when they are named in C. virginica "* homolog"
# Need to mention in methods that I rmemoved "like" from the end of names
# C_ virginica: lipopolysaccharide-induced tumor necrosis factor-alpha factor homolog, C_gig: lipopolysaccharide-induced tumor necrosis factor-alpha
# C_vir: putative transcription factor p65 homolog, C_gig: transcription factor p65 homolog
# C_vir: macrophage migration inhibitory factor homolog, C_gig: macrophage migration inhibitory factor

# Genes in C vir and not in C gig (all the C_gig_gene_LOC NA's)
c_vir_not_c_gig <- combined_gene_name_yes_no_table_unique %>% filter(is.na(C_gig_gene_LOC))
# Genes in C gig and not in C vir (all the C_vir_gene_LOC NA's)
c_gig_not_c_vir <- combined_gene_name_yes_no_table_unique %>% filter(is.na(C_vir_gene_LOC))
# shared in both 
shared_apoptosis_gene_names <- na.omit(combined_gene_name_yes_no_table_unique)

# Load gene pathway key with protein aliases curated in excel
gene_name_pathway_key_merged <- read.csv("Gene_name_pathway_key.csv", head=TRUE)



## Combine raw counts data frames from C. virginica
# All the rows should be in the same order because I used the same apoptosis data frame to join them
all(rownames(Dermo_Susceptible_counts_apop) == rownames(Dermo_Tolerant_counts_apop)) # TRUE
all(rownames(Dermo_Susceptible_counts_apop) %in% rownames(Dermo_Tolerant_counts_apop)) # TRUE

# Check probiotic table order
all(rownames(Probiotic_counts_apop) == rownames(Dermo_Tolerant_counts_apop)) # FALSE
all(rownames(Probiotic_counts_apop) %in% rownames(Dermo_Tolerant_counts_apop)) # TRUE
Probiotic_counts_apop <- Probiotic_counts_apop[row.names(Dermo_Tolerant_counts_apop),]
all(rownames(Probiotic_counts_apop) == rownames(Dermo_Tolerant_counts_apop)) # TRUE

# Check ROD order and change if necessary 
all(rownames(ROD_Resistant_counts_apop) == rownames(Dermo_Tolerant_counts_apop)) # FALSE
all(rownames(ROD_Resistant_counts_apop) %in% rownames(Dermo_Tolerant_counts_apop)) # TRUE
ROD_Resistant_counts_apop <- ROD_Resistant_counts_apop[row.names(Dermo_Tolerant_counts_apop),]
ROD_Susceptible_counts_apop <-  ROD_Susceptible_counts_apop[row.names(Dermo_Tolerant_counts_apop),]
all(rownames(ROD_Resistant_counts_apop) == rownames(Dermo_Tolerant_counts_apop)) # TRUE
all(rownames(ROD_Susceptible_counts_apop) == rownames(Dermo_Tolerant_counts_apop)) # TRUE
all(rownames(ROD_Susceptible_counts_apop) == rownames(ROD_Resistant_counts_apop)) # TRUE

# Colbind all C. virginica tables
C_virginica_apop_counts <- cbind(Dermo_Susceptible_counts_apop,Dermo_Tolerant_counts_apop,
                                 Probiotic_counts_apop,ROD_Resistant_counts_apop,
                                 ROD_Susceptible_counts_apop)

# Set equal the rownames and colnames of the coldata and count data
all(rownames(C_vir_coldata ) %in% colnames(C_virginica_apop_counts))  #Should return TRUE
# returns TRUE
all(colnames(C_virginica_apop_counts) %in% rownames(C_vir_coldata  ))  
# returns TRUE
all(rownames(C_vir_coldata ) == colnames(C_virginica_apop_counts)) # FALSE

# Fix the order
C_virginica_apop_counts <- C_virginica_apop_counts[,row.names(C_vir_coldata)]

all(rownames(C_vir_coldata ) %in% colnames(C_virginica_apop_counts))  #Should return TRUE
# returns TRUE
all(colnames(C_virginica_apop_counts) %in% rownames(C_vir_coldata  ))  
# returns TRUE
all(rownames(C_vir_coldata ) == colnames(C_virginica_apop_counts))  # TRUE

# Make DEseq data set from matrix so that the coldata gets attached
C_virginica_apop_counts_dds <- DESeqDataSetFromMatrix(countData = C_virginica_apop_counts,
                                                      colData= C_vir_coldata,
                                                      design = ~Condition)
# Collapse technical replicates 
C_virginica_apop_counts_dds <- collapseReplicates(C_virginica_apop_counts_dds, C_virginica_apop_counts_dds$Sample, C_virginica_apop_counts_dds$TechRep)

# Calculate the vst
C_virginica_apop_counts_vst <- varianceStabilizingTransformation(C_virginica_apop_counts_dds)

## Combine C_gig data frame
# Check row order before combining
all(rownames(Zhang_counts_apop) == rownames(Rubio_counts_apop)) # FALSE
all(rownames(Zhang_counts_apop) %in% rownames(Rubio_counts_apop)) # TRUE
Zhang_counts_apop <- Zhang_counts_apop[row.names(Rubio_counts_apop),]
all(rownames(Zhang_counts_apop) == rownames(Rubio_counts_apop)) # TRUE

deLorgeril_Susceptible_counts_apop <- deLorgeril_Susceptible_counts_apop[row.names(Rubio_counts_apop),]
deLorgeril_Resistant_counts_apop <- deLorgeril_Resistant_counts_apop[row.names(Rubio_counts_apop),]
He_counts_apop <- He_counts_apop[row.names(Rubio_counts_apop),]
all(rownames(deLorgeril_Susceptible_counts_apop) == rownames(Zhang_counts_apop)) # TRUE
all(rownames(deLorgeril_Resistant_counts_apop) == rownames(Zhang_counts_apop)) # TRUE
all(rownames(He_counts_apop) == rownames(Zhang_counts_apop)) # TRUE

C_gigas_apop_counts <- cbind(Zhang_counts_apop,Rubio_counts_apop,deLorgeril_Susceptible_counts_apop,
                             deLorgeril_Resistant_counts_apop,He_counts_apop)
# Set equal the rownames and colnames of the coldata and count data
all(rownames(C_gig_coldata ) %in% colnames(C_gigas_apop_counts))  #Should return TRUE
# returns TRUE
all(colnames(C_gigas_apop_counts) %in% rownames(C_gig_coldata  ))  
# returns TRUE
all(rownames(C_gig_coldata ) == colnames(C_gigas_apop_counts)) # FALSE

# Fix the order
C_gigas_apop_counts <- C_gigas_apop_counts[,row.names(C_gig_coldata)]

all(rownames(C_gig_coldata ) %in% colnames(C_gigas_apop_counts))  #Should return TRUE
# returns TRUE
all(colnames(C_gigas_apop_counts) %in% rownames(C_gig_coldata  ))  
# returns TRUE
all(rownames(C_gig_coldata ) == colnames(C_gigas_apop_counts))  # TRUE

# Make DEseq data set from matrix so that the coldata gets attached
C_gigas_apop_dds <- DESeqDataSetFromMatrix(countData = C_gigas_apop_counts ,
                                           colData = C_gig_coldata,
                                           design = ~Condition)
# Calculate the vst
C_gigas_apop_counts_vst <- varianceStabilizingTransformation(C_gigas_apop_dds)

## Remove Batch effects from experiment for C_vir
plotPCA(C_virginica_apop_counts_vst, "Experiment") # grouping by experiment
mat_C_vir <- assay(C_virginica_apop_counts_vst)
mat_C_vir <- limma::removeBatchEffect(mat_C_vir, C_virginica_apop_counts_vst$Experiment)
assay(C_virginica_apop_counts_vst) <- mat_C_vir
plotPCA(C_virginica_apop_counts_vst, "Experiment") # Probiotic and ROD now cluster together
plotPCA(C_virginica_apop_counts_vst, "Sample")
plotPCA(C_virginica_apop_counts_vst, "Time") # no clustering by time
plotPCA(C_virginica_apop_counts_vst, "Family")

## Remove Batch effects from experiment for C_vir
plotPCA(C_gigas_apop_counts_vst, "Experiment") # grouping by experiment, but He and Zhang also cluster
mat_C_gig <- assay(C_gigas_apop_counts_vst)
mat_C_gig <- limma::removeBatchEffect(mat_C_gig, C_gigas_apop_counts_vst$Experiment)
assay(C_gigas_apop_counts_vst) <- mat_C_gig
plotPCA(C_gigas_apop_counts_vst, "Experiment") # He and deLorgeril cluster some, Rubio and Zhang definitely cluster
plotPCA(C_gigas_apop_counts_vst, "Sample")
plotPCA(C_gigas_apop_counts_vst, "Time") # no clustering by time
plotPCA(C_gigas_apop_counts_vst, "Family")

### Heatmaps of apoptosis gene vsts ###

# heatmap of all apoptosis genes 
C_virginica_apop_counts_assay <-  assay(C_virginica_apop_counts_vst)[,]
C_virginica_apop_counts_assay_mat <- C_virginica_apop_counts_assay - rowMeans(C_virginica_apop_counts_assay)
C_virginica_apop_counts_assay_anno <- as.data.frame(colData(C_virginica_apop_counts_vst)[, c("Condition","Time","Experiment")])
C_virginica_apop_counts_assay_heatmap <- pheatmap(C_virginica_apop_counts_assay_mat  , annotation_col = C_virginica_apop_counts_assay_anno)

C_gigas_apop_counts_assay <-  assay(C_gigas_apop_counts_vst)[,]
C_gigas_apop_counts_assay_mat <- C_gigas_apop_counts_assay - rowMeans(C_gigas_apop_counts_assay)
C_gigas_apop_counts_assay_anno <- as.data.frame(colData(C_gigas_apop_counts_vst)[, c("Condition","Experiment")])
C_gigas_apop_counts_assay_heatmap <- pheatmap(C_gigas_apop_counts_assay_mat  , annotation_col = C_gigas_apop_counts_assay_anno)
# clustering of delorgeril and HE

# heatmap of most variable apoptosis genes for C_vir (this selects genes with the greatest variance in the sample)
topVarGenes_C_virginica_apop_assay <-  head(order(rowVars(assay(C_virginica_apop_counts_vst)), decreasing = TRUE), 100) 
top_Var_C_virginica_apop_assay_mat<- assay(C_virginica_apop_counts_vst)[topVarGenes_C_virginica_apop_assay,]
top_Var_C_virginica_apop_assay_mat <- top_Var_C_virginica_apop_assay_mat - rowMeans(top_Var_C_virginica_apop_assay_mat)
top_Var_C_virginica_apop_assay_anno <- as.data.frame(colData(C_virginica_apop_counts_vst)[, c("Experiment","Condition")])
top_Var_C_virginica_apop_assay_heatmap <- pheatmap(top_Var_C_virginica_apop_assay_mat  , annotation_col = top_Var_C_virginica_apop_assay_anno)
head(top_Var_C_virginica_apop_assay_mat ) # some ROD and probiotic clustering..dermo samples mostly cluster together

topVarGenes_C_gigas_apop_assay <-  head(order(rowVars(assay(C_gigas_apop_counts_vst)), decreasing = TRUE), 100) 
top_Var_C_gigas_apop_assay_mat<- assay(C_gigas_apop_counts_vst)[topVarGenes_C_gigas_apop_assay,]
top_Var_C_gigas_apop_assay_mat <- top_Var_C_gigas_apop_assay_mat - rowMeans(top_Var_C_gigas_apop_assay_mat)
top_Var_C_gigas_apop_assay_anno <- as.data.frame(colData(C_gigas_apop_counts_vst)[, c("Family", "Experiment")])
top_Var_C_gigas_apop_assay_heatmap <- pheatmap(top_Var_C_gigas_apop_assay_mat  , annotation_col = top_Var_C_gigas_apop_assay_anno)
head(top_Var_C_gigas_apop_assay_mat ) # OsHV1 susceptible clusters well with HE susceptible

#### DIFFERENTIAL EXPRESSION ANALYSIS ####


### GENE CLUSTERING ANALYSIS HEATMAPS  
# Extract genes with the highest variance across samples for each comparison using either vst or rlog transformed data
# This heatmap rather than plotting absolute expression strength plot the amount by which each gene deviates in a specific sample from the gene’s average across all samples. 
# example codes from RNAseq workflow: https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#other-comparisons

Pro_RE22_dds_deseq_Challenge_res_LFC_sig_assay <-  head(order(rowVars(assay(Pro_RE22_dds_rlog  )), decreasing = TRUE), 200)
family_Pro_RE22_broken_mat <- assay(Pro_RE22_dds_rlog )[Pro_RE22_dds_deseq_Challenge_res_LFC_sig_assay ,]
family_Pro_RE22_broken_mat <- family_Pro_RE22_broken_mat - rowMeans(family_Pro_RE22_broken_mat)
family_Pro_RE22_broken_anno <- as.data.frame(colData(Pro_RE22_dds_rlog )[, c("Condition","Time")])
family_Pro_RE22_broken_heatmap <- pheatmap(family_Pro_RE22_broken_mat , annotation_col = family_Pro_RE22_broken_anno)
head(family_Pro_RE22_broken_mat) # 

# Gene clustering heatmap with only apoptosis genes #
# Search original Probiotic_counts for apoptosis genes and do rlog on just these
Pro_RE22_counts_apop <- Pro_RE22_counts[row.names(Pro_RE22_counts) %in% C_vir_rtracklayer_apop_product_final_ID,]
nrow(Pro_RE22_counts_apop) #1026
head(Pro_RE22_counts_apop)
Pro_RE22_counts_apop_dds <- DESeqDataSetFromMatrix(countData = Pro_RE22_counts_apop,
                                                   colData = Pro_RE22_coldata,
                                                   design = ~Time + Condition) # add time to control for larval age effect 
# Prefiltering the data and running rlog
Pro_RE22_counts_apop_dds<- Pro_RE22_counts_apop_dds[ rowSums(counts(Pro_RE22_counts_apop_dds)) > 10, ]
Pro_RE22_counts_apop_dds_rlog <- rlog(Pro_RE22_counts_apop_dds, blind=TRUE)

## PCA plot of rlog transformed counts for apoptosis
plotPCA(Pro_RE22_counts_apop_dds_rlog , intgroup="Time") # still overall clustering by time and not condition

# heatmap of all apoptosis genes 
Pro_RE22_counts_apop_assay <-  assay(Pro_RE22_counts_apop_dds_rlog)[,]
Pro_RE22_counts_apop_assay_mat <- Pro_RE22_counts_apop_assay - rowMeans(Pro_RE22_counts_apop_assay)
Pro_RE22_counts_apop_assay_anno <- as.data.frame(colData(Pro_RE22_counts_apop_dds_rlog )[, c("Condition","Sample")])
Pro_RE22_counts_apop_assay_heatmap <- pheatmap(Pro_RE22_counts_apop_assay_mat  , annotation_col = Pro_RE22_counts_apop_assay_anno)
head(Pro_RE22_counts_apop_assay_mat ) # 

# heatmap of most variable apoptosis genes (this selects genes with the greatest variance in the sample)
topVarGenes_Pro_RE22_counts_apop_assay <-  head(order(rowVars(assay(Pro_RE22_counts_apop_dds_rlog)), decreasing = TRUE), 100) 
top_Var_Pro_RE22_counts_apop_assay_mat<- assay(Pro_RE22_counts_apop_dds_rlog)[topVarGenes_Pro_RE22_counts_apop_assay,]
top_Var_Pro_RE22_counts_apop_assay_mat <- top_Var_Pro_RE22_counts_apop_assay_mat - rowMeans(top_Var_Pro_RE22_counts_apop_assay_mat)
top_Var_Pro_RE22_counts_apop_assay_anno <- as.data.frame(colData(Pro_RE22_counts_apop_dds_rlog)[, c("Condition","Time")])
top_Var_Pro_RE22_counts_apop_assay_heatmap <- pheatmap(top_Var_Pro_RE22_counts_apop_assay_mat  , annotation_col = top_Var_Pro_RE22_counts_apop_assay_anno)
head(top_Var_Pro_RE22_counts_apop_assay_mat )
# some clustering patterns here in signature

# annotate the top 100 most variable genes  
# reorder annotation table to match ordering in heatmap 
top_Var_Pro_RE22_counts_apop_assay_heatmap_reorder <-rownames(top_Var_Pro_RE22_counts_apop_assay_mat[top_Var_Pro_RE22_counts_apop_assay_heatmap$tree_row[["order"]],])
# annotate the row.names
top_Var_Pro_RE22_counts_apop_assay_prot <- as.data.frame(top_Var_Pro_RE22_counts_apop_assay_heatmap_reorder)
colnames(top_Var_Pro_RE22_counts_apop_assay_prot)[1] <- "ID"
top_Var_Pro_RE22_counts_apop_assay_prot_annot <- left_join(top_Var_Pro_RE22_counts_apop_assay_prot, select(C_vir_rtracklayer_apop_product_final, ID, product, gene), by = "ID")

### Extract list of significant Apoptosis Genes (not less than or greater than 1 LFC) using merge
Pro_RE22_dds_deseq_Challenge_res_LFC_sig_APOP <- merge(Pro_RE22_dds_deseq_Challenge_res_LFC_sig, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_dds_deseq_Challenge_res_LFC_sig_APOP_arranged <- arrange(Pro_RE22_dds_deseq_Challenge_res_LFC_sig_APOP, -log2FoldChange) 
nrow(Pro_RE22_dds_deseq_Challenge_res_LFC_sig_APOP) # 

Pro_RE22_dds_deseq_Challenge_res_LFC_sig_APOP_plot <- ggplot(Pro_RE22_dds_deseq_Challenge_res_LFC_sig_APOP , aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("") +
  ylab("Log2 Fold Change")
