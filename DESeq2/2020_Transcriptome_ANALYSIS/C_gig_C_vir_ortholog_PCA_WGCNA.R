# Script to Analyze C. virginica and C. gigas orthologs involved in apoptosis
# Orthologs for apoptosis will be input into PCA and WGCNA
# Erin Roberts, PhD Candidate University of Rhode Island 
# 2/26/2020


# Helpful tutorial for OrthoFinder
# OrthoFinder Manual: https://github.com/davidemms/OrthoFinder/blob/master/README.md
# http://arken.nmbu.no/~larssn/teach/bin310/week8.html

# Load packages 
library(tidyverse)

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
head(Pro_RE22_gene_counts)
head(Rubio_gene_counts) # header need the "_1" removed
head(deLorgeril_gene_counts) # header need the "_1" removed
head(He_gene_counts )
head(ROD_gene_counts)
head(Zhang_gene_counts )

# # colnames all have "_1", remove this. It was an artifact of the PE sample names
colnames(Probiotic_gene_counts ) <- sub('\\_[^_]+$', '', colnames(Probiotic_gene_counts))
colnames(Rubio_gene_counts ) <- sub('\\_[^_]+$', '', colnames(Rubio_gene_counts))
colnames(deLorgeril_gene_counts) <- sub('\\_[^_]+$', '', colnames(deLorgeril_gene_counts))

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

Pro_RE22_gene_counts_single_copy <- Pro_RE22_gene_counts[Pro_RE22_gene_counts$gene %in% C_vir_orthogroups_single_copy_unique$gene]
Pro_RE22_gene_counts_single_copy <- left_join(Pro_RE22_gene_counts_single_copy, C_vir_orthogroups_single_copy_unique)
nrow(Pro_RE22_gene_counts_single_copy) #771

ROD_gene_counts_single_copy <- ROD_gene_counts[ROD_gene_counts$gene %in% C_vir_orthogroups_single_copy_unique$gene]
ROD_gene_counts_single_copy <- left_join(ROD_gene_counts_single_copy, C_vir_orthogroups_single_copy_unique)
nrow(ROD_gene_counts_single_copy) #771

Rubio_gene_counts_single_copy  <- ROD_gene_counts[ROD_gene_counts$gene %in% C_gig_orthogroups_single_copy_unique$gene]
Rubio_gene_counts_single_copy  <- left_join(Rubio_gene_counts_single_copy, C_gig_orthogroups_single_copy_unique)
nrow(Rubio_gene_counts_single_copy ) #771

deLorgeril_gene_counts_single_copy <- deLorgeril_gene_counts[deLorgeril_gene_counts$gene %in% C_gig_orthogroups_single_copy_unique$gene]
deLorgeril_gene_counts_single_copy <- left_join(deLorgeril_gene_counts_single_copy, C_gig_orthogroups_single_copy_unique)
nrow(deLorgeril_gene_counts_single_copy ) #771

He_gene_counts_single_copy <- He_gene_counts[He_gene_counts$gene %in% C_gig_orthogroups_single_copy_unique$gene]
He_gene_counts_single_copy <- left_join(He_gene_counts_single_copy, C_gig_orthogroups_single_copy_unique)
nrow(He_gene_counts_single_copy ) #771

Zhang_gene_counts_single_copy <- Zhang_gene_counts[Zhang_gene_counts$gene %in% C_gig_orthogroups_single_copy_unique$gene]
Zhang_gene_counts_single_copy <- left_join(Zhang_gene_counts_single_copy, C_gig_orthogroups_single_copy_unique)
nrow(Zhang_gene_counts_single_copy ) #771

# Join within species based on gene and orthogroup, starting with one that has the most 
C_vir_ortholog_gene_counts <- left_join(Probiotic_gene_counts_single_copy,ROD_gene_counts_single_copy, by =c("gene","Orthogroup"))

Dermo_gene_counts_single_copy
Probiotic_gene_counts_single_copy
Pro_RE22_gene_counts_single_copy
ROD_gene_counts_single_copy
Rubio_gene_counts_single_copy 
deLorgeril_gene_counts_single_copy 
He_gene_counts_single_copy 
Zhang_gene_counts_single_copy 

# merge starting with the largest first 
Cvir

nrow(Dermo_gene_counts_single_copy )
nrow(Probiotic_gene_counts_single_copy)
nrow(ROD_gene_counts_single_copy)
nrow(Rubio_gene_counts_single_copy )
nrow(deLorgeril_gene_counts_single_copy)
nrow(He_gene_counts_single_copy )
nrow(Zhang_gene_counts_single_copy )

## Merge data frames based on matching orthologroup IDs to get full table of counts 







#### COMPARING APOPTOSIS GENE EXPRESSION BETWEEN EXPERIMENTS PCA HEATMAPS FULL THEN SUBSET ####


# split product 

C_vir_full_gene_counts <- C_vir_full_gene_counts[,-7] # remove rownames to allow for vst 
colnames(C_vir_full_gene_counts)
head(C_vir_full_gene_counts)


C_gig_full_gene_counts <- left_join(Zhang_gene_counts,He_gene_counts, by ="rownames")
C_gig_full_gene_counts <- left_join(C_gig_full_gene_counts,Rubio_gene_counts, by = "rownames")
C_gig_full_gene_counts <- left_join(C_gig_full_gene_counts,deLorgeril_gene_counts, by = "rownames")
colnames(C_gig_full_gene_counts)
row.names(C_gig_full_gene_counts) <- C_gig_full_gene_counts$rownames
C_gig_full_gene_counts <- C_gig_full_gene_counts[,-10] # remove rownames to allow for vst 

# Combine dataframes starting with largest first (probiotic)
C_vir_apop_gene_counts <- left_join(Probiotic_gene_counts_apop, ROD_gene_counts_apop, by = "Parent")
C_vir_apop_gene_counts <- left_join(C_vir_apop_gene_counts,Dermo_gene_counts_apop, by = "Parent")
class(C_vir_apop_gene_counts$Parent) # character
class(C_vir_rtracklayer_apop_product_final$Parent) # AsIs
C_vir_rtracklayer_apop_product_final$Parent <- as.character(C_vir_rtracklayer_apop_product_final$Parent) # change the class to character
View(C_vir_apop_gene_counts)
nrow(C_vir_apop_gene_counts) # 198


# Match colnames for each to their product names in order to combine by gene name

# For C_gig: the LOC is the gene column C_gig: Search counts tables for apoptosis genes using the Parent column from genome
# vector C_gig_apop transcript IDs
C_gig_rtracklayer_apop_product_final_gene <- unique(C_gig_rtracklayer_apop_product_final$gene)

nrow(Rubio_gene_counts) #6888
nrow(deLorgeril_gene_counts) # 5023
nrow(He_gene_counts) # 7888
nrow(Zhang_gene_counts)



# Set equal the rownames and colnames of the coldata and count data
all(rownames(C_vir_coldata ) %in% colnames(C_vir_full_gene_counts))  #Should return TRUE
# returns TRUE
all(colnames(C_vir_full_gene_counts) %in% rownames(C_vir_coldata  ))  
# returns TRUE
all(rownames(C_vir_coldata ) == colnames(C_vir_full_gene_counts)) # FALSE

# Fix the order (already in correct order)
C_vir_full_gene_counts <- C_vir_full_gene_counts[,row.names(C_vir_coldata)]

# C_gig
all(rownames(C_gig_coldata ) %in% colnames(C_gig_full_counts ))  #Should return TRUE
# returns TRUE
all(colnames(C_gig_full_counts ) %in% rownames(C_gig_coldata  ))  
# returns TRUE
all(rownames(C_gig_coldata ) == colnames(C_gig_full_counts )) # TRUE

# Fix the order
C_gig_full_gene_counts  <- C_gig_full_gene_counts[,row.names(C_gig_coldata)]
all(rownames(C_gig_coldata ) %in% colnames(C_gig_full_gene_counts ))  #Should return TRUE
# returns TRUE
all(colnames(C_gig_full_gene_counts ) %in% rownames(C_gig_coldata  ))  
# returns TRUE
all(rownames(C_gig_coldata ) == colnames(C_gig_full_gene_counts )) # TRUE

# Change NA values to be zero for both 
C_vir_full_gene_counts[is.na(C_vir_full_gene_counts)] <- 0
C_gig_full_gene_counts[is.na(C_gig_full_gene_counts)] <- 0




# Make DEseq data set from matrix so that the coldata gets attached
C_vir_full_gene_counts_dds <- DESeqDataSetFromMatrix(countData = C_vir_full_gene_counts,
                                                     colData= C_vir_coldata,
                                                     design = ~Condition)
# Collapse technical replicates 
C_vir_full_counts_dds <- collapseReplicates(C_vir_full_counts_dds, C_vir_full_counts_dds$Sample, C_vir_full_counts_dds$TechRep)

# Calculate the vst
C_vir_full_counts_vst <- varianceStabilizingTransformation(C_vir_full_counts_dds)

# Make DEseq data set from matrix so that the coldata gets attached
C_gig_full_counts_dds <- DESeqDataSetFromMatrix(countData = C_gig_full_counts ,
                                                colData = C_gig_coldata,
                                                design = ~Condition)
# Calculate the vst
C_gig_full_counts_vst <- varianceStabilizingTransformation(C_gig_full_counts_dds)

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

## LOAD DATA
Pro_RE22_counts <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/Probiotic_RE22_transcript_count_matrix.csv", header=TRUE,
                            row.names = "transcript_id")
head(Pro_RE22_counts )
colnames(Pro_RE22_counts )

# colnames all have "_1", remove this. It was an artifact of the PE sample names
colnames(Pro_RE22_counts  ) <- sub('\\_[^_]+$', '', colnames(Pro_RE22_counts ))
colnames(Pro_RE22_counts )

# remove MSTRG novel transcript lines (can assess these later if necessary)
Pro_RE22_counts  <- Pro_RE22_counts [!grepl("MSTRG", row.names(Pro_RE22_counts )),]
row.names(Pro_RE22_counts ) <- remove_rna(row.names(Pro_RE22_counts ))
head(Pro_RE22_counts )

#Load in sample metadata
Pro_RE22_coldata <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/Probiotic_coldata.csv", row.names = 1 )
View(Pro_RE22_coldata)  
nrow(Pro_RE22_coldata) 

# Make sure the columns of the count matrix and rows of the column data (sample metadata) are in the same order. 
all(rownames(Pro_RE22_coldata) %in% colnames(Pro_RE22_counts))  #Should return TRUE
# returns TRUE
all(colnames(Pro_RE22_counts) %in% rownames(Pro_RE22_coldata))  
# returns TRUE
all(rownames(Pro_RE22_coldata) == colnames(Pro_RE22_counts))
# returns FALSE

# Fix the order
Pro_RE22_counts <-Pro_RE22_counts[,row.names(Pro_RE22_coldata)]
row.names(Pro_RE22_coldata)

all(rownames(Pro_RE22_coldata) %in% colnames(Pro_RE22_counts))  #Should return TRUE
# returns TRUE
all(colnames(Pro_RE22_counts) %in% rownames(Pro_RE22_coldata))  
# returns TRUE
all(rownames(Pro_RE22_coldata) == colnames(Pro_RE22_counts))
# returns TRUE

### DATA QC PCA PLOT 
# rlog transform data is recommended over vst for small data sets 
# PCA plots of data (https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/02_Preprocessing_Data.nb.html#count-distribution-boxplots)
Pro_RE22_counts_matrix <- as.matrix(Pro_RE22_counts)
Pro_RE22rlogcounts <- rlog(Pro_RE22_counts_matrix, blind =TRUE)

# run PCA
pcPro_RE22 <- prcomp(t(Pro_RE22rlogcounts))

# Plot PCA
autoplot(pcPro_RE22,
         data = Pro_RE22_coldata, 
         colour="Condition", 
         size=5) 

# Plot PCA 2 and 3 for comparison
autoplot(pcPro_RE22,
         data = Pro_RE22_coldata, 
         colour = "Condition", 
         size = 5,
         x = 2,
         y = 3) 
# Plot PCA 4 and 4 for comparison
autoplot(pcPro_RE22,
         data = Pro_RE22_coldata, 
         colour = "Condition", 
         size = 5,
         x = 3,
         y = 4)  # this PCA axis shows some clustering of samples by treatment

## MAKE DESEQ DATA SET FROM MATRIX
# This object specifies the count data and metadata you will work with. The design piece is critical.
# Correct for batch effects if necessary in this original formula: see this thread https://support.bioconductor.org/p/121408/
# do not correct counts using the removeBatchEffects from limma based on thread above 

## Check levels 
# It is prefered in R that the first level of a factor be the reference level for comparison
# (e.g. control, or untreated samples), so we can relevel the factor like so
# Check factor levels, set it so that comparison group is the first
levels(Pro_RE22_coldata$Condition) # 
Pro_RE22_coldata$Condition <- factor(Pro_RE22_coldata$Condition , levels = c())
levels(Pro_RE22_coldata$Condition)
levels(Pro_RE22_coldata$Time) # Time this is what I will control for
Pro_RE22_coldata$Time <- factor(Pro_RE22_coldata$Time , levels = c())
levels(Pro_RE22_coldata$Time)

## Creating deseq data set from matrix, controlling for the effect of time
Pro_RE22_dds <- DESeqDataSetFromMatrix(countData = Pro_RE22_counts,
                                       colData = Pro_RE22_coldata,
                                       design = ~ Time + Condition ) 

## Prefiltering the data
# Data prefiltering helps decrease the size of the data set and get rid of
# rows with no data or very minimal data (<10). Apply a minimal filtering here as more stringent filtering will be applied later
Pro_RE22_dds <- Pro_RE22_dds [ rowSums(counts(Pro_RE22_dds )) > 10, ]

## DATA TRANSFORMATION AND VISUALIZATION
# Assess sample clustering after setting initial formula for comparison
Pro_RE22_dds_rlog <- rlog(Pro_RE22_dds , blind = TRUE) # keep blind = true before deseq function has been run

## PCA plot visualization of individuals in the family 
plotPCA(Pro_RE22_dds_rlog, intgroup=c("Sample", "Condition")) # clustering by time is not as tight

### DIFFERENTIAL EXPRESSION ANALYSIS
# run pipeline with single command because the formula has already been specified
# Steps: estimation of size factors (controlling for differences in the sequencing depth of the samples), 
# the estimation of dispersion values for each gene,and fitting a generalized linear model.
Pro_RE22_dds_deseq <- DESeq(Pro_RE22_dds) 

## Check the resultsNames object of each to look at the available coefficients for use in lfcShrink command
resultsNames(Pro_RE22_dds_deseq) #

## BUILD THE RESULTS OBJECT
# Examining the results object, change alpha to p <0.05, looking at object metadata
# use mcols to look at metadata for each table
mcols(Pro_RE22_dds_deseq)
Pro_RE22_dds_deseq_Challenge_res <- results(Pro_RE22_dds_deseq, alpha=0.05, name= "")
head(Pro_RE22_dds_deseq_Challenge_res) # 

### Perform LFC Shrinkage with apeglm
## NOTES 
# Before plotting we need to apply an LFC shrinkage, set the coef as the specific comparison in the ResultsNames function of the deseq object
# Issue that the specific comparisons I set in my results formulas are not available in the ResultsNames coef list.
# More detailed notes about LFC Shrinkage are in the code for the Zhang Vibrio

## DECISION: USE SAME RES OBJECT TO KEEP ALPHA ADJUSTMENT, and use LFCShrink apeglm
Pro_RE22_dds_deseq_Challenge_res_LFC<- lfcShrink(Pro_RE22_dds_deseq, coef="", type="apeglm", res= Pro_RE22_dds_deseq_Challenge_res)

## EXPLORATORY PLOTTING OF RESULTS 
## MA Plotting
plotMA(Pro_RE22_dds_deseq_Challenge_res_LFC, ylim = c(-5, 5))

## Histogram of P values 
# exclude genes with very small counts to avoid spikes and plot using the LFCshrinkage
hist(Pro_RE22_dds_deseq_Challenge_res_LFC$padj[Pro_RE22_dds_deseq_Challenge_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

### Subsetting Significant Genes by padj < 0.05
# again, only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
Pro_RE22_dds_deseq_Challenge_res_LFC_sig <-  subset(Pro_RE22_dds_deseq_Challenge_res_LFC , padj < 0.05)
Pro_RE22_dds_deseq_Challenge_res_LFC_sig$ID<- row.names(Pro_RE22_dds_deseq_Challenge_res_LFC_sig)
Pro_RE22_dds_deseq_Challenge_res_LFC_sig <- as.data.frame(Pro_RE22_dds_deseq_Challenge_res_LFC_sig)
nrow(Pro_RE22_dds_deseq_Challenge_res_LFC_sig) # 1762

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
