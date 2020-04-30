## R script to run TOMsimilarityFromExpr on the bluewaves cluster 
## Erin Michele Roberts, PhD Candidate

library(WGCNA)

options(stringsAsFactors = FALSE) # run every time
allowWGCNAThreads()
cor <- WGCNA::cor

#### LOAD SAVED DATA ####
Pro_RE22_dds_rlog_matrix_RE22 <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_dds_rlog_matrix_RE22.table")
#print(class(Pro_RE22_dds_rlog_matrix_RE22))
#print(head(Pro_RE22_dds_rlog_matrix_RE22))
#
Dermo_Tolerant_dds_vst_matrix <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Dermo_Tolerant_dds_vst_matrix.table")
#print(class(Dermo_Tolerant_dds_vst_matrix))
#print(head(Dermo_Tolerant_dds_vst_matrix))

#Dermo_Tol_full_TOM = TOMsimilarityFromExpr(Dermo_Tolerant_dds_vst_matrix, power = 3, TOMType = "signed", networkType= "signed hybrid", corType = "bicor") 
#
#save(Dermo_Tol_full_TOM, file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Dermo_Tol_full_TOM.RData")
#
#Pro_RE22_RE22_full_TOM = TOMsimilarityFromExpr(Pro_RE22_dds_rlog_matrix_RE22, power = 9, TOMType = "signed", networkType= "signed hybrid", corType ="bicor" )
#
#save(Pro_RE22_RE22_full_TOM, file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_RE22_full_TOM.RData" )
#


###### RUNNING THE FOLLOWING CODE IN AN INTERACTIVE R SESSION ON BLUEWAVES ####

# $ pwd 
# /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA
# $ interactive 
# $ module load R/3.6.0-intel-2019a
# $ R

library(WGCNA)

# Read in the annotation files 
load(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/C_gig_C_vir_annotations.RData")

# Check they loaded
nrow(C_vir_rtracklayer)

# Read in the TOM files as tables
load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Dermo_Tol_full_TOM.RData")
load(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_RE22_full_TOM.RData")

# Read in the matrix files
Pro_RE22_dds_rlog_matrix_RE22 <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_dds_rlog_matrix_RE22.table")
Dermo_Tolerant_dds_vst_matrix <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Dermo_Tolerant_dds_vst_matrix.table")

# Read in the module colors files
# export moduleColors file for use in cluster
load(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Dermo_Tol_full_moduleColors.RData")
load(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_RE22_full_net_moduleColors.RData")

class(Dermo_Tol_full_moduleColors) # character

## Run the exportnetwork to cytoscape command
Dermo_Tol_full_modules = c("darkslateblue", "turquoise", "greenyellow",   "skyblue3" , "cyan"  ,"red"  , "tan" )
# Select module probes
Dermo_Tol_full_probes = colnames(Dermo_Tolerant_dds_vst_matrix)
Dermo_Tol_full_inModule = is.finite(match(Dermo_Tol_full_moduleColors, Dermo_Tol_full_modules))
Dermo_Tol_full_modProbes = Dermo_Tol_full_probes[Dermo_Tol_full_inModule]
Dermo_Tol_full_modGenes = C_vir_rtracklayer$ID[match(Dermo_Tol_full_modProbes, C_vir_rtracklayer$ID)]
# Select the corresponding Topological Overlap
Dermo_Tol_full_modTOM = Dermo_Tol_full_TOM[Dermo_Tol_full_inModule, Dermo_Tol_full_inModule]
dimnames(Dermo_Tol_full_modTOM) = list(Dermo_Tol_full_modProbes, Dermo_Tol_full_modProbes)
# Export the network into edge and node list files Cytoscape can read
Dermo_Tol_full_cyt = exportNetworkToCytoscape(Dermo_Tol_full_modTOM,
  edgeFile = paste("CytoscapeInput-edges-Dermo_Tol_full", paste(Dermo_Tol_full_modules, collapse="-"), ".txt", sep=""),
                                              nodeFile = paste("CytoscapeInput-nodes-Dermo_Tol_full", paste(Dermo_Tol_full_modules, collapse="-"), ".txt", sep=""),
                                              weighted = TRUE,
                                              threshold = 0.00, # using 0 threshold so no genes are subset out 
                                              nodeNames = Dermo_Tol_full_modProbes,
                                              altNodeNames = Dermo_Tol_full_modGenes,
                                              nodeAttr = Dermo_Tol_full_moduleColors[Dermo_Tol_full_inModule])
str(Dermo_Tol_full_cyt)
#List of 2
#$ edgeData:'data.frame':	7285542 obs. of  6 variables:
#  ..$ fromNode   : Factor w/ 6469 levels "gene10065","gene10438",..: 2072 2072 2072 2072 2072 2072 2072 2072 2072 2072 ...
#..$ toNode     : Factor w/ 6467 levels "gene10065","gene10438",..: 2069 2068 5470 2071 1108 5198 3243 2843 5466 5352 ...
#..$ weight     : num [1:7285542] 0.027 0.0295 0.0224 0.074 0.0242 ...
#..$ direction  : Factor w/ 1 level "undirected": 1 1 1 1 1 1 1 1 1 1 ...
#..$ fromAltName: Factor w/ 6469 levels "gene10065","gene10438",..: 2072 2072 2072 2072 2072 2072 2072 2072 2072 2072 ...
#..$ toAltName  : Factor w/ 6467 levels "gene10065","gene10438",..: 2069 2068 5470 2071 1108 5198 3243 2843 5466 5352 ...
#$ nodeData:'data.frame':	6478 obs. of  3 variables:
#  ..$ nodeName                : Factor w/ 6478 levels "gene10065","gene10438",..: 2073 2071 2070 2069 5480 2075 2074 1108 1107 1109 ...
#..$ altName                 : Factor w/ 6478 levels "gene10065","gene10438",..: 2073 2071 2070 2069 5480 2075 2074 1108 1107 1109 ...
#..$ nodeAttr[nodesPresent, ]: Factor w/ 7 levels "cyan","darkslateblue",..: 7 7 7 2 7 7 7 7 7 4 ...

# The command writes results to file automatically 

#### Repeat for the Pro_RE22_RE22 modules ####

# Select modules
Pro_RE22_RE22_full_modules = c("darkgreen",  "tan" ,       "turquoise",  "darkorange" )
# Select module probes
Pro_RE22_RE22_full_probes = colnames(Pro_RE22_dds_rlog_matrix_RE22)
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
                                                  threshold = 0.00, # use zero threshold so no genes are subset out 
                                                  nodeNames = Pro_RE22_RE22_full_modProbes,
                                                  altNodeNames = Pro_RE22_RE22_full_modGenes,
                                                  nodeAttr = Pro_RE22_RE22_full_net_moduleColors[Pro_RE22_RE22_full_inModule])

# The command writes results to file automatically 

sessionInfo()
#sessionInfo()
#R version 3.6.0 (2019-04-26)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: CentOS release 6.5 (Final)

#Matrix products: default
#BLAS:   /net/clusterhn.cluster.com/opt/software/R/3.6.0-intel-2019a/lib64/R/lib/libR.so
#LAPACK: /net/clusterhn.cluster.com/opt/software/R/3.6.0-intel-2019a/lib64/R/modules/lapack.so

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#[3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#[5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#[7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#[9] LC_ADDRESS=C               LC_TELEPHONE=C            
#[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#
#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     
#
#other attached packages:
#  [1] WGCNA_1.69            fastcluster_1.1.25    dynamicTreeCut_1.63-1
#
#loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.1            lattice_0.20-38       GO.db_3.10.0         
#[4] assertthat_0.2.1      digest_0.6.19         foreach_1.4.4        
#[7] R6_2.4.0              plyr_1.8.4            backports_1.1.4      
#[10] acepack_1.4.1         stats4_3.6.0          RSQLite_2.1.1        
#[13] ggplot2_3.1.1         pillar_1.4.1          rlang_0.3.4          
#[16] lazyeval_0.2.2        rstudioapi_0.10       data.table_1.12.2    
#[19] blob_1.1.1            S4Vectors_0.24.4      rpart_4.1-15         
#[22] Matrix_1.2-17         preprocessCore_1.48.0 checkmate_1.9.3      
#[25] splines_3.6.0         stringr_1.4.0         foreign_0.8-71       
#[28] htmlwidgets_1.3       bit_1.1-14            munsell_0.5.0        
#[31] compiler_3.6.0        xfun_0.7              pkgconfig_2.0.2      
#[34] BiocGenerics_0.32.0   base64enc_0.1-3       htmltools_0.3.6      
#[37] nnet_7.3-12           tidyselect_0.2.5      tibble_2.1.3         
#[40] gridExtra_2.3         htmlTable_1.13.1      Hmisc_4.2-0          
#[43] IRanges_2.20.2        codetools_0.2-16      matrixStats_0.54.0   
#[46] crayon_1.3.4          dplyr_0.8.1           grid_3.6.0           
#[49] gtable_0.3.0          DBI_1.0.0             magrittr_1.5         
#[52] scales_1.0.0          stringi_1.4.3         impute_1.60.0        
#[55] doParallel_1.0.14     latticeExtra_0.6-28   Formula_1.2-3        
#[58] RColorBrewer_1.1-2    iterators_1.0.10      tools_3.6.0          
#[61] bit64_0.9-7           Biobase_2.46.0        glue_1.3.1           
#[64] purrr_0.3.2           parallel_3.6.0        survival_2.44-1.1    
#[67] AnnotationDbi_1.48.0  colorspace_1.4-1      cluster_2.0.9        
#[70] memoise_1.1.0         knitr_1.23  #