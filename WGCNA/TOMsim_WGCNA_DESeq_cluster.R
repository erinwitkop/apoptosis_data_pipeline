## R script to run TOMsimilarityFromExpr on the bluewaves cluster 
## Erin Michele Roberts, PhD Candidate

# Clear global workspace
rm(list = ls())

# Load libraries
library(WGCNA)

options(stringsAsFactors = FALSE) # run every time
allowWGCNAThreads()
cor <- WGCNA::cor

## LIST OF MODULES TO EXPORT ##
# Which modules? - updated Sept. 21st, 2020 to add 4 modules 
# deLorg_Res: MEturquoise
# deLorg_Sus: MEturquoise
# Dermo_Tol: MEturquoise
# Dermo_Sus: MElightpink4 
# Pro_RE22_Pro_RI:  MEskyblue3,  MEsteelblue
# Pro_RE22_Pro_S4:  MEwhite, MEsteelblue
# He:  MEpurple, MEyellow
# Rubio: MEmagenta, MEturquoise, MEblue, MEbrown, MEblack   
# Zhang: MEblack

## Load saved matrices one by one and then unload to free up memory 

#### CALCULATE TOM ####
# power taken from original WGCNA run 

# redoing the Dermo TOM because data updated
#Dermo_Tolerant_dds_vst_matrix <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Dermo_Tolerant_dds_vst_matrix.table")
#Dermo_Tol_full_TOM = TOMsimilarityFromExpr(Dermo_Tolerant_dds_vst_matrix, power = 3, TOMType = "signed", networkType= "signed hybrid", corType = "bicor") 
#save(Dermo_Tol_full_TOM, file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Dermo_Tol_full_TOM.RData")

# Dermo Sus
#Dermo_Susceptible_dds_vst_matrix <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Dermo_Susceptible_dds_vst_matrix.table")
#Dermo_Sus_full_TOM = TOMsimilarityFromExpr(Dermo_Susceptible_dds_vst_matrix, power = 6, TOMType = "signed", networkType= "signed hybrid", corType = "bicor") 
#save(Dermo_Sus_full_TOM, file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Dermo_Sus_full_TOM.RData")
#rm(Dermo_Sus_full_TOM)
#rm(Dermo_Susceptible_dds_vst_matrix)

# C. virginica experiments
#Pro_RE22_dds_rlog_matrix_Pro <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_dds_rlog_matrix_Pro.table")
#Pro_RE22_Pro_full_TOM = TOMsimilarityFromExpr(Pro_RE22_dds_rlog_matrix_Pro, power = 4, TOMType = "signed", networkType= "signed hybrid", corType = "bicor") 
#save(Pro_RE22_Pro_full_TOM , file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_Pro_full_TOM.RData")
#rm(Pro_RE22_Pro_full_TOM) # remove from workspace to free up memory
#rm(Pro_RE22_dds_rlog_matrix_Pro)

## C. gigas experiments
#Zhang_dds_rlog_matrix <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Zhang_dds_rlog_matrix.table")
#Zhang_full_TOM = TOMsimilarityFromExpr(Zhang_dds_rlog_matrix, power = 4, TOMType = "signed", networkType= "signed hybrid", corType = "bicor") 
#save(Zhang_full_TOM , file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Zhang_full_TOM.RData")
#rm(Zhang_full_TOM)
#rm(Zhang_dds_rlog_matrix)
#
#Rubio_dds_rlog_matrix <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Rubio_dds_rlog_matrix.table")
#Rubio_full_TOM = TOMsimilarityFromExpr(Rubio_dds_rlog_matrix, power = 7, TOMType = "signed", networkType= "signed hybrid", corType = "bicor") 
#save(Rubio_full_TOM, file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Rubio_full_TOM.RData")
#rm(Rubio_full_TOM)
#rm(Rubio_dds_rlog_matrix)
#
#deLorgeril_Resistant_dds_vst_matrix <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/deLorgeril_Resistant_dds_vst_matrix.table")
#deLorg_Res_full_TOM = TOMsimilarityFromExpr(deLorgeril_Resistant_dds_vst_matrix, power = 8,TOMType = "signed", networkType= "signed hybrid", corType = "bicor") 
#save(deLorg_Res_full_TOM, file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/deLorg_Res_full_TOM.RData")
#rm(deLorgeril_Resistant_dds_vst_matrix)
#rm(deLorg_Res_full_TOM)
#
#deLorgeril_Susceptible_dds_vst_matrix <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/deLorgeril_Susceptible_dds_vst_matrix.table")
#deLorg_Sus_full_TOM = TOMsimilarityFromExpr(deLorgeril_Susceptible_dds_vst_matrix, power = 3, TOMType = "signed", networkType= "signed hybrid", corType = "bicor") 
#save(deLorg_Sus_full_TOM, file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/deLorg_Sus_full_TOM.RData")
#rm(deLorgeril_Susceptible_dds_vst_matrix)
#rm(deLorg_Sus_full_TOM)
#
#He_dds_vst_matrix <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/He_dds_vst_matrix.table")
#He_full_TOM = TOMsimilarityFromExpr(He_dds_vst_matrix, power = 5,TOMType = "signed", networkType= "signed hybrid", corType = "bicor") 
#save(He_full_TOM, file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/He_full_TOM.RData")
#rm(He_dds_vst_matrix)
#rm(He_full_TOM)

###### CAN RUN BELOW EITHER INTERACTIVE OR AS A SCRIPT ####
## EXPORTING EACH AS SEPARATE MODULES BECAUSE THE INDIVIDUAL NETWORKS ARE TOO LARGE 
# Read in the TOM files as tables, read in matrix files, read in moduleColor files and remove from workspace after exporting modules for each

# Read in the annotation files for all (updated)
load(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/C_gig_C_vir_annotations.RData")

##### Dermo_Tol modules ####
##(already done in bluewaves)
#load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Dermo_Tol_full_TOM.RData") # already done
#load(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Dermo_Tol_full_moduleColors.RData")
#
## TURQUOISE 
#Dermo_Tol_full_modules = "turquoise"
## Select module probes
#Dermo_Tol_full_probes = colnames(Dermo_Tolerant_dds_vst_matrix)
#Dermo_Tol_full_inModule = is.finite(match(Dermo_Tol_full_moduleColors, Dermo_Tol_full_modules))
#Dermo_Tol_full_modProbes = Dermo_Tol_full_probes[Dermo_Tol_full_inModule]
#Dermo_Tol_full_modGenes = C_vir_rtracklayer$ID[match(Dermo_Tol_full_modProbes, C_vir_rtracklayer$ID)]
## Select the corresponding Topological Overlap
#Dermo_Tol_full_modTOM = Dermo_Tol_full_TOM[Dermo_Tol_full_inModule, Dermo_Tol_full_inModule]
#dimnames(Dermo_Tol_full_modTOM) = list(Dermo_Tol_full_modProbes, Dermo_Tol_full_modProbes)
## Export the network into edge and node list files Cytoscape can read
#Dermo_Tol_full_cyt = exportNetworkToCytoscape(Dermo_Tol_full_modTOM,
#                                              edgeFile = paste("CytoscapeInput-edges-Dermo_Tol_full", paste(Dermo_Tol_full_modules, collapse="-"), ".txt", sep=""),
#                                              nodeFile = paste("CytoscapeInput-nodes-Dermo_Tol_full", paste(Dermo_Tol_full_modules, collapse="-"), ".txt", sep=""),
#                                              weighted = TRUE,
#                                              threshold = 0.00, # using 0 threshold so no genes are subset out 
#                                              nodeNames = Dermo_Tol_full_modProbes,
#                                              altNodeNames = Dermo_Tol_full_modGenes,
#                                              nodeAttr = Dermo_Tol_full_moduleColors[Dermo_Tol_full_inModule])
#rm(Dermo_Tol_full_TOM)
#rm(Dermo_Tol_full_moduleColors)

#### Dermo_Sus modules ####
#(already done in bluewaves)
load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Dermo_Sus_full_TOM.RData") # already done
load(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Dermo_Sus_full_moduleColors.RData")
Dermo_Susceptible_dds_vst_matrix <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Dermo_Susceptible_dds_vst_matrix.table")


# TURQUOISE 
Dermo_Sus_full_modules = "lightpink4"
# Select module probes
Dermo_Sus_full_probes = colnames(Dermo_Susceptible_dds_vst_matrix)
Dermo_Sus_full_inModule = is.finite(match(Dermo_Sus_full_moduleColors, Dermo_Sus_full_modules))
Dermo_Sus_full_modProbes = Dermo_Sus_full_probes[Dermo_Sus_full_inModule]
Dermo_Sus_full_modGenes = C_vir_rtracklayer$ID[match(Dermo_Sus_full_modProbes, C_vir_rtracklayer$ID)]
# Select the corresponding Topological Overlap
Dermo_Sus_full_modTOM = Dermo_Sus_full_TOM[Dermo_Sus_full_inModule, Dermo_Sus_full_inModule]
dimnames(Dermo_Sus_full_modTOM) = list(Dermo_Sus_full_modProbes, Dermo_Sus_full_modProbes)
# Export the network into edge and node list files Cytoscape can read
Dermo_Sus_full_cyt = exportNetworkToCytoscape(Dermo_Sus_full_modTOM,
                                              edgeFile = paste("CytoscapeInput-edges-Dermo_Sus_full", paste(Dermo_Sus_full_modules, collapse="-"), ".txt", sep=""),
                                              nodeFile = paste("CytoscapeInput-nodes-Dermo_Sus_full", paste(Dermo_Sus_full_modules, collapse="-"), ".txt", sep=""),
                                              weighted = TRUE,
                                              threshold = 0.00, # using 0 threshold so no genes are subset out 
                                              nodeNames = Dermo_Sus_full_modProbes,
                                              altNodeNames = Dermo_Sus_full_modGenes,
                                              nodeAttr = Dermo_Sus_full_moduleColors[Dermo_Sus_full_inModule])
rm(Dermo_Sus_full_TOM)
rm(Dermo_Sus_full_moduleColors)
rm(Dermo_Susceptible_dds_vst_matrix)

#### Pro_RE22_Pro modules ####
# Pro_RE22_Pro_RI:MEturquoise    , MEdarkslateblue, MEsteelblue    
# Pro_RE22_Pro_S4: MEroyalblue, MEturquoise, MEsteelblue

load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_Pro_full_TOM.RData")
load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_Pro_full_moduleColors.RData")
Pro_RE22_dds_rlog_matrix_Pro <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_dds_rlog_matrix_Pro.table")

## MEskyblue3,  MEwhite
Pro_RE22_Pro_full_modules = "skyblue3"
# Select module probes
Pro_RE22_Pro_full_probes = colnames(Pro_RE22_dds_rlog_matrix_Pro)
Pro_RE22_Pro_full_inModule = is.finite(match(Pro_RE22_Pro_full_moduleColors, Pro_RE22_Pro_full_modules))
Pro_RE22_Pro_full_modProbes = Pro_RE22_Pro_full_probes[Pro_RE22_Pro_full_inModule]
Pro_RE22_Pro_full_modGenes = C_vir_rtracklayer$ID[match(Pro_RE22_Pro_full_modProbes, C_vir_rtracklayer$ID)]
# Select the corresponding Topological Overlap
Pro_RE22_Pro_full_modTOM = Pro_RE22_Pro_full_TOM[Pro_RE22_Pro_full_inModule, Pro_RE22_Pro_full_inModule]
dimnames(Pro_RE22_Pro_full_modTOM ) = list(Pro_RE22_Pro_full_modProbes, Pro_RE22_Pro_full_modProbes)
# Export the network into edge and node list files Cytoscape can read
Pro_RE22_Pro_full_cyt = exportNetworkToCytoscape(Pro_RE22_Pro_full_modTOM,
                                                  edgeFile = paste("CytoscapeInput-edges-Pro_RE22_Pro_full", paste(Pro_RE22_Pro_full_modules, collapse="-"), ".txt", sep=""),
                                                  nodeFile = paste("CytoscapeInput-nodes-Pro_RE22_Pro_full", paste(Pro_RE22_Pro_full_modules, collapse="-"), ".txt", sep=""),
                                                  weighted = TRUE,
                                                  threshold = 0.00, # use zero threshold so no genes are subset out 
                                                  nodeNames = Pro_RE22_Pro_full_modProbes,
                                                  altNodeNames = Pro_RE22_Pro_full_modGenes,
                                                  nodeAttr = Pro_RE22_Pro_full_moduleColors[Pro_RE22_Pro_full_inModule])
## MEwhite
Pro_RE22_Pro_full_modules = "white"
# Select module probes
Pro_RE22_Pro_full_probes = colnames(Pro_RE22_dds_rlog_matrix_Pro)
Pro_RE22_Pro_full_inModule = is.finite(match(Pro_RE22_Pro_full_moduleColors, Pro_RE22_Pro_full_modules))
Pro_RE22_Pro_full_modProbes = Pro_RE22_Pro_full_probes[Pro_RE22_Pro_full_inModule]
Pro_RE22_Pro_full_modGenes = C_vir_rtracklayer$ID[match(Pro_RE22_Pro_full_modProbes, C_vir_rtracklayer$ID)]
# Select the corresponding Topological Overlap
Pro_RE22_Pro_full_modTOM = Pro_RE22_Pro_full_TOM[Pro_RE22_Pro_full_inModule, Pro_RE22_Pro_full_inModule]
dimnames(Pro_RE22_Pro_full_modTOM ) = list(Pro_RE22_Pro_full_modProbes, Pro_RE22_Pro_full_modProbes)
# Export the network into edge and node list files Cytoscape can read
Pro_RE22_Pro_full_cyt = exportNetworkToCytoscape(Pro_RE22_Pro_full_modTOM,
                                                 edgeFile = paste("CytoscapeInput-edges-Pro_RE22_Pro_full", paste(Pro_RE22_Pro_full_modules, collapse="-"), ".txt", sep=""),
                                                 nodeFile = paste("CytoscapeInput-nodes-Pro_RE22_Pro_full", paste(Pro_RE22_Pro_full_modules, collapse="-"), ".txt", sep=""),
                                                 weighted = TRUE,
                                                 threshold = 0.00, # use zero threshold so no genes are subset out 
                                                 nodeNames = Pro_RE22_Pro_full_modProbes,
                                                 altNodeNames = Pro_RE22_Pro_full_modGenes,
                                                 nodeAttr = Pro_RE22_Pro_full_moduleColors[Pro_RE22_Pro_full_inModule])


## TURQUOISE
#Pro_RE22_Pro_full_modules = "turquoise"
## Select module probes
#Pro_RE22_Pro_full_probes = colnames(Pro_RE22_dds_rlog_matrix_Pro)
#Pro_RE22_Pro_full_inModule = is.finite(match(Pro_RE22_Pro_full_moduleColors, Pro_RE22_Pro_full_modules))
#Pro_RE22_Pro_full_modProbes = Pro_RE22_Pro_full_probes[Pro_RE22_Pro_full_inModule]
#Pro_RE22_Pro_full_modGenes = C_vir_rtracklayer$ID[match(Pro_RE22_Pro_full_modProbes, C_vir_rtracklayer$ID)]
## Select the corresponding Topological Overlap
#Pro_RE22_Pro_full_modTOM = Pro_RE22_Pro_full_TOM[Pro_RE22_Pro_full_inModule, Pro_RE22_Pro_full_inModule]
#dimnames(Pro_RE22_Pro_full_modTOM ) = list(Pro_RE22_Pro_full_modProbes, Pro_RE22_Pro_full_modProbes)
## Export the network into edge and node list files Cytoscape can read
#Pro_RE22_Pro_full_cyt = exportNetworkToCytoscape(Pro_RE22_Pro_full_modTOM,
#                                                  edgeFile = paste("CytoscapeInput-edges-Pro_RE22_Pro_full", paste(Pro_RE22_Pro_full_modules, collapse="-"), ".txt", sep=""),
#                                                  nodeFile = paste("CytoscapeInput-nodes-Pro_RE22_Pro_full", paste(Pro_RE22_Pro_full_modules, collapse="-"), ".txt", sep=""),
#                                                  weighted = TRUE,
#                                                  threshold = 0.00, # use zero threshold so no genes are subset out 
#                                                  nodeNames = Pro_RE22_Pro_full_modProbes,
#                                                  altNodeNames = Pro_RE22_Pro_full_modGenes,
#                                                  nodeAttr = Pro_RE22_Pro_full_moduleColors[Pro_RE22_Pro_full_inModule])
#
# darkslateblue
#Pro_RE22_Pro_full_modules = "darkslateblue"
## Select module probes
#Pro_RE22_Pro_full_probes = colnames(Pro_RE22_dds_rlog_matrix_Pro)
#Pro_RE22_Pro_full_inModule = is.finite(match(Pro_RE22_Pro_full_moduleColors, Pro_RE22_Pro_full_modules))
#Pro_RE22_Pro_full_modProbes = Pro_RE22_Pro_full_probes[Pro_RE22_Pro_full_inModule]
#Pro_RE22_Pro_full_modGenes = C_vir_rtracklayer$ID[match(Pro_RE22_Pro_full_modProbes, C_vir_rtracklayer$ID)]
## Select the corresponding Topological Overlap
#Pro_RE22_Pro_full_modTOM = Pro_RE22_Pro_full_TOM[Pro_RE22_Pro_full_inModule, Pro_RE22_Pro_full_inModule]
#dimnames(Pro_RE22_Pro_full_modTOM ) = list(Pro_RE22_Pro_full_modProbes, Pro_RE22_Pro_full_modProbes)
## Export the network into edge and node list files Cytoscape can read
#Pro_RE22_Pro_full_cyt = exportNetworkToCytoscape(Pro_RE22_Pro_full_modTOM,
#                                                 edgeFile = paste("CytoscapeInput-edges-Pro_RE22_Pro_full", paste(Pro_RE22_Pro_full_modules, collapse="-"), ".txt", sep=""),
#                                                 nodeFile = paste("CytoscapeInput-nodes-Pro_RE22_Pro_full", paste(Pro_RE22_Pro_full_modules, collapse="-"), ".txt", sep=""),
#                                                 weighted = TRUE,
#                                                 threshold = 0.00, # use zero threshold so no genes are subset out 
#                                                 nodeNames = Pro_RE22_Pro_full_modProbes,
#                                                 altNodeNames = Pro_RE22_Pro_full_modGenes,
#                                                 nodeAttr = Pro_RE22_Pro_full_moduleColors[Pro_RE22_Pro_full_inModule])# steelblue 
#Pro_RE22_Pro_full_modules = "steelblue"
## Select module probes
#Pro_RE22_Pro_full_probes = colnames(Pro_RE22_dds_rlog_matrix_Pro)
#Pro_RE22_Pro_full_inModule = is.finite(match(Pro_RE22_Pro_full_moduleColors, Pro_RE22_Pro_full_modules))
#Pro_RE22_Pro_full_modProbes = Pro_RE22_Pro_full_probes[Pro_RE22_Pro_full_inModule]
#Pro_RE22_Pro_full_modGenes = C_vir_rtracklayer$ID[match(Pro_RE22_Pro_full_modProbes, C_vir_rtracklayer$ID)]
## Select the corresponding Topological Overlap
#Pro_RE22_Pro_full_modTOM = Pro_RE22_Pro_full_TOM[Pro_RE22_Pro_full_inModule, Pro_RE22_Pro_full_inModule]
#dimnames(Pro_RE22_Pro_full_modTOM ) = list(Pro_RE22_Pro_full_modProbes, Pro_RE22_Pro_full_modProbes)
## Export the network into edge and node list files Cytoscape can read
#Pro_RE22_Pro_full_cyt = exportNetworkToCytoscape(Pro_RE22_Pro_full_modTOM,
#                                                 edgeFile = paste("CytoscapeInput-edges-Pro_RE22_Pro_full", paste(Pro_RE22_Pro_full_modules, collapse="-"), ".txt", sep=""),
#                                                 nodeFile = paste("CytoscapeInput-nodes-Pro_RE22_Pro_full", paste(Pro_RE22_Pro_full_modules, collapse="-"), ".txt", sep=""),
#                                                 weighted = TRUE,
#                                                 threshold = 0.00, # use zero threshold so no genes are subset out 
#                                                 nodeNames = Pro_RE22_Pro_full_modProbes,
#                                                 altNodeNames = Pro_RE22_Pro_full_modGenes,
#                                                 nodeAttr = Pro_RE22_Pro_full_moduleColors[Pro_RE22_Pro_full_inModule])# royalblue
#Pro_RE22_Pro_full_modules = "royalblue"
## Select module probes
#Pro_RE22_Pro_full_probes = colnames(Pro_RE22_dds_rlog_matrix_Pro)
#Pro_RE22_Pro_full_inModule = is.finite(match(Pro_RE22_Pro_full_moduleColors, Pro_RE22_Pro_full_modules))
#Pro_RE22_Pro_full_modProbes = Pro_RE22_Pro_full_probes[Pro_RE22_Pro_full_inModule]
#Pro_RE22_Pro_full_modGenes = C_vir_rtracklayer$ID[match(Pro_RE22_Pro_full_modProbes, C_vir_rtracklayer$ID)]
## Select the corresponding Topological Overlap
#Pro_RE22_Pro_full_modTOM = Pro_RE22_Pro_full_TOM[Pro_RE22_Pro_full_inModule, Pro_RE22_Pro_full_inModule]
#dimnames(Pro_RE22_Pro_full_modTOM ) = list(Pro_RE22_Pro_full_modProbes, Pro_RE22_Pro_full_modProbes)
## Export the network into edge and node list files Cytoscape can read
#Pro_RE22_Pro_full_cyt = exportNetworkToCytoscape(Pro_RE22_Pro_full_modTOM,
#                                                 edgeFile = paste("CytoscapeInput-edges-Pro_RE22_Pro_full", paste(Pro_RE22_Pro_full_modules, collapse="-"), ".txt", sep=""),
#                                                 nodeFile = paste("CytoscapeInput-nodes-Pro_RE22_Pro_full", paste(Pro_RE22_Pro_full_modules, collapse="-"), ".txt", sep=""),
#                                                 weighted = TRUE,
#                                                 threshold = 0.00, # use zero threshold so no genes are subset out 
#                                                 nodeNames = Pro_RE22_Pro_full_modProbes,
#                                                 altNodeNames = Pro_RE22_Pro_full_modGenes,
#                                                 nodeAttr = Pro_RE22_Pro_full_moduleColors[Pro_RE22_Pro_full_inModule])
rm(Pro_RE22_Pro_full_TOM)
rm(Pro_RE22_Pro_full_moduleColors)
rm(Pro_RE22_dds_rlog_matrix_Pro)
#
##### ZHANG ####
## Zhang: MEblack
#load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Zhang_full_TOM.RData")
#Zhang_dds_rlog_matrix <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Zhang_dds_rlog_matrix.table")
#load(file=  "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Zhang_full_moduleColors.RData") 
#
#Zhang_full_modules = "black"
## Select module probes
#Zhang_full_probes = colnames(Zhang_dds_rlog_matrix)
#Zhang_full_inModule = is.finite(match(Zhang_full_moduleColors, Zhang_full_modules))
#Zhang_full_modProbes = Zhang_full_probes[Zhang_full_inModule]
#Zhang_full_modGenes = C_vir_rtracklayer$ID[match(Zhang_full_modProbes, C_vir_rtracklayer$ID)]
## Select the corresponding Topological Overlap
#Zhang_full_modTOM = Zhang_full_TOM[Zhang_full_inModule,  Zhang_full_inModule]
#dimnames(Zhang_full_modTOM ) = list(Zhang_full_modProbes, Zhang_full_modProbes)
## Export the network into edge and node list files Cytoscape can read
#Zhang_full_cyt = exportNetworkToCytoscape(Zhang_full_modTOM,
#                                          edgeFile = paste("CytoscapeInput-edges-Zhang_full", paste(Zhang_full_modules, collapse="-"), ".txt", sep=""),
#                                          nodeFile = paste("CytoscapeInput-nodes-Zhang_full", paste(Zhang_full_modules, collapse="-"), ".txt", sep=""),
#                                          weighted = TRUE,
#                                          threshold = 0.00, # use zero threshold so no genes are subset out 
#                                          nodeNames = Zhang_full_modProbes,
#                                          altNodeNames = Zhang_full_modGenes,
#                                          nodeAttr = Zhang_full_moduleColors[Zhang_full_inModule])
#
#rm(Zhang_full_TOM)
#rm(Zhang_dds_rlog_matrix)
#rm(Zhang_full_moduleColors)
#
##### RUBIO ####
## Rubio: MEmagenta, MEturquoise, MEblue, MEbrown, MEblack   
load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Rubio_full_TOM.RData")
Rubio_dds_rlog_matrix <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Rubio_dds_rlog_matrix.table")
load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Rubio_full_moduleColors.RData")

Rubio_full_modules = "black"
# Select module probes
Rubio_full_probes = colnames(Rubio_dds_rlog_matrix)
Rubio_full_inModule = is.finite(match(Rubio_full_moduleColors, Rubio_full_modules))
Rubio_full_modProbes = Rubio_full_probes[Rubio_full_inModule]
Rubio_full_modGenes = C_vir_rtracklayer$ID[match(Rubio_full_modProbes, C_vir_rtracklayer$ID)]
# Select the corresponding Topological Overlap
Rubio_full_modTOM = Rubio_full_TOM[Rubio_full_inModule,  Rubio_full_inModule]
dimnames(Rubio_full_modTOM ) = list(Rubio_full_modProbes, Rubio_full_modProbes)
# Export the network into edge and node list files Cytoscape can read
Rubio_full_cyt = exportNetworkToCytoscape(Rubio_full_modTOM,
                                               edgeFile = paste("CytoscapeInput-edges-Rubio_full", paste(Rubio_full_modules, collapse="-"), ".txt", sep=""),
                                               nodeFile = paste("CytoscapeInput-nodes-Rubio_full", paste(Rubio_full_modules, collapse="-"), ".txt", sep=""),
                                               weighted = TRUE,
                                               threshold = 0.00, # use zero threshold so no genes are subset out 
                                               nodeNames = Rubio_full_modProbes,
                                               altNodeNames = Rubio_full_modGenes,
                                               nodeAttr = Rubio_full_moduleColors[Rubio_full_inModule])



#Rubio_full_modules = "magenta"
## Select module probes
#Rubio_full_probes = colnames(Rubio_dds_rlog_matrix)
#Rubio_full_inModule = is.finite(match(Rubio_full_moduleColors, Rubio_full_modules))
#Rubio_full_modProbes = Rubio_full_probes[Rubio_full_inModule]
#Rubio_full_modGenes = C_vir_rtracklayer$ID[match(Rubio_full_modProbes, C_vir_rtracklayer$ID)]
## Select the corresponding Topological Overlap
#Rubio_full_modTOM = Rubio_full_TOM[Rubio_full_inModule,  Rubio_full_inModule]
#dimnames(Rubio_full_modTOM ) = list(Rubio_full_modProbes, Rubio_full_modProbes)
## Export the network into edge and node list files Cytoscape can read
#Rubio_full_cyt = exportNetworkToCytoscape(Rubio_full_modTOM,
#                                               edgeFile = paste("CytoscapeInput-edges-Rubio_full", paste(Rubio_full_modules, collapse="-"), ".txt", sep=""),
#                                               nodeFile = paste("CytoscapeInput-nodes-Rubio_full", paste(Rubio_full_modules, collapse="-"), ".txt", sep=""),
#                                               weighted = TRUE,
#                                               threshold = 0.00, # use zero threshold so no genes are subset out 
#                                               nodeNames = Rubio_full_modProbes,
#                                               altNodeNames = Rubio_full_modGenes,
#                                               nodeAttr = Rubio_full_moduleColors[Rubio_full_inModule])
#
#Rubio_full_modules = "turquoise"
## Select module probes
#Rubio_full_probes = colnames(Rubio_dds_rlog_matrix)
#Rubio_full_inModule = is.finite(match(Rubio_full_moduleColors, Rubio_full_modules))
#Rubio_full_modProbes = Rubio_full_probes[Rubio_full_inModule]
#Rubio_full_modGenes = C_vir_rtracklayer$ID[match(Rubio_full_modProbes, C_vir_rtracklayer$ID)]
## Select the corresponding Topological Overlap
#Rubio_full_modTOM = Rubio_full_TOM[Rubio_full_inModule,  Rubio_full_inModule]
#dimnames(Rubio_full_modTOM ) = list(Rubio_full_modProbes, Rubio_full_modProbes)
## Export the network into edge and node list files Cytoscape can read
#Rubio_full_cyt = exportNetworkToCytoscape(Rubio_full_modTOM,
#                                          edgeFile = paste("CytoscapeInput-edges-Rubio_full", paste(Rubio_full_modules, collapse="-"), ".txt", sep=""),
#                                          nodeFile = paste("CytoscapeInput-nodes-Rubio_full", paste(Rubio_full_modules, collapse="-"), ".txt", sep=""),
#                                          weighted = TRUE,
#                                          threshold = 0.00, # use zero threshold so no genes are subset out 
#                                          nodeNames = Rubio_full_modProbes,
#                                          altNodeNames = Rubio_full_modGenes,
#                                          nodeAttr = Rubio_full_moduleColors[Rubio_full_inModule])
#
#Rubio_full_modules = "blue"
## Select module probes
#Rubio_full_probes = colnames(Rubio_dds_rlog_matrix)
#Rubio_full_inModule = is.finite(match(Rubio_full_moduleColors, Rubio_full_modules))
#Rubio_full_modProbes = Rubio_full_probes[Rubio_full_inModule]
#Rubio_full_modGenes = C_vir_rtracklayer$ID[match(Rubio_full_modProbes, C_vir_rtracklayer$ID)]
## Select the corresponding Topological Overlap
#Rubio_full_modTOM = Rubio_full_TOM[Rubio_full_inModule,  Rubio_full_inModule]
#dimnames(Rubio_full_modTOM ) = list(Rubio_full_modProbes, Rubio_full_modProbes)
## Export the network into edge and node list files Cytoscape can read
#Rubio_full_cyt = exportNetworkToCytoscape(Rubio_full_modTOM,
#                                          edgeFile = paste("CytoscapeInput-edges-Rubio_full", paste(Rubio_full_modules, collapse="-"), ".txt", sep=""),
#                                          nodeFile = paste("CytoscapeInput-nodes-Rubio_full", paste(Rubio_full_modules, collapse="-"), ".txt", sep=""),
#                                          weighted = TRUE,
#                                          threshold = 0.00, # use zero threshold so no genes are subset out 
#                                          nodeNames = Rubio_full_modProbes,
#                                          altNodeNames = Rubio_full_modGenes,
#                                          nodeAttr = Rubio_full_moduleColors[Rubio_full_inModule])
#
#Rubio_full_modules = "brown"
## Select module probes
#Rubio_full_probes = colnames(Rubio_dds_rlog_matrix)
#Rubio_full_inModule = is.finite(match(Rubio_full_moduleColors, Rubio_full_modules))
#Rubio_full_modProbes = Rubio_full_probes[Rubio_full_inModule]
#Rubio_full_modGenes = C_vir_rtracklayer$ID[match(Rubio_full_modProbes, C_vir_rtracklayer$ID)]
## Select the corresponding Topological Overlap
#Rubio_full_modTOM = Rubio_full_TOM[Rubio_full_inModule,  Rubio_full_inModule]
#dimnames(Rubio_full_modTOM ) = list(Rubio_full_modProbes, Rubio_full_modProbes)
## Export the network into edge and node list files Cytoscape can read
#Rubio_full_cyt = exportNetworkToCytoscape(Rubio_full_modTOM,
#                                          edgeFile = paste("CytoscapeInput-edges-Rubio_full", paste(Rubio_full_modules, collapse="-"), ".txt", sep=""),
#                                          nodeFile = paste("CytoscapeInput-nodes-Rubio_full", paste(Rubio_full_modules, collapse="-"), ".txt", sep=""),
#                                          weighted = TRUE,
#                                          threshold = 0.00, # use zero threshold so no genes are subset out 
#                                          nodeNames = Rubio_full_modProbes,
#                                          altNodeNames = Rubio_full_modGenes,
#                                          nodeAttr = Rubio_full_moduleColors[Rubio_full_inModule])
#
rm(Rubio_full_TOM)
rm(Rubio_dds_rlog_matrix)
rm(Rubio_full_moduleColors)
#
##### DELORG RES ####
## deLorg_Res: MEturquoise
#load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/deLorg_Res_full_TOM.RData")
#load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/deLorg_Res_full_moduleColors.RData")
#deLorgeril_Resistant_dds_vst_matrix <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/deLorgeril_Resistant_dds_vst_matrix.table")
#
#deLorg_Res_full_modules = "turquoise"
## Select module probes
#deLorg_Res_full_probes = colnames(deLorgeril_Resistant_dds_vst_matrix)
#deLorg_Res_full_inModule = is.finite(match(deLorg_Res_full_moduleColors, deLorg_Res_full_modules))
#deLorg_Res_full_modProbes = deLorg_Res_full_probes[deLorg_Res_full_inModule]
#deLorg_Res_full_modGenes = C_vir_rtracklayer$ID[match(deLorg_Res_full_modProbes, C_vir_rtracklayer$ID)]
## Select the corresponding Topological Overlap
#deLorg_Res_full_modTOM = deLorg_Res_full_TOM[deLorg_Res_full_inModule,  deLorg_Res_full_inModule]
#dimnames(deLorg_Res_full_modTOM ) = list(deLorg_Res_full_modProbes, deLorg_Res_full_modProbes)
## Export the network into edge and node list files Cytoscape can read
#deLorg_Res_full_cyt = exportNetworkToCytoscape(deLorg_Res_full_modTOM,
#                                                 edgeFile = paste("CytoscapeInput-edges-deLorg_Res_full", paste(deLorg_Res_full_modules, collapse="-"), ".txt", sep=""),
#                                                 nodeFile = paste("CytoscapeInput-nodes-deLorg_Res_full", paste(deLorg_Res_full_modules, collapse="-"), ".txt", sep=""),
#                                                 weighted = TRUE,
#                                                 threshold = 0.00, # use zero threshold so no genes are subset out 
#                                                 nodeNames = deLorg_Res_full_modProbes,
#                                                 altNodeNames = deLorg_Res_full_modGenes,
#                                                 nodeAttr = deLorg_Res_full_moduleColors[deLorg_Res_full_inModule])
#
#rm(deLorg_Res_full_TOM)
#rm(deLorg_Res_full_moduleColors)
#rm(deLorgeril_Resistant_dds_vst_matrix)
#
##### DELORG SUS ####
## deLorg_Sus: MEturquoise
#load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/deLorg_Sus_full_TOM.RData")
#deLorgeril_Susceptible_dds_vst_matrix <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/deLorgeril_Susceptible_dds_vst_matrix.table")
#load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/deLorg_Sus_full_moduleColors.RData")
#
#deLorg_Sus_full_modules = "turquoise"
## Select module probes
#deLorg_Sus_full_probes = colnames(deLorgeril_Susceptible_dds_vst_matrix)
#deLorg_Sus_full_inModule = is.finite(match(deLorg_Sus_full_moduleColors, deLorg_Sus_full_modules))
#deLorg_Sus_full_modProbes = deLorg_Sus_full_probes[deLorg_Sus_full_inModule]
#deLorg_Sus_full_modGenes = C_vir_rtracklayer$ID[match(deLorg_Sus_full_modProbes, C_vir_rtracklayer$ID)]
## Select the corresponding Topological Overlap
#deLorg_Sus_full_modTOM = deLorg_Sus_full_TOM[deLorg_Sus_full_inModule,  deLorg_Sus_full_inModule]
#dimnames(deLorg_Sus_full_modTOM ) = list(deLorg_Sus_full_modProbes, deLorg_Sus_full_modProbes)
## Export the network into edge and node list files Cytoscape can read
#deLorg_Sus_full_cyt = exportNetworkToCytoscape(deLorg_Res_full_modTOM,
#                                               edgeFile = paste("CytoscapeInput-edges-deLorg_Sus_full", paste(deLorg_Sus_full_modules, collapse="-"), ".txt", sep=""),
#                                               nodeFile = paste("CytoscapeInput-nodes-deLorg_Sus_full", paste(deLorg_Sus_full_modules, collapse="-"), ".txt", sep=""),
#                                               weighted = TRUE,
#                                               threshold = 0.00, # use zero threshold so no genes are subset out 
#                                               nodeNames = deLorg_Sus_full_modProbes,
#                                               altNodeNames = deLorg_Sus_full_modGenes,
#                                               nodeAttr = deLorg_Sus_full_moduleColors[deLorg_Sus_full_inModule])
#
#rm(deLorg_Sus_full_TOM)
#rm(deLorgeril_Susceptible_dds_vst_matrix)
#rm(deLorg_Sus_full_moduleColors)
#
##### HE ####
## He:  MEpurple, MEyellow
#load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/He_full_TOM.RData")
#He_dds_vst_matrix <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/He_dds_vst_matrix.table")
#load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/He_full_moduleColors.RData")
#
#He_full_modules = "purple"
## Select module probes
#He_full_probes = colnames(He_dds_vst_matrix)
#He_full_inModule = is.finite(match(He_full_moduleColors, He_full_modules))
#He_full_modProbes = He_full_probes[He_full_inModule]
#He_full_modGenes = C_vir_rtracklayer$ID[match(He_full_modProbes, C_vir_rtracklayer$ID)]
## Select the corresponding Topological Overlap
#He_full_modTOM = He_full_TOM[He_full_inModule,  He_full_inModule]
#dimnames(He_full_modTOM ) = list(He_full_modProbes, He_full_modProbes)
## Export the network into edge and node list files Cytoscape can read
#He_full_cyt = exportNetworkToCytoscape(He_full_modTOM,
#                                               edgeFile = paste("CytoscapeInput-edges-He_full", paste(He_full_modules, collapse="-"), ".txt", sep=""),
#                                               nodeFile = paste("CytoscapeInput-nodes-He_full", paste(He_full_modules, collapse="-"), ".txt", sep=""),
#                                               weighted = TRUE,
#                                               threshold = 0.00, # use zero threshold so no genes are subset out 
#                                               nodeNames = He_full_modProbes,
#                                               altNodeNames = He_full_modGenes,
#                                               nodeAttr = He_full_moduleColors[He_full_inModule])
#
#He_full_modules = "yellow"
## Select module probes
#He_full_probes = colnames(He_dds_vst_matrix)
#He_full_inModule = is.finite(match(He_full_moduleColors, He_full_modules))
#He_full_modProbes = He_full_probes[He_full_inModule]
#He_full_modGenes = C_vir_rtracklayer$ID[match(He_full_modProbes, C_vir_rtracklayer$ID)]
## Select the corresponding Topological Overlap
#He_full_modTOM = He_full_TOM[He_full_inModule,  He_full_inModule]
#dimnames(He_full_modTOM ) = list(He_full_modProbes, He_full_modProbes)
## Export the network into edge and node list files Cytoscape can read
#He_full_cyt = exportNetworkToCytoscape(He_full_modTOM,
#                                       edgeFile = paste("CytoscapeInput-edges-He_full", paste(He_full_modules, collapse="-"), ".txt", sep=""),
#                                       nodeFile = paste("CytoscapeInput-nodes-He_full", paste(He_full_modules, collapse="-"), ".txt", sep=""),
#                                       weighted = TRUE,
#                                       threshold = 0.00, # use zero threshold so no genes are subset out 
#                                       nodeNames = He_full_modProbes,
#                                       altNodeNames = He_full_modGenes,
#                                       nodeAttr = He_full_moduleColors[He_full_inModule])
#rm(He_full_TOM)
#rm(He_dds_vst_matrix)
#rm(He_full_moduleColors)
#
sessionInfo()
