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
## Which modules? 
# Which modules? 
# deLorg_Res: MEturquoise
# deLorg_Sus: MEturquoise
# Dermo: MEturquoise
# Pro_RE22_Pro_RI:MEturquoise    , MEdarkslateblue, MEsteelblue    
# Pro_RE22_Pro_S4: MEroyalblue, MEturquoise, MEsteelblue
# He:  MEpurple, MEyellow
# Rubio: MEmagenta, MEturquoise, MEblue, MEbrown    
# Zhang: MEblack

## Load saved matrices one by one and then unload to free up memory 

#### CALCULATE TOM ####

#Dermo_Tolerant_dds_vst_matrix <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Dermo_Tolerant_dds_vst_matrix.table")
#Dermo_Tol_full_TOM = TOMsimilarityFromExpr(Dermo_Tolerant_dds_vst_matrix, power = 3, TOMType = "signed", networkType= "signed hybrid", corType = "bicor") 
#save(Dermo_Tol_full_TOM, file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Dermo_Tol_full_TOM.RData")

# power taken from original WGCNA run 

# C. virginica experiments
#Pro_RE22_dds_rlog_matrix_Pro <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_dds_rlog_matrix_Pro.table")
#Pro_RE22_Pro_full_TOM = TOMsimilarityFromExpr(Pro_RE22_dds_rlog_matrix_Pro, power = 4, TOMType = "signed", networkType= "signed hybrid", corType = "bicor") 
#save(Pro_RE22_Pro_full_TOM , file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_Pro_full_TOM.RData")
#rm(Pro_RE22_Pro_full_TOM) # remove from workspace to free up memory
#rm(Pro_RE22_dds_rlog_matrix_Pro)
#
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

#### Dermo_Tol modules ####
#(already done in bluewaves)
# load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Dermo_Tol_full_TOM.RData") # already done
#load(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Dermo_Tol_full_moduleColors.RData")

# TURQUOISE 
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

#### Pro_RE22_Pro modules ####
# Pro_RE22_Pro_RI:MEturquoise    , MEdarkslateblue, MEsteelblue    
# Pro_RE22_Pro_S4: MEroyalblue, MEturquoise, MEsteelblue

load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_Pro_full_TOM.RData")
load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_Pro_full_moduleColors.RData")
Pro_RE22_dds_rlog_matrix_Pro <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_dds_rlog_matrix_Pro.table")

# TURQUOISE
Pro_RE22_Pro_full_modules = "turquoise"
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

# darkslateblue
Pro_RE22_Pro_full_modules = "darkslateblue"
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
                                                 nodeAttr = Pro_RE22_Pro_full_moduleColors[Pro_RE22_Pro_full_inModule])# steelblue 
Pro_RE22_Pro_full_modules = "steelblue"
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
                                                 nodeAttr = Pro_RE22_Pro_full_moduleColors[Pro_RE22_Pro_full_inModule])# royalblue
Pro_RE22_Pro_full_modules = "royalblue"
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
rm(Pro_RE22_Pro_full_TOM)
rm(Pro_RE22_Pro_full_moduleColors)
rm(Pro_RE22_dds_rlog_matrix_Pro)

#### ZHANG ####
# Zhang: MEblack
load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Zhang_full_TOM.RData")
Zhang_dds_rlog_matrix <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Zhang_dds_rlog_matrix.table")
load(file=  "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Zhang_full_moduleColors.RData") 

Zhang_full_modules = "black"
# Select module probes
Zhang_full_probes = colnames(Zhang_dds_rlog_matrix)
Zhang_full_inModule = is.finite(match(Zhang_full_moduleColors, Zhang_full_modules))
Zhang_full_modProbes = Zhang_full_probes[Zhang_full_inModule]
Zhang_full_modGenes = C_vir_rtracklayer$ID[match(Zhang_full_modProbes, C_vir_rtracklayer$ID)]
# Select the corresponding Topological Overlap
Zhang_full_modTOM = Zhang_full_TOM[Zhang_full_inModule,  Zhang_full_inModule]
dimnames(Zhang_full_modTOM ) = list(Zhang_full_modProbes, Zhang_full_modProbes)
# Export the network into edge and node list files Cytoscape can read
Zhang_full_cyt = exportNetworkToCytoscape(Zhang_full_modTOM,
                                          edgeFile = paste("CytoscapeInput-edges-Zhang_full", paste(Zhang_full_modules, collapse="-"), ".txt", sep=""),
                                          nodeFile = paste("CytoscapeInput-nodes-Zhang_full", paste(Zhang_full_modules, collapse="-"), ".txt", sep=""),
                                          weighted = TRUE,
                                          threshold = 0.00, # use zero threshold so no genes are subset out 
                                          nodeNames = Zhang_full_modProbes,
                                          altNodeNames = Zhang_full_modGenes,
                                          nodeAttr = Zhang_full_moduleColors[Zhang_full_inModule])

rm(Zhang_full_TOM)
rm(Zhang_dds_rlog_matrix)
rm(Zhang_full_moduleColors)

#### RUBIO ####
# Rubio: MEmagenta, MEturquoise, MEblue, MEbrown    
load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Rubio_full_TOM.RData")
Rubio_dds_rlog_matrix <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Rubio_dds_rlog_matrix.table")
load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Rubio_full_moduleColors.RData")

Rubio_full_modules = "magenta"
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

Rubio_full_modules = "turquoise"
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

Rubio_full_modules = "blue"
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

Rubio_full_modules = "brown"
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

rm(Rubio_full_TOM)
rm(Rubio_dds_rlog_matrix)
rm(Rubio_full_moduleColors)

#### DELORG RES ####
# deLorg_Res: MEturquoise
load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/deLorg_Res_full_TOM.RData")
load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/deLorg_Res_full_moduleColors.RData")
deLorgeril_Resistant_dds_vst_matrix <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/deLorgeril_Resistant_dds_vst_matrix.table")

deLorg_Res_full_modules = "turquoise"
# Select module probes
deLorg_Res_full_probes = colnames(deLorgeril_Resistant_dds_vst_matrix)
deLorg_Res_full_inModule = is.finite(match(deLorg_Res_full_moduleColors, deLorg_Res_full_modules))
deLorg_Res_full_modProbes = deLorg_Res_full_probes[deLorg_Res_full_inModule]
deLorg_Res_full_modGenes = C_vir_rtracklayer$ID[match(deLorg_Res_full_modProbes, C_vir_rtracklayer$ID)]
# Select the corresponding Topological Overlap
deLorg_Res_full_modTOM = deLorg_Res_full_TOM[deLorg_Res_full_inModule,  deLorg_Res_full_inModule]
dimnames(deLorg_Res_full_modTOM ) = list(deLorg_Res_full_modProbes, deLorg_Res_full_modProbes)
# Export the network into edge and node list files Cytoscape can read
deLorg_Res_full_cyt = exportNetworkToCytoscape(deLorg_Res_full_modTOM,
                                                 edgeFile = paste("CytoscapeInput-edges-deLorg_Res_full", paste(deLorg_Res_full_modules, collapse="-"), ".txt", sep=""),
                                                 nodeFile = paste("CytoscapeInput-nodes-deLorg_Res_full", paste(deLorg_Res_full_modules, collapse="-"), ".txt", sep=""),
                                                 weighted = TRUE,
                                                 threshold = 0.00, # use zero threshold so no genes are subset out 
                                                 nodeNames = deLorg_Res_full_modProbes,
                                                 altNodeNames = deLorg_Res_full_modGenes,
                                                 nodeAttr = deLorg_Res_full_moduleColors[deLorg_Res_full_inModule])

rm(deLorg_Res_full_TOM)
rm(deLorg_Res_full_moduleColors)
rm(deLorgeril_Resistant_dds_vst_matrix)

#### DELORG SUS ####
# deLorg_Sus: MEturquoise
load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/deLorg_Sus_full_TOM.RData")
deLorgeril_Susceptible_dds_vst_matrix <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/deLorgeril_Susceptible_dds_vst_matrix.table")
load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/deLorg_Sus_full_moduleColors.RData")

deLorg_Sus_full_modules = "turquoise"
# Select module probes
deLorg_Sus_full_probes = colnames(deLorgeril_Susceptible_dds_vst_matrix)
deLorg_Sus_full_inModule = is.finite(match(deLorg_Sus_full_moduleColors, deLorg_Sus_full_modules))
deLorg_Sus_full_modProbes = deLorg_Sus_full_probes[deLorg_Sus_full_inModule]
deLorg_Sus_full_modGenes = C_vir_rtracklayer$ID[match(deLorg_Sus_full_modProbes, C_vir_rtracklayer$ID)]
# Select the corresponding Topological Overlap
deLorg_Sus_full_modTOM = deLorg_Sus_full_TOM[deLorg_Sus_full_inModule,  deLorg_Sus_full_inModule]
dimnames(deLorg_Sus_full_modTOM ) = list(deLorg_Sus_full_modProbes, deLorg_Sus_full_modProbes)
# Export the network into edge and node list files Cytoscape can read
deLorg_Sus_full_cyt = exportNetworkToCytoscape(deLorg_Res_full_modTOM,
                                               edgeFile = paste("CytoscapeInput-edges-deLorg_Sus_full", paste(deLorg_Sus_full_modules, collapse="-"), ".txt", sep=""),
                                               nodeFile = paste("CytoscapeInput-nodes-deLorg_Sus_full", paste(deLorg_Sus_full_modules, collapse="-"), ".txt", sep=""),
                                               weighted = TRUE,
                                               threshold = 0.00, # use zero threshold so no genes are subset out 
                                               nodeNames = deLorg_Sus_full_modProbes,
                                               altNodeNames = deLorg_Sus_full_modGenes,
                                               nodeAttr = deLorg_Sus_full_moduleColors[deLorg_Sus_full_inModule])

rm(deLorg_Sus_full_TOM)
rm(deLorgeril_Susceptible_dds_vst_matrix)
rm(deLorg_Sus_full_moduleColors)

#### HE ####
# He:  MEpurple, MEyellow
load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/He_full_TOM.RData")
He_dds_vst_matrix <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/He_dds_vst_matrix.table")
load(file = "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/He_full_moduleColors.RData")

He_full_modules = "purple"
# Select module probes
He_full_probes = colnames(He_dds_vst_matrix)
He_full_inModule = is.finite(match(He_full_moduleColors, He_full_modules))
He_full_modProbes = He_full_probes[He_full_inModule]
He_full_modGenes = C_vir_rtracklayer$ID[match(He_full_modProbes, C_vir_rtracklayer$ID)]
# Select the corresponding Topological Overlap
He_full_modTOM = He_full_TOM[He_full_inModule,  He_full_inModule]
dimnames(He_full_modTOM ) = list(He_full_modProbes, He_full_modProbes)
# Export the network into edge and node list files Cytoscape can read
He_full_cyt = exportNetworkToCytoscape(He_full_modTOM,
                                               edgeFile = paste("CytoscapeInput-edges-He_full", paste(He_full_modules, collapse="-"), ".txt", sep=""),
                                               nodeFile = paste("CytoscapeInput-nodes-He_full", paste(He_full_modules, collapse="-"), ".txt", sep=""),
                                               weighted = TRUE,
                                               threshold = 0.00, # use zero threshold so no genes are subset out 
                                               nodeNames = He_full_modProbes,
                                               altNodeNames = He_full_modGenes,
                                               nodeAttr = He_full_moduleColors[He_full_inModule])

He_full_modules = "yellow"
# Select module probes
He_full_probes = colnames(He_dds_vst_matrix)
He_full_inModule = is.finite(match(He_full_moduleColors, He_full_modules))
He_full_modProbes = He_full_probes[He_full_inModule]
He_full_modGenes = C_vir_rtracklayer$ID[match(He_full_modProbes, C_vir_rtracklayer$ID)]
# Select the corresponding Topological Overlap
He_full_modTOM = He_full_TOM[He_full_inModule,  He_full_inModule]
dimnames(He_full_modTOM ) = list(He_full_modProbes, He_full_modProbes)
# Export the network into edge and node list files Cytoscape can read
He_full_cyt = exportNetworkToCytoscape(He_full_modTOM,
                                       edgeFile = paste("CytoscapeInput-edges-He_full", paste(He_full_modules, collapse="-"), ".txt", sep=""),
                                       nodeFile = paste("CytoscapeInput-nodes-He_full", paste(He_full_modules, collapse="-"), ".txt", sep=""),
                                       weighted = TRUE,
                                       threshold = 0.00, # use zero threshold so no genes are subset out 
                                       nodeNames = He_full_modProbes,
                                       altNodeNames = He_full_modGenes,
                                       nodeAttr = He_full_moduleColors[He_full_inModule])
rm(He_full_TOM)
rm(He_dds_vst_matrix)
rm(He_full_moduleColors)

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