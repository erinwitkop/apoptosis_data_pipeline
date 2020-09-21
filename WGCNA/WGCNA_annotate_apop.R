# September 21st, 2020
# annotate apoptosis transcripts for all WGCNA subset modules in bluewaves

# Clear global workspace
rm(list = ls())

# Load libraries
library(tidyverse)

# Read in the annotation files for all (updated)
load(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/C_gig_C_vir_annotations.RData")

# read in the apoptosis text files
#Dermo_Tol_fullturquoise_apop_hits <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Dermo_Tol_fullturquoise_apop_hits.txt", col.names = c("fromNode","toNode","weight","direction", "fromAltName","toAltName"))
Dermo_Sus_fulllightpink4_apop_hits <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Dermo_Sus_fulllightpink4_apop_hits.txt", col.names = c("fromNode","toNode","weight","direction", "fromAltName","toAltName"))

#deLorg_Res_fullturquoise_apop_hits <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/deLorg_Res_fullturquoise_apop_hits.txt", col.names = c("fromNode","toNode","weight","direction", "fromAltName","toAltName"))
#deLorg_Sus_fullturquoise_apop_hits <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/deLorg_Sus_fullturquoise_apop_hits.txt", col.names = c("fromNode","toNode","weight","direction", "fromAltName","toAltName"))
#
#He_fullpurple_apop_hits <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/He_fullpurple_apop_hits.txt", col.names = c("fromNode","toNode","weight","direction", "fromAltName","toAltName"))
#He_fullyellow_apop_hits <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/He_fullyellow_apop_hits.txt", col.names = c("fromNode","toNode","weight","direction", "fromAltName","toAltName"))
#
#Zhang_fullblack_apop_hits <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Zhang_fullblack_apop_hits.txt", col.names = c("fromNode","toNode","weight","direction", "fromAltName","toAltName"))
#
#Rubio_fullblue_apop_hits <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Rubio_fullblue_apop_hits.txt", col.names = c("fromNode","toNode","weight","direction", "fromAltName","toAltName"))
#Rubio_fullbrown_apop_hits <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Rubio_fullbrown_apop_hits.txt", col.names = c("fromNode","toNode","weight","direction", "fromAltName","toAltName"))
#Rubio_fullmagenta_apop_hits <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Rubio_fullmagenta_apop_hits.txt", col.names = c("fromNode","toNode","weight","direction", "fromAltName","toAltName"))
#Rubio_fullturquoise_apop_hits <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Rubio_fullturquoise_apop_hits.txt", col.names = c("fromNode","toNode","weight","direction", "fromAltName","toAltName"))
Rubio_fullblack_apop_hits <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Rubio_fullblack_apop_hits.txt", col.names = c("fromNode","toNode","weight","direction", "fromAltName","toAltName"))


#Pro_RE22_Pro_fulldarkslateblue_apop_hits  <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_Pro_fulldarkslateblue_apop_hits.txt", col.names = c("fromNode","toNode","weight","direction", "fromAltName","toAltName"))
#Pro_RE22_Pro_fullroyalblue_apop_hits <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_Pro_fullroyalblue_apop_hits.txt", col.names = c("fromNode","toNode","weight","direction", "fromAltName","toAltName"))
#Pro_RE22_Pro_fullsteelblue_apop_hits <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_Pro_fullsteelblue_apop_hits.txt", col.names = c("fromNode","toNode","weight","direction", "fromAltName","toAltName"))
#Pro_RE22_Pro_fullturquoise_apop_hits <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_Pro_fullturquoise_apop_hits.txt", col.names = c("fromNode","toNode","weight","direction", "fromAltName","toAltName"))
Pro_RE22_Pro_fullskyblue3_apop_hits <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_Pro_fullskyblue3_apop_hits.txt", col.names = c("fromNode","toNode","weight","direction", "fromAltName","toAltName"))
Pro_RE22_Pro_fullwhite_apop_hits <- read.table(file="/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_Pro_fullwhite_apop_hits.txt", col.names = c("fromNode","toNode","weight","direction", "fromAltName","toAltName"))

# join annotations with the to and from nodes 
#Dermo_Tol_fullturquoise_apop_hits_annot <- Dermo_Tol_fullturquoise_apop_hits %>% dplyr::rename(ID = fromNode) %>% 
#  left_join(., unique(C_vir_rtracklayer[,c("ID","product")])) %>% dplyr::rename(fromNode = ID) %>% dplyr::rename(fromNodeproduct = product) %>%
#  dplyr::rename(ID = toNode) %>% left_join(., unique(C_vir_rtracklayer[,c("ID","product")])) %>%
#  dplyr::rename(toNode = ID) %>% dplyr::rename(toNodeproduct = product)

Dermo_Sus_fulllightpink4_apop_hits_annot <- Dermo_Sus_fulllightpink4_apop_hits %>% dplyr::rename(ID = fromNode) %>% 
  left_join(., unique(C_vir_rtracklayer[,c("ID","product")])) %>% dplyr::rename(fromNode = ID) %>% dplyr::rename(fromNodeproduct = product) %>%
  dplyr::rename(ID = toNode) %>% left_join(., unique(C_vir_rtracklayer[,c("ID","product")])) %>%
  dplyr::rename(toNode = ID) %>% dplyr::rename(toNodeproduct = product)

#deLorg_Res_fullturquoise_apop_hits_annot <- deLorg_Res_fullturquoise_apop_hits %>% dplyr::rename(transcript_id = fromNode) %>% 
#  left_join(., unique(C_gig_rtracklayer[,c("transcript_id","product")])) %>% dplyr::rename(fromNode = transcript_id) %>% dplyr::rename(fromNodeproduct = product) %>%
#dplyr::rename(transcript_id = toNode) %>% left_join(., unique(C_gig_rtracklayer[,c("transcript_id","product")])) %>%
#  dplyr::rename(toNode = transcript_id) %>% dplyr::rename(toNodeproduct = product)
#
#deLorg_Sus_fullturquoise_apop_hits_annot <- deLorg_Sus_fullturquoise_apop_hits %>% dplyr::rename(transcript_id = fromNode) %>% 
#  left_join(., unique(C_gig_rtracklayer[,c("transcript_id","product")])) %>% dplyr::rename(fromNode = transcript_id) %>% dplyr::rename(fromNodeproduct = product) %>%
#dplyr::rename(transcript_id = toNode) %>% left_join(., unique(C_gig_rtracklayer[,c("transcript_id","product")])) %>%
#  dplyr::rename(toNode = transcript_id) %>% dplyr::rename(toNodeproduct = product)
#
#He_fullpurple_apop_hits_annot <- He_fullpurple_apop_hits  %>% dplyr::rename(transcript_id = fromNode) %>% 
#  left_join(., unique(C_gig_rtracklayer[,c("transcript_id","product")])) %>% dplyr::rename(fromNode = transcript_id) %>% dplyr::rename(fromNodeproduct = product) %>%
#dplyr::rename(transcript_id = toNode) %>% left_join(., unique(C_gig_rtracklayer[,c("transcript_id","product")])) %>%
#  dplyr::rename(toNode = transcript_id) %>% dplyr::rename(toNodeproduct = product)
#
#He_fullyellow_apop_hits_annot <- He_fullyellow_apop_hits %>% dplyr::rename(transcript_id = fromNode) %>% 
#  left_join(., unique(C_gig_rtracklayer[,c("transcript_id","product")])) %>% dplyr::rename(fromNode = transcript_id) %>% dplyr::rename(fromNodeproduct = product) %>%
#dplyr::rename(transcript_id = toNode) %>% left_join(., unique(C_gig_rtracklayer[,c("transcript_id","product")])) %>%
#  dplyr::rename(toNode = transcript_id) %>% dplyr::rename(toNodeproduct = product)
#
#Zhang_fullblack_apop_hits_annot <- Zhang_fullblack_apop_hits %>% dplyr::rename(transcript_id = fromNode) %>% 
#  left_join(., unique(C_gig_rtracklayer[,c("transcript_id","product")])) %>% dplyr::rename(fromNode = transcript_id) %>% dplyr::rename(fromNodeproduct = product) %>%
#dplyr::rename(transcript_id = toNode) %>% left_join(., unique(C_gig_rtracklayer[,c("transcript_id","product")])) %>%
#  dplyr::rename(toNode = transcript_id) %>% dplyr::rename(toNodeproduct = product)
#
#Rubio_fullblue_apop_hits_annot <- Rubio_fullblue_apop_hits %>% dplyr::rename(transcript_id = fromNode) %>% 
#  left_join(., unique(C_gig_rtracklayer[,c("transcript_id","product")])) %>% dplyr::rename(fromNode = transcript_id) %>% dplyr::rename(fromNodeproduct = product) %>%
#dplyr::rename(transcript_id = toNode) %>% left_join(., unique(C_gig_rtracklayer[,c("transcript_id","product")])) %>%
#  dplyr::rename(toNode = transcript_id) %>% dplyr::rename(toNodeproduct = product)
#
#Rubio_fullbrown_apop_hits_annot <- Rubio_fullbrown_apop_hits %>% dplyr::rename(transcript_id = fromNode) %>% 
#  left_join(., unique(C_gig_rtracklayer[,c("transcript_id","product")])) %>% dplyr::rename(fromNode = transcript_id) %>% dplyr::rename(fromNodeproduct = product) %>%
#dplyr::rename(transcript_id = toNode) %>% left_join(., unique(C_gig_rtracklayer[,c("transcript_id","product")])) %>%
#  dplyr::rename(toNode = transcript_id) %>% dplyr::rename(toNodeproduct = product)
#
#Rubio_fullmagenta_apop_hits_annot <- Rubio_fullmagenta_apop_hits %>% dplyr::rename(transcript_id = fromNode) %>% 
#  left_join(., unique(C_gig_rtracklayer[,c("transcript_id","product")])) %>% dplyr::rename(fromNode = transcript_id) %>% dplyr::rename(fromNodeproduct = product) %>%
#dplyr::rename(transcript_id = toNode) %>% left_join(., unique(C_gig_rtracklayer[,c("transcript_id","product")])) %>%
#  dplyr::rename(toNode = transcript_id) %>% dplyr::rename(toNodeproduct = product)
#
#Rubio_fullturquoise_apop_hits_annot <- Rubio_fullturquoise_apop_hits %>% dplyr::rename(transcript_id = fromNode) %>% 
#  left_join(., unique(C_gig_rtracklayer[,c("transcript_id","product")])) %>% dplyr::rename(fromNode = transcript_id) %>% dplyr::rename(fromNodeproduct = product) %>%
#dplyr::rename(transcript_id = toNode) %>% left_join(., unique(C_gig_rtracklayer[,c("transcript_id","product")])) %>%
#  dplyr::rename(toNode = transcript_id) %>% dplyr::rename(toNodeproduct = product)

Rubio_fullblack_apop_hits_annot <- Rubio_fullblack_apop_hits %>% dplyr::rename(transcript_id = fromNode) %>% 
  left_join(., unique(C_gig_rtracklayer[,c("transcript_id","product")])) %>% dplyr::rename(fromNode = transcript_id) %>% dplyr::rename(fromNodeproduct = product) %>%
  dplyr::rename(transcript_id = toNode) %>% left_join(., unique(C_gig_rtracklayer[,c("transcript_id","product")])) %>%
  dplyr::rename(toNode = transcript_id) %>% dplyr::rename(toNodeproduct = product)

#Pro_RE22_Pro_fulldarkslateblue_apop_hits_annot <- Pro_RE22_Pro_fulldarkslateblue_apop_hits  %>% dplyr::rename(ID = fromNode) %>% 
#  left_join(., unique(C_vir_rtracklayer[,c("ID","product")])) %>% dplyr::rename(fromNode = ID) %>% dplyr::rename(fromNodeproduct = product) %>%
#dplyr::rename(ID = toNode) %>% left_join(., unique(C_vir_rtracklayer[,c("ID","product")])) %>%
#  dplyr::rename(toNode = ID) %>% dplyr::rename(toNodeproduct = product)
#
#Pro_RE22_Pro_fullroyalblue_apop_hits_annot <- Pro_RE22_Pro_fullroyalblue_apop_hits  %>% dplyr::rename(ID = fromNode) %>% 
#  left_join(., unique(C_vir_rtracklayer[,c("ID","product")])) %>% dplyr::rename(fromNode = ID) %>% dplyr::rename(fromNodeproduct = product) %>%
#dplyr::rename(ID = toNode) %>% left_join(., unique(C_vir_rtracklayer[,c("ID","product")])) %>%
#  dplyr::rename(toNode = ID) %>% dplyr::rename(toNodeproduct = product)
#
#Pro_RE22_Pro_fullsteelblue_apop_hits_annot <- Pro_RE22_Pro_fullsteelblue_apop_hits %>% dplyr::rename(ID = fromNode) %>% 
#  left_join(., unique(C_vir_rtracklayer[,c("ID","product")])) %>% dplyr::rename(fromNode = ID) %>% dplyr::rename(fromNodeproduct = product) %>%
#dplyr::rename(ID = toNode) %>% left_join(., unique(C_vir_rtracklayer[,c("ID","product")])) %>%
#  dplyr::rename(toNode = ID) %>% dplyr::rename(toNodeproduct = product)
#
#Pro_RE22_Pro_fullturquoise_apop_hits_annot <- Pro_RE22_Pro_fullturquoise_apop_hits  %>% dplyr::rename(ID = fromNode) %>% 
#  left_join(., unique(C_vir_rtracklayer[,c("ID","product")])) %>% dplyr::rename(fromNode = ID) %>% dplyr::rename(fromNodeproduct = product) %>%
#dplyr::rename(ID = toNode) %>% left_join(., unique(C_vir_rtracklayer[,c("ID","product")])) %>%
#  dplyr::rename(toNode = ID) %>% dplyr::rename(toNodeproduct = product)

Pro_RE22_Pro_fullskyblue3_apop_hits_annot <- Pro_RE22_Pro_fullskyblue3_apop_hits %>% dplyr::rename(ID = fromNode) %>% 
  left_join(., unique(C_vir_rtracklayer[,c("ID","product")])) %>% dplyr::rename(fromNode = ID) %>% dplyr::rename(fromNodeproduct = product) %>%
  dplyr::rename(ID = toNode) %>% left_join(., unique(C_vir_rtracklayer[,c("ID","product")])) %>%
  dplyr::rename(toNode = ID) %>% dplyr::rename(toNodeproduct = product)

Pro_RE22_Pro_fullwhite_apop_hits_annot <- Pro_RE22_Pro_fullwhite_apop_hits  %>% dplyr::rename(ID = fromNode) %>% 
  left_join(., unique(C_vir_rtracklayer[,c("ID","product")])) %>% dplyr::rename(fromNode = ID) %>% dplyr::rename(fromNodeproduct = product) %>%
  dplyr::rename(ID = toNode) %>% left_join(., unique(C_vir_rtracklayer[,c("ID","product")])) %>%
  dplyr::rename(toNode = ID) %>% dplyr::rename(toNodeproduct = product)

# Export the annotated lists
#write.table(Dermo_Tol_fullturquoise_apop_hits_annot , file= "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Dermo_Tol_fullturquoise_apop_hits_annot.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(Dermo_Sus_fulllightpink4_apop_hits_annot , file= "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Dermo_Sus_fulllightpink4_apop_hits_annot.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

#write.table(deLorg_Res_fullturquoise_apop_hits_annot , file= "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/deLorg_Res_fullturquoise_apop_hits_annot.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
#write.table(deLorg_Sus_fullturquoise_apop_hits_annot , file= "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/deLorg_Sus_fullturquoise_apop_hits_annot.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
#write.table(He_fullpurple_apop_hits_annot , file= "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/He_fullpurple_apop_hits_annot.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
#write.table(He_fullyellow_apop_hits_annot , file= "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/He_fullyellow_apop_hits_annot.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
#write.table(Zhang_fullblack_apop_hits_annot , file= "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Zhang_fullblack_apop_hits_annot.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
#write.table(Rubio_fullblue_apop_hits_annot , file= "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Rubio_fullblue_apop_hits_annot.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
#write.table(Rubio_fullbrown_apop_hits_annot , file= "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Rubio_fullbrown_apop_hits_annot.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
#write.table(Rubio_fullmagenta_apop_hits_annot , file= "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Rubio_fullmagenta_apop_hits_annot.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
#write.table(Rubio_fullturquoise_apop_hits_annot , file= "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Rubio_fullturquoise_apop_hits_annot.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(Rubio_fullblack_apop_hits_annot, file= "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Rubio_fullblack_apop_hits_annot.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

#write.table(Pro_RE22_Pro_fulldarkslateblue_apop_hits_annot , file= "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_Pro_fulldarkslateblue_apop_hits_annot.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
#write.table(Pro_RE22_Pro_fullroyalblue_apop_hits_annot , file= "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_Pro_fullroyalblue_apop_hits_annot.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
#write.table(Pro_RE22_Pro_fullsteelblue_apop_hits_annot , file= "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_Pro_fullsteelblue_apop_hits_annot.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
#write.table(Pro_RE22_Pro_fullturquoise_apop_hits_annot , file= "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_Pro_fullturquoise_apop_hits_annot.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(Pro_RE22_Pro_fullskyblue3_apop_hits_annot , file= "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_Pro_fullskyblue3_apop_hits_annot.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(Pro_RE22_Pro_fullwhite_apop_hits_annot , file= "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/Pro_RE22_Pro_fullwhite_apop_hits_annot.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)


# Clear global workspace
rm(list = ls())
