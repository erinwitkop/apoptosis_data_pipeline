# Script to Analyze Differential Expression of Transcriptomes using DESeq2
# Erin Roberts, PhD Candidate University of Rhode Island 
# 2/13/2020

# This script will calculate differential expression of apoptosis genes from C. gigas and C. virginica, for each experiment separately.
# Comparisons and formulas used to calculate differential expression for each experiment will be unique to that experiment, specifically
# tailored to when the infection was most acute. Challenge group samlpes will always be compared to their own control.

# This script contains all the multiple testing performed by Erin Roberts during DEG analysis, while C_gig_C_vir_transcriptome_DE_GO_Plotting.R 
# Contains the final DEG analysis presented in publication. 

#### LOAD PACKAGES ####

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
library(ComplexHeatmap)
library(reshape2)
library(plyr)
library(Repitools)
library(purrr)
library(tibble)
library(ggfortify)
library(ggpubr)
library(limma)

# VERSIONS (see sessionInfo at bottom of script for full information)
# R version 3.6.1 (2019-07-05)
# DESeq2_1.24.0 

## Resources used to help build code
# Updated 2020 vignette: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html 
# Updated 2018 Vignette of DESeq2 : https://bioconductor.org/packages/3.7/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
# Updated Workflow for DESeq2: https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#the-deseqdataset-object-sample-information-and-the-design-formula
# Oct. 2019 workflow: https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#the-deseqdataset-object-sample-information-and-the-design-formula 
# Deseq2 Manual updated July 18th, 2019: https://www.bioconductor.org/packages/devel/bioc/manuals/DESeq2/man/DESeq2.pdf
# the results section of this is particularly useful

##### LOADING GENOME DATA AND GO TERMS ##### 

#### C. virginica 

# Import gff file with rtracklayer
C_vir_rtracklayer <- import("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/ref_C_virginica-3.0_top_level.gff3")
C_vir_rtracklayer <- as.data.frame(C_vir_rtracklayer)

# Load finished C_vir_rtracklayer_transcripts_GO.csv with GO terms mapped using code and processed below. Do not repeat every time.
C_vir_rtracklayer_transcripts_GO <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/C_vir_rtracklayer_transcripts_GO.csv", header=TRUE)

# Isolate transcript lines in GFF annotation so that I can use Batch Entrez lookup for their parent proteins
#C_vir_rtracklayer_transcripts <- filter(C_vir_rtracklayer, grepl("XM", transcript_id))
#C_vir_rtracklayer_transcripts_unique <- subset(C_vir_rtracklayer_transcripts, !duplicated(transcript_id))
#write.table(file="./Dermo_2015_Analysis/OUTPUT/C_vir_unique_transcripts.txt", C_vir_rtracklayer_transcripts_unique$transcript_id) 
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
#C_vir_XM_with_XP <- import("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/C_VIRGINICA_PIPELINE/Dermo_2015_Analysis/OUTPUT/combined_batch_lookup_cds.gff")
#C_vir_XM_with_XP <- as.data.frame(C_vir_XM_with_XP)
#C_vir_XM_with_XP_index <- C_vir_XM_with_XP[,c("seqnames","protein_id")]
#colnames(C_vir_XM_with_XP_index)[1] <- "transcript_id"

# Load in Interproscan GO annotation from LSU Kevin using Rtracklayer
# cat all the edited header removed files from Interproscan Kevin LSU files Import gff file with rtracklayer
# for i in edited*.gff3; do cat $i >> combined_CV_prot_id_interproscan.gff3 ; done

# In terminal, combine all the gff files and remove the comment lines from the script 
#Cvir_Interproscan1 <- import("/Users/erinroberts/Documents/PhD_Research/GOMEZCHIARRI2_FILES_DAILYNOTES/GENOME_TRANSCRIPTOME_SEQUENCES/IP5_out/protein_id/combined_CV_prot_id_interproscan.gff3")
#class(Cvir_Interproscan1)
# convert to dataframe object using annoGR2DF (repitools)
#Cvir_Interproscan_DF <- annoGR2DF(Cvir_Interproscan1)
#class(Cvir_Interproscan_DF) #data.frame

# Interproscan file has multiple lines per entry, need to subset for the lines that do have a GO entry because not all of them do, then subset for unique protein lines
#class(Cvir_Interproscan_DF$Ontology_term)
#Cvir_Interproscan_DF$Ontology_term <- as.character(Cvir_Interproscan_DF$Ontology_term)
#C_vir_Interproscan_GO <- Cvir_Interproscan_DF %>% filter(Ontology_term !="character(0)")
#head(C_vir_Interproscan_GO)

# keep lines with unique protein and GO terms 
#C_vir_Interproscan_GO_unique <- C_vir_Interproscan_GO[!duplicated(C_vir_Interproscan_GO[,c("chr","Ontology_term")]),]

# Format GO term column correctly 
#C_vir_Interproscan_GO_unique$Ontology_term <- gsub("[^[:alnum:][:blank:]+:/\\-]", "", C_vir_Interproscan_GO_unique$Ontology_term )
#C_vir_Interproscan_GO_unique$Ontology_term <- gsub(C_vir_Interproscan_GO_unique$Ontology_term, pattern="c\\", replacement="", fixed=TRUE)
#C_vir_Interproscan_GO_unique$Ontology_term <- gsub(C_vir_Interproscan_GO_unique$Ontology_term, pattern="\\", replacement="", fixed=TRUE)
#C_vir_Interproscan_GO_unique$Ontology_term <- gsub(C_vir_Interproscan_GO_unique$Ontology_term, pattern=" ",replacement =",", fixed=TRUE)
#class(C_vir_Interproscan_GO_unique)

# merge GO terms with the same protein name so there aren't multiple lines for a transcript
#C_vir_Interproscan_GO_unique <- ddply(C_vir_Interproscan_GO_unique, "chr", summarize, Combined_ontology = toString(Ontology_term))

# Remove duplicate GO strings in the same column and create new column
#C_vir_Interproscan_GO_unique$Combined_ontology <- gsub(" ", "", C_vir_Interproscan_GO_unique$Combined_ontology)
#C_vir_Interproscan_GO_unique <- mutate(C_vir_Interproscan_GO_unique, unique_go = map(str_split(Combined_ontology, ","), unique))

# make new column that removes the "_ORF" from the end of the seqnames columns
#C_vir_Interproscan_GO_unique$protein_id <- str_remove(C_vir_Interproscan_GO_unique$chr, "_ORF")

# merge XM and XP list with Interproscan list
#C_vir_Interproscan_GO_unique_XM_merged <- C_vir_Interproscan_GO_unique %>% left_join(select(C_vir_XM_with_XP_index, "transcript_id","protein_id"), by = "protein_id")

# Merge Interproscan list onto transcript annotation using left_join of the TRANSCRIPT ID so that it joins correctly at the level of transcripts
# non unique file
#C_vir_rtracklayer_GO <-  C_vir_rtracklayer %>% left_join(select(C_vir_Interproscan_GO_unique_XM_merged,"unique_go","transcript_id"), by = "transcript_id")

#unique file with one line per transcript
#C_vir_rtracklayer_transcripts_GO <- C_vir_rtracklayer_transcripts_unique %>% left_join(select(C_vir_Interproscan_GO_unique_XM_merged,"unique_go","transcript_id"), by = "transcript_id")

# are they all NULL?
#C_vir_rtracklayer_GO %>% filter(unique_go !="NULL") # nope they are not!!

# Export for use in cluster later
#class(C_vir_rtracklayer_transcripts_GO)
# coerce data frame to be all character
# First coerce the data.frame to all-character
#C_vir_rtracklayer_transcripts_GO_2 = data.frame(lapply(C_vir_rtracklayer_transcripts_GO, as.character), stringsAsFactors=FALSE)
#head(C_vir_rtracklayer_transcripts_GO_2)
# now write to csv
#write.csv(C_vir_rtracklayer_transcripts_GO_2, file="C_vir_rtracklayer_transcripts_GO.csv")

### C. gigas
# C. gigas genome already has GO terms mapped to it. 

# Import gff file, using new version of genome annotation
C_gig_rtracklayer <- import("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/GCF_000297895.1_oyster_v9_genomic.gff")
C_gig_rtracklayer <- as.data.frame(C_gig_rtracklayer)
C_gig_rtracklayer_transcripts <-  C_gig_rtracklayer %>% filter(type =="mRNA")

#### IMPORT APOPTOSIS GENE NAMES LISTS FOR EACH SPECIES AND MAP ####
## NOTE: THIS CODE IS THE SAME AS THAT TO IDENTIFY APOPTOSIS GENES IN THE GENOME FOR THE APOPTOSIS PATHWAY ANNOTATION PART OF THE PAPER FROM Identify_count_annotated_apop_genes.R

Apoptosis_names_list <- c('bcl-2-related protein A1',
                          'apoptosis-inducing factor 1',
                          'akt1',
                          'RAC-alpha serine/threonine-protein kinase-',
                          'RAC-gamma serine/threonine-protein kinase-',
                          'methylthioribulose-1-phosphate dehydratase',
                          'tumor necrosis factor ligand superfamily member 10',
                          'tumor necrosis factor superfamily member 12',
                          'cell death regulator Aven',
                          'BCL2 associated agonist of cell death',
                          'BAG family molecular chaperone regulator',
                          'bcl-2 homologous antagonist/killer',
                          'apoptosis regulator BAX',
                          'bcl-2-like protein 2',
                          'bcl-2-like protein 1',
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
                          'caspase-1-like',
                          'caspase-2',
                          'caspase-3-like',
                          'caspase-7',
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
                          'caspase activity and apoptosis inhibitor 1',
                          'baculoviral IAP repeat-containing protein',
                          'cAMP-responsive element',
                          'cytochrome c-like',
                          'death-associated inhibitor of apoptosis',
                          'tumor necrosis factor receptor superfamily member 25',
                          'tumor necrosis factor receptor superfamily member 10A',
                          'tumor necrosis factor receptor superfamily member 10B',
                          'endonuclease G, mitochondrial',
                          'FAS-associated death domain protein',
                          'fas apoptotic inhibitory molecule 1',
                          'tumor necrosis factor receptor superfamily member 6',
                          'fas cell surface death receptor',
                          'GTPase IMAP family member',
                          'harakiri',
                          'baculoviral IAP repeat-containing protein',
                          'DNA fragmentation factor subunit alpha',
                          'interferon-induced protein 44',
                          'NF-kappa-B inhibitor alpha',
                          'NF-kappa-B inhibitor epsilon',
                          'inositol 1,4,5-trisphosphate receptor',
                          'stress-activated protein kinase JNK',
                          'lipopolysaccharide-induced tumor necrosis factor-alpha',
                          'induced myeloid leukemia cell differentiation protein Mcl-1-like',
                          'mitogen-activated protein kinase kinase kinase 1-like',
                          'mitogen-activated protein kinase 1',
                          'dual specificity mitogen-activated protein kinase kinase 1',
                          'mitogen-activated protein kinase kinase kinase 7-like',
                          'transcriptional regulator Myc-A',
                          'myeloid differentiation primary response protein MyD88',
                          'phorbol-12-myristate-13-acetate-induced protein 1',
                          'nuclear factor NF-kappa-B p105 subunit',
                          'nuclear factor NF-kappa-B p100 subunit',
                          'transcription factor p65',
                          'RELB proto-oncogene',
                          'NF-kB subunit',
                          'reticuloendotheliosis oncogene',
                          'anti-apoptotic protein NR13',
                          'nuclear mitotic apparatus protein 1',
                          'dynamin-like 120 kDa protein',
                          'cyclin-dependent kinase 5 activator 1',
                          'cellular tumor antigen p53',
                          'programmed cell death protein',
                          'p53 and DNA damage-regulated protein 1',
                          'phosphatidylinositol 3-kinase',
                          'putative inhibitor of apoptosis',
                          'cAMP-dependent protein kinase',
                          'protein kinase C delta type',
                          'protein kinase C iota type',
                          'BCL2 binding component 3',
                          'cdc42 homolog',
                          'ras-like GTP-binding protein rho',
                          'rho-related GTP-binding protein RhoE-like',
                          'ras-related C3 botulinum toxin substrate 1',
                          'rho-related protein racA',
                          'mitochondrial Rho GTPase 1',
                          'receptor-interacting serine/threonine-protein kinase 1',
                          'receptor-interacting serine/threonine-protein kinase 4',
                          'diablo homolog, mitochondrial',
                          'toll-like receptor',
                          'tumor necrosis factor',
                          'lymphotoxin-alpha',
                          'CD40 ligand',
                          'tumor necrosis factor receptor superfamily member',
                          'TNFRSF1A associated via death domain',
                          'TNF receptor-associated factor',
                          'E3 ubiquitin-protein ligase XIAP',
                          'netrin receptor',
                          'neurotrophic receptor tyrosine kinase 1',
                          'sonic hedgehog receptor',
                          'receptor-interacting serine/threonine-protein kinase',
                          'mixed lineage kinase domain',
                          'heat shock protein',
                          'E3 ubiquitin-protein ligase CHIP',
                          'tumor necrosis factor alpha-induced protein 3',
                          'protein phosphatase 1B',
                          'aurora kinase A',
                          'glutathione peroxidase 4',
                          'gasdermin',
                          'poly \\[ADP-ribose]\\ polymerase 1-like',
                          'macrophage migration inhibitory factor',
                          'hexokinase-1',
                          'Raf-1 protooncogene serine/threonine kinase',
                          'elastase, neutrophil expressed',
                          'cathepsin',
                          'PRKC apoptosis WT1 regulator protein',
                          'apoptosis-stimulating of p53 protein 1',
                          'apoptosis-stimulating of p53 protein 2',
                          'apoptosis inhibitor 5',
                          'apoptotic chromatin condensation inducer in the nucleus',
                          "high mobility group box 1",
                          'ceramide synthase',
                          'cyclic AMP-responsive element-binding protein',
                          'cell death-inducing p53-target protein 1',
                          'TP53-binding protein 1',
                          'p53-induced death domain-containing protein 1',
                          'death domain-containing protein CRADD',
                          'p63',
                          'p73')

#### Grep Apoptosis protein names in genome files ####

# C virginica 
C_vir_rtracklayer_mRNA <- C_vir_rtracklayer %>% filter(type == "mRNA")
C_vir_rtracklayer_apop_product <- C_vir_rtracklayer_mRNA[grepl(paste(Apoptosis_names_list,collapse="|"), 
                                                               C_vir_rtracklayer_mRNA$product, ignore.case = TRUE),]
nrow(C_vir_rtracklayer_apop_product) #1106

# Terms to remove
# remove complement C1q proteins, dual specificity protein phosphatase 1B-like, remove kunitz-type, and NOT other kDA protein names so I can keep all heat shock proteins
C_vir_rtracklayer_apop_product_final <-C_vir_rtracklayer_apop_product[!grepl("complement C1q", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) & 
                                                                         !grepl("dual specificity protein phosphatase 1B-like", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
                                                                         !grepl("kunitz-type", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
                                                                         !grepl("mannosyl", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
                                                                         !grepl("S-acyl", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
                                                                         !grepl("dynamin-like",C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
                                                                         !grepl("activator of 90 kDa", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
                                                                         !grepl("zinc finger protein", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) & 
                                                                         !grepl("phosphatidylinositol 3-kinase 1", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE)  & 
                                                                         !grepl("phosphatidylinositol 3-kinase 2", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE)  & 
                                                                         !grepl("DDB_G0272098",C_vir_rtracklayer_apop_product$product, ignore.case = TRUE)  & 
                                                                         !grepl("transcriptional regulator Myc-A",C_vir_rtracklayer_apop_product$product, ignore.case = TRUE)  & 
                                                                         !grepl("caspase-14", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
                                                                         !grepl("WD repeat-containing protein WRAP73", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
                                                                         !grepl("tumor protein p63-regulated gene 1-like protein", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE),]
nrow(C_vir_rtracklayer_apop_product_final) #1026

## Identify CG apoptosis genes 
Apoptosis_names_list_CG <- c('bcl-2-related protein A1',
                             'apoptosis-inducing factor 1',
                             'akt1',
                             'RAC-alpha serine/threonine-protein kinase-',
                             'RAC-gamma serine/threonine-protein kinase-',
                             'methylthioribulose-1-phosphate dehydratase',
                             'tumor necrosis factor ligand superfamily member 10',
                             'tumor necrosis factor superfamily member 12',
                             'cell death regulator Aven',
                             'BCL2 associated agonist of cell death',
                             'BAG family molecular chaperone regulator',
                             'bcl-2 homologous antagonist/killer',
                             'apoptosis regulator BAX',
                             'bcl-2-like protein 2',
                             'bcl-2-like protein 1',
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
                             'caspase-1-like',
                             'caspase-2',
                             'caspase-3-like',
                             'caspase-7',
                             'caspase-8',
                             'caspase-9',
                             'caspase-10',
                             'caspase-11',
                             'caspase-6',
                             'caspase-4',
                             'caspase-5',
                             'cell division cycle and apoptosis regulator protein 1',
                             'CD151 antigen',
                             'B-cell translocation gene 1',
                             'interferon alpha',
                             'caspase activity and apoptosis inhibitor 1',
                             'baculoviral IAP repeat-containing protein',
                             'cAMP-responsive element',
                             'cytochrome c',
                             'death-associated inhibitor of apoptosis',
                             'tumor necrosis factor receptor superfamily member 25',
                             'tumor necrosis factor receptor superfamily member 10A',
                             'tumor necrosis factor receptor superfamily member 10B',
                             'endonuclease G, mitochondrial',
                             'FAS-associated death domain protein',
                             'fas apoptotic inhibitory molecule 1',
                             'tumor necrosis factor receptor superfamily member 6',
                             'fas cell surface death receptor',
                             'GTPase IMAP family member',
                             'harakiri',
                             'baculoviral IAP repeat-containing protein',
                             'DNA fragmentation factor subunit alpha',
                             'interferon-induced protein 44',
                             'NF-kappa-B inhibitor alpha',
                             'NF-kappa-B inhibitor epsilon',
                             'inositol 1,4,5-trisphosphate receptor',
                             'stress-activated protein kinase JNK',
                             'lipopolysaccharide-induced tumor necrosis factor-alpha',
                             'induced myeloid leukemia cell differentiation protein Mcl-1',
                             'mitogen-activated protein kinase kinase kinase 1,',
                             'mitogen-activated protein kinase 1',
                             'dual specificity mitogen-activated protein kinase kinase 1',
                             'mitogen-activated protein kinase kinase kinase 7',
                             'transcriptional regulator Myc-A',
                             'myeloid differentiation primary response protein MyD88',
                             'phorbol-12-myristate-13-acetate-induced protein 1',
                             'nuclear factor NF-kappa-B p105 subunit',
                             'nuclear factor NF-kappa-B p100 subunit',
                             'transcription factor p65',
                             'RELB proto-oncogene',
                             'NF-kB subunit',
                             'reticuloendotheliosis oncogene',
                             'anti-apoptotic protein NR13',
                             'nuclear mitotic apparatus protein 1',
                             'dynamin-like 120 kDa protein',
                             'cyclin-dependent kinase 5 activator 1',
                             'cellular tumor antigen p53',
                             'programmed cell death protein',
                             'p53 and DNA damage-regulated protein 1',
                             'phosphatidylinositol 3-kinase',
                             'putative inhibitor of apoptosis',
                             'cAMP-dependent protein kinase',
                             'protein kinase C delta type',
                             'protein kinase C iota type',
                             'BCL2 binding component 3',
                             'cdc42 homolog',
                             'ras-like GTP-binding protein rho',
                             'rho-related GTP-binding protein RhoE',
                             'ras-related C3 botulinum toxin substrate 1',
                             'rho-related protein racA',
                             'mitochondrial Rho GTPase 1',
                             'receptor-interacting serine/threonine-protein kinase 1',
                             'receptor-interacting serine/threonine-protein kinase 4',
                             'diablo homolog, mitochondrial',
                             'toll-like receptor',
                             'tumor necrosis factor',
                             'lymphotoxin-alpha',
                             'CD40 ligand',
                             'tumor necrosis factor receptor superfamily member',
                             'TNFRSF1A associated via death domain',
                             'TNF receptor-associated factor',
                             'E3 ubiquitin-protein ligase XIAP',
                             'netrin receptor',
                             'neurotrophic receptor tyrosine kinase 1',
                             'sonic hedgehog receptor',
                             'receptor-interacting serine/threonine-protein kinase',
                             'mixed lineage kinase domain',
                             'heat shock protein',
                             'E3 ubiquitin-protein ligase CHIP',
                             'tumor necrosis factor alpha-induced protein 3',
                             'protein phosphatase 1B',
                             'aurora kinase A',
                             'glutathione peroxidase 4',
                             'gasdermin',
                             'poly \\[ADP-ribose]\\ polymerase 1 ',
                             'poly \\[ADP-ribose]\\ polymerase 1-',
                             'macrophage migration inhibitory factor',
                             'hexokinase-1',
                             'Raf-1 protooncogene serine/threonine kinase',
                             'elastase, neutrophil expressed',
                             'cathepsin',
                             'PRKC apoptosis WT1 regulator protein',
                             'apoptosis-stimulating of p53 protein 1',
                             'apoptosis-stimulating of p53 protein 2',
                             'apoptosis inhibitor 5',
                             'apoptotic chromatin condensation inducer in the nucleus',
                             "high mobility group box 1",
                             'ceramide synthase',
                             'cyclic AMP-responsive element-binding protein',
                             'cell death-inducing p53-target protein 1',
                             'TP53-binding protein 1',
                             'p53-induced death domain-containing protein 1',
                             'death domain-containing protein CRADD',
                             'p63',
                             'p73')

#### Grep Apoptosis protein names in genome files CG ####
## Note unlike C. vir file, C. gigas file has gene names in the "product" column, and still includes comma with transcript info for some
# than a comman and the transcript name
# All transcripts from the same gene share the LOC ID in the "gene" column, though the Name column differs
# C gigas, filter for rows that have NA gene, this will also keep all the lines with transcript information. Filter out gbkey CDS 
# because I only care about the number of transcripts and genes
C_gig_rtracklayer_filtered <- C_gig_rtracklayer %>% filter(type =="mRNA")

# only genes have the description, the mRNA id's will match with the ID
C_gig_rtracklayer_apop_product <- C_gig_rtracklayer_filtered[grepl(paste(Apoptosis_names_list_CG,collapse="|"), 
                                                                   C_gig_rtracklayer_filtered$product, ignore.case = TRUE),]
nrow(C_gig_rtracklayer_apop_product) #797

# Terms to remove
# remove complement C1q proteins, dual specificity protein phosphatase 1B-like, remove kunitz-type, and other NOT kDA protein names so I can keep heat shock proteins in both
C_gig_rtracklayer_apop_product_final <- C_gig_rtracklayer_apop_product[!grepl("complement C1q", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
                                                                         # !grepl("dual specificity protein phosphatase 1B-like", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) & not present in list
                                                                         !grepl("kunitz-type", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
                                                                         !grepl("mannosyl", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
                                                                         !grepl("S-acyl", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
                                                                         !grepl("cytochrome c oxidase", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
                                                                         !grepl("cytochrome c-type heme lyase", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
                                                                         !grepl("cytochrome c1", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
                                                                         !grepl("interferon alpha-inducible protein 27-like protein 2",C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
                                                                         !grepl("dynamin-like",C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
                                                                         !grepl("activator of 90 kDa", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE)  & 
                                                                         !grepl("phosphatidylinositol 3-kinase 1", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE)  & 
                                                                         !grepl("phosphatidylinositol 3-kinase 2", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE)  & 
                                                                         !grepl("DDB_G0272098",C_gig_rtracklayer_apop_product$product, ignore.case = TRUE)  & 
                                                                         !grepl("transcriptional regulator Myc-A",C_gig_rtracklayer_apop_product$product, ignore.case = TRUE)  & 
                                                                         !grepl("caspase-14",C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
                                                                         !grepl("WD repeat-containing protein WRAP73", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
                                                                         !grepl("tumor protein p63-regulated gene 1-like protein", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE),]
nrow(C_gig_rtracklayer_apop_product_final) #660

#### ZHANG VIBRIO TRANSCRIPTOME ANALYSIS ####

## LOAD DATA
Zhang_counts <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/Zhang_transcript_count_matrix.csv", header=TRUE,
                         row.names = "transcript_id")
head(Zhang_counts)
colnames(Zhang_counts)

# remove MSTRG novel transcript lines (can assess these later if necessary)
Zhang_counts <- Zhang_counts[!grepl("MSTRG", row.names(Zhang_counts)),]

# Cute the "rna-" from the beginning of rownames
remove_rna = function(x){
  return(gsub("rna-","",x))
}
row.names(Zhang_counts) <- remove_rna(row.names(Zhang_counts))
head(Zhang_counts)

#Load in sample metadata
Zhang_coldata <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/Zhang_coldata2.csv", row.names = 1 )
View(Zhang_coldata)  
nrow(Zhang_coldata) 

# Make sure the columns of the count matrix and rows of the column data (sample metadata) are in the same order. 
all(rownames(Zhang_coldata) %in% colnames(Zhang_counts))  #Should return TRUE
# returns TRUE
all(colnames(Zhang_counts) %in% rownames(Zhang_coldata))  
# returns TRUE
all(rownames(Zhang_coldata) == colnames(Zhang_counts))    # should return TRUE
# returns TRUE


### DATA QC PCA PLOT 
# rlog transform data is recommended over vst for small data sets 
# PCA plots of data (https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/02_Preprocessing_Data.nb.html#count-distribution-boxplots)
Zhang_counts_matrix <- as.matrix(Zhang_counts)
Zhangrlogcounts <- rlog(Zhang_counts_matrix, blind =TRUE)
# run PCA
pcZhang <- prcomp(t(Zhangrlogcounts))
# plot PCA
autoplot(pcZhang)

# Lets add colour to look at the clustering for Status
autoplot(pcZhang,
         data = Zhang_coldata, 
         colour="condition", 
         size=5) # Vibrio tubiashi and LPS cluster together (somewhat)
# with MSTRG removed, V aes and V alg2 cluster, Mlut, LPS and PBS cluster
# control and PBS clustered together, LPS and M lut clustered together (see this in downstream PCAs too)
autoplot(pcZhang,
         data = Zhang_coldata, 
         colour="group_by_sim", 
         size=5) 

# Is there a time effect? 
autoplot(pcZhang,
         data = Zhang_coldata, 
         colour="time", 
         size=5)  # time may have some effect but because there are no replicates at this time point it is difficult to tell. 
                  #Going to remove this from the DESeq2 formula
# Where do the different conditions cluster?
autoplot(pcZhang,
         data = Zhang_coldata, 
         colour="condition", 
         size=5)
# Plot PCA 2 and 3 for comparison
autoplot(pcZhang,
         data = Zhang_coldata, 
         colour = "condition", 
         size = 5,
         x = 2,
         y = 3) # control and PBS clustering here, LPS and V. tubiashi still clustering
  # with MSTRG removed, V tub and V aes are closest cluster, they also cluster with LPS
autoplot(pcZhang,
         data = Zhang_coldata, 
         colour = "group", 
         size = 5,
         x = 2,
         y = 3)
autoplot(pcZhang,
         data = Zhang_coldata, 
         colour = "condition", 
         size = 5,
         x = 3,
         y = 4) # LPS and V. tubiashi still clustering
# with MSTRG removed, V ang and V alg1 cluster the closest 

## MAKE DESEQ DATA SET FROM MATRIX
# This object specifies the count data and metadata you will work with. The design piece is critical.
# Design here is to combine control and PBS as control and the others as challenge
# Correct for batch effects if necessary in this original formula: see this thread https://support.bioconductor.org/p/121408/
# do not correct counts using the removeBatchEffects from limma based on thread above 

## Creating three here so I can compare the results
Zhang_dds <- DESeqDataSetFromMatrix(countData = Zhang_counts,
                                    colData = Zhang_coldata,
                                    design = ~time + group_by_sim) # add time to control for injection and time effect
#Zhang_dds_path <- DESeqDataSetFromMatrix(countData = Zhang_counts,
#                                                 colData = Zhang_coldata,
#                                                 design = ~  path) # add time to control for injection and time effect

# USE THIS ONE IN FUTURE OTHERS WERE FOR TESTING
#Zhang_dds_broken_group <- DESeqDataSetFromMatrix(countData = Zhang_counts,
#                                                 colData = Zhang_coldata,
#                                                 design = ~time+ group_by_sim) # add time to control for injection and time effect

## Prefiltering the data
# Data prefiltering helps decrease the size of the data set and get rid of
# rows with no data or very minimal data (<10). Apply a minimal filtering here as more stringent filtering will be applied later
Zhang_dds <- Zhang_dds[ rowSums(counts(Zhang_dds)) > 10, ]
#Zhang_dds_broken_group <- Zhang_dds_broken_group[ rowSums(counts(Zhang_dds_broken_group)) > 10, ]
#Zhang_dds_path <- Zhang_dds_path[ rowSums(counts(Zhang_dds_path)) > 10, ]

## Check levels 
# It is prefered in R that the first level of a factor be the reference level for comparison
# (e.g. control, or untreated samples), so we can relevel the factor like so
# Check factor levels, set it so that comparison group is the first
levels(Zhang_coldata$condition) #"Control" "LPS"     "M_lut"   "PBS"     "V_aes"   "V_alg_1" "V_alg_2" "V_ang"   "V_tub"  
levels(Zhang_coldata$group_by_sim) # "control"             "LPS_M_lut"           "V_aes_V_alg1_V_alg2" "V_tub_V_ang"    
levels(Zhang_coldata$time) # 12hr, no injection
Zhang_dds$time <- factor(Zhang_dds$time , levels = c("No_injection","12h"))

## DATA TRANSFORMATION AND VISUALIZATION
# Assess sample clustering after setting initial formula for comparison
Zhang_dds_rlog <- rlog(Zhang_dds, blind = TRUE) # keep blind = true before deseq function has been run
#Zhang_dds_broken_group_rlog <- rlog(Zhang_dds_broken_group, blind = TRUE)
#Zhang_dds_path_rlog <- rlog(Zhang_dds_path, blind = TRUE)

## PCA plot visualization of individuals in the family 
plotPCA(Zhang_dds_rlog, intgroup=c("group_by_sim", "condition"))
  # control and PBS still clustering, LPS M. lut and V. aes clustering
  # with MSTRG removed, the LPS and M. lut cluster most closely 

### DIFFERENTIAL EXPRESSION ANALYSIS
# run pipeline with single command because the formula has already been specified
# Steps: estimation of size factors (controlling for differences in the sequencing depth of the samples), 
# the estimation of dispersion values for each gene, 
# and fitting a generalized linear model.

Zhang_dds_deseq <- DESeq(Zhang_dds) 

#Zhang_dds_broken_group_deseq <- DESeq(Zhang_dds_broken_group) # USE THIS ONE IN FUTURE OTHERS WERE FOR TESTING
#Zhang_dds_path_deseq <- DESeq(Zhang_dds_path)

## Check the resultsNames object of each to look at the available coefficients for use in lfcShrink command
resultsNames(Zhang_dds_deseq) # [1] "Intercept"                                   "time_12h_vs_No_injection"                   
#[3] "group_by_sim_LPS_M_lut_vs_control"           "group_by_sim_V_aes_V_alg1_V_alg2_vs_control"
#[5] "group_by_sim_V_tub_V_ang_vs_control" 
resultsNames(Zhang_dds_broken_group_deseq) 
#[1] "Intercept"                                   "time_No_injection_vs_12h"                    "group_by_sim_LPS_M_lut_vs_control"          
#[4] "group_by_sim_V_aes_V_alg1_V_alg2_vs_control" "group_by_sim_V_tub_V_ang_vs_control" 

#resultsNames(Zhang_dds_path_deseq) 
# [1] "Intercept"               "path_nonpath_vs_control" "path_path_vs_control"   

## BUILD THE RESULTS OBJECT
# Examining the results object, change alpha to p <0.05, looking at object metadata
# use mcols to look at metadata for each table

mcols(Zhang_dds_broked_group_deseq)
mcols(Zhang_dds_deseq)
Zhang_dds_deseq_res <- results(Zhang_dds_deseq, alpha=0.05, name="group_challenge_vs_control")
Zhang_dds_deseq_res_V_alg1 <- results(Zhang_dds_broken_group_deseq, alpha=0.05, name = "group_by_sim_V_aes_V_alg1_V_alg2_vs_control"  )
Zhang_dds_deseq_res_V_tub <- results(Zhang_dds_broken_group_deseq, alpha=0.05, name= "group_by_sim_V_tub_V_ang_vs_control" )
Zhang_dds_deseq_res_LPS <- results(Zhang_dds_broken_group_deseq, alpha=0.05, name= "group_by_sim_LPS_M_lut_vs_control")

# Zhang_dds_path_deseq_res_path_nonpath <- results(Zhang_dds_path_deseq, alpha=0.05, name=  "path_nonpath_vs_control")
# Zhang_dds_path_deseq_res_path_path <- results(Zhang_dds_path_deseq, alpha=0.05, name=  "path_path_vs_control")

head(Zhang_dds_deseq_res) # group challenge vs control 
head(Zhang_dds_deseq_res_V_alg1) #  group by sim Vibrio vs control 
head(Zhang_dds_deseq_res_V_tub) # group by sim LPS M lut Vtub vs control 
head(Zhang_dds_deseq_res_LPS) # group by sim V aes v alg2 vs control 
head(Zhang_dds_path_deseq_res_path_nonpath) # : path nonpath vs control 
head(Zhang_dds_path_deseq_res_path_path) #  path path vs control 

### Perform LFC Shrinkage with apeglm
## NOTES 
# Before plotting we need to apply an LFC shrinkage, set the coef as the specific comparison in the ResultsNames function of the deseq object
# Issue that the specific comparisons I set in my results formulas are not available in the ResultsNames coef list.

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

# Notes on setting up coefficients for apeglm, https://support.bioconductor.org/p/115435/ , https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#extended-section-on-shrinkage-estimators
# Although apeglm cannot be used with contrast, we note that many designs can be easily rearranged such that what was a contrast becomes its own coefficient.
# In this case, the dispersion does not have to be estimated again, as the designs are equivalent, up to the meaning of the coefficients. 
# Instead, one need only run nbinomWaldTest to re-estimate MLE coefficients – these are necessary for apeglm – and then run lfcShrink specifying 
# the coefficient of interest in resultsNames(dds)
# The user would for example, either change the levels of dds$condition or replace the design using design(dds)<-, then run nbinomWaldTest followed by lfcShrink

# For each LFCshrink I can pass to it my res object for each so that I can keep my alpha setting at 0.05. Doing this procedure will 
# keep the p-values and padj from the results() call, and simply update the LFCs so they are posterior estimates.

## DECISION: USE SAME RES OBJECT TO KEEP ALPHA ADJUSTMENT, and use LFCShrink apeglm

#Zhang_dds_deseq_res_LFC <- lfcShrink(Zhang_dds_deseq, coef="group_challenge_vs_control", type= "apeglm", res=Zhang_dds_deseq_res)
# Review results object summary
#summary(Zhang_dds_deseq_res_LFC) # SHOWS NUMBER OF SIGNIFICANT GENES
# MSTRG removed 
#out of 32394 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 649, 2%
#LFC < 0 (down)     : 146, 0.45%
#outliers [1]       : 0, 0%
#low counts [2]     : 4786, 15%
#(mean count < 3)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

Zhang_dds_deseq_res_V_alg1_LFC <- lfcShrink(Zhang_dds_deseq, coef="group_by_sim_V_aes_V_alg1_V_alg2_vs_control" , type= "apeglm", res=Zhang_dds_deseq_res_V_alg1)
# Review results object summary
summary(Zhang_dds_deseq_res_V_alg1_LFC)
#out of 33418 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 1156, 3.5%
#LFC < 0 (down)     : 172, 0.51%
#outliers [1]       : 0, 0%
#low counts [2]     : 6479, 19%
#(mean count < 5)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

Zhang_dds_deseq_res_V_tub_LFC <- lfcShrink(Zhang_dds_deseq, coef="group_by_sim_V_tub_V_ang_vs_control" , type= "apeglm", res=Zhang_dds_deseq_res_V_tub)
summary(Zhang_dds_deseq_res_V_tub_LFC)
#out of 33418 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 666, 2%
#LFC < 0 (down)     : 216, 0.65%
#outliers [1]       : 0, 0%
#low counts [2]     : 11014, 33%
#(mean count < 12)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

Zhang_dds_deseq_res_LFC_LPS <- lfcShrink(Zhang_dds_deseq, coef="group_by_sim_LPS_M_lut_vs_control", type= "apeglm", res=Zhang_dds_deseq_res_LPS)
summary(Zhang_dds_deseq_res_LFC_LPS)
#out of 33418 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 612, 1.8%
#LFC < 0 (down)     : 210, 0.63%
#outliers [1]       : 0, 0%
#low counts [2]     : 11662, 35%
#(mean count < 13)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

#Zhang_dds_path_deseq_res_path_nonpath_LFC <- lfcShrink(Zhang_dds_path_deseq, coef="path_nonpath_vs_control", type= "apeglm", res=Zhang_dds_path_deseq_res_path_nonpath)
#summary(Zhang_dds_path_deseq_res_path_nonpath_LFC)
#out of 33418 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 446, 1.3%
#LFC < 0 (down)     : 137, 0.41%
#outliers [1]       : 1598, 4.8%
#low counts [2]     : 3888, 12%
#(mean count < 3)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

#Zhang_dds_path_deseq_res_path_path_LFC <- lfcShrink(Zhang_dds_path_deseq, coef="path_path_vs_control", type= "apeglm", res=Zhang_dds_path_deseq_res_path_path)
#summary(Zhang_dds_path_deseq_res_path_path_LFC)
#out of 33418 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 289, 0.86%
#LFC < 0 (down)     : 228, 0.68%
#outliers [1]       : 1598, 4.8%
#low counts [2]     : 3888, 12%
#(mean count < 3)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

## EXPLORATORY PLOTTING OF RESULTS 
## MA Plotting
plotMA(Zhang_dds_deseq_res_LFC, ylim = c(-5, 5))
## Histogram of P values 
# exclude genes with very small counts to avoid spikes and plot using the LFCshrinkage
hist(Zhang_dds_deseq_res_LFC$padj[Zhang_dds_deseq_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

### Subsetting Significant Genes by padj < 0.05
# again, only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
Zhang_dds_deseq_res_LFC_sig <- subset(Zhang_dds_deseq_res_LFC, padj < 0.05)
Zhang_dds_deseq_res_LFC_sig$transcript_id <- row.names(Zhang_dds_deseq_res_LFC_sig)
Zhang_dds_deseq_res_LFC_sig <- as.data.frame(Zhang_dds_deseq_res_LFC_sig)
nrow(Zhang_dds_deseq_res_LFC_sig)  #795 with MSTRG removed 

Zhang_dds_deseq_res_V_alg1_LFC_sig <- subset(Zhang_dds_deseq_res_V_alg1_LFC, padj < 0.05)
Zhang_dds_deseq_res_V_alg1_LFC_sig $transcript_id <- row.names(Zhang_dds_deseq_res_V_alg1_LFC_sig )
Zhang_dds_deseq_res_V_alg1_LFC_sig  <- as.data.frame(Zhang_dds_deseq_res_V_alg1_LFC_sig )
nrow(Zhang_dds_deseq_res_V_alg1_LFC_sig )  #1328

Zhang_dds_deseq_res_V_tub_LFC_sig <- subset(Zhang_dds_deseq_res_V_tub_LFC, padj < 0.05)
Zhang_dds_deseq_res_V_tub_LFC_sig  $transcript_id <- row.names(Zhang_dds_deseq_res_V_tub_LFC_sig  )
Zhang_dds_deseq_res_V_tub_LFC_sig   <- as.data.frame(Zhang_dds_deseq_res_V_tub_LFC_sig )
nrow(Zhang_dds_deseq_res_V_tub_LFC_sig  )  #882

Zhang_dds_deseq_res_LFC_LPS_sig <- subset(Zhang_dds_deseq_res_LFC_LPS, padj < 0.05)
Zhang_dds_deseq_res_LFC_LPS_sig $transcript_id <- row.names(Zhang_dds_deseq_res_LFC_LPS_sig )
Zhang_dds_deseq_res_LFC_LPS_sig  <- as.data.frame(Zhang_dds_deseq_res_LFC_LPS_sig)
nrow(Zhang_dds_deseq_res_LFC_LPS_sig)  #822

# Zhang_dds_path_deseq_res_path_nonpath_LFC_sig <-subset(Zhang_dds_path_deseq_res_path_nonpath_LFC, padj < 0.05)
# Zhang_dds_path_deseq_res_path_nonpath_LFC_sig$transcript_id <- row.names(Zhang_dds_path_deseq_res_path_nonpath_LFC_sig )
# Zhang_dds_path_deseq_res_path_nonpath_LFC_sig <- as.data.frame(Zhang_dds_path_deseq_res_path_nonpath_LFC_sig)
# nrow(Zhang_dds_path_deseq_res_path_nonpath_LFC_sig )  #583
# 
# Zhang_dds_path_deseq_res_path_path_LFC_sig  <-subset(Zhang_dds_path_deseq_res_path_path_LFC , padj < 0.05)
# Zhang_dds_path_deseq_res_path_path_LFC_sig $transcript_id <- row.names(Zhang_dds_path_deseq_res_path_path_LFC_sig  )
# Zhang_dds_path_deseq_res_path_path_LFC_sig  <- as.data.frame(Zhang_dds_path_deseq_res_path_path_LFC_sig )
# nrow(Zhang_dds_path_deseq_res_path_path_LFC_sig  )  #517

### GENE CLUSTERING ANALYSIS HEATMAPS  
# Extract genes with the highest variance across samples for each comparison using either vst or rlog transformed data
# This heatmap rather than plotting absolute expression strength plot the amount by which each gene deviates in a specific sample from the gene’s average across all samples. 
# example codes from RNAseq workflow: https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#other-comparisons
# topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 40)
# mat  <- assay(vsd)[ topVarGenes, ]
# mat  <- mat - rowMeans(mat)
# anno <- as.data.frame(colData(vsd)[, c("cell","dex")])
# pheatmap(mat, annotation_col = anno)

topVarGenes_Zhang_dds_rlog <-  head(order(rowVars(assay(Zhang_dds_rlog )), decreasing = TRUE), 100)
family_Zhang_mat <- assay(Zhang_dds_rlog)[topVarGenes_Zhang_dds_rlog,]
family_Zhang_mat <- family_Zhang_mat - rowMeans(family_Zhang_mat)
family_Zhang_anno <- as.data.frame(colData(Zhang_dds_rlog)[, c("condition","group")])
family_Zhang_heatmap <- pheatmap(family_Zhang_mat , annotation_col = family_Zhang_anno)
head(family_Zhang_mat)
# With MSTRG removed, the control PBS and control control cluster together, then the LPS, M lut and V alg1 cluster
  # and the other Vibrio challenges cluster 

topVarGenes_Zhang_dds_broken_rlog <-  head(order(rowVars(assay(Zhang_dds_broked_group_rlog)), decreasing = TRUE), 200)
family_Zhang_broken_mat <- assay(Zhang_dds_broked_group_rlog)[topVarGenes_Zhang_dds_broken_rlog,]
family_Zhang_broken_mat <- family_Zhang_broken_mat - rowMeans(family_Zhang_broken_mat)
family_Zhang_broken_anno <- as.data.frame(colData(Zhang_dds_broked_group_rlog)[, c("condition","group_by_sim")])
family_Zhang_broken_heatmap <- pheatmap(family_Zhang_broken_mat , annotation_col = family_Zhang_broken_anno)
head(family_Zhang_broken_mat)
# the V_aes and V_alg2 cluster is supported, the LPS V alg 1 and M lut cluster, V ang and V tub cluster, though LPS shares the most with control and PBS
# I'm going to re-do the groups so that V_aes, V_alg2 together, LPS, PBS and control together, and V alg 1 and Mtub together
# or can do the non pathogenic to oyster V aes, V alg 1 and  V alg 2 together, and the V. ang and V. tub together, and LPS, PBS, and control together
# with 200 genes: PBS, V alg1, M lut cluster,control and LPS cluster, V aes and V alg 2 cluster, V ang and V tub clustered

topVarGenes_Zhang_dds_path_rlog <-  head(order(rowVars(assay(Zhang_dds_path_rlog )), decreasing = TRUE), 200)
family_Zhang_path_mat <- assay(Zhang_dds_path_rlog )[topVarGenes_Zhang_dds_path_rlog,]
family_Zhang_path_mat <- family_Zhang_path_mat - rowMeans(family_Zhang_path_mat)
family_Zhang_path_anno <- as.data.frame(colData(Zhang_dds_path_rlog )[, c("condition","path")])
family_Zhang_path_heatmap <- pheatmap(family_Zhang_path_mat , annotation_col = family_Zhang_path_anno)
head(family_Zhang_path_mat)
# with 200 genes: V ang and V tub still clustering, V aes V alg2 clutered, M lut and V alg 1 clustered with PBS

# reorder annotation table to match ordering in heatmap 
family_Zhang_heatmap_reorder <-rownames(family_Zhang_mat[family_Zhang_heatmap$tree_row[["order"]],])
# annotate the row.names
family_Zhang_mat_prot <- as.data.frame(family_Zhang_heatmap_reorder)
colnames(family_Zhang_mat_prot)[1] <- "transcript_id"
family_Zhang_mat_prot_annot <- left_join(family_Zhang_mat_prot, select(C_gig_rtracklayer_transcripts, transcript_id, product, gene), by = "transcript_id")

# Gene clustering heatmap with only apoptosis genes #
# vector C_gig_apop transcript IDs
C_gig_rtracklayer_apop_product_final_transcript_id <- as.vector(C_gig_rtracklayer_apop_product_final$transcript_id)
# Search original Zhang_counts for apoptosis genes and do rlog on just these
Zhang_counts_apop <- Zhang_counts[row.names(Zhang_counts) %in% C_gig_rtracklayer_apop_product_final_transcript_id,]
nrow(Zhang_counts_apop) #659
head(Zhang_counts_apop)
Zhang_counts_apop_dds <- DESeqDataSetFromMatrix(countData = Zhang_counts_apop,
                                         colData = Zhang_coldata,
                                         design = ~time + group_by_sim) # add time to control for injection and time effect
# Prefiltering the data and running rlog
Zhang_counts_apop_dds <- Zhang_counts_apop_dds[ rowSums(counts(Zhang_counts_apop_dds)) > 10, ]
Zhang_counts_apop_rlog <- rlog(Zhang_counts_apop_dds, blind=TRUE)

## PCA plot of rlog transformed counts for apoptosis
plotPCA(Zhang_counts_apop_rlog , intgroup=c("group", "condition"))

# heatmap of all apoptosis genes 
Zhang_counts_apop_assay <-  assay(Zhang_counts_apop_rlog)[,]
Zhang_counts_apop_assay_mat <- Zhang_counts_apop_assay - rowMeans(Zhang_counts_apop_assay)
Zhang_counts_apop_assay_anno <- as.data.frame(colData(Zhang_counts_apop_rlog )[, c("condition","group_by_sim")])
Zhang_counts_apop_assay_heatmap <- pheatmap(Zhang_counts_apop_assay_mat  , annotation_col = Zhang_counts_apop_assay_anno)
head(Zhang_counts_apop_assay_mat )
# control and PBS grouping, V_aes and V_ alg2 grouping, V_ang V alg 1 and V tub clustering

# heatmap of most variable apoptosis genes (this selects genes with the greatest variance in the sample)
topVarGenes_Zhang_counts_apop_assay <-  head(order(rowVars(assay(Zhang_counts_apop_rlog)), decreasing = TRUE), 100) 
top_Var_Zhang_counts_apop_assay_mat<- assay(Zhang_counts_apop_rlog )[topVarGenes_Zhang_counts_apop_assay,]
top_Var_Zhang_counts_apop_assay_mat <- top_Var_Zhang_counts_apop_assay_mat - rowMeans(top_Var_Zhang_counts_apop_assay_mat)
top_Var_Zhang_counts_apop_assay_anno <- as.data.frame(colData(Zhang_counts_apop_rlog )[, c("condition","group_by_sim")])
top_Var_Zhang_counts_apop_assay_heatmap <- pheatmap(top_Var_Zhang_counts_apop_assay_mat  , annotation_col = top_Var_Zhang_counts_apop_assay_anno)
head(top_Var_Zhang_counts_apop_assay_mat )
# same grouping as above with top 200, top 100 changes grouping M lut LPS V alg1 and V alg2 group together in a cluster, while V tub, V ang and V aes group with control

# annotate the top 100 most variable genes  
# reorder annotation table to match ordering in heatmap 
top_Var_Zhang_counts_apop_assay_heatmap_reorder <-rownames(top_Var_Zhang_counts_apop_assay_mat[top_Var_Zhang_counts_apop_assay_heatmap $tree_row[["order"]],])
# annotate the row.names
top_Var_Zhang_counts_apop_assay_prot <- as.data.frame(top_Var_Zhang_counts_apop_assay_heatmap_reorder)
colnames(top_Var_Zhang_counts_apop_assay_prot)[1] <- "transcript_id"
top_Var_Zhang_counts_apop_assay_prot_annot <- left_join(top_Var_Zhang_counts_apop_assay_prot, select(C_gig_rtracklayer_transcripts, transcript_id, product, gene), by = "transcript_id")
#isolate interesting clusters

### Extract list of significant Apoptosis Genes (not less than or greater than 1 LFC) using merge
Zhang_dds_deseq_res_LFC_sig_APOP <- merge(Zhang_dds_deseq_res_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")
Zhang_dds_deseq_res_LFC_sig_APOP_arranged <- arrange(Zhang_dds_deseq_res_LFC_sig_APOP, -log2FoldChange) 
nrow(Zhang_dds_deseq_res_LFC_sig_APOP ) # 13
# only 13 apoptotis genes when looking at the results object as a whole with only control vs. challenged

Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP <- merge(Zhang_dds_deseq_res_V_alg1_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")
Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP_arranged <- arrange(Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP , -log2FoldChange) 
nrow(Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP) #19

Zhang_dds_deseq_res_V_tub_LFC_sig_APOP <- merge(Zhang_dds_deseq_res_V_tub_LFC_sig , C_gig_rtracklayer_apop_product_final, by = "transcript_id")
Zhang_dds_deseq_res_V_tub_LFC_sig_APOP_arranged <- arrange(Zhang_dds_deseq_res_V_tub_LFC_sig_APOP  , -log2FoldChange) 
nrow(Zhang_dds_deseq_res_V_tub_LFC_sig_APOP) # 16

Zhang_dds_deseq_res_LFC_LPS_sig_APOP <- merge(Zhang_dds_deseq_res_LFC_LPS_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")
Zhang_dds_deseq_res_LFC_LPS_sig_APOP_arranged <- arrange(Zhang_dds_deseq_res_LFC_LPS_sig_APOP , -log2FoldChange) 
nrow(Zhang_dds_deseq_res_LFC_LPS_sig_APOP )  #14

# Zhang_dds_path_deseq_res_path_nonpath_LFC_sig_APOP <- merge(Zhang_dds_path_deseq_res_path_nonpath_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")
# Zhang_dds_path_deseq_res_path_nonpath_LFC_sig_APOP_arranged <- arrange(Zhang_dds_path_deseq_res_path_nonpath_LFC_sig_APOP , -log2FoldChange) 
# nrow(Zhang_dds_path_deseq_res_path_nonpath_LFC_sig_APOP ) # 9
# 
# Zhang_dds_path_deseq_res_path_path_LFC_sig_APOP <- merge(Zhang_dds_path_deseq_res_path_path_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")
# Zhang_dds_path_deseq_res_path_path_LFC_sig_APOP_arranged <- arrange(Zhang_dds_path_deseq_res_path_path_LFC_sig_APOP , -log2FoldChange) 
# nrow(Zhang_dds_path_deseq_res_path_path_LFC_sig_APOP ) # 6

# Compare apoptosis genes between group_by_sim groups
Zhang_dds_deseq_res_V_alg1_APOP_short <- Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP[,c(1:6,28)]
Zhang_dds_deseq_res_V_alg1_APOP_short$group_by_sim <- "V_aes_V_alg1_V_alg2"
Zhang_dds_deseq_res_V_tub_APOP_short <- Zhang_dds_deseq_res_V_tub_LFC_sig_APOP[,c(1:6,28)] 
Zhang_dds_deseq_res_V_tub_APOP_short$group_by_sim <- "V_tub_V_ang"
Zhang_dds_deseq_res_LPS_APOP_short <- Zhang_dds_deseq_res_LFC_LPS_sig_APOP[,c(1:6,28)] 
Zhang_dds_deseq_res_LPS_APOP_short$group_by_sim <- "LPS_M_lut"

# combine data frames 
Zhang_upset_all_sig_APOP <- rbind(Zhang_dds_deseq_res_V_alg1_APOP_short[,c("product","group_by_sim","log2FoldChange")],
                                  Zhang_dds_deseq_res_V_tub_APOP_short[,c("product","group_by_sim","log2FoldChange")],
                                  Zhang_dds_deseq_res_LPS_APOP_short[,c("product","group_by_sim","log2FoldChange")] )


# Convert into wide format using reshape
Zhang_upset_all_sig_APOP_tally <- Zhang_upset_all_sig_APOP %>% group_by(product) %>% tally() 
Zhang_upset_all_sig_APOP_upset <- Zhang_upset_all_sig_APOP %>% group_by(product) %>% mutate(value=1) %>% spread(group_by_sim, value, fill =0 )
Zhang_upset_all_sig_APOP_upset <- as.matrix(Zhang_upset_all_sig_APOP_upset)

# Make plot
Zhang_full_LFC_plot <- ggplot(Zhang_upset_all_sig_APOP, aes(x=product,y=log2FoldChange, fill=group_by_sim )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()

Zhang_dds_deseq_res_V_alg1_APOP_short_plot <- ggplot(Zhang_dds_deseq_res_V_alg1_APOP_short, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Zhang V. alg, V. aes vs Control") +
  ylab("Log2 Fold Change")
Zhang_dds_deseq_res_V_tub_APOP_short_plot <- ggplot(Zhang_dds_deseq_res_V_tub_APOP_short, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("V. tub V. ang vs Control") +
  ylab("Log2 Fold Change")
Zhang_dds_deseq_res_LPS_APOP_short_plot <- ggplot(Zhang_dds_deseq_res_LPS_APOP_short, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Zhang LPS and M Lut LFC vs Control") +
  ylab("Log2 Fold Change")

ggarrange(Zhang_dds_deseq_res_V_alg1_APOP_short_plot , Zhang_dds_deseq_res_V_tub_APOP_short_plot,Zhang_dds_deseq_res_LPS_APOP_short_plot)


## compe back to finishing this analysis after I have done the other dataset analyses

#### RUBIO VIBRIO TRANSCRIPTOME ANALYSIS ####

## LOAD DATA
Rubio_counts <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/Rubio_transcript_count_matrix.csv", header=TRUE,
                         row.names = "transcript_id")
head(Rubio_counts)
colnames(Rubio_counts)
# colnames all have "_1", remove this. It was an artifact of the PE sample names
colnames(Rubio_counts) <- sub('\\_[^_]+$', '', colnames(Rubio_counts))
colnames(Rubio_counts)

# remove MSTRG novel transcript lines (can assess these later if necessary)
Rubio_counts <- Rubio_counts[!grepl("MSTRG", row.names(Rubio_counts)),]

# Cute the "rna-" from the beginning of rownames
remove_rna = function(x){
  return(gsub("rna-","",x))
}
row.names(Rubio_counts) <- remove_rna(row.names(Rubio_counts))
head(Rubio_counts)

#Load in sample metadata
Rubio_coldata <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/Rubio_coldata.csv", row.names = 1 )
View(Rubio_coldata)  
nrow(Rubio_coldata) 

# Make sure the columns of the count matrix and rows of the column data (sample metadata) are in the same order. 
all(rownames(Rubio_coldata) %in% colnames(Rubio_counts ))  #Should return TRUE
# returns TRUE
all(colnames(Rubio_counts ) %in% rownames(Rubio_coldata))  
# returns TRUE
all(rownames(Rubio_coldata) == colnames(Rubio_counts ))    # should return TRUE
# returns FALSE

# Fix the order
Rubio_counts <-Rubio_counts[,row.names(Rubio_coldata)]
row.names(Rubio_coldata)

all(rownames(Rubio_coldata) %in% colnames(Rubio_counts ))  #Should return TRUE
# returns TRUE
all(colnames(Rubio_counts ) %in% rownames(Rubio_coldata))  
# returns TRUE
all(rownames(Rubio_coldata) == colnames(Rubio_counts ))    # should return TRUE
# returns TRUE

### DATA QC PCA PLOT 
# rlog transform data is recommended over vst for small data sets 
# PCA plots of data (https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/02_Preprocessing_Data.nb.html#count-distribution-boxplots)
Rubio_counts_matrix <- as.matrix(Rubio_counts)
Rubiorlogcounts <- rlog(Rubio_counts_matrix, blind =TRUE)
  #-- note: fitType='parametric', but the dispersion trend was not well captured by the
  #function: y = a/x + b, and a local regression fit was automatically substituted.
  #specify fitType='local' or 'mean' to avoid this message next time.

# run PCA
pcRubio <- prcomp(t(Rubiorlogcounts))

# Plot PCA
autoplot(pcRubio,
         data = Rubio_coldata, 
         colour="Condition", 
         size=5) # PCA axes explain very little of the variation! (6 and 7 percent). There is a high amount of variablility between samples. 
        # control untreated has the least variability however

# Plot PCA 2 and 3 for comparison
autoplot(pcRubio,
         data = Rubio_coldata, 
         colour = "Condition", 
         size = 5,
         x = 2,
         y = 3) # no new trends emerge
# Plot PCA 4 and 4 for comparison
autoplot(pcRubio,
         data = Rubio_coldata, 
         colour = "Condition", 
         size = 5,
         x = 3,
         y = 4) # extremely spread, still high variability 

## MAKE DESEQ DATA SET FROM MATRIX
# This object specifies the count data and metadata you will work with. The design piece is critical.
# Correct for batch effects if necessary in this original formula: see this thread https://support.bioconductor.org/p/121408/
# do not correct counts using the removeBatchEffects from limma based on thread above 

## Check levels 
# It is prefered in R that the first level of a factor be the reference level for comparison
# (e.g. control, or untreated samples), so we can relevel the factor like so
# Check factor levels, set it so that comparison group is the first
levels(Rubio_coldata$Condition) # "Control_anesthesis" "Control_untreated"  "Vcrass_J2_8"        "Vcrass_J2_9"        "Vtasma_LGP32"       "Vtasma_LMG20012T"  
Rubio_coldata$Condition <- factor(Rubio_coldata$Condition , levels = c("Control_untreated","Control_anesthesis", "Vcrass_J2_8", "Vcrass_J2_9","Vtasma_LGP32", "Vtasma_LMG20012T"  ))
levels(Rubio_coldata$Condition)
levels(Rubio_coldata$Group) # "Control"      "Non_virulent" "Virulent"  

## Creating three here so I can compare the results
Rubio_dds <- DESeqDataSetFromMatrix(countData = Rubio_counts,
                                    colData = Rubio_coldata,
                                    design = ~ Condition) 


## Prefiltering the data
# Data prefiltering helps decrease the size of the data set and get rid of
# rows with no data or very minimal data (<10). Apply a minimal filtering here as more stringent filtering will be applied later
Rubio_dds <- Rubio_dds [ rowSums(counts(Rubio_dds )) > 10, ]

## DATA TRANSFORMATION AND VISUALIZATION
# Assess sample clustering after setting initial formula for comparison
Rubio_dds_rlog <- rlog(Rubio_dds, blind = TRUE) # keep blind = true before deseq function has been run
#-- note: fitType='parametric', but the dispersion trend was not well captured by the
#function: y = a/x + b, and a local regression fit was automatically substituted.
#specify fitType='local' or 'mean' to avoid this message next time.

## PCA plot visualization of individuals in the family 
plotPCA(Rubio_dds_rlog, intgroup=c("Sample", "Condition"))
# Still extremely high variation 

### DIFFERENTIAL EXPRESSION ANALYSIS
# run pipeline with single command because the formula has already been specified
# Steps: estimation of size factors (controlling for differences in the sequencing depth of the samples), 
# the estimation of dispersion values for each gene,and fitting a generalized linear model.
Rubio_dds_deseq <- DESeq(Rubio_dds) 
#-- note: fitType='parametric', but the dispersion trend was not well captured by the
#function: y = a/x + b, and a local regression fit was automatically substituted.
#specify fitType='local' or 'mean' to avoid this message next time.

## Check the resultsNames object of each to look at the available coefficients for use in lfcShrink command
resultsNames(Rubio_dds_deseq) # [1] "Intercept", "Condition_Control_anesthesis_vs_Control_untreated", "Condition_Vcrass_J2_8_vs_Control_untreated"       
# [4] "Condition_Vcrass_J2_9_vs_Control_untreated"        "Condition_Vtasma_LGP32_vs_Control_untreated"       "Condition_Vtasma_LMG20012T_vs_Control_untreated"   


## BUILD THE RESULTS OBJECT
# Examining the results object, change alpha to p <0.05, looking at object metadata
# use mcols to look at metadata for each table

mcols(Rubio_dds_deseq)
Rubio_dds_deseq_J2_8_res <- results(Rubio_dds_deseq, alpha=0.05, name="Condition_Vcrass_J2_8_vs_Control_untreated" )
Rubio_dds_deseq_J2_9_res <- results(Rubio_dds_deseq, alpha=0.05, name="Condition_Vcrass_J2_9_vs_Control_untreated"  )
Rubio_dds_deseq_LGP32_res <- results(Rubio_dds_deseq, alpha=0.05, name="Condition_Vtasma_LGP32_vs_Control_untreated"   )
Rubio_dds_deseq_LMG20012T_res <- results(Rubio_dds_deseq, alpha=0.05, name="Condition_Vtasma_LMG20012T_vs_Control_untreated"  )

head(Rubio_dds_deseq_J2_8_res) # Condition Vcrass J2 8 vs Control untreated 
head(Rubio_dds_deseq_J2_9_res) #  Condition Vcrass J2 9 vs Control untreated
head(Rubio_dds_deseq_LGP32_res) # Condition Vtasma LGP32 vs Control untreated 
head(Rubio_dds_deseq_LMG20012T_res) # Condition Vtasma LMG20012T vs Control untreated 

### Perform LFC Shrinkage with apeglm
## NOTES 
# Before plotting we need to apply an LFC shrinkage, set the coef as the specific comparison in the ResultsNames function of the deseq object
# Issue that the specific comparisons I set in my results formulas are not available in the ResultsNames coef list.
# More detailed notes about LFC Shrinkage are in the code for the Zhang Vibrio

## DECISION: USE SAME RES OBJECT TO KEEP ALPHA ADJUSTMENT, and use LFCShrink apeglm
Rubio_dds_deseq_J2_8_res_LFC<- lfcShrink(Rubio_dds_deseq, coef="Condition_Vcrass_J2_8_vs_Control_untreated", type="apeglm", res= Rubio_dds_deseq_J2_8_res)
Rubio_dds_deseq_J2_9_res_LFC  <- lfcShrink(Rubio_dds_deseq, coef="Condition_Vcrass_J2_9_vs_Control_untreated" , type= "apeglm", res= Rubio_dds_deseq_J2_9_res   )
Rubio_dds_deseq_LGP32_res_LFC<- lfcShrink(Rubio_dds_deseq, coef="Condition_Vtasma_LGP32_vs_Control_untreated" , type= "apeglm", res= Rubio_dds_deseq_LGP32_res   )
Rubio_dds_deseq_LMG20012T_res_LFC<- lfcShrink(Rubio_dds_deseq, coef="Condition_Vtasma_LMG20012T_vs_Control_untreated", type= "apeglm", res= Rubio_dds_deseq_LMG20012T_res )

## EXPLORATORY PLOTTING OF RESULTS 
## MA Plotting
plotMA(Rubio_dds_deseq_J2_8_res_LFC, ylim = c(-5, 5))
plotMA(Rubio_dds_deseq_J2_9_res_LFC, ylim = c(-5, 5))
plotMA(Rubio_dds_deseq_LGP32_res_LFC, ylim = c(-5, 5))
plotMA(Rubio_dds_deseq_LMG20012T_res_LFC, ylim = c(-5, 5))
## Histogram of P values 
# exclude genes with very small counts to avoid spikes and plot using the LFCshrinkage
hist(Rubio_dds_deseq_J2_8_res_LFC$padj[Rubio_dds_deseq_J2_8_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
hist(Rubio_dds_deseq_J2_9_res_LFC$padj[Rubio_dds_deseq_J2_9_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
hist(Rubio_dds_deseq_LGP32_res_LFC$padj[Rubio_dds_deseq_LGP32_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
hist(Rubio_dds_deseq_LMG20012T_res_LFC$padj[Rubio_dds_deseq_LMG20012T_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

### Subsetting Significant Genes by padj < 0.05
# again, only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
Rubio_dds_deseq_J2_8_res_LFC_sig <-  subset(Rubio_dds_deseq_J2_8_res_LFC , padj < 0.05)
Rubio_dds_deseq_J2_9_res_LFC_sig <-  subset(Rubio_dds_deseq_J2_9_res_LFC , padj < 0.05)
Rubio_dds_deseq_LGP32_res_LFC_sig <-  subset(Rubio_dds_deseq_LGP32_res_LFC , padj < 0.05)
Rubio_dds_deseq_LMG20012T_res_LFC_sig <-  subset(Rubio_dds_deseq_LMG20012T_res_LFC, padj < 0.05)

Rubio_dds_deseq_J2_8_res_LFC_sig$transcript_id <- row.names(Rubio_dds_deseq_J2_8_res_LFC_sig)
Rubio_dds_deseq_J2_9_res_LFC_sig$transcript_id <- row.names(Rubio_dds_deseq_J2_9_res_LFC_sig)
Rubio_dds_deseq_LGP32_res_LFC_sig$transcript_id <- row.names(Rubio_dds_deseq_LGP32_res_LFC_sig)
Rubio_dds_deseq_LMG20012T_res_LFC_sig$transcript_id <- row.names(Rubio_dds_deseq_LMG20012T_res_LFC_sig)

Rubio_dds_deseq_J2_8_res_LFC_sig <- as.data.frame(Rubio_dds_deseq_J2_8_res_LFC_sig)
Rubio_dds_deseq_J2_9_res_LFC_sig <- as.data.frame(Rubio_dds_deseq_J2_9_res_LFC_sig)
Rubio_dds_deseq_LGP32_res_LFC_sig <- as.data.frame(Rubio_dds_deseq_LGP32_res_LFC_sig)
Rubio_dds_deseq_LMG20012T_res_LFC_sig<- as.data.frame(Rubio_dds_deseq_LMG20012T_res_LFC_sig)

nrow(Rubio_dds_deseq_J2_8_res_LFC_sig) # 3532
nrow(Rubio_dds_deseq_J2_9_res_LFC_sig) #3719
nrow(Rubio_dds_deseq_LGP32_res_LFC_sig) # 3783
nrow(Rubio_dds_deseq_LMG20012T_res_LFC_sig) # 3571

### GENE CLUSTERING ANALYSIS HEATMAPS  
# Extract genes with the highest variance across samples for each comparison using either vst or rlog transformed data
# This heatmap rather than plotting absolute expression strength plot the amount by which each gene deviates in a specific sample from the gene’s average across all samples. 
# example codes from RNAseq workflow: https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#other-comparisons

Rubio_dds_deseq_J2_8_res_LFC_sig_rlog <-  head(order(rowVars(assay(Rubio_dds_rlog )), decreasing = TRUE), 200)
family_Rubio_broken_mat <- assay(Rubio_dds_rlog)[Rubio_dds_deseq_J2_8_res_LFC_sig_rlog,]
family_Rubio_broken_mat <- family_Rubio_broken_mat - rowMeans(family_Rubio_broken_mat)
family_Rubio_broken_anno <- as.data.frame(colData(Rubio_dds_rlog)[, c("Condition","Sample")])
family_Rubio_broken_heatmap <- pheatmap(family_Rubio_broken_mat , annotation_col = family_Rubio_broken_anno)
head(family_Rubio_broken_mat)
# Control untreated group together, two samples in Vtasm LPG32 cluster, two V2crass J2_8 cluster, but all other samples do not have a clear clustering pattern for these 
  # most variable genes

# reorder annotation table to match ordering in heatmap 
family_Rubio_broken_heatmap_reorder <-rownames(family_Rubio_broken_mat[family_Rubio_broken_heatmap$tree_row[["order"]],])
# annotate the row.names
family_Rubio_broken_mat_prot <- as.data.frame(family_Rubio_broken_heatmap_reorder )
colnames(family_Rubio_broken_mat_prot)[1] <- "transcript_id"
family_Rubio_broken_mat_prot_annot <- left_join(family_Rubio_broken_mat_prot, select(C_gig_rtracklayer_transcripts, transcript_id, product, gene), by = "transcript_id")
# transcription factor AP1 is the most variable gene across all samples 

# Gene clustering heatmap with only apoptosis genes #
# vector C_gig_apop transcript IDs
C_gig_rtracklayer_apop_product_final_transcript_id 
# Search original Rubio_counts for apoptosis genes and do rlog on just these
Rubio_counts_apop <- Rubio_counts[row.names(Rubio_counts) %in% C_gig_rtracklayer_apop_product_final_transcript_id,]
nrow(Rubio_counts_apop) #659
head(Rubio_counts_apop)
Rubio_counts_apop_dds <- DESeqDataSetFromMatrix(countData = Rubio_counts_apop,
                                                colData = Rubio_coldata,
                                                design = ~ Condition) 
# Prefiltering the data and running rlog
Rubio_counts_apop_dds <- Rubio_counts_apop_dds[ rowSums(counts(Rubio_counts_apop_dds)) > 10, ]
Rubio_counts_apop_dds_rlog <- rlog(Rubio_counts_apop_dds, blind=TRUE)

## PCA plot of rlog transformed counts for apoptosis
plotPCA(Rubio_counts_apop_dds_rlog, intgroup="Condition") # some clustering by condition

# heatmap of all apoptosis genes 
Rubio_counts_apop_assay <-  assay(Rubio_counts_apop_dds_rlog)[,]
Rubio_counts_apop_assay_mat <- Rubio_counts_apop_assay - rowMeans(Rubio_counts_apop_assay)
Rubio_counts_apop_assay_anno <- as.data.frame(colData(Rubio_counts_apop_dds_rlog )[, c("Condition","Sample")])
Rubio_counts_apop_assay_heatmap <- pheatmap(Rubio_counts_apop_assay_mat  , annotation_col = Rubio_counts_apop_assay_anno)
head(Rubio_counts_apop_assay_mat ) 
# clustering of V crass and control anesthesia

# heatmap of most variable apoptosis genes (this selects genes with the greatest variance in the sample)
topVarGenes_Rubio_counts_apop_assay <-  head(order(rowVars(assay(Rubio_counts_apop_dds_rlog)), decreasing = TRUE), 100) 
top_Var_Rubio_counts_apop_assay_mat<- assay(Rubio_counts_apop_dds_rlog)[topVarGenes_Rubio_counts_apop_assay,]
top_Var_Rubio_counts_apop_assay_mat <- top_Var_Rubio_counts_apop_assay_mat - rowMeans(top_Var_Rubio_counts_apop_assay_mat)
top_Var_Rubio_counts_apop_assay_anno <- as.data.frame(colData(Rubio_counts_apop_dds_rlog)[, c("Condition","Sample")])
top_Var_Rubio_counts_apop_assay_heatmap <- pheatmap(top_Var_Rubio_counts_apop_assay_mat  , annotation_col = top_Var_Rubio_counts_apop_assay_anno)
head(top_Var_Rubio_counts_apop_assay_mat )
# some clustering patterns here in signature

# annotate the top 100 most variable genes  
# reorder annotation table to match ordering in heatmap 
top_Var_Rubio_counts_apop_assay_heatmap_reorder <-rownames(top_Var_Rubio_counts_apop_assay_mat[top_Var_Rubio_counts_apop_assay_heatmap $tree_row[["order"]],])
# annotate the row.names
top_Var_Rubio_counts_apop_assay_prot <- as.data.frame(top_Var_Rubio_counts_apop_assay_heatmap_reorder)
colnames(top_Var_Rubio_counts_apop_assay_prot)[1] <- "transcript_id"
top_Var_Rubio_counts_apop_assay_prot_annot <- left_join(top_Var_Rubio_counts_apop_assay_prot, select(C_gig_rtracklayer_transcripts, transcript_id, product, gene), by = "transcript_id")

### Extract list of significant Apoptosis Genes (not less than or greater than 1 LFC) using merge
Rubio_dds_deseq_J2_8_res_LFC_sig_APOP <-    merge(Rubio_dds_deseq_J2_8_res_LFC_sig , C_gig_rtracklayer_apop_product_final, by = "transcript_id")
Rubio_dds_deseq_J2_9_res_LFC_sig_APOP <-    merge(Rubio_dds_deseq_J2_9_res_LFC_sig , C_gig_rtracklayer_apop_product_final, by = "transcript_id")
Rubio_dds_deseq_LGP32_res_LFC_sig_APOP <-   merge(Rubio_dds_deseq_LGP32_res_LFC_sig , C_gig_rtracklayer_apop_product_final, by = "transcript_id")
Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP<-merge(Rubio_dds_deseq_LMG20012T_res_LFC_sig , C_gig_rtracklayer_apop_product_final, by = "transcript_id")

Rubio_dds_deseq_J2_8_res_LFC_sig_APOP_arranged <-    arrange(Rubio_dds_deseq_J2_8_res_LFC_sig_APOP, -log2FoldChange)
Rubio_dds_deseq_J2_9_res_LFC_sig_APOP_arranged <-    arrange(Rubio_dds_deseq_J2_9_res_LFC_sig_APOP, -log2FoldChange)
Rubio_dds_deseq_LGP32_res_LFC_sig_APOP_arranged <-   arrange(Rubio_dds_deseq_LGP32_res_LFC_sig_APOP, -log2FoldChange)
Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP_arranged<-arrange(Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP, -log2FoldChange)

nrow(Rubio_dds_deseq_J2_8_res_LFC_sig_APOP_arranged ) # 63
nrow(Rubio_dds_deseq_J2_9_res_LFC_sig_APOP_arranged  ) # 66
nrow(Rubio_dds_deseq_LGP32_res_LFC_sig_APOP_arranged  ) #68
nrow(Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP_arranged) # 62

View(Rubio_dds_deseq_J2_8_res_LFC_sig_APOP_arranged)
View(Rubio_dds_deseq_J2_9_res_LFC_sig_APOP_arranged)
View(Rubio_dds_deseq_LGP32_res_LFC_sig_APOP_arranged)
View(Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP_arranged)

# Compare apoptosis genes between group_by_sim groups
Rubio_dds_deseq_J2_8_res_LFC_sig_APOP_short <- Rubio_dds_deseq_J2_8_res_LFC_sig_APOP[,c(1:6,28)]
Rubio_dds_deseq_J2_8_res_LFC_sig_APOP_short$group_by_sim <- "J2-8 non-vir"
Rubio_dds_deseq_J2_9_res_LFC_sig_APOP_short <- Rubio_dds_deseq_J2_9_res_LFC_sig_APOP[,c(1:6,28)] 
Rubio_dds_deseq_J2_9_res_LFC_sig_APOP_short$group_by_sim <- "J2-9 vir"
Rubio_dds_deseq_LGP32_res_LFC_sig_APOP_short <-Rubio_dds_deseq_LGP32_res_LFC_sig_APOP[,c(1:6,28)] 
Rubio_dds_deseq_LGP32_res_LFC_sig_APOP_short$group_by_sim <- "LGP32 vir"
Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP_short <-Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP[,c(1:6,28)] 
Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP_short$group_by_sim <- "LMG20012T non-vir"

# combine data frames 
Rubio_all_sig_APOP <- rbind(Rubio_dds_deseq_J2_8_res_LFC_sig_APOP_short[,c("product","group_by_sim","log2FoldChange")],
                            Rubio_dds_deseq_J2_9_res_LFC_sig_APOP_short[,c("product","group_by_sim","log2FoldChange")],
                            Rubio_dds_deseq_LGP32_res_LFC_sig_APOP_short[,c("product","group_by_sim","log2FoldChange")],
                            Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP_short[,c("product","group_by_sim","log2FoldChange")])
# Make plot
Rubio_full_LFC_plot <- ggplot(Rubio_all_sig_APOP , aes(x=product,y=log2FoldChange, fill=group_by_sim )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()


Rubio_dds_deseq_J2_8_res_LFC_sig_APOP_short_plot <- ggplot(Rubio_dds_deseq_J2_8_res_LFC_sig_APOP_short , aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("J2-8 Non-Virulent vs Control") +
  ylab("Log2 Fold Change")
Rubio_dds_deseq_J2_9_res_LFC_sig_APOP_short_plot <- ggplot(Rubio_dds_deseq_J2_9_res_LFC_sig_APOP_short , aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("J2-9 Virulent vs Control") +
  ylab("Log2 Fold Change")
Rubio_dds_deseq_LGP32_res_LFC_sig_APOP_short_plot <- ggplot(Rubio_dds_deseq_LGP32_res_LFC_sig_APOP_short, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("LGP32 Virulent vs Control") +
  ylab("Log2 Fold Change")
Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP_short  <- ggplot(Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP_short , aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("LMG20012T Non-Virulent vs Control") +
  ylab("Log2 Fold Change")
ggarrange(Zhang_dds_deseq_res_V_alg1_APOP_short_plot , Zhang_dds_deseq_res_V_tub_APOP_short_plot,Zhang_dds_deseq_res_LPS_APOP_short_plot)
 

#### DELORGERIL OSHV1 TRANSCRIPTOME ANALYSIS ####

## LOAD DATA
deLorgeril_counts <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/deLorgeril_transcript_count_matrix.csv", header=TRUE,
                      row.names = "transcript_id")
head(deLorgeril_counts)
colnames(deLorgeril_counts)

# colnames all have "_1", remove this. It was an artifact of the PE sample names
colnames(deLorgeril_counts) <- sub('\\_[^_]+$', '', colnames(deLorgeril_counts))
colnames(deLorgeril_counts)

# remove MSTRG novel transcript lines (can assess these later if necessary)
deLorgeril_counts<- deLorgeril_counts[!grepl("MSTRG", row.names(deLorgeril_counts)),]

# Cute the "rna-" from the beginning of rownames
remove_rna = function(x){
  return(gsub("rna-","",x))
}
row.names(deLorgeril_counts) <- remove_rna(row.names(deLorgeril_counts))
head(deLorgeril_counts)

#Load in sample metadata
deLorgeril_coldata <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/deLorgeril_coldata.csv", row.names = 1 )
View(deLorgeril_coldata)  
nrow(deLorgeril_coldata) 

# Make sure the columns of the count matrix and rows of the column data (sample metadata) are in the same order. 
all(rownames(deLorgeril_coldata) %in% colnames(deLorgeril_counts))  #Should return TRUE
# returns TRUE
all(colnames(deLorgeril_counts ) %in% rownames(deLorgeril_coldata))  
# returns TRUE
all(rownames(deLorgeril_coldata) == colnames(deLorgeril_counts ))  # FALSE
# returns true

# Change order
deLorgeril_counts <- deLorgeril_counts[,rownames(deLorgeril_coldata)]
all(rownames(deLorgeril_coldata) %in% colnames(deLorgeril_counts))  #Should return TRUE
# returns TRUE
all(colnames(deLorgeril_counts) %in% rownames(deLorgeril_coldata))  
# returns TRUE
all(rownames(deLorgeril_coldata) == colnames(deLorgeril_counts ))  # TRUE

# split up counts and coldata into the resistant and sucsceptible families (since comparing families is not what I want)
deLorgeril_Resistant_coldata <- deLorgeril_coldata %>% subset(Condition == "AF21_Resistant" | Condition == "AF21_Resistant_control")
deLorgeril_Resistant_counts <- deLorgeril_counts[,row.names(deLorgeril_Resistant_coldata)]
colnames(deLorgeril_Resistant_counts)

deLorgeril_Susceptible_coldata <- deLorgeril_coldata %>% subset(Condition == "AF11_Susceptible" | Condition == "AF11_Susceptible_control")
deLorgeril_Susceptible_counts <- deLorgeril_counts[,row.names(deLorgeril_Susceptible_coldata)]
colnames(deLorgeril_Susceptible_counts)

### DATA QC PCA PLOT 
# rlog transform data is recommended over vst for small data sets 
# PCA plots of data (https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/02_Preprocessing_Data.nb.html#count-distribution-boxplots)
deLorgeril_Resistant_counts_matrix <- as.matrix(deLorgeril_Resistant_counts)
deLorgeril_Susceptible_counts_matrix <- as.matrix(deLorgeril_Susceptible_counts)

deLorgeril_Resistant_vstcounts <- vst(deLorgeril_Resistant_counts_matrix, blind =TRUE)
deLorgeril_Susceptible_vstcounts <- vst(deLorgeril_Susceptible_counts_matrix, blind =TRUE)

# run PCA
pcdeLorgeril <- prcomp(t(deLorgeril_Resistant_vstcounts))
pcdeLorgeril_susceptible <- prcomp(t(deLorgeril_Susceptible_vstcounts ))

# Plot PCA
autoplot(pcdeLorgeril,
         data = deLorgeril_Resistant_coldata, 
         colour="Condition", 
         size=5)# ~19% of variance explained, family is the main separation
autoplot(pcdeLorgeril,
         data = deLorgeril_Resistant_coldata, 
         colour="Time", 
         size=5) # some clustering by time
autoplot(pcdeLorgeril_susceptible,
         data = deLorgeril_Susceptible_coldata, 
         colour="Condition", 
         size=5) # ~25% of the variance explained 
autoplot(pcdeLorgeril_susceptible,
         data = deLorgeril_Susceptible_coldata, 
         colour="Time", 
         size=5)
# Plot PCA 2 and 3 for comparison
autoplot(pcdeLorgeril,
         data = deLorgeril_Resistant_coldata, 
         colour = "Condition", 
         size = 5,
         x = 2,
         y = 3) # greater clustering by treatment
autoplot(pcdeLorgeril_susceptible,
         data = deLorgeril_Susceptible_coldata, 
         colour = "Time", 
         size = 5,
         x = 2,
         y = 3) # clustering mostly by time 

## MAKE DESEQ DATA SET FROM MATRIX
# This object specifies the count data and metadata you will work with. The design piece is critical.
# Correct for batch effects if necessary in this original formula: see this thread https://support.bioconductor.org/p/121408/
# do not correct counts using the removeBatchEffects from limma based on thread above 

## Check levels 
# It is prefered in R that the first level of a factor be the reference level for comparison
# (e.g. control, or untreated samples), so we can relevel the factor like so
# Check factor levels, set it so that comparison group is the first
levels(deLorgeril_Resistant_coldata$Condition)# "AF11_Susceptible"         "AF11_Susceptible_control" "AF21_Resistant"           "AF21_Resistant_control" 
deLorgeril_Resistant_coldata$Condition <- droplevels(deLorgeril_Resistant_coldata$Condition)
levels(deLorgeril_Resistant_coldata$Condition) #  "AF21_Resistant"         "AF21_Resistant_control"       
deLorgeril_Resistant_coldata$Condition <- factor(deLorgeril_Resistant_coldata$Condition, levels=c("AF21_Resistant_control", "AF21_Resistant"))
levels(deLorgeril_Susceptible_coldata$Condition)  
deLorgeril_Susceptible_coldata$Condition <- droplevels(deLorgeril_Susceptible_coldata$Condition)
levels(deLorgeril_Susceptible_coldata$Condition)   #    "AF11_Susceptible"         "AF11_Susceptible_control"    
deLorgeril_Susceptible_coldata$Condition <- factor(deLorgeril_Susceptible_coldata$Condition, levels=c("AF11_Susceptible_control","AF11_Susceptible" ))


levels(deLorgeril_Resistant_coldata$Time) # "0h"  "12h" "24h" "48h" "60h" "6h"  "72h"
deLorgeril_Resistant_coldata$Time <- factor(deLorgeril_Resistant_coldata$Time, levels=c("0h","6h","12h", "24h", "48h", "60h", "72h" ))
levels(deLorgeril_Susceptible_coldata$Time) # "0h"  "12h" "24h" "48h" "60h" "6h"  "72h"
deLorgeril_Susceptible_coldata$Time <- factor(deLorgeril_Susceptible_coldata$Time, levels=c("0h","6h","12h", "24h", "48h", "60h", "72h" ))

## Creating three here so I can compare the results
deLorgeril_Resistant_dds <- DESeqDataSetFromMatrix(countData = deLorgeril_Resistant_counts,
                                 colData = deLorgeril_Resistant_coldata,
                                 design = ~ Time) 
                                  # with both time and condition included, 'Model matrix not full rank', getting rid of accounting for time

deLorgeril_Susceptible_dds <- DESeqDataSetFromMatrix(countData = deLorgeril_Susceptible_counts,
                                 colData = deLorgeril_Susceptible_coldata,
                                 design = ~ Time)  # with both time and position included, 'Model matrix not full rank'. 
                                #Looking at it this way allows me to look for the acute response

## Prefiltering the data
# Data prefiltering helps decrease the size of the data set and get rid of
# rows with no data or very minimal data (<10). Apply a minimal filtering here as more stringent filtering will be applied later
deLorgeril_Resistant_dds <- deLorgeril_Resistant_dds[ rowSums(counts(deLorgeril_Resistant_dds)) > 10, ]
deLorgeril_Susceptible_dds <- deLorgeril_Susceptible_dds[ rowSums(counts(deLorgeril_Susceptible_dds)) > 10, ]

## DATA TRANSFORMATION AND VISUALIZATION
# Assess sample clustering after setting initial formula for comparison
deLorgeril_Resistant_dds_vst <- vst(deLorgeril_Resistant_dds, blind = TRUE) # keep blind = true before deseq function has been run
deLorgeril_Susceptible_dds_vst <- vst(deLorgeril_Susceptible_dds , blind = TRUE)

## PCA plot visualization of individuals in the family, use vst because greater than 30 samples 
plotPCA(deLorgeril_Resistant_dds_vst , intgroup=c("Time", "Condition")) # highly variable, a bit more of variation explained than before
plotPCA(deLorgeril_Susceptible_dds_vst , intgroup=c("Time", "Condition"))

### DIFFERENTIAL EXPRESSION ANALYSIS
# run pipeline with single command because the formula has already been specified
# Steps: estimation of size factors (controlling for differences in the sequencing depth of the samples), 
# the estimation of dispersion values for each gene,and fitting a generalized linear model.
deLorgeril_Resistant_dds_deseq <- DESeq(deLorgeril_Resistant_dds) 
#-- note: fitType='parametric', but the dispersion trend was not well captured by the
#function: y = a/x + b, and a local regression fit was automatically substituted.
#specify fitType='local' or 'mean' to avoid this message next time.
deLorgeril_Susceptible_dds_deseq <- DESeq(deLorgeril_Susceptible_dds) 

## Check the resultsNames object of each to look at the available coefficients for use in lfcShrink command
resultsNames(deLorgeril_Resistant_dds_deseq ) #"Intercept""Time_6h_vs_0h"  "Time_12h_vs_0h" "Time_24h_vs_0h" "Time_48h_vs_0h" "Time_60h_vs_0h" "Time_72h_vs_0h"
resultsNames(deLorgeril_Susceptible_dds_deseq ) #"Intercept"      "Time_6h_vs_0h"  "Time_12h_vs_0h" "Time_24h_vs_0h" "Time_48h_vs_0h" "Time_60h_vs_0h" "Time_72h_vs_0h"

## BUILD THE RESULTS OBJECT
# Examining the results object, change alpha to p <0.05, looking at object metadata
# use mcols to look at metadata for each table, creating three results objects to look find the acute response. The paper observed acute response around 48hr
mcols(deLorgeril_Resistant_dds_deseq)
deLorgeril_Resistant_dds_res_6 <- results(deLorgeril_Resistant_dds_deseq , alpha=0.05, name= "Time_6h_vs_0h" )
deLorgeril_Resistant_dds_res_12 <- results(deLorgeril_Resistant_dds_deseq , alpha=0.05, name= "Time_12h_vs_0h" )
deLorgeril_Resistant_dds_res_24 <- results(deLorgeril_Resistant_dds_deseq , alpha=0.05, name= "Time_24h_vs_0h" )
deLorgeril_Resistant_dds_res_48 <- results(deLorgeril_Resistant_dds_deseq , alpha=0.05, name= "Time_48h_vs_0h" )
deLorgeril_Resistant_dds_res_60 <- results(deLorgeril_Resistant_dds_deseq , alpha=0.05, name= "Time_60h_vs_0h" )
deLorgeril_Resistant_dds_res_72 <- results(deLorgeril_Resistant_dds_deseq , alpha=0.05, name= "Time_72h_vs_0h" )

mcols(deLorgeril_Susceptible_dds_deseq)
deLorgeril_Susceptible_dds_res_6 <- results(deLorgeril_Susceptible_dds_deseq , alpha=0.05, name= "Time_6h_vs_0h" )
deLorgeril_Susceptible_dds_res_12 <- results(deLorgeril_Susceptible_dds_deseq , alpha=0.05, name= "Time_12h_vs_0h" )
deLorgeril_Susceptible_dds_res_24 <- results(deLorgeril_Susceptible_dds_deseq , alpha=0.05, name= "Time_24h_vs_0h" )
deLorgeril_Susceptible_dds_res_48 <- results(deLorgeril_Susceptible_dds_deseq , alpha=0.05, name= "Time_48h_vs_0h" )
deLorgeril_Susceptible_dds_res_60 <- results(deLorgeril_Susceptible_dds_deseq , alpha=0.05, name= "Time_60h_vs_0h" )
deLorgeril_Susceptible_dds_res_72 <- results(deLorgeril_Susceptible_dds_deseq , alpha=0.05, name= "Time_72h_vs_0h" )

### Perform LFC Shrinkage with apeglm
## NOTES 
# Before plotting we need to apply an LFC shrinkage, set the coef as the specific comparison in the ResultsNames function of the deseq object
# Issue that the specific comparisons I set in my results formulas are not available in the ResultsNames coef list.
# More detailed notes about LFC Shrinkage are in the code for the Zhang Vibrio

## DECISION: USE SAME RES OBJECT TO KEEP ALPHA ADJUSTMENT, and use LFCShrink apeglm
deLorgeril_Resistant_dds_res_6_LFC <- lfcShrink(deLorgeril_Resistant_dds_deseq , coef="Time_6h_vs_0h", type="apeglm",res=deLorgeril_Resistant_dds_res_6)
deLorgeril_Resistant_dds_res_12_LFC <- lfcShrink(deLorgeril_Resistant_dds_deseq , coef="Time_12h_vs_0h", type="apeglm",res=deLorgeril_Resistant_dds_res_12)
deLorgeril_Resistant_dds_res_24_LFC <- lfcShrink(deLorgeril_Resistant_dds_deseq , coef="Time_24h_vs_0h", type="apeglm",res=deLorgeril_Resistant_dds_res_24)
deLorgeril_Resistant_dds_res_48_LFC <- lfcShrink(deLorgeril_Resistant_dds_deseq , coef="Time_48h_vs_0h", type="apeglm",res=deLorgeril_Resistant_dds_res_48)
deLorgeril_Resistant_dds_res_60_LFC <- lfcShrink(deLorgeril_Resistant_dds_deseq , coef="Time_60h_vs_0h", type="apeglm",res=deLorgeril_Resistant_dds_res_60)
deLorgeril_Resistant_dds_res_72_LFC <- lfcShrink(deLorgeril_Resistant_dds_deseq , coef="Time_72h_vs_0h", type="apeglm",res=deLorgeril_Resistant_dds_res_72)

deLorgeril_Susceptible_dds_res_6_LFC <- lfcShrink(deLorgeril_Susceptible_dds_deseq , coef= "Time_6h_vs_0h", type="apeglm", res=deLorgeril_Susceptible_dds_res_6)
deLorgeril_Susceptible_dds_res_12_LFC <- lfcShrink(deLorgeril_Susceptible_dds_deseq , coef= "Time_12h_vs_0h", type="apeglm", res=deLorgeril_Susceptible_dds_res_12)
deLorgeril_Susceptible_dds_res_24_LFC <- lfcShrink(deLorgeril_Susceptible_dds_deseq , coef= "Time_24h_vs_0h", type="apeglm", res=deLorgeril_Susceptible_dds_res_24)
deLorgeril_Susceptible_dds_res_48_LFC <- lfcShrink(deLorgeril_Susceptible_dds_deseq , coef= "Time_48h_vs_0h", type="apeglm", res=deLorgeril_Susceptible_dds_res_48)
deLorgeril_Susceptible_dds_res_60_LFC <- lfcShrink(deLorgeril_Susceptible_dds_deseq , coef= "Time_60h_vs_0h", type="apeglm", res=deLorgeril_Susceptible_dds_res_60)
deLorgeril_Susceptible_dds_res_72_LFC <- lfcShrink(deLorgeril_Susceptible_dds_deseq , coef= "Time_72h_vs_0h", type="apeglm", res=deLorgeril_Susceptible_dds_res_72)

## EXPLORATORY PLOTTING OF RESULTS 
## MA Plotting
plotMA(deLorgeril_Resistant_dds_res_6_LFC , ylim = c(-5, 5))
plotMA(deLorgeril_Resistant_dds_res_12_LFC , ylim = c(-5, 5))
plotMA(deLorgeril_Resistant_dds_res_24_LFC , ylim = c(-5, 5))
plotMA(deLorgeril_Resistant_dds_res_48_LFC , ylim = c(-5, 5))
plotMA(deLorgeril_Resistant_dds_res_60_LFC , ylim = c(-5, 5)) # 60hr has more significant genes
plotMA(deLorgeril_Resistant_dds_res_72_LFC , ylim = c(-5, 5)) # 72 hrs has less
plotMA(deLorgeril_Susceptible_dds_res_6_LFC, ylim = c(-5, 5))
plotMA(deLorgeril_Susceptible_dds_res_12_LFC, ylim = c(-5, 5))
plotMA(deLorgeril_Susceptible_dds_res_24_LFC, ylim = c(-5, 5))
plotMA(deLorgeril_Susceptible_dds_res_48_LFC, ylim = c(-5, 5))
plotMA(deLorgeril_Susceptible_dds_res_60_LFC, ylim = c(-5, 5)) # 60hr has many more significant genes
plotMA(deLorgeril_Susceptible_dds_res_72_LFC, ylim = c(-5, 5)) # 72hr has less

## Histogram of P values 
# exclude genes with very small counts to avoid spikes and plot using the LFCshrinkage
hist(deLorgeril_Resistant_dds_res_6_LFC $padj[deLorgeril_Resistant_dds_res_6_LFC$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white")
hist(deLorgeril_Resistant_dds_res_12_LFC $padj[deLorgeril_Resistant_dds_res_12_LFC$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white")
hist(deLorgeril_Resistant_dds_res_24_LFC $padj[deLorgeril_Resistant_dds_res_24_LFC$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white")
hist(deLorgeril_Resistant_dds_res_48_LFC $padj[deLorgeril_Resistant_dds_res_48_LFC$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white")
hist(deLorgeril_Resistant_dds_res_60_LFC $padj[deLorgeril_Resistant_dds_res_60_LFC$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white")
hist(deLorgeril_Resistant_dds_res_72_LFC $padj[deLorgeril_Resistant_dds_res_72_LFC$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white")

hist(deLorgeril_Susceptible_dds_res_6_LFC$padj[deLorgeril_Susceptible_dds_res_6_LFC$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white")
hist(deLorgeril_Susceptible_dds_res_12_LFC$padj[deLorgeril_Susceptible_dds_res_12_LFC$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white")
hist(deLorgeril_Susceptible_dds_res_24_LFC$padj[deLorgeril_Susceptible_dds_res_24_LFC$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white")
hist(deLorgeril_Susceptible_dds_res_48_LFC$padj[deLorgeril_Susceptible_dds_res_48_LFC$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white")
hist(deLorgeril_Susceptible_dds_res_60_LFC$padj[deLorgeril_Susceptible_dds_res_60_LFC$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white")
hist(deLorgeril_Susceptible_dds_res_72_LFC$padj[deLorgeril_Susceptible_dds_res_72_LFC$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white")
     
### Subsetting Significant Genes by padj < 0.05
# again, only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
deLorgeril_Resistant_dds_res_6_LFC_sig <- subset(  deLorgeril_Resistant_dds_res_6_LFC, padj < 0.05)
deLorgeril_Resistant_dds_res_12_LFC_sig <- subset(  deLorgeril_Resistant_dds_res_12_LFC, padj < 0.05)
deLorgeril_Resistant_dds_res_24_LFC_sig <- subset(  deLorgeril_Resistant_dds_res_24_LFC, padj < 0.05)
deLorgeril_Resistant_dds_res_48_LFC_sig <- subset(  deLorgeril_Resistant_dds_res_48_LFC, padj < 0.05)
deLorgeril_Resistant_dds_res_60_LFC_sig <- subset(  deLorgeril_Resistant_dds_res_60_LFC, padj < 0.05)
deLorgeril_Resistant_dds_res_72_LFC_sig <- subset(  deLorgeril_Resistant_dds_res_72_LFC, padj < 0.05)

deLorgeril_Susceptible_dds_res_6_LFC_sig <- subset(deLorgeril_Susceptible_dds_res_6_LFC, padj < 0.05)
deLorgeril_Susceptible_dds_res_12_LFC_sig <- subset(deLorgeril_Susceptible_dds_res_12_LFC, padj < 0.05)
deLorgeril_Susceptible_dds_res_24_LFC_sig <- subset(deLorgeril_Susceptible_dds_res_24_LFC, padj < 0.05)
deLorgeril_Susceptible_dds_res_48_LFC_sig <- subset(deLorgeril_Susceptible_dds_res_48_LFC, padj < 0.05)
deLorgeril_Susceptible_dds_res_60_LFC_sig <- subset(deLorgeril_Susceptible_dds_res_60_LFC, padj < 0.05)
deLorgeril_Susceptible_dds_res_72_LFC_sig <- subset(deLorgeril_Susceptible_dds_res_72_LFC, padj < 0.05)
 
deLorgeril_Resistant_dds_res_6_LFC_sig$transcript_id <- row.names(  deLorgeril_Resistant_dds_res_6_LFC_sig) 
deLorgeril_Resistant_dds_res_12_LFC_sig$transcript_id <- row.names(  deLorgeril_Resistant_dds_res_12_LFC_sig) 
deLorgeril_Resistant_dds_res_24_LFC_sig$transcript_id <- row.names(  deLorgeril_Resistant_dds_res_24_LFC_sig) 
deLorgeril_Resistant_dds_res_48_LFC_sig$transcript_id <- row.names(  deLorgeril_Resistant_dds_res_48_LFC_sig) 
deLorgeril_Resistant_dds_res_60_LFC_sig$transcript_id <- row.names(  deLorgeril_Resistant_dds_res_60_LFC_sig) 
deLorgeril_Resistant_dds_res_72_LFC_sig$transcript_id <- row.names(  deLorgeril_Resistant_dds_res_72_LFC_sig) 

deLorgeril_Susceptible_dds_res_6_LFC_sig$transcript_id <- row.names(deLorgeril_Susceptible_dds_res_6_LFC_sig) 
deLorgeril_Susceptible_dds_res_12_LFC_sig$transcript_id <- row.names(deLorgeril_Susceptible_dds_res_12_LFC_sig) 
deLorgeril_Susceptible_dds_res_24_LFC_sig$transcript_id <- row.names(deLorgeril_Susceptible_dds_res_24_LFC_sig) 
deLorgeril_Susceptible_dds_res_48_LFC_sig$transcript_id <- row.names(deLorgeril_Susceptible_dds_res_48_LFC_sig) 
deLorgeril_Susceptible_dds_res_60_LFC_sig$transcript_id <- row.names(deLorgeril_Susceptible_dds_res_60_LFC_sig) 
deLorgeril_Susceptible_dds_res_72_LFC_sig$transcript_id <- row.names(deLorgeril_Susceptible_dds_res_72_LFC_sig) 

deLorgeril_Resistant_dds_res_6_LFC_sig   <-as.data.frame(deLorgeril_Resistant_dds_res_6_LFC_sig)
deLorgeril_Resistant_dds_res_12_LFC_sig   <-as.data.frame(deLorgeril_Resistant_dds_res_12_LFC_sig)
deLorgeril_Resistant_dds_res_24_LFC_sig   <-as.data.frame(deLorgeril_Resistant_dds_res_24_LFC_sig)
deLorgeril_Resistant_dds_res_48_LFC_sig   <-as.data.frame(deLorgeril_Resistant_dds_res_48_LFC_sig)
deLorgeril_Resistant_dds_res_60_LFC_sig   <-as.data.frame(deLorgeril_Resistant_dds_res_60_LFC_sig)
deLorgeril_Resistant_dds_res_72_LFC_sig   <-as.data.frame(deLorgeril_Resistant_dds_res_72_LFC_sig)

deLorgeril_Susceptible_dds_res_6_LFC_sig<-as.data.frame(deLorgeril_Susceptible_dds_res_6_LFC_sig)
deLorgeril_Susceptible_dds_res_12_LFC_sig<-as.data.frame(deLorgeril_Susceptible_dds_res_12_LFC_sig)
deLorgeril_Susceptible_dds_res_24_LFC_sig<-as.data.frame(deLorgeril_Susceptible_dds_res_24_LFC_sig)
deLorgeril_Susceptible_dds_res_48_LFC_sig<-as.data.frame(deLorgeril_Susceptible_dds_res_48_LFC_sig)
deLorgeril_Susceptible_dds_res_60_LFC_sig<-as.data.frame(deLorgeril_Susceptible_dds_res_60_LFC_sig)
deLorgeril_Susceptible_dds_res_72_LFC_sig<-as.data.frame(deLorgeril_Susceptible_dds_res_72_LFC_sig)

nrow(deLorgeril_Resistant_dds_res_6_LFC_sig) #1644
nrow(deLorgeril_Resistant_dds_res_12_LFC_sig) #1644
nrow(deLorgeril_Resistant_dds_res_24_LFC_sig) #3775
nrow(deLorgeril_Resistant_dds_res_48_LFC_sig) #1593
nrow(deLorgeril_Resistant_dds_res_60_LFC_sig) # 3403
nrow(deLorgeril_Resistant_dds_res_72_LFC_sig) # 2309
nrow(deLorgeril_Susceptible_dds_res_6_LFC_sig) # 1445
nrow(deLorgeril_Susceptible_dds_res_12_LFC_sig) # 3435
nrow(deLorgeril_Susceptible_dds_res_24_LFC_sig) # 8298
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
deLorgeril_Resistant_dds_res_6_LFC_sig_APOP <-  merge(deLorgeril_Resistant_dds_res_6_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id" )
deLorgeril_Resistant_dds_res_12_LFC_sig_APOP <- merge(deLorgeril_Resistant_dds_res_12_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id" )
deLorgeril_Resistant_dds_res_24_LFC_sig_APOP <- merge(deLorgeril_Resistant_dds_res_24_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id" )
deLorgeril_Resistant_dds_res_48_LFC_sig_APOP <- merge(deLorgeril_Resistant_dds_res_48_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id" )
deLorgeril_Resistant_dds_res_60_LFC_sig_APOP <- merge(deLorgeril_Resistant_dds_res_60_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id" )
deLorgeril_Resistant_dds_res_72_LFC_sig_APOP <- merge(deLorgeril_Resistant_dds_res_72_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id" )

deLorgeril_Susceptible_dds_res_6_LFC_sig_APOP <-  merge(deLorgeril_Susceptible_dds_res_6_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")
deLorgeril_Susceptible_dds_res_12_LFC_sig_APOP <- merge(deLorgeril_Susceptible_dds_res_12_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")
deLorgeril_Susceptible_dds_res_24_LFC_sig_APOP <- merge(deLorgeril_Susceptible_dds_res_24_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")
deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP <- merge(deLorgeril_Susceptible_dds_res_48_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")
deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP <- merge(deLorgeril_Susceptible_dds_res_60_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")
deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP <- merge(deLorgeril_Susceptible_dds_res_72_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")

deLorgeril_Resistant_dds_res_6_LFC_sig_APOP  <- arrange(deLorgeril_Resistant_dds_res_6_LFC_sig_APOP , -log2FoldChange)
deLorgeril_Resistant_dds_res_12_LFC_sig_APOP  <- arrange(deLorgeril_Resistant_dds_res_12_LFC_sig_APOP , -log2FoldChange)
deLorgeril_Resistant_dds_res_24_LFC_sig_APOP  <- arrange(deLorgeril_Resistant_dds_res_24_LFC_sig_APOP , -log2FoldChange)
deLorgeril_Resistant_dds_res_48_LFC_sig_APOP  <- arrange(deLorgeril_Resistant_dds_res_48_LFC_sig_APOP , -log2FoldChange)
deLorgeril_Resistant_dds_res_60_LFC_sig_APOP  <- arrange(deLorgeril_Resistant_dds_res_60_LFC_sig_APOP , -log2FoldChange)
deLorgeril_Resistant_dds_res_72_LFC_sig_APOP  <- arrange(deLorgeril_Resistant_dds_res_72_LFC_sig_APOP , -log2FoldChange)

deLorgeril_Susceptible_dds_res_6_LFC_sig_APOP <- arrange(deLorgeril_Susceptible_dds_res_6_LFC_sig_APOP, -log2FoldChange)
deLorgeril_Susceptible_dds_res_12_LFC_sig_APOP <- arrange(deLorgeril_Susceptible_dds_res_12_LFC_sig_APOP, -log2FoldChange)
deLorgeril_Susceptible_dds_res_24_LFC_sig_APOP <- arrange(deLorgeril_Susceptible_dds_res_24_LFC_sig_APOP, -log2FoldChange)
deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP <- arrange(deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP, -log2FoldChange)
deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP <- arrange(deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP, -log2FoldChange)
deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP <- arrange(deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP, -log2FoldChange)

nrow(deLorgeril_Resistant_dds_res_6_LFC_sig_APOP) #21
nrow(deLorgeril_Resistant_dds_res_12_LFC_sig_APOP) #50
nrow(deLorgeril_Resistant_dds_res_24_LFC_sig_APOP) #73
nrow(deLorgeril_Resistant_dds_res_48_LFC_sig_APOP) #25
nrow(deLorgeril_Resistant_dds_res_60_LFC_sig_APOP) # 54
nrow(deLorgeril_Resistant_dds_res_72_LFC_sig_APOP) # 33

nrow(deLorgeril_Susceptible_dds_res_6_LFC_sig_APOP) #28
nrow(deLorgeril_Susceptible_dds_res_12_LFC_sig_APOP) #85
nrow(deLorgeril_Susceptible_dds_res_24_LFC_sig_APOP) #154
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
deLorgeril_Resistant_dds_res_6_LFC_sig_APOP  $group_by_sim <- "Resistant_dds_res_6"
deLorgeril_Resistant_dds_res_12_LFC_sig_APOP  $group_by_sim <- "Resistant_dds_res_12"
deLorgeril_Resistant_dds_res_24_LFC_sig_APOP  $group_by_sim <- "Resistant_dds_res_24"
deLorgeril_Resistant_dds_res_48_LFC_sig_APOP  $group_by_sim <- "Resistant_dds_res_48"
deLorgeril_Resistant_dds_res_60_LFC_sig_APOP  $group_by_sim <- "Resistant_dds_res_60"
deLorgeril_Resistant_dds_res_72_LFC_sig_APOP  $group_by_sim <- "Resistant_dds_res_72"

deLorgeril_Susceptible_dds_res_6_LFC_sig_APOP$group_by_sim <- "Susceptible_dds_res_6"
deLorgeril_Susceptible_dds_res_12_LFC_sig_APOP$group_by_sim <- "Susceptible_dds_res_12"
deLorgeril_Susceptible_dds_res_24_LFC_sig_APOP$group_by_sim <- "Susceptible_dds_res_24"
deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP$group_by_sim <- "Susceptible_dds_res_48"
deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP$group_by_sim <- "Susceptible_dds_res_60"
deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP$group_by_sim <- "Susceptible_dds_res_72"

# combine data frames 
deLorgeril_all_sig_APOP <- rbind( 
  deLorgeril_Resistant_dds_res_6_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
  deLorgeril_Resistant_dds_res_12_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
  deLorgeril_Resistant_dds_res_24_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
deLorgeril_Resistant_dds_res_48_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
deLorgeril_Resistant_dds_res_60_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
deLorgeril_Resistant_dds_res_72_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
deLorgeril_Susceptible_dds_res_6_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
deLorgeril_Susceptible_dds_res_12_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
deLorgeril_Susceptible_dds_res_24_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
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

#### HE OSHV1 TRANSCRIPTOME ANALYSIS ####
## LOAD DATA
He_counts <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/He_transcript_count_matrix.csv", header=TRUE,
                         row.names = "transcript_id")
head(He_counts )
colnames(He_counts)

# remove MSTRG novel transcript lines (can assess these later if necessary)
He_counts <- He_counts[!grepl("MSTRG", row.names(He_counts)),]

# Cut the "rna-" from the beginning of rownames
remove_rna = function(x){
  return(gsub("rna-","",x))
}
row.names(He_counts ) <- remove_rna(row.names(He_counts ))
head(He_counts )

#Load in sample metadata
He_coldata <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/He_coldata.csv", row.names = 1 )
View(He_coldata)  
nrow(He_coldata) 

# Make sure the columns of the count matrix and rows of the column data (sample metadata) are in the same order. 
all(rownames(He_coldata) %in% colnames(He_counts ))  #Should return TRUE
# returns TRUE
all(colnames(He_counts ) %in% rownames(He_coldata))  
# returns TRUE
all(rownames(He_coldata) == colnames(He_counts ))  # should return TRUE
# returns true

### DATA QC PCA PLOT 
# rlog transform data is recommended over vst for small data sets 
# PCA plots of data (https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/02_Preprocessing_Data.nb.html#count-distribution-boxplots)
He_counts_matrix <- as.matrix(He_counts)
Hevstcounts <- vst(He_counts_matrix, blind =TRUE)

# run PCA
pcHe <- prcomp(t(Hevstcounts))

# Plot PCA
autoplot(pcHe,
         data = He_coldata, 
         colour="Condition", 
         size=5) # clustering by condition explains more of the varaince than clustering by time, however very little of the variance can
        # be explained in the PCA (~12%)
autoplot(pcHe,
         data = He_coldata, 
         colour="Time", 
         size=5) # replicates cluster by time within treatment, not between control and treated
# Plot PCA 2 and 3 for comparison
autoplot(pcHe,
         data = He_coldata, 
         colour = "Condition", 
         size = 5,
         x = 2,
         y = 3) 
# Plot PCA 4 and 5 for comparison
autoplot(pcHe,
         data = He_coldata, 
         colour = "Condition", 
         size = 5,
         x = 3,
         y = 4) 

## MAKE DESEQ DATA SET FROM MATRIX
# This object specifies the count data and metadata you will work with. The design piece is critical.
# Correct for batch effects if necessary in this original formula: see this thread https://support.bioconductor.org/p/121408/
# do not correct counts using the removeBatchEffects from limma based on thread above 

## Check levels 
# It is prefered in R that the first level of a factor be the reference level for comparison
# (e.g. control, or untreated samples), so we can relevel the factor like so
# Check factor levels, set it so that comparison group is the first
levels(He_coldata$Condition) # "control" "OsHV1"  
levels(He_coldata$Time) # "120hr" "12h"   "24h"   "48h"   "6h"    "Time0"
He_coldata$Time <- factor(He_coldata$Time, levels=c("Time0","6h", "12h", "24h" , "48h" ,"120hr" ))

## Creating three here so I can compare the results
He_dds <- DESeqDataSetFromMatrix(countData = He_counts,
                                    colData = He_coldata,
                                    design = ~ Condition + Time)  # again accounting for time in my formula, keeping condition as the last

## Prefiltering the data
# Data prefiltering helps decrease the size of the data set and get rid of
# rows with no data or very minimal data (<10). Apply a minimal filtering here as more stringent filtering will be applied later
He_dds <- He_dds [ rowSums(counts(He_dds )) > 10, ]

## DATA TRANSFORMATION AND VISUALIZATION
# Assess sample clustering after setting initial formula for comparison
He_dds_vst <- vst(He_dds, blind = TRUE) # keep blind = true before deseq function has been run

## PCA plot visualization of individuals in the family, use vst because greater than 30 samples 
plotPCA(He_dds_vst, intgroup=c("Time", "Condition")) # highly variable, a bit more of variation explained than before

### DIFFERENTIAL EXPRESSION ANALYSIS
# run pipeline with single command because the formula has already been specified
# Steps: estimation of size factors (controlling for differences in the sequencing depth of the samples), 
# the estimation of dispersion values for each gene,and fitting a generalized linear model.
He_dds_deseq <- DESeq(He_dds) 

## Check the resultsNames object of each to look at the available coefficients for use in lfcShrink command
resultsNames(He_dds_deseq) # [1] "Intercept", "Time_6h_vs_Time0","Time_12h_vs_Time0","Time_24h_vs_Time0","Time_48h_vs_Time0","Time_120hr_vs_Time0" ,
# "Condition_OsHV1_vs_control"

## BUILD THE RESULTS OBJECT
# Examining the results object, change alpha to p <0.05, looking at object metadata
# use mcols to look at metadata for each table
mcols(He_dds_deseq)
He_dds_res <- results(He_dds_deseq, alpha=0.05, name= "Condition_OsHV1_vs_control")
He_dds_res_6hr <- results(He_dds_deseq, alpha=0.05, name= "Time_6h_vs_Time0" )
He_dds_res_12hr <- results(He_dds_deseq, alpha=0.05, name= "Time_12h_vs_Time0")
He_dds_res_24hr <- results(He_dds_deseq, alpha=0.05, name=  "Time_24h_vs_Time0")
He_dds_res_48hr <- results(He_dds_deseq, alpha=0.05, name= "Time_48h_vs_Time0")
He_dds_res_120hr <- results(He_dds_deseq, alpha=0.05, name= "Time_120hr_vs_Time0")
He_dds_res <- results(He_dds_deseq, alpha=0.05, name= "Condition_OsHV1_vs_control")
head(He_dds_res) 

summary(He_dds_res_6hr )
summary(He_dds_res_12hr )
summary(He_dds_res_24hr )
summary(He_dds_res_48hr )
summary(He_dds_res_120hr)


### Perform LFC Shrinkage with apeglm
## NOTES 
# Before plotting we need to apply an LFC shrinkage, set the coef as the specific comparison in the ResultsNames function of the deseq object
# Issue that the specific comparisons I set in my results formulas are not available in the ResultsNames coef list.
# More detailed notes about LFC Shrinkage are in the code for the Zhang Vibrio

## DECISION: USE SAME RES OBJECT TO KEEP ALPHA ADJUSTMENT, and use LFCShrink apeglm
He_dds_res_LFC<- lfcShrink(He_dds_deseq, coef="Condition_OsHV1_vs_control", type="apeglm", res= He_dds_res)
He_dds_res_6hr  <- lfcShrink(He_dds_deseq, coef="Time_6h_vs_Time0"   , type= "apeglm", res=He_dds_res_6hr)
He_dds_res_12hr <- lfcShrink(He_dds_deseq, coef="Time_12h_vs_Time0"  , type= "apeglm", res=He_dds_res_12hr )
He_dds_res_24hr <- lfcShrink(He_dds_deseq, coef= "Time_24h_vs_Time0" , type= "apeglm", res=He_dds_res_24hr )
He_dds_res_48hr <- lfcShrink(He_dds_deseq, coef="Time_48h_vs_Time0"  , type= "apeglm", res=He_dds_res_48hr )
He_dds_res_120hr <- lfcShrink(He_dds_deseq, coef="Time_120hr_vs_Time0", type= "apeglm", res=He_dds_res_120hr)

## EXPLORATORY PLOTTING OF RESULTS 
## MA Plotting
plotMA(He_dds_res_LFC, ylim = c(-5, 5))
plotMA(He_dds_res_6hr  , ylim=c(-5,5))
plotMA(He_dds_res_12hr , ylim=c(-5,5))
plotMA(He_dds_res_24hr , ylim=c(-5,5))
plotMA(He_dds_res_48hr , ylim=c(-5,5))
plotMA(He_dds_res_120hr, ylim=c(-5,5))

## Histogram of P values 
# exclude genes with very small counts to avoid spikes and plot using the LFCshrinkage
hist(He_dds_res_LFC$padj[He_dds_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

### Subsetting Significant Genes by padj < 0.05
# again, only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
He_dds_res_LFC_sig <-  subset(He_dds_res_LFC , padj < 0.05)
He_dds_res_LFC_sig $transcript_id <- row.names(He_dds_res_LFC_sig )
He_dds_res_LFC_sig  <- as.data.frame(He_dds_res_LFC_sig )
nrow(He_dds_res_LFC_sig ) # 3322
summary(He_dds_res_LFC_sig)

He_dds_res_6hr_sig <- subset(He_dds_res_6hr, padj < 0.05)
He_dds_res_12hr_sig <- subset(He_dds_res_12hr, padj < 0.05)
He_dds_res_24hr_sig <- subset(He_dds_res_24hr, padj < 0.05)
He_dds_res_48hr_sig <- subset(He_dds_res_48hr, padj < 0.05)
He_dds_res_120hr_sig  <- subset(He_dds_res_120hr, padj < 0.05)

He_dds_res_6hr_sig $transcript_id <- row.names(He_dds_res_6hr_sig )
He_dds_res_12hr_sig $transcript_id <- row.names(He_dds_res_12hr_sig )
He_dds_res_24hr_sig $transcript_id <- row.names(He_dds_res_24hr_sig )
He_dds_res_48hr_sig $transcript_id <- row.names(He_dds_res_48hr_sig )
He_dds_res_120hr_sig$transcript_id <- row.names(He_dds_res_120hr_sig)

He_dds_res_6hr_sig <- as.data.frame(He_dds_res_6hr_sig)
He_dds_res_12hr_sig <- as.data.frame(He_dds_res_12hr_sig)
He_dds_res_24hr_sig <- as.data.frame(He_dds_res_24hr_sig)
He_dds_res_48hr_sig <- as.data.frame(He_dds_res_48hr_sig)
He_dds_res_120hr_sig <- as.data.frame(He_dds_res_120hr_sig)

nrow(He_dds_res_6hr_sig) #742
nrow(He_dds_res_12hr_sig) # 525
nrow(He_dds_res_24hr_sig) # 699
nrow(He_dds_res_48hr_sig) # 641
nrow(He_dds_res_120hr_sig) # 520 

### GENE CLUSTERING ANALYSIS HEATMAPS  
# Extract genes with the highest variance across samples for each comparison using either vst or rlog transformed data
# This heatmap rather than plotting absolute expression strength plot the amount by which each gene deviates in a specific sample from the gene’s average across all samples. 
# example codes from RNAseq workflow: https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#other-comparisons
He_dds_res_LFC_sig_vst <-  head(order(rowVars(assay(He_dds_vst )), decreasing = TRUE), 200)
He_mat <- assay(He_dds_vst)[He_dds_res_LFC_sig_vst,]
He_mat <- He_mat - rowMeans(He_mat)
He_anno <- as.data.frame(colData(He_dds_vst)[, c("Condition","Time")])
He_heatmap <- pheatmap(He_mat, annotation_col = He_anno)
head(He_mat) # osHV1 and control cluster pretty strongly 

# reorder annotation table to match ordering in heatmap 
He_heatmap_reorder <-rownames(He_mat[He_heatmap$tree_row[["order"]],])
# annotate the row.names
He_mat_prot <- as.data.frame(He_heatmap_reorder )
colnames(He_mat_prot)[1] <- "transcript_id"
He_mat_prot_annot <- left_join(He_mat_prot, select(C_gig_rtracklayer_transcripts, transcript_id, product, gene), by = "transcript_id")

# Gene clustering heatmap with only apoptosis genes #
# vector C_gig_apop transcript IDs
C_gig_rtracklayer_apop_product_final_transcript_id 
# Search original Rubio_counts for apoptosis genes and do rlog on just these
He_counts_apop <- He_counts[row.names(He_counts) %in% C_gig_rtracklayer_apop_product_final_transcript_id,]
nrow(He_counts_apop ) #659
He_counts_apop_dds <- DESeqDataSetFromMatrix(countData = He_counts_apop,
                                                colData = He_coldata,
                                                design = ~ Time + Condition) # using same formula as before
# Prefiltering the data and running vst
He_counts_apop_dds <- He_counts_apop_dds[ rowSums(counts(He_counts_apop_dds)) > 10, ]
He_counts_apop_dds_vst <- varianceStabilizingTransformation(He_counts_apop_dds, blind=TRUE)

## PCA plot of rlog transformed counts for apoptosis
plotPCA(He_counts_apop_dds_vst , intgroup="Condition") # good separation of control and challenge

# heatmap of all apoptosis genes 
He_counts_apop_assay <-  assay(He_counts_apop_dds_vst)[,]
He_counts_apop_assay_mat <- He_counts_apop_assay - rowMeans(He_counts_apop_assay)
He_counts_apop_assay_anno <- as.data.frame(colData(He_counts_apop_dds_vst)[, c("Condition","Time")])
He_counts_apop_assay_heatmap <- pheatmap(He_counts_apop_assay_mat  , annotation_col = He_counts_apop_assay_anno)
head(He_counts_apop_assay_mat ) 
# good clustering still by control and OsHV1 

# heatmap of most variable apoptosis genes (this selects genes with the greatest variance in the sample)
topVarGenes_He_counts_apop_assay <-  head(order(rowVars(assay(He_counts_apop_dds_vst )), decreasing = TRUE), 100) 
top_Var_He_counts_apop_assay_mat<- assay(He_counts_apop_dds_vst)[topVarGenes_He_counts_apop_assay,]
top_Var_He_counts_apop_assay_mat <- top_Var_He_counts_apop_assay_mat - rowMeans(top_Var_He_counts_apop_assay_mat)
top_Var_He_counts_apop_assay_anno <- as.data.frame(colData(He_counts_apop_dds_vst)[, c("Condition","Time")])
top_Var_He_counts_apop_assay_heatmap <- pheatmap(top_Var_He_counts_apop_assay_mat  , annotation_col = top_Var_He_counts_apop_assay_anno)
head(top_Var_He_counts_apop_assay_mat )
# clustering not as good as looking at heatmap of all genes 

# annotate the top 100 most variable genes  
# reorder annotation table to match ordering in heatmap 
top_Var_He_counts_apop_assay_heatmap_reorder <-rownames(top_Var_He_counts_apop_assay_mat[top_Var_He_counts_apop_assay_heatmap$tree_row[["order"]],])
# annotate the row.names
top_Var_He_counts_apop_assay_prot <- as.data.frame(top_Var_He_counts_apop_assay_heatmap_reorder)
colnames(top_Var_He_counts_apop_assay_prot)[1] <- "transcript_id"
top_Var_He_counts_apop_assay_prot_annot <- left_join(top_Var_He_counts_apop_assay_prot, select(C_gig_rtracklayer_transcripts, transcript_id, product, gene), by = "transcript_id")

### Extract list of significant Apoptosis Genes (not less than or greater than 1 LFC) using merge
He_dds_res_LFC_sig_APOP <- merge(He_dds_res_LFC_sig , C_gig_rtracklayer_apop_product_final, by = "transcript_id")
He_dds_res_LFC_sig_arranged <-arrange(He_dds_res_LFC_sig_APOP, -log2FoldChange)
nrow(He_dds_res_LFC_sig_arranged ) #66
View(He_dds_res_LFC_sig_arranged)

He_dds_res_6hr_sig_APOP <- merge(He_dds_res_6hr_sig,    C_gig_rtracklayer_apop_product_final, by = "transcript_id")
He_dds_res_12hr_sig_APOP <- merge(He_dds_res_12hr_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")
He_dds_res_24hr_sig_APOP <- merge(He_dds_res_24hr_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")
He_dds_res_48hr_sig_APOP <- merge(He_dds_res_48hr_sig , C_gig_rtracklayer_apop_product_final, by = "transcript_id")
He_dds_res_120hr_sig_APOP <- merge(He_dds_res_120hr_sig ,C_gig_rtracklayer_apop_product_final, by = "transcript_id" )

He_dds_res_6hr_sig_arranged <- arrange(He_dds_res_6hr_sig_APOP , -log2FoldChange)
He_dds_res_12hr_sig_arranged <- arrange(He_dds_res_12hr_sig_APOP , -log2FoldChange)
He_dds_res_24hr_sig_arranged <- arrange(He_dds_res_24hr_sig_APOP , -log2FoldChange)
He_dds_res_48hr_sig_arranged <- arrange(He_dds_res_48hr_sig_APOP , -log2FoldChange)
He_dds_res_120hr_sig_arranged <- arrange(He_dds_res_120hr_sig_APOP, -log2FoldChange)

nrow(He_dds_res_6hr_sig_arranged ) # 20
nrow(He_dds_res_12hr_sig_arranged) # 14
nrow(He_dds_res_24hr_sig_arranged) # 19
nrow(He_dds_res_48hr_sig_arranged) # 15
nrow(He_dds_res_120hr_sig_arranged) # 9

He_dds_res_6hr_sig_APOP$group_by_sim <-"He_dds_res_6hr_sig_APOP"
He_dds_res_12hr_sig_APOP$group_by_sim <-"He_dds_res_12hr_sig_APOP"
He_dds_res_24hr_sig_APOP$group_by_sim <-"He_dds_res_24hr_sig_APOP"
He_dds_res_48hr_sig_APOP$group_by_sim <-"He_dds_res_48hr_sig_APOP"
He_dds_res_120hr_sig_APOP$group_by_sim <-"He_dds_res_120hr_sig_APOP"

# combine data frames 
He_all_sig_APOP <- rbind(He_dds_res_6hr_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
He_dds_res_12hr_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
He_dds_res_24hr_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
He_dds_res_48hr_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
He_dds_res_120hr_sig_APOP[,c("product","group_by_sim","log2FoldChange")])

# Make plot
He_full_LFC_plot <- ggplot(He_all_sig_APOP , aes(x=product,y=log2FoldChange, fill=group_by_sim )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()


# Compare apoptosis genes between group_by_sim groups
He_dds_res_LFC_sig_APOP_plot <- ggplot(He_dds_res_LFC_sig_APOP , aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("OsHV1 vs Control") +
  ylab("Log2 Fold Change")
# mostly upregulation

#### PROBIOTIC TRANSCRIPTOME ANALYSIS ####

## LOAD DATA
Probiotic_counts <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/Probiotic_transcript_count_matrix.csv", header=TRUE,
                         row.names = "transcript_id")
head(Probiotic_counts)
colnames(Probiotic_counts)

# colnames all have "_1", remove this. It was an artifact of the PE sample names
colnames(Probiotic_counts ) <- sub('\\_[^_]+$', '', colnames(Probiotic_counts))
colnames(Probiotic_counts )

# remove MSTRG novel transcript lines (can assess these later if necessary)
Probiotic_counts <- Probiotic_counts[!grepl("MSTRG", row.names(Probiotic_counts)),]
row.names(Probiotic_counts) <- remove_rna(row.names(Probiotic_counts))
head(Probiotic_counts)

#Load in sample metadata
Probiotic_coldata <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/Probiotic_coldata.csv", row.names = 1 )
View(Probiotic_coldata)  
nrow(Probiotic_coldata) 

# Make sure the columns of the count matrix and rows of the column data (sample metadata) are in the same order. 
all(rownames(Probiotic_coldata) %in% colnames(Probiotic_counts))  #Should return TRUE
# returns TRUE
all(colnames(Probiotic_counts) %in% rownames(Probiotic_coldata))  
# returns TRUE
all(rownames(Probiotic_coldata) == colnames(Probiotic_counts))
# returns FALSE

# Fix the order
Probiotic_counts <-Probiotic_counts[,row.names(Probiotic_coldata)]
row.names(Probiotic_coldata)

all(rownames(Probiotic_coldata) %in% colnames(Probiotic_counts))  #Should return TRUE
# returns TRUE
all(colnames(Probiotic_counts) %in% rownames(Probiotic_coldata))  
# returns TRUE
all(rownames(Probiotic_coldata) == colnames(Probiotic_counts))
# returns TRUE

### DATA QC PCA PLOT 
# rlog transform data is recommended over vst for small data sets 
# PCA plots of data (https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/02_Preprocessing_Data.nb.html#count-distribution-boxplots)
Probiotic_counts_matrix <- as.matrix(Probiotic_counts)
Probioticrlogcounts <- rlog(Probiotic_counts_matrix, blind =TRUE)

# run PCA
pcProbiotic <- prcomp(t(Probioticrlogcounts))

# Plot PCA
autoplot(pcProbiotic,
         data = Probiotic_coldata, 
         colour="Time", 
         size=5) # strong clustering by Day!!, ~50% of the variation explained by these first two PCA axes. 
autoplot(pcProbiotic,
         data = Probiotic_coldata, 
         colour="Condition", 
         size=5)
# Plot PCA 2 and 3 for comparison
autoplot(pcProbiotic,
         data = Probiotic_coldata, 
         colour = "Time", 
         size = 5,
         x = 2,
         y = 3) 
# Plot PCA 4 and 4 for comparison
autoplot(pcProbiotic,
         data = Probiotic_coldata, 
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
levels(Probiotic_coldata$Condition) # "Bacillus_pumilus_RI0695" "Untreated_control" 
Probiotic_coldata$Condition <- factor(Probiotic_coldata$Condition , levels = c("Untreated_control","Bacillus_pumilus_RI0695"))
levels(Probiotic_coldata$Condition)
levels(Probiotic_coldata$Time) # "12_d" "16_d" "5_d" 
Probiotic_coldata$Time <- factor(Probiotic_coldata$Time , levels = c("5_d","12_d", "16_d"))
levels(Probiotic_coldata$Time)

## Creating deseq data set from matrix, controlling for the effect of time
Probiotic_dds <- DESeqDataSetFromMatrix(countData = Probiotic_counts,
                                    colData = Probiotic_coldata,
                                    design = ~ Time + Condition ) 

## Prefiltering the data
# Data prefiltering helps decrease the size of the data set and get rid of
# rows with no data or very minimal data (<10). Apply a minimal filtering here as more stringent filtering will be applied later
Probiotic_dds <- Probiotic_dds [ rowSums(counts(Probiotic_dds )) > 10, ]

## DATA TRANSFORMATION AND VISUALIZATION
# Assess sample clustering after setting initial formula for comparison
Probiotic_dds_rlog <- rlog(Probiotic_dds , blind = TRUE) # keep blind = true before deseq function has been run

## PCA plot visualization of individuals in the family 
plotPCA(Probiotic_dds_rlog, intgroup=c("Sample", "Condition")) # clustering by time is not as tight

### DIFFERENTIAL EXPRESSION ANALYSIS
# run pipeline with single command because the formula has already been specified
# Steps: estimation of size factors (controlling for differences in the sequencing depth of the samples), 
# the estimation of dispersion values for each gene,and fitting a generalized linear model.
Probiotic_dds_deseq <- DESeq(Probiotic_dds) 
#-- note: fitType='parametric', but the dispersion trend was not well captured by the
#function: y = a/x + b, and a local regression fit was automatically substituted.
#specify fitType='local' or 'mean' to avoid this message next time.

## Check the resultsNames object of each to look at the available coefficients for use in lfcShrink command
resultsNames(Probiotic_dds_deseq) # [1] "Intercept", "Time_12_d_vs_5_d", "Time_16_d_vs_5_d", "Condition_Bacillus_pumilus_RI0695_vs_Untreated_control"

## BUILD THE RESULTS OBJECT
# Examining the results object, change alpha to p <0.05, looking at object metadata
# use mcols to look at metadata for each table
mcols(Probiotic_dds_deseq)
Probiotic_dds_deseq_Challenge_res <- results(Probiotic_dds_deseq, alpha=0.05, name= "Condition_Bacillus_pumilus_RI0695_vs_Untreated_control")
head(Probiotic_dds_deseq_Challenge_res) # Condition Bacillus pumilus RI0695 vs Untreated control 

### Perform LFC Shrinkage with apeglm
## NOTES 
# Before plotting we need to apply an LFC shrinkage, set the coef as the specific comparison in the ResultsNames function of the deseq object
# Issue that the specific comparisons I set in my results formulas are not available in the ResultsNames coef list.
# More detailed notes about LFC Shrinkage are in the code for the Zhang Vibrio

## DECISION: USE SAME RES OBJECT TO KEEP ALPHA ADJUSTMENT, and use LFCShrink apeglm
Probiotic_dds_deseq_Challenge_res_LFC<- lfcShrink(Probiotic_dds_deseq, coef="Condition_Bacillus_pumilus_RI0695_vs_Untreated_control", type="apeglm", res= Probiotic_dds_deseq_Challenge_res)

## EXPLORATORY PLOTTING OF RESULTS 
## MA Plotting
plotMA(Probiotic_dds_deseq_Challenge_res_LFC, ylim = c(-5, 5))

## Histogram of P values 
# exclude genes with very small counts to avoid spikes and plot using the LFCshrinkage
hist(Probiotic_dds_deseq_Challenge_res_LFC$padj[Probiotic_dds_deseq_Challenge_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

### Subsetting Significant Genes by padj < 0.05
# again, only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
Probiotic_dds_deseq_Challenge_res_LFC_sig <-  subset(Probiotic_dds_deseq_Challenge_res_LFC , padj < 0.05)
Probiotic_dds_deseq_Challenge_res_LFC_sig$ID<- row.names(Probiotic_dds_deseq_Challenge_res_LFC_sig)
Probiotic_dds_deseq_Challenge_res_LFC_sig <- as.data.frame(Probiotic_dds_deseq_Challenge_res_LFC_sig)
nrow(Probiotic_dds_deseq_Challenge_res_LFC_sig) # 1762

### GENE CLUSTERING ANALYSIS HEATMAPS  
# Extract genes with the highest variance across samples for each comparison using either vst or rlog transformed data
# This heatmap rather than plotting absolute expression strength plot the amount by which each gene deviates in a specific sample from the gene’s average across all samples. 
# example codes from RNAseq workflow: https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#other-comparisons

Probiotic_dds_deseq_Challenge_res_LFC_sig_assay <-  head(order(rowVars(assay(Probiotic_dds_rlog  )), decreasing = TRUE), 200)
family_Probiotic_broken_mat <- assay(Probiotic_dds_rlog )[Probiotic_dds_deseq_Challenge_res_LFC_sig_assay ,]
family_Probiotic_broken_mat <- family_Probiotic_broken_mat - rowMeans(family_Probiotic_broken_mat)
family_Probiotic_broken_anno <- as.data.frame(colData(Probiotic_dds_rlog )[, c("Condition","Time")])
family_Probiotic_broken_heatmap <- pheatmap(family_Probiotic_broken_mat , annotation_col = family_Probiotic_broken_anno)
head(familyProbiotic_broken_mat) # some clustering by bacillus, still overall clustering by day

# Gene clustering heatmap with only apoptosis genes #
# vector C_vir_apop transcript IDs
C_vir_rtracklayer_apop_product_final_ID <- C_vir_rtracklayer_apop_product_final$ID
# Search original Probiotic_counts for apoptosis genes and do rlog on just these
Probiotic_counts_apop <- Probiotic_counts[row.names(Probiotic_counts) %in% C_vir_rtracklayer_apop_product_final_ID,]
nrow(Probiotic_counts_apop) #1026
head(Probiotic_counts_apop)
Probiotic_counts_apop_dds <- DESeqDataSetFromMatrix(countData = Probiotic_counts_apop,
                                                colData = Probiotic_coldata,
                                                design = ~Time + Condition) # add time to control for injection and time effect
# Prefiltering the data and running rlog
Probiotic_counts_apop_dds<- Probiotic_counts_apop_dds[ rowSums(counts(Probiotic_counts_apop_dds)) > 10, ]
Probiotic_counts_apop_dds_rlog <- rlog(Probiotic_counts_apop_dds, blind=TRUE)

## PCA plot of rlog transformed counts for apoptosis
plotPCA(Probiotic_counts_apop_dds_rlog , intgroup="Time") # still overall clustering by time and not condition

# heatmap of all apoptosis genes 
Probiotic_counts_apop_assay <-  assay(Probiotic_counts_apop_dds_rlog)[,]
Probiotic_counts_apop_assay_mat <- Probiotic_counts_apop_assay - rowMeans(Probiotic_counts_apop_assay)
Probiotic_counts_apop_assay_anno <- as.data.frame(colData(Probiotic_counts_apop_dds_rlog )[, c("Condition","Sample")])
Probiotic_counts_apop_assay_heatmap <- pheatmap(Probiotic_counts_apop_assay_mat  , annotation_col = Probiotic_counts_apop_assay_anno)
head(Probiotic_counts_apop_assay_mat ) # the untreated controls are clustering more now

# heatmap of most variable apoptosis genes (this selects genes with the greatest variance in the sample)
topVarGenes_Probiotic_counts_apop_assay <-  head(order(rowVars(assay(Probiotic_counts_apop_dds_rlog)), decreasing = TRUE), 100) 
top_Var_Probiotic_counts_apop_assay_mat<- assay(Probiotic_counts_apop_dds_rlog)[topVarGenes_Probiotic_counts_apop_assay,]
top_Var_Probiotic_counts_apop_assay_mat <- top_Var_Probiotic_counts_apop_assay_mat - rowMeans(top_Var_Probiotic_counts_apop_assay_mat)
top_Var_Probiotic_counts_apop_assay_anno <- as.data.frame(colData(Probiotic_counts_apop_dds_rlog)[, c("Condition","Time")])
top_Var_Probiotic_counts_apop_assay_heatmap <- pheatmap(top_Var_Probiotic_counts_apop_assay_mat  , annotation_col = top_Var_Probiotic_counts_apop_assay_anno)
head(top_Var_Probiotic_counts_apop_assay_mat )
# some clustering patterns here in signature

# annotate the top 100 most variable genes  
# reorder annotation table to match ordering in heatmap 
top_Var_Probiotic_counts_apop_assay_heatmap_reorder <-rownames(top_Var_Probiotic_counts_apop_assay_mat[top_Var_Probiotic_counts_apop_assay_heatmap$tree_row[["order"]],])
# annotate the row.names
top_Var_Probiotic_counts_apop_assay_prot <- as.data.frame(top_Var_Probiotic_counts_apop_assay_heatmap_reorder)
colnames(top_Var_Probiotic_counts_apop_assay_prot)[1] <- "ID"
top_Var_Probiotic_counts_apop_assay_prot_annot <- left_join(top_Var_Probiotic_counts_apop_assay_prot, select(C_vir_rtracklayer_apop_product_final, ID, product, gene), by = "ID")

### Extract list of significant Apoptosis Genes (not less than or greater than 1 LFC) using merge
Probiotic_dds_deseq_Challenge_res_LFC_sig_APOP <- merge(Probiotic_dds_deseq_Challenge_res_LFC_sig, C_vir_rtracklayer_apop_product_final, by = "ID")
Probiotic_dds_deseq_Challenge_res_LFC_sig_APOP_arranged <- arrange(Probiotic_dds_deseq_Challenge_res_LFC_sig_APOP, -log2FoldChange) 
nrow(Probiotic_dds_deseq_Challenge_res_LFC_sig_APOP) # 22

Probiotic_dds_deseq_Challenge_res_LFC_sig_APOP_plot <- ggplot(Probiotic_dds_deseq_Challenge_res_LFC_sig_APOP , aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Bacillus pumilus vs Control") +
  ylab("Log2 Fold Change")
 # mostly downregulation of transcripts 

#### MODAK PROBIOTIC RE22 TRANSCRIPTOME ANALYSIS ####

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
Pro_RE22_coldata <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/Modak_Pro_RE22_coldata.csv", row.names = 1 )
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
autoplot(pcPro_RE22,
         data = Pro_RE22_coldata, 
         colour="Family", # Family (source of material_) explains the largest amount of variance
         size=5) 
autoplot(pcPro_RE22,
         data = Pro_RE22_coldata, 
         colour="Time", #same exact pattern as the Family 
         size=5) 

# Plot PCA 2 and 3 for comparison
autoplot(pcPro_RE22,
         data = Pro_RE22_coldata, 
         colour = "Condition", 
         size = 5,
         x = 2,
         y = 3)
autoplot(pcPro_RE22,
         data = Pro_RE22_coldata, 
         colour = "Family", 
         size = 5,
         x = 2,
         y = 3) 
# Plot PCA 4 and 4 for comparison
autoplot(pcPro_RE22,
         data = Pro_RE22_coldata, 
         colour = "Condition", 
         size = 5,
         x = 3,
         y = 4)  
# Overall little to no clustering by treatment 

## MAKE DESEQ DATA SET FROM MATRIX
# This object specifies the count data and metadata you will work with. The design piece is critical.
# Correct for batch effects if necessary in this original formula: see this thread https://support.bioconductor.org/p/121408/
# do not correct counts using the removeBatchEffects from limma based on thread above 

## Check levels 
# It is prefered in R that the first level of a factor be the reference level for comparison
# (e.g. control, or untreated samples), so we can relevel the factor like so
# Check factor levels, set it so that comparison group is the first
levels(Pro_RE22_coldata$Condition) # 
#"Bacillus_pumilus_RI06_95_exposure_24h"   "Bacillus_pumilus_RI06_95_exposure_6h"    "Control_no_treatment"                   
#"Phaeobacter_inhibens_S4_exposure_24h"    "Phaeobacter_inhibens_S4_exposure_6h"     "Vibrio_coralliilyticus_RE22_exposure_6h"

Pro_RE22_coldata$Condition <- factor(Pro_RE22_coldata$Condition , levels = c("Control_no_treatment"  , "Bacillus_pumilus_RI06_95_exposure_6h" ,"Bacillus_pumilus_RI06_95_exposure_24h"    ,                   
"Phaeobacter_inhibens_S4_exposure_6h" ,"Phaeobacter_inhibens_S4_exposure_24h" ,  "Vibrio_coralliilyticus_RE22_exposure_6h"
))

levels(Pro_RE22_coldata$Condition)
levels(Pro_RE22_coldata$Time) # Time this is what I will control for

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
plotPCA(Pro_RE22_dds_rlog, intgroup="Condition") # less clustering now by family, some clustering of different treatments

### DIFFERENTIAL EXPRESSION ANALYSIS
# run pipeline with single command because the formula has already been specified
# Steps: estimation of size factors (controlling for differences in the sequencing depth of the samples), 
# the estimation of dispersion values for each gene,and fitting a generalized linear model.
Pro_RE22_dds_deseq <- DESeq(Pro_RE22_dds) 

## Check the resultsNames object of each to look at the available coefficients for use in lfcShrink command
resultsNames(Pro_RE22_dds_deseq) #[1] "Intercept"                                                                
# [2] "Time_6d_vs_10d"                                                           
# [3] "Time_7d_vs_10d"                                                           
# [4] "Condition_Bacillus_pumilus_RI06_95_exposure_6h_vs_Control_no_treatment"   
# [5] "Condition_Bacillus_pumilus_RI06_95_exposure_24h_vs_Control_no_treatment"  
# [6] "Condition_Phaeobacter_inhibens_S4_exposure_6h_vs_Control_no_treatment"    
# [7] "Condition_Phaeobacter_inhibens_S4_exposure_24h_vs_Control_no_treatment"   
# [8] "Condition_Vibrio_coralliilyticus_RE22_exposure_6h_vs_Control_no_treatment"

## BUILD THE RESULTS OBJECT
# Examining the results object, change alpha to p <0.05, looking at object metadata
# use mcols to look at metadata for each table
mcols(Pro_RE22_dds_deseq)
Pro_RE22_dds_deseq_res_RI_6h <- results(Pro_RE22_dds_deseq, alpha=0.05, name= "Condition_Bacillus_pumilus_RI06_95_exposure_6h_vs_Control_no_treatment")
Pro_RE22_dds_deseq_res_RI_24h <- results(Pro_RE22_dds_deseq, alpha=0.05, name= "Condition_Bacillus_pumilus_RI06_95_exposure_24h_vs_Control_no_treatment" )
Pro_RE22_dds_deseq_res_S4_6h <- results(Pro_RE22_dds_deseq, alpha=0.05, name= "Condition_Phaeobacter_inhibens_S4_exposure_6h_vs_Control_no_treatment"  )
Pro_RE22_dds_deseq_res_S4_24h <- results(Pro_RE22_dds_deseq, alpha=0.05, name= "Condition_Phaeobacter_inhibens_S4_exposure_24h_vs_Control_no_treatment"   )
Pro_RE22_dds_deseq_res_RE22 <- results(Pro_RE22_dds_deseq, alpha=0.05, name= "Condition_Vibrio_coralliilyticus_RE22_exposure_6h_vs_Control_no_treatment"  )

### Perform LFC Shrinkage with apeglm
## NOTES 
# Before plotting we need to apply an LFC shrinkage, set the coef as the specific comparison in the ResultsNames function of the deseq object
# Issue that the specific comparisons I set in my results formulas are not available in the ResultsNames coef list.
# More detailed notes about LFC Shrinkage are in the code for the Zhang Vibrio

## DECISION: USE SAME RES OBJECT TO KEEP ALPHA ADJUSTMENT, and use LFCShrink apeglm
Pro_RE22_dds_deseq_res_RI_6h_LFC <- lfcShrink(Pro_RE22_dds_deseq, coef="Condition_Bacillus_pumilus_RI06_95_exposure_6h_vs_Control_no_treatment", type="apeglm", res= Pro_RE22_dds_deseq_res_RI_6h )
Pro_RE22_dds_deseq_res_RI_24h_LFC <- lfcShrink(Pro_RE22_dds_deseq, coef= "Condition_Bacillus_pumilus_RI06_95_exposure_24h_vs_Control_no_treatment" , type="apeglm", res= Pro_RE22_dds_deseq_res_RI_24h)
Pro_RE22_dds_deseq_res_S4_6h_LFC <- lfcShrink(Pro_RE22_dds_deseq, coef="Condition_Phaeobacter_inhibens_S4_exposure_6h_vs_Control_no_treatment", type="apeglm", res= Pro_RE22_dds_deseq_res_S4_6h )
Pro_RE22_dds_deseq_res_S4_24h_LFC <- lfcShrink(Pro_RE22_dds_deseq, coef= "Condition_Phaeobacter_inhibens_S4_exposure_24h_vs_Control_no_treatment"  , type="apeglm", res= Pro_RE22_dds_deseq_res_S4_24h)
Pro_RE22_dds_deseq_res_RE22_LFC <- lfcShrink(Pro_RE22_dds_deseq, coef= "Condition_Vibrio_coralliilyticus_RE22_exposure_6h_vs_Control_no_treatment" , type="apeglm", res= Pro_RE22_dds_deseq_res_RE22 )

## EXPLORATORY PLOTTING OF RESULTS 
## MA Plotting
plotMA(Pro_RE22_dds_deseq_res_RI_6h_LFC , ylim=c(-5,5))
plotMA(Pro_RE22_dds_deseq_res_RI_24h_LFC , ylim=c(-5,5))
plotMA(Pro_RE22_dds_deseq_res_S4_6h_LFC , ylim=c(-5,5))
plotMA(Pro_RE22_dds_deseq_res_S4_24h_LFC , ylim=c(-5,5))
plotMA(Pro_RE22_dds_deseq_res_RE22_LFC , ylim=c(-5,5))

## Histogram of P values 
# exclude genes with very small counts to avoid spikes and plot using the LFCshrinkage
hist(Pro_RE22_dds_deseq_res_RI_6h_LFC$padj[Pro_RE22_dds_deseq_res_RI_6h_LFC$baseMean >1], breaks=0:20/20, col="grey50", border="white")
hist(Pro_RE22_dds_deseq_res_RI_24h_LFC$padj[Pro_RE22_dds_deseq_res_RI_24h_LFC$baseMean >1], breaks=0:20/20, col="grey50", border="white")
hist(Pro_RE22_dds_deseq_res_S4_6h_LFC$padj[Pro_RE22_dds_deseq_res_S4_6h_LFC$baseMean >1], breaks=0:20/20, col="grey50", border="white")
hist(Pro_RE22_dds_deseq_res_S4_24h_LFC$padj[Pro_RE22_dds_deseq_res_S4_24h_LFC$baseMean >1], breaks=0:20/20, col="grey50", border="white")
hist(Pro_RE22_dds_deseq_res_RE22_LFC$padj[Pro_RE22_dds_deseq_res_RE22_LFC$baseMean >1], breaks=0:20/20, col="grey50", border="white")

### Subsetting Significant Genes by padj < 0.05
# again, only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering

Pro_RE22_dds_deseq_res_RI_6h_LFC_sig <- subset(Pro_RE22_dds_deseq_res_RI_6h_LFC , padj < 0.05)
Pro_RE22_dds_deseq_res_RI_24h_LFC_sig <- subset(Pro_RE22_dds_deseq_res_RI_24h_LFC, padj < 0.05)
Pro_RE22_dds_deseq_res_S4_6h_LFC_sig <-      subset(Pro_RE22_dds_deseq_res_S4_6h_LFC , padj < 0.05) 
Pro_RE22_dds_deseq_res_S4_24h_LFC_sig <- subset(Pro_RE22_dds_deseq_res_S4_24h_LFC, padj < 0.05)
Pro_RE22_dds_deseq_res_RE22_LFC_sig <-  subset(Pro_RE22_dds_deseq_res_RE22_LFC, padj < 0.05)

Pro_RE22_dds_deseq_res_RI_6h_LFC_sig$ID <- row.names(Pro_RE22_dds_deseq_res_RI_6h_LFC_sig)
Pro_RE22_dds_deseq_res_RI_24h_LFC_sig$ID <- row.names(Pro_RE22_dds_deseq_res_RI_24h_LFC_sig)
Pro_RE22_dds_deseq_res_S4_6h_LFC_sig$ID <- row.names(Pro_RE22_dds_deseq_res_S4_6h_LFC_sig)
Pro_RE22_dds_deseq_res_S4_24h_LFC_sig$ID <- row.names(Pro_RE22_dds_deseq_res_S4_24h_LFC_sig)
Pro_RE22_dds_deseq_res_RE22_LFC_sig$ID <- row.names(Pro_RE22_dds_deseq_res_RE22_LFC_sig)

Pro_RE22_dds_deseq_res_RI_6h_LFC_sig  <- as.data.frame(Pro_RE22_dds_deseq_res_RI_6h_LFC_sig)
Pro_RE22_dds_deseq_res_RI_24h_LFC_sig <- as.data.frame(Pro_RE22_dds_deseq_res_RI_24h_LFC_sig)
Pro_RE22_dds_deseq_res_S4_6h_LFC_sig  <- as.data.frame(Pro_RE22_dds_deseq_res_S4_6h_LFC_sig)
Pro_RE22_dds_deseq_res_S4_24h_LFC_sig <- as.data.frame(Pro_RE22_dds_deseq_res_S4_24h_LFC_sig)
Pro_RE22_dds_deseq_res_RE22_LFC_sig   <- as.data.frame(Pro_RE22_dds_deseq_res_RE22_LFC_sig)

nrow(Pro_RE22_dds_deseq_res_RI_6h_LFC_sig ) # 1795
nrow(Pro_RE22_dds_deseq_res_RI_24h_LFC_sig) # 2570
nrow(Pro_RE22_dds_deseq_res_S4_6h_LFC_sig ) # 2424
nrow(Pro_RE22_dds_deseq_res_S4_24h_LFC_sig) # 3683
nrow(Pro_RE22_dds_deseq_res_RE22_LFC_sig  ) # 2005

### GENE CLUSTERING ANALYSIS HEATMAPS  
# Extract genes with the highest variance across samples for each comparison using either vst or rlog transformed data
# This heatmap rather than plotting absolute expression strength plot the amount by which each gene deviates in a specific sample from the gene’s average across all samples. 
# example codes from RNAseq workflow: https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#other-comparisons

Pro_RE22_dds_deseq_Challenge_res_LFC_sig_assay <-  head(order(rowVars(assay(Pro_RE22_dds_rlog  )), decreasing = TRUE), 200)
family_Pro_RE22_broken_mat <- assay(Pro_RE22_dds_rlog )[Pro_RE22_dds_deseq_Challenge_res_LFC_sig_assay ,]
family_Pro_RE22_broken_mat <- family_Pro_RE22_broken_mat - rowMeans(family_Pro_RE22_broken_mat)
family_Pro_RE22_broken_anno <- as.data.frame(colData(Pro_RE22_dds_rlog )[, c("Condition","Family")])
family_Pro_RE22_broken_heatmap <- pheatmap(family_Pro_RE22_broken_mat , annotation_col = family_Pro_RE22_broken_anno)
head(family_Pro_RE22_broken_mat) # mostly clustering by treatment and still some by family 

# Gene clustering heatmap with only apoptosis genes #
# Search original Probiotic_counts for apoptosis genes and do rlog on just these
Pro_RE22_counts_apop <- Pro_RE22_counts[row.names(Pro_RE22_counts) %in% C_vir_rtracklayer_apop_product_final_ID,]
nrow(Pro_RE22_counts_apop) #1026
head(Pro_RE22_counts_apop)
Pro_RE22_counts_apop_dds <- DESeqDataSetFromMatrix(countData = Pro_RE22_counts_apop,
                                                    colData = Pro_RE22_coldata,
                                                    design = ~Time + Condition) # Control for larval age and source of larvae 
# Prefiltering the data and running rlog
Pro_RE22_counts_apop_dds<- Pro_RE22_counts_apop_dds[ rowSums(counts(Pro_RE22_counts_apop_dds)) > 10, ]
Pro_RE22_counts_apop_dds_rlog <- rlog(Pro_RE22_counts_apop_dds, blind=TRUE)

## PCA plot of rlog transformed counts for apoptosis
plotPCA(Pro_RE22_counts_apop_dds_rlog , intgroup="Family") # little clustering by treatment overall still clustering by family 

# heatmap of all apoptosis genes 
Pro_RE22_counts_apop_assay <-  assay(Pro_RE22_counts_apop_dds_rlog)[,]
Pro_RE22_counts_apop_assay_mat <- Pro_RE22_counts_apop_assay - rowMeans(Pro_RE22_counts_apop_assay)
Pro_RE22_counts_apop_assay_anno <- as.data.frame(colData(Pro_RE22_counts_apop_dds_rlog )[, c("Condition","Family")])
Pro_RE22_counts_apop_assay_heatmap <- pheatmap(Pro_RE22_counts_apop_assay_mat  , annotation_col = Pro_RE22_counts_apop_assay_anno)
head(Pro_RE22_counts_apop_assay_mat ) # more clustering by larvae source than by disease response 

# heatmap of most variable apoptosis genes (this selects genes with the greatest variance in the sample)
topVarGenes_Pro_RE22_counts_apop_assay <-  head(order(rowVars(assay(Pro_RE22_counts_apop_dds_rlog)), decreasing = TRUE), 100) 
top_Var_Pro_RE22_counts_apop_assay_mat<- assay(Pro_RE22_counts_apop_dds_rlog)[topVarGenes_Pro_RE22_counts_apop_assay,]
top_Var_Pro_RE22_counts_apop_assay_mat <- top_Var_Pro_RE22_counts_apop_assay_mat - rowMeans(top_Var_Pro_RE22_counts_apop_assay_mat)
top_Var_Pro_RE22_counts_apop_assay_anno <- as.data.frame(colData(Pro_RE22_counts_apop_dds_rlog)[, c("Condition","Family")])
top_Var_Pro_RE22_counts_apop_assay_heatmap <- pheatmap(top_Var_Pro_RE22_counts_apop_assay_mat  , annotation_col = top_Var_Pro_RE22_counts_apop_assay_anno)
head(top_Var_Pro_RE22_counts_apop_assay_mat ) # very high amoung of variation

# annotate the top 100 most variable genes  
# reorder annotation table to match ordering in heatmap 
top_Var_Pro_RE22_counts_apop_assay_heatmap_reorder <-rownames(top_Var_Pro_RE22_counts_apop_assay_mat[top_Var_Pro_RE22_counts_apop_assay_heatmap$tree_row[["order"]],])
# annotate the row.names
top_Var_Pro_RE22_counts_apop_assay_prot <- as.data.frame(top_Var_Pro_RE22_counts_apop_assay_heatmap_reorder)
colnames(top_Var_Pro_RE22_counts_apop_assay_prot)[1] <- "ID"
top_Var_Pro_RE22_counts_apop_assay_prot_annot <- left_join(top_Var_Pro_RE22_counts_apop_assay_prot, select(C_vir_rtracklayer_apop_product_final, ID, product, gene), by = "ID")

### Extract list of significant Apoptosis Genes (not less than or greater than 1 LFC) using merge
Pro_RE22_dds_deseq_res_RI_6h_LFC_sig_APOP <- merge(Pro_RE22_dds_deseq_res_RI_6h_LFC_sig, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_dds_deseq_res_RI_24h_LFC_sig_APOP <- merge(Pro_RE22_dds_deseq_res_RI_24h_LFC_sig, C_vir_rtracklayer_apop_product_final, by =  "ID")
Pro_RE22_dds_deseq_res_S4_6h_LFC_sig_APOP <- merge(Pro_RE22_dds_deseq_res_S4_6h_LFC_sig, C_vir_rtracklayer_apop_product_final, by = "ID")
Pro_RE22_dds_deseq_res_S4_24h_LFC_sig_APOP <- merge(Pro_RE22_dds_deseq_res_S4_24h_LFC_sig, C_vir_rtracklayer_apop_product_final, by =  "ID")
Pro_RE22_dds_deseq_res_RE22_LFC_sig_APOP <- merge(Pro_RE22_dds_deseq_res_RE22_LFC_sig, C_vir_rtracklayer_apop_product_final, by = "ID")

nrow(Pro_RE22_dds_deseq_res_RI_6h_LFC_sig_APOP ) # 17
nrow(Pro_RE22_dds_deseq_res_RI_24h_LFC_sig_APOP) # 32
nrow(Pro_RE22_dds_deseq_res_S4_6h_LFC_sig_APOP ) # 40
nrow(Pro_RE22_dds_deseq_res_S4_24h_LFC_sig_APOP) # 41
nrow(Pro_RE22_dds_deseq_res_RE22_LFC_sig_APOP ) # 27 

# Compare apoptosis genes between group_by_sim groups
Pro_RE22_dds_deseq_res_RI_6h_LFC_sig_APOP $group_by_sim <- "RI_6h"
Pro_RE22_dds_deseq_res_RI_24h_LFC_sig_APOP$group_by_sim <- "RI_24h"
Pro_RE22_dds_deseq_res_S4_6h_LFC_sig_APOP $group_by_sim <- "S4_6h"
Pro_RE22_dds_deseq_res_S4_24h_LFC_sig_APOP$group_by_sim <- "S4_24h"
Pro_RE22_dds_deseq_res_RE22_LFC_sig_APOP $group_by_sim <-  "RE22"

# combine data frames 
Pro_RE22_all_sig_APOP <- rbind( 
  Pro_RE22_dds_deseq_res_RI_6h_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
  Pro_RE22_dds_deseq_res_RI_24h_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
  Pro_RE22_dds_deseq_res_S4_6h_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
  Pro_RE22_dds_deseq_res_S4_24h_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
  Pro_RE22_dds_deseq_res_RE22_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")])

# Make plot or up and downregulated
Pro_RE22_all_sig_APOP_downregulated <- Pro_RE22_all_sig_APOP  %>% filter(log2FoldChange <= 0)
Pro_RE22_all_sig_APOP_upregulated <-   Pro_RE22_all_sig_APOP %>% filter(log2FoldChange > 0)

Pro_RE22_all_sig_APOP_downregulated_plot <- ggplot(Pro_RE22_all_sig_APOP_downregulated  , aes(x=product,y=log2FoldChange, fill=group_by_sim )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()
Pro_RE22_all_sig_APOP_upregulated_plot <- ggplot(Pro_RE22_all_sig_APOP_upregulated, aes(x=product,y=log2FoldChange, fill=group_by_sim )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()


Pro_RE22_dds_deseq_res_RI_6h_LFC_sig_APOP_plot <- ggplot(Pro_RE22_dds_deseq_res_RI_6h_LFC_sig_APOP , aes(x=product, y = log2FoldChange, fill= log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("RI 6h") 
Pro_RE22_dds_deseq_res_RI_24h_LFC_sig_APOP_plot <- ggplot(Pro_RE22_dds_deseq_res_RI_24h_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill= log2FoldChange))+ geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("RI 24h")
Pro_RE22_dds_deseq_res_S4_6h_LFC_sig_APOP_plot <- ggplot(Pro_RE22_dds_deseq_res_S4_6h_LFC_sig_APOP , aes(x=product, y = log2FoldChange, fill= log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("S4 6h")
Pro_RE22_dds_deseq_res_S4_24h_LFC_sig_APOP_plot <- ggplot(Pro_RE22_dds_deseq_res_S4_24h_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill= log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("S4 24h")
Pro_RE22_dds_deseq_res_RE22_LFC_sig_APOP_plot <- ggplot(Pro_RE22_dds_deseq_res_RE22_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill= log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("RE22 ")



#### ROD TRANSCRIPTOME ANALYSIS #### 

ROD_counts <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/ROD_transcript_count_matrix.csv", header=TRUE,
                             row.names = "transcript_id")
head(ROD_counts)
colnames(ROD_counts)

# remove MSTRG novel transcript lines (can assess these later if necessary)
ROD_counts <- ROD_counts[!grepl("MSTRG", row.names(ROD_counts)),]
row.names(ROD_counts) <- remove_rna(row.names(ROD_counts))
head(ROD_counts)

#Load in sample metadata
ROD_coldata <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/ROD_coldata.csv", row.names = 1 )
View(ROD_coldata )  
nrow(ROD_coldata ) 

# Make sure the columns of the count matrix and rows of the column data (sample metadata) are in the same order. 
all(rownames(ROD_coldata ) %in% colnames(ROD_counts))  #Should return TRUE
# returns TRUE
all(colnames(ROD_counts) %in% rownames(ROD_coldata ))  
# returns TRUE
all(rownames(ROD_coldata ) == colnames(ROD_counts))
# returns FALSE

# Fix the order
ROD_counts <- ROD_counts[,row.names(ROD_coldata)]

all(rownames(ROD_coldata ) %in% colnames(ROD_counts))  #Should return TRUE
# returns TRUE
all(colnames(ROD_counts) %in% rownames(ROD_coldata ))  
# returns TRUE
all(rownames(ROD_coldata ) == colnames(ROD_counts))
# returns TRUE

# going to split the analysis into the Susceptible 1d and 5d vs. Susceptible 5d 15, and then the control and resistant going to be controlled 

ROD_Susceptible_coldata <- ROD_coldata %>%  subset(Breed == "F3L") # subset() keeps rownames
row.names(ROD_Susceptible_coldata )
ROD_Susceptible_counts <- ROD_counts[,row.names(ROD_Susceptible_coldata )]
colnames(ROD_Susceptible_counts)
ROD_Resistant_coldata <- ROD_coldata %>%  subset(Breed == "GX09") # subset() keeps rownames
ROD_Resistant_counts <- ROD_counts[,row.names(ROD_Resistant_coldata)]

### DATA QC PCA PLOT 
# rlog transform data is recommended over vst for small data sets 
# PCA plots of data (https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/02_Preprocessing_Data.nb.html#count-distribution-boxplots)
ROD_Susceptible_counts_matrix <- as.matrix(ROD_Susceptible_counts)
ROD_Resistant_counts_matrix <- as.matrix(ROD_Resistant_counts)

ROD_Susceptible_rlogcounts <- rlog(ROD_Susceptible_counts_matrix , blind =TRUE)
ROD_Resistant_rlogcounts <- rlog(ROD_Resistant_counts_matrix , blind =TRUE)

# run PCA
pcROD_Susceptible <- prcomp(t(ROD_Susceptible_rlogcounts))
pcROD_Resistant <- prcomp(t(ROD_Resistant_rlogcounts))

# Plot PCA
autoplot(pcROD_Susceptible,
         data = ROD_Susceptible_coldata, 
         colour="Time", 
         size=5) # 1d and 15 days actually cluster most closely
autoplot(pcROD_Susceptible,
         data = ROD_Susceptible_coldata, 
         colour="Disease_stage", 
         size=5) # the sample with no ROD signs is a strong outlier, could use disease stage as the basis for comparison here, day 15 and 30 cluster closely together
autoplot(pcROD_Resistant,
         data = ROD_Resistant_coldata, 
         colour="Time", 
         size=5) # day 15 and 30 cluster most closely 
autoplot(pcROD_Resistant,
         data = ROD_Resistant_coldata, 
         colour="Condition", 
         size=5) # time explains most of the variation 
# Plot PCA 2 and 3 for comparison
autoplot(pcROD_Susceptible,
         data = ROD_Susceptible_coldata, 
         colour = "Time", 
         size = 5,
         x = 2,
         y = 3) # very spread
autoplot(pcROD_Resistant,
         data = ROD_Resistant_coldata, 
         colour = "Condition", 
         size = 5,
         x = 2,
         y = 3) # some clustering by day and condition

## MAKE DESEQ DATA SET FROM MATRIX
# This object specifies the count data and metadata you will work with. The design piece is critical.
# Correct for batch effects if necessary in this original formula: see this thread https://support.bioconductor.org/p/121408/
# do not correct counts using the removeBatchEffects from limma based on thread above 

## Check levels 
# It is prefered in R that the first level of a factor be the reference level for comparison
# (e.g. control, or untreated samples), so we can relevel the factor like so
# Check factor levels, set it so that comparison group is the first
levels(ROD_Susceptible_coldata$Condition) #"Control_Resistant"     "Resistant_Challenge"   "Susceptible_Challenge"
ROD_Susceptible_coldata$Condition <- droplevels(ROD_Susceptible_coldata$Condition)
levels(ROD_Susceptible_coldata$Condition)
levels(ROD_Susceptible_coldata$Disease_stage) # [1] "No _signs_of_ROD", "Signs_of_ROD_15%_cumulative_percent_mortality",  "Signs_of_ROD_30%_cumulative_percent_mortality"
#[4] "Signs_of_ROD_5%_cumulative_percent_mortality" 

levels(ROD_Resistant_coldata$Condition) # "Control_Resistant"   "Early_Susceptible"   "Late_Susecptible"    "Resistant_Challenge" 
ROD_Resistant_coldata$Condition <- droplevels(ROD_Resistant_coldata$Condition)
levels(ROD_Resistant_coldata$Condition)
levels(ROD_Resistant_coldata$Time ) # "15d" "1d"  "30d" "5d" 
ROD_Resistant_coldata$Time<- factor(ROD_Resistant_coldata$Time, levels = c( "1d","5d","15d",  "30d" ))
 
## Creating two data set from matrix, one for each family  
ROD_Resistant_dds <- DESeqDataSetFromMatrix(countData = ROD_Resistant_counts,
                                        colData = ROD_Resistant_coldata,
                                        design = ~Time + Condition) # control for time effect here to get at challenge response only 
ROD_Susceptible_dds <- DESeqDataSetFromMatrix(countData = ROD_Susceptible_counts,
                                            colData = ROD_Susceptible_coldata,
                                            design = ~ Condition) # not adding in time here because basically condition here is representing early vs late response 

## Prefiltering the data
# Data prefiltering helps decrease the size of the data set and get rid of
# rows with no data or very minimal data (<10). Apply a minimal filtering here as more stringent filtering will be applied later
ROD_Resistant_dds <- ROD_Resistant_dds[ rowSums(counts(ROD_Resistant_dds )) > 10, ]
ROD_Susceptible_dds <- ROD_Susceptible_dds[rowSums(counts(ROD_Susceptible_dds )) > 10, ]

## DATA TRANSFORMATION AND VISUALIZATION
# Assess sample clustering after setting initial formula for comparison
ROD_Resistant_dds_rlog <- rlog(ROD_Resistant_dds , blind = TRUE) # keep blind = true before deseq function has been run
ROD_Susceptible_dds_rlog <- rlog(ROD_Susceptible_dds , blind = TRUE) # keep blind = true before deseq function has been run

## PCA plot visualization of individuals in the family 
plotPCA(ROD_Resistant_dds_rlog, intgroup= "Time") # no clustering by treatment
plotPCA(ROD_Susceptible_dds_rlog , intgroup="Condition") # some clustering

### DIFFERENTIAL EXPRESSION ANALYSIS
# run pipeline with single command because the formula has already been specified
# Steps: estimation of size factors (controlling for differences in the sequencing depth of the samples), 
# the estimation of dispersion values for each gene,and fitting a generalized linear model.
ROD_Resistant_dds_deseq <- DESeq(ROD_Resistant_dds) 
ROD_Susceptible_dds_deseq <- DESeq(ROD_Susceptible_dds) 

## Check the resultsNames object of each to look at the available coefficients for use in lfcShrink command
resultsNames(ROD_Resistant_dds_deseq) # [1] "Intercept" , "Time_5d_vs_1d" , "Time_15d_vs_1d", "Time_30d_vs_1d", "Condition_Resistant_Challenge_vs_Control_Resistant"
resultsNames(ROD_Susceptible_dds_deseq) # [1] "Intercept"  "Condition_Late_Susecptible_vs_Early_Susceptible"

## BUILD THE RESULTS OBJECT
# Examining the results object, change alpha to p <0.05, looking at object metadata
# use mcols to look at metadata for each table
mcols(ROD_Resistant_dds_deseq)
mcols(ROD_Susceptible_dds_deseq)

ROD_Resistant_dds_res <- results(ROD_Resistant_dds_deseq, alpha=0.05, name= "Condition_Resistant_Challenge_vs_Control_Resistant")
head(ROD_Resistant_dds_res)  # : Condition Resistant Challenge vs Control Resistant 
ROD_Susceptible_dds_res <- results(ROD_Susceptible_dds_deseq, alpha=0.05, name= "Condition_Late_Susecptible_vs_Early_Susceptible")
head(ROD_Susceptible_dds_res)  # Condition Late Susecptible vs Early Susceptible 

### Perform LFC Shrinkage with apeglm
## NOTES 
# Before plotting we need to apply an LFC shrinkage, set the coef as the specific comparison in the ResultsNames function of the deseq object
# Issue that the specific comparisons I set in my results formulas are not available in the ResultsNames coef list.
# More detailed notes about LFC Shrinkage are in the code for the Zhang Vibrio

## DECISION: USE SAME RES OBJECT TO KEEP ALPHA ADJUSTMENT, and use LFCShrink apeglm
ROD_Resistant_dds_res_LFC<- lfcShrink(ROD_Resistant_dds_deseq, coef="Condition_Resistant_Challenge_vs_Control_Resistant", type="apeglm", res= ROD_Resistant_dds_res)
ROD_Susceptible_dds_res_LFC <- lfcShrink(ROD_Susceptible_dds_deseq, coef="Condition_Late_Susecptible_vs_Early_Susceptible", type="apeglm", res=ROD_Susceptible_dds_res)

## EXPLORATORY PLOTTING OF RESULTS 
## MA Plotting
plotMA(ROD_Resistant_dds_res_LFC, ylim = c(-5, 5))
plotMA(ROD_Susceptible_dds_res_LFC, ylim = c(-5, 5))

## Histogram of P values 
# exclude genes with very small counts to avoid spikes and plot using the LFCshrinkage
hist(ROD_Resistant_dds_res_LFC$padj[ROD_Resistant_dds_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
hist(ROD_Susceptible_dds_res_LFC$padj[ROD_Susceptible_dds_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

### Subsetting Significant Genes by padj < 0.05
# again, only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
ROD_Resistant_dds_res_LFC_sig <-  subset(ROD_Resistant_dds_res_LFC , padj < 0.05)
ROD_Resistant_dds_res_LFC_sig$ID<- row.names(ROD_Resistant_dds_res_LFC_sig)
ROD_Resistant_dds_res_LFC_sig <- as.data.frame(ROD_Resistant_dds_res_LFC_sig)
nrow(ROD_Resistant_dds_res_LFC_sig) # 68

ROD_Susceptible_dds_res_LFC_sig <-  subset(ROD_Susceptible_dds_res_LFC , padj < 0.05)
ROD_Susceptible_dds_res_LFC_sig$ID<- row.names(ROD_Susceptible_dds_res_LFC_sig)
ROD_Susceptible_dds_res_LFC_sig <- as.data.frame(ROD_Susceptible_dds_res_LFC_sig)
nrow(ROD_Susceptible_dds_res_LFC_sig)  #2020

### GENE CLUSTERING ANALYSIS HEATMAPS  
# Extract genes with the highest variance across samples for each comparison using either vst or rlog transformed data
# This heatmap rather than plotting absolute expression strength plot the amount by which each gene deviates in a specific sample from the gene’s average across all samples. 
# example codes from RNAseq workflow: https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#other-comparisons
ROD_Resistant_dds_rlog
ROD_Susceptible_dds_rlog

ROD_Resistant_dds_LFC_assay <-  head(order(rowVars(assay(ROD_Resistant_dds_rlog  )), decreasing = TRUE), 200)
ROD_Resistant_mat <- assay(ROD_Resistant_dds_rlog )[ROD_Resistant_dds_LFC_assay,]
ROD_Resistant_mat <- ROD_Resistant_mat - rowMeans(ROD_Resistant_mat)
ROD_Resistant_anno <- as.data.frame(colData(ROD_Resistant_dds_rlog )[, c("Condition","Time")])
ROD_Resistant_heatmap <- pheatmap(ROD_Resistant_mat , annotation_col = ROD_Resistant_anno)
head(ROD_Resistant_mat) # some clustering by bacillus, still overall clustering by day

ROD_Susceptible_dds_LFC_assay <-  head(order(rowVars(assay(ROD_Susceptible_dds_rlog  )), decreasing = TRUE), 200)
ROD_Susceptible_mat <- assay(ROD_Susceptible_dds_rlog )[ROD_Susceptible_dds_LFC_assay,]
ROD_Susceptible_mat <- ROD_Susceptible_mat - rowMeans(ROD_Susceptible_mat)
ROD_Susceptible_anno <- as.data.frame(colData(ROD_Susceptible_dds_rlog )[, c("Condition","Time")])
ROD_Susceptible_heatmap <- pheatmap(ROD_Susceptible_mat , annotation_col = ROD_Susceptible_anno)
head(ROD_Susceptible_mat) # some clustering by bacillus, still overall clustering by day

# Gene clustering heatmap with only apoptosis genes #
# vector C_vir_apop transcript IDs
# Search original counts for apoptosis genes and do rlog on just these
ROD_Resistant_counts_apop <- ROD_Resistant_counts[row.names(ROD_Resistant_counts) %in% C_vir_rtracklayer_apop_product_final_ID,]
ROD_Susceptible_counts_apop <- ROD_Susceptible_counts[row.names(ROD_Susceptible_counts) %in% C_vir_rtracklayer_apop_product_final_ID,]
nrow(ROD_Resistant_counts_apop) #1026
nrow( ROD_Susceptible_counts_apop) # 1026

ROD_Resistant_counts_apop_dds <- DESeqDataSetFromMatrix(countData = ROD_Resistant_counts_apop,
                                                    colData = ROD_Resistant_coldata,
                                                    design = ~Time + Condition) # add time to control for time effect
ROD_Susceptible_counts_apop_dds <- DESeqDataSetFromMatrix(countData = ROD_Susceptible_counts_apop,
                                                        colData = ROD_Susceptible_coldata,
                                                        design = ~Condition) # add time to control for time effect

# Prefiltering the data and running rlog
ROD_Resistant_counts_apop_dds<- ROD_Resistant_counts_apop_dds[ rowSums(counts(ROD_Resistant_counts_apop_dds)) > 10, ]
ROD_Susceptible_counts_apop_dds<- ROD_Susceptible_counts_apop_dds[ rowSums(counts(ROD_Susceptible_counts_apop_dds)) > 10, ]

ROD_Resistant_counts_apop_dds_rlog <- rlog(ROD_Resistant_counts_apop_dds, blind=TRUE)
ROD_Susceptible_counts_apop_dds_rlog <- rlog(ROD_Susceptible_counts_apop_dds, blind=TRUE)


## PCA plot of rlog transformed counts for apoptosis
plotPCA(ROD_Resistant_counts_apop_dds_rlog  , intgroup="Time") # still overall clustering by time and not condition
plotPCA(ROD_Susceptible_counts_apop_dds_rlog  , intgroup="Condition") # still overall clustering by time and not condition

# heatmap of all apoptosis genes 
ROD_Resistant_apop_assay <-  assay(ROD_Resistant_counts_apop_dds_rlog)[,]
ROD_Resistant_apop_assay_mat <- ROD_Resistant_apop_assay - rowMeans(ROD_Resistant_apop_assay)
ROD_Resistant_apop_assay_anno <- as.data.frame(colData(ROD_Resistant_counts_apop_dds_rlog )[, c("Condition","Time")])
ROD_Resistant_apop_assay_heatmap <- pheatmap(ROD_Resistant_apop_assay_mat  , annotation_col = ROD_Resistant_apop_assay_anno)
head(ROD_Resistant_apop_assay_mat ) # tighter grouping of challenge and control

ROD_Susceptible_apop_assay <-  assay(ROD_Susceptible_counts_apop_dds_rlog)[,]
ROD_Susceptible_apop_assay_mat <- ROD_Susceptible_apop_assay - rowMeans(ROD_Susceptible_apop_assay)
ROD_Susceptible_apop_assay_anno <- as.data.frame(colData(ROD_Susceptible_counts_apop_dds_rlog )[, c("Condition","Time")])
ROD_Susceptible_apop_assay_heatmap <- pheatmap(ROD_Susceptible_apop_assay_mat  , annotation_col = ROD_Susceptible_apop_assay_anno)
head(ROD_Susceptible_apop_assay_mat ) # early and late cluster perfectly

# heatmap of most variable apoptosis genes for Resistant family (this selects genes with the greatest variance in the sample)
topVarGenes_ROD_Resistant_apop_assay <-  head(order(rowVars(assay(ROD_Resistant_counts_apop_dds_rlog)), decreasing = TRUE), 100) 
top_Var_ROD_Resistant_apop_assay_mat<- assay(ROD_Resistant_counts_apop_dds_rlog)[topVarGenes_ROD_Resistant_apop_assay,]
top_Var_ROD_Resistant_apop_assay_mat <- top_Var_ROD_Resistant_apop_assay_mat - rowMeans(top_Var_ROD_Resistant_apop_assay_mat)
top_Var_ROD_Resistant_apop_assay_anno <- as.data.frame(colData(ROD_Resistant_counts_apop_dds_rlog)[, c("Condition","Time")])
top_Var_ROD_Resistant_apop_assay_heatmap <- pheatmap(top_Var_ROD_Resistant_apop_assay_mat  , annotation_col = top_Var_ROD_Resistant_apop_assay_anno)
head(top_Var_ROD_Resistant_apop_assay_mat )

# annotate the top 100 most variable genes  
# reorder annotation table to match ordering in heatmap 
top_Var_ROD_Resistant_apop_assay_heatmap_reorder <-rownames(top_Var_ROD_Resistant_apop_assay_mat[top_Var_ROD_Resistant_apop_assay_heatmap$tree_row[["order"]],])
# annotate the row.names
top_Var_ROD_Resistant_apop_assay_prot <- as.data.frame(top_Var_ROD_Resistant_apop_assay_heatmap_reorder)
colnames(top_Var_ROD_Resistant_apop_assay_prot)[1] <- "ID"
top_Var_ROD_Resistant_apop_assay_prot_annot <- left_join(top_Var_ROD_Resistant_apop_assay_prot, select(C_vir_rtracklayer_apop_product_final, ID, product, gene), by = "ID")

# heatmap of most variable apoptosis genes for Susceptible family (this selects genes with the greatest variance in the sample)
topVarGenes_ROD_Susceptible_apop_assay <-  head(order(rowVars(assay(ROD_Susceptible_counts_apop_dds_rlog)), decreasing = TRUE), 100) 
top_Var_ROD_Susceptible_apop_assay_mat<- assay(ROD_Susceptible_counts_apop_dds_rlog)[topVarGenes_ROD_Susceptible_apop_assay,]
top_Var_ROD_Susceptible_apop_assay_mat <- top_Var_ROD_Susceptible_apop_assay_mat - rowMeans(top_Var_ROD_Susceptible_apop_assay_mat)
top_Var_ROD_Susceptible_apop_assay_anno <- as.data.frame(colData(ROD_Susceptible_counts_apop_dds_rlog)[, c("Condition","Time")])
top_Var_ROD_Susceptible_apop_assay_heatmap <- pheatmap(top_Var_ROD_Susceptible_apop_assay_mat  , annotation_col = top_Var_ROD_Susceptible_apop_assay_anno)
head(top_Var_ROD_Susceptible_apop_assay_mat )

# annotate the top 100 most variable genes  
# reorder annotation table to match ordering in heatmap 
top_Var_ROD_Susceptible_apop_assay_heatmap_reorder <-rownames(top_Var_ROD_Susceptible_apop_assay_mat[top_Var_ROD_Susceptible_apop_assay_heatmap$tree_row[["order"]],])
# annotate the row.names
top_Var_ROD_Susceptible_apop_assay_prot <- as.data.frame(top_Var_ROD_Susceptible_apop_assay_heatmap_reorder)
colnames(top_Var_ROD_Susceptible_apop_assay_prot)[1] <- "ID"
top_Var_ROD_Susceptible_apop_assay_prot_annot <- left_join(top_Var_ROD_Susceptible_apop_assay_prot, select(C_vir_rtracklayer_apop_product_final, ID, product, gene), by = "ID")

### Extract list of significant Apoptosis Genes (not less than or greater than 1 LFC) using merge
ROD_Susceptible_dds_res_LFC_sig_APOP <- merge(ROD_Susceptible_dds_res_LFC_sig, C_vir_rtracklayer_apop_product_final, by = "ID")
ROD_Susceptible_dds_res_LFC_sig_APOP_arranged <- arrange(ROD_Susceptible_dds_res_LFC_sig_APOP, -log2FoldChange) 
nrow(ROD_Susceptible_dds_res_LFC_sig_APOP) # 39

ROD_Resistant_dds_res_LFC_sig_APOP <- merge(ROD_Resistant_dds_res_LFC_sig, C_vir_rtracklayer_apop_product_final, by = "ID")
ROD_Resistant_dds_res_LFC_sig_APOP_arranged <- arrange(ROD_Resistant_dds_res_LFC_sig_APOP, -log2FoldChange) 
nrow(ROD_Resistant_dds_res_LFC_sig_APOP) # 3

ROD_Susceptible_dds_res_LFC_sig_APOP_plot <- ggplot(ROD_Susceptible_dds_res_LFC_sig_APOP , aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Susceptible Early vs Late") +
  ylab("Log2 Fold Change")

ROD_Resistant_dds_res_LFC_sig_APOP_plot <- ggplot(ROD_Resistant_dds_res_LFC_sig_APOP , aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Susceptible Early vs Late") +
  ylab("Log2 Fold Change")

###### DERMO PROESTOU TRANSCRIPTOME ANALYSIS ####

# DA is susceptible and LB is tolerant

Dermo_counts <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/Dermo_transcript_count_matrix.csv", header=TRUE,
                       row.names = "transcript_id")
head(Dermo_counts)
colnames(Dermo_counts)
ncol(Dermo_counts)

# remove MSTRG novel transcript lines (can assess these later if necessary)
Dermo_counts <- Dermo_counts[!grepl("MSTRG", row.names(Dermo_counts)),]
row.names(Dermo_counts) <- remove_rna(row.names(Dermo_counts))
head(Dermo_counts)

#Load in sample metadata
Dermo_coldata <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/Dermo_coldata.csv", row.names = 1 )
View(Dermo_coldata )  
nrow(Dermo_coldata ) 

# Make sure the columns of the count matrix and rows of the column data (sample metadata) are in the same order. 
all(rownames(Dermo_coldata ) %in% colnames(Dermo_counts))  #Should return TRUE
# returns TRUE
all(colnames(Dermo_counts) %in% rownames(Dermo_coldata ))  
# returns TRUE
all(rownames(Dermo_coldata ) == colnames(Dermo_counts))

# Fix the order
Dermo_counts <- Dermo_counts[,row.names(Dermo_coldata)]

all(rownames(Dermo_coldata ) %in% colnames(Dermo_counts))  #Should return TRUE
# returns TRUE
all(colnames(Dermo_counts) %in% rownames(Dermo_coldata ))  
# returns TRUE
all(rownames(Dermo_coldata ) == colnames(Dermo_counts))
# returns TRUE

# going to split the analysis into DA and LB 
Dermo_Susceptible_coldata <- Dermo_coldata %>%  subset(Family == "DA") # subset() keeps rownames
row.names(Dermo_Susceptible_coldata )
Dermo_Susceptible_counts <- Dermo_counts[,row.names(Dermo_Susceptible_coldata )]
colnames(Dermo_Susceptible_counts)
Dermo_Tolerant_coldata <- Dermo_coldata %>%  subset(Family == "LB") # subset() keeps rownames
Dermo_Tolerant_counts <- Dermo_counts[,row.names(Dermo_Tolerant_coldata)]

### DATA QC PCA PLOT 
# rlog transform data is recommended over vst for small data sets 
# PCA plots of data (https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/02_Preprocessing_Data.nb.html#count-distribution-boxplots)
Dermo_Susceptible_counts_matrix <- as.matrix(Dermo_Susceptible_counts)
Dermo_Tolerant_counts_matrix <-    as.matrix(Dermo_Tolerant_counts)

Dermo_Susceptible_vstcounts <- vst(Dermo_Susceptible_counts_matrix , blind =TRUE)
Dermo_Tolerant_vstcounts <-    vst(Dermo_Tolerant_counts_matrix, blind =TRUE)

# run PCA
pcDermo_Susceptible <- prcomp(t(Dermo_Susceptible_vstcounts))
pcDermo_Tolerant <-   prcomp(t(Dermo_Tolerant_vstcounts))

# Plot PCA
autoplot(pcDermo_Susceptible,
         data = Dermo_Susceptible_coldata, 
         colour="Time", 
         size=5) # PCA axes don't explain much of the variation, there are three clusters but they don't really have a lot in common 
autoplot(pcDermo_Susceptible,
         data = Dermo_Susceptible_coldata, 
         colour="LogConc", 
         size=5)
autoplot(pcDermo_Susceptible,
         data = Dermo_Susceptible_coldata, 
         colour="Condition", 
         size=5)
autoplot(pcDermo_Susceptible,
         data = Dermo_Susceptible_coldata, 
         colour="Tech_rep", 
         size=5) # technical replicates track very closely 
autoplot(pcDermo_Susceptible,
         data = Dermo_Susceptible_coldata, 
         colour="Lib_prep_date", # batches are totally based on library prep date 
         size=5)
autoplot(pcDermo_Susceptible,
         data = Dermo_Susceptible_coldata, 
         colour="LogConc", 
         size=5) # the sample with no ROD signs is a strong outlier, could use disease stage as the basis for comparison here, day 15 and 30 cluster closely together
autoplot(pcDermo_Tolerant,
         data = Dermo_Tolerant_coldata, 
         colour="Time", 
         size=5) #
autoplot(pcDermo_Tolerant,
         data = Dermo_Tolerant_coldata, 
         colour="Condition", 
         size=5) #
autoplot(pcDermo_Tolerant,
         data = Dermo_Tolerant_coldata, 
         colour="Lib_prep_date", 
         size=5) #large separation by lib prep date 

## MAKE DESEQ DATA SET FROM MATRIX
# This object specifies the count data and metadata you will work with. The design piece is critical.
# Correct for batch effects if necessary in this original formula: see this thread https://support.bioconductor.org/p/121408/
# do not correct counts using the removeBatchEffects from limma based on thread above 

## Check levels 
# It is prefered in R that the first level of a factor be the reference level for comparison
# (e.g. control, or untreated samples), so we can relevel the factor like so
# Check factor levels, set it so that comparison group is the first
levels(Dermo_Susceptible_coldata$Condition) #"Control"  "Injected"
levels(Dermo_Tolerant_coldata$Condition) #"Control"  "Injected"
levels(Dermo_Susceptible_coldata$Time) # "28d" "36h" "7d" 
levels(Dermo_Tolerant_coldata$Time) # "28d" "36h" "7d"
levels(Dermo_Tolerant_coldata$Lib_prep_date)
levels(Dermo_Susceptible_coldata$Lib_prep_date)
Dermo_Susceptible_coldata$Time <- factor(Dermo_Susceptible_coldata$Time, levels= c("36h","7d","28d"))
Dermo_Tolerant_coldata$Time <- factor(Dermo_Tolerant_coldata$Time, levels= c("36h","7d","28d"))

## Creating two data set from matrix, one for each family  
Dermo_Tolerant_dds <- DESeqDataSetFromMatrix(countData =Dermo_Tolerant_counts,
                                            colData = Dermo_Tolerant_coldata,
                                            design = ~Lib_prep_date + Condition + Time) # keep time at end so I can look for condition effect through time
Dermo_Susceptible_dds <- DESeqDataSetFromMatrix(countData = Dermo_Susceptible_counts,
                                              colData = Dermo_Susceptible_coldata,
                                              design = ~Lib_prep_date + Condition + Time) # keep time at end so I can look for condition effect through time
## Collapse technical replicates
Dermo_Tolerant_dds <- collapseReplicates(Dermo_Tolerant_dds, Dermo_Tolerant_dds$Sample_ID, Dermo_Tolerant_dds$Tech_rep)
Dermo_Dermo_Susceptible_dds <- collapseReplicates(Dermo_Susceptible_dds, Dermo_Susceptible_dds$Sample_ID, Dermo_Susceptible_dds$Tech_rep)

## Prefiltering the data
# Data prefiltering helps decrease the size of the data set and get rid of
# rows with no data or very minimal data (<10). Apply a minimal filtering here as more stringent filtering will be applied later
Dermo_Tolerant_dds <- Dermo_Tolerant_dds[ rowSums(counts(Dermo_Tolerant_dds )) > 10, ]
Dermo_Susceptible_dds <- Dermo_Susceptible_dds[rowSums(counts(Dermo_Susceptible_dds )) > 10, ]

## DATA TRANSFORMATION AND VISUALIZATION
# Assess sample clustering after setting initial formula for comparison
Dermo_Tolerant_dds_vst <- vst(Dermo_Tolerant_dds , blind = TRUE) # keep blind = true before deseq function has been run
Dermo_Susceptible_dds_vst <- vst(Dermo_Susceptible_dds , blind = TRUE) # keep blind = true before deseq function has been run

## PCA plot visualization of individuals in the family 
plotPCA(Dermo_Tolerant_dds_vst, intgroup= c("Time","Condition")) # no some clustering by treatment
plotPCA(Dermo_Susceptible_dds_vst , intgroup=c("Time","Condition")) # some clustering

### DIFFERENTIAL EXPRESSION ANALYSIS
# run pipeline with single command because the formula has already been specified
# Steps: estimation of size factors (controlling for differences in the sequencing depth of the samples), 
# the estimation of dispersion values for each gene,and fitting a generalized linear model.
Dermo_Tolerant_dds_deseq <- DESeq(Dermo_Tolerant_dds) 
Dermo_Susceptible_dds_deseq <- DESeq(Dermo_Susceptible_dds) 

## Check the resultsNames object of each to look at the available coefficients for use in lfcShrink command
resultsNames(Dermo_Tolerant_dds_deseq) # [1] "Intercept" "Lib_prep_date_12_17_2020_vs_12_15_2020" "Condition_Injected_vs_Control"         
  # [4] "Time_7d_vs_36h"    "Time_28d_vs_36h"        
resultsNames(Dermo_Susceptible_dds_deseq ) #[1] "Intercept" "Lib_prep_date_12_17_2020_vs_12_15_2020" "Condition_Injected_vs_Control"         
  # [4] "Time_7d_vs_36h"   "Time_28d_vs_36h" 

## BUILD THE RESULTS OBJECT
# Examining the results object, change alpha to p <0.05, looking at object metadata
# use mcols to look at metadata for each table
mcols(Dermo_Tolerant_dds_deseq)
mcols(Dermo_Susceptible_dds_deseq)

Dermo_Tolerant_dds_7d_res <- results(Dermo_Tolerant_dds_deseq, alpha=0.05, name="Time_7d_vs_36h" )
Dermo_Tolerant_dds_28d_res <- results(Dermo_Tolerant_dds_deseq, alpha=0.05, name="Time_28d_vs_36h"  )
Dermo_Susceptible_dds_7d_res <- results(Dermo_Susceptible_dds_deseq, alpha=0.05, name="Time_7d_vs_36h" )
Dermo_Susceptible_dds_28d_res <- results(Dermo_Susceptible_dds_deseq, alpha=0.05, name="Time_28d_vs_36h" )

### Perform LFC Shrinkage with apeglm
## NOTES 
# Before plotting we need to apply an LFC shrinkage, set the coef as the specific comparison in the ResultsNames function of the deseq object
# Issue that the specific comparisons I set in my results formulas are not available in the ResultsNames coef list.
# More detailed notes about LFC Shrinkage are in the code for the Zhang Vibrio

## DECISION: USE SAME RES OBJECT TO KEEP ALPHA ADJUSTMENT, and use LFCShrink apeglm
Dermo_Tolerant_dds_7d_res_LFC <- lfcShrink(Dermo_Tolerant_dds_deseq, coef="Time_7d_vs_36h", type="apeglm", res= Dermo_Tolerant_dds_7d_res)
Dermo_Tolerant_dds_28d_res_LFC <- lfcShrink(Dermo_Tolerant_dds_deseq, coef="Time_28d_vs_36h", type="apeglm", res= Dermo_Tolerant_dds_28d_res)
Dermo_Susceptible_dds_7d_res_LFC <- lfcShrink(Dermo_Susceptible_dds_deseq, coef="Time_7d_vs_36h", type="apeglm", res=Dermo_Susceptible_dds_7d_res)
Dermo_Susceptible_dds_28d_res_LFC <- lfcShrink(Dermo_Susceptible_dds_deseq, coef="Time_28d_vs_36h", type="apeglm", res=Dermo_Susceptible_dds_28d_res)

## EXPLORATORY PLOTTING OF RESULTS 
## MA Plotting
plotMA(Dermo_Tolerant_dds_7d_res_LFC, ylim = c(-5, 5))
plotMA(Dermo_Tolerant_dds_28d_res_LFC, ylim = c(-5, 5))
plotMA(Dermo_Susceptible_dds_7d_res_LFC, ylim = c(-5, 5))
plotMA(Dermo_Susceptible_dds_28d_res_LFC, ylim = c(-5, 5))

## Histogram of P values 
# exclude genes with very small counts to avoid spikes and plot using the LFCshrinkage
hist(Dermo_Tolerant_dds_7d_res_LFC$padj[Dermo_Tolerant_dds_7d_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
hist(Dermo_Tolerant_dds_28d_res_LFC$padj[Dermo_Tolerant_dds_28d_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
hist(Dermo_Susceptible_dds_7d_res_LFC$padj[Dermo_Susceptible_dds_7d_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
hist(Dermo_Susceptible_dds_28d_res_LFC$padj[Dermo_Susceptible_dds_28d_res_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

### Subsetting Significant Genes by padj < 0.05
# again, only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
Dermo_Tolerant_dds_7d_res_LFC_sig <-  subset(Dermo_Tolerant_dds_7d_res_LFC , padj < 0.05)
Dermo_Tolerant_dds_7d_res_LFC_sig$ID<- row.names(Dermo_Tolerant_dds_7d_res_LFC_sig)
Dermo_Tolerant_dds_7d_res_LFC_sig <- as.data.frame(Dermo_Tolerant_dds_7d_res_LFC_sig)
nrow(Dermo_Tolerant_dds_7d_res_LFC_sig) #1390

Dermo_Tolerant_dds_28d_res_LFC_sig <-  subset(Dermo_Tolerant_dds_28d_res_LFC , padj < 0.05)
Dermo_Tolerant_dds_28d_res_LFC_sig$ID<- row.names(Dermo_Tolerant_dds_28d_res_LFC_sig)
Dermo_Tolerant_dds_28d_res_LFC_sig <- as.data.frame(Dermo_Tolerant_dds_28d_res_LFC_sig)
nrow(Dermo_Tolerant_dds_28d_res_LFC_sig) # 1848

Dermo_Susceptible_dds_7d_res_LFC_sig <-  subset(Dermo_Susceptible_dds_7d_res_LFC , padj < 0.05)
Dermo_Susceptible_dds_7d_res_LFC_sig$ID<- row.names(Dermo_Susceptible_dds_7d_res_LFC_sig)
Dermo_Susceptible_dds_7d_res_LFC_sig <- as.data.frame(Dermo_Susceptible_dds_7d_res_LFC_sig)
nrow(Dermo_Susceptible_dds_7d_res_LFC_sig)  #2283

Dermo_Susceptible_dds_28d_res_LFC_sig <-  subset(Dermo_Susceptible_dds_28d_res_LFC , padj < 0.05)
Dermo_Susceptible_dds_28d_res_LFC_sig$ID<- row.names(Dermo_Susceptible_dds_28d_res_LFC_sig)
Dermo_Susceptible_dds_28d_res_LFC_sig <- as.data.frame(Dermo_Susceptible_dds_28d_res_LFC_sig)
nrow(Dermo_Susceptible_dds_28d_res_LFC_sig)  #2477

### GENE CLUSTERING ANALYSIS HEATMAPS  
# Extract genes with the highest variance across samples for each comparison using either vst or rlog transformed data
# This heatmap rather than plotting absolute expression strength plot the amount by which each gene deviates in a specific sample from the gene’s average across all samples. 
# example codes from RNAseq workflow: https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#other-comparisons
Dermo_Tolerant_dds_vst
Dermo_Susceptible_dds_vst

Dermo_Tolerant_dds_LFC_assay <-  head(order(rowVars(assay(Dermo_Tolerant_dds_vst  )), decreasing = TRUE), 200)
Dermo_Tolerant_mat <- assay(Dermo_Tolerant_dds_vst )[Dermo_Tolerant_dds_LFC_assay,]
Dermo_Tolerant_mat <- Dermo_Tolerant_mat - rowMeans(Dermo_Tolerant_mat)
Dermo_Tolerant_anno <- as.data.frame(colData(Dermo_Tolerant_dds_vst )[, c("Condition","Time")])
Dermo_Tolerant_heatmap <- pheatmap(Dermo_Tolerant_mat , annotation_col =Dermo_Tolerant_anno)
head(Dermo_Tolerant_mat) # clustering mostly by control and injected 

Dermo_Susceptible_dds_LFC_assay <-  head(order(rowVars(assay(Dermo_Susceptible_dds_vst  )), decreasing = TRUE), 200)
Dermo_Susceptible_mat <- assay(Dermo_Susceptible_dds_vst)[Dermo_Susceptible_dds_LFC_assay,]
Dermo_Susceptible_mat <- Dermo_Susceptible_mat - rowMeans(Dermo_Susceptible_mat)
Dermo_Susceptible_anno <- as.data.frame(colData(Dermo_Susceptible_dds_vst)[, c("Condition","Time")])
Dermo_Susceptible_heatmap <- pheatmap(Dermo_Susceptible_mat , annotation_col =Dermo_Susceptible_anno)
head(Dermo_Susceptible_mat) # clustering by control and injected

# Gene clustering heatmap with only apoptosis genes #
# vector C_vir_apop transcript IDs
# Search original counts for apoptosis genes and do rlog on just these
Dermo_Tolerant_counts_apop <- Dermo_Tolerant_counts[row.names(Dermo_Tolerant_counts) %in% C_vir_rtracklayer_apop_product_final_ID,]
Dermo_Susceptible_counts_apop <- Dermo_Susceptible_counts[row.names(Dermo_Susceptible_counts) %in% C_vir_rtracklayer_apop_product_final_ID,]
nrow(Dermo_Tolerant_counts_apop) #1026
nrow( Dermo_Susceptible_counts_apop) # 1026

Dermo_Tolerant_counts_apop_dds <- DESeqDataSetFromMatrix(countData = Dermo_Tolerant_counts_apop,
                                                        colData = Dermo_Tolerant_coldata,
                                                        design = ~Lib_prep_date + Condition + Time) # same as before
Dermo_Susceptible_counts_apop_dds <- DESeqDataSetFromMatrix(countData = Dermo_Susceptible_counts_apop,
                                                          colData =Dermo_Susceptible_coldata,
                                                          design =  ~Lib_prep_date + Condition + Time) # same as before

# Prefiltering the data and running rlog
Dermo_Tolerant_counts_apop_dds<- Dermo_Tolerant_counts_apop_dds[ rowSums(counts(Dermo_Tolerant_counts_apop_dds)) > 10, ]
Dermo_Susceptible_counts_apop_dds<- Dermo_Susceptible_counts_apop_dds[ rowSums(counts(Dermo_Susceptible_counts_apop_dds)) > 10, ]

Dermo_Tolerant_counts_apop_dds_vst <- varianceStabilizingTransformation(Dermo_Tolerant_counts_apop_dds, blind=TRUE)
Dermo_Susceptible_counts_apop_dds_vst <- varianceStabilizingTransformation(Dermo_Susceptible_counts_apop_dds, blind=TRUE)

## PCA plot of rlog transformed counts for apoptosis
plotPCA(Dermo_Tolerant_counts_apop_dds_vst,  intgroup="Condition") # mostly clustering by condition
plotPCA(Dermo_Susceptible_counts_apop_dds_vst  , intgroup="Time") # some clustering by condition

# heatmap of all apoptosis genes 
Dermo_Tolerant_apop_assay <-  assay(Dermo_Tolerant_counts_apop_dds_vst)[,]
Dermo_Tolerant_apop_assay_mat <- Dermo_Tolerant_apop_assay - rowMeans(Dermo_Tolerant_apop_assay)
Dermo_Tolerant_apop_assay_anno <- as.data.frame(colData(Dermo_Tolerant_counts_apop_dds_vst )[, c("Condition","Time")])
Dermo_Tolerant_apop_assay_heatmap <- pheatmap(Dermo_Tolerant_apop_assay_mat  , annotation_col = Dermo_Tolerant_apop_assay_anno)
head(Dermo_Tolerant_apop_assay_mat ) # some clustering by treatment and condition

Dermo_Susceptible_apop_assay <-  assay(Dermo_Susceptible_counts_apop_dds_vst)[,]
Dermo_Susceptible_apop_assay_mat <- Dermo_Susceptible_apop_assay - rowMeans(Dermo_Susceptible_apop_assay)
Dermo_Susceptible_apop_assay_anno <- as.data.frame(colData(Dermo_Susceptible_counts_apop_dds_vst )[, c("Condition","Time")])
Dermo_Susceptible_apop_assay_heatmap <- pheatmap(Dermo_Susceptible_apop_assay_mat  , annotation_col = Dermo_Susceptible_apop_assay_anno)
head(Dermo_Susceptible_apop_assay_mat ) #some clustering by treatment and condition

# heatmap of most variable apoptosis genes for Resistant family (this selects genes with the greatest variance in the sample)
topVarGenes_Dermo_Tolerant_apop_assay <-  head(order(rowVars(assay(Dermo_Tolerant_counts_apop_dds_vst)), decreasing = TRUE), 100) 
top_Var_Dermo_Tolerant_apop_assay_mat<- assay(Dermo_Tolerant_counts_apop_dds_vst)[topVarGenes_Dermo_Tolerant_apop_assay,]
top_Var_Dermo_Tolerant_apop_assay_mat <- top_Var_Dermo_Tolerant_apop_assay_mat - rowMeans(top_Var_Dermo_Tolerant_apop_assay_mat)
top_Var_Dermo_Tolerant_apop_assay_anno <- as.data.frame(colData(Dermo_Tolerant_counts_apop_dds_vst)[, c("Condition","Time")])
top_Var_Dermo_Tolerant_apop_assay_heatmap <- pheatmap(top_Var_Dermo_Tolerant_apop_assay_mat  , annotation_col = top_Var_Dermo_Tolerant_apop_assay_anno)
head(top_Var_Dermo_Tolerant_apop_assay_mat ) # some clustering by condition, but not great

# annotate the top 100 most variable genes  
# reorder annotation table to match ordering in heatmap 
top_Var_Dermo_Tolerant_apop_assay_heatmap_reorder <-rownames(top_Var_Dermo_Tolerant_apop_assay_mat[top_Var_Dermo_Tolerant_apop_assay_heatmap$tree_row[["order"]],])
# annotate the row.names
top_Var_Dermo_Tolerant_apop_assay_prot <- as.data.frame(top_Var_Dermo_Tolerant_apop_assay_heatmap_reorder)
colnames(top_Var_Dermo_Tolerant_apop_assay_prot)[1] <- "ID"
top_Var_Dermo_Tolerant_apop_assay_prot_annot <- left_join(top_Var_Dermo_Tolerant_apop_assay_prot, select(C_vir_rtracklayer_apop_product_final, ID, product, gene), by = "ID")

# heatmap of most variable apoptosis genes for Susceptible family (this selects genes with the greatest variance in the sample)
topVarGenes_Dermo_Susceptible_apop_assay <-  head(order(rowVars(assay(Dermo_Susceptible_counts_apop_dds_vst)), decreasing = TRUE), 100) 
top_Var_Dermo_Susceptible_apop_assay_mat<- assay(Dermo_Susceptible_counts_apop_dds_vst)[topVarGenes_Dermo_Susceptible_apop_assay,]
top_Var_Dermo_Susceptible_apop_assay_mat <- top_Var_Dermo_Susceptible_apop_assay_mat - rowMeans(top_Var_Dermo_Susceptible_apop_assay_mat)
top_Var_Dermo_Susceptible_apop_assay_anno <- as.data.frame(colData(Dermo_Susceptible_counts_apop_dds_vst)[, c("Condition","Time")])
top_Var_Dermo_Susceptible_apop_assay_heatmap <- pheatmap(top_Var_Dermo_Susceptible_apop_assay_mat  , annotation_col = top_Var_Dermo_Susceptible_apop_assay_anno)
head(top_Var_Dermo_Susceptible_apop_assay_mat )

# annotate the top 100 most variable genes  
# reorder annotation table to match ordering in heatmap 
top_Var_Dermo_Susceptible_apop_assay_heatmap_reorder <-rownames(top_Var_Dermo_Susceptible_apop_assay_mat[top_Var_Dermo_Susceptible_apop_assay_heatmap$tree_row[["order"]],])
# annotate the row.names
top_Var_Dermo_Susceptible_apop_assay_prot <- as.data.frame(top_Var_Dermo_Susceptible_apop_assay_heatmap_reorder)
colnames(top_Var_Dermo_Susceptible_apop_assay_prot)[1] <- "ID"
top_Var_Dermo_Susceptible_apop_assay_prot_annot <- left_join(top_Var_Dermo_Susceptible_apop_assay_prot, select(C_vir_rtracklayer_apop_product_final, ID, product, gene), by = "ID")

### Extract list of significant Apoptosis Genes (not less than or greater than 1 LFC) using merge
Dermo_Susceptible_dds_7d_res_LFC_sig_APOP <- merge(Dermo_Susceptible_dds_7d_res_LFC_sig, C_vir_rtracklayer_apop_product_final, by = "ID")
Dermo_Susceptible_dds_7d_res_LFC_sig_APOP_arranged <- arrange(Dermo_Susceptible_dds_7d_res_LFC_sig_APOP, -log2FoldChange) 
nrow(Dermo_Susceptible_dds_7d_res_LFC_sig_APOP) # 40

Dermo_Susceptible_dds_28d_res_LFC_sig_APOP <- merge(Dermo_Susceptible_dds_28d_res_LFC_sig, C_vir_rtracklayer_apop_product_final, by = "ID")
Dermo_Susceptible_dds_28d_res_LFC_sig_APOP_arranged <- arrange(Dermo_Susceptible_dds_28d_res_LFC_sig_APOP, -log2FoldChange) 
nrow(Dermo_Susceptible_dds_28d_res_LFC_sig_APOP) # 40

Dermo_Tolerant_dds_7d_res_LFC_sig_APOP <- merge(Dermo_Tolerant_dds_7d_res_LFC_sig, C_vir_rtracklayer_apop_product_final, by = "ID")
Dermo_Tolerant_dds_7d_res_LFC_sig_APOP_arranged <- arrange(Dermo_Tolerant_dds_7d_res_LFC_sig_APOP, -log2FoldChange) 
nrow(Dermo_Tolerant_dds_7d_res_LFC_sig_APOP) # 17

Dermo_Tolerant_dds_28d_res_LFC_sig_APOP <- merge(Dermo_Tolerant_dds_28d_res_LFC_sig, C_vir_rtracklayer_apop_product_final, by = "ID")
Dermo_Tolerant_dds_28d_res_LFC_sig_APOP_arranged <- arrange(Dermo_Tolerant_dds_28d_res_LFC_sig_APOP, -log2FoldChange) 
nrow(Dermo_Tolerant_dds_28d_res_LFC_sig_APOP) # 31

# Combined LFC plot
Dermo_Susceptible_dds_7d_res_LFC_sig_APOP
Dermo_Susceptible_dds_28d_res_LFC_sig_APOP
Dermo_Tolerant_dds_7d_res_LFC_sig_APOP
Dermo_Tolerant_dds_28d_res_LFC_sig_APOP

# Compare apoptosis genes between group_by_sim groups
Dermo_Susceptible_dds_7d_res_LFC_sig_APOP$group_by_sim <- "Susceptible_dds_7d"
Dermo_Susceptible_dds_28d_res_LFC_sig_APOP$group_by_sim <- "Susceptible_dds_28d"
Dermo_Tolerant_dds_7d_res_LFC_sig_APOP$group_by_sim <- "Tolerant_dds_7d"
Dermo_Tolerant_dds_28d_res_LFC_sig_APOP$group_by_sim <- "Tolerant_dds_28d"

# combine data frames 
Dermo_all_sig_APOP <- rbind( 
  Dermo_Susceptible_dds_7d_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
  Dermo_Susceptible_dds_28d_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
  Dermo_Tolerant_dds_7d_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
  Dermo_Tolerant_dds_28d_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")]) 
  
# Make plot or up and downregulated
Dermo_all_sig_APOP_downregulated <- Dermo_all_sig_APOP%>% filter(log2FoldChange <= 0)
Dermo_all_sig_APOP_upregulated <-   Dermo_all_sig_APOP%>% filter(log2FoldChange > 0)

Dermo_all_sig_APOP_downregulated_plot <- ggplot(Dermo_all_sig_APOP_downregulated , aes(x=product,y=log2FoldChange, fill=group_by_sim )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()
Dermo_all_sig_APOP_upregulated_plot <- ggplot(Dermo_all_sig_APOP_upregulated, aes(x=product,y=log2FoldChange, fill=group_by_sim )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()


Dermo_Susceptible_dds_7d_res_LFC_sig_APOP_plot <- ggplot(Dermo_Susceptible_dds_7d_res_LFC_sig_APOP , aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Dermo Susceptible 7d vs 36hr") +
  ylab("Log2 Fold Change")
Dermo_Susceptible_dds_28d_res_LFC_sig_APOP_plot <- ggplot(Dermo_Susceptible_dds_28d_res_LFC_sig_APOP , aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Dermo Susceptible 28d vs 36hr") +
  ylab("Log2 Fold Change")
Dermo_Tolerant_dds_7d_res_LFC_sig_APOP_plot <- ggplot(Dermo_Tolerant_dds_7d_res_LFC_sig_APOP , aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Dermo Tolerant 7d vs 36hr") +
  ylab("Log2 Fold Change")
Dermo_Tolerant_dds_28d_res_LFC_sig_APOP_plot <- ggplot(Dermo_Tolerant_dds_28d_res_LFC_sig_APOP , aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Dermo Tolerant 28d vs 36hr") +
  ylab("Log2 Fold Change")


#### LFC PLOTS BY SPECIES ####

ROD_Susceptible_dds_res_LFC_sig_APOP$group_by_sim <- "ROD_susceptible"
ROD_Susceptible_dds_res_LFC_sig_APOP$experiment <- "ROD"
ROD_Resistant_dds_res_LFC_sig_APOP$group_by_sim <- "ROD_resistant"
ROD_Resistant_dds_res_LFC_sig_APOP$experiment <- "ROD"
Probiotic_dds_deseq_Challenge_res_LFC_sig_APOP$group_by_sim <- "Probiotic"
Probiotic_dds_deseq_Challenge_res_LFC_sig_APOP$experiment <- "Probiotic"
Dermo_Susceptible_dds_7d_res_LFC_sig_APOP$experiment <- "Dermo"
Dermo_Susceptible_dds_28d_res_LFC_sig_APOP$experiment <- "Dermo"
Dermo_Tolerant_dds_7d_res_LFC_sig_APOP$experiment <- "Dermo"
Dermo_Tolerant_dds_28d_res_LFC_sig_APOP$experiment <- "Dermo"
Pro_RE22_dds_deseq_res_RI_6h_LFC_sig_APOP$experiment <- "Pro_RE22"
Pro_RE22_dds_deseq_res_RI_24h_LFC_sig_APOP$experiment <- "Pro_RE22"
Pro_RE22_dds_deseq_res_S4_6h_LFC_sig_APOP$experiment <- "Pro_RE22"
Pro_RE22_dds_deseq_res_S4_24h_LFC_sig_APOP$experiment <- "Pro_RE22"
Pro_RE22_dds_deseq_res_RE22_LFC_sig_APOP$experiment <- "RE22"

# combine all data 
C_vir_apop_LFC <- rbind(Dermo_Susceptible_dds_7d_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange", "experiment")],
Dermo_Susceptible_dds_28d_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
Dermo_Tolerant_dds_7d_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
Dermo_Tolerant_dds_28d_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
ROD_Susceptible_dds_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")], 
Probiotic_dds_deseq_Challenge_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
Pro_RE22_dds_deseq_res_RI_6h_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
Pro_RE22_dds_deseq_res_RI_24h_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
Pro_RE22_dds_deseq_res_S4_6h_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
Pro_RE22_dds_deseq_res_S4_24h_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
Pro_RE22_dds_deseq_res_RE22_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")])

# Make plot for up and downregulated
C_vir_apop_APOP_downregulated <- C_vir_apop_LFC %>% filter(log2FoldChange <= 0)
C_vir_apop_APOP_upregulated <- C_vir_apop_LFC %>% filter(log2FoldChange > 0)

# Make plot
C_vir_apop_APOP_downregulated_plot <- ggplot(C_vir_apop_APOP_downregulated , aes(x=product,y=log2FoldChange, fill=experiment )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()
C_vir_apop_APOP_upregulated_plot <- ggplot(C_vir_apop_APOP_upregulated , aes(x=product,y=log2FoldChange, fill=experiment )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()

# C_gig
Zhang_dds_deseq_res_V_alg1_APOP_short$experiment <- "Zhang"
Zhang_dds_deseq_res_V_tub_APOP_short$experiment <- "Zhang"
Zhang_dds_deseq_res_LPS_APOP_short$experiment <- "Zhang"
Rubio_dds_deseq_J2_8_res_LFC_sig_APOP_short$experiment <- "Rubio"
Rubio_dds_deseq_J2_9_res_LFC_sig_APOP_short$experiment <- "Rubio"
Rubio_dds_deseq_LGP32_res_LFC_sig_APOP_short$experiment <- "Rubio"
Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP_short$experiment <- "Rubio"
He_dds_res_6hr_sig_APOP$experiment <- "He"
He_dds_res_12hr_sig_APOP$experiment <- "He"
He_dds_res_24hr_sig_APOP$experiment <- "He"
He_dds_res_48hr_sig_APOP$experiment <- "He"
He_dds_res_120hr_sig_APOP$experiment <- "He"
deLorgeril_Resistant_dds_res_6_LFC_sig_APOP$experiment <- "deLorgeril"
deLorgeril_Resistant_dds_res_12_LFC_sig_APOP$experiment <- "deLorgeril"
deLorgeril_Resistant_dds_res_24_LFC_sig_APOP$experiment <- "deLorgeril"
deLorgeril_Resistant_dds_res_48_LFC_sig_APOP$experiment <- "deLorgeril"
deLorgeril_Resistant_dds_res_60_LFC_sig_APOP$experiment <- "deLorgeril"
deLorgeril_Resistant_dds_res_72_LFC_sig_APOP$experiment <- "deLorgeril"
deLorgeril_Susceptible_dds_res_6_LFC_sig_APOP$experiment <- "deLorgeril"
deLorgeril_Susceptible_dds_res_12_LFC_sig_APOP$experiment <- "deLorgeril"
deLorgeril_Susceptible_dds_res_24_LFC_sig_APOP$experiment <- "deLorgeril"
deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP$experiment <- "deLorgeril"
deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP$experiment <- "deLorgeril"
deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP$experiment <- "deLorgeril"

C_gig_apop_LFC <- rbind(Zhang_dds_deseq_res_V_alg1_APOP_short[,c("product","group_by_sim","log2FoldChange","experiment")],
                        Zhang_dds_deseq_res_V_tub_APOP_short[,c("product","group_by_sim","log2FoldChange","experiment")],
                        Zhang_dds_deseq_res_LPS_APOP_short[,c("product","group_by_sim","log2FoldChange","experiment")],
                        Rubio_dds_deseq_J2_8_res_LFC_sig_APOP_short[,c("product","group_by_sim","log2FoldChange","experiment")],
                        Rubio_dds_deseq_J2_9_res_LFC_sig_APOP_short[,c("product","group_by_sim","log2FoldChange","experiment")],
                        Rubio_dds_deseq_LGP32_res_LFC_sig_APOP_short[,c("product","group_by_sim","log2FoldChange","experiment")],
                        Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP_short[,c("product","group_by_sim","log2FoldChange","experiment")],
                        He_dds_res_6hr_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
                        He_dds_res_12hr_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
                        He_dds_res_24hr_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
                        He_dds_res_48hr_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
                        He_dds_res_120hr_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
                        deLorgeril_Resistant_dds_res_6_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
                        deLorgeril_Resistant_dds_res_12_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
                        deLorgeril_Resistant_dds_res_24_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
                        deLorgeril_Resistant_dds_res_48_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
                        deLorgeril_Resistant_dds_res_60_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
                        deLorgeril_Resistant_dds_res_72_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
                        deLorgeril_Susceptible_dds_res_6_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
                        deLorgeril_Susceptible_dds_res_12_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
                        deLorgeril_Susceptible_dds_res_24_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
                        deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
                        deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")],
                        deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment")])

# Make plot for up and downregulated
C_gig_apop_APOP_downregulated <- C_gig_apop_LFC %>% filter(log2FoldChange <= 0)
C_gig_apop_APOP_upregulated <- C_gig_apop_LFC %>% filter(log2FoldChange > 0)

# Make plot
C_gig_apop_APOP_downregulated_plot <- ggplot(C_gig_apop_APOP_downregulated , aes(x=product,y=log2FoldChange, fill=experiment )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()
C_gig_apop_APOP_upregulated_plot <- ggplot(C_gig_apop_APOP_upregulated , aes(x=product,y=log2FoldChange, fill=experiment )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()

#### TRANSCRIPT UPSET PLOT ####
# helpful link: https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
# Using UpsetR NEED TO FIX THIS BELOW
# Input data can be a list of sets where each set is a vector: https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html#input-data
Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP_vector <- as.vector(Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP$transcript_id)
Zhang_dds_deseq_res_V_tub_LFC_sig_APOP_vector <- as.vector(Zhang_dds_deseq_res_V_tub_LFC_sig_APOP$transcript_id)
Zhang_dds_deseq_res_LFC_LPS_sig_APOP_vector <- as.vector(Zhang_dds_deseq_res_LFC_LPS_sig_APOP$transcript_id)
Rubio_dds_deseq_J2_8_res_LFC_sig_APOP_vector <- as.vector(Rubio_dds_deseq_J2_8_res_LFC_sig_APOP$transcript_id)
Rubio_dds_deseq_J2_9_res_LFC_sig_APOP_vector <- as.vector(Rubio_dds_deseq_J2_9_res_LFC_sig_APOP$transcript_id)
Rubio_dds_deseq_LGP32_res_LFC_sig_APOP_vector <- as.vector(Rubio_dds_deseq_LGP32_res_LFC_sig_APOP$transcript_id)
Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP_vector <- as.vector(Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP$transcript_id)
He_dds_res_6hr_sig_APOP_vector <- as.vector(He_dds_res_6hr_sig_APOP$transcript_id)
He_dds_res_12hr_sig_APOP_vector <- as.vector(He_dds_res_12hr_sig_APOP$transcript_id)
He_dds_res_24hr_sig_APOP_vector <- as.vector(He_dds_res_24hr_sig_APOP$transcript_id)
He_dds_res_48hr_sig_APOP_vector <- as.vector(He_dds_res_48hr_sig_APOP$transcript_id)
He_dds_res_120hr_sig_APOP_vector <- as.vector(He_dds_res_120hr_sig_APOP$transcript_id)
deLorgeril_Resistant_dds_res_6_LFC_sig_APOP_vector <- as.vector(deLorgeril_Resistant_dds_res_6_LFC_sig_APOP$transcript_id)
deLorgeril_Resistant_dds_res_12_LFC_sig_APOP_vector <- as.vector(deLorgeril_Resistant_dds_res_12_LFC_sig_APOP$transcript_id)
deLorgeril_Resistant_dds_res_24_LFC_sig_APOP_vector <- as.vector(deLorgeril_Resistant_dds_res_24_LFC_sig_APOP$transcript_id)
deLorgeril_Resistant_dds_res_48_LFC_sig_APOP_vector <- as.vector(deLorgeril_Resistant_dds_res_48_LFC_sig_APOP$transcript_id)
deLorgeril_Resistant_dds_res_60_LFC_sig_APOP_vector <- as.vector(deLorgeril_Resistant_dds_res_60_LFC_sig_APOP$transcript_id)
deLorgeril_Resistant_dds_res_72_LFC_sig_APOP_vector <- as.vector(deLorgeril_Resistant_dds_res_72_LFC_sig_APOP$transcript_id)
deLorgeril_Susceptible_dds_res_6_LFC_sig_APOP_vector <- as.vector(deLorgeril_Susceptible_dds_res_6_LFC_sig_APOP$transcript_id)
deLorgeril_Susceptible_dds_res_12_LFC_sig_APOP_vector <- as.vector(deLorgeril_Susceptible_dds_res_12_LFC_sig_APOP$transcript_id)
deLorgeril_Susceptible_dds_res_24_LFC_sig_APOP_vector <- as.vector(deLorgeril_Susceptible_dds_res_24_LFC_sig_APOP$transcript_id)
deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP_vector <- as.vector(deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP$transcript_id)
deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP_vector <- as.vector(deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP$transcript_id)
deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP_vector <- as.vector(deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP$transcript_id)

# make list of vectors
C_gig_transcript_list <- 
  list(Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP_vector = Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP_vector,
       Zhang_dds_deseq_res_V_tub_LFC_sig_APOP_vector = Zhang_dds_deseq_res_V_tub_LFC_sig_APOP_vector,
       Zhang_dds_deseq_res_LFC_LPS_sig_APOP_vector = Zhang_dds_deseq_res_LFC_LPS_sig_APOP_vector,
       Rubio_dds_deseq_J2_8_res_LFC_sig_APOP_vector = Rubio_dds_deseq_J2_8_res_LFC_sig_APOP_vector,
       Rubio_dds_deseq_J2_9_res_LFC_sig_APOP_vector = Rubio_dds_deseq_J2_9_res_LFC_sig_APOP_vector,
       Rubio_dds_deseq_LGP32_res_LFC_sig_APOP_vector = Rubio_dds_deseq_LGP32_res_LFC_sig_APOP_vector,
       Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP_vector = Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP_vector,
       He_dds_res_6hr_sig_APOP_vector = He_dds_res_6hr_sig_APOP_vector,
       He_dds_res_12hr_sig_APOP_vector = He_dds_res_12hr_sig_APOP_vector,
       He_dds_res_24hr_sig_APOP_vector = He_dds_res_24hr_sig_APOP_vector,
       He_dds_res_48hr_sig_APOP_vector = He_dds_res_48hr_sig_APOP_vector,
       He_dds_res_120hr_sig_APOP_vector = He_dds_res_120hr_sig_APOP_vector,
       deLorgeril_Resistant_dds_res_6_LFC_sig_APOP_vector = deLorgeril_Resistant_dds_res_6_LFC_sig_APOP_vector,
       deLorgeril_Resistant_dds_res_12_LFC_sig_APOP_vector = deLorgeril_Resistant_dds_res_12_LFC_sig_APOP_vector,
       deLorgeril_Resistant_dds_res_24_LFC_sig_APOP_vector = deLorgeril_Resistant_dds_res_24_LFC_sig_APOP_vector,
       deLorgeril_Resistant_dds_res_48_LFC_sig_APOP_vector = deLorgeril_Resistant_dds_res_48_LFC_sig_APOP_vector,
       deLorgeril_Resistant_dds_res_60_LFC_sig_APOP_vector = deLorgeril_Resistant_dds_res_60_LFC_sig_APOP_vector,
       deLorgeril_Resistant_dds_res_72_LFC_sig_APOP_vector = deLorgeril_Resistant_dds_res_72_LFC_sig_APOP_vector,
       deLorgeril_Susceptible_dds_res_6_LFC_sig_APOP_vector = deLorgeril_Susceptible_dds_res_6_LFC_sig_APOP_vector,
       deLorgeril_Susceptible_dds_res_12_LFC_sig_APOP_vector = deLorgeril_Susceptible_dds_res_12_LFC_sig_APOP_vector,
       deLorgeril_Susceptible_dds_res_24_LFC_sig_APOP_vector = deLorgeril_Susceptible_dds_res_24_LFC_sig_APOP_vector,
       deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP_vector = deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP_vector,
       deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP_vector = deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP_vector,
       deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP_vector = deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP_vector)

# Make combination matrix in intersect mode with the list 
C_gig_transcript_list_matrix_comb <- make_comb_mat(C_gig_transcript_list, mode = "intersect")


# plot in intersect mode 
C_gig_transcript_upset_plot <- UpSet(C_gig_transcript_list_matrix_comb)


#### COMPARING APOPTOSIS TRANSCRIPT EXPRESSION BETWEEN EXPERIMENTS PCA HEATMAPS VST ON APOP SUBSET ALONE ####
# some helpful forum posts on the topic: https://www.biostars.org/p/364768/
# Suggest combining, using limma to remove batch effects for each experiment, and then calculate the rlog all together
# Could combine pvalues from DEseq for comparison using metaRNAseq R package to compare p values with a Fisher method

# Load allcoldata.csv 
All_coldata <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/All_coldata.csv", row.names = 1 )
View(All_coldata)  
C_vir_coldata <- subset(All_coldata, Species =="C_vir")
C_gig_coldata <-  subset(All_coldata, Species =="C_gig")

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



#### COMPARING APOPTOSIS TRANSCRIPT EXPRESSION BETWEEN EXPERIMENTS PCA HEATMAPS VST ON FULL THEN SUBSET ####
# some helpful forum posts on the topic: https://www.biostars.org/p/364768/
# Suggest combining, using limma to remove batch effects for each experiment, and then calculate the rlog all together
# Could combine pvalues from DEseq for comparison using metaRNAseq R package to compare p values with a Fisher method

# Load allcoldata.csv 
All_coldata <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/All_coldata.csv", row.names = 1 )
View(All_coldata)  
C_vir_coldata <- subset(All_coldata, Species =="C_vir")
C_gig_coldata <-  subset(All_coldata, Species =="C_gig")

# Combine transcript count matrices by rownames
nrow(Dermo_counts) # 67868
nrow(Probiotic_counts) # 67876
nrow(ROD_counts) # 67870
Dermo_counts$rownames<- row.names(Dermo_counts)
Probiotic_counts$rownames<- row.names(Probiotic_counts)
ROD_counts$rownames<- row.names(ROD_counts)
 
nrow(deLorgeril_counts) # 53701
nrow(He_counts) # 86859
nrow(Zhang_counts ) # 53705
nrow(Rubio_counts) # 53705
deLorgeril_counts$rownames <- row.names(deLorgeril_counts)
He_counts$rownames <- row.names(He_counts)
Zhang_counts $rownames <- row.names(Zhang_counts )
Rubio_counts$rownames <- row.names(Rubio_counts)

# merge based on rownames (then delete rownames), starting with largest first
C_vir_full_counts <- left_join(Probiotic_counts,ROD_counts, by ="rownames")
C_vir_full_counts <- left_join(C_vir_full_counts,Dermo_counts, by = "rownames")
colnames(C_vir_full_counts)
row.names(C_vir_full_counts) <- C_vir_full_counts$rownames
head(C_vir_full_counts)
C_vir_full_counts <- C_vir_full_counts[,-7] # remove rownames to allow for vst 

C_gig_full_counts <- left_join(He_counts,Zhang_counts, by ="rownames")
C_gig_full_counts <- left_join(C_gig_full_counts,Rubio_counts, by = "rownames")
C_gig_full_counts <- left_join(C_gig_full_counts,deLorgeril_counts, by = "rownames")
colnames(C_gig_full_counts)
row.names(C_gig_full_counts) <- C_gig_full_counts$rownames
C_gig_full_counts <- C_gig_full_counts[,-33] # remove rownames to allow for vst 

# Set equal the rownames and colnames of the coldata and count data
all(rownames(C_vir_coldata ) %in% colnames(C_vir_full_counts))  #Should return TRUE
# returns TRUE
all(colnames(C_vir_full_counts) %in% rownames(C_vir_coldata  ))  
# returns TRUE
all(rownames(C_vir_coldata ) == colnames(C_vir_full_counts)) # TRUE

# Fix the order (already in correct order)
# C_vir_full_counts <- C_vir_full_counts[,row.names(C_vir_coldata)]

# C_gig
all(rownames(C_gig_coldata ) %in% colnames(C_gig_full_counts ))  #Should return TRUE
# returns TRUE
all(colnames(C_gig_full_counts ) %in% rownames(C_gig_coldata  ))  
# returns TRUE
all(rownames(C_gig_coldata ) == colnames(C_gig_full_counts )) # FALSE

# Fix the order
C_gig_full_counts  <- C_gig_full_counts [,row.names(C_gig_coldata)]
all(rownames(C_gig_coldata ) %in% colnames(C_gig_full_counts ))  #Should return TRUE
# returns TRUE
all(colnames(C_gig_full_counts ) %in% rownames(C_gig_coldata  ))  
# returns TRUE
all(rownames(C_gig_coldata ) == colnames(C_gig_full_counts )) # TRUE

# Change NA values to be zero for both 
C_vir_full_counts[is.na(C_vir_full_counts)] <- 0
C_gig_full_counts[is.na(C_gig_full_counts)] <- 0

# Make DEseq data set from matrix so that the coldata gets attached
C_vir_full_counts_dds <- DESeqDataSetFromMatrix(countData = C_vir_full_counts,
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


#### SESSION INFO FOR RUNNING SCRIPTS FEB 2020 ####

sessionInfo()
# R version 3.6.1 (2019-07-05)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Mojave 10.14
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] forcats_0.4.0               readr_1.3.1                 tidyverse_1.2.1             limma_3.40.6                ggpubr_0.2.4               
# [6] ggfortify_0.4.8             tibble_2.1.3                purrr_0.3.3                 Repitools_1.30.0            plyr_1.8.5                 
# [11] reshape2_1.4.3              UpSetR_1.4.0                rtracklayer_1.44.3          stringr_1.4.0               tidyr_1.0.0                
# [16] fission_1.4.0               genefilter_1.66.0           apeglm_1.6.0                questionr_0.7.0             RColorBrewer_1.1-2         
# [21] pheatmap_1.0.12             dplyr_0.8.3                 magrittr_1.5                ggplot2_3.2.1               DESeq2_1.24.0              
# [26] SummarizedExperiment_1.14.1 DelayedArray_0.10.0         BiocParallel_1.18.1         matrixStats_0.54.0          Biobase_2.44.0             
# [31] GenomicRanges_1.36.0        GenomeInfoDb_1.20.0         IRanges_2.18.2              S4Vectors_0.22.0            BiocGenerics_0.30.0        
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1             backports_1.1.5          Hmisc_4.2-0              R.devices_2.16.0         aroma.light_3.14.0      
# [6] R.rsp_0.43.1             lazyeval_0.2.2           splines_3.6.1            listenv_0.7.0            digest_0.6.23           
# [11] htmltools_0.3.6          fansi_0.4.1              gdata_2.18.0             Rsolnp_1.16              checkmate_1.9.4         
# [16] memoise_1.1.0            BSgenome_1.52.0          aroma.apd_0.6.0          cluster_2.1.0            globals_0.12.4          
# [21] Biostrings_2.52.0        annotate_1.62.0          modelr_0.1.5             R.utils_2.9.0            aroma.core_3.2.0        
# [26] colorspace_1.4-1         rvest_0.3.4              blob_1.2.0               haven_2.1.1              xfun_0.9                
# [31] jsonlite_1.6             crayon_1.3.4             RCurl_1.95-4.12          zeallot_0.1.0            survival_2.44-1.1       
# [36] glue_1.3.1               R.huge_0.9.0             gtable_0.3.0             zlibbioc_1.30.0          XVector_0.24.0          
# [41] R.cache_0.13.0           scales_1.1.0             vsn_3.52.0               DBI_1.0.0                edgeR_3.26.8            
# [46] miniUI_0.1.1.1           Rcpp_1.0.3               xtable_1.8-4             emdbook_1.3.11           htmlTable_1.13.1        
# [51] foreign_0.8-72           bit_1.1-14               preprocessCore_1.46.0    Formula_1.2-3            truncnorm_1.0-8         
# [56] httr_1.4.1               htmlwidgets_1.3          gplots_3.0.1.1           acepack_1.4.1            farver_2.0.3            
# [61] pkgconfig_2.0.3          XML_3.98-1.20            R.methodsS3_1.7.1        nnet_7.3-12              locfit_1.5-9.1          
# [66] DNAcopy_1.58.0           labeling_0.3             tidyselect_0.2.5         rlang_0.4.2              later_0.8.0             
# [71] AnnotationDbi_1.46.1     cellranger_1.1.0         munsell_0.5.0            tools_3.6.1              cli_2.0.1               
# [76] generics_0.0.2           RSQLite_2.1.2            broom_0.5.2              yaml_2.2.0               knitr_1.24              
# [81] bit64_0.9-7              caTools_1.17.1.2         nlme_3.1-141             future_1.14.0            mime_0.7                
# [86] R.oo_1.22.0              xml2_1.2.2               compiler_3.6.1           rstudioapi_0.10          ggsignif_0.6.0          
# [91] affyio_1.54.0            geneplotter_1.62.0       stringi_1.4.5            highr_0.8                gsmoothr_0.1.7          
# [96] lattice_0.20-38          Matrix_1.2-17            vctrs_0.2.1              pillar_1.4.3             lifecycle_0.1.0         
# [101] BiocManager_1.30.4       data.table_1.12.8        bitops_1.0-6             httpuv_1.5.1             R.filesets_2.13.0       
# [106] R6_2.4.1                 latticeExtra_0.6-28      affy_1.62.0              promises_1.0.1           KernSmooth_2.23-15      
# [111] gridExtra_2.3            aroma.affymetrix_3.2.0   codetools_0.2-16         Ringo_1.48.0             MASS_7.3-51.4           
# [116] gtools_3.8.1             assertthat_0.2.1         withr_2.1.2              GenomicAlignments_1.20.1 Rsamtools_2.0.0         
# [121] GenomeInfoDbData_1.2.1   hms_0.5.1                grid_3.6.1               rpart_4.1-15             coda_0.19-3             
# [126] PSCBS_0.65.0             bbmle_1.0.20             lubridate_1.7.4          numDeriv_2016.8-1.1      shiny_1.3.2             
# [131] base64enc_0.1-3         
