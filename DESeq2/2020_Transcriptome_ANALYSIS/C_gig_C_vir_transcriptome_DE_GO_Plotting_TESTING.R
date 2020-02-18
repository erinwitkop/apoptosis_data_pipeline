# Script to Analyze Differential Expression of Transcriptomes using DESeq2
# Erin Roberts, PhD Candidate University of Rhode Island 
# 2/13/2020

# This script will calculate differential expression of apoptosis genes from C. gigas and C. virginica, for each experiment separately.
# Comparisons and formulas used to calculate differential expression for each experiment will be unique to that experiment, specifically
# tailored to when the infection was most acute. Challenge group samlpes will always be compared to their own control.

# This script contains all the multiple testing performed by Erin Roberts during analysis, while C_gig_C_vir_transcriptome_DE_GO_Plotting.R is the finished code
# only including code necessary to recreate the analysis and create plots used for comparison and publication 

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
library(reshape2)
library(plyr)
library(Repitools)
library(purrr)
library(tibble)
library(ggfortify)
library(ggpubr)

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
C_vir_rtracklayer_transcripts <- filter(C_vir_rtracklayer, grepl("XM", transcript_id))

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
                                                          C_vir_rtracklayer$product, ignore.case = TRUE),]

# Terms to remove
# remove complement C1q proteins, dual specificity protein phosphatase 1B-like, remove kunitz-type, and NOT other kDA protein names so I can keep all heat shock proteins
C_vir_rtracklayer_apop_product_final <- C_vir_rtracklayer_apop_product[!grepl("complement C1q", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) & 
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

#### DELORGERIL OSHV1 TRANSCRIPTOME ANALYSIS ####

#### HE OSHV1 TRANSCRIPTOME ANALYSIS ####



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

# Fix the order using code format below if necessary 
# Dermo_counts <- Dermo_counts[,colnames(Dermo_counts)]

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
                                    design = ~ group) # add time to control for injection and time effect
Zhang_dds_path <- DESeqDataSetFromMatrix(countData = Zhang_counts,
                                                 colData = Zhang_coldata,
                                                 design = ~  path) # add time to control for injection and time effect
Zhang_dds_broken_group <- DESeqDataSetFromMatrix(countData = Zhang_counts,
                                                 colData = Zhang_coldata,
                                                 design = ~time+ group_by_sim) # add time to control for injection and time effect

## Prefiltering the data
# Data prefiltering helps decrease the size of the data set and get rid of
# rows with no data or very minimal data (<10). Apply a minimal filtering here as more stringent filtering will be applied later
Zhang_dds <- Zhang_dds[ rowSums(counts(Zhang_dds)) > 10, ]
Zhang_dds_broken_group <- Zhang_dds_broken_group[ rowSums(counts(Zhang_dds_broken_group)) > 10, ]
Zhang_dds_path <- Zhang_dds_path[ rowSums(counts(Zhang_dds_path)) > 10, ]

## Check levels 
# It is prefered in R that the first level of a factor be the reference level for comparison
# (e.g. control, or untreated samples), so we can relevel the factor like so
# Check factor levels, set it so that comparison group is the first
levels(Zhang_coldata$condition) #"Control" "LPS"     "M_lut"   "PBS"     "V_aes"   "V_alg_1" "V_alg_2" "V_ang"   "V_tub"  
levels(Zhang_coldata$group)# "challenge" "control"    
Zhang_dds$group <- factor(Zhang_dds$group , levels = c("control","challenge"))
levels(Zhang_dds$group)
levels(Zhang_coldata$group_by_sim) # "control"             "LPS_M_lut"           "V_aes_V_alg1_V_alg2" "V_tub_V_ang"   
levels(Zhang_coldata$path)  # [1] "control" "nonpath" "path"    don't need to relevel
levels(Zhang_coldata$time) # "12h"          "No_injection"
Zhang_dds$time <- factor(Zhang_dds$group , levels = c("12h"))

## DATA TRANSFORMATION AND VISUALIZATION
# Assess sample clustering after setting initial formula for comparison
Zhang_dds_rlog <- rlog(Zhang_dds, blind = TRUE) # keep blind = true before deseq function has been run
Zhang_dds_broken_group_rlog <- rlog(Zhang_dds_broken_group, blind = TRUE)
Zhang_dds_path_rlog <- rlog(Zhang_dds_path, blind = TRUE)

## PCA plot visualization of individuals in the family 
plotPCA(Zhang_dds_rlog, intgroup=c("group", "condition"))
  # control and PBS still clustering, LPS M. lut and V. aes clustering
  # with MSTRG removed, the LPS and M. lut cluster most closely 

### DIFFERENTIAL EXPRESSION ANALYSIS
# run pipeline with single command because the formula has already been specified
# Steps: estimation of size factors (controlling for differences in the sequencing depth of the samples), 
# the estimation of dispersion values for each gene, 
# and fitting a generalized linear model.

Zhang_dds_deseq <- DESeq(Zhang_dds) 
Zhang_dds_broken_group_deseq <- DESeq(Zhang_dds_broken_group)
Zhang_dds_path_deseq <- DESeq(Zhang_dds_path)

## Check the resultsNames object of each to look at the available coefficients for use in lfcShrink command
resultsNames(Zhang_dds_deseq) # [1] "Intercept"                  "group_challenge_vs_control"
resultsNames(Zhang_dds_broken_group_deseq) 
#[1] "Intercept"                                   "time_No_injection_vs_12h"                    "group_by_sim_LPS_M_lut_vs_control"          
#[4] "group_by_sim_V_aes_V_alg1_V_alg2_vs_control" "group_by_sim_V_tub_V_ang_vs_control" 

resultsNames(Zhang_dds_path_deseq) 
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

Zhang_dds_path_deseq_res_path_nonpath <- results(Zhang_dds_path_deseq, alpha=0.05, name=  "path_nonpath_vs_control")
Zhang_dds_path_deseq_res_path_path <- results(Zhang_dds_path_deseq, alpha=0.05, name=  "path_path_vs_control")

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

## DECISION: USE SAME RES OBJECT TO KEEP ALPHA ADJUSTMENT, and use LFCShrink ashr

Zhang_dds_deseq_res_LFC <- lfcShrink(Zhang_dds_deseq, coef="group_challenge_vs_control", type= "apeglm", res=Zhang_dds_deseq_res)
# Review results object summary
summary(Zhang_dds_deseq_res_LFC) # SHOWS NUMBER OF SIGNIFICANT GENES
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

Zhang_dds_deseq_res_V_alg1_LFC <- lfcShrink(Zhang_dds_broken_group_deseq, coef="group_by_sim_V_aes_V_alg1_V_alg2_vs_control" , type= "apeglm", res=Zhang_dds_deseq_res_V_alg1)
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

Zhang_dds_deseq_res_V_tub_LFC <- lfcShrink(Zhang_dds_broken_group_deseq, coef="group_by_sim_V_tub_V_ang_vs_control" , type= "apeglm", res=Zhang_dds_deseq_res_V_tub)
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

Zhang_dds_deseq_res_LFC_LPS <- lfcShrink(Zhang_dds_broken_group_deseq, coef="group_by_sim_LPS_M_lut_vs_control", type= "apeglm", res=Zhang_dds_deseq_res_LPS)
summary(Zhang_dds_deseq_res_LFC_LPS)
#out of 33418 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 615, 1.8%
#LFC < 0 (down)     : 207, 0.62%
#outliers [1]       : 0, 0%
#low counts [2]     : 11662, 35%
#(mean count < 13)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

Zhang_dds_path_deseq_res_path_nonpath_LFC <- lfcShrink(Zhang_dds_path_deseq, coef="path_nonpath_vs_control", type= "apeglm", res=Zhang_dds_path_deseq_res_path_nonpath)
summary(Zhang_dds_path_deseq_res_path_nonpath_LFC)
#out of 33418 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 446, 1.3%
#LFC < 0 (down)     : 137, 0.41%
#outliers [1]       : 1598, 4.8%
#low counts [2]     : 3888, 12%
#(mean count < 3)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

Zhang_dds_path_deseq_res_path_path_LFC <- lfcShrink(Zhang_dds_path_deseq, coef="path_path_vs_control", type= "apeglm", res=Zhang_dds_path_deseq_res_path_path)
summary(Zhang_dds_path_deseq_res_path_path_LFC)
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

Zhang_dds_path_deseq_res_path_nonpath_LFC_sig <-subset(Zhang_dds_path_deseq_res_path_nonpath_LFC, padj < 0.05)
Zhang_dds_path_deseq_res_path_nonpath_LFC_sig$transcript_id <- row.names(Zhang_dds_path_deseq_res_path_nonpath_LFC_sig )
Zhang_dds_path_deseq_res_path_nonpath_LFC_sig <- as.data.frame(Zhang_dds_path_deseq_res_path_nonpath_LFC_sig)
nrow(Zhang_dds_path_deseq_res_path_nonpath_LFC_sig )  #583

Zhang_dds_path_deseq_res_path_path_LFC_sig  <-subset(Zhang_dds_path_deseq_res_path_path_LFC , padj < 0.05)
Zhang_dds_path_deseq_res_path_path_LFC_sig $transcript_id <- row.names(Zhang_dds_path_deseq_res_path_path_LFC_sig  )
Zhang_dds_path_deseq_res_path_path_LFC_sig  <- as.data.frame(Zhang_dds_path_deseq_res_path_path_LFC_sig )
nrow(Zhang_dds_path_deseq_res_path_path_LFC_sig  )  #517

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

# heatmap of all apoptosis genes 
Zhang_counts_apop_assay <-  assay(Zhang_counts_apop_rlog)[,]
Zhang_counts_apop_assay_mat <- Zhang_counts_apop_assay - rowMeans(Zhang_counts_apop_assay)
Zhang_counts_apop_assay_anno <- as.data.frame(colData(Zhang_counts_apop_rlog )[, c("condition","group_by_sim")])
Zhang_counts_apop_assay_heatmap <- pheatmap(Zhang_counts_apop_assay_mat  , annotation_col = Zhang_counts_apop_assay_anno)
head(Zhang_counts_apop_assay_mat )
# control and PBS grouping, V_aes and V_ alg2 grouping, V_ang V alg 1 and V tub clustering

# heatmap of most variable apoptosis genes (this selects genes with the greatest variance in the sample)
topVarGenes_Zhang_counts_apop_assay <-  head(order(rowVars(assay(Zhang_counts_apop_rlog)), decreasing = TRUE), 50) 
top_Var_Zhang_counts_apop_assay_mat<- assay(Zhang_counts_apop_rlog )[topVarGenes_Zhang_counts_apop_assay,]
top_Var_Zhang_counts_apop_assay_mat <- top_Var_Zhang_counts_apop_assay_mat - rowMeans(top_Var_Zhang_counts_apop_assay_mat)
top_Var_Zhang_counts_apop_assay_anno <- as.data.frame(colData(Zhang_counts_apop_rlog )[, c("condition","group_by_sim")])
top_Var_Zhang_counts_apop_assay_heatmap <- pheatmap(top_Var_Zhang_counts_apop_assay_mat  , annotation_col = top_Var_Zhang_counts_apop_assay_anno)
head(top_Var_Zhang_counts_apop_assay_mat )
# same grouping as above with top 200, top 100 changes grouping M lut LPS V alg1 and V alg2 group together in a cluster, while V tub, V ang and V aes group with control

# annotate the top 100 genes 
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

Zhang_dds_path_deseq_res_path_nonpath_LFC_sig_APOP <- merge(Zhang_dds_path_deseq_res_path_nonpath_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")
Zhang_dds_path_deseq_res_path_nonpath_LFC_sig_APOP_arranged <- arrange(Zhang_dds_path_deseq_res_path_nonpath_LFC_sig_APOP , -log2FoldChange) 
nrow(Zhang_dds_path_deseq_res_path_nonpath_LFC_sig_APOP ) # 9

Zhang_dds_path_deseq_res_path_path_LFC_sig_APOP <- merge(Zhang_dds_path_deseq_res_path_path_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")
Zhang_dds_path_deseq_res_path_path_LFC_sig_APOP_arranged <- arrange(Zhang_dds_path_deseq_res_path_path_LFC_sig_APOP , -log2FoldChange) 
nrow(Zhang_dds_path_deseq_res_path_path_LFC_sig_APOP ) # 6

# Compare apoptosis genes between group_by_sim groups
Zhang_dds_deseq_res_V_alg1_APOP_short <- Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP[,c(1:6,28)]
Zhang_dds_deseq_res_V_alg1_APOP_short$group_by_sim <- "V_aes_V_alg1_V_alg2"
Zhang_dds_deseq_res_V_tub_APOP_short <- Zhang_dds_deseq_res_V_tub_LFC_sig_APOP[,c(1:6,28)] 
Zhang_dds_deseq_res_V_tub_APOP_short$group_by_sim <- "V_tub_V_ang"
Zhang_dds_deseq_res_LPS_APOP_short <- Zhang_dds_deseq_res_LFC_LPS_sig_APOP[,c(1:6,28)] 
Zhang_dds_deseq_res_LPS_APOP_short$group_by_sim <- "LPS_M_lut"

Zhang_dds_deseq_res_V_alg1_APOP_short_plot <- ggplot(Zhang_dds_deseq_res_V_alg1_APOP_short, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Zhang Non-pathogenic Vibrio vs Control") +
  ylab("Log2 Fold Change")
Zhang_dds_deseq_res_V_tub_APOP_short_plot <- ggplot(Zhang_dds_deseq_res_V_tub_APOP_short, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Zhang Pathogenic Vibrio vs Control") +
  ylab("Log2 Fold Change")
Zhang_dds_deseq_res_LPS_APOP_short_plot <- ggplot(Zhang_dds_deseq_res_LPS_APOP_short, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Zhang LPS and M Lut LFC vs Control") +
  ylab("Log2 Fold Change")

ggarrange(Zhang_dds_deseq_res_V_alg1_APOP_short_plot , Zhang_dds_deseq_res_V_tub_APOP_short_plot,Zhang_dds_deseq_res_LPS_APOP_short_plot)

# still working on this 
Zhang_apop_group_by_sim_comparison <- rbind(Zhang_dds_deseq_res_V_alg1_APOP_short, Zhang_dds_deseq_res_V_tub_APOP_short,Zhang_dds_deseq_res_LPS_APOP_short)
Zhang_apop_group_by_sim_comparison_LFC_plot <-  ggplot(Zhang_apop_group_by_sim_comparison, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Zhang Pathogenic Vibrio, Non-pathogenic Vibrio, LPS and M Lut LFC vs Control") +
  ylab("Log2 Fold Change")


#### RUBIO VIBRIO TRANSCRIPTOME ANALYSIS ####


 
#### PROBIOTIC TRANSCRIPTOME ANALYSIS ####

#### ROD TRANSCRIPTOME ANALYSIS #### 

###### PROESTOU DERMO TRANSCRIPTOME ANALYSIS ####
#Load in TRANSCRIPT expression data as count matrix
Dermo_counts <- read.csv("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS", 
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
all(colnames(Dermo_counts) %in% rownames(Dermo_coldata))  
# returns TRUE
all(rownames(Dermo_coldata) == colnames(Dermo_counts))    # should return TRUE
# returns FALSE
# Fix the order
Dermo_counts <- Dermo_counts[,colnames(Dermo_counts)]
colnames(Dermo_counts)
rownames(Dermo_coldata)
# recheck the order
all(rownames(Dermo_coldata) == colnames(Dermo_counts)) # FALSE


# add column that is specifically rownames
Dermo_coldata$rownames <- rownames(Dermo_coldata)
head(Dermo_coldata)

####### Check levels #######
# It is prefered in R that the first level of a factor be the reference level for comparison
# (e.g. control, or untreated samples), so we can relevel the dex factor like so

# Check factor levels, set it so that comparison group is the first
levels(Dermo_coldata$FamCode) #"DA" "LB"

# LB all timepoints
levels(Dermo_coldata_LB_all$Timepoint) # #"6h"  "36h" "7d"
Dermo_coldata_LB_all$Timepoint <- factor(Dermo_coldata_LB_all$Timepoint, levels=c("6h","36h","7d"))

# DA all timepoints
levels(Dermo_coldata_DA_all$Timepoint) # #"6h"  "36h" "7d" 
Dermo_coldata_DA_all$Timepoint <- factor(Dermo_coldata_DA_all$Timepoint, levels=c("6h","36h","7d"))

#### Build DESeqDataSetFromMatrix ####

#This object specifies the count data and metadata you will work with. The design piece is critical.
#For this, I will include the library prep date because Mary said there are significant batch effects present.
#Unlike their analysis, I will set my second design variable as the Family code as it first
# differences between families that I care about. 
# dds object with all data combined for comparing all sample clustering in exploratory analysis

# Batch effects corrected in the original formula: see this thread https://support.bioconductor.org/p/121408/
# do not correct counts using the removeBatchEffects from limma based on thread above 

Dermo_dds_family <- DESeqDataSetFromMatrix(countData = Dermo_counts,
                                           colData = Dermo_coldata,
                                           design = ~Library_Prep_Date + Treat +FamCode)
##### Prefiltering the data #####

# Data prefiltering helps decrease the size of the data set and get rid of
# rows with no data or very minimal data (<10). Apply a minimal filtering here as more stringent filtering will be applied later

Dermo_dds_family <- Dermo_dds_family[ rowSums(counts(Dermo_dds_family)) > 10, ]

#### VST Transformation of count data for visualization #### 

# VST Transforming the count data for visualization
# this is only going to be performed for the full data set and not each of the individual data sets
# rlog transformation is recommended for small datasets (n<30) while VST is recommended for larger data sets
nrow(Dermo_coldata)
# there are 53 samples, VST is recommended

Dermo_dds_family_vsd <- vst(Dermo_dds_family, blind=FALSE)

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
plotPCA(Dermo_dds_family_vsd, intgroup=c("FamCode","Treat"))
#the DA family samples are closely clustered, while the LB family is less close together,
# but still falling out in the same location 

##### Differential Expression Analysis ####

## RUN DESEQ PIPELINE WITH DESEQ COMMAND

# run pipeline with single command because the formula has already been specified
# Steps: estimation of size factors (controlling for differences in the sequencing depth of the samples), 
# the estimation of dispersion values for each gene, 
# and fitting a generalized linear model.

Dermo_dds_family_deseq <- DESeq(Dermo_dds_family) 
## Check the resultsNames object of each to look at the available coef

resultsNames(Dermo_dds_family_deseq) # "Intercept" "Library_Prep_Date_17_Dec_vs_15_Dec" "FamCode_LB_vs_DA"

## BUILD THE RESULTS OBJECT

# Examining the results object, change alpha to p <0.05, looking at object metadata
# use mcols to look at metadata for each table

# Full Family DESeq
Dermo_dds_family_res <- results(Dermo_dds_family_deseq, alpha=0.05, name="FamCode_LB_vs_DA")
Dermo_dds_family_res # comparison is log2 fold change (MLE): FamCode LB vs DA

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

## REVIEW RESULTS OBJECT SUMMARY ####
# SHOWS NUMBER OF SIGNIFICANT GENES

summary(Dermo_dds_timecourse_family_res) # NO LFC performed for this one yet as it may not end up being used
#### Exploratory Plotting of Results ####

## MA Plotting

plotMA(Dermo_dds_LB_time_res_early_late_LFC , ylim = c(-5, 5))
## Histogram of P values ##
# exclude genes with very small counts to avoid spikes and plot using the LFCshrinkage

hist(Dermo_dds_LB_time_res_early_late_LFC$padj[Dermo_dds_LB_time_res_early_late_LFC$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
#### Subsetting Significant Genes by padj < 0.05 and L2FC of greater than 1.0 or less than -1.0 #####

# again, only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
Dermo_dds_family_res_LFC_sig <- subset(Dermo_dds_family_res_LFC, padj < 0.05)
Dermo_dds_family_res_LFC_sig$transcript_id <- row.names(Dermo_dds_family_res_LFC_sig)
Dermo_dds_family_res_LFC_sig <- as.data.frame(Dermo_dds_family_res_LFC_sig)
Dermo_dds_family_res_LFC_sig <- Dermo_dds_family_res_LFC_sig %>% filter(abs(log2FoldChange) >= 1.0)
nrow(Dermo_dds_family_res_LFC_sig) #3636
nrow(Dermo_dds_family_res_LFC) # 53379

#### Which transcripts are different between timepoint comparisons? #####

## Comparison of transcripts within families
# Transcripts in DA in the 6 vs 7d comparison not in the 6 vs 36r comparison
DA_early_late_not_early_middle <- Dermo_dds_DA_time_res_early_late_LFC_sig$transcript_id[!(Dermo_dds_DA_time_res_early_late_LFC_sig$transcript_id %in% Dermo_dds_DA_time_res_early_middle_LFC_sig$transcript_id)] 
length(DA_early_late_not_early_middle) #1561 , DA_early_late has 4828
#### Upset plots of overall gene expression changes between samples ####
# helpful tutorial for doing this: http://genomespot.blogspot.com/2017/09/upset-plots-as-replacement-to-venn.html
# http://crazyhottommy.blogspot.com/2016/01/upset-plot-for-overlapping-chip-seq.html
# UpsetR vignette: https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
# first extract gene list for each set
Dermo_dds_family_res_LFC_sig_id <- Dermo_dds_family_res_LFC_sig$transcript_id
Dermo_dds_6h_res_LFC_sig_id  <-Dermo_dds_6h_res_LFC_sig$transcript_id

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

#### Extract list of significant Apoptosis Genes ####
Dermo_dds_family_res_LFC_sig_annot_APOP <- Dermo_dds_family_res_LFC_sig_annot[grepl(paste(Apoptosis_names,collapse="|"), 
                                                                                    Dermo_dds_family_res_LFC_sig_annot$product, ignore.case = TRUE),]
arrange(Dermo_dds_family_res_LFC_sig_annot_APOP, -log2FoldChange) 
nrow(Dermo_dds_family_res_LFC_sig_annot_APOP) # 78





#### SESSION INFO FOR RUNNING SCRIPTS FEB 2020 ####

sessionInfo()

# R version 3.6.1 (2019-07-05)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Mojave 10.14

# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

#locale:
  # [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
  # [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] tibble_2.1.3                purrr_0.3.3                 Repitools_1.30.0            plyr_1.8.5                 
# [5] reshape2_1.4.3              UpSetR_1.4.0                rtracklayer_1.44.3          stringr_1.4.0              
# [9] tidyr_1.0.0                 fission_1.4.0               genefilter_1.66.0           apeglm_1.6.0               
# [13] questionr_0.7.0             RColorBrewer_1.1-2          pheatmap_1.0.12             dplyr_0.8.3                
# [17] magrittr_1.5                ggplot2_3.2.1               DESeq2_1.24.0               SummarizedExperiment_1.14.1
# [21] DelayedArray_0.10.0         BiocParallel_1.18.1         matrixStats_0.54.0          Biobase_2.44.0             
# [25] GenomicRanges_1.36.0        GenomeInfoDb_1.20.0         IRanges_2.18.2              S4Vectors_0.22.0           
# [29] BiocGenerics_0.30.0        

# loaded via a namespace (and not attached):
# [1] backports_1.1.5          Hmisc_4.2-0              R.devices_2.16.0         aroma.light_3.14.0       R.rsp_0.43.1            
# [6] lazyeval_0.2.2           splines_3.6.1            listenv_0.7.0            digest_0.6.23            htmltools_0.3.6         
# [11] gdata_2.18.0             Rsolnp_1.16              checkmate_1.9.4          memoise_1.1.0            BSgenome_1.52.0         
# [16] aroma.apd_0.6.0          cluster_2.1.0            limma_3.40.6             globals_0.12.4           Biostrings_2.52.0       
# [21] annotate_1.62.0          R.utils_2.9.0            aroma.core_3.2.0         colorspace_1.4-1         blob_1.2.0              
# [26] xfun_0.9                 crayon_1.3.4             RCurl_1.95-4.12          zeallot_0.1.0            survival_2.44-1.1       
# [31] glue_1.3.1               R.huge_0.9.0             gtable_0.3.0             zlibbioc_1.30.0          XVector_0.24.0          
# [36] R.cache_0.13.0           scales_1.1.0             vsn_3.52.0               DBI_1.0.0                edgeR_3.26.8            
# [41] miniUI_0.1.1.1           Rcpp_1.0.3               xtable_1.8-4             emdbook_1.3.11           htmlTable_1.13.1        
# [46] foreign_0.8-72           bit_1.1-14               preprocessCore_1.46.0    Formula_1.2-3            truncnorm_1.0-8         
# [51] htmlwidgets_1.3          gplots_3.0.1.1           acepack_1.4.1            pkgconfig_2.0.3          XML_3.98-1.20           
# [56] R.methodsS3_1.7.1        nnet_7.3-12              locfit_1.5-9.1           DNAcopy_1.58.0           tidyselect_0.2.5        
# [61] rlang_0.4.2              later_0.8.0              AnnotationDbi_1.46.1     munsell_0.5.0            tools_3.6.1             
# [66] RSQLite_2.1.2            yaml_2.2.0               knitr_1.24               bit64_0.9-7              caTools_1.17.1.2        
# [71] future_1.14.0            mime_0.7                 R.oo_1.22.0              compiler_3.6.1           rstudioapi_0.10         
# [76] affyio_1.54.0            geneplotter_1.62.0       stringi_1.4.5            highr_0.8                gsmoothr_0.1.7          
# [81] lattice_0.20-38          Matrix_1.2-17            vctrs_0.2.1              pillar_1.4.3             lifecycle_0.1.0         
# [86] BiocManager_1.30.4       data.table_1.12.8        bitops_1.0-6             httpuv_1.5.1             R.filesets_2.13.0       
# [91] R6_2.4.1                 latticeExtra_0.6-28      affy_1.62.0              promises_1.0.1           KernSmooth_2.23-15      
# [96] gridExtra_2.3            aroma.affymetrix_3.2.0   codetools_0.2-16         Ringo_1.48.0             MASS_7.3-51.4           
# [101] gtools_3.8.1             assertthat_0.2.1         withr_2.1.2              GenomicAlignments_1.20.1 Rsamtools_2.0.0         
# [106] GenomeInfoDbData_1.2.1   grid_3.6.1               rpart_4.1-15             coda_0.19-3              PSCBS_0.65.0            
# [111] bbmle_1.0.20             numDeriv_2016.8-1.1      shiny_1.3.2              base64enc_0.1-3         

