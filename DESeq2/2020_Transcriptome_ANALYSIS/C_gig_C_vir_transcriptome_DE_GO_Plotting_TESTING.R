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

##### LOADING GENOME DATA ##### 

#### C. virginica 
# Import gff file with rtracklayer
C_vir_rtracklayer <- readGFF("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/ref_C_virginica-3.0_top_level.gff3")
C_vir_rtracklayer <- as.data.frame(C_vir_rtracklayer)

#### C. gigas
# Import gff file, using new version of genome annotation
C_gig_rtracklayer <- readGFF("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS/GCF_000297895.1_oyster_v9_genomic.gff")
C_gig_rtracklayer <- as.data.frame(C_gig_rtracklayer)

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
                          'tumor necrosis factor receptor superfamily member',
                          'endonuclease G, mitochondrial',
                          'FAS-associated death domain protein',
                          'fas apoptotic inhibitory molecule 1',
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
                          'mitogen-activated protein kinase 1',
                          'dual specificity mitogen-activated protein kinase kinase 1',
                          'mitogen-activated protein kinase kinase kinase',
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
                          'high mobility group box 1',
                          'ceramide synthase',
                          'cyclic AMP-responsive element-binding protein',
                          'cell death-inducing p53-target protein 1',
                          'TP53-binding protein 1',
                          'p53-induced death domain-containing protein 1',
                          'death domain-containing protein CRADD',
                          'p63',
                          'p73',
                          'interferon regulatory factor',
                          'stimulator of interferon genes',
                          'interleukin 17-like protein',
                          'interleukin-17',
                          'interleukin-1 receptor-associated kinase 4',
                          'interleukin-1 receptor-associated kinase 1-binding protein',
                          'TGF-beta-activated kinase 1 and MAP3K7-binding protein',
                          'sterile alpha and TIR motif-containing protein',
                          'tyrosine-protein kinase JAK2',
                          'signal transducer and activator of transcription',
                          'serine-protein kinase ATM',
                          'MAP kinase-activating death domain protein',
                          'death domain-associated protein 6',
                          'leucine-rich repeat and death domain-containing protein',
                          'serine/threonine-protein kinase/endoribonuclease IRE1',
                          'eukaryotic translation initiation factor 2-alpha kinase 3',
                          'growth arrest and DNA damage-inducible protein GADD45',
                          'calpain',
                          'tyrosine-protein phosphatase non-receptor type 13',
                          'HTRA2',
                          '3-phosphoinositide-dependent protein kinase 1',
                          'dronc',
                          'pyrin',
                          'proto-oncogene c-Rel',
                          'leukocyte elastase inhibitor',
                          'protein patched homolog 1',
                          'cyclic AMP-dependent transcription factor ATF-4',
                          'Siva')

#### Grep Apoptosis protein names in genome files ####

# C virginica 
C_vir_rtracklayer <- C_vir_rtracklayer %>% filter(type == "mRNA")
C_vir_rtracklayer_apop_product <- C_vir_rtracklayer[grepl(paste(Apoptosis_names_list,collapse="|"), 
                                                          C_vir_rtracklayer$product, ignore.case = TRUE),]

# Terms to remove
# remove complement C1q proteins, dual specificity protein phosphatase 1B-like, remove kunitz-type, and NOT other kDA protein names so I can keep all heat shock proteins
C_vir_rtracklayer_apop_product_final <- C_vir_rtracklayer_apop_product[ 
  !grepl("dual specificity protein phosphatase 1B-like", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("kunitz-type", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("complement C1q", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("mannosyl", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("S-acyl", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("activator of 90 kDa", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("zinc finger protein", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) & 
    !grepl("DDB_G0272098",C_vir_rtracklayer_apop_product$product, ignore.case = TRUE)  & 
    !grepl("WD repeat-containing protein WRAP73", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("tumor protein p63-regulated gene 1-like protein", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("mitogen-activated protein kinase kinase kinase kinase 3", C_vir_rtracklayer_apop_product$product, ignore.case = TRUE),] # get rid of MAP4K

nrow(C_vir_rtracklayer_apop_product_final) # 1245
#View(C_vir_rtracklayer_apop_product_final[26][!duplicated(C_vir_rtracklayer_apop_product_final$product),])

# remove duplicated rows
C_vir_rtracklayer_apop_product_final <- C_vir_rtracklayer_apop_product_final[!duplicated(C_vir_rtracklayer_apop_product_final$product),]
nrow(C_vir_rtracklayer_apop_product_final) # 601 unique transcripts
#View(C_vir_rtracklayer_apop_product_final)

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
                             'tumor necrosis factor receptor superfamily member',
                             'endonuclease G, mitochondrial',
                             'FAS-associated death domain protein',
                             'fas apoptotic inhibitory molecule 1',
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
                             'mitogen-activated protein kinase kinase kinase',
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
                             'p73',
                             'interferon regulatory factor',
                             'stimulator of interferon genes',
                             'interleukin 17-like protein',
                             'interleukin-17',
                             'interleukin-1 receptor-associated kinase 4',
                             'interleukin-1 receptor-associated kinase 1-binding protein',
                             'TGF-beta-activated kinase 1 and MAP3K7-binding protein',
                             'sterile alpha and TIR motif-containing protein',
                             'tyrosine-protein kinase JAK2',
                             'signal transducer and activator of transcription',
                             'serine-protein kinase ATM',
                             'MAP kinase-activating death domain protein',
                             'death domain-associated protein 6',
                             'leucine-rich repeat and death domain-containing protein',
                             'serine/threonine-protein kinase/endoribonuclease IRE1',
                             'eukaryotic translation initiation factor 2-alpha kinase 3',
                             'growth arrest and DNA damage-inducible protein GADD45',
                             'calpain',
                             'tyrosine-protein phosphatase non-receptor type 13',
                             'HTRA2',
                             '3-phosphoinositide-dependent protein kinase 1',
                             'dronc',
                             'pyrin',
                             'proto-oncogene c-Rel',
                             'leukocyte elastase inhibitor',
                             'protein patched homolog 1',
                             'cyclic AMP-dependent transcription factor ATF-4',
                             'Siva')

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
nrow(C_gig_rtracklayer_apop_product) #937

# Terms to remove
# remove complement C1q proteins, dual specificity protein phosphatase 1B-like, remove kunitz-type, and other NOT kDA protein names so I can keep heat shock proteins in both
C_gig_rtracklayer_apop_product_final <- C_gig_rtracklayer_apop_product[
  # !grepl("dual specificity protein phosphatase 1B-like", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) & not present in list
  !grepl("kunitz-type", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("mannosyl", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("S-acyl", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("complement C1q", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("cytochrome c oxidase", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("cytochrome c-type heme lyase", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("cytochrome c1", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("activator of 90 kDa", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE)  & 
    !grepl("phosphatidylinositol 3-kinase 1", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE)  & 
    !grepl("phosphatidylinositol 3-kinase 2", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE)  & 
    !grepl("DDB_G0272098",C_gig_rtracklayer_apop_product$product, ignore.case = TRUE)  & 
    !grepl("WD repeat-containing protein WRAP73", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("tumor protein p63-regulated gene 1-like protein", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("mitogen-activated protein kinase kinase kinase kinase 3", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE) &
    !grepl("mitogen-activated protein kinase kinase kinase kinase 4", C_gig_rtracklayer_apop_product$product, ignore.case = TRUE),]
nrow(C_gig_rtracklayer_apop_product_final) #834

# remove duplicated rows
C_gig_rtracklayer_apop_product_final <- C_gig_rtracklayer_apop_product_final[!duplicated(C_gig_rtracklayer_apop_product_final$product),]
nrow(C_gig_rtracklayer_apop_product_final) # 564 unique transcripts

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
         colour = "group_by_sim", 
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

## Prefiltering the data
# Data prefiltering helps decrease the size of the data set and get rid of
# rows with no data or very minimal data (<10). Apply a minimal filtering here as more stringent filtering will be applied later
Zhang_dds <- Zhang_dds[ rowSums(counts(Zhang_dds)) > 10, ]

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

## Check the resultsNames object of each to look at the available coefficients for use in lfcShrink command
resultsNames(Zhang_dds_deseq) # [1] "Intercept"                                   "time_12h_vs_No_injection"                   
#[3] "group_by_sim_LPS_M_lut_vs_control"           "group_by_sim_V_aes_V_alg1_V_alg2_vs_control"
#[5] "group_by_sim_V_tub_V_ang_vs_control" 

## BUILD THE RESULTS OBJECT
# Examining the results object, change alpha to p <0.05, looking at object metadata
# use mcols to look at metadata for each table
Zhang_dds_deseq_res_V_alg1 <- results(Zhang_dds_deseq, alpha=0.05, name = "group_by_sim_V_aes_V_alg1_V_alg2_vs_control"  )
Zhang_dds_deseq_res_V_tub <- results(Zhang_dds_deseq, alpha=0.05, name= "group_by_sim_V_tub_V_ang_vs_control" )
Zhang_dds_deseq_res_LPS <- results(Zhang_dds_deseq, alpha=0.05, name= "group_by_sim_LPS_M_lut_vs_control")

head(Zhang_dds_deseq_res_V_alg1) #  group by sim Vibrio vs control 
head(Zhang_dds_deseq_res_V_tub) # group by sim LPS M lut Vtub vs control 
head(Zhang_dds_deseq_res_LPS) # group by sim V aes v alg2 vs control 

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

Zhang_dds_deseq_res_V_alg1_LFC <- lfcShrink(Zhang_dds_deseq, coef="group_by_sim_V_aes_V_alg1_V_alg2_vs_control" , type= "apeglm", res=Zhang_dds_deseq_res_V_alg1)
# Review results object summary
summary(Zhang_dds_deseq_res_V_alg1_LFC)
#out of 33418 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 1149, 3.4%
#LFC < 0 (down)     : 173, 0.52%
#outliers [1]       : 0, 0%
#low counts [2]     : 5831, 17%
#(mean count < 5)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

Zhang_dds_deseq_res_V_tub_LFC <- lfcShrink(Zhang_dds_deseq, coef="group_by_sim_V_tub_V_ang_vs_control" , type= "apeglm", res=Zhang_dds_deseq_res_V_tub)
summary(Zhang_dds_deseq_res_V_tub_LFC)
#out of 33418 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 665, 2%
#LFC < 0 (down)     : 217, 0.65%
#outliers [1]       : 0, 0%
#low counts [2]     : 9719, 29%
#(mean count < 10)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

Zhang_dds_deseq_res_LFC_LPS <- lfcShrink(Zhang_dds_deseq, coef="group_by_sim_LPS_M_lut_vs_control", type= "apeglm", res=Zhang_dds_deseq_res_LPS)
summary(Zhang_dds_deseq_res_LFC_LPS)
#out of 33418 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 620, 1.9%
#LFC < 0 (down)     : 213, 0.64%
#outliers [1]       : 0, 0%
#low counts [2]     : 11662, 35%
#(mean count < 13)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

## EXPLORATORY PLOTTING OF RESULTS 
## MA Plotting
#plotMA(Zhang_dds_deseq_res_LFC, ylim = c(-5, 5))
## Histogram of P values 
# exclude genes with very small counts to avoid spikes and plot using the LFCshrinkage
#hist(Zhang_dds_deseq_res_LFC$padj[Zhang_dds_deseq_res_LFC$baseMean > 1], breaks = 0:20/20,
  #   col = "grey50", border = "white")

### Subsetting Significant Genes by padj < 0.05
# again, only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
Zhang_dds_deseq_res_V_alg1_LFC_sig <- subset(Zhang_dds_deseq_res_V_alg1_LFC, padj < 0.05)
Zhang_dds_deseq_res_V_alg1_LFC_sig $transcript_id <- row.names(Zhang_dds_deseq_res_V_alg1_LFC_sig )
Zhang_dds_deseq_res_V_alg1_LFC_sig  <- as.data.frame(Zhang_dds_deseq_res_V_alg1_LFC_sig )
nrow(Zhang_dds_deseq_res_V_alg1_LFC_sig )  #1322

Zhang_dds_deseq_res_V_tub_LFC_sig <- subset(Zhang_dds_deseq_res_V_tub_LFC, padj < 0.05)
Zhang_dds_deseq_res_V_tub_LFC_sig  $transcript_id <- row.names(Zhang_dds_deseq_res_V_tub_LFC_sig  )
Zhang_dds_deseq_res_V_tub_LFC_sig   <- as.data.frame(Zhang_dds_deseq_res_V_tub_LFC_sig )
nrow(Zhang_dds_deseq_res_V_tub_LFC_sig  )  #882

Zhang_dds_deseq_res_LFC_LPS_sig <- subset(Zhang_dds_deseq_res_LFC_LPS, padj < 0.05)
Zhang_dds_deseq_res_LFC_LPS_sig $transcript_id <- row.names(Zhang_dds_deseq_res_LFC_LPS_sig )
Zhang_dds_deseq_res_LFC_LPS_sig  <- as.data.frame(Zhang_dds_deseq_res_LFC_LPS_sig)
nrow(Zhang_dds_deseq_res_LFC_LPS_sig)  #833

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
family_Zhang_anno <- as.data.frame(colData(Zhang_dds_rlog)[, c("condition","group_by_sim")])
family_Zhang_heatmap <- pheatmap(family_Zhang_mat , annotation_col = family_Zhang_anno)
head(family_Zhang_mat)
# With MSTRG removed, the control PBS and control control cluster together, then the LPS, M lut and V alg1 cluster
  # and the other Vibrio challenges cluster 

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
nrow(Zhang_counts_apop) #683
head(Zhang_counts_apop)

Zhang_counts_apop_dds <- DESeqDataSetFromMatrix(countData = Zhang_counts_apop,
                                         colData = Zhang_coldata,
                                         design = ~time + group_by_sim) # add time to control for injection and time effect
# Prefiltering the data and running rlog
Zhang_counts_apop_dds <- Zhang_counts_apop_dds[ rowSums(counts(Zhang_counts_apop_dds)) > 10, ]
Zhang_counts_apop_rlog <- rlog(Zhang_counts_apop_dds, blind=TRUE)

## PCA plot of rlog transformed counts for apoptosis
plotPCA(Zhang_counts_apop_rlog , intgroup=c("group_by_sim", "condition"))

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
Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP <- merge(Zhang_dds_deseq_res_V_alg1_LFC_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")
Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP_arranged <- arrange(Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP , -log2FoldChange) 
nrow(Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP) #23

Zhang_dds_deseq_res_V_tub_LFC_sig_APOP <- merge(Zhang_dds_deseq_res_V_tub_LFC_sig , C_gig_rtracklayer_apop_product_final, by = "transcript_id")
Zhang_dds_deseq_res_V_tub_LFC_sig_APOP_arranged <- arrange(Zhang_dds_deseq_res_V_tub_LFC_sig_APOP  , -log2FoldChange) 
nrow(Zhang_dds_deseq_res_V_tub_LFC_sig_APOP) # 15

Zhang_dds_deseq_res_LFC_LPS_sig_APOP <- merge(Zhang_dds_deseq_res_LFC_LPS_sig, C_gig_rtracklayer_apop_product_final, by = "transcript_id")
Zhang_dds_deseq_res_LFC_LPS_sig_APOP_arranged <- arrange(Zhang_dds_deseq_res_LFC_LPS_sig_APOP , -log2FoldChange) 
nrow(Zhang_dds_deseq_res_LFC_LPS_sig_APOP )  #15

# Compare apoptosis genes between group_by_sim groups
Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP$group_by_sim <- "V_aes_V_alg1_V_alg2"
Zhang_dds_deseq_res_V_tub_LFC_sig_APOP$group_by_sim <- "V_tub_V_ang"
Zhang_dds_deseq_res_LFC_LPS_sig_APOP$group_by_sim <- "LPS_M_lut"

# combine data frames 
Zhang_upset_all_sig_APOP <- rbind(Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
                                  Zhang_dds_deseq_res_V_tub_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
                                  Zhang_dds_deseq_res_LFC_LPS_sig_APOP[,c("product","group_by_sim","log2FoldChange")])

Zhang_upset_all_sig_APOP_joined <- full_join( Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP[,c("product","group_by_sim")], Zhang_dds_deseq_res_V_tub_LFC_sig_APOP[,c("product","group_by_sim")], by ="product")
Zhang_upset_all_sig_APOP_joined <- full_join(Zhang_upset_all_sig_APOP_joined , Zhang_dds_deseq_res_LFC_LPS_sig_APOP[,c("product","group_by_sim")], by ="product")
Zhang_upset_all_sig_APOP_joined <- na.omit(Zhang_upset_all_sig_APOP_joined)
nrow(Zhang_upset_all_sig_APOP_joined) # 11 shared between all

  # Convert into wide format using reshape
Zhang_upset_all_sig_APOP_tally <- Zhang_upset_all_sig_APOP %>% group_by(product) %>% tally() 
Zhang_upset_all_sig_APOP_upset <- Zhang_upset_all_sig_APOP %>% group_by(product) %>% mutate(value=1) %>% spread(group_by_sim, value, fill =0 )
Zhang_upset_all_sig_APOP_upset <- as.matrix(Zhang_upset_all_sig_APOP_upset)

# Make plot
Zhang_full_LFC_plot <- ggplot(Zhang_upset_all_sig_APOP, aes(x=product,y=log2FoldChange, fill=group_by_sim )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()

Zhang_dds_deseq_res_V_alg1_APOP_plot <- ggplot(Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Zhang V. alg, V. aes vs Control") +
  ylab("Log2 Fold Change")
Zhang_dds_deseq_res_V_tub_APOP_plot <- ggplot(Zhang_dds_deseq_res_V_tub_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("V. tub V. ang vs Control") +
  ylab("Log2 Fold Change")
Zhang_dds_deseq_res_LPS_APOP_plot <- ggplot(Zhang_dds_deseq_res_LFC_LPS_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Zhang LPS and M Lut LFC vs Control") +
  ylab("Log2 Fold Change")

ggarrange(Zhang_dds_deseq_res_V_alg1_APOP_plot , Zhang_dds_deseq_res_V_tub_APOP_plot,Zhang_dds_deseq_res_LPS_APOP_plot)


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
#plotMA(Rubio_dds_deseq_J2_8_res_LFC, ylim = c(-5, 5))
#plotMA(Rubio_dds_deseq_J2_9_res_LFC, ylim = c(-5, 5))
#plotMA(Rubio_dds_deseq_LGP32_res_LFC, ylim = c(-5, 5))
#plotMA(Rubio_dds_deseq_LMG20012T_res_LFC, ylim = c(-5, 5))
### Histogram of P values 
## exclude genes with very small counts to avoid spikes and plot using the LFCshrinkage
#hist(Rubio_dds_deseq_J2_8_res_LFC$padj[Rubio_dds_deseq_J2_8_res_LFC$baseMean > 1], breaks = 0:20/20,
#     col = "grey50", border = "white")
#hist(Rubio_dds_deseq_J2_9_res_LFC$padj[Rubio_dds_deseq_J2_9_res_LFC$baseMean > 1], breaks = 0:20/20,
#     col = "grey50", border = "white")
#hist(Rubio_dds_deseq_LGP32_res_LFC$padj[Rubio_dds_deseq_LGP32_res_LFC$baseMean > 1], breaks = 0:20/20,
#     col = "grey50", border = "white")
#hist(Rubio_dds_deseq_LMG20012T_res_LFC$padj[Rubio_dds_deseq_LMG20012T_res_LFC$baseMean > 1], breaks = 0:20/20,
#     col = "grey50", border = "white")

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
nrow(Rubio_counts_apop) #683
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

nrow(Rubio_dds_deseq_J2_8_res_LFC_sig_APOP_arranged ) # 41
nrow(Rubio_dds_deseq_J2_9_res_LFC_sig_APOP_arranged  ) # 50
nrow(Rubio_dds_deseq_LGP32_res_LFC_sig_APOP_arranged  ) #46
nrow(Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP_arranged) # 45

#View(Rubio_dds_deseq_J2_8_res_LFC_sig_APOP_arranged)
#View(Rubio_dds_deseq_J2_9_res_LFC_sig_APOP_arranged)
#View(Rubio_dds_deseq_LGP32_res_LFC_sig_APOP_arranged)
#View(Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP_arranged)

# Compare apoptosis genes between group_by_sim groups
Rubio_dds_deseq_J2_8_res_LFC_sig_APOP$group_by_sim <- "J2-8 non-vir"
Rubio_dds_deseq_J2_9_res_LFC_sig_APOP$group_by_sim <- "J2-9 vir"
Rubio_dds_deseq_LGP32_res_LFC_sig_APOP$group_by_sim <- "LGP32 vir"
Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP$group_by_sim <- "LMG20012T non-vir"

# combine data frames 
Rubio_all_sig_APOP <- rbind(Rubio_dds_deseq_J2_8_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
                            Rubio_dds_deseq_J2_9_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
                            Rubio_dds_deseq_LGP32_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
                            Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")])
# Make plot
Rubio_full_LFC_plot <- ggplot(Rubio_all_sig_APOP , aes(x=product,y=log2FoldChange, fill=group_by_sim )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()


Rubio_dds_deseq_J2_8_res_LFC_sig_APOP_plot <- ggplot(Rubio_dds_deseq_J2_8_res_LFC_sig_APOP , aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("J2-8 Non-Virulent vs Control") +
  ylab("Log2 Fold Change")
Rubio_dds_deseq_J2_9_res_LFC_sig_APOP_plot <- ggplot(Rubio_dds_deseq_J2_9_res_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("J2-9 Virulent vs Control") +
  ylab("Log2 Fold Change")
Rubio_dds_deseq_LGP32_res_LFC_sig_APOP_plot <- ggplot(Rubio_dds_deseq_LGP32_res_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("LGP32 Virulent vs Control") +
  ylab("Log2 Fold Change")
Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP_plot  <- ggplot(Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP , aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("LMG20012T Non-Virulent vs Control") +
  ylab("Log2 Fold Change")

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

nrow(deLorgeril_Resistant_dds_res_6_LFC_sig_APOP) #19
nrow(deLorgeril_Resistant_dds_res_12_LFC_sig_APOP) #48
nrow(deLorgeril_Resistant_dds_res_24_LFC_sig_APOP) #66
nrow(deLorgeril_Resistant_dds_res_48_LFC_sig_APOP) #18
nrow(deLorgeril_Resistant_dds_res_60_LFC_sig_APOP) # 45
nrow(deLorgeril_Resistant_dds_res_72_LFC_sig_APOP) # 23

nrow(deLorgeril_Susceptible_dds_res_6_LFC_sig_APOP) #24
nrow(deLorgeril_Susceptible_dds_res_12_LFC_sig_APOP) #75
nrow(deLorgeril_Susceptible_dds_res_24_LFC_sig_APOP) #129
nrow(deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP) #30
nrow(deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP) # 149
nrow(deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP) # 40

#View(deLorgeril_Resistant_dds_res_48_LFC_sig_APOP  $product)
#View(deLorgeril_Resistant_dds_res_60_LFC_sig_APOP  $product)
#View(deLorgeril_Resistant_dds_res_72_LFC_sig_APOP  $product)
#View(deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP$product)
#View(deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP$product) 
#View(deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP$product)

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

deLorgeril_Resistant_dds_res_48_LFC_sig_APOP_plot <- ggplot(deLorgeril_Resistant_dds_res_48_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Resistant 48hr vs Control") +
  ylab("Log2 Fold Change")
deLorgeril_Resistant_dds_res_60_LFC_sig_APOP_plot   <- ggplot(deLorgeril_Resistant_dds_res_60_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Resistant 60hr vs Control") +
  ylab("Log2 Fold Change") # mostly upregulation
deLorgeril_Resistant_dds_res_72_LFC_sig_APOP_plot   <- ggplot(deLorgeril_Resistant_dds_res_72_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Resistant 72hr vs Control") +
  ylab("Log2 Fold Change") # downregulation of cathepsins

deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP_plot <- ggplot(deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Susceptible 48hr vs Control") +
  ylab("Log2 Fold Change")
deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP_plot <- ggplot(deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Susceptible 60hr vs Control") +
  ylab("Log2 Fold Change") # TLR downregulation, a lot of upregulation
deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP_plot <- ggplot(deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP, aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
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

### Perform LFC Shrinkage with apeglm
## NOTES 
# Before plotting we need to apply an LFC shrinkage, set the coef as the specific comparison in the ResultsNames function of the deseq object
# Issue that the specific comparisons I set in my results formulas are not available in the ResultsNames coef list.
# More detailed notes about LFC Shrinkage are in the code for the Zhang Vibrio

## DECISION: USE SAME RES OBJECT TO KEEP ALPHA ADJUSTMENT, and use LFCShrink apeglm
He_dds_res_LFC<- lfcShrink(He_dds_deseq, coef= 'Condition_OsHV1_vs_control', type='apeglm', res= He_dds_res)

## EXPLORATORY PLOTTING OF RESULTS 
## MA Plotting
plotMA(He_dds_res_LFC, ylim = c(-5, 5))

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

### HE SEPARATED FOR EACH TIMEPOINT ###
# Taking same approach to splitting up the data as I did with the Dermo samples 
He_coldata_6h <- He_coldata %>% subset(Time == "6h")
He_coldata_12h <- He_coldata %>% subset(Time == "12h")
He_coldata_24h <- He_coldata %>% subset(Time == "24h")
He_coldata_48h <- He_coldata %>% subset(Time == "48h")
He_coldata_120h <- He_coldata %>% subset(Time == "120hr")

He_counts_6h <- He_counts[, row.names(He_coldata_6h)]
He_counts_12h <- He_counts[, row.names(He_coldata_12h)]
He_counts_24h <- He_counts[, row.names(He_coldata_24h)]
He_counts_48h <- He_counts[, row.names(He_coldata_48h)]
He_counts_120h <- He_counts[, row.names(He_coldata_120h)]

levels(He_coldata_6h $Time)
levels(He_coldata_12h $Time)
levels(He_coldata_24h $Time)
levels(He_coldata_48h $Time)
levels(He_coldata_120h$Time)

He_coldata_6h $Time <- droplevels(He_coldata_6h $Time)
He_coldata_12h $Time <- droplevels(He_coldata_12h $Time)
He_coldata_24h $Time <- droplevels(He_coldata_24h $Time)
He_coldata_48h $Time <- droplevels(He_coldata_48h $Time)
He_coldata_120h$Time <- droplevels(He_coldata_120h$Time)

levels(He_coldata_6h $Time)
levels(He_coldata_12h $Time)
levels(He_coldata_24h $Time)
levels(He_coldata_48h $Time)
levels(He_coldata_120h$Time)

## Creating three here so I can compare the results
He_dds_6h <- DESeqDataSetFromMatrix(countData = He_counts_6h,
                                 colData = He_coldata_6h,
                                 design = ~ Condition) 
He_dds_12h <- DESeqDataSetFromMatrix(countData = He_counts_12h,
                                 colData = He_coldata_12h,
                                 design = ~ Condition) 
He_dds_24h <- DESeqDataSetFromMatrix(countData = He_counts_24h,
                                 colData = He_coldata_24h,
                                 design = ~ Condition) 
He_dds_48h <- DESeqDataSetFromMatrix(countData = He_counts_48h,
                                 colData = He_coldata_48h,
                                 design = ~ Condition) 
He_dds_120h <- DESeqDataSetFromMatrix(countData = He_counts_120h,
                                 colData = He_coldata_120h,
                                 design = ~ Condition) 

## Prefiltering the data
He_dds_6h <- He_dds_6h [ rowSums(counts(He_dds_6h )) > 10, ]
He_dds_12h <- He_dds_12h [ rowSums(counts(He_dds_12h )) > 10, ]
He_dds_24h <- He_dds_24h [ rowSums(counts(He_dds_24h )) > 10, ]
He_dds_48h <- He_dds_48h [ rowSums(counts(He_dds_48h )) > 10, ]
He_dds_120h <- He_dds_120h [ rowSums(counts(He_dds_120h )) > 10, ]

### DIFFERENTIAL EXPRESSION ANALYSIS
He_dds_6h_deseq  <- DESeq(He_dds_6h )
He_dds_12h_deseq <- DESeq(He_dds_12h )
He_dds_24h_deseq <- DESeq(He_dds_24h )
He_dds_48h_deseq <- DESeq(He_dds_48h )
He_dds_120h_deseq <- DESeq(He_dds_120h)

## Check the resultsNames object of each to look at the available coefficients for use in lfcShrink command
resultsNames(He_dds_6h_deseq  ) # [1] "Intercept"                  "Condition_OsHV1_vs_control"
resultsNames(He_dds_12h_deseq ) # [1] "Intercept"                  "Condition_OsHV1_vs_control"
resultsNames(He_dds_24h_deseq ) # [1] "Intercept"                  "Condition_OsHV1_vs_control"
resultsNames(He_dds_48h_deseq ) # [1] "Intercept"                  "Condition_OsHV1_vs_control"
resultsNames(He_dds_120h_deseq) # [1] "Intercept"                  "Condition_OsHV1_vs_control"

## BUILD THE RESULTS OBJECT

He_dds_6h_res  <- results(He_dds_6h_deseq, alpha=0.05, name = "Condition_OsHV1_vs_control" )
He_dds_12h_res <- results(He_dds_12h_deseq, alpha=0.05, name = "Condition_OsHV1_vs_control" )
He_dds_24h_res <- results(He_dds_24h_deseq, alpha=0.05, name = "Condition_OsHV1_vs_control" )
He_dds_48h_res <- results(He_dds_48h_deseq, alpha=0.05, name = "Condition_OsHV1_vs_control" )
He_dds_120h_res <- results(He_dds_120h_deseq, alpha=0.05, name = "Condition_OsHV1_vs_control" ) 

summary(He_dds_6h_res  )
#out of 31228 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 392, 1.3%
#LFC < 0 (down)     : 391, 1.3%
#outliers [1]       : 2333, 7.5%
#low counts [2]     : 1211, 3.9%
#(mean count < 2)

summary(He_dds_12h_res )
#out of 27478 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 666, 2.4%
#LFC < 0 (down)     : 270, 0.98%
#outliers [1]       : 186, 0.68%
#low counts [2]     : 2130, 7.8%
#(mean count < 3)

summary(He_dds_24h_res )
#out of 29221 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 1392, 4.8%
#LFC < 0 (down)     : 1503, 5.1%
#outliers [1]       : 1653, 5.7%
#low counts [2]     : 1699, 5.8%
#(mean count < 2)

summary(He_dds_48h_res )
#out of 26840 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 432, 1.6%
#LFC < 0 (down)     : 404, 1.5%
#outliers [1]       : 1495, 5.6%
#low counts [2]     : 2082, 7.8%
#(mean count < 3)

summary(He_dds_120h_res)
#out of 29656 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 506, 1.7%
#LFC < 0 (down)     : 485, 1.6%
#outliers [1]       : 1899, 6.4%
#low counts [2]     : 575, 1.9%
#(mean count < 2)

### Perform LFC Shrinkage with apeglm
He_dds_6h_res_LFC  <- lfcShrink(He_dds_6h_deseq, coef = "Condition_OsHV1_vs_control", type = "apeglm", res = He_dds_6h_res)
He_dds_12h_res_LFC <- lfcShrink(He_dds_12h_deseq, coef = "Condition_OsHV1_vs_control", type = "apeglm", res = He_dds_12h_res )
He_dds_24h_res_LFC <- lfcShrink(He_dds_24h_deseq, coef = "Condition_OsHV1_vs_control", type = "apeglm", res = He_dds_24h_res )
He_dds_48h_res_LFC <- lfcShrink(He_dds_48h_deseq, coef = "Condition_OsHV1_vs_control", type = "apeglm", res = He_dds_48h_res )
He_dds_120h_res_LFC <- lfcShrink(He_dds_120h_deseq, coef = "Condition_OsHV1_vs_control", type = "apeglm", res = He_dds_120h_res)

## EXPLORATORY PLOTTING OF RESULTS 

## Histogram of P values 
hist(He_dds_6h_res_LFC  $padj[He_dds_6h_res_LFC  $baseMean >1], breaks = 0:20/20, col = "grey50", border = "white")
hist(He_dds_12h_res_LFC $padj[He_dds_12h_res_LFC $baseMean >1], breaks = 0:20/20, col = "grey50", border = "white")    
hist(He_dds_24h_res_LFC $padj[He_dds_24h_res_LFC $baseMean >1], breaks = 0:20/20, col = "grey50", border = "white")
hist(He_dds_48h_res_LFC $padj[He_dds_48h_res_LFC $baseMean >1], breaks = 0:20/20, col = "grey50", border = "white")
hist(He_dds_120h_res_LFC$padj[He_dds_120h_res_LFC$baseMean >1], breaks = 0:20/20, col = "grey50", border = "white")

### Subsetting Significant Genes by padj < 0.05
He_dds_res_6hr_sig <- subset(He_dds_6h_res_LFC, padj < 0.05)
He_dds_res_12hr_sig <- subset(He_dds_12h_res_LFC, padj < 0.05)
He_dds_res_24hr_sig <- subset(He_dds_24h_res_LFC, padj < 0.05)
He_dds_res_48hr_sig <- subset(He_dds_48h_res_LFC, padj < 0.05)
He_dds_res_120hr_sig  <- subset(He_dds_120h_res_LFC, padj < 0.05)

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

nrow(He_dds_res_6hr_sig) #783
nrow(He_dds_res_12hr_sig) # 936
nrow(He_dds_res_24hr_sig) # 2895
nrow(He_dds_res_48hr_sig) # 836
nrow(He_dds_res_120hr_sig) # 991

### Extract list of significant Apoptosis Genes (not less than or greater than 1 LFC) using merge
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

nrow(He_dds_res_6hr_sig_arranged ) # 16
nrow(He_dds_res_12hr_sig_arranged) # 20
nrow(He_dds_res_24hr_sig_arranged) # 45
nrow(He_dds_res_48hr_sig_arranged) # 18
nrow(He_dds_res_120hr_sig_arranged) # 14

He_dds_res_6hr_sig_APOP$group_by_sim <-"He_dds_res_6hr_sig_APOP"
He_dds_res_12hr_sig_APOP$group_by_sim <-"He_dds_res_12hr_sig_APOP"
He_dds_res_24hr_sig_APOP$group_by_sim <-"He_dds_res_24hr_sig_APOP"
He_dds_res_48hr_sig_APOP$group_by_sim <-"He_dds_res_48hr_sig_APOP"
He_dds_res_120hr_sig_APOP$group_by_sim <-"He_dds_res_120hr_sig_APOP"

# combine dataframes 
He_all_sig_APOP <- rbind(He_dds_res_6hr_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
                         He_dds_res_12hr_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
                         He_dds_res_24hr_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
                         He_dds_res_48hr_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
                         He_dds_res_120hr_sig_APOP[,c("product","group_by_sim","log2FoldChange")])

# Make plot
He_full_LFC_plot <- ggplot(He_all_sig_APOP , aes(x=product,y=log2FoldChange, fill=group_by_sim )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()

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
nrow(Probiotic_dds_deseq_Challenge_res_LFC_sig_APOP) # 17

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

## Creating deseq data set from matrix, controlling for the effect of time (which is the same as effect of different larvae sources )
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

nrow(Pro_RE22_dds_deseq_res_RI_6h_LFC_sig_APOP ) # 18
nrow(Pro_RE22_dds_deseq_res_RI_24h_LFC_sig_APOP) # 35
nrow(Pro_RE22_dds_deseq_res_S4_6h_LFC_sig_APOP ) # 30
nrow(Pro_RE22_dds_deseq_res_S4_24h_LFC_sig_APOP) # 42
nrow(Pro_RE22_dds_deseq_res_RE22_LFC_sig_APOP ) # 20 

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
nrow(ROD_Susceptible_dds_res_LFC_sig_APOP) # 29

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


## Alternative testing of looking at the condition effect, controlling for time  = this gets the same results as the previous formula
Dermo_Tolerant_condition_dds <- DESeqDataSetFromMatrix(countData =Dermo_Tolerant_counts,
                                             colData = Dermo_Tolerant_coldata,
                                             design = ~Lib_prep_date + Time + Condition) 
Dermo_Susceptible_condition_dds <- DESeqDataSetFromMatrix(countData = Dermo_Susceptible_counts,
                                                colData = Dermo_Susceptible_coldata,
                                                design = ~Lib_prep_date + Time + Condition)

## Collapse technical replicates
Dermo_Tolerant_dds <- collapseReplicates(Dermo_Tolerant_dds, Dermo_Tolerant_dds$Sample_ID, Dermo_Tolerant_dds$Tech_rep)
Dermo_Susceptible_dds <- collapseReplicates(Dermo_Susceptible_dds, Dermo_Susceptible_dds$Sample_ID, Dermo_Susceptible_dds$Tech_rep)

Dermo_Tolerant_condition_dds <- collapseReplicates(Dermo_Tolerant_condition_dds, Dermo_Tolerant_condition_dds$Sample_ID, Dermo_Tolerant_condition_dds$Tech_rep)
Dermo_Susceptible_condition_dds <- collapseReplicates(Dermo_Susceptible_condition_dds, Dermo_Susceptible_condition_dds$Sample_ID, Dermo_Susceptible_condition_dds$Tech_rep)


## Prefiltering the data
# Data prefiltering helps decrease the size of the data set and get rid of
# rows with no data or very minimal data (<10). Apply a minimal filtering here as more stringent filtering will be applied later
Dermo_Tolerant_dds <- Dermo_Tolerant_dds[ rowSums(counts(Dermo_Tolerant_dds )) > 10, ]
Dermo_Susceptible_dds <- Dermo_Susceptible_dds[rowSums(counts(Dermo_Susceptible_dds )) > 10, ]

Dermo_Tolerant_condition_dds <- Dermo_Tolerant_condition_dds[ rowSums(counts(Dermo_Tolerant_condition_dds )) > 10, ]
Dermo_Susceptible_condition_dds <- Dermo_Susceptible_condition_dds[rowSums(counts(Dermo_Susceptible_condition_dds )) > 10, ]

## DATA TRANSFORMATION AND VISUALIZATION
# VST - not affected by the design formula
# Assess sample clustering after setting initial formula for comparison
Dermo_Tolerant_dds_vst <- vst(Dermo_Tolerant_dds , blind = TRUE) # keep blind = true before deseq function has been run
Dermo_Susceptible_dds_vst <- vst(Dermo_Susceptible_dds , blind = TRUE) # keep blind = true before deseq function has been run

Dermo_Tolerant_condition_dds_vst <- vst(Dermo_Tolerant_condition_dds , blind = TRUE) # keep blind = true before deseq function has been run
Dermo_Susceptible_condition_dds_vst <- vst(Dermo_Susceptible_condition_dds , blind = TRUE) # keep blind = true before deseq function has been run

## PCA plot visualization of individuals in the family 
plotPCA(Dermo_Tolerant_dds_vst, intgroup= c("Time","Condition")) # no some clustering by treatment
plotPCA(Dermo_Susceptible_dds_vst , intgroup=c("Time","Condition")) # some clustering

### DIFFERENTIAL EXPRESSION ANALYSIS
# run pipeline with single command because the formula has already been specified
# Steps: estimation of size factors (controlling for differences in the sequencing depth of the samples), 
# the estimation of dispersion values for each gene,and fitting a generalized linear model.
Dermo_Tolerant_dds_deseq <- DESeq(Dermo_Tolerant_dds) 
Dermo_Susceptible_dds_deseq <- DESeq(Dermo_Susceptible_dds)
Dermo_Tolerant_condition_dds_deseq <- DESeq(Dermo_Tolerant_condition_dds) 
Dermo_Susceptible_condition_dds_deseq <- DESeq(Dermo_Susceptible_condition_dds)

## Check the resultsNames object of each to look at the available coefficients for use in lfcShrink command
resultsNames(Dermo_Tolerant_dds_deseq) # [1] "Intercept" "Lib_prep_date_12_17_2020_vs_12_15_2020" "Condition_Injected_vs_Control"         
  # [4] "Time_7d_vs_36h"    "Time_28d_vs_36h"        
resultsNames(Dermo_Susceptible_condition_dds_deseq ) #[1] "Intercept" "Lib_prep_date_12_17_2020_vs_12_15_2020" "Condition_Injected_vs_Control"         
  # [4] "Time_7d_vs_36h"   "Time_28d_vs_36h" 

resultsNames(Dermo_Tolerant_condition_dds_deseq)    
#[1] "Intercept"                              "Lib_prep_date_12_17_2020_vs_12_15_2020" "Time_7d_vs_36h"                        
#[4] "Time_28d_vs_36h"                        "Condition_Injected_vs_Control" 
resultsNames(Dermo_Susceptible_condition_dds_deseq ) 
#[1] "Intercept"                              "Lib_prep_date_12_17_2020_vs_12_15_2020" "Time_7d_vs_36h"                        
#[4] "Time_28d_vs_36h"                        "Condition_Injected_vs_Control"         

## BUILD THE RESULTS OBJECT
# Examining the results object, change alpha to p <0.05, looking at object metadata
# use mcols to look at metadata for each table
mcols(Dermo_Tolerant_dds_deseq)
mcols(Dermo_Susceptible_dds_deseq)

Dermo_Tolerant_condition_dds_res <- results(Dermo_Tolerant_condition_dds_deseq, alpha=0.05, name="Condition_Injected_vs_Control" )
Dermo_Susceptible_condition_dds_res <- results(Dermo_Susceptible_condition_dds_deseq, alpha=0.05, name="Condition_Injected_vs_Control")

Dermo_Tolerant_dds_res <- results(Dermo_Tolerant_dds_deseq, alpha=0.05, name="Condition_Injected_vs_Control" )
Dermo_Susceptible_dds_res <- results(Dermo_Susceptible_dds_deseq, alpha=0.05, name="Condition_Injected_vs_Control")

### Perform LFC Shrinkage with apeglm
## NOTES 
# Before plotting we need to apply an LFC shrinkage, set the coef as the specific comparison in the ResultsNames function of the deseq object
# Issue that the specific comparisons I set in my results formulas are not available in the ResultsNames coef list.
# More detailed notes about LFC Shrinkage are in the code for the Zhang Vibrio

## DECISION: USE SAME RES OBJECT TO KEEP ALPHA ADJUSTMENT, and use LFCShrink apeglm
Dermo_Tolerant_dds_res_LFC <- lfcShrink(Dermo_Tolerant_dds_deseq, coef="Condition_Injected_vs_Control", type="apeglm", res= Dermo_Tolerant_dds_res)
Dermo_Susceptible_dds_res_LFC <- lfcShrink(Dermo_Susceptible_dds_deseq, coef="Condition_Injected_vs_Control", type="apeglm", res=Dermo_Susceptible_dds_res)
Dermo_Tolerant_condition_dds_res_LFC <- lfcShrink(Dermo_Tolerant_condition_dds_deseq, coef="Condition_Injected_vs_Control", type="apeglm", res= Dermo_Tolerant_condition_dds_res)
Dermo_Susceptible_condition_dds_res_LFC <- lfcShrink(Dermo_Susceptible_condition_dds_deseq, coef="Condition_Injected_vs_Control", type="apeglm", res=Dermo_Susceptible_condition_dds_res)

## EXPLORATORY PLOTTING OF RESULTS 

### Subsetting Significant Genes by padj < 0.05
# again, only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
Dermo_Tolerant_condition_dds_res_LFC_sig <-  subset(Dermo_Tolerant_condition_dds_res_LFC , padj < 0.05)
Dermo_Tolerant_condition_dds_res_LFC_sig$ID<- row.names(Dermo_Tolerant_condition_dds_res_LFC_sig)
Dermo_Tolerant_condition_dds_res_LFC_sig <- as.data.frame(Dermo_Tolerant_condition_dds_res_LFC_sig)
nrow(Dermo_Tolerant_condition_dds_res_LFC_sig) # 1200

Dermo_Susceptible_condition_dds_res_LFC_sig <-  subset(Dermo_Susceptible_condition_dds_res_LFC , padj < 0.05)
Dermo_Susceptible_condition_dds_res_LFC_sig$ID<- row.names(Dermo_Susceptible_condition_dds_res_LFC_sig)
Dermo_Susceptible_condition_dds_res_LFC_sig <- as.data.frame(Dermo_Susceptible_condition_dds_res_LFC_sig)
nrow(Dermo_Susceptible_condition_dds_res_LFC_sig)  # 659

Dermo_Tolerant_dds_res_LFC_sig <-  subset(Dermo_Tolerant_dds_res_LFC , padj < 0.05)
Dermo_Tolerant_dds_res_LFC_sig$ID<- row.names(Dermo_Tolerant_dds_res_LFC_sig)
Dermo_Tolerant_dds_res_LFC_sig <- as.data.frame(Dermo_Tolerant_dds_res_LFC_sig)
nrow(Dermo_Tolerant_dds_res_LFC_sig) #1200

Dermo_Susceptible_dds_res_LFC_sig <-  subset(Dermo_Susceptible_dds_res_LFC , padj < 0.05)
Dermo_Susceptible_dds_res_LFC_sig$ID<- row.names(Dermo_Susceptible_dds_res_LFC_sig)
Dermo_Susceptible_dds_res_LFC_sig <- as.data.frame(Dermo_Susceptible_dds_res_LFC_sig)
nrow(Dermo_Susceptible_dds_res_LFC_sig)  #659

### GENE CLUSTERING ANALYSIS HEATMAPS  
# Extract genes with the highest variance across samples for each comparison using either vst or rlog transformed data
# This heatmap rather than plotting absolute expression strength plot the amount by which each gene deviates in a specific sample from the gene’s average across all samples. 
# example codes from RNAseq workflow: https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#other-comparisons

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
Dermo_Tolerant_condition_dds_res_LFC_sig_APOP <- merge(Dermo_Tolerant_condition_dds_res_LFC_sig, C_vir_rtracklayer_apop_product_final, by = "ID")
Dermo_Tolerant_condition_dds_res_LFC_sig_APOP_arranged <- arrange(Dermo_Tolerant_condition_dds_res_LFC_sig_APOP, -log2FoldChange) 
nrow(Dermo_Tolerant_condition_dds_res_LFC_sig_APOP) # 22

Dermo_Susceptible_condition_dds_res_LFC_sig_APOP <- merge(Dermo_Susceptible_condition_dds_res_LFC_sig, C_vir_rtracklayer_apop_product_final, by = "ID")
Dermo_Susceptible_condition_dds_res_LFC_sig_APOP_arranged <- arrange(Dermo_Susceptible_condition_dds_res_LFC_sig_APOP, -log2FoldChange) 
nrow(Dermo_Susceptible_condition_dds_res_LFC_sig_APOP) # 12 

Dermo_Susceptible_dds_res_LFC_sig_APOP <- merge(Dermo_Susceptible_dds_res_LFC_sig, C_vir_rtracklayer_apop_product_final, by = "ID")
Dermo_Susceptible_dds_res_LFC_sig_APOP_arranged <- arrange(Dermo_Susceptible_dds_res_LFC_sig_APOP, -log2FoldChange) 
nrow(Dermo_Susceptible_dds_res_LFC_sig_APOP) # 12

Dermo_Tolerant_dds_res_LFC_sig_APOP <- merge(Dermo_Tolerant_dds_res_LFC_sig, C_vir_rtracklayer_apop_product_final, by = "ID")
Dermo_Tolerant_dds_res_LFC_sig_APOP_arranged <- arrange(Dermo_Tolerant_dds_res_LFC_sig_APOP, -log2FoldChange) 
nrow(Dermo_Tolerant_dds_res_LFC_sig_APOP) # 22

# Combined LFC plot
Dermo_Susceptible_condition_dds_res_LFC_sig_APOP
Dermo_Tolerant_condition_dds_dres_LFC_sig_APOP 

Dermo_Susceptible_dds_res_LFC_sig_APOP
Dermo_Tolerant_dds_res_LFC_sig_APOP 

# Compare apoptosis genes between group_by_sim groups
Dermo_Susceptible_condition_dds_res_LFC_sig_APOP$group_by_sim <- "Susceptible_Control_vs_injected"
Dermo_Tolerant_condition_dds_res_LFC_sig_APOP $group_by_sim <-  "Tolerant_Control_vs_injected"
Dermo_Susceptible_dds_res_LFC_sig_APOP$group_by_sim <- "Susceptible_Control_vs_injected"
Dermo_Tolerant_dds_res_LFC_sig_APOP $group_by_sim <- "Tolerant_Control_vs_injected"

# combine data frames 
Dermo_all_sig_APOP <- rbind( 
  Dermo_Susceptible_condition_dds_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
  Dermo_Tolerant_condition_dds_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
  Dermo_Susceptible_dds_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
  Dermo_Tolerant_dds_res_LFC_sig_APOP [,c("product","group_by_sim","log2FoldChange")])

# Make plot or up and downregulated
Dermo_all_sig_APOP_downregulated <- Dermo_all_sig_APOP%>% filter(log2FoldChange <= 0)
Dermo_all_sig_APOP_upregulated <-   Dermo_all_sig_APOP%>% filter(log2FoldChange > 0)

Dermo_all_sig_APOP_downregulated_plot <- ggplot(Dermo_all_sig_APOP_downregulated , aes(x=product,y=log2FoldChange, fill=group_by_sim )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()
Dermo_all_sig_APOP_upregulated_plot <- ggplot(Dermo_all_sig_APOP_upregulated, aes(x=product,y=log2FoldChange, fill=group_by_sim )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()

Dermo_Susceptible_condition_dds_res_LFC_sig_APOP_plot <- ggplot(Dermo_Susceptible_condition_dds_res_LFC_sig_APOP , aes(x=product, y = log2FoldChange, fill=log2FoldChange)) + geom_col(position="dodge") +
  coord_flip() + scale_fill_gradient2(low="purple",mid = "grey", high="darkgreen") + ggtitle("Dermo Susceptible Control vs. Injected") +
  ylab("Log2 Fold Change")

### DERMO ANALYSIS AFTER SUBSETTING CONTROL AND TREATED FOR EACH TIMEPOINT USE THESE DATA SETS AND NOT THOSE ABOVE###
# going to split the analysis into DA and LB and split into each timepoint
Dermo_Susceptible_36hr_coldata <- Dermo_coldata %>%  subset(Family == "DA" & Time == "36h") 
Dermo_Susceptible_7d_coldata <- Dermo_coldata %>%  subset(Family == "DA" & Time == "7d") 
Dermo_Susceptible_28d_coldata <- Dermo_coldata %>%  subset(Family == "DA" & Time == "28d") 

Dermo_Susceptible_36hr_counts <- Dermo_counts[,row.names(Dermo_Susceptible_36hr_coldata )]
Dermo_Susceptible_7d_counts <- Dermo_counts[,row.names(Dermo_Susceptible_7d_coldata )]
Dermo_Susceptible_28d_counts <- Dermo_counts[,row.names(Dermo_Susceptible_28d_coldata )]

Dermo_Tolerant_36hr_coldata <- Dermo_coldata %>%  subset(Family == "LB" & Time == "36h") # subset() keeps rownames
Dermo_Tolerant_7d_coldata <- Dermo_coldata %>%  subset(Family == "LB" & Time == "7d") # subset() keeps rownames
Dermo_Tolerant_28d_coldata <- Dermo_coldata %>%  subset(Family == "LB" & Time == "28d") # subset() keeps rownames

Dermo_Tolerant_36hr_counts <- Dermo_counts[,row.names(Dermo_Tolerant_36hr_coldata)]
Dermo_Tolerant_7d_counts <- Dermo_counts[,row.names(Dermo_Tolerant_7d_coldata)]
Dermo_Tolerant_28d_counts <- Dermo_counts[,row.names(Dermo_Tolerant_28d_coldata)]

colnames(Dermo_Tolerant_36hr_counts)
colnames(Dermo_Tolerant_7d_counts )
colnames(Dermo_Tolerant_28d_counts )

## Check levels 
levels(Dermo_Tolerant_36hr_coldata$Condition)
levels(Dermo_Tolerant_36hr_coldata$Time)

## Creating two data set from matrix, one for each family  
Dermo_Tolerant_36hr_dds <- DESeqDataSetFromMatrix(countData =Dermo_Tolerant_36hr_counts,
                                             colData = Dermo_Tolerant_36hr_coldata,
                                             design = ~Lib_prep_date + Condition) 
Dermo_Tolerant_7d_dds <- DESeqDataSetFromMatrix(countData =Dermo_Tolerant_7d_counts,
                                                  colData = Dermo_Tolerant_7d_coldata,
                                                  design = ~Lib_prep_date + Condition) 
Dermo_Tolerant_28d_dds <- DESeqDataSetFromMatrix(countData =Dermo_Tolerant_28d_counts,
                                                colData = Dermo_Tolerant_28d_coldata,
                                                design = ~Lib_prep_date + Condition) 

Dermo_Susceptible_36hr_dds <- DESeqDataSetFromMatrix(countData =Dermo_Susceptible_36hr_counts,
                                                  colData = Dermo_Susceptible_36hr_coldata,
                                                  design = ~Lib_prep_date + Condition) 
Dermo_Susceptible_7d_dds <- DESeqDataSetFromMatrix(countData =Dermo_Susceptible_7d_counts,
                                                colData = Dermo_Susceptible_7d_coldata,
                                                design = ~Lib_prep_date + Condition) 
Dermo_Susceptible_28d_dds <- DESeqDataSetFromMatrix(countData =Dermo_Susceptible_28d_counts,
                                                 colData = Dermo_Susceptible_28d_coldata,
                                                 design = ~Lib_prep_date + Condition) 

## Collapse technical replicates
Dermo_Tolerant_36hr_dds <- collapseReplicates(Dermo_Tolerant_36hr_dds, Dermo_Tolerant_36hr_dds$Sample_ID, Dermo_Tolerant_36hr_dds$Tech_rep)
Dermo_Tolerant_7d_dds <- collapseReplicates(Dermo_Tolerant_7d_dds, Dermo_Tolerant_7d_dds$Sample_ID, Dermo_Tolerant_7d_dds$Tech_rep)
Dermo_Tolerant_28d_dds <- collapseReplicates(Dermo_Tolerant_28d_dds, Dermo_Tolerant_28d_dds$Sample_ID, Dermo_Tolerant_28d_dds$Tech_rep)

Dermo_Susceptible_36hr_dds <- collapseReplicates(Dermo_Susceptible_36hr_dds, Dermo_Susceptible_36hr_dds$Sample_ID, Dermo_Susceptible_36hr_dds$Tech_rep)
Dermo_Susceptible_7d_dds <- collapseReplicates(Dermo_Susceptible_7d_dds, Dermo_Susceptible_7d_dds$Sample_ID, Dermo_Susceptible_7d_dds$Tech_rep)
Dermo_Susceptible_28d_dds <- collapseReplicates(Dermo_Susceptible_28d_dds, Dermo_Susceptible_28d_dds$Sample_ID, Dermo_Susceptible_28d_dds$Tech_rep)

## Prefiltering the data
Dermo_Tolerant_36hr_dds <- Dermo_Tolerant_36hr_dds [rowSums(counts(Dermo_Tolerant_36hr_dds ))>10,]
Dermo_Tolerant_7d_dds<- Dermo_Tolerant_7d_dds[rowSums(counts(Dermo_Tolerant_7d_dds))>10,]
Dermo_Tolerant_28d_dds <- Dermo_Tolerant_28d_dds [rowSums(counts(Dermo_Tolerant_28d_dds ))>10,]
Dermo_Susceptible_36hr_dds<- Dermo_Susceptible_36hr_dds[rowSums(counts(Dermo_Susceptible_36hr_dds))>10,]
Dermo_Susceptible_7d_dds <- Dermo_Susceptible_7d_dds [rowSums(counts(Dermo_Susceptible_7d_dds ))>10,]
Dermo_Susceptible_28d_dds <- Dermo_Susceptible_28d_dds [rowSums(counts(Dermo_Susceptible_28d_dds ))>10,]

### DIFFERENTIAL EXPRESSION ANALYSIS
Dermo_Tolerant_36hr_dds_deseq <- DESeq(Dermo_Tolerant_36hr_dds )
Dermo_Tolerant_7d_dds_deseq <- DESeq(Dermo_Tolerant_7d_dds)
Dermo_Tolerant_28d_dds_deseq <- DESeq(Dermo_Tolerant_28d_dds )
Dermo_Susceptible_36hr_dds_deseq <- DESeq(Dermo_Susceptible_36hr_dds)
Dermo_Susceptible_7d_dds_deseq <- DESeq(Dermo_Susceptible_7d_dds )
Dermo_Susceptible_28d_dds_deseq <- DESeq(Dermo_Susceptible_28d_dds )

## Check the resultsNames object of each to look at the available coefficients for use in lfcShrink command
resultsNames(Dermo_Tolerant_36hr_dds_deseq)
resultsNames(Dermo_Tolerant_7d_dds_deseq)
resultsNames(Dermo_Tolerant_28d_dds_deseq)
resultsNames(Dermo_Susceptible_36hr_dds_deseq)
resultsNames(Dermo_Susceptible_7d_dds_deseq)
resultsNames(Dermo_Susceptible_28d_dds_deseq)      

## BUILD THE RESULTS OBJECT
Dermo_Tolerant_36hr_dds_res    <- results(Dermo_Tolerant_36hr_dds_deseq, alpha=0.05, name = "Condition_Injected_vs_Control")
Dermo_Tolerant_7d_dds_res      <- results(Dermo_Tolerant_7d_dds_deseq, alpha=0.05, name = "Condition_Injected_vs_Control")
Dermo_Tolerant_28d_dds_res     <- results(Dermo_Tolerant_28d_dds_deseq, alpha=0.05, name = "Condition_Injected_vs_Control")
Dermo_Susceptible_36hr_dds_res <- results(Dermo_Susceptible_36hr_dds_deseq, alpha=0.05, name = "Condition_Injected_vs_Control")
Dermo_Susceptible_7d_dds_res   <- results(Dermo_Susceptible_7d_dds_deseq, alpha=0.05, name = "Condition_Injected_vs_Control")
Dermo_Susceptible_28d_dds_res  <- results(Dermo_Susceptible_28d_dds_deseq, alpha=0.05, name = "Condition_Injected_vs_Control")

##  LFCShrink apeglm
Dermo_Tolerant_36hr_dds_res_LFC <- lfcShrink(Dermo_Tolerant_36hr_dds_deseq      , coef="Condition_Injected_vs_Control", type="apeglm", res=Dermo_Tolerant_36hr_dds_res    )
Dermo_Tolerant_7d_dds_res_LFC <- lfcShrink(Dermo_Tolerant_7d_dds_deseq          , coef="Condition_Injected_vs_Control", type="apeglm", res=Dermo_Tolerant_7d_dds_res      )
Dermo_Tolerant_28d_dds_res_LFC <- lfcShrink(Dermo_Tolerant_28d_dds_deseq        , coef="Condition_Injected_vs_Control", type="apeglm", res=Dermo_Tolerant_28d_dds_res     )
Dermo_Susceptible_36hr_dds_res_LFC <-lfcShrink(Dermo_Susceptible_36hr_dds_deseq , coef="Condition_Injected_vs_Control", type="apeglm", res=Dermo_Susceptible_36hr_dds_res )
Dermo_Susceptible_7d_dds_res_LFC <- lfcShrink(Dermo_Susceptible_7d_dds_deseq    , coef="Condition_Injected_vs_Control", type="apeglm", res=Dermo_Susceptible_7d_dds_res   )
Dermo_Susceptible_28d_dds_res_LFC <- lfcShrink(Dermo_Susceptible_28d_dds_deseq  , coef="Condition_Injected_vs_Control", type="apeglm", res=Dermo_Susceptible_28d_dds_res  )

### Subsetting Significant Genes by padj < 0.05
# again, only working with the LFCshrinkage adjusted log fold changes, and with the BH adjusted p-value
# first make sure to make the rownames with the transcript ID as a new column, then make it a dataframe for filtering
Dermo_Tolerant_36hr_dds_res_LFC_sig <-  subset(Dermo_Tolerant_36hr_dds_res_LFC , padj < 0.05)
Dermo_Tolerant_36hr_dds_res_LFC_sig$ID<- row.names(Dermo_Tolerant_36hr_dds_res_LFC_sig)
Dermo_Tolerant_36hr_dds_res_LFC_sig <- as.data.frame(Dermo_Tolerant_36hr_dds_res_LFC_sig)
nrow(Dermo_Tolerant_36hr_dds_res_LFC_sig) # 747

Dermo_Tolerant_7d_dds_res_LFC_sig <-  subset(Dermo_Tolerant_7d_dds_res_LFC , padj < 0.05)
Dermo_Tolerant_7d_dds_res_LFC_sig$ID<- row.names(Dermo_Tolerant_7d_dds_res_LFC_sig)
Dermo_Tolerant_7d_dds_res_LFC_sig <- as.data.frame(Dermo_Tolerant_7d_dds_res_LFC_sig)
nrow(Dermo_Tolerant_7d_dds_res_LFC_sig) # 996

Dermo_Tolerant_28d_dds_res_LFC_sig <-  subset(Dermo_Tolerant_28d_dds_res_LFC , padj < 0.05)
Dermo_Tolerant_28d_dds_res_LFC_sig$ID<- row.names(Dermo_Tolerant_28d_dds_res_LFC_sig)
Dermo_Tolerant_28d_dds_res_LFC_sig <- as.data.frame(Dermo_Tolerant_28d_dds_res_LFC_sig)
nrow(Dermo_Tolerant_28d_dds_res_LFC_sig) # 437

Dermo_Susceptible_36hr_dds_res_LFC_sig <-  subset(Dermo_Susceptible_36hr_dds_res_LFC , padj < 0.05)
Dermo_Susceptible_36hr_dds_res_LFC_sig$ID<- row.names(Dermo_Susceptible_36hr_dds_res_LFC_sig)
Dermo_Susceptible_36hr_dds_res_LFC_sig <- as.data.frame(Dermo_Susceptible_36hr_dds_res_LFC_sig)
nrow(Dermo_Susceptible_36hr_dds_res_LFC_sig) # 617

Dermo_Susceptible_7d_dds_res_LFC_sig <-  subset(Dermo_Susceptible_7d_dds_res_LFC , padj < 0.05)
Dermo_Susceptible_7d_dds_res_LFC_sig$ID<- row.names(Dermo_Susceptible_7d_dds_res_LFC_sig)
Dermo_Susceptible_7d_dds_res_LFC_sig <- as.data.frame(Dermo_Susceptible_7d_dds_res_LFC_sig)
nrow(Dermo_Susceptible_7d_dds_res_LFC_sig) # 378

Dermo_Susceptible_28d_dds_res_LFC_sig <-  subset(Dermo_Susceptible_28d_dds_res_LFC , padj < 0.05)
Dermo_Susceptible_28d_dds_res_LFC_sig$ID<- row.names(Dermo_Susceptible_28d_dds_res_LFC_sig)
Dermo_Susceptible_28d_dds_res_LFC_sig <- as.data.frame(Dermo_Susceptible_28d_dds_res_LFC_sig)
nrow(Dermo_Susceptible_28d_dds_res_LFC_sig) # 2462

### Extract list of significant Apoptosis Genes (not less than or greater than 1 LFC) using merge
Dermo_Susceptible_36hr_dds_res_LFC_sig_APOP <- merge(Dermo_Susceptible_36hr_dds_res_LFC_sig, C_vir_rtracklayer_apop_product_final, by = "ID")
Dermo_Susceptible_36hr_dds_res_LFC_sig_APOP_arranged <- arrange(Dermo_Susceptible_36hr_dds_res_LFC_sig_APOP, -log2FoldChange) 
nrow(Dermo_Susceptible_36hr_dds_res_LFC_sig_APOP) # 6

Dermo_Susceptible_7d_dds_res_LFC_sig_APOP <- merge(Dermo_Susceptible_7d_dds_res_LFC_sig, C_vir_rtracklayer_apop_product_final, by = "ID")
Dermo_Susceptible_7d_dds_res_LFC_sig_APOP_arranged <- arrange(Dermo_Susceptible_7d_dds_res_LFC_sig_APOP, -log2FoldChange) 
nrow(Dermo_Susceptible_7d_dds_res_LFC_sig_APOP) # 0

Dermo_Susceptible_28d_dds_res_LFC_sig_APOP <- merge(Dermo_Susceptible_28d_dds_res_LFC_sig, C_vir_rtracklayer_apop_product_final, by = "ID")
Dermo_Susceptible_28d_dds_res_LFC_sig_APOP_arranged <- arrange(Dermo_Susceptible_28d_dds_res_LFC_sig_APOP, -log2FoldChange) 
nrow(Dermo_Susceptible_28d_dds_res_LFC_sig_APOP) # 15

Dermo_Tolerant_36hr_dds_res_LFC_sig_APOP <- merge(Dermo_Tolerant_36hr_dds_res_LFC_sig, C_vir_rtracklayer_apop_product_final, by = "ID")
Dermo_Tolerant_36hr_dds_res_LFC_sig_APOP_arranged <- arrange(Dermo_Tolerant_36hr_dds_res_LFC_sig_APOP, -log2FoldChange) 
nrow(Dermo_Tolerant_36hr_dds_res_LFC_sig_APOP)  # 9

Dermo_Tolerant_7d_dds_res_LFC_sig_APOP <- merge(Dermo_Tolerant_7d_dds_res_LFC_sig, C_vir_rtracklayer_apop_product_final, by = "ID")
Dermo_Tolerant_7d_dds_res_LFC_sig_APOP_arranged <- arrange(Dermo_Tolerant_7d_dds_res_LFC_sig_APOP, -log2FoldChange) 
nrow(Dermo_Tolerant_7d_dds_res_LFC_sig_APOP) # 6

Dermo_Tolerant_28d_dds_res_LFC_sig_APOP <- merge(Dermo_Tolerant_28d_dds_res_LFC_sig, C_vir_rtracklayer_apop_product_final, by = "ID")
Dermo_Tolerant_28d_dds_res_LFC_sig_APOP_arranged <- arrange(Dermo_Tolerant_28d_dds_res_LFC_sig_APOP, -log2FoldChange) 
nrow(Dermo_Tolerant_28d_dds_res_LFC_sig_APOP) # 10
 
# Compare apoptosis genes between group_by_sim groups
Dermo_Susceptible_36hr_dds_res_LFC_sig_APOP$group_by_sim <- "Susceptible_36hr"
#Dermo_Susceptible_7d_dds_res_LFC_sig_APOP$group_by_sim <- "Susceptible_7d"
Dermo_Susceptible_28d_dds_res_LFC_sig_APOP$group_by_sim <- "Susceptible_28d"
Dermo_Tolerant_36hr_dds_res_LFC_sig_APOP$group_by_sim <- "Tolerant_36hr"
Dermo_Tolerant_7d_dds_res_LFC_sig_APOP $group_by_sim <- "Tolerant_7d"
Dermo_Tolerant_28d_dds_res_LFC_sig_APOP $group_by_sim <- "Tolerant_28d"

# combine data frames 
Dermo_separated_all_sig_APOP <- rbind( 
  Dermo_Susceptible_36hr_dds_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
  #Dermo_Susceptible_7d_dds_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
  Dermo_Susceptible_28d_dds_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
    Dermo_Tolerant_36hr_dds_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
    Dermo_Tolerant_7d_dds_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange")],
  Dermo_Tolerant_28d_dds_res_LFC_sig_APOP [,c("product","group_by_sim","log2FoldChange")])

# Make plot
Dermo_separated_all_sig_APOP_plot <- ggplot(Dermo_separated_all_sig_APOP, aes(x=product,y=log2FoldChange, fill=group_by_sim )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()

#### LFC PLOTS BY SPECIES ####

ROD_Susceptible_dds_res_LFC_sig_APOP$group_by_sim <- "ROD_susceptible"
ROD_Resistant_dds_res_LFC_sig_APOP$group_by_sim <- "ROD_resistant"
ROD_Susceptible_dds_res_LFC_sig_APOP$experiment <- "ROD"
ROD_Resistant_dds_res_LFC_sig_APOP$experiment <- "ROD"

Probiotic_dds_deseq_Challenge_res_LFC_sig_APOP$group_by_sim <- "Probiotic"
Probiotic_dds_deseq_Challenge_res_LFC_sig_APOP$experiment <- "Probiotic"

Dermo_Susceptible_36hr_dds_res_LFC_sig_APOP $experiment <- "Dermo"
Dermo_Susceptible_7d_dds_res_LFC_sig_APOP$experiment <- "Dermo"
Dermo_Susceptible_28d_dds_res_LFC_sig_APOP$experiment <- "Dermo"
Dermo_Tolerant_36hr_dds_res_LFC_sig_APOP$experiment <- "Dermo"
Dermo_Tolerant_7d_dds_res_LFC_sig_APOP$experiment <- "Dermo"
Dermo_Tolerant_28d_dds_res_LFC_sig_APOP $experiment <- "Dermo"
Dermo_Susceptible_36hr_dds_res_LFC_sig_APOP $group_by_sim <- "Susceptible_36hr"
Dermo_Susceptible_7d_dds_res_LFC_sig_APOP   $group_by_sim <- "Susceptible_7d"
Dermo_Susceptible_28d_dds_res_LFC_sig_APOP  $group_by_sim <- "Susceptible_28d"
Dermo_Tolerant_36hr_dds_res_LFC_sig_APOP    $group_by_sim <- "Tolerant_36hr"
Dermo_Tolerant_7d_dds_res_LFC_sig_APOP      $group_by_sim <- "Tolerant_7d"
Dermo_Tolerant_28d_dds_res_LFC_sig_APOP     $group_by_sim <- "Tolerant_28d"

Pro_RE22_dds_deseq_res_RI_6h_LFC_sig_APOP$experiment <- "Pro_RE22"
Pro_RE22_dds_deseq_res_RI_24h_LFC_sig_APOP$experiment <- "Pro_RE22"
Pro_RE22_dds_deseq_res_S4_6h_LFC_sig_APOP$experiment <- "Pro_RE22"
Pro_RE22_dds_deseq_res_S4_24h_LFC_sig_APOP$experiment <- "Pro_RE22"
Pro_RE22_dds_deseq_res_RE22_LFC_sig_APOP$experiment <- "Pro_RE22"

# combine all data , add in the gene and XM info. The XM lines don't also contain the XP
C_vir_apop_LFC <- rbind(
  ROD_Susceptible_dds_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "gene","transcript_id","padj")],
  ROD_Resistant_dds_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "gene","transcript_id","padj")],
  Probiotic_dds_deseq_Challenge_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "gene","transcript_id","padj")],
  Dermo_Susceptible_36hr_dds_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "gene","transcript_id","padj")],
  #Dermo_Susceptible_7d_dds_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "gene","transcript_id","padj")],
  Dermo_Susceptible_28d_dds_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "gene","transcript_id","padj")],
  Dermo_Tolerant_36hr_dds_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "gene","transcript_id","padj")],
  Dermo_Tolerant_7d_dds_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "gene","transcript_id","padj")],
  Dermo_Tolerant_28d_dds_res_LFC_sig_APOP [,c("product","group_by_sim","log2FoldChange","experiment", "gene","transcript_id","padj")],
  Pro_RE22_dds_deseq_res_RI_6h_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "gene","transcript_id","padj")],
  Pro_RE22_dds_deseq_res_RI_24h_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "gene","transcript_id","padj")],
  Pro_RE22_dds_deseq_res_S4_6h_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "gene","transcript_id","padj")],
  Pro_RE22_dds_deseq_res_S4_24h_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "gene","transcript_id","padj")],
  Pro_RE22_dds_deseq_res_RE22_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "gene","transcript_id","padj")])
  
View(arrange(C_vir_apop_LFC, -log2FoldChange) )

# assign to gene families
Apoptosis_names_df <- data.frame(product=c(# removing this because duplicated with bcl-2 'bcl-2-related protein A1',
  '^apoptosis-inducing factor 1',
  'nuclear apoptosis-inducing factor 1',
  'akt1',
  'RAC-alpha',
  'RAC-gamma',
  'methylthioribulose-1-phosphate dehydratase',
  '^tumor necrosis factor receptor',
  '^tumor necrosis factor ligand',
  'tumor necrosis factor alpha-induced protein',
  'lipopolysaccharide-induced tumor necrosis factor-alpha',
  'cell death regulator Aven',
  'bcl-2',
  'BAX',
  'BAG family molecular chaperone regulator',
  'BH3 interacting domain death agonist',
  'Bik-like killer protein',
  'CASP8 and FADD like apoptosis regulator',
  'transcription factor AP-1',
  'adenylate cyclase',
  'caspase-',
  'caspase activity and apoptosis inhibitor 1',
  'cell division cycle and apoptosis regulator protein 1',
  'CD151 antigen',
  'protein BTG1',
  'baculoviral IAP',
  'E3 ubiquitin-protein ligase XIAP',
  'cAMP-responsive element',
  'cytochrome c',
  'endonuclease G, mitochondrial',
  'FAS-associated death domain protein',
  'fas apoptotic inhibitory molecule 1',
  'fas cell surface death receptor',
  'GTPase IMAP family member',
  'harakiri',
  'DNA fragmentation factor',
  'interferon-induced protein 44',
  'interferon alpha-inducible protein 27',
  'NF-kappa-B',
  'inositol 1,4,5-trisphosphate receptor',
  'stress-activated protein kinase JNK',
  'induced myeloid leukemia cell differentiation protein Mcl-1',
  'mitogen-activated protein kinase kinase kinase',
  'mitogen-activated protein kinase-',
  'transcriptional regulator Myc-A',
  'myeloid differentiation primary response protein MyD88',
  'phorbol-12-myristate-13-acetate-induced protein 1',
  'transcription factor p65',
  'RELB proto-oncogene',
  'reticuloendotheliosis oncogene',
  'anti-apoptotic protein NR13',
  'nuclear mitotic apparatus protein 1',
  'dynamin-like 120 kDa protein',
  'cyclin-dependent kinase 5 activator 1',
  'cellular tumor antigen p53',
  'programmed cell death protein',
  'p53 and DNA damage-regulated protein 1',
  'phosphatidylinositol 3-kinase',
  'cAMP-dependent protein kinase',
  'protein kinase C',
  'BCL2 binding component 3',
  'cdc42 homolog',
  'ras-related',
  'rho-related',
  'ras-like',
  'receptor-interacting serine/threonine-protein kinase',
  'diablo homolog, mitochondrial',
  'toll-like receptor',
  'lymphotoxin-alpha',
  'CD40 ligand',
  'TNFRSF1A associated via death domain',
  'TNF receptor-associated factor',
  'netrin receptor',
  'neurotrophic receptor tyrosine kinase 1',
  'sonic hedgehog receptor',
  'mixed lineage kinase domain',
  'heat shock protein',
  'E3 ubiquitin-protein ligase CHIP',
  'protein phosphatase 1B',
  'aurora kinase A',
  'glutathione peroxidase 4',
  'gasdermin',
  'poly \\[ADP-ribose]\\ polymerase 1',
  'macrophage migration inhibitory factor',
  'hexokinase-1',
  'Raf-1 protooncogene serine/threonine kinase',
  'elastase, neutrophil expressed',
  'cathepsin',
  'PRKC apoptosis WT1 regulator protein',
  'apoptosis-stimulating of p53 protein',
  'apoptosis inhibitor 5',
  'apoptotic chromatin condensation inducer in the nucleus',
  "high mobility group box 1",
  "ceramide synthase",
  'inhibitor of apoptosis',
  'cyclic AMP-responsive element-binding protein',
  'cell death-inducing p53-target protein 1',
  'TP53-binding protein 1',
  'p53-induced death domain-containing protein 1',
  'death domain-containing protein CRADD',
  'p63',
  'p73',
  'interferon regulatory factor',
  'stimulator of interferon genes',
  'interleukin 17-like protein',
  'interleukin-17',
  'interleukin-1 receptor-associated kinase 4',
  'interleukin-1 receptor-associated kinase 1-binding protein',
  'TGF-beta-activated kinase 1 and MAP3K7-binding protein',
  'sterile alpha and TIR motif-containing protein',
  'tyrosine-protein kinase JAK2',
  'signal transducer and activator of transcription',
  'serine-protein kinase ATM',
  'MAP kinase-activating death domain protein',
  'death domain-associated protein 6',
  'leucine-rich repeat and death domain-containing protein',
  'serine/threonine-protein kinase/endoribonuclease IRE1',
  'eukaryotic translation initiation factor 2-alpha kinase 3',
  'growth arrest and DNA damage-inducible protein GADD45',
  'calpain',
  'tyrosine-protein phosphatase non-receptor type 13',
  'HTRA2',
  '3-phosphoinositide-dependent protein kinase 1',
  'dronc',
  'pyrin',
  'proto-oncogene c-Rel',
  'leukocyte elastase inhibitor',
  'protein patched homolog 1',
  'cyclic AMP-dependent transcription factor ATF-4',
  'dual specificity mitogen-activated protein kinase kinase 1',
  'mitogen-activated protein kinase 1',
  'mitochondrial Rho GTPase 1',
  'Siva'))

Apoptosis_names_df$rownames <- rownames(Apoptosis_names_df)
# Make new index
idx2 <- sapply(Apoptosis_names_df$product, grep, C_vir_apop_LFC$product)
idx1 <- sapply(seq_along(idx2), function(i) rep(i, length(idx2[[i]])))

# create df_merged which has duplicated product column for gene name hits, first column is the index name
df_merged <- cbind(Apoptosis_names_df[unlist(idx1),,drop=F], C_vir_apop_LFC[unlist(idx2),,drop=F])
head(df_merged)
df_merged[is.na(df_merged),]
nrow(df_merged) # 240 (none duplicated)
View(df_merged)
# change column names
colnames(df_merged)[1] <- "apoptosis_names_query" 
# subset rownames column
df_merged <- df_merged[,-2]
head(df_merged)

# Check df_merged apoptosis gene family name assignments here
setdiff(C_vir_apop_LFC$transcript_id,df_merged$transcript_id) # none missing
setdiff(df_merged$transcript_id,C_vir_apop_LFC$transcript_id)
df_merged[duplicated(df_merged$transcript_id),]

# HOW MANY DIFFERENTIALLY EXPRESSED IN FAMILIES
df_merged_family_count <- df_merged %>% group_by(apoptosis_names_query) %>% mutate(count=n())
df_merged_family_count_unique <- df_merged_family_count[!duplicated(df_merged_family_count),]
df_merged_family_count_unique <- unique(df_merged_family_count_unique[,c("apoptosis_names_query","count")])

# which is most frequent across experiments (total number in each experiment, divided by total transcript number for that family)
df_merged_family_exp <-  df_merged %>% group_by(experiment, apoptosis_names_query) %>% mutate(count_per_exp=n())
df_merged_family_exp_unique <- unique(df_merged_family_exp[,c("apoptosis_names_query","experiment", "count_per_exp")])
ggplot(df_merged_family_exp_unique, aes(x=apoptosis_names_query, y= count_per_exp, fill= experiment)) + geom_col() + coord_flip()
# looked at bars to see which were common across all experiments 

# how many genes are being used across the significant groups 
df_merged_family_exp_genes <- df_merged_family_exp %>% group_by(apoptosis_names_query, experiment, gene) %>% mutate(number_genes = n())
df_merged_family_exp_genes_unique <- unique(df_merged_family_exp_genes[,c("apoptosis_names_query","experiment", "number_genes")])
ggplot(df_merged_family_exp_genes_unique , aes(x=apoptosis_names_query, y=number_genes , fill= experiment)) + geom_col() + coord_flip()

# graph sd to see which ones are most variable between experiments
df_merged_family_count_sd <- df_merged_family_exp %>% group_by(apoptosis_names_query, experiment) %>% mutate(sd=sd(log2FoldChange))
df_merged_family_count_sd <- df_merged_family_count_sd %>% filter(!is.na(sd))
df_merged_family_count_sd_unique<- unique(df_merged_family_count_sd[,c("apoptosis_names_query","experiment", "sd")])
View(df_merged_family_count_sd_unique)
ggplot(df_merged_family_count_sd_unique, aes(x=apoptosis_names_query, y=sd , fill= experiment)) + geom_col() + coord_flip()

# Make plot for up and downregulated
C_vir_apop_APOP_downregulated <- C_vir_apop_LFC %>% filter(log2FoldChange <= 0)
C_vir_apop_APOP_upregulated <- C_vir_apop_LFC %>% filter(log2FoldChange > 0)

# Make plot
C_vir_apop_APOP_downregulated_plot <- ggplot(C_vir_apop_APOP_downregulated , aes(x=product,y=log2FoldChange, fill=experiment )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()
C_vir_apop_APOP_upregulated_plot <- ggplot(C_vir_apop_APOP_upregulated , aes(x=product,y=log2FoldChange, fill=experiment )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()
C_vir_apop_APOP_plot <- ggplot(C_vir_apop_LFC, aes(x=product,y=log2FoldChange, fill=experiment )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()


# C_gig
Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP$experiment <- "Zhang"
Zhang_dds_deseq_res_V_tub_LFC_sig_APOP$experiment <- "Zhang"
Zhang_dds_deseq_res_LFC_LPS_sig_APOP$experiment <- "Zhang"

Rubio_dds_deseq_J2_8_res_LFC_sig_APOP$experiment <- "Rubio"
Rubio_dds_deseq_J2_9_res_LFC_sig_APOP$experiment <- "Rubio"
Rubio_dds_deseq_LGP32_res_LFC_sig_APOP$experiment <- "Rubio"
Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP$experiment <- "Rubio"

He_dds_res_6hr_sig_APOP$experiment <- "He"
He_dds_res_12hr_sig_APOP$experiment <- "He"
He_dds_res_24hr_sig_APOP$experiment <- "He"
He_dds_res_48hr_sig_APOP$experiment <- "He"
He_dds_res_120hr_sig_APOP$experiment <- "He"

deLorgeril_Resistant_dds_res_6_LFC_sig_APOP  $experiment <- "deLorgeril"
deLorgeril_Resistant_dds_res_12_LFC_sig_APOP  $experiment <- "deLorgeril"
deLorgeril_Resistant_dds_res_24_LFC_sig_APOP  $experiment <- "deLorgeril"
deLorgeril_Resistant_dds_res_48_LFC_sig_APOP  $experiment <- "deLorgeril"
deLorgeril_Resistant_dds_res_60_LFC_sig_APOP  $experiment <- "deLorgeril"
deLorgeril_Resistant_dds_res_72_LFC_sig_APOP  $experiment <- "deLorgeril"
deLorgeril_Susceptible_dds_res_6_LFC_sig_APOP$experiment <- "deLorgeril"
deLorgeril_Susceptible_dds_res_12_LFC_sig_APOP$experiment <- "deLorgeril"
deLorgeril_Susceptible_dds_res_24_LFC_sig_APOP$experiment <- "deLorgeril"
deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP$experiment <- "deLorgeril"
deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP$experiment <- "deLorgeril"
deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP$experiment <- "deLorgeril"

# add in protein, gene and XM info for merging with ortholog information 
C_gig_apop_LFC <- rbind(Zhang_dds_deseq_res_V_alg1_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        Zhang_dds_deseq_res_V_tub_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        Zhang_dds_deseq_res_LFC_LPS_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        Rubio_dds_deseq_J2_8_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        Rubio_dds_deseq_J2_9_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        Rubio_dds_deseq_LGP32_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        Rubio_dds_deseq_LMG20012T_res_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        He_dds_res_6hr_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        He_dds_res_12hr_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        He_dds_res_24hr_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        He_dds_res_48hr_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        He_dds_res_120hr_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        deLorgeril_Resistant_dds_res_6_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        deLorgeril_Resistant_dds_res_12_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        deLorgeril_Resistant_dds_res_24_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        deLorgeril_Resistant_dds_res_48_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        deLorgeril_Resistant_dds_res_60_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        deLorgeril_Resistant_dds_res_72_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        deLorgeril_Susceptible_dds_res_6_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        deLorgeril_Susceptible_dds_res_12_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        deLorgeril_Susceptible_dds_res_24_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        deLorgeril_Susceptible_dds_res_48_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        deLorgeril_Susceptible_dds_res_60_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")],
                        deLorgeril_Susceptible_dds_res_72_LFC_sig_APOP[,c("product","group_by_sim","log2FoldChange","experiment", "seqnames","Name","gene","padj")])
nrow(C_gig_apop_LFC) # 1014

# assign to gene families 
Apoptosis_names_df_CG <- data.frame(product=c(     # removing this because duplicated with bcl-2 'bcl-2-related protein A1',
  '^apoptosis-inducing factor 1',
  'nuclear apoptosis-inducing factor 1',
  'akt1',
  'RAC-alpha',
  'RAC-gamma',
  'methylthioribulose-1-phosphate dehydratase',
  '^tumor necrosis factor receptor',
  '^tumor necrosis factor ligand',
  '^tumor necrosis factor$',# dollar sign indicates the end of the string
  'tumor necrosis factor alpha-induced protein',
  'lipopolysaccharide-induced tumor necrosis factor-alpha',
  'cell death regulator Aven',
  'bcl-2',
  'BAX',
  'BAG family molecular chaperone regulator',
  'BH3 interacting domain death agonist',
  'Bik-like killer protein',
  'CASP8 and FADD like apoptosis regulator',
  'transcription factor AP-1',
  'adenylate cyclase',
  'caspase-',
  'caspase activity and apoptosis inhibitor 1',
  'cell division cycle and apoptosis regulator protein 1',
  'CD151 antigen',
  'B-cell translocation gene 1',
  'baculoviral IAP',
  'E3 ubiquitin-protein ligase XIAP',
  'cAMP-responsive element',
  'cytochrome c',
  'endonuclease G, mitochondrial',
  'FAS-associated death domain protein',
  'fas-associated death domain protein', #as with Heat shock protein beta will need to merge these
  'fas apoptotic inhibitory molecule 1',
  'fas cell surface death receptor',
  'GTPase IMAP family member',
  'harakiri',
  'DNA fragmentation factor',
  'interferon-induced protein 44',
  'interferon alpha-inducible protein 27',
  'NF-kappa-B',
  'inositol 1,4,5-trisphosphate receptor',
  'stress-activated protein kinase JNK',
  'induced myeloid leukemia cell differentiation protein Mcl-1',
  'mitogen-activated protein kinase kinase kinase',
  'mitogen-activated protein kinase-',
  'transcriptional regulator Myc-A',
  'myeloid differentiation primary response protein MyD88',
  'phorbol-12-myristate-13-acetate-induced protein 1',
  'transcription factor p65',
  'RELB proto-oncogene',
  'reticuloendotheliosis oncogene',
  'anti-apoptotic protein NR13',
  'nuclear mitotic apparatus protein 1',
  'dynamin-like 120 kDa protein',
  'cyclin-dependent kinase 5 activator 1',
  'cellular tumor antigen p53',
  'programmed cell death protein',
  'p53 and DNA damage-regulated protein 1',
  'phosphatidylinositol 3-kinase',
  'cAMP-dependent protein kinase',
  'protein kinase C',
  'BCL2 binding component 3',
  'cdc42 homolog',
  'ras-related',
  'rho-related',
  'ras-like',
  'receptor-interacting serine/threonine-protein kinase',
  'diablo homolog, mitochondrial',
  'toll-like receptor',
  'lymphotoxin-alpha',
  'CD40 ligand',
  'TNFRSF1A associated via death domain',
  'TNF receptor-associated factor',
  'netrin receptor',
  'neurotrophic receptor tyrosine kinase 1',
  'sonic hedgehog receptor',
  'mixed lineage kinase domain',
  'heat shock protein',
  'Heat shock protein', #will need to recode these later to be in the same group, case sensitivity caused problem
  'E3 ubiquitin-protein ligase CHIP',
  'protein phosphatase 1B',
  'aurora kinase A',
  'glutathione peroxidase 4',
  'gasdermin',
  'poly \\[ADP-ribose]\\ polymerase 1',
  'macrophage migration inhibitory factor',
  'hexokinase-1',
  'Raf-1 protooncogene serine/threonine kinase',
  'elastase, neutrophil expressed',
  'cathepsin',
  'PRKC apoptosis WT1 regulator protein',
  'apoptosis-stimulating of p53 protein',
  'apoptosis inhibitor 5',
  'apoptotic chromatin condensation inducer in the nucleus',
  "high mobility group box 1",
  "ceramide synthase",
  'inhibitor of apoptosis',
  'cyclic AMP-responsive element-binding protein',
  'cell death-inducing p53-target protein 1',
  'TP53-binding protein 1',
  'p53-induced death domain-containing protein 1',
  'death domain-containing protein CRADD',
  'p63',
  'p73',
  'interferon regulatory factor',
  'stimulator of interferon genes',
  'interleukin 17-like protein',
  'interleukin-17',
  'interleukin-1 receptor-associated kinase 4',
  'interleukin-1 receptor-associated kinase 1-binding protein',
  'TGF-beta-activated kinase 1 and MAP3K7-binding protein',
  'sterile alpha and TIR motif-containing protein',
  'tyrosine-protein kinase JAK2',
  'signal transducer and activator of transcription',
  'serine-protein kinase ATM',
  'MAP kinase-activating death domain protein',
  'death domain-associated protein 6',
  'leucine-rich repeat and death domain-containing protein',
  'serine/threonine-protein kinase/endoribonuclease IRE1',
  'eukaryotic translation initiation factor 2-alpha kinase 3',
  'growth arrest and DNA damage-inducible protein GADD45',
  'calpain',
  'tyrosine-protein phosphatase non-receptor type 13',
  'HTRA2',
  '3-phosphoinositide-dependent protein kinase 1',
  'Dronc',
  'pyrin',
  'proto-oncogene c-Rel',
  'leukocyte elastase inhibitor',
  'protein patched homolog 1',
  'cyclic AMP-dependent transcription factor ATF-4',
  'dual specificity mitogen-activated protein kinase kinase 1',
  'mitogen-activated protein kinase 1',
  'mitochondrial Rho GTPase 1',
  'Siva'))

Apoptosis_names_df_CG$rownames <- rownames(Apoptosis_names_df_CG)

# Make new index
idx2 <- sapply(Apoptosis_names_df_CG$product, grep, C_gig_apop_LFC$product)
idx1 <- sapply(seq_along(idx2), function(i) rep(i, length(idx2[[i]])))

# create df_merged which has duplicated product column for gene name hits, first column is the index name
df_merged_CG <- cbind(Apoptosis_names_df_CG[unlist(idx1),,drop=F], C_gig_apop_LFC[unlist(idx2),,drop=F])
head(df_merged_CG)
# change column names
colnames(df_merged_CG)[1] <- "apoptosis_names_query" 
# subset rownames column
df_merged_CG <- df_merged_CG[,-2]
head(df_merged_CG)
df_merged_CG[is.na(df_merged_CG),]
nrow(df_merged_CG) # 1014 (none duplicated, missing two)
View(df_merged_CG)

# HOW MANY DIFFERENTIALLY EXPRESSED IN FAMILIES
df_merged_family_count_CG <- df_merged_CG %>% group_by(apoptosis_names_query) %>% mutate(count=n())
df_merged_family_count_unique_CG <- df_merged_family_count_CG[!duplicated(df_merged_family_count_CG),]
df_merged_family_count_unique_CG <- unique(df_merged_family_count_unique_CG[,c("apoptosis_names_query","count")])

# which is most frequent across experiments (total number in each experiment, divided by total transcript number for that family)
df_merged_family_exp_CG <-  df_merged_CG %>% group_by(experiment, apoptosis_names_query) %>% mutate(count_per_exp=n())
df_merged_family_exp_unique_CG <- unique(df_merged_family_exp_CG[,c("apoptosis_names_query","experiment", "count_per_exp")])
ggplot(df_merged_family_exp_unique_CG, aes(x=apoptosis_names_query, y= count_per_exp, fill= experiment)) + geom_col() + coord_flip()

# how many genes are being used across the significant groups 
df_merged_family_CG_exp_genes <- df_merged_family_exp_CG %>% group_by(apoptosis_names_query, group_by_sim, gene) %>% mutate(number_genes = n())
df_merged_family_CG_exp_genes_unique <- unique(df_merged_family_CG_exp_genes[,c("apoptosis_names_query","experiment", "number_genes")])
ggplot(df_merged_family_CG_exp_genes_unique , aes(x=apoptosis_names_query, y=number_genes , fill= experiment)) + geom_col() + coord_flip()

# graph sd to see which ones are most variable between experiments
df_merged_family_count_CG_sd <- df_merged_family_exp_CG %>% group_by(apoptosis_names_query, experiment) %>% mutate(sd=sd(log2FoldChange))
df_merged_family_count_CG_sd <- df_merged_family_count_CG_sd %>% filter(!is.na(sd))
df_merged_family_count_CG_sd_unique<- unique(df_merged_family_count_CG_sd[,c("apoptosis_names_query","experiment", "sd")])
View(df_merged_family_count_CG_sd_unique)

ggplot(df_merged_family_count_CG_sd_unique, aes(x=apoptosis_names_query, y=sd , fill= experiment)) + geom_col() + coord_flip()

# Make plot for up and downregulated
C_gig_apop_APOP_downregulated <- C_gig_apop_LFC %>% filter(log2FoldChange <= 0)
C_gig_apop_APOP_upregulated <- C_gig_apop_LFC %>% filter(log2FoldChange > 0)

# Make plot
C_gig_apop_APOP_downregulated_plot <- ggplot(C_gig_apop_APOP_downregulated , aes(x=product,y=log2FoldChange, fill=experiment )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()
C_gig_apop_APOP_upregulated_plot <- ggplot(C_gig_apop_APOP_upregulated , aes(x=product,y=log2FoldChange, fill=experiment )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()
C_gig_apop_APOP_plot <- ggplot(C_gig_apop_LFC  , aes(x=product,y=log2FoldChange, fill=experiment )) + geom_col(position="dodge") + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) + coord_flip()


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
all(rownames(Probiotic_counts_apop) == rownames(Dermo_Tolerant_counts_apop)) # TRUE
all(rownames(Probiotic_counts_apop) %in% rownames(Dermo_Tolerant_counts_apop)) # TRUE
Probiotic_counts_apop <- Probiotic_counts_apop[row.names(Dermo_Tolerant_counts_apop),]
all(rownames(Probiotic_counts_apop) == rownames(Dermo_Tolerant_counts_apop)) # TRUE

# Check ROD order and change if necessary 
all(rownames(ROD_Resistant_counts_apop) == rownames(Dermo_Tolerant_counts_apop)) # TRUE
all(rownames(ROD_Resistant_counts_apop) %in% rownames(Dermo_Tolerant_counts_apop)) # TRUE
ROD_Resistant_counts_apop <- ROD_Resistant_counts_apop[row.names(Dermo_Tolerant_counts_apop),]
ROD_Susceptible_counts_apop <-  ROD_Susceptible_counts_apop[row.names(Dermo_Tolerant_counts_apop),]
all(rownames(ROD_Resistant_counts_apop) == rownames(Dermo_Tolerant_counts_apop)) # TRUE
all(rownames(ROD_Susceptible_counts_apop) == rownames(Dermo_Tolerant_counts_apop)) # TRUE
all(rownames(ROD_Susceptible_counts_apop) == rownames(ROD_Resistant_counts_apop)) # TRUE

# Colbind all C. virginica tables
C_virginica_apop_counts <- cbind(Dermo_Susceptible_counts_apop,Dermo_Tolerant_counts_apop,
                                 Probiotic_counts_apop,Pro_RE22_counts_apop, ROD_Resistant_counts_apop,
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
all(rownames(Zhang_counts_apop) == rownames(Rubio_counts_apop)) # TRUE
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

# Combine transcript count matrices by rownames (change this so rowname not added down here to same dataframe so data doesnt need to be reloaded
nrow(Dermo_counts) # 67868
nrow(Probiotic_counts) # 67876
nrow(ROD_counts) # 67870
nrow(Pro_RE22_counts) # 67873

Dermo_counts_PCA <- Dermo_counts
Probiotic_counts_PCA <- Probiotic_counts
ROD_counts_PCA <- ROD_counts
Pro_RE22_counts_PCA <- Pro_RE22_counts

Dermo_counts_PCA$rownames<- row.names(Dermo_counts_PCA)
Probiotic_counts_PCA$rownames<- row.names(Probiotic_counts_PCA)
ROD_counts_PCA$rownames<- row.names(ROD_counts_PCA)
Pro_RE22_counts_PCA$rownames <- row.names(Pro_RE22_counts_PCA)
 
nrow(deLorgeril_counts) # 53701
nrow(He_counts) # 86859
nrow(Zhang_counts ) # 53705
nrow(Rubio_counts) # 53705

deLorgeril_counts_PCA <- deLorgeril_counts
He_counts_PCA <- He_counts
Zhang_counts_PCA <- Zhang_counts
Rubio_counts_PCA <- Rubio_counts

deLorgeril_counts_PCA $rownames <- row.names(deLorgeril_counts_PCA )
He_counts_PCA$rownames <- row.names(He_counts_PCA )
Zhang_counts_PCA $rownames <- row.names(Zhang_counts_PCA  )
Rubio_counts_PCA $rownames <- row.names(Rubio_counts_PCA )

# merge based on rownames (then delete rownames), starting with largest first
C_vir_full_counts <- left_join(Probiotic_counts_PCA , Pro_RE22_counts_PCA ,by ="rownames")
C_vir_full_counts <- left_join(C_vir_full_counts ,ROD_counts_PCA , by = "rownames")
C_vir_full_counts <- left_join(C_vir_full_counts ,Dermo_counts_PCA , by = "rownames")
colnames(C_vir_full_counts)
row.names(C_vir_full_counts) <- C_vir_full_counts$rownames
C_vir_common_transcripts <- C_vir_full_counts$rownames
length(C_vir_common_transcripts) #67876
 
head(C_vir_full_counts)
C_vir_full_counts <- C_vir_full_counts[,-7] # remove rownames to allow for vst 

C_gig_full_counts <- left_join(He_counts_PCA ,Zhang_counts_PCA , by ="rownames")
C_gig_full_counts <- left_join(C_gig_full_counts,Rubio_counts_PCA , by = "rownames")
C_gig_full_counts <- left_join(C_gig_full_counts,deLorgeril_counts_PCA , by = "rownames")
colnames(C_gig_full_counts)
row.names(C_gig_full_counts) <- C_gig_full_counts$rownames
C_gig_common_transcripts <- C_gig_full_counts$rownames
C_gig_full_counts <- C_gig_full_counts[,-33] # remove rownames to allow for vst 

# Set equal the rownames and colnames of the coldata and count data
all(rownames(C_vir_coldata ) %in% colnames(C_vir_full_counts))  #Should return TRUE
# returns TRUE
all(colnames(C_vir_full_counts) %in% rownames(C_vir_coldata  ))  
# returns TRUE
all(rownames(C_vir_coldata ) == colnames(C_vir_full_counts)) # TRUE

# Fix the order (already in correct order)
 C_vir_full_counts <- C_vir_full_counts[,row.names(C_vir_coldata)]

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
                                                      design = ~Experiment)
# Collapse technical replicates 
C_vir_full_counts_dds <- collapseReplicates(C_vir_full_counts_dds, C_vir_full_counts_dds$Sample, C_vir_full_counts_dds$TechRep)
colnames(C_vir_full_counts_dds)
# Calculate the vst
C_vir_full_counts_vst <- varianceStabilizingTransformation(C_vir_full_counts_dds)

# Make DEseq data set from matrix so that the coldata gets attached
C_gig_full_counts_dds <- DESeqDataSetFromMatrix(countData = C_gig_full_counts ,
                                           colData = C_gig_coldata,
                                           design = ~Experiment)
# Calculate the vst
C_gig_full_counts_vst <- varianceStabilizingTransformation(C_gig_full_counts_dds)

## Remove Batch effects from experiment for C_vir
plotPCA(C_vir_full_counts_vst, "Experiment") # grouping by experiment, ROD and probiotic cluster more closely
mat_C_vir <- assay(C_vir_full_counts_vst)
mat_C_vir <- limma::removeBatchEffect(mat_C_vir, C_vir_full_counts_vst$Experiment) # remove batch effects for experiment
mat_C_vir <- limma::removeBatchEffect(mat_C_vir, C_vir_full_counts_vst$Species)
#Coefficients not estimable: batch1 batch3 batch5 batch6 
#Warning message:
#  Partial NA coefficients for 67876 probe(s) 
colnames(C_vir_full_counts_vst) # colnames got changed from SRA ID to the sample name
assay(C_vir_full_counts_vst) <- mat_C_vir
plotPCA(C_vir_full_counts_vst, "Experiment") # Probiotic and ROD now cluster together
plotPCA(C_vir_full_counts_vst, "Condition")
plotPCA(C_vir_full_counts_vst, "Sample")
plotPCA(C_vir_full_counts_vst, "Time") # no clustering by time
plotPCA(C_vir_full_counts_vst, "Family")

plotPCA(C_gig_full_counts_vst, "Experiment") # grouping by experiment, Rubio delorgeril cluster, HE and Zhang far apart
mat_C_gig <- assay(C_gig_full_counts_vst)
mat_C_gig <- limma::removeBatchEffect(mat_C_gig, C_gig_full_counts_vst$Experiment)
mat_C_gig <- limma::removeBatchEffect(mat_C_gig, C_gig_full_counts_vst$Species)
#Coefficients not estimable: batch4 batch5 batch6 
#Warning message:
#  Partial NA coefficients for 86859 probe(s) 
assay(C_gig_full_counts_vst) <- mat_C_gig
plotPCA(C_gig_full_counts_vst, "Experiment") # He, Rubio and Zhang cluster closely 
plotPCA(C_gig_full_counts_vst, "Condition")
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

#### DESCRIPTIVE TABLES ####

### Core set of transcripts and products across each challenge in each species
### Apr 6th, assess whether the code below is correct 

head(C_vir_apop_LFC)
C_vir_apop_LFC_core <- Reduce(intersect, split(C_vir_apop_LFC$transcript_id, C_vir_apop_LFC$experiment))
  # "XM_022474766.1"
C_vir_apop_LFC_core_product <- Reduce(intersect, split(C_vir_apop_LFC$product, C_vir_apop_LFC$experiment))
  # "tyrosine-protein kinase JAK2-like, transcript variant X5"
C_gig_apop_LFC_core <- Reduce(intersect, split(C_gig_apop_LFC$Name, C_gig_apop_LFC$experiment))
  # "XM_011456852.2" "XM_020072623.1" "XM_011451291.2"
C_gig_apop_LFC_product <- Reduce(intersect, split(C_gig_apop_LFC$product, C_gig_apop_LFC$experiment))
  # [1] "tyrosine-protein phosphatase non-receptor type 13, transcript variant X4" "transcription factor AP-1, transcript variant X2"                        
  # "interferon regulatory factor 1, transcript variant X2"      

# Investigate C_vir bacterial response
C_vir_apop_LFC_bac <- C_vir_apop_LFC %>% filter(experiment == "Probiotic" | experiment == "Pro_RE22" | experiment == "ROD" )
C_vir_apop_LFC_bac_core <- Reduce(intersect, split(C_vir_apop_LFC_bac$transcript_id, C_vir_apop_LFC_bac$experiment))
  # "XM_022474766.1"
C_vir_apop_LFC_bac_core_product <- Reduce(intersect, split(C_vir_apop_LFC_bac$product, C_vir_apop_LFC_bac$experiment))
  # "tyrosine-protein kinase JAK2-like, transcript variant X5"
C_vir_apop_LFC_bac_core <- Reduce(intersect, split(C_vir_apop_LFC_bac$transcript_id, C_vir_apop_LFC_bac$experiment))
# "XM_022474766.1"

# group RE22 and all the RI and S$ into a group to see which 5 interesect, based on the upset plot
C_vir_apop_LFC_bac_group_by_sim <- C_vir_apop_LFC %>% filter(experiment == "Probiotic" | experiment == "Pro_RE22" )
C_vir_apop_LFC_bac_core_group_by_sim <- Reduce(intersect, split(C_vir_apop_LFC_bac_group_by_sim$product, C_vir_apop_LFC_bac_group_by_sim$group_by_sim))
# "tyrosine-protein kinase JAK2-like, transcript variant X5"


C_vir_apop_LFC_BAC_path <- C_vir_apop_LFC_bac %>% filter(group_by_sim == "RE22" | group_by_sim == "ROD_susceptible")
C_vir_apop_LFC_BAC_path_core <- Reduce(intersect, split(C_vir_apop_LFC_BAC_path$transcript_id, C_vir_apop_LFC_BAC_path$group_by_sim))
  # "XM_022460059.1" "XM_022467669.1" "XM_022474529.1" "XM_022463142.1"
C_vir_apop_LFC_BAC_path_core_product <- Reduce(intersect, split(C_vir_apop_LFC_BAC_path$product, C_vir_apop_LFC_BAC_path$group_by_sim))
#"mitochondrial Rho GTPase 1-like, transcript variant X3"                          "mitogen-activated protein kinase kinase kinase 2-like"                          
# "signal transducer and activator of transcription 5B-like, transcript variant X3" "stress-activated protein kinase JNK-like, transcript variant X4" 

C_vir_apop_LFC_BAC_pro <- C_vir_apop_LFC_bac %>% filter(group_by_sim == "Probiotic" | group_by_sim == "RI_6h" | group_by_sim == "RI_24h" | group_by_sim == "S4_6h"| group_by_sim == "S4_24h")
C_vir_apop_LFC_BAC_pro_core <- Reduce(intersect, split(C_vir_apop_LFC_BAC_pro $transcript_id, C_vir_apop_LFC_BAC_pro$group_by_sim))
  # 0 
C_vir_apop_LFC_BAC_pro_core_product <- Reduce(intersect, split(C_vir_apop_LFC_BAC_pro $product, C_vir_apop_LFC_BAC_pro$group_by_sim))
  # 0

C_vir_apop_LFC_BAC_pro_RI <- C_vir_apop_LFC_bac %>% filter(group_by_sim == "Probiotic" | group_by_sim == "RI_6h" | group_by_sim == "RI_24h")
C_vir_apop_LFC_BAC_pro_RI_core <- Reduce(intersect, split(C_vir_apop_LFC_BAC_pro_RI$transcript_id, C_vir_apop_LFC_BAC_pro_RI$group_by_sim))
  # "XM_022437239.1" = mitogen-activated protein kinase kinase kinase 7-like, transcript variant X2
C_vir_apop_LFC_BAC_pro_RI_core_product <- Reduce(intersect, split(C_vir_apop_LFC_BAC_pro_RI$product, C_vir_apop_LFC_BAC_pro_RI$group_by_sim))
  # [1] "mitogen-activated protein kinase kinase kinase 7-like, transcript variant X2"

C_vir_apop_LFC_BAC_pro_S4 <- C_vir_apop_LFC_bac %>% filter(group_by_sim == "S4_6h" | group_by_sim == "S4_24h")
C_vir_apop_LFC_BAC_pro_S4_core <- Reduce(intersect, split(C_vir_apop_LFC_BAC_pro_S4$transcript_id,C_vir_apop_LFC_BAC_pro_S4$group_by_sim))
C_vir_apop_LFC_BAC_pro_S4_core <- as.data.frame(C_vir_apop_LFC_BAC_pro_S4_core )
colnames(C_vir_apop_LFC_BAC_pro_S4_core )[1] <- "transcript_id"  
C_vir_apop_LFC_BAC_pro_S4_core_name <- left_join(C_vir_apop_LFC_BAC_pro_S4_core, C_vir_apop_LFC_bac, by = "transcript_id")
View(unique(C_vir_apop_LFC_BAC_pro_S4_core_name$product))
# caspase-3-like, transcript variant X1
# adenylate cyclase type 9-like
# baculoviral IAP repeat-containing protein 6-like, transcript variant X5
# calpain-B-like, transcript variant X19
# calpain-B-like, transcript variant X23
# eukaryotic translation initiation factor 2-alpha kinase 3-like, transcript variant X2
# tumor necrosis factor alpha-induced protein 3-like, transcript variant X3
# signal transducer and activator of transcription 5B-like, transcript variant X3
# baculoviral IAP repeat-containing protein 7-like, transcript variant X4
# ceramide synthase 5-like, transcript variant X3
# interferon alpha-inducible protein 27-like protein 2B, transcript variant X1
# rho-related protein racA-like
# lipopolysaccharide-induced tumor necrosis factor-alpha factor-like, transcript variant X1
# myeloid differentiation primary response protein MyD88-like, transcript variant X1

C_vir_apop_LFC_BAC_pro_S4_core_product <- Reduce(intersect, split(C_vir_apop_LFC_BAC_pro_S4$product,C_vir_apop_LFC_BAC_pro_S4$group_by_sim))
 # same as above

# Investigate C_gig bacterial response
C_gig_apop_LFC_bac <- C_gig_apop_LFC %>% filter(experiment == "Zhang" | experiment == "Rubio" )
C_gig_apop_LFC_bac_core <- Reduce(intersect, split(C_gig_apop_LFC_bac$Name, C_gig_apop_LFC_bac$experiment))
  # "XM_011414513.2" "XM_011425781.1" "XM_011456317.2" "XM_011456757.2" "XM_020072623.1" "XM_011447078.2" "XM_020070534.1" "XM_011426791.2" "XM_011451291.2" "XM_011456852.2" 
C_gig_apop_LFC_bac_core_product <- Reduce(intersect, split(C_gig_apop_LFC_bac$product, C_gig_apop_LFC_bac$experiment))
# hexokinase-1
# cdc42 homolog
# myeloid differentiation primary response protein MyD88-like, transcript variant X2
# GTPase IMAP family member 4, transcript variant X3
# transcription factor AP-1, transcript variant X2
# baculoviral IAP repeat-containing protein 3-like, transcript variant X1
# adenylate cyclase type 2, transcript variant X2
# cyclic AMP-responsive element-binding protein 3-like protein 3-B, transcript variant X2
# interferon regulatory factor 1, transcript variant X2
# tyrosine-protein phosphatase non-receptor type 13, transcript variant X4

C_gig_apop_LFC_BAC_path <- C_gig_apop_LFC_bac %>% filter(group_by_sim == "V_tub_V_ang" | group_by_sim == "V_aes_V_alg1_V_alg2" | group_by_sim == "J2-9 vir" | group_by_sim == "LGP32 vir")
C_gig_apop_LFC_BAC_path_core <- Reduce(intersect, split(C_gig_apop_LFC_BAC_path $Name, C_gig_apop_LFC_BAC_path $group_by_sim))
  # "XM_011425781.1" "XM_020070534.1" "XM_020072623.1"
C_gig_apop_LFC_BAC_path_core_product <- Reduce(intersect, split(C_gig_apop_LFC_BAC_path $product, C_gig_apop_LFC_BAC_path $group_by_sim))
  # [1] "cdc42 homolog" "caspase-7-like" "adenylate cyclase type 2, transcript variant X2" "transcription factor AP-1, transcript variant X2"

C_gig_apop_LFC_BAC_nonpath <- C_gig_apop_LFC_bac %>% filter(group_by_sim == "LPS_M_lut" | group_by_sim == "J2-8 non-vir" | group_by_sim == "LMG20012T non-vir")
C_gig_apop_LFC_BAC_nonpath_core <- Reduce(intersect, split(C_gig_apop_LFC_BAC_nonpath$Name, C_gig_apop_LFC_BAC_nonpath$group_by_sim))
  # "XM_011425781.1" "XM_020072623.1"
C_gig_apop_LFC_BAC_nonpath_core_product <- Reduce(intersect, split(C_gig_apop_LFC_BAC_nonpath$product, C_gig_apop_LFC_BAC_nonpath$group_by_sim))
  #  "cdc42 homolog" "transcription factor AP-1, transcript variant X2"

# The lists of shared transcript IDs and shared product names is more similar in C. gigas than in C. virginica...maybe there is just higher transcript diversity in C. virginica?

# Bacterial list upset plot 
# helpful tutorial for doing this: http://genomespot.blogspot.com/2017/09/upset-plots-as-replacement-to-venn.html
# http://crazyhottommy.blogspot.com/2016/01/upset-plot-for-overlapping-chip-seq.html
# https://www.littlemissdata.com/blog/set-analysis
# UpsetR vignette: https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html
C_vir_apop_LFC_upset <- C_vir_apop_LFC[,c("transcript_id","group_by_sim")]
C_gig_apop_LFC_upset <- C_gig_apop_LFC[,c("Name","group_by_sim")]
C_vir_apop_LFC_bac_upset <- C_vir_apop_LFC_bac[,c("transcript_id","group_by_sim")]
C_gig_apop_LFC_bac_upset <- C_gig_apop_LFC_bac[,c("Name","group_by_sim")]

# Convert into wide format using reshape
C_vir_apop_LFC_upset_wide <- C_vir_apop_LFC_upset %>% mutate(value=1) %>% spread(group_by_sim, value, fill =0 )
head(C_vir_apop_LFC_upset_wide)
C_vir_apop_LFC_bac_upset_wide <- C_vir_apop_LFC_bac_upset %>% mutate(value=1) %>% spread(group_by_sim, value, fill =0 )
head(C_vir_apop_LFC_bac_upset_wide)

C_gig_apop_LFC_upset_wide <- C_gig_apop_LFC_upset %>% mutate(value=1) %>% spread(group_by_sim, value, fill =0 )
head(C_gig_apop_LFC_upset_wide)
C_gig_apop_LFC_bac_upset_wide <- C_gig_apop_LFC_bac_upset %>% mutate(value=1) %>% spread(group_by_sim, value, fill =0 )
head(C_gig_apop_LFC_bac_upset_wide)

# Make upset plots
C_vir_apop_LFC_upset_wide_GROUP <- upset(C_vir_apop_LFC_upset_wide,  mainbar.y.label = "Transcript id Intersections", 
      sets.x.label = "# Significantly Differentially Expressed Apoptosis Transcripts", text.scale = c(1.3, 1.3, 1.3, 1.3, 2, 1.0),
      sets= c("Probiotic","RE22","RI_24h","RI_6h","ROD_resistant","ROD_susceptible","S4_24h", "S4_6h","Susceptible_28d",
              "Susceptible_36hr", "Tolerant_28d","Tolerant_36hr","Tolerant_7d"),
       order.by="freq")
C_vir_apop_LFC_bac_upset_wide_GROUP <- upset(C_vir_apop_LFC_bac_upset_wide,  mainbar.y.label = "Transcript id Intersections", 
                                         sets.x.label = "# Significantly Differentially Expressed Apoptosis Transcripts", text.scale = c(1.3, 1.3, 1.3, 1.3, 2, 1.0),
                                         sets= c("Probiotic","RE22","RI_24h","RI_6h","ROD_resistant","ROD_susceptible",
                                                 "S4_24h","S4_6h"),order.by="freq")
C_gig_apop_LFC_upset_wide_GROUP <- upset(C_gig_apop_LFC_upset_wide,  mainbar.y.label = "Transcript id Intersections", 
                                         sets.x.label = "# Significantly Differentially Expressed Apoptosis Transcripts", text.scale = c(1.3, 1.3, 1.3, 1.3, 2, 1.0),
                                         sets= c("He_dds_res_120hr_sig_APOP" ,"He_dds_res_12hr_sig_APOP",  "He_dds_res_24hr_sig_APOP",  "He_dds_res_48hr_sig_APOP",
                                                 "He_dds_res_6hr_sig_APOP"   ,"J2-8 non-vir"            ,  "J2-9 vir"                ,  "LGP32 vir"               , "LMG20012T non-vir"       , 
                                                     "Resistant_dds_res_12"   ,   "Resistant_dds_res_24"   ,   "Resistant_dds_res_48"   ,  "Resistant_dds_res_6"    ,  
                                                  "Resistant_dds_res_60"     , "Resistant_dds_res_72"   ,   "Susceptible_dds_res_12" ,   "Susceptible_dds_res_24" ,  "Susceptible_dds_res_48" ,  
                                                  "Susceptible_dds_res_6"    , "Susceptible_dds_res_60" ,   "Susceptible_dds_res_72" ,   "V_aes_V_alg1_V_alg2"    ,  "V_tub_V_ang", "LPS_M_lut"   ),
                                         order.by="freq")
                # Most unique sets belong to the early delorgeril response
                # greatest between experiment interaction is between the Susceptible delorgeril 24 and 60 hr, and all the Rubio samples, 5 intersections 
C_gig_apop_LFC_bac_upset_wide_GROUP <- upset(C_gig_apop_LFC_bac_upset_wide,  mainbar.y.label = "Transcript id Intersections", 
                                             sets.x.label = "# Significantly Differentially Expressed Apoptosis Transcripts", text.scale = c(1.3, 1.3, 1.3, 1.3, 2, 1.0),
                                             sets= c("J2-8 non-vir","J2-9 vir" ,"LGP32 vir","LMG20012T non-vir","LPS_M_lut","V_aes_V_alg1_V_alg2",
                                                     "V_tub_V_ang" ),order.by="freq")
                # greatest sharing of all transcripts is within experiments, most unique sets in the J2-8, J2-9 and LGP32 

## Extract products that appear only once 
C_vir_apop_LFC_notshared  <- C_vir_apop_LFC[!(duplicated(C_vir_apop_LFC$transcript_id) | duplicated(C_vir_apop_LFC$transcript_id, fromLast = TRUE)), ]
C_gig_apop_LFC_notshared <- C_gig_apop_LFC[!(duplicated(C_gig_apop_LFC$Name) | duplicated(C_gig_apop_LFC$Name, fromLast = TRUE)), ]


### Comparative tables of apoptosis gene family expression ###
### Apr. 6th, 2020 - need to check the frequency table code below, not sure if it's working as expected

C_vir_apop_LFC_summmary_apop_transcripts_table <- C_vir_apop_LFC %>% group_by(experiment, group_by_sim) %>% mutate(number_apop_transcripts=n())
C_gig_apop_LFC_summmary_apop_transcripts_table <- C_gig_apop_LFC %>% group_by(experiment, group_by_sim) %>% mutate(number_apop_transcripts=n())

Combined_summmary_apop_transcripts_table <- rbind(C_vir_apop_LFC_summmary_apop_transcripts_table,C_gig_apop_LFC_summmary_apop_transcripts_table)

C_vir_apop_LFC_per_group <- C_vir_apop_LFC %>% group_by(experiment, group_by_sim) %>% mutate(transcripts_per_group = n())
C_gig_apop_LFC_per_group <- C_gig_apop_LFC %>% group_by(experiment, group_by_sim) %>% mutate(transcripts_per_group = n())

# Number GIMAP per challenge
C_vir_apop_LFC_GIMAP <- C_vir_apop_LFC[grepl("IMAP", C_vir_apop_LFC$product, ignore.case = TRUE),]
C_gig_apop_LFC_GIMAP <- C_gig_apop_LFC[grepl("IMAP", C_gig_apop_LFC$product, ignore.case = TRUE),]
# per condition frequency
# C_vir_apop_LFC_GIMAP_freq <- C_vir_apop_LFC_GIMAP %>% group_by(experiment, group_by_sim) %>% filter(grepl("IMAP", product,ignore.case = TRUE)) %>% mutate(GIMAP_number = n())
# C_gig_apop_LFC_GIMAP_freq <- C_gig_apop_LFC_GIMAP %>% group_by(experiment, group_by_sim) %>% filter(grepl("IMAP", product,ignore.case = TRUE)) %>% mutate(GIMAP_number = n())
# C_vir_apop_LFC_GIMAP_freq <- left_join(C_vir_apop_LFC_GIMAP_freq,C_vir_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(GIMAP_freq = (GIMAP_number/transcripts_per_group))
# C_gig_apop_LFC_GIMAP_freq <- left_join(C_gig_apop_LFC_GIMAP_freq,C_gig_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(GIMAP_freq = (GIMAP_number/transcripts_per_group))
# # compare frequency
# GIMAP_freq_combined <- rbind(C_vir_apop_LFC_GIMAP_freq ,C_gig_apop_LFC_GIMAP_freq )

# Number IAP per challenge
C_vir_apop_LFC_IAP <- C_vir_apop_LFC[grepl("IAP", C_vir_apop_LFC$product, ignore.case = TRUE),]
C_gig_apop_LFC_IAP <- C_gig_apop_LFC[grepl("IAP", C_gig_apop_LFC$product, ignore.case = TRUE),]
# per condition frequency
# C_vir_apop_LFC_IAP_freq <- C_vir_apop_LFC_IAP %>% group_by(experiment, group_by_sim) %>% filter(grepl("IAP", product,ignore.case = TRUE)) %>% mutate(IAP_number = n())
# C_gig_apop_LFC_IAP_freq <- C_gig_apop_LFC_IAP %>% group_by(experiment, group_by_sim) %>% filter(grepl("IAP", product,ignore.case = TRUE)) %>% mutate(IAP_number = n())
# C_vir_apop_LFC_IAP_freq <- left_join(C_vir_apop_LFC_IAP_freq,C_vir_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(IAP_freq = (IAP_number/transcripts_per_group))
# C_gig_apop_LFC_IAP_freq <- left_join(C_gig_apop_LFC_IAP_freq,C_gig_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(IAP_freq = (IAP_number/transcripts_per_group))
# # compare frequency
# IAP_freq_combined <- rbind(C_vir_apop_LFC_IAP_freq ,C_gig_apop_LFC_IAP_freq )

# Number caspase per challenge
C_vir_apop_LFC_caspase <- C_vir_apop_LFC[grepl("caspase", C_vir_apop_LFC$product, ignore.case = TRUE),]
C_gig_apop_LFC_caspase <- C_gig_apop_LFC[grepl("caspase", C_gig_apop_LFC$product, ignore.case = TRUE),]
# per condition frequency
#C_vir_apop_LFC_caspase_freq <- C_vir_apop_LFC_caspase %>% group_by(experiment, group_by_sim) %>% filter(grepl("caspase", product,ignore.case = TRUE)) %>% mutate(caspase_number = n())
#C_gig_apop_LFC_caspase_freq <- C_gig_apop_LFC_caspase %>% group_by(experiment, group_by_sim) %>% filter(grepl("caspase", product,ignore.case = TRUE)) %>% mutate(caspase_number = n())
#C_vir_apop_LFC_caspase_freq <- left_join(C_vir_apop_LFC_caspase_freq,C_vir_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(caspase_freq = (caspase_number/transcripts_per_group))
#C_gig_apop_LFC_caspase_freq <- left_join(C_gig_apop_LFC_caspase_freq,C_gig_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(caspase_freq = (caspase_number/transcripts_per_group))
## compare frequency
#caspase_freq_combined <- rbind(C_vir_apop_LFC_caspase_freq ,C_gig_apop_LFC_caspase_freq )

# Number TNF family per challenge
C_vir_apop_LFC_TNF <- C_vir_apop_LFC[grepl("TNF", C_vir_apop_LFC$product, ignore.case = TRUE) | grepl("tumor", C_vir_apop_LFC$product, ignore.case = TRUE),]
C_gig_apop_LFC_TNF <- C_gig_apop_LFC[grepl("TNF", C_gig_apop_LFC$product, ignore.case = TRUE) | grepl("tumor", C_gig_apop_LFC$product, ignore.case = TRUE),]
#per condition frequency
#C_vir_apop_LFC_TNF_freq <- C_vir_apop_LFC_TNF %>% group_by(experiment, group_by_sim) %>% filter(grepl("TNF", product,ignore.case = TRUE) |grepl("tumor", product,ignore.case = TRUE) ) %>% mutate(TNF_number = n())
#C_gig_apop_LFC_TNF_freq <- C_gig_apop_LFC_TNF %>% group_by(experiment, group_by_sim) %>% filter(grepl("TNF", product,ignore.case = TRUE) |grepl("tumor", product,ignore.case = TRUE) ) %>% mutate(TNF_number = n())
#C_vir_apop_LFC_TNF_freq <- left_join(C_vir_apop_LFC_TNF_freq,C_vir_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(tnf_freq = (TNF_number/transcripts_per_group))
#C_gig_apop_LFC_TNF_freq <- left_join(C_gig_apop_LFC_TNF_freq,C_gig_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(tnf_freq = (TNF_number/transcripts_per_group))
## compare frequency
#TNF_freq_combined <- rbind(C_vir_apop_LFC_TNF_freq ,C_gig_apop_LFC_TNF_freq )

# Number Toll family per challenge
C_vir_apop_LFC_Toll <- C_vir_apop_LFC[grepl("toll", C_vir_apop_LFC$product, ignore.case = TRUE),]
C_gig_apop_LFC_Toll <- C_gig_apop_LFC[grepl("toll", C_gig_apop_LFC$product, ignore.case = TRUE),]
#per condition frequency
# C_vir_apop_LFC_toll_freq <- C_vir_apop_LFC_Toll %>% group_by(experiment, group_by_sim) %>% filter(grepl("toll", product,ignore.case = TRUE)) %>% mutate(toll_number = n())
# C_gig_apop_LFC_toll_freq <- C_gig_apop_LFC_Toll %>% group_by(experiment, group_by_sim) %>% filter(grepl("toll", product,ignore.case = TRUE)) %>% mutate(toll_number = n())
# C_vir_apop_LFC_toll_freq <- left_join(C_vir_apop_LFC_toll_freq,C_vir_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(toll_freq = (toll_number/transcripts_per_group))
# C_gig_apop_LFC_toll_freq <- left_join(C_gig_apop_LFC_toll_freq,C_gig_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(toll_freq = (toll_number/transcripts_per_group))
# # compare frequency
# Toll_freq_combined <- rbind(C_vir_apop_LFC_toll_freq ,C_gig_apop_LFC_toll_freq )

# Number interferon per challenge
C_vir_apop_LFC_interferon <- C_vir_apop_LFC[grepl("interferon", C_vir_apop_LFC$product, ignore.case = TRUE),]
C_gig_apop_LFC_interferon <- C_gig_apop_LFC[grepl("interferon", C_gig_apop_LFC$product, ignore.case = TRUE),]
#per condition frequency
# C_vir_apop_LFC_interferon_freq <- C_vir_apop_LFC_interferon %>% group_by(experiment, group_by_sim) %>% filter(grepl("interferon", product,ignore.case = TRUE)) %>% mutate(interferon_number = n())
# C_gig_apop_LFC_interferon_freq <- C_gig_apop_LFC_interferon %>% group_by(experiment, group_by_sim) %>% filter(grepl("interferon", product,ignore.case = TRUE)) %>% mutate(interferon_number = n())
# C_vir_apop_LFC_interferon_freq <- left_join(C_vir_apop_LFC_interferon_freq,C_vir_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(interferon_freq = (interferon_number/transcripts_per_group))
# C_gig_apop_LFC_interferon_freq <- left_join(C_gig_apop_LFC_interferon_freq,C_gig_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(interferon_freq = (interferon_number/transcripts_per_group))
# # compare frequency
# interferon_freq_combined <- rbind(C_vir_apop_LFC_interferon_freq ,C_gig_apop_LFC_interferon_freq )

# Number heatshock per challenge
C_vir_apop_LFC_hsp <- C_vir_apop_LFC[grepl("heat", C_vir_apop_LFC$product, ignore.case = TRUE),]
C_gig_apop_LFC_hsp <- C_gig_apop_LFC[grepl("heat", C_gig_apop_LFC$product, ignore.case = TRUE),]
#per condition frequency
# C_vir_apop_LFC_hsp_freq <- C_vir_apop_LFC_hsp %>% group_by(experiment, group_by_sim) %>% filter(grepl("heat", product,ignore.case = TRUE)) %>% mutate(hsp_number = n())
# C_gig_apop_LFC_hsp_freq <- C_gig_apop_LFC_hsp %>% group_by(experiment, group_by_sim) %>% filter(grepl("heat", product,ignore.case = TRUE)) %>% mutate(hsp_number = n())
# C_vir_apop_LFC_hsp_freq <- left_join(C_vir_apop_LFC_hsp_freq,C_vir_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(hsp_freq = (hsp_number/transcripts_per_group))
# C_gig_apop_LFC_hsp_freq <- left_join(C_gig_apop_LFC_hsp_freq,C_gig_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(hsp_freq = (hsp_number/transcripts_per_group))
# # compare frequency
# Hsp_freq_combined <- rbind(C_vir_apop_LFC_hsp_freq ,C_gig_apop_LFC_hsp_freq )

# Number cathepsin per challenge
C_vir_apop_LFC_cathepsin <- C_vir_apop_LFC[grepl("cathepsin", C_vir_apop_LFC$product, ignore.case = TRUE),]
C_gig_apop_LFC_cathepsin <- C_gig_apop_LFC[grepl("cathepsin", C_gig_apop_LFC$product, ignore.case = TRUE),]
#per condition frequency
# C_vir_apop_LFC_cathepsin_freq <- C_vir_apop_LFC_cathepsin %>% group_by(experiment, group_by_sim) %>% filter(grepl("cathepsin", product,ignore.case = TRUE)) %>% mutate(cathepsin_number = n())
# C_gig_apop_LFC_cathepsin_freq <- C_gig_apop_LFC_cathepsin %>% group_by(experiment, group_by_sim) %>% filter(grepl("cathepsin", product,ignore.case = TRUE)) %>% mutate(cathepsin_number = n())
# C_vir_apop_LFC_cathepsin_freq <- left_join(C_vir_apop_LFC_cathepsin_freq,C_vir_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(cathepsin_freq = (cathepsin_number/transcripts_per_group))
# C_gig_apop_LFC_cathepsin_freq <- left_join(C_gig_apop_LFC_cathepsin_freq,C_gig_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(cathepsin_freq = (cathepsin_number/transcripts_per_group))
# # compare frequency
# Cathepsin_freq_combined <- rbind(C_vir_apop_LFC_cathepsin_freq ,C_gig_apop_LFC_cathepsin_freq )

# Number programmed cell death protein per challenge
C_vir_apop_LFC_PCDP <- C_vir_apop_LFC[grepl("programmed", C_vir_apop_LFC$product, ignore.case = TRUE),]
C_gig_apop_LFC_PCDP <- C_gig_apop_LFC[grepl("programmed", C_gig_apop_LFC$product, ignore.case = TRUE),]
#per condition frequency
# C_vir_apop_LFC_PCDP_freq <- C_vir_apop_LFC_PCDP %>% group_by(experiment, group_by_sim) %>% filter(grepl("programmed", product,ignore.case = TRUE)) %>% mutate(PCDP_number = n())
# C_gig_apop_LFC_PCDP_freq <- C_gig_apop_LFC_PCDP %>% group_by(experiment, group_by_sim) %>% filter(grepl("programmed", product,ignore.case = TRUE)) %>% mutate(PCDP_number = n())
# C_vir_apop_LFC_PCDP_freq <- left_join(C_vir_apop_LFC_PCDP_freq,C_vir_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(PCDP_freq = (PCDP_number/transcripts_per_group))
# C_gig_apop_LFC_PCDP_freq <- left_join(C_gig_apop_LFC_PCDP_freq,C_gig_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(PCDP_freq = (PCDP_number/transcripts_per_group))
# # compare frequency
# PCDP_freq_combined <- rbind(C_vir_apop_LFC_PCDP_freq ,C_gig_apop_LFC_PCDP_freq )

# Number rho-related
C_vir_apop_LFC_rho <- C_vir_apop_LFC[grepl("rho", C_vir_apop_LFC$product, ignore.case = TRUE) | grepl("ras", C_vir_apop_LFC$product, ignore.case = TRUE) & 
                                       !grepl("polymerase",C_vir_apop_LFC$product, ignore.case = TRUE) ,]
C_gig_apop_LFC_rho <- C_gig_apop_LFC[grepl("rho", C_gig_apop_LFC$product, ignore.case = TRUE) | grepl("ras", C_gig_apop_LFC$product, ignore.case = TRUE) &
                                       !grepl("polymerase",C_gig_apop_LFC$product, ignore.case = TRUE) ,]
#per condition frequency
# C_vir_apop_LFC_rho_freq <- C_vir_apop_LFC_rho %>% group_by(experiment, group_by_sim) %>% filter(grepl("rho", product,ignore.case = TRUE) |grepl("ras", product,ignore.case = TRUE) & !grepl("polymerase", product, ignore.case = TRUE)) %>% mutate(rho_number = n())
# C_gig_apop_LFC_rho_freq <- C_gig_apop_LFC_rho %>% group_by(experiment, group_by_sim) %>% filter(grepl("rho", product,ignore.case = TRUE) |grepl("ras", product,ignore.case = TRUE) & !grepl("polymerase", product, ignore.case = TRUE)) %>% mutate(rho_number = n())
# C_vir_apop_LFC_rho_freq <- left_join(C_vir_apop_LFC_rho_freq,C_vir_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(rho_freq = (rho_number/transcripts_per_group))
# C_gig_apop_LFC_rho_freq <- left_join(C_gig_apop_LFC_rho_freq,C_gig_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(rho_freq = (rho_number/transcripts_per_group))
# # compare frequency
# Rho_freq_combined <- rbind(C_vir_apop_LFC_rho_freq ,C_gig_apop_LFC_rho_freq )

# Number  myD88
C_vir_apop_LFC_myD88 <- C_vir_apop_LFC[grepl("myD88", C_vir_apop_LFC$product, ignore.case = TRUE),]
C_gig_apop_LFC_myD88 <- C_gig_apop_LFC[grepl("myD88", C_gig_apop_LFC$product, ignore.case = TRUE),]
#per condition frequency
# C_vir_apop_LFC_myD88_freq <- C_vir_apop_LFC_myD88 %>% group_by(experiment, group_by_sim) %>% filter(grepl("myD88", product,ignore.case = TRUE)) %>% mutate(myD88_number = n())
# C_gig_apop_LFC_myD88_freq <- C_gig_apop_LFC_myD88 %>% group_by(experiment, group_by_sim) %>% filter(grepl("myD88", product,ignore.case = TRUE)) %>% mutate(myD88_number = n())
# C_vir_apop_LFC_myD88_freq <- left_join(C_vir_apop_LFC_myD88_freq,C_vir_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(myD88_freq = (myD88_number/transcripts_per_group))
# C_gig_apop_LFC_myD88_freq <- left_join(C_gig_apop_LFC_myD88_freq,C_gig_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(myD88_freq = (myD88_number/transcripts_per_group))
# # compare frequency
# myD88_freq_combined <- rbind(C_vir_apop_LFC_myD88_freq ,C_gig_apop_LFC_myD88_freq )

# Number inositol
C_vir_apop_LFC_inositol <- C_vir_apop_LFC[grepl("inositol", C_vir_apop_LFC$product, ignore.case = TRUE) & !grepl("kinase", C_vir_apop_LFC$product, ignore.case = TRUE),]
C_gig_apop_LFC_inositol <- C_gig_apop_LFC[grepl("inositol", C_gig_apop_LFC$product, ignore.case = TRUE) & !grepl("kinase",C_gig_apop_LFC$product, ignore.case = TRUE),]
#per condition frequency
# C_vir_apop_LFC_inositol_freq <- C_vir_apop_LFC_inositol %>% group_by(experiment, group_by_sim) %>% filter(grepl("inositol", product,ignore.case = TRUE) & !grepl("kinase",product, ignore.case=TRUE)) %>% mutate(inositol_number = n())
# C_gig_apop_LFC_inositol_freq <- C_gig_apop_LFC_inositol %>% group_by(experiment, group_by_sim) %>% filter(grepl("inositol", product,ignore.case = TRUE) & !grepl("kinase",product, ignore.case=TRUE)) %>% mutate(inositol_number = n())
# C_vir_apop_LFC_inositol_freq <- left_join(C_vir_apop_LFC_inositol_freq,C_vir_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(inositol_freq = (inositol_number/transcripts_per_group))
# C_gig_apop_LFC_inositol_freq <- left_join(C_gig_apop_LFC_inositol_freq,C_gig_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(inositol_freq = (inositol_number/transcripts_per_group))
# # compare frequency
# inositol_freq_combined <- rbind(C_vir_apop_LFC_inositol_freq ,C_gig_apop_LFC_inositol_freq )

# Number adenylate cyclase
C_vir_apop_LFC_adenylate <- C_vir_apop_LFC[grepl("adenylate", C_vir_apop_LFC$product, ignore.case = TRUE),]
C_gig_apop_LFC_adenylate <- C_gig_apop_LFC[grepl("adenylate", C_gig_apop_LFC$product, ignore.case = TRUE),]
#per condition frequency
# C_vir_apop_LFC_adenylate_freq <- C_vir_apop_LFC_adenylate %>% group_by(experiment, group_by_sim) %>% filter(grepl("adenylate", product,ignore.case = TRUE)) %>% mutate(adenylate_number = n())
# C_gig_apop_LFC_adenylate_freq <- C_gig_apop_LFC_adenylate %>% group_by(experiment, group_by_sim) %>% filter(grepl("adenylate", product,ignore.case = TRUE)) %>% mutate(adenylate_number = n())
# C_vir_apop_LFC_adenylate_freq <- left_join(C_vir_apop_LFC_adenylate_freq,C_vir_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(adenylate_freq = (adenylate_number/transcripts_per_group))
# C_gig_apop_LFC_adenylate_freq <- left_join(C_gig_apop_LFC_adenylate_freq,C_gig_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(adenylate_freq = (adenylate_number/transcripts_per_group))
# # compare frequency
# adenylate_freq_combined <- rbind(C_vir_apop_LFC_adenylate_freq ,C_gig_apop_LFC_adenylate_freq )

# Number apoptosis inhibitor
C_vir_apop_LFC_ApopI <- C_vir_apop_LFC[grepl("apoptosis inhibitor", C_vir_apop_LFC$product, ignore.case = TRUE),]
C_gig_apop_LFC_ApopI <- C_gig_apop_LFC[grepl("apoptosis inhibitor", C_gig_apop_LFC$product, ignore.case = TRUE),]
#per condition frequency
# C_vir_apop_LFC_ApopI_freq <- C_vir_apop_LFC_ApopI %>% group_by(experiment, group_by_sim) %>% filter(grepl("apoptosis inhibitor", product,ignore.case = TRUE)) %>% mutate(ApopI_number = n())
# C_gig_apop_LFC_ApopI_freq <- C_gig_apop_LFC_ApopI %>% group_by(experiment, group_by_sim) %>% filter(grepl("apoptosis inhibitor", product,ignore.case = TRUE)) %>% mutate(ApopI_number = n())
# C_vir_apop_LFC_ApopI_freq <- left_join(C_vir_apop_LFC_ApopI_freq,C_vir_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(ApopI_freq = (ApopI_number/transcripts_per_group))
# C_gig_apop_LFC_ApopI_freq <- left_join(C_gig_apop_LFC_ApopI_freq,C_gig_apop_LFC_per_group, by = c("experiment","group_by_sim")) %>% mutate(ApopI_freq = (ApopI_number/transcripts_per_group))
# # compare frequency
# ApopI_freq_combined <- rbind(C_vir_apop_LFC_ApopI_freq ,C_gig_apop_LFC_ApopI_freq )

## Join all gene family frequency tables 
#Gene_family_frequency <- full_join(GIMAP_freq_combined[,c("experiment","group_by_sim","GIMAP_freq")], IAP_freq_combined[,c("experiment","group_by_sim", "IAP_freq")], 
#                                   by = c("experiment","group_by_sim"))  
#Gene_family_frequency <- full_join(Gene_family_frequency, caspase_freq_combined[,c("experiment","group_by_sim","caspase_freq")], 
#                                   by = c("experiment","group_by_sim")) 
#Gene_family_frequency <- full_join(Gene_family_frequency, TNF_freq_combined[,c("experiment","group_by_sim","tnf_freq")], 
#                                   by = c("experiment","group_by_sim")) 
#Gene_family_frequency <- full_join(Gene_family_frequency, Toll_freq_combined [,c("experiment","group_by_sim","toll_freq")], 
#                                   by = c("experiment","group_by_sim")) 
#Gene_family_frequency <- full_join(Gene_family_frequency, interferon_freq_combined[,c("experiment","group_by_sim","interferon_freq")], 
#                                   by = c("experiment","group_by_sim")) 
#Gene_family_frequency <- full_join(Gene_family_frequency, Hsp_freq_combined[,c("experiment","group_by_sim","hsp_freq")], 
#                                   by = c("experiment","group_by_sim"))  
#Gene_family_frequency <- full_join(Gene_family_frequency, Cathepsin_freq_combined[,c("experiment","group_by_sim","cathepsin_freq")], 
#                                   by = c("experiment","group_by_sim"))  
#Gene_family_frequency <- full_join(Gene_family_frequency, PCDP_freq_combined[,c("experiment","group_by_sim","PCDP_freq")], 
#                                   by = c("experiment","group_by_sim"))  
#Gene_family_frequency <- full_join(Gene_family_frequency, Rho_freq_combined  [,c("experiment","group_by_sim","rho_freq")], 
#                                   by = c("experiment","group_by_sim"))  
#Gene_family_frequency <- full_join(Gene_family_frequency, myD88_freq_combined[,c("experiment","group_by_sim","myD88_freq")], 
#                                   by = c("experiment","group_by_sim"))  
#Gene_family_frequency <- full_join(Gene_family_frequency, inositol_freq_combined[,c("experiment","group_by_sim","inositol_freq")], 
#                                   by = c("experiment","group_by_sim")) 
#Gene_family_frequency <- full_join(Gene_family_frequency, adenylate_freq_combined[,c("experiment","group_by_sim","adenylate_freq")], 
#                                   by = c("experiment","group_by_sim")) 
#Gene_family_frequency <- full_join(Gene_family_frequency, ApopI_freq_combined [,c("experiment","group_by_sim","ApopI_freq")], 
#                                   by = c("experiment","group_by_sim")) 
#Gene_family_frequency_long <- melt(Gene_family_frequency, id.vars= c("experiment","group_by_sim"), measure.vars=c(3:16))
#is.na(Gene_family_frequency_long) <- "0"
#Gene_family_frequency_long_plot <- ggplot(Gene_family_frequency_long, aes(x=experiment, y = value, fill=experiment)) + 
#  geom_col(position="dodge") + facet_grid(.~variable) + theme(axis.text.x = element_text(angle = 90, hjust = 1))


#### SESSION INFO FOR RUNNING SCRIPTS Spring 2020 ####

sessionInfo()
