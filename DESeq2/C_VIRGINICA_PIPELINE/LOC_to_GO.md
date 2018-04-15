# *C. virginica* Add GO ID and gene info
## April 15, 2018
## By: Erin Roberts

### GOAL:
1. extract transcript info from GFF 
2. merge LOC ID from gff3 to LOC ID info from merged.gtf file from merged gtf file 
3. use referencen (XP or XM) transcript info in this new file to map the interproscan GO IDs onto the transcript ids 

Desired Output file = all the transcripts and exons found with HISAT/Stringtie using Cvir transcriptomes are mapped to the GFF file information and have the InterProScan information added

### Step 1: Load reference GFF3 (/data3/marine_diseases_lab/shared/ref_C_virginica-3.0_top_level.gff3) 

```
### Load Libraries
library(dplyr)
library(tidyr)
install.packages("gtools")
library(gtools)
library(stringr)

#read in reference GFF3 file, with the header lines removed from file (sed '1,9d' ref_C_virginica-3.0_top_level.gff3 > edited_ref_C_virginica-3.0_top_level.gff3)
Cvir_reference <- read.csv("edited_ref_C_virginica-3.0_top_level.gff3", sep="\t", header=FALSE)

# tell R what the column names should be
colnames(Cvir_reference) <- c('seqname', 'source', 'feature','start','end','score','strand','frame', 'attribute')

# subset just the transcript lines and not the exon lines 
Cvir_reference_transcripts_exons <- Cvir_reference[Cvir_reference$feature %in% c("transcript","exon"),]

# separate attributes column with the ; separator in tidyr
# use the extra  argument to control what happens when every row doesn't split into the same number of pieces
Cvir_reference_transcripts_exons_separated <- Cvir_reference_transcripts_exons  %>% separate(attribute, c("ID","Parent_RNA","Dbxref_GeneID", "gbkey","geneLOC", "product","transcript_id"), ";", extra="merge")   

# separate gene name into Dbxref and Genbank
Cvir_reference_transcripts_exons_separated_gene_name <- Cvir_reference_transcripts_exons_separated %>% separate(Dbxref_GeneID, c("Dbxref", "Genbank"), ",", extra="merge")

# separate "gene=" from the geneLOC column in order to merge with merged file
Cvir_reference_transcripts_exons_separated_gene_name_LOC <-Cvir_reference_transcripts_exons_separated_gene_name %>% separate(geneLOC, c("gene=", "LOC"), "=")

# clean up white space in the data
Cvir_reference_transcripts_exons_separated_clean <- data.frame(lapply(Cvir_reference_transcripts_exons_separated_gene_name_LOC, trimws),stringsAsFactors = FALSE)
```

### Step 2: Load the Stringtie -merged GTF file from bluewaves 

```
# Read in stringtie merged file
C_vir_stringtie <- read.csv(file= "C_Vir_stringtie_merged_RERUN.gtf", sep="\t", header = FALSE)

# set colnames for the attributes
colnames(C_vir_stringtie) <- c('seqname', 'source', 'feature','start','end','score','strand','frame', 'attribute')

# subset for transcript and exons
C_vir_stringtie_transcripts_exons <- C_vir_stringtie[C_vir_stringtie$feature %in% c("transcript","exon"),]

# grep out only the lines with a "LOC" value
C_vir_stringtie_transcripts_LOC <- C_vir_stringtie_transcripts_exons[grepl("LOC", C_vir_stringtie_transcripts_exons$attribute),]

# make into df
C_vir_stringtie_transcripts_LOC <- as.data.frame(C_vir_stringtie_transcripts_LOC)

# separate attributes column with the ; separator in tidyr, separate only the transcripts first and then the exons, and then merge together
# do this because the transcripts have a different number of attributes and the columns get off kilter if you don't separate them separately
C_vir_stringtie_transcripts_LOC_separated_2 <- C_vir_stringtie_transcripts_LOC[C_vir_stringtie_transcripts_LOC$feature %in% "transcript",] %>% separate(attribute, c("gene_id","transcript_id","gene_name","ref_gene_id"), ";")   
Cvir_stringtie_exons_LOC_separated_2 <- C_vir_stringtie_transcripts_LOC[C_vir_stringtie_transcripts_LOC$feature %in% "exon",] %>% separate(attribute, c("gene_id","transcript_id","exon_number","gene_name", "ref_gene_id"), ";")   

# merge the two dfs using smartbind, which will know to put NAs into the exon_number column for the transcripts
Cvir_stringtie_combined <- smartbind(C_vir_stringtie_transcripts_LOC_separated_2, Cvir_stringtie_exons_LOC_separated_2)

# clean up white space in the data
Cvir_stringtie_combined_clean <- data.frame(lapply(Cvir_stringtie_combined, trimws),stringsAsFactors = FALSE)

# separate gene name column and exon number column in both into separate columns with similar names so they can be merged 
Cvir_stringtie_combined_LOC <- Cvir_stringtie_combined_clean %>% separate(gene_name, into= c("gene", "LOC"), sep = " ")
```

### Step 3: Find the entries in the Stringtie file that match to the reference GFF file 

```
# Merge columns from gtf with columns from GFF3 based on match in the "LOC" column
# make sure the entries being merged really are the same by merging with the LOC, feature, start, and end columns
Cvir_stringtie_merged_GFF_combined <- merge(Cvir_reference_transcripts_exons_separated_clean, Cvir_stringtie_combined_LOC, by=c("LOC","feature","start", "end"))

#split the genbank column apart to get the XM info alone
Cvir_stringtie_merged_GFF_combined_split <- Cvir_stringtie_merged_GFF_combined %>% separate(Genbank, c("gen", "Acc_ID"), sep=":", extra="merge")

#remove "gen" column 
Cvir_stringtie_merged_GFF_combined_split <- Cvir_stringtie_merged_GFF_combined_split[,-13]
```

### Step 4: LOAD in Interproscan files (from IP5 out) and Merge BLAST2GO information from Interproscan file

```
# 1. remove all lines starting with # in bash using `for i in *.gff3; do sed '/^#/ d' $i > edited.$i`; done in bash
Cvir_Interproscan1 <- read.csv("/Users/erinroberts/Documents/PhD_Research/GOMEZCHIARRI2_FILES_DAILYNOTES/GENOME_TRANSCRIPTOME_SEQUENCES/IP5_out/protein_id/edited.CV_prot_id_aa1.fsa.gff3", sep="\t", header=FALSE)
Cvir_Interproscan2 <- read.csv("/Users/erinroberts/Documents/PhD_Research/GOMEZCHIARRI2_FILES_DAILYNOTES/GENOME_TRANSCRIPTOME_SEQUENCES/IP5_out/protein_id/edited.CV_prot_id_aa2.fsa.gff3", sep="\t", header=FALSE)
Cvir_Interproscan3 <- read.csv("/Users/erinroberts/Documents/PhD_Research/GOMEZCHIARRI2_FILES_DAILYNOTES/GENOME_TRANSCRIPTOME_SEQUENCES/IP5_out/protein_id/edited.CV_prot_id_aa3.fsa.gff3", sep="\t", header=FALSE)
Cvir_Interproscan4 <- read.csv("/Users/erinroberts/Documents/PhD_Research/GOMEZCHIARRI2_FILES_DAILYNOTES/GENOME_TRANSCRIPTOME_SEQUENCES/IP5_out/protein_id/edited.CV_prot_id_aa4.fsa.gff3", sep="\t", header=FALSE)
Cvir_Interproscan5 <- read.csv("/Users/erinroberts/Documents/PhD_Research/GOMEZCHIARRI2_FILES_DAILYNOTES/GENOME_TRANSCRIPTOME_SEQUENCES/IP5_out/protein_id/edited.CV_prot_id_aa5.fsa.gff3", sep="\t", header=FALSE)
Cvir_Interproscan6 <- read.csv("/Users/erinroberts/Documents/PhD_Research/GOMEZCHIARRI2_FILES_DAILYNOTES/GENOME_TRANSCRIPTOME_SEQUENCES/IP5_out/protein_id/edited.CV_prot_id_aa6.fsa.gff3", sep="\t", header=FALSE)

# Combine all the files (they were split up to decrease processing time)
Cvir_Interproscan_combined <- rbind(Cvir_Interproscan1, Cvir_Interproscan2, Cvir_Interproscan3, Cvir_Interproscan4, Cvir_Interproscan5, Cvir_Interproscan6)

# tell R what the column names should be, separate attributes column
colnames(Cvir_Interproscan_combined) <- c('Acc_ID', 'source', 'feature','start','end','score','strand','frame', 'attribute')

#remove lines with no match, these are called "polypeptide" in the feature column
Cvir_Interproscan_combined_only_match <- Cvir_Interproscan_combined[!grepl("polypeptide", Cvir_Interproscan_combined$feature),]

#extract only thos lines with mapped Ontology term
Cvir_Interproscan_combined_only_match_GO <- Cvir_Interproscan_combined_only_match[grepl("Ontology_term", Cvir_Interproscan_combined_only_match$attribute),]

#separate attributes column 
Cvir_Interproscan_combined_separated <- Cvir_Interproscan_combined_only_match_GO %>% separate(attribute, c("date", "target", "Ontology_term", "ID","signature_desc", "Name", "Dbxref"), sep=";", extra="merge")

# Separate and remove "ORF" from the seqname column using the function str_sub
Cvir_Interproscan_combined_separated_no_ORF <- str_sub(Cvir_Interproscan_combined_separated$Acc_ID, start=1, end = -5)

#Change all the XPs to XM so that they will merge 
substr(Cvir_Interproscan_combined_separated_no_ORF, 2, 2) <- "M"

# Switch out old Acc_ID colummn for the newly edited colum 
Cvir_Interproscan_combined_separated$Acc_ID <- Cvir_Interproscan_combined_separated_no_ORF

# clean up white space in the data
Cvir_Interproscan_combined_separated_clean <- data.frame(lapply(Cvir_Interproscan_combined_separated, trimws),stringsAsFactors = FALSE)

# Merge Acc_ID column from InterProScan file with the merged GFF and GTF file 
index <- match(Cvir_stringtie_merged_GFF_combined_split$Acc_ID, Cvir_Interproscan_combined_separated_clean$Acc_ID)
Cvir_stringtie_merged_GFF_combined_split[na.omit(index),]
