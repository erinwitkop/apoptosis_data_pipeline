## Script to match GO terms to XM IDs for C. virginica
# Erin Roberts, PhD Candidate University of Rhode Island
# Marta Gomez-Chiarri Lab
# 2019

## Load libraries
library(Repitools)
library(reshape2)
library(tidyverse)


#### Import Genome Annotation and add GO terms ####

# Import gff file with rtracklayer
C_vir_rtracklayer <- import("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/C_VIRGINICA_PIPELINE/ref_C_virginica-3.0_top_level.gff3")
C_vir_rtracklayer <- as.data.frame(C_vir_rtracklayer)

### CODE BELOW IS HOW ORIGINAL MATCHING TO GO TERMS WAS COMPLETED

# Isolate transcript lines in GFF annotation so that I can use Batch Entrez lookup for their parent proteins
C_vir_rtracklayer_transcripts <- filter(C_vir_rtracklayer, grepl("XM", transcript_id))
C_vir_rtracklayer_transcripts_unique <- subset(C_vir_rtracklayer_transcripts, !duplicated(transcript_id))
write.table(file="./Dermo_2015_Analysis/OUTPUT/C_vir_unique_transcripts.txt", C_vir_rtracklayer_transcripts_unique$transcript_id) 
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
C_vir_XM_with_XP <- import("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/C_VIRGINICA_PIPELINE/Dermo_2015_Analysis/OUTPUT/combined_batch_lookup_cds.gff")
C_vir_XM_with_XP <- as.data.frame(C_vir_XM_with_XP)
C_vir_XM_with_XP_index <- C_vir_XM_with_XP[,c("seqnames","protein_id")]
colnames(C_vir_XM_with_XP_index)[1] <- "transcript_id"

# Load in Interproscan GO annotation from LSU Kevin using Rtracklayer
# cat all the edited header removed files from Interproscan Kevin LSU files Import gff file with rtracklayer
# for i in edited*.gff3; do cat $i >> combined_CV_prot_id_interproscan.gff3 ; done

# In terminal, combine all the gff files and remove the comment lines from the script 
Cvir_Interproscan1 <- import("/Users/erinroberts/Documents/PhD_Research/GOMEZCHIARRI2_FILES_DAILYNOTES/GENOME_TRANSCRIPTOME_SEQUENCES/IP5_out/protein_id/combined_CV_prot_id_interproscan.gff3")
class(Cvir_Interproscan1)
# convert to dataframe object using annoGR2DF (repitools)
Cvir_Interproscan_DF <- annoGR2DF(Cvir_Interproscan1)
class(Cvir_Interproscan_DF) #data.frame

# Interproscan file has multiple lines per entry, need to subset for the lines that do have a GO entry because not all of them do, then subset for unique protein lines
class(Cvir_Interproscan_DF$Ontology_term)
Cvir_Interproscan_DF$Ontology_term <- as.character(Cvir_Interproscan_DF$Ontology_term)
C_vir_Interproscan_GO <- Cvir_Interproscan_DF %>% filter(Ontology_term !="character(0)")
head(C_vir_Interproscan_GO)

# keep lines with unique protein and GO terms 
C_vir_Interproscan_GO_unique <- C_vir_Interproscan_GO[!duplicated(C_vir_Interproscan_GO[,c("chr","Ontology_term")]),]

# Format GO term column correctly 
C_vir_Interproscan_GO_unique$Ontology_term <- gsub("[^[:alnum:][:blank:]+:/\\-]", "", C_vir_Interproscan_GO_unique$Ontology_term )
C_vir_Interproscan_GO_unique$Ontology_term <- gsub(C_vir_Interproscan_GO_unique$Ontology_term, pattern="c\\", replacement="", fixed=TRUE)
C_vir_Interproscan_GO_unique$Ontology_term <- gsub(C_vir_Interproscan_GO_unique$Ontology_term, pattern="\\", replacement="", fixed=TRUE)
C_vir_Interproscan_GO_unique$Ontology_term <- gsub(C_vir_Interproscan_GO_unique$Ontology_term, pattern=" ",replacement =",", fixed=TRUE)
class(C_vir_Interproscan_GO_unique)

# merge GO terms with the same protein name so there aren't multiple lines for a transcript
C_vir_Interproscan_GO_unique <- ddply(C_vir_Interproscan_GO_unique, "chr", summarize, Combined_ontology = toString(Ontology_term))

# Remove duplicate GO strings in the same column and create new column
C_vir_Interproscan_GO_unique$Combined_ontology <- gsub(" ", "", C_vir_Interproscan_GO_unique$Combined_ontology)
C_vir_Interproscan_GO_unique <- mutate(C_vir_Interproscan_GO_unique, unique_go = map(str_split(Combined_ontology, ","), unique))

# make new column that removes the "_ORF" from the end of the seqnames columns
C_vir_Interproscan_GO_unique$protein_id <- str_remove(C_vir_Interproscan_GO_unique$chr, "_ORF")

# merge XM and XP list with Interproscan list
C_vir_Interproscan_GO_unique_XM_merged <- C_vir_Interproscan_GO_unique %>% left_join(select(C_vir_XM_with_XP_index, "transcript_id","protein_id"), by = "protein_id")

# Merge Interproscan list onto transcript annotation using left_join of the TRANSCRIPT ID so that it joins correctly at the level of transcripts
# non unique file
C_vir_rtracklayer_GO <-  C_vir_rtracklayer %>% left_join(select(C_vir_Interproscan_GO_unique_XM_merged,"unique_go","transcript_id"), by = "transcript_id")

#unique file with one line per transcript
C_vir_rtracklayer_transcripts_GO <- C_vir_rtracklayer_transcripts_unique %>% left_join(select(C_vir_Interproscan_GO_unique_XM_merged,"unique_go","transcript_id"), by = "transcript_id")

# are they all NULL?
C_vir_rtracklayer_GO %>% filter(unique_go !="NULL") # nope they are not!!

# Export for use in cluster later
class(C_vir_rtracklayer_transcripts_GO)
# coerce data frame to be all character
# First coerce the data.frame to all-character
C_vir_rtracklayer_transcripts_GO_2 = data.frame(lapply(C_vir_rtracklayer_transcripts_GO, as.character), stringsAsFactors=FALSE)
head(C_vir_rtracklayer_transcripts_GO_2)
# now write to csv
write.csv(C_vir_rtracklayer_transcripts_GO_2, file="C_vir_rtracklayer_transcripts_GO.csv")
                                            