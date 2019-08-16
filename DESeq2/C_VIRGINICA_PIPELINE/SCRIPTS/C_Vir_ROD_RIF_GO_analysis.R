#C_Vir_ROD_RIF_GO_analysis

#Load packages
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("topGO","ALL","Rgraphviz"))
library(topGO)
#install.packages("tm")
library(tm)
library(ALL)
library(GO.db)
library(tidyr)
library(dplyr)
library(plyr)
library(genefilter)
library(Rgraphviz)
#install.packages("org.Hs.egGO")
library()

#Using GO slim and my subset of GO terms. Go slims give a broad overview of content without the fine grained terms
#"PANTHER GO-slim" uses a selected set of terms from the Gene Ontology TM for classifications by molecular function, 
# biological process and cellular component. The PANTHER Protein Class ontology was adapted from the PANTHER/X
# molecular function ontology, and includes commonly used classes of protein functions, many of which are not 
# covered by GO molecular function. Download the classes and relationship information
#Molecular function is the function that a protein performs on its direct molecular targets; 
# e.g. the insulin receptor has transmembrane receptor tyrosine protein kinase activity (GO:0004714), 
#which means it catalyzes the reaction that adds a phosphate group to a tyrosine in another protein (its target).
#Cellular component is the location where the protein performs its molecular function; 
#e.g. the insulin receptor is located in the plasma membrane (GO:0005886).
#Biological process covers the biological systems to which a protein contributes; 
#e.g. the insulin receptor is involved in regulation of carbohydrate metabolic process (GO:0006109).


#Package to get the GO terms 
#org.Hs.egGO
#https://stuff.mit.edu/afs/athena/software/r/current/lib/R/library/org.Hs.eg.db/html/org.Hs.egGO.html

#example code
#Load packages


####Data preparation####
#input is a table that contains the genes, the GO mappings, and the p-values
#OsHV1
head(resoshv1Tran_05_df)
OsHV1_topGOsubset <- 
  resoshv1Tran_05_df[grep("transcript:", rownames(resoshv1Tran_05_df)), ]
head(OsHV1_topGOsubset)
nrow(OsHV1_topGOsubset) #rows = 12112
#make a new column with the rownames, called Ensembl_genomes_transcript
OsHV1_topGOsubset$Ensembl_Genomes_Transcript<-rownames(OsHV1_topGOsubset)  
##needs to be made a data frame so that can be joined later by plyr
OsHV1_topGOsubset_df <- as.data.frame(OsHV1_topGOsubset)  
OsHV1_topGOtranscript_subset <- as.data.frame(rownames(OsHV1_topGOsubset))
head(OsHV1_topGOtranscript_subset)
colnames(OsHV1_topGOtranscript_subset) <- "transcript"
OsHV1transcriptIDdf <- OsHV1_topGOtranscript_subset %>%
  separate(transcript, c("transcript", "EKC"), ":")
#write to string
OsHV1_topGO_transcriptIDstring <- toString(OsHV1transcriptIDdf[,2], sep=',')
write(OsHV1_topGO_transcriptIDstring, "OsHV1_topGO_transcriptIDstring", sep = ",")
#cut off commas with sed 's/,//g' INFILE > outfile

#Load annotation file
OsHV1_topGO_annotations <- read.csv("OsHV1_topGO_GOIDS.csv", sep=",", header=TRUE) #12092 rows
nrow(OsHV1_topGO_annotations)
head(OsHV1_topGO_annotations)
#add "transcript: before every string" in Ensembl Genomes Transcript
OsHV1_topGO_annotations$Ensembl_Genomes_Transcript = paste('transcript', OsHV1_topGO_annotations$Ensembl_Genomes_Transcript, sep=':')
#concatentate gene expression and uniprot table
#Perform join by plyr
OsHV1_topGO_annotations
OsHV1_topGO_annotations_FULL <- join(OsHV1_topGO_annotations, OsHV1_topGOsubset_df, by="Ensembl_Genomes_Transcript")
nrow(OsHV1_topGO_annotations_FULL) #12092 This is what I wanted 

#Bac
head(resBacTran_05_df)
BAC_topGOsubset <- 
  resBacTran_05_df[grep("transcript:", rownames(resBacTran_05_df)), ]
head(BAC_topGOsubset)
#make a new column with the rownames, called Ensembl_genomes_transcript
BAC_topGOsubset$Ensembl_Genomes_Transcript<-rownames(BAC_topGOsubset)  
##needs to be made a data frame so that can be joined later by plyr
BAC_topGOsubset_df <- as.data.frame(BAC_topGOsubset)  #use this for the plyr join!
BAC_topGOtranscript_subset <- as.data.frame(rownames(BAC_topGOsubset))
head(BAC_topGOtranscript_subset)
colnames(BAC_topGOtranscript_subset) <- "transcript"
BACtranscriptIDdf <- BAC_topGOtranscript_subset %>%
  separate(transcript, c("transcript", "EKC"), ":") #edited to make rownames
nrow(BACtranscriptIDdf) #8075
#write to string so I can do lookup
BAC_topGO_transcriptIDstring <- toString(BACtranscriptIDdf[,2], sep=',')
write(BAC_topGO_transcriptIDstring, "BAC_topGO_transcriptIDstring", sep = ",")

#Load annotation file
BAC_topGO_annotations <- read.csv("Bac_topGO_transcript_GOIDs.csv", sep=",", header=TRUE) #include GO IDs
#8065 out of 8075 were mapped to 8092
nrow(BAC_topGO_annotations) #8092
head(BAC_topGO_annotations)
#add "transcript: before every string" in Ensembl Genomes Transcript
BAC_topGO_annotations$Ensembl_Genomes_Transcript = paste('transcript', BAC_topGO_annotations$Ensembl_Genomes_Transcript, sep=':')
#concatentate gene expression and uniprot table
#Perform join by plyr
BAC_topGO_annotations_FULL <- join(BAC_topGO_annotations, BAC_topGOsubset_df, by="Ensembl_Genomes_Transcript")
nrow(BAC_topGO_annotations_FULL) #8092


#####Prepare topGO data ####
#filtering step using genefilter?

#subset for correct rows
OsHV1_data <- OsHV1_topGO_annotations_FULL[,c("Ensembl_Genomes_Transcript", "Gene_names", "Protein_names","pvalue","padj", "Gene_ontology_IDs")]
nrow(OsHV1_data) #12092
Bac_data <- BAC_topGO_annotations_FULL[,c("Ensembl_Genomes_Transcript", "Gene_names", "Protein_names","pvalue","padj", "Gene_ontology_IDs")]
nrow(Bac_data) #8092
#remove genes with no mapping
OsHV1_data_noblanks <- OsHV1_data[!(OsHV1_data$Gene_ontology_IDs == ""), ] #7958 rows
Bac_data_noblanks <- Bac_data[!(Bac_data$Gene_ontology_IDs == ""),] #5518 rows

#Create topGO object
#This object will contain all information necessary for the GO analysis, 
# namely the list of genes, the list of interesting genes, the gene
#scores (if available) and the part of the GO ontology (the GO graph) 
#which needs to be used in the analysis. The topGOdata object will be the
#input of the testing procedures, the evaluation and visualisation functions.
#Can input custom annotations already found, need a gene to GO mapping approach
#parse the annotations file with function readMappings()
#file format required by readMappings = gene_ID<TAB>GO_ID1, GO_ID2, GO_ID3,

#format GO mapping files to readMapping format
OsHV1_readMapping <- OsHV1_data_noblanks[,c("Ensembl_Genomes_Transcript", "Gene_ontology_IDs")]
Bac_readMapping <- Bac_data_noblanks[,c("Ensembl_Genomes_Transcript", "Gene_ontology_IDs")]
write.table(OsHV1_readMapping, file="OsHV1_readMapping", row.names=FALSE, sep= "\t")
write.table(Bac_readMapping, file="Bac_readMapping", row.names=FALSE, sep= "\t") 
#must be tab delimited
#reformat using bash commands (remove quotes, change semicolon to commas, remove "transcript:")
# sed -e 's/\;/,/g' Bac_readMapping > Bac_readMapping_commas
#sed 's/\"//g' Bac_readMapping_commas > Bac_readMapping_finished.txt
# sed 's/\"//g' OsHv1_readMapping_commas > OsHV1_readMapping_finished.txt
#sed 's/^.\{11\}//g' Bac_readMapping_finished.txt > Bac_readMapping_finished2.txt
#remove first line: tail -n +2 Bac_readMapping_finished2.txt > Bac_readMapping_finished3.txt

#Use readMappings()
OsHV1geneID2GO <- readMappings(file = "OsHV1_readMapping_finished3.txt")
str(head(OsHV1geneID2GO))
BacgeneID2GO <- readMappings(file= "Bac_readMapping_finished3.txt")
str(head(BacgeneID2GO))

#Defining your list of genes of interest, and the 'gene universe' to compare it to
#Define genes of interest: the number of genes with a padj less than 0.05
#OsHV1 
OsHV1geneUniverse <- names(OsHV1geneID2GO) 
OsHV1_topDiff <- OsHV1_data_noblanks[ which( OsHV1_data_noblanks$padj < 0.05 ), ]
OsHV1_geneofInterest <- OsHV1_topDiff[,1] #need to remove "transcript:
stopwords="transcript:"
OsHV1_geneofInterest <- removeWords(OsHV1_geneofInterest,stopwords)
#tell TopGO where these interesting genes appear in the 'geneUniverse' vector:
#geneList tells topGO which genes are of interest
OsHV1geneList <- factor(as.integer(OsHV1geneUniverse %in% OsHV1_geneofInterest))
names(OsHV1geneList) <- OsHV1geneUniverse

#Build topGO object for biological process
OsHV1_GOdata <- new("topGOdata", description = "OsHV1 Gene Enrichment", ontology = "BP",
                    allGenes = OsHV1geneList,
                    nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = OsHV1geneID2GO)
OsHV1_GOdata
#I care about BP, biological process
#nodeSize=used to prune the GO hierarchy from the terms which have less than 5 annotated genes
#annFUN.gene2GO = this function is used when the annotations are provided as a gene-to-GOs mapping.

#The list of genes of interest can be accessed using the method sigGenes():
#the specified genes of interest are the ones called significant
OsHV1_sg <- sigGenes(OsHV1_GOdata)
str(OsHV1_sg)
numSigGenes(OsHV1_GOdata) 

#Bacterial challenge genes
BacgeneUniverse <- names(BacgeneID2GO) 
Bac_topDiff <- Bac_data_noblanks[which( Bac_data_noblanks$padj < 0.05), ]
Bac_geneofInterest <- Bac_topDiff[,1] #need to remove "transcript:
stopwords="transcript:"
Bac_geneofInterest <- removeWords(Bac_geneofInterest,stopwords)
#tell TopGO where these interesting genes appear in the 'geneUniverse' vector:
#geneList tells topGO which genes are of interest
BacgeneList <- factor(as.integer(BacgeneUniverse %in% Bac_geneofInterest))
names(BacgeneList) <- BacgeneUniverse

#Build topGO object for biological process
Bac_GOdata <- new("topGOdata", description = "Bac Gene Enrichment", ontology = "BP",
                  allGenes = BacgeneList,
                  nodeSize = 5,  annot = annFUN.gene2GO, gene2GO = BacgeneID2GO)


#### OsHV1 Enrichment Analysis ####
#Fisher's Exact test
#Run all (though weight01 us default)
#to take GO hierarchy into account use algorithm='weight01', classic doesn't take into account
OsHV1Fisher <- runTest(OsHV1_GOdata, algorithm = "classic", statistic = "fisher")
OsHV1Fisher_weight <- runTest(OsHV1_GOdata, algorithm="weight01", statistic="fisher")
OsHV1FisherElim <- runTest(OsHV1_GOdata, algorithm="elim", statistic="fisher")
OsHV1FisherParentchild <- runTest(OsHV1_GOdata, algorithm="parentchild", statistic="fisher")

# see how many results we get where weight01 gives a P-value <= 0.05:
OsHV1summary <- summary(attributes(OsHV1Fisher_weight)$score <= 0.05)
OsHV1_numsignif <- as.integer(OsHV1summary[[3]]) # how many terms is it true that P <= 0.05
OsHV1_numsignif: 6!
  
  #print out the top 'numsignif' results:
  OsHV1_numsignif_Res <- GenTable(OsHV1_GOdata, classicFisher = OsHV1Fisher,
                                  elimFisher = OsHV1FisherElim, topgoFisher = OsHV1Fisher_weight,
                                  parentchildFisher = OsHV1FisherParentchild, orderBy = "topgoFisher", 
                                  ranksOf = "classicFisher", topNodes = OsHV1_numsignif)
write.csv(OsHV1_numsignif_Res, file = "OsHV1_numsignif_Res.csv")

# print a graph (to a pdf file) with the top 'numsignif' results from weight01 analysis:
printGraph(OsHV1_GOdata, OsHV1Fisher_weight, firstSigNodes = OsHV1_numsignif, useInfo = "all", pdfSW = TRUE)
#rename your graph! 

#Find top significant nodes 
OsHV1_topRes <- GenTable(OsHV1_GOdata, classicFisher = OsHV1Fisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 20)
write.csv(OsHV1_topRes, file= "OsHV1_topGO_topSignificantNodes.csv")

# print out the genes that are annotated with the significantly enriched GO terms:
OsHV1_terms <- OsHV1_numsignif_Res$GO.ID
OsHV1_genes <- genesInTerm(OsHV1_GOdata, OsHV1_terms)
for (i in 1:length(OsHV1_terms))
{
  OsHV1_GOterm <- OsHV1_terms[i]
  OsHV1genesforterm <- OsHV1_genes[OsHV1_GOterm][[1]] 
  OsHV1_factor <- OsHV1genesforterm %in% OsHV1_geneofInterest # find the genes that are in the list of genes of interest
  OsHV1genesforterm2 <- OsHV1genesforterm[OsHV1_factor == TRUE] 
  OsHV1genesforterm2 <- paste(OsHV1genesforterm2, collapse=',')
  print(paste("Term",OsHV1_GOterm,"transcripts:",OsHV1genesforterm2))
}

#Find gene name for these genes, listed in OsHV1_data
#need to remove "transcript:", add back new column
OsHV1_data_edited1 <- removeWords(OsHV1_data$Ensembl_Genomes_Transcript, stopwords)
OsHV1_data_edited2<- cbind(OsHV1_data, OsHV1_data_edited1)
#Lookup table with match(), return indices of match in second arg to first arg:
OsHV1_sigGO_name_lookup <- read.csv("OsHV1_sigGO_transcriptIDs.csv", header=TRUE, sep=",")
OsHV1_sigGO_gene_names_index <- OsHV1_data_edited2$OsHV1_data_edited1 %in% OsHV1_sigGO_name_lookup$OsHV1_data_edited1
!is.na(match(OsHV1_data_edited2$OsHV1_data_edited1,OsHV1_sigGO_name_lookup$OsHV1_data_edited1))
OsHV1_sigGO_gene_names <- OsHV1_data_edited2[OsHV1_sigGO_gene_names_index, ]
write.csv(OsHV1_sigGO_gene_names, file="OsHV1_sigGO_gene_names.csv")

#### Bac Enrichment Analysis ####
#Fisher's Exact test
#Run all (though weight01 us default)
#to take GO hierarchy into account use algorithm='weight01', classic doesn't take into account
BacFisher <- runTest(Bac_GOdata, algorithm = "classic", statistic = "fisher")
BacFisher_weight <- runTest(Bac_GOdata, algorithm="weight01", statistic="fisher")
BacFisherElim <- runTest(Bac_GOdata, algorithm="elim", statistic="fisher")
BacFisherParentchild <- runTest(Bac_GOdata, algorithm="parentchild", statistic="fisher")

# see how many results we get where weight01 gives a P-value <= 0.05:
Bacsummary <- summary(attributes(BacFisher_weight)$score <= 0.05)
Bac_numsignif <- as.integer(Bacsummary[[3]]) # how many terms is it true that P <= 0.05
Bac_numsignif # = 4 

#print out the top 'numsignif' results:
Bac_numsignif_Res <- GenTable(Bac_GOdata, classicFisher = BacFisher,
                              elimFisher = BacFisherElim, topgoFisher = BacFisher_weight,
                              parentchildFisher = BacFisherParentchild, orderBy = "topgoFisher", 
                              ranksOf = "classicFisher", topNodes = Bac_numsignif)
write.csv(Bac_numsignif_Res, file= "Bac_numsignif_Res.csv")

# print a graph (to a pdf file) with the top 'numsignif' results from weight01 analysis:
printGraph(Bac_GOdata, BacFisher_weight, firstSigNodes = Bac_numsignif, useInfo = "all", pdfSW = TRUE)
#no. nodes = 38
#rename your graph! 

#Find top significant nodes 
Bac_topRes <- GenTable(Bac_GOdata, classicFisher = BacFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 20)
write.csv(Bac_topRes, file="Bac_topGO_topSignificantNodes.csv")

# print out the genes that are annotated with the significantly enriched GO terms:
Bac_terms <- Bac_numsignif_Res$GO.ID
Bac_genes <- genesInTerm(Bac_GOdata, Bac_terms)
for (i in 1:length(Bac_terms))
{
  Bac_GOterm <- Bac_terms[i]
  Bacgenesforterm <- Bac_genes[Bac_GOterm][[1]] 
  Bac_factor <- Bacgenesforterm %in% Bac_geneofInterest # find the genes that are in the list of genes of interest
  Bacgenesforterm2 <- Bacgenesforterm[Bac_factor == TRUE] 
  Bacgenesforterm2 <- paste(Bacgenesforterm2, collapse=',')
  print(paste("Term",Bac_GOterm,"transcripts:",Bacgenesforterm2))
}

#Find gene name for these genes, listed in OsHV1_data
#need to remove "transcript:", add back new column
Bac_data_edited1 <- removeWords(Bac_data$Ensembl_Genomes_Transcript, stopwords)
Bac_data_edited2<- cbind(Bac_data, Bac_data_edited1)
#Lookup table with match(), return indices of match in second arg to first arg:
Bac_sigGO_name_lookup <- read.csv("Bac_sigGO_transcriptIDs.csv", header=TRUE, sep=",")
Bac_sigGO_gene_names_index <- Bac_data_edited2$Bac_data_edited1 %in% Bac_sigGO_name_lookup$Bac_data_edited1
!is.na(match(Bac_data_edited2$Bac_data_edited1,Bac_sigGO_name_lookup$Bac_data_edited1))
Bac_sigGO_gene_names <- Bac_data_edited2[Bac_sigGO_gene_names_index, ]
write.csv(Bac_sigGO_gene_names, file="Bac_sigGO_gene_names.csv")

####COMPILED RESULTS TABLES ####
Bac_numsignif_Res["Challenge"] <- "Bac"
OsHV1_numsignif_Res["Challenge"] <- "OsHV1"
COMBINED_TOPGO_RES_NUMSIGNIF <- rbind(Bac_numsignif_Res,OsHV1_numsignif_Res)
write.csv(COMBINED_TOPGO_RES_NUMSIGNIF, file="COMBINED_TOPGO_NUMSIGNIF_RES.csv")
Bac_topRes["Challenge"] <- "Bac"
OsHV1_topRes["Challenge"] <- "OsHV1"
COMBINED_TOPGO_RES_TOPSIGNIFICANTNODES <- rbind(Bac_topRes,OsHV1_topRes)
write.csv(COMBINED_TOPGO_RES_TOPSIGNIFICANTNODES, file="COMBINED_TOPGO_RES_TOPSIGNIFICANTNODES.csv")
Bac_sigGO_gene_names["Challenge"] <- "Bac"
OsHV1_sigGO_gene_names["Challenge"] <- "OsHV1"
#get rid of column with unmatching names
Bac_sigGO_gene_names <- Bac_sigGO_gene_names[,c(1,2,3,4,5,6,8)]
OsHV1_sigGO_gene_names <- OsHV1_sigGO_gene_names[,c(1,2,3,4,5,6,8)]
COMBINED_TOPGO_RES_SIGGO_GENE_NAMES <- rbind(Bac_sigGO_gene_names,OsHV1_sigGO_gene_names)
COMBINED_TOPGO_RES_SIGGO_GENE_NAMES <- COMBINED_TOPGO_RES_SIGGO_GENE_NAMES %>%
  filter(Protein_names !="Uncharacterized protein")
write.csv(COMBINED_TOPGO_RES_SIGGO_GENE_NAMES, file="COMBINED_TOPGO_RES_SIGGO_GENE_NAMES.csv")

#References: Gene set enrichment analysis with topGO
#Adrian Alexa, Jorg Rahnenfuhrer, April 24, 2017, http://www.mpi-sb.mpg.de/alexa
#http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html
