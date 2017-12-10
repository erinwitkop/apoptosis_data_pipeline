#05_2_Bac_Viral_DESeq2_OSHV1_Bac_Challenge_FUNCTIONAL_ENRICHMENT.R

#this script utilized topGO to perform Gene Set Ennrichment Analysis (GSEA)
#downloaded the Generic GO slim from the GO consortium, http://www.geneontology.org/page/go-slim-and-subset-guide

#Load packages
source("http://bioconductor.org/biocLite.R")
biocLite(c("topGO","ALL","Rgraphviz"))
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
#see what ontologies are supported

ls("package:GO.db")

#Only using those genes that already have annotated terms

####Data preparation####
#input is a table that contains the genes, the GO mappings, and the p-values
#OsHV1
head(resoshvTran_05_df_FULL)
colnames(resoshvTran_05_df_FULL)
oshv_GO_table <- resoshvTran_05_df_FULL[,c("ID","padj", "Protein.names", "Gene.ontology.IDs")]
head(oshv_GO_table)

#Bac
colnames(resBacTran_05_df_FULL)
Bac_GO_table <- resBacTran_05_df_FULL[,c("ID","padj", "Protein.names", "Gene.ontology.IDs")]
head(Bac_GO_table)

#####Prepare topGO data ####

#remove genes with no mapping
OsHV1_data_noblanks <- oshv_GO_table[!(oshv_GO_table$ Gene.ontology.IDs == ""), ] #7958 rows
Bac_data_noblanks <- Bac_GO_table[!(Bac_GO_table$ Gene.ontology.IDs == ""), ] #5518 rows


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
OsHV1_readMapping <- OsHV1_data_noblanks[,c("ID", "Gene.ontology.IDs")]
Bac_readMapping <- Bac_data_noblanks[,c("ID", "Gene.ontology.IDs")]
write.table(OsHV1_readMapping, file="OsHV1_readMapping", row.names=FALSE, sep= "\t")
write.table(Bac_readMapping, file="Bac_readMapping", row.names=FALSE, sep= "\t") 
#must be tab delimited
#reformat using bash commands (remove quotes, change semicolon to commas, remove "transcript:")
sed 's/\"//g' Bac_readMapping > Bac_readMapping_quotes.txt
sed 's/\"//g' OsHV1_readMapping > OsHV1_readMapping_quotes.txt
 sed -e 's/\;/,/g' Bac_readMapping_quotes.txt > Bac_readMapping_commas.txt
sed -e 's/\;/,/g' OsHV1_readMapping_quotes.txt > OsHV1_readMapping_commas.txt
tail -n +2 Bac_readMapping_commas.txt > Bac_readMapping_finished.txt
tail -n +2 OsHV1_readMapping_commas.txt > OsHV1_readMapping_finished.txt

  #sed 's/^.\{12\}//g' Bac_readMapping_finished.txt > Bac_readMapping_finished2.txt
 # sed 's/^.\{12\}//g' OsHV1_readMapping_finished.txt > OsHV1_readMapping_finished2.txt
 
#Use readMappings()
OsHV1geneID2GO <- readMappings(file = "OsHV1_readMapping_finished.txt")
str(head(OsHV1geneID2GO))
BacgeneID2GO <- readMappings(file= "Bac_readMapping_finished.txt")
str(head(BacgeneID2GO))

#how to use a GO slim ontology with topGO

#Defining your list of genes of interest, and the 'gene universe' to compare it to
#Define genes of interest: the number of genes with a padj less than 0.05
#OsHV1 
OsHV1geneUniverse <- names(OsHV1geneID2GO) 
OsHV1_topDiff <- OsHV1_data_noblanks[ which( OsHV1_data_noblanks$padj < 0.05 ), ]
head(OsHV1_topDiff)
OsHV1_geneofInterest <- OsHV1_topDiff[,1] #need to remove "transcript:
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
numSigGenes(OsHV1_GOdata) #70

#Bacterial challenge genes
BacgeneUniverse <- names(BacgeneID2GO) 
Bac_topDiff <- Bac_data_noblanks[which( Bac_data_noblanks$padj < 0.05), ]
Bac_geneofInterest <- Bac_topDiff[,1] 
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
OsHV1_numsignif: 7
  
#print out the top 'numsignif' results:
OsHV1_numsignif_Res <- GenTable(OsHV1_GOdata, classicFisher = OsHV1Fisher,
        elimFisher = OsHV1FisherElim, topgoFisher = OsHV1Fisher_weight,
        parentchildFisher = OsHV1FisherParentchild, orderBy = "topgoFisher", 
        ranksOf = "classicFisher", topNodes = OsHV1_numsignif)
write.csv(OsHV1_numsignif_Res, file = "OsHV1_numsignif_Res_NEW.csv")

# print a graph (to a pdf file) with the top 'numsignif' results from weight01 analysis:
printGraph(OsHV1_GOdata, OsHV1Fisher_weight, firstSigNodes = OsHV1_numsignif, useInfo = "all", pdfSW = TRUE)
#rename your graph! 

#Find top significant nodes 
OsHV1_topRes <- GenTable(OsHV1_GOdata, classicFisher = OsHV1Fisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 20)
write.csv(OsHV1_topRes, file= "OsHV1_topGO_topSignificantNodes_NEW.csv")

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
#[1] "Term GO:0043248 transcripts: EKC34293,EKC42394"
#[1] "Term GO:0006298 transcripts: EKC25057,EKC26394"
#[1] "Term GO:0007154 transcripts: EKC17376,EKC28428,EKC33131,EKC35537"
#[1] "Term GO:0015992 transcripts: EKC18664,EKC37019"
#[1] "Term GO:0007010 transcripts: EKC17965,EKC18223,EKC21714,EKC35202,EKC35203"
#[1] "Term GO:0045454 transcripts: EKC32832,EKC37672,EKC39509"
#[1] "Term GO:0006040 transcripts: EKC27126,EKC40852"


#Find gene name for these genes, listed in OsHV1_data
#Lookup table with match(), return indices of match in second arg to first arg:
OsHV1_sigGO_name_lookup <- read.csv("OsHV1_sigGO_transcriptIDs_NEW.csv", header=TRUE, sep=",")
OsHV1_sigGO_gene_names_index <- oshv_GO_table$ID %in% OsHV1_sigGO_name_lookup$ID
!is.na(match(oshv_GO_table$ID,oshv_GO_table$ID))
OsHV1_sigGO_gene_names <- oshv_GO_table[OsHV1_sigGO_gene_names_index, ]
write.csv(OsHV1_sigGO_gene_names, file="OsHV1_sigGO_gene_names_NEW.csv")

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
Bac_numsignif # = 2 
  
#print out the top 'numsignif' results:
Bac_numsignif_Res <- GenTable(Bac_GOdata, classicFisher = BacFisher,
                                  elimFisher = BacFisherElim, topgoFisher = BacFisher_weight,
                                  parentchildFisher = BacFisherParentchild, orderBy = "topgoFisher", 
                                  ranksOf = "classicFisher", topNodes = Bac_numsignif)
write.csv(Bac_numsignif_Res, file= "Bac_numsignif_Res_NEW.csv")

# print a graph (to a pdf file) with the top 'numsignif' results from weight01 analysis:
printGraph(Bac_GOdata, BacFisher_weight, firstSigNodes = Bac_numsignif, useInfo = "all", pdfSW = TRUE)
#rename your graph! 

#Find top significant nodes 
Bac_topRes <- GenTable(Bac_GOdata, classicFisher = BacFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 20)
write.csv(Bac_topRes, file="Bac_topGO_topSignificantNodes_NEW.csv")

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
#[1] "Term GO:0007264 transcripts: EKC20004,EKC21485,EKC24479,EKC26517,EKC28675,EKC39204"
#[1] "Term GO:0006518 transcripts: EKC17829,EKC18290,EKC19096,EKC25781,EKC28353,EKC35385,EKC37643,EKC39112,EKC42018"


#Find gene name for these genes, listed in OsHV1_data
#need to remove "transcript:", add back new column
#Find gene name for these genes, listed in OsHV1_data
#Lookup table with match(), return indices of match in second arg to first arg:
Bac_sigGO_name_lookup <- read.csv("Bac_sigGO_transcriptIDs_NEW.csv", header=TRUE, sep=",")
Bac_sigGO_gene_names_index <- Bac_GO_table$ID %in% Bac_sigGO_name_lookup$ID
!is.na(match(Bac_GO_table$ID,Bac_GO_table$ID))
Bac_sigGO_gene_names <- Bac_GO_table[Bac_sigGO_gene_names_index, ]
write.csv(Bac_sigGO_gene_names, file="Bac_sigGO_gene_names_NEW.csv")


#### list to use for REVIGO analysis ####
OsHV1_REVIGO <- score(OsHV1Fisher_weight)
write.csv(OsHV1_REVIGO, file="OsHV1_REVIGO.csv")
Bac_REVIGO <- score(BacFisher_weight)
write.csv(Bac_REVIGO, file="Bac_REVIGO.csv")

####COMPILED RESULTS TABLES ####
Bac_numsignif_Res["Challenge"] <- "Bac"
OsHV1_numsignif_Res["Challenge"] <- "OsHV1"
COMBINED_TOPGO_RES_NUMSIGNIF <- rbind(Bac_numsignif_Res,OsHV1_numsignif_Res)
write.csv(COMBINED_TOPGO_RES_NUMSIGNIF, file="COMBINED_TOPGO_NUMSIGNIF_RES_NEW.csv")
Bac_topRes["Challenge"] <- "Bac"
OsHV1_topRes["Challenge"] <- "OsHV1"
COMBINED_TOPGO_RES_TOPSIGNIFICANTNODES <- rbind(Bac_topRes,OsHV1_topRes)
write.csv(COMBINED_TOPGO_RES_TOPSIGNIFICANTNODES, file="COMBINED_TOPGO_RES_TOPSIGNIFICANTNODES_NEW.csv")
Bac_sigGO_gene_names["Challenge"] <- "Bac"
OsHV1_sigGO_gene_names["Challenge"] <- "OsHV1"
COMBINED_TOPGO_RES_SIGGO_GENE_NAMES <- rbind(Bac_sigGO_gene_names,OsHV1_sigGO_gene_names)
COMBINED_TOPGO_RES_SIGGO_GENE_NAMES <- COMBINED_TOPGO_RES_SIGGO_GENE_NAMES %>%
  filter(Protein.names !="Uncharacterized protein")
write.csv(COMBINED_TOPGO_RES_SIGGO_GENE_NAMES, file="COMBINED_TOPGO_RES_SIGGO_GENE_NAMES_NEW.csv")

#graph results with topGO

#References: Gene set enrichment analysis with topGO
#Adrian Alexa, Jorg Rahnenfuhrer, April 24, 2017, http://www.mpi-sb.mpg.de/alexa
#http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html
