## Formatting Perkinsus Transcriptomes for input into DESeq2 ##

# Aug. 16th 2019

# Load packages
library(tidyverse)
library(readtext)
library(RSQLite)

setwd("/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/RAW DATA AND INFO/2015 Dermo Challenge /2015_transcriptomes_CLC_PMG_mapping/")


# #import the individual sample count data and rename column with sample name
# List of DP samples to add 
DP22 <- read.csv("DP22_GTGGCC.R1_GTGGCC (paired) trimmed (paired) (TE).csv")
DP23 <- read.csv("DP23_TTAGGC.R1_TTAGGC (paired) trimmed (paired) (TE).csv")
DP25<- read.csv("DP25_GTTTCG.R1_GTTTCG (paired) trimmed (paired) (TE).csv")
DP26<- read.csv("DP26_GTGGCC.R1_GTGGCC (paired) trimmed (paired) (TE).csv")
DP29<- read.csv("DP29_GTTTCG.R1_GTTTCG (paired) trimmed (paired) (TE).csv")
DP30<- read.csv("DP30_GTGGCC.R1_GTGGCC (paired) trimmed (paired) (TE).csv")
DP31<- read.csv("DP31_TTAGGC.R1_TTAGGC (paired) trimmed (paired) (TE).csv")
DP32<- read.csv("DP32_GTTTCG.R1_GTTTCG (paired) trimmed (paired) (TE).csv")
DP45<- read.csv("DP45_TAGCTT.R1_TAGCTT (paired) trimmed (paired) (TE).csv")
DP46<- read.csv("DP46_ACTGAT.R1_ACTGAT (paired) trimmed (paired) (TE).csv")
DP49<- read.csv("DP49_ATTCCT.R1_ATTCCT (paired) trimmed (paired) (TE).csv")
DP50<- read.csv("DP50_ACTGAT.R1_ACTGAT (paired) trimmed (paired) (TE).csv")
DP51<- read.csv("DP51_GGCTAC.R1_GGCTAC (paired) trimmed (paired) (TE).csv")
DP52<- read.csv("DP52_ATTCCT.R1_ATTCCT (paired) trimmed (paired) (TE).csv")
DP54<- read.csv("DP54_ACTGAT.R1_ACTGAT (paired) trimmed (paired) (TE).csv")
DP55<- read.csv("DP55_GGCTAC.R1_GGCTAC (paired) trimmed (paired) (TE).csv")
DP56<- read.csv("DP56_ATTCCT.R1_ATTCCT (paired) trimmed (paired) (TE).csv")
DP58<- read.csv("DP58_AGTCAA.R1_AGTCAA (paired) trimmed (paired) (TE).csv")
DP60<- read.csv("DP60_AGTTCC.R1_AGTTCC (paired) trimmed (paired) (TE).csv")
DP63<- read.csv("DP63_TGACCA.R1_TGACCA (paired) trimmed (paired) (TE).csv")
DP64<- read.csv("DP64_AGTTCC.R1_AGTTCC (paired) trimmed (paired) (TE).csv")
DP66<- read.csv("DP66_AGTCAA.R1_AGTCAA (paired) trimmed (paired) (TE).csv")
DP67<- read.csv("DP67_TGACCA.R1_TGACCA (paired) trimmed (paired) (TE).csv")
DP83<- read.csv("DP83_CTTGTA.R1_CTTGTA (paired) trimmed (paired) (TE).csv")
DP84<- read.csv("DP84_GTGAAA.R1_GTGAAA (paired) trimmed (paired) (TE).csv")
DP86<- read.csv("DP86_GTCCGC.R1_GTCCGC (paired) trimmed (paired) (TE).csv")
DP87<- read.csv("DP87_CTTGTA.R1_CTTGTA (paired) trimmed (paired) (TE).csv")
DP90<- read.csv("DP90_GTCCGC.R1_GTCCGC (paired) trimmed (paired) (TE).csv")
DP91<- read.csv("DP91_CTTGTA.R1_CTTGTA (paired) trimmed (paired) (TE).csv")
DP59<- read.csv("DP59_TGACCA.R1_TGACCA (paired) trimmed (paired) (TE).csv")
DP61<- read.csv("DP61_CGATGT.R1_CGATGT (paired) trimmed (paired) (TE).csv")
DP65<- read.csv("DP65_CGATGT.R1_CGATGT (paired) trimmed (paired) (TE).csv")
DP81<- read.csv("DP81_CAGATC.R1_CAGATC (paired) trimmed (paired) (TE).csv")
DP89<- read.csv("DP89_GTGAAA.R1_GTGAAA (paired) trimmed (paired) (TE).csv")
DP92<- read.csv("DP92_GTGAAA.R1_GTGAAA (paired) trimmed (paired) (TE).csv")
DA159 <- read.csv("DA159_AGTCAA.R1 (paired) trimmed (paired) (TE).csv")
DA170 <- read.csv("DA170_AGTCAA.R1 (paired) trimmed (paired) (TE).csv")
DA176 <- read.csv("DA176_CGATGT.R1 (paired) trimmed (paired) (TE).csv")
DA179 <- read.csv("DA179_AGTTCC.R1 (paired) trimmed (paired) (TE).csv")
DA180 <- read.csv("DA180_AGTCAA.R1 (paired) trimmed (paired) (TE).csv")
DA53 <- read.csv("DA53_AGTCAA.R1 (paired) trimmed (paired) (TE).csv")
DA56 <- read.csv("DA56_CGATGT.R1 (paired) trimmed (paired) (TE).csv")
DA59 <- read.csv("DA59_CGATGT.R1 (paired) trimmed (paired) (TE).csv")
DA62 <- read.csv("DA62_TGACCA.R1 (paired) trimmed (paired) (TE).csv")
DA69 <- read.csv("DA69_CGATGT.R1 (paired) trimmed (paired) (TE).csv")
LB165 <- read.csv("LB165_ACAGTG.R1 (paired) trimmed (paired) (TE).csv")
LB173 <- read.csv("LB173_ATGTCA.R1 (paired) trimmed (paired) (TE).csv")
LB179 <- read.csv("LB179_ACAGTG.R1 (paired) trimmed (paired) (TE).csv")
LB180 <- read.csv("LB180_GCCAAT.R1 (paired) trimmed (paired) (TE).csv")
LB46 <- read.csv("LB46_ACAGTG.R1 (paired) trimmed (paired) (TE).csv")
LB62 <- read.csv("LB62_GCCAAT.R1 (paired) trimmed (paired) (TE).csv")
LB67 <- read.csv("LB67_ATGTCA.R1 (paired) trimmed (paired) (TE).csv")
LB69 <- read.csv("LB69_AGTTCC.R1 (paired) trimmed (paired) (TE).csv")

# Make df list in the same order as the coldata metadata

df_list <- list(DP46 , DP50 , DP56,  DP22,  DP26 , DP32 , DP83 , 
                DP87 , DP90 , DP60 , DP64 , DP66 , DP45 , DP51 , DP54 , DA176 ,DA180 ,DP23 ,
DP30 , DA56 , DA59,  DA69,  DP49 , DP52 , DP55 , DA159, DA170, DA179, DP25 , DP29 , DP31 , DA53 , DA62, DP84  ,DP86  ,DP91 ,
LB173, LB179, DP58,  DP63,  DP67 , LB46 , LB69 , DP81 , DP89 , DP92 , LB165, LB180, DP59 , DP61 , DP65, LB62  ,LB67)

# Subset rows with a non empty Transcript.ID, select Expression.value and Transcript.ID columns 
format <- lapply(df_list, function(df) {
  df <- subset(df, Transcript.ID !="") # select only rows that have a transcript ID
  df <- df[,c("Transcript.ID", "Expression.value")] #select only the Expression value column
  })
# change names of the saved function, MUST BE IN THE SAME ORDER AS DF LIST 
names(format) <- c("DP46_1",  "DP50_1",  "DP56_1" , "DP22_1",  "DP26_1",  "DP32_1",
                   "DP83_1",  "DP87_1",  "DP90_1" , "DP60_1" , "DP64_1",  "DP66_1",  
                   "DP45_1" , "DP51_1" , "DP54_1" , "DA176_1" ,"DA180_1","DP23_1" ,
                  "DP30_1",  "DA56_1" , "DA59_1",  "DA69_1" , "DP49_1",  "DP52_1" , "DP55_1",
                   "DA159_1", "DA170_1", "DA179_1" ,"DP25_1",  "DP29_1" , "DP31_1" , "DA53_1",  "DA62_1" , "DP84_1",  "DP86_1",
                   "DP91_1" ,"LB173_1", "LB179_1", "DP58_1",  "DP63_1",  "DP67_1",  "LB46_1" , "LB69_1",  "DP81_1",
                   "DP89_1",  "DP92_1",  "LB165_1" ,"LB180_1", "DP59_1",  "DP61_1" , "DP65_1" , "LB62_1",  "LB67_1")

# create a dataframe with the list of each name
list2env(format, envir = .GlobalEnv)

#change name of Expression value of each column 
  colnames(DP22_1)[2] <- "DP22"
  colnames(DP23_1)[2] <-"DP23"
  colnames(DP25_1)[2]<-"DP25"
  colnames(DP26_1)[2]<-"DP26"
  colnames(DP29_1)[2]<-"DP29"
  colnames(DP30_1)[2]<-"DP30"
  colnames(DP31_1)[2]<-"DP31"
  colnames(DP32_1)[2]<-"DP32"
  colnames(DP45_1)[2]<-"DP45"
  colnames(DP46_1)[2]<-"DP46"
  colnames(DP49_1)[2]<-"DP49"
  colnames(DP50_1)[2]<-"DP50"
  colnames(DP51_1)[2]<-"DP51"
  colnames(DP52_1)[2]<-"DP52"
  colnames(DP54_1)[2]<-"DP54"
  colnames(DP55_1)[2]<-"DP55"
  colnames(DP56_1)[2]<-"DP56"
  colnames(DP58_1)[2]<-"DP58"
  colnames(DP60_1)[2]<-"DP60"
  colnames(DP63_1)[2]<-"DP63"
  colnames(DP64_1)[2]<-"DP64"
  colnames(DP66_1)[2]<-"DP66"
  colnames(DP67_1)[2]<-"DP67"
  colnames(DP83_1)[2]<-"DP83"
  colnames(DP84_1)[2]<-"DP86"
  colnames(DP86_1)[2]<-"DP86"
  colnames(DP87_1)[2]<-"DP87"
  colnames(DP90_1)[2]<-"DP90"
  colnames(DP91_1)[2]<-"DP91"
  colnames(DP59_1)[2]<-"DP59"
  colnames(DP61_1)[2]<-"DP61"
  colnames(DP65_1)[2]<-"DP65"
  colnames(DP81_1)[2]<-"DP81"
  colnames(DP89_1)[2]<-"DP89"
  colnames(DP92_1)[2]<-"DP92"
  colnames(DA159_1)[2]<-"DA159"
  colnames(DA170_1)[2]<-"DA170"
  colnames(DA176_1)[2]<-"DA176"
  colnames(DA179_1)[2]<-"DA179"
  colnames(DA180_1)[2]<-"DA180"
  colnames(DA53_1)[2]<- "DA53"
  colnames(DA56_1)[2]<- "DA56"
  colnames(DA59_1)[2]<- "DA59"
  colnames(DA62_1)[2]<- "DA62"
  colnames(DA69_1)[2]<- "DA69"
  colnames(LB165_1)[2]<-"LB165"
  colnames(LB173_1)[2]<-"LB173"
  colnames(LB179_1)[2]<-"LB179"
  colnames(LB180_1)[2]<-"LB180"
  colnames(LB46_1)[2]<-"LB46"
  colnames(LB62_1)[2]<-"LB62"
  colnames(LB67_1)[2]<-"LB67"
  colnames(LB69_1)[2]<-"LB69"


# Merge all of the data frames together by Transcript.ID using the reduce function from the purr package
# Merging them all at once caused the computer to crash, split into groups of ten and try 
df_list_1 <- list(DP46_1, DP50_1 , DP56_1,  DP22_1,  DP26_1 , DP32_1 , DP83_1 , 
                  DP87_1 , DP90_1 , DP60_1)
df_list_2 <- list(DP64_1, DP66_1 , DP45_1 , DP51_1 , DP54_1 , DA176_1 ,DA180_1 ,DP23_1 ,
                   DP30_1 , DA56_1)
df_list_3 <- list(DA59_1,  DA69_1,  DP49_1 , DP52_1 , DP55_1 , DA159_1, DA170_1, DA179_1, DP25_1 , DP29_1)
df_list_4 <- list(DP31_1 , DA53_1 , DA62_1, DP84_1  ,DP86_1  ,DP91_1 ,
                   LB173_1, LB179_1, DP58_1,  DP63_1)
df_list_5 <- list(DP67_1 , LB46_1 , LB69_1 , DP81_1 , DP89_1 , DP92_1 , LB165_1, LB180_1, DP59_1 , DP61_1)
df_list_6 <- list(DP65_1, LB62_1  ,LB67_1)

merge1 <- df_list_1  %>% reduce(left_join, by = "Transcript.ID")
merge2 <- df_list_2  %>% reduce(left_join, by = "Transcript.ID")
merge3 <- df_list_3  %>% reduce(left_join, by = "Transcript.ID")
merge4 <- df_list_4  %>% reduce(left_join, by = "Transcript.ID")
merge5 <- df_list_5  %>% reduce(left_join, by = "Transcript.ID")
merge6 <- df_list_6  %>% reduce(left_join, by = "Transcript.ID")

# Merge all the lists into one
merge_list <- list(merge1, merge2, merge3, merge4, merge5, merge6)

PMG_counts_matrix <- merge_list %>% reduce(left_join, by ="Transcript.ID")
head(PMG_counts_matrix)

# Make the Transcript.ID the rownames
rownames(PMG_counts_matrix) <- PMG_counts_matrix$Transcript.ID
head(PMG_counts_matrix)
PMG_counts_matrix <- PMG_counts_matrix[,-1] # remove the Transcript.ID column
head(PMG_counts_matrix)

write.csv(PMG_counts_matrix, file="Perkinsus_mapping_TE_Counts_Matrix.csv")
