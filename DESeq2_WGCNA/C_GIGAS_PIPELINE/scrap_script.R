#Unusued/incorrect code taken from the pipeline_just_pathogen_challenges file


#LPS
LPSCountData <- gramNegCountData[ ,c(1,2)]
head(LPSCountData)
LPSColData <- gramNegColData[c(1,2),c(1,2) ]
head(LPSColData)
ddsLPS <- DESeqDataSetFromMatrix(countData = LPSCountData, 
                                 colData = LPSColData, 
                                 design = ~ condition)
head(ddsLPS) #review LPS data set

#V. aes
V_aesCountData <- gramNegCountData[ ,c(1,3)]
head(V_aesCountData)
V_aesColData <- gramNegColData[c(1,3), ]
head(V_aesColData)
ddsV_aes <- DESeqDataSetFromMatrix(countData = V_aesCountData, 
                                   colData = V_aesColData, 
                                   design = ~ condition)
head(ddsV_aes) #review V_aes data set

#V. ang
V_angCountData <- gramNegCountData[ ,c(1,4)]
head(V_angCountData)
V_angColData <- gramNegColData[c(1,4), ]
head(V_angColData)
ddsV_ang <- DESeqDataSetFromMatrix(countData = V_angCountData, 
                                   colData = V_angColData, 
                                   design = ~ condition)
head(ddsV_ang) #review V_ang data set

#V. alg1

V_alg1CountData <- gramNegCountData[ ,c(1,5)]
head(V_alg1CountData)
V_alg1ColData <- gramNegColData[c(1,5), ]
head(V_alg1ColData)
ddsV_alg1 <- DESeqDataSetFromMatrix(countData = V_alg1CountData, 
                                    colData = V_alg1ColData, 
                                    design = ~ condition)
head(ddsV_alg1) #review V_alg1 data set

#V. alg2
V_alg2CountData <- gramNegCountData[ ,c(1,6)]
head(V_alg2CountData)
V_alg2ColData <- gramNegColData[c(1,6), ]
head(V_alg2ColData)
ddsV_alg2 <- DESeqDataSetFromMatrix(countData = V_alg2CountData, 
                                    colData = V_alg2ColData, 
                                    design = ~ condition)
head(ddsV_alg2) #review V_alg2 data set

#V. tub 
V_tubCountData <- gramNegCountData[ ,c(1,7)]
head(V_tubCountData)
V_tubColData <- gramNegColData[c(1,7), ]
head(V_tubColData)
ddsV_tub <- DESeqDataSetFromMatrix(countData = V_tubCountData, 
                                   colData = V_tubColData, 
                                   design = ~ condition)
head(ddsV_tub) #review V_aes data set


as.data.frame( colData(ddsV_aes) )
as.data.frame( colData(ddsV_ang) )
as.data.frame( colData(ddsV_alg1) )
as.data.frame( colData(ddsV_alg2) )
as.data.frame( colData(ddsV_tub) )



#Running the DEG pipeline for no replication (similar structure and code order as sr320 "SCRIPT_DESeq_LT_no replication.R"
#LPS

ddsLPS_DE <- DESeq(ddsLPS)
ddsV_aes_DE <- DESeq(ddsV_aes)
ddsV_ang_DE <- DESeq(ddsV_ang)
ddsV_alg1_DE <- DESeq(ddsV_alg1)
ddsV_alg2_DE <- DESeq(ddsV_alg2)
ddsV_tub_DE <- DESeq(ddsV_tub)

#Inspect results table
#to extract just log2fold change and p values
resLPS <- results(ddsLPS_DE)
head(resLPS)
resV_aes <- results(ddsV_aes_DE)
head(resV_aes)
resV_ang <- results(ddsV_ang_DE)
head(resV_ang)
resV_alg1 <- results(ddsV_alg1_DE)
head(resV_alg1)
resV_alg2 <- results(ddsV_alg2_DE)
head(resV_alg2)
resV_tub <- results(ddsV_tub_DE)
head(resV_tub)


####Experimenting with Alternative Method used in DESeq with no replicates (Adapted from Stephen Roberts github script "SCRIPT_DESeq_CLAM_no replication.R" )
ddsLPS <- DESeqDataSetFromMatrix(countData = LPSCountData, 
                                 colData = LPSColData, 
                                 design = ~ condition)
colData(ddsLPS)$condition<-factor(colData(ddsLPS)$condition,levels=c("control
                                                                     ","treatment"))  #set factors
ddsLPS <- estimateSizeFactors(ddsLPS)
sizeFactors(ddsLPS)

ddsLPS <- estimateDispersions( ddsLPS, fitType="local") 



#12 hr challenge
resoshv1_12_dfSig <- resoshv1_12_df[ which(resoshv1_12_df$padj < 0.1 ), ]
head( resoshv1_12_dfSig[ order( resoshv1_12_dfSig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resoshv1_12_dfSig[ order( resoshv1_12_dfSig$log2FoldChange ), ] ) #tail for strongest up regulation

#24hr challenge
resoshv1_24_dfSig <- resoshv1_24_df[ which(resoshv1_24_df$padj < 0.1 ), ]
head( resoshv1_24_dfSig[ order( resoshv1_24_dfSig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resoshv1_24_dfSig[ order( resoshv1_24_dfSig$log2FoldChange ), ] ) #tail for strongest up regulation

#48 hr challenge
resoshv1_48_dfSig <- resoshv1_48_df[ which(resoshv1_48_df$padj < 0.1 ), ]
head( resoshv1_48_dfSig[ order( resoshv1_48_dfSig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resoshv1_48_dfSig[ order( resoshv1_48_dfSig$log2FoldChange ), ] ) #tail for strongest up regulation

#120 hr challenge
resoshv1_120_dfSig <- resoshv1_120_df[ which(resoshv1_120_df$padj < 0.1 ), ]
head( resoshv1_120_dfSig[ order( resoshv1_120_dfSig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resoshv1_120_dfSig[ order( resoshv1_120_dfSig$log2FoldChange ), ] ) #tail for strongest up regulation


write.csv( as.data.frame(resoshv1_12_df), file="resoshv1_12_df.csv")
write.csv( as.data.frame(resoshv1_12_dfSig), file="resoshv1_12_dfSig.csv")
write.csv( as.data.frame(resoshv1_24_df), file="resoshv1_24_df.csv")
write.csv( as.data.frame(resoshv1_24_dfSig), file="resoshv1_24_dfSig.csv")
write.csv( as.data.frame(resoshv1_48_df), file="resoshv1_48_df.csv")
write.csv( as.data.frame(resoshv1_48_dfSig), file="resoshv1_48_dfSig.csv")
write.csv( as.data.frame(resoshv1_120_df), file="resoshv1_120_df.csv")
write.csv( as.data.frame(resoshv1_120_dfSig), file="resoshv1_120_dfSig.csv")
