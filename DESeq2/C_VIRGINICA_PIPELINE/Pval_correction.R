#Script to perform p-value correction 

####p-value correction for both sets of results ####
#First Visualize histograms
#histogram of P- values to visualize any "hills" or "U shape"
# hill means variance of the null distribution too high, U shape means variance assumed too low
hist(resRODTran_05$pvalue, breaks= 20, col = "grey")
hist(resRODTran_05_Sig$pvalue, breaks = 20, col = "grey") #hill

#remove filtered out genes by independent filtering, they have NA adj. pvals
resRODTran_05_df <- resRODTran_05[ !is.na(resRODTran_05$padj), ]

#remove genes with NA pvals (outliers)
resRODTran_05_df <- resRODTran_05_df[ !is.na(resRODTran_05_df$pvalue), ]

#remove adjsuted pvalues, since we add the fdrtool results later on (based on the correct p-values)
resRODTran_05_df <- resRODTran_05_df[, -which(names(resRODTran_05_df) == "padj")]

#use z-scores as input to FDRtool to re-estimate the p-value
FDR.resRODTran_05_df <- fdrtool(resRODTran_05_df$stat, statistic= "normal", plot = T)

#add values to the results data frame, also ad new BH- adjusted p-values
resRODTran_05_df[,"padj"] <- p.adjust(FDR.resRODTran_05_df$pval, method = "BH")

#replot corrected p-values 
hist(FDR.resRODTran_05_df$pval, col = "royalblue4",
     main = "Correct null model ROD Transcript Count", xlab = "CORRECTED p-values")

#Check how many genes have BH adjusted p values of less than 0.05 after P-value correction?
sum( resRODTran_05_df$padj < 0.05, na.rm=TRUE ) #3138

#Subset the results table to the differentially expressed genes under FDR 0.1, order the Log2FC table first by strongest down regulation
resRODTran_05_dfSig <- resRODTran_05_df[ which(resRODTran_05_df$padj < 0.05 ), ]
head( resRODTran_05_dfSig[ order( resRODTran_05_dfSig$log2FoldChange ), ] ) #head for strongest downregulation
tail( resRODTran_05_dfSig[ order( resRODTran_05_dfSig$log2FoldChange ), ] ) #tail for strongest up regulation
summary(resRODTran_05_dfSig)
resRODTran_05_df_non_Sig <- resRODTran_05_df[ which(resRODTran_05_df$padj > 0.05 ), ]
summary(resRODTran_05_df_non_Sig)