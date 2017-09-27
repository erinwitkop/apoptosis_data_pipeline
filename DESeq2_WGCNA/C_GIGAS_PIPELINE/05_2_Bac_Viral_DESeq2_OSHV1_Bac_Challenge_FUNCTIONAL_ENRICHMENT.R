#05_2_Bac_Viral_DESeq2_OSHV1_Bac_Challenge_FUNCTIONAL_ENRICHMENT.R

####BLAST2GO the MSTRG SIG genes (AND MSTRG NON SIG GENES), then perform gene set enrichment on these and pull out more specific apoptosis genes ####
##subset out only the genes that have "MSTRG ID"
#lookup these lines in the stringtie.merge file to find the sequence for them (SEE 06_MSTRG_isolate_getfasta.sh), then put them into BLAST2GO
#Sig MSTRG's from OsHV1
resoshv1Tran_05_dfSig_Transcript_MSTRG <- resoshv1Tran_05_dfSig[grep("MSTRG", rownames(resoshv1Tran_05_dfSig)), ] #for Significant genes
OsHv1_MSTRGID_tran_Sig <- as.data.frame(rownames(resoshv1Tran_05_dfSig_Transcript_MSTRG))
head(OsHv1_MSTRGID_tran_Sig)
write.table(OsHv1_MSTRGID_tran_Sig[,1], file="OsHv1_MSTRGID_tran_Sig", sep= "\t")

#NON SI MSTRGs OsHV1
resoshv1Tran_05_df_non_Sig_Transcript_MSTRG <- resoshv1Tran_05_df_non_Sig[grep("MSTRG", rownames(resoshv1Tran_05_df_non_Sig)),]
OsHV1_MSTRGID_tran_non_Sig <- as.data.frame(rownames(resoshv1Tran_05_df_non_Sig_Transcript_MSTRG))
head(OsHV1_MSTRGID_tran_non_Sig)
write.table(OsHV1_MSTRGID_tran_non_Sig[,1], file="OsHv1_MSTRGID_tran_non_Sig", sep= "\t")

#Sig MSTRGs from Bacterial Challenge
resBacTran_05_dfSig_Transcript_MSTRG <- resBacTran_05_dfSig[grep("MSTRG", rownames(resBacTran_05_dfSig)), ]
Bac_MSTRGID_tran_Sig <- as.data.frame(rownames(resoshv1Tran_05_dfSig_Transcript_MSTRG))
head(Bac_MSTRGID_tran_Sig)
write.table(Bac_MSTRGID_tran_Sig[,1], file="Bac_MSTRGID_tran_Sig", sep= "\t")

#NON SIG MSTRGs from Bacterial Challenge
resBacTran_05_df_non_Sig_MSTRG <- resBacTran_05_df_non_Sig[grep("MSTRG", rownames(resBacTran_05_df_non_Sig)),]
Bac_MSTRG_tran_non_Sig <- as.data.frame(rownames(resBacTran_05_df_non_Sig_MSTRG))
head(Bac_MSTRG_tran_non_Sig)
write.table(write.table(Bac_MSTRG_tran_non_Sig[,1], file="Bac_MSTRGID_tran_non_Sig", sep= "\t"))

#USE 06_MSTRG_isolate_getfasta.sh script to get the sequence of genes while using the bam files of each and the sequence
#INSERT SIG SEQUENCES SEQUENCES INTO BLAST2GO

#TOPGO and REVIGO.R to get functional enrichment

# Cytoscape word cloud ?
