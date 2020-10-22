#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --export=NONE
#SBATCH -o /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/SCRIPTS/Script_out_error_files/module_apop_10_22_2020_out
#SBATCH -e /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA/SCRIPTS/Script_out_error_files/module_apop_10_22_2020_error

echo "START $(date)"

cd /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/WGCNA

# Modules that were already exported and searched for are commented out

# Dermo_Tol: MEturquoise
#grep -f C_vir_rtracklayer_apop_product_final_rnaID.txt CytoscapeInput-edges-Dermo_Tol_fullturquoise.txt > Dermo_Tol_fullturquoise_apop_hits.txt

# Dermo_Sus: MElightpink4
grep -f C_vir_rtracklayer_apop_product_final_rnaID.txt CytoscapeInput-edges-Dermo_Sus_fulllightpink4.txt > Dermo_Sus_fulllightpink4_apop_hits.txt

# deLorg Res
#grep -f C_gig_rtracklayer_apop_product_final_transcript_ID.txt CytoscapeInput-edges-deLorg_Res_fullturquoise.txt > deLorg_Res_fullturquoise_apop_hits.txt

# de Lorg sus
#grep -f C_gig_rtracklayer_apop_product_final_transcript_ID.txt CytoscapeInput-edges-deLorg_Sus_fullturquoise.txt > deLorg_Sus_fullturquoise_apop_hits.txt

# He:  MEpurple, MEyellow
#grep -f C_gig_rtracklayer_apop_product_final_transcript_ID.txt CytoscapeInput-edges-He_fullpurple.txt > He_fullpurple_apop_hits.txt
#grep -f C_gig_rtracklayer_apop_product_final_transcript_ID.txt CytoscapeInput-edges-He_fullyellow.txt > He_fullyellow_apop_hits.txt

# Zhang: MEblack
#grep -f C_gig_rtracklayer_apop_product_final_transcript_ID.txt CytoscapeInput-edges-Zhang_fullblack.txt > Zhang_fullblack_apop_hits.txt

# Rubio: MEmagenta, MEturquoise, MEblue, MEbrown, MEblack
#grep -f C_gig_rtracklayer_apop_product_final_transcript_ID.txt CytoscapeInput-edges-Rubio_fullblue.txt > Rubio_fullblue_apop_hits.txt
#grep -f C_gig_rtracklayer_apop_product_final_transcript_ID.txt CytoscapeInput-edges-Rubio_fullbrown.txt > Rubio_fullbrown_apop_hits.txt
#grep -f C_gig_rtracklayer_apop_product_final_transcript_ID.txt CytoscapeInput-edges-Rubio_fullmagenta.txt > Rubio_fullmagenta_apop_hits.txt
#grep -f C_gig_rtracklayer_apop_product_final_transcript_ID.txt CytoscapeInput-edges-Rubio_fullturquoise.txt > Rubio_fullturquoise_apop_hits.txt
grep -f C_gig_rtracklayer_apop_product_final_transcript_ID.txt CytoscapeInput-edges-Rubio_fullblack.txt > Rubio_fullblack_apop_hits.txt

## Pro_RE22_Pro_RI:MEturquoise    , MEdarkslateblue, MEsteelblue , MEroyalblue
#grep -f C_vir_rtracklayer_apop_product_final_rnaID.txt CytoscapeInput-edges-Pro_RE22_Pro_fulldarkslateblue.txt > Pro_RE22_Pro_fulldarkslateblue_apop_hits.txt
#grep -f C_vir_rtracklayer_apop_product_final_rnaID.txt CytoscapeInput-edges-Pro_RE22_Pro_fullroyalblue.txt > Pro_RE22_Pro_fullroyalblue_apop_hits.txt
#grep -f C_vir_rtracklayer_apop_product_final_rnaID.txt CytoscapeInput-edges-Pro_RE22_Pro_fullsteelblue.txt > Pro_RE22_Pro_fullsteelblue_apop_hits.txt
#grep -f C_vir_rtracklayer_apop_product_final_rnaID.txt CytoscapeInput-edges-Pro_RE22_Pro_fullturquoise.txt > Pro_RE22_Pro_fullturquoise_apop_hits.txt
grep -f C_vir_rtracklayer_apop_product_final_rnaID.txt CytoscapeInput-edges-Pro_RE22_Pro_fullwhite.txt > Pro_RE22_Pro_fullwhite_apop_hits.txt
grep -f C_vir_rtracklayer_apop_product_final_rnaID.txt CytoscapeInput-edges-Pro_RE22_Pro_fullskyblue3.txt > Pro_RE22_Pro_fullskyblue3_apop_hits.txt

echo "DONE $(date)"
