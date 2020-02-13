# Script to Analyze Differential Expression of Transcriptomes using DESeq2
# Erin Roberts, PhD Candidate University of Rhode Island 
# 2/13/2020

# This script will calculate differential expression of apoptosis genes from C. gigas and C. virginica, for each experiment separately.
# Comparisons and formulas used to calculate differential expression for each experiment will be unique to that experiment, specifically
# tailored to when the infection was most acute. Challenge group samlpes will always be compared to their own control.

#### Load Packages ####

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
library(reshape2)
library(plyr)
library(Repitools)
library(purrr)
library(tibble)

# VERSIONS (see sessionInfo at bottom of script for full information)
# R version 3.6.1 (2019-07-05)
# DESeq2_1.24.0 

######





#### SESSION INFO FOR RUNNING SCRIPTS ####

sessionInfo()

# R version 3.6.1 (2019-07-05)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Mojave 10.14

# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

#locale:
  # [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
  # [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] tibble_2.1.3                purrr_0.3.3                 Repitools_1.30.0            plyr_1.8.5                 
# [5] reshape2_1.4.3              UpSetR_1.4.0                rtracklayer_1.44.3          stringr_1.4.0              
# [9] tidyr_1.0.0                 fission_1.4.0               genefilter_1.66.0           apeglm_1.6.0               
# [13] questionr_0.7.0             RColorBrewer_1.1-2          pheatmap_1.0.12             dplyr_0.8.3                
# [17] magrittr_1.5                ggplot2_3.2.1               DESeq2_1.24.0               SummarizedExperiment_1.14.1
# [21] DelayedArray_0.10.0         BiocParallel_1.18.1         matrixStats_0.54.0          Biobase_2.44.0             
# [25] GenomicRanges_1.36.0        GenomeInfoDb_1.20.0         IRanges_2.18.2              S4Vectors_0.22.0           
# [29] BiocGenerics_0.30.0        

# loaded via a namespace (and not attached):
# [1] backports_1.1.5          Hmisc_4.2-0              R.devices_2.16.0         aroma.light_3.14.0       R.rsp_0.43.1            
# [6] lazyeval_0.2.2           splines_3.6.1            listenv_0.7.0            digest_0.6.23            htmltools_0.3.6         
# [11] gdata_2.18.0             Rsolnp_1.16              checkmate_1.9.4          memoise_1.1.0            BSgenome_1.52.0         
# [16] aroma.apd_0.6.0          cluster_2.1.0            limma_3.40.6             globals_0.12.4           Biostrings_2.52.0       
# [21] annotate_1.62.0          R.utils_2.9.0            aroma.core_3.2.0         colorspace_1.4-1         blob_1.2.0              
# [26] xfun_0.9                 crayon_1.3.4             RCurl_1.95-4.12          zeallot_0.1.0            survival_2.44-1.1       
# [31] glue_1.3.1               R.huge_0.9.0             gtable_0.3.0             zlibbioc_1.30.0          XVector_0.24.0          
# [36] R.cache_0.13.0           scales_1.1.0             vsn_3.52.0               DBI_1.0.0                edgeR_3.26.8            
# [41] miniUI_0.1.1.1           Rcpp_1.0.3               xtable_1.8-4             emdbook_1.3.11           htmlTable_1.13.1        
# [46] foreign_0.8-72           bit_1.1-14               preprocessCore_1.46.0    Formula_1.2-3            truncnorm_1.0-8         
# [51] htmlwidgets_1.3          gplots_3.0.1.1           acepack_1.4.1            pkgconfig_2.0.3          XML_3.98-1.20           
# [56] R.methodsS3_1.7.1        nnet_7.3-12              locfit_1.5-9.1           DNAcopy_1.58.0           tidyselect_0.2.5        
# [61] rlang_0.4.2              later_0.8.0              AnnotationDbi_1.46.1     munsell_0.5.0            tools_3.6.1             
# [66] RSQLite_2.1.2            yaml_2.2.0               knitr_1.24               bit64_0.9-7              caTools_1.17.1.2        
# [71] future_1.14.0            mime_0.7                 R.oo_1.22.0              compiler_3.6.1           rstudioapi_0.10         
# [76] affyio_1.54.0            geneplotter_1.62.0       stringi_1.4.5            highr_0.8                gsmoothr_0.1.7          
# [81] lattice_0.20-38          Matrix_1.2-17            vctrs_0.2.1              pillar_1.4.3             lifecycle_0.1.0         
# [86] BiocManager_1.30.4       data.table_1.12.8        bitops_1.0-6             httpuv_1.5.1             R.filesets_2.13.0       
# [91] R6_2.4.1                 latticeExtra_0.6-28      affy_1.62.0              promises_1.0.1           KernSmooth_2.23-15      
# [96] gridExtra_2.3            aroma.affymetrix_3.2.0   codetools_0.2-16         Ringo_1.48.0             MASS_7.3-51.4           
# [101] gtools_3.8.1             assertthat_0.2.1         withr_2.1.2              GenomicAlignments_1.20.1 Rsamtools_2.0.0         
# [106] GenomeInfoDbData_1.2.1   grid_3.6.1               rpart_4.1-15             coda_0.19-3              PSCBS_0.65.0            
# [111] bbmle_1.0.20             numDeriv_2016.8-1.1      shiny_1.3.2              base64enc_0.1-3         

