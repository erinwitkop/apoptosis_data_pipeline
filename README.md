# Pipeline to perform differential transcript analysis with transcriptomes
## By: Erin M. Roberts
## 1/29/2020

This pipeline includes scripts I generated to analyze transcriptomes for differential expression with samples from *C. gigas* and *C. virginica*.

For each species they perfom the following analysis
1. Adapter trimming, quality trimming using BBTools
2. Mapping with HISAT2 to the reference genome
3. Alignment with Stringtie to the reference annotation
4. Differential Transcript analysis of count data with DESeq2

The following script is only available for *C. gigas* transcriptomes.
1. Gene set enrichment analysis (GSEA) with topGO

The /SCRIPTS folder contains bash scripts to be executed from a cluster computing environment. The /Bac_Viral_Subset
folder contains scripts to process *C. gigas* transcriptomes. The /C_Virginica_Subset folder contains scripts to process
*C. virginica* transcriptomes. 

The DESeq2 folder contains R scripts to perform DESeq2 differential transcript analysis and GSEA (Gene Set Enrichment Analysis).
Scripts again are separated into separate folders by species. 

Adidtionally, I have provided a folder called "Streamlined Pipeline Tutorial" where I provide a Markdown file with resources for where to find each tool, step by step instructions of how to use the pipeline, and rationale for how I set up the pipeline. I also include relevant papers as well as a powerpoint providing an introduction to reference-based RNA-seq in general.

To cite this work: Roberts, E.M. 2020. "Pipeline to perform differential expression analysis with transcriptomes".
https://github.com/erinroberts/apoptosis_data_pipeline. The University of Rhode Island. 
 


