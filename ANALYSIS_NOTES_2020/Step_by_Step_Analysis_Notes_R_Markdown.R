# Step-by-Step Analysis Notes R Markdown
## by Erin Roberts, PhD Candidate University of Rhode Island,
### Jan 2020 -

This R markdown project notebook details my full analysis steps for performing comparative transcriptomics to
compare apoptosis gene expression between disease challenges.

# 1/29/2020-1/30/2020 Data Preparation, Data Storage Plan

  ## 1. Review Which SRAs to download and analyze from each full experiment

      - C. virginica Probiotic = Download and analyze all available SRAs
      - C. virginica ROD = Download and analyze all available SRAs
      - C. virginica Dermo challenge = use all avialable transcriptomes (provided by D. Proestou)
      but merge technical replicates from samples that had two technical replicates (see below)
      - C. gigas He et al. 2015 OsHV1 = Download and analyze all available SRAs
      - C. gigas Zhang et al., 2015 Vibrio = Download and analyze all available SRAs (only the
      transcriptomes that were challenged with individual strains of Vibrio, LPS, M. lut and PBS on NCBI)
      - C. gigas Rubio et al. 2019 Vibrio = Download and analyze all available SRAs
      - C. gigas de Lorgeril et al., 2017 = Download and analyze only the samples that
      were injected with OsHV-1?

  ## 2. Create New folders on bluewaves cluster where data will be housed.


  ## 3.  Merge technical replicates from Proestou et al. 2015 transcriptomes

  ## 4. C


  ## 5. Make an individual File for each experiment with the SRA's for each (not for Dermo transcriptomes)

        $ cd /Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/Transcriptome_Bluewaves_SRA_Data_Analysis_2020/Bio_projects_Sample_Metadata
        $ ls *SRA_ID.txt
        $ C_gig_He_2015_OsHV1_SRA_ID.txt
          C_vir_Probiotic_SRA_ID.txt
          C_gig_Rubio_Vibrio_SRA_ID.txt
          C_vir_ROD_SRA_ID.txt
          C_gig_Zhang_Vibrio_SRA_ID.txt

  ## Add metadata for all the samples to "Organized_SRA_info.xlsx" spreadsheet for reference later.


  ## Create Full_SRA_PE_list.txt and Full_SRA_SE_list.txt that has full SRA list with each species

# 1/30/2020 Data Download from SRA database (except for Dermo transcriptomes from Dina)
