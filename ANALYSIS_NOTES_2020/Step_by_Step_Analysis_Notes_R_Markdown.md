# Step-by-Step Analysis Notes R Markdown
## by Erin Roberts, PhD Candidate University of Rhode Island,
### Jan 2020 -

This R markdown project notebook details my full analysis steps for performing comparative transcriptomics to
compare apoptosis gene expression between disease challenges.

## 1/29/2020-1/30/2020 Data Preparation, Data Storage Plan

  ### 1. Review Which SRAs to download and analyze from each full experiment:

      - C. virginica Probiotic = Download and analyze all available SRAs
      - C. virginica ROD = Download and analyze all available SRAs
      - C. virginica Dermo challenge = use all avialable transcriptomes (provided by D. Proestou)
      but merge technical replicates from samples that had two technical replicates (see below)
      - C. gigas He et al. 2015 OsHV1 = Download and analyze all available SRAs
      - C. gigas Zhang et al., 2015 Vibrio = Download and analyze all available SRAs (only the
      transcriptomes that were challenged with individual strains of Vibrio, LPS, M. lut and PBS on NCBI)
      - C. gigas Rubio et al. 2019 Vibrio = Download and analyze all available SRAs
      - C. gigas de Lorgeril et al., 2017 = Download and analyze only the samples from families 11 and
      21 that were infected with OsHV-1 during a "natural" infection and sequenced. These were all sequenced
      paired end on hiseq 500. 42 total transcriptomes

  ### 2. Create New folders on bluewaves cluster where data will be housed.

      - C. virginica Raw data will be in the following folder, with a separate directory per experiment:
      ```
      /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/2020_Raw_Transcriptome_Data
      ```
      -C. gigas Raw data will be in the following folder, with a separate directory per experiment:
      ```
      /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/Bac_Viral_subset/2020_Raw_Transcriptome_Data
      ```

  ### 3. Make an individual text file for each experiment with the SRA's for each (not for Dermo transcriptomes),
  and Create Full_SRA_PE_list.txt and Full_SRA_SE_list.txt that has full SRA list (with species combined)

      ```
        $ pwd
        /Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/Transcriptome_Bluewaves_Cluster_Analysis_2020/Bio_projects_Sample_Metadata        $ ls *SRA_ID.txt

        $ ls *.txt
            C_gig_He_2015_OsHV1_SRA_ID.txt		C_gig_deLorgeril_OsHV1_SRA_ID.txt	Full_SRA_PE_list.txt
            C_gig_Rubio_Vibrio_SRA_ID.txt		C_vir_Probiotic_SRA_ID.txt		Full_SRA_SE_list.txt
            C_gig_Zhang_Vibrio_SRA_ID.txt		C_vir_ROD_SRA_ID.txt
      ```

  ### 6. Add metadata for all the samples to "Organized_SRA_info.xlsx" spreadsheet for reference later.

        - Completed

## 1/31/2020 Data Download from SRA database (except for Dermo transcriptomes from Dina)

  ### 1. Write script to download each SRA in a separate loop and place files in correct folder.

    - Script added to github.
      Path: /Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/Transcriptome_Bluewaves_Cluster_Analysis_2020/2020_SCRIPTS/fetch_all_SRA_2020.sh
    - Helpful link regarding use of SRA-toolkit: https://reneshbedre.github.io/blog/fqutil.html
    - Version of SRA toolkit used: SRA-Toolkit/2.9.0-centos_linux64
    - Out and error files in `/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_files/`
          a. Outfile named `fetch_SRA_output_1_31_2020`
          b. Error file named `fetch_SRA_error_1_31_2020`  
    - After download, check that all SRA's were successfully loaded onto the cluster.
    - The prefetch command downloads a `*.sra` version of each file in `$HOME/ncbi/public/sra`. You can change the default location of where this data is dowloaded by using `vdb-config`, however changing the config file for a shared cluster resource is not recommended. I'm going to download data for each experiment separately, run a check sum using `vdb-validate`, and then delete the SRA files from home, and then start the next group of transcriptomes in my script and go down the list.
      -Error messages received during download regarding timeout. Disregard message as long as you receive read and written output.
        `2020-02-03T15:54:05 prefetch.2.9.0 sys: timeout exhausted while reading file within network system module - mbedtls_ssl_read returned -76 ( NET - Reading information from the socket failed )
        2020-02-03T15:54:05 prefetch.2.9.0 int: timeout exhausted while reading file within network system module - ?-]: Cannot KStreamRe
        `
      -As long as every file receives an output report like below, the download was successful
        `Read 11992053 spots for SRR2002962
        Written 11992053 spots for SRR2002962`
    -Check that the number of samples read and written in the output file and the number of total files match
        `$ grep 'Read' fetch_SRA_output_He_OsHV1_2_3_2020 | wc -l`
    -Added `vdb-validate --option-file *SRA_ID.txt` for each file in the script. This file is saved in each folder.

    - Experiments finished downloading and reviewing checksum:
      1. HE OsHV1 - downloaded and checksum confirmed
      2. Zhang Vibrio - downloaded and checksum confirmed
      3. ROD - downloaded and checksum confirmed
      4. de lorgeril OsHV1- downloaded and confirmed
      5. Rubio Vibrio - downloading
      6. Probiotic - downloading

    - Download data as several separate scripts. There are individual output and error files for each download for future reference:


  ### 2. Merge technical replicates from Proestou et al. 2015 transcriptomes. (decided to do in DESeq2)

    - Options:
      a. Merge raw SRA fastq files for samples with technical replicates
      b. Combined in DESeq2 using collapseReplicates() function. See site for details: https://support.bioconductor.org/p/85536/
        1. see the collapseReplicates section of this DESeq tutorial: https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf for code. Thus function merges the count tables before proceeeding to calculate DeSeq2.

## 2/1/2020-2/4/2020 Adapter Trimming, Quality Filtering of data

  ### 1. Combined trimming scripts for both C_gig and C_vir into single script. In order to maintain folder organization of files from each experiment, I'm going to run a loop through every filtering and trimming command on all the files for each experiment separately.
      a. Though I have them all in a single script. I'm commenting out parts of the script and running them individually to check. A separate output file will exist for each. Deleting intermediate files between each step.
      b. Starting with the Dermo transcriptomes that I know are fully downloaded, while the other SRAs are finishing being downloaded.  
      c. Running Dermo transcriptome script separately (because file name format is different). Script called "01_SRA_Trim_Filter_Dermo_only.sh"
          - finished successfully
          - Checked that all files had been pre-processed. 194 total raw files and processed files.
          - moved raw data into one lower folder `Dermo_Raw_Transcriptomes`
          - Compressed the raw data folder in interactive mode
              `$tar -zcvf Dermo_Raw_Transcriptomes.archive.tar.gz Dermo_Raw_Transcriptomes`
      d. Running HE OsHV1 preprocessing
          - Completed
          - Checked all files had been processed
             `$ ls *.filter.gz | wc -l ` 32
             `$ ls ./He_OsHV1_Raw_Transcriptomes/* | wc -l` 32
          - Moved raw data into `He_OsHV1_Raw_Transcriptomes`
          - Compressed raw data folder
      e. Starting Zhang Vibrio Pre-processing
          - Completed
          - Checked all files had been processed
          `[erin_roberts@bluewaves Zhang_Vibrio_Raw_Transcriptomes]$ ls * | wc -l ` 9
          `[erin_roberts@bluewaves C_gig_Zhang_Vibrio_SRA]$ ls *.gz | wc -l ` 9
          - Moved raw data into `Zhang_Vibrio_Raw_Transcriptomes`
          - Compressed Raw data
            `$ tar -zcvf Zhang_Vibrio_Raw_Transcriptomes.archive.tar.gz Zhang_Vibrio_Raw_Transcriptomes`
      f. Starting ROD Pre-processing
      g. Starting deLorgeril OsHV1 pre-processing
      h. Starting Rubio Vibrio pre-processing
      i. Starting Probiotic pre-processing

      - All scripts have unique output and error files. Ended up running pre-processing on data subsets one at a time in order to expedite process.

## 2/4/2020 - 2/5/2020 Building HISAT2 Script and Running HISAT2

  1. Building new HISAT index.
      a. Note, a new version of HISAT is now available on the cluster. I'm going to use this version (HISAT2/2.1.0-foss-2016b)
          - Notes also that the code I'm adding to github is all combined into one script, though I ran the same code in the cluster separated into multiple scripts to aid in running multiple scripts at once. The code used however was not changed.
      a. Checking how HISAT2 genome indexes we made. Creating new indexes for both C_vir and C_gig using updated software version

          - Checking C_virginica genome file on hand which was edited from `ref_C_virginica-3.0_top_level.gff3` to remove a space from header that was making it incompatible with Stringtie annotation. Checking also that mitochondrial DNA was added.

            `$ grep '^>' cvir_edited.fa
              >NC_035780.1 Crassostrea virginica isolate RU13XGHG1-28 chromosome 1, C_virginica-3.0, whole genome shotgun sequence
              >NC_035781.1 Crassostrea virginica isolate RU13XGHG1-28 chromosome 2, C_virginica-3.0, whole genome shotgun sequence
              >NC_035782.1 Crassostrea virginica isolate RU13XGHG1-28 chromosome 3, C_virginica-3.0, whole genome shotgun sequence
              >NC_035783.1 Crassostrea virginica isolate RU13XGHG1-28 chromosome 4, C_virginica-3.0, whole genome shotgun sequence
              >NC_035784.1 Crassostrea virginica isolate RU13XGHG1-28 chromosome 5, C_virginica-3.0, whole genome shotgun sequence
              >NC_035785.1 Crassostrea virginica isolate RU13XGHG1-28 chromosome 6, C_virginica-3.0, whole genome shotgun sequence
              >NC_035786.1 Crassostrea virginica isolate RU13XGHG1-28 chromosome 7, C_virginica-3.0, whole genome shotgun sequence
              >NC_035787.1 Crassostrea virginica isolate RU13XGHG1-28 chromosome 8, C_virginica-3.0, whole genome shotgun sequence
              >NC_035788.1 Crassostrea virginica isolate RU13XGHG1-28 chromosome 9, C_virginica-3.0, whole genome shotgun sequence
              >NC_035789.1 Crassostrea virginica isolate RU13XGHG1-28 chromosome 10, C_virginica-3.0, whole genome shotgun sequence
              >NC_007175.2 Crassostrea virginica mitochondrion, complete genome`
            - Checking C_gigas genome file. Only on the cluster was `Crassostrea_gigas.gff`. A newer version of C_gigas genome is avaiable on NCBI now. Downloading that version and will use it to create HISAT2 index for C_gigas sequences.
                - The C. gigas genome file contains all the genomic scaffolds (it has never been assembled to chromosome level like oysters) and the mitochondrial genome.

            - Creating new indexes using code in `02_Build_Hisat_Indexes.sh`. Indices created!
  2. Creating script to perform HISAT 2 alignment and SAMtools sorting on files for each experiment.
