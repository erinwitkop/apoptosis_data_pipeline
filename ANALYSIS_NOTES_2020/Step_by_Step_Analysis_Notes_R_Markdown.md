# Step-by-Step Analysis Notes R Markdown
## by Erin Roberts, PhD Candidate University of Rhode Island,
### Jan 2020 -

This R markdown project notebook details my full analysis steps for performing comparative transcriptomics to
compare apoptosis gene expression between disease challenges.

## 1/29/2020-1/30/2020 Data Preparation, Data Storage Plan

### 1. Review Which SRAs to download and analyze from each full experiment:

* C. virginica Probiotic = Download and analyze all available SRAs
* C. virginica ROD = Download and analyze all available SRAs
* C. virginica Dermo challenge = use all avialable transcriptomes (provided by D. Proestou) but merge technical replicates from samples that had two technical replicates (see below)
* C. gigas He et al. 2015 OsHV1 = Download and analyze all available SRAs
* C. gigas Zhang et al., 2015 Vibrio = Download and analyze all available SRAs (only the transcriptomes that were challenged with individual strains of Vibrio, LPS, M. lut and PBS on NCBI)
* C. gigas Rubio et al. 2019 Vibrio = Download and analyze all available SRAs
* C. gigas de Lorgeril et al., 2017 = Download and analyze only the samples from families 11 and 21 that were infected with OsHV-1 during a "natural" infection and sequenced. These were all sequenced paired end on hiseq 500. 42 total transcriptomes

* Counting total number of Transcriptome samples (to report in methods) 181 total SRA samples

         `# deLorgeril OsHV1 = 42 samples
         erin_roberts@bluewaves 2020_Raw_Transcriptome_Data]$ ls ./C_gig_deLorgeril_OsHV1_SRA/*.filter.gz | wc -l # 84
         erin_roberts@bluewaves C_gig_deLorgeril_OsHV1_SRA]$ cat C_gig_deLorgeril_OsHV1_SRA_ID.txt | wc -l # 42

         # Rubio Vibrio 18
         $ cat C_gig_Rubio_Vibrio_SRA_ID.txt | wc -l # 17 (actually 18)
         erin_roberts@bluewaves C_gig_Rubio_Vibrio_SRA]$ ls *_1.fastq.gz | wc -l #18

         # Zhang Vibrio 9
         $ cat C_gig_Zhang_Vibrio_SRA_ID.txt | wc -l # 9
         $ ls *.filter.gz | wc -l # 9

         # He OsHV1 32
         $ ls *.filter.gz | wc -l # 32
         $ cat C_gig_He_2015_OsHV1_SRA_ID.txt | wc -l # 32

         # ROD 12
         $ cat C_vir_ROD_SRA_ID.txt | wc -l # 12
         $ ls *.filter.gz | wc -l # 12

         # Probiotic 6
         $  ls *_1.fastq.gz | wc -l # 6
         $ cat C_vir_Probiotic_SRA_ID.txt | wc -l # 6

         # Dermo 62
          # count ones with Technical replicates first
          $ ls *TechRep1.R1*filter.gz | wc -l # 35
          # then count other samples
          $ ls DCS2015*R1*filter.gz | wc -l # 27


### 2. Create New folders on bluewaves cluster where data will be housed.

* C. virginica Raw data will be in the following folder, with a separate directory per experiment:

      `/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/2020_Raw_Transcriptome_Data`

* C. gigas Raw data will be in the following folder, with a separate directory per experiment:

      `/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/Bac_Viral_subset/2020_Raw_Transcriptome_Data`

### 3. Make an individual text file for each experiment with the SRA's for each (not for Dermo transcriptomes), and Create Full_SRA_PE_list.txt and Full_SRA_SE_list.txt that has full SRA list (with species combined)

      ```
        $ pwd
        /Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/Transcriptome_Bluewaves_Cluster_Analysis_2020/Bio_projects_Sample_Metadata        $ ls *SRA_ID.txt

        $ ls *.txt
            C_gig_He_2015_OsHV1_SRA_ID.txt		C_gig_deLorgeril_OsHV1_SRA_ID.txt	Full_SRA_PE_list.txt
            C_gig_Rubio_Vibrio_SRA_ID.txt		C_vir_Probiotic_SRA_ID.txt		Full_SRA_SE_list.txt
            C_gig_Zhang_Vibrio_SRA_ID.txt		C_vir_ROD_SRA_ID.txt
      ```

### 6. Add metadata for all the samples to "Organized_SRA_info.xlsx" spreadsheet for reference later.

* Completed

## 1/31/2020 Data Download from SRA database (except for Dermo transcriptomes from Dina)

### 1. Write script to download each SRA in a separate loop and place files in correct folder.

* Script added to github.
      Path: /Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/Transcriptome_Bluewaves_Cluster_Analysis_2020/2020_SCRIPTS/fetch_all_SRA_2020.sh
* Helpful link regarding use of SRA-toolkit: https://reneshbedre.github.io/blog/fqutil.html
* Version of SRA toolkit used: SRA-Toolkit/2.9.0-centos_linux64
* Out and error files in    

      `/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts/Script_out_error_file/`
          a. Outfile named `fetch_SRA_output_1_31_2020`
          b. Error file named `fetch_SRA_error_1_31_2020`  

* After download, check that all SRA's were successfully loaded onto the cluster.
* The prefetch command downloads a `*.sra` version of each file in `$HOME/ncbi/public/sra`. You can change the default location of where this data is dowloaded by using `vdb-config`, however changing the config file for a shared cluster resource is not recommended. I'm going to download data for each experiment separately, run a check sum using `vdb-validate`, and then delete the SRA files from home, and then start the next group of transcriptomes in my script and go down the list.
       * Error messages received during download regarding timeout. Disregard message as long as you receive read and written output.
        `2020-02-03T15:54:05 prefetch.2.9.0 sys: timeout exhausted while reading file within network system module - mbedtls_ssl_read returned -76 ( NET - Reading information from the socket failed )
        2020-02-03T15:54:05 prefetch.2.9.0 int: timeout exhausted while reading file within network system module - ?-]: Cannot KStreamRe
        `
    * As long as every file receives an output report like below, the download was successful
        `Read 11992053 spots for SRR2002962
        Written 11992053 spots for SRR2002962`
    * Check that the number of samples read and written in the output file and the number of total files match
        `$ grep 'Read' fetch_SRA_output_He_OsHV1_2_3_2020 | wc -l`
    * Added `vdb-validate --option-file *SRA_ID.txt` for each file in the script. This file is saved in each folder.

    * Experiments finished downloading and reviewing checksum:
      1. HE OsHV1 - downloaded and checksum confirmed
      2. Zhang Vibrio - downloaded and checksum confirmed
      3. ROD - downloaded and checksum confirmed
      4. de lorgeril OsHV1- downloaded and confirmed
      5. Rubio Vibrio - downloaded and confirmed. Correct number of files present
      6. Probiotic - downloaded and confirmed. correct number of files present

* Download data as several separate scripts. There are individual output and error files for each download for future reference:


### 2. Merge technical replicates from Proestou et al. 2015 transcriptomes. (decided to do in DESeq2)

* Options:
      * Merge raw SRA fastq files for samples with technical replicates
      * Combined in DESeq2 using collapseReplicates() function. See site for details: https://support.bioconductor.org/p/85536/
          * see the collapseReplicates section of this DESeq tutorial: https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf for code. Thus function merges the count tables before proceeeding to calculate DeSeq2.

## 2/1/2020-2/4/2020 Adapter Trimming, Quality Filtering of data

### 1. Combined trimming scripts for both C_gig and C_vir into single script. In order to maintain folder organization of files from each experiment, I'm going to run a loop through every filtering and trimming command on all the files for each experiment separately.

* Though I have them all in a single script I'm commenting out parts of the script and running them individually to check. A separate output file will exist for each. Deleting intermediate files between each step.
* Starting with the Dermo transcriptomes that I know are fully downloaded, while the other SRAs are finishing being downloaded.  
* Running Dermo transcriptome script separately (because file name format is different). Script called "01_SRA_Trim_Filter_Dermo_only.sh"
  * finished successfully
  * Checked that all files had been pre-processed. 194 total raw files and processed files.
  * moved raw data into one lower folder `Dermo_Raw_Transcriptomes`
  * Compressed the raw data folder in interactive mode
          `$tar -zcvf Dermo_Raw_Transcriptomes.archive.tar.gz Dermo_Raw_Transcriptomes`
* Running HE OsHV1 preprocessing
  * Completed
  * Checked all files had been processed
          `$ ls *.filter.gz | wc -l ` 32
          `$ ls ./He_OsHV1_Raw_Transcriptomes/* | wc -l` 32
  * Moved raw data into `He_OsHV1_Raw_Transcriptomes`
  * Compressed raw data folder
* Starting Zhang Vibrio Pre-processing
  * Completed
  * Checked all files had been processed
          `[erin_roberts@bluewaves Zhang_Vibrio_Raw_Transcriptomes]$ ls * | wc -l ` 9
          `[erin_roberts@bluewaves C_gig_Zhang_Vibrio_SRA]$ ls *.gz | wc -l ` 9
  * Moved raw data into `Zhang_Vibrio_Raw_Transcriptomes`
  * Compressed Raw data
            `$ tar -zcvf Zhang_Vibrio_Raw_Transcriptomes.archive.tar.gz Zhang_Vibrio_Raw_Transcriptomes`
* Starting ROD Pre-processing
  * completed
  * Checking number of files
        `ls *.gz | wc -l 12` , 12 in both places
  * compressed raw data and put in `ROD_Raw_Transcriptomes`
        `tar -zcvf ROD_Raw_Transcriptomes.archive.tar.gz ROD_Raw_Transcriptomes/`
* Starting deLorgeril OsHV1 pre-processing
  * completed
  * Checking number of files
          `$ ls *_1*filter.gz | wc -l # 42`
  * compressed raw data

* Starting Rubio Vibrio pre-processing
  * completed
  * Checking number of files
          `$ ls *_1*.filter.gz | wc -l #18 `
          `$ ls ./Rubio_Raw_Transcriptomes/*_1*fastq.gz | wc -l #18`
  * compressing raw data
          `$ tar -zcvf Rubio_Raw_Transcriptomes.archive.tar.gz Rubio_Raw_Transcriptomes/`

* Starting Probiotic pre-processing
  * completed
  * correct number of files present
  * Compressing raw data
        `$ tar -zcvf Probiotic_Raw_Transcriptomes.archive.tar.gz Probiotic_Raw_Transcriptomes/`


* All scripts have unique output and error files. Ended up running pre-processing on data subsets one at a time in order to expedite process.

## 2/4/2020 - 2/5/2020 Building HISAT2 Script and Running HISAT2

### 1. Building new HISAT index.

* Note, a new version of HISAT is now available on the cluster. I'm going to use this version (HISAT2/2.1.0-foss-2016b)
  *  Note also that the code I'm adding to github is all combined into one script, though I ran the same   code in the cluster separated into multiple scripts to aid in running multiple scripts at once. The code used however was not changed.
* Checking how HISAT2 genome indexes we made. Creating new indexes for both C_vir and C_gig using updated software version

  * Checking C_virginica genome file on hand which was edited from `ref_C_virginica-3.0_top_level.gff3` to remove a space from header that was making it incompatible with Stringtie annotation. Checking also that mitochondrial DNA was added.

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

    * Checking C_gigas genome file. Only on the cluster was `Crassostrea_gigas.gff`. A newer version of C_gigas genome is avaiable on NCBI now. Downloading that version and will use it to create HISAT2 index for C_gigas sequences.
    * The C. gigas genome file contains all the genomic scaffolds (it has never been assembled to chromosome level like oysters) and the mitochondrial genome.

* Creating new indexes using code in `02_Build_Hisat_Indexes.sh`. Indices created!

### 2. Creating script to perform HISAT 2 alignment and SAMtools sorting on files for each experiment.

* Dermo transcriptomes, He OsHV1, Zhang Vibrio and ROD transcriptomes are currently finished being preprocessed. Starting HISAT on these transcriptomes. Checked script to make sure correct indexes and paths were used.

* Received error running the script: Checking with kevin about how to resolve conflict

        `foss/2016b(13):ERROR:150: Module 'foss/2016b' conflicts with the currently loaded module(s) 'foss/2018b'
         foss/2016b(13):ERROR:102: Tcl command execution failed: conflict foss

         GCCcore/5.4.0(13):ERROR:150: Module 'GCCcore/5.4.0' conflicts with the currently loaded module(s) 'GCCcore/7.3.0'
         GCCcore/5.4.0(13):ERROR:102: Tcl command execution failed: conflict GCCcore`

* Resolved error. Cannot download simultaneously packages that have different foss toolchains. Kevin downloaded the HISAT2 foss-2018b so that is is compatible with the SAMtools foss-2018b.

* Dermo transcriptomes finished HISAT2
  * Compressing .sam files and putting HISAT2 stats in a separate folder. There are the same number of .sam and .bam files.

* Ran all other transcriptomes in script. Checking through error and output files to make sure all worked correctly. Checking correct number of sorted bam files for each before Stringtie.

* ROD, Dermo, HE, Zhang completed HISAT
  * Deleting SAM files from each
  * deleting clean.trim.filter files also after checking that correct bam files present
    * ROD - deleted filter files
  `$ ls *.bam | wc -l #12`
  `$ ls *.filter.gz | wc -l # 12`
   * Dermo - deleted filter files
   `$ ls *.bam | wc -l`
   `$ ls *.filter.gz | wc -l # 194`
   * He - deleted filter files
   `$ ls *filter.gz | wc -l # 32 `
   `$ ls *bam | wc -l # 32`
   * Zhang - deleted filter files
    `$ ls *filter.gz | wc -l # 9`
    `$ ls *.bam | wc -l # 9`

* Probiotic, deLorgeril, and Rubio STOPPED with error and stopped half way through deLorgeril:
      `File write failed
samtools sort: failed to create temporary file "./samtools.15357.2415.tmp.0007.bam": No space left on device
/var/spool/slurmd/job1035890/slurm_script: line 52: /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_deLorgeril_OsHV1_SRA/SRR6679068_1.fastq.gz.clean.trim.filter.gz.bam.stats: Disk quota exceeded
`
    * Deleting the bam files already done for deLorgeril because I'm not sure where they stopped
    * Deleting sam files for those experiments that have already finished running
* Rerunning Probiotic, deLorgeril and Rubio HISAT2 script but combining script with Stringtie at the end so that it will run over the weekend. See how the Stringtie script was developed below

## 2/7/2020 Build Stringtie Script and run

### 1. Setting up stringtie script

* After HISAT2 there are no longer two files for paired end reads, so I can run the same script for all.
* Creating mergelist for each experiment
* Ran script `03_Stringtie_Assembly_Quantify.sh` with the ROD_Dermo_HE_Zhang all together
* Ran script for performing HISAT mapping and Stringtie quantification for Rubio, Probiotic, and deLorgeril all together `02_HISAT2_samtools_sort_Rubio_Pro_deLorgeril.sh`

## 2/10/2020 Checking Output of Stringtie to check that it worked

* Checking output and error files resulting from the scripts above
      `$ cat  Stringtie_ROD_Dermo_HE_Zhang_out_2_7_2020`
      * Dermo DONE : assembled, then merged, then gffcompared, then transcript abundance re-estimated. BUT DIDN'T WORK DUE TO MERGELIST.TXT ERROR.
      * HE DONE : assembled, then merged, then gffcompared, then transcript abundance re-estimated. BUT DIDN'T WORK DUE TO MERGELIST.TXT ERROR.
      * Zhang DONE, with assembly, gffcompare, merge, and re-abundnace but it seems like then the script re-assembled one of the files and then stopped before moving on to the ROD files.
      * ROD not merged at all.
`
      * Looking at the Stringtie_ROD_Dermo_HE_Zhang_error_2_7_2020 file. Shows that the mergelist.txt files were not created correctly so no transcripts were found and no valid reference transcripts were found either. Need to fix merglist and re-try the merging step. Also need to check why a syntax error occurred for array 7 on line 185.

      `Error: no transcripts were found in input file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/2020_Raw_Transcriptome_Data/Dermo_2015_transcriptomes/Fastq/Dermo_mergelist.txt
  67891 reference transcripts loaded.
  38 duplicate reference transcripts discarded.
  0 query transfrags loaded.
Error: could not any valid reference transcripts in /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/2020_Raw_Transcriptome_Data/Dermo_2015_transcriptomes/Fastq/Dermo_stringtie_merged.gtf (invalid GTF/GFF file?)
Error: could not any valid reference transcripts in /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/2020_Raw_Transcriptome_Data/Dermo_2015_transcriptomes/Fastq/Dermo_stringtie_merged.gtf (invalid GTF/GFF file?)

/var/spool/slurmd/job1036269/slurm_script: line 185: syntax error near unexpected token `;'
/var/spool/slurmd/job1036269/slurm_script: line 185: `for i in ${array7[@]}; do'
WARNING: no reference transcripts were found for the genomic sequences where reads were mapped!
Please make sure the -G annotation file uses the same naming convention for the genome sequences.
`
* Syntax error due to accidentally deleting closeout parenthesis when creating array7 :`array7=($(ls $CR/*.bam)`.

* Same Stringtie merge error due to the mergelisted happen for the deLorgeril, Rubio, and Probiotic file (`HISAT2_Stringtie_deLorgeril_Rubio_Pro_error_2_7_2020`). Though the HISAT2 alignment finished correctly. Need to redo Stringtie mergelist onward for these transcriptomes as well.
        ` $ cat HISAT2_Stringtie_deLorgeril_Rubio_Pro_error_2_7_2020` # this files holds the alignment rate statistics for the HISAT mapping  

* Before fixing Stringtie mergelist error, going to delete the sam files and trimmed and filtered files from deLorgeril_Rubio_Pro that finished alignment with HISAT2.
        `# Checking number of sam files and bam files to see that they match
        $ pwd  /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_deLorgeril_OsHV1_SRA
        $ ls *.bam | wc -l # 42
        $ ls *.sam | wc -l # 42
        # Checking out a Stringtie mapping file before deleting .sam in case I get the old error I used to get where the annotation file was not compatible with the original genome headers that the mapping was done with
        $ less SRR6679052_1.fastq.gz.clean.trim.filter.gz.bam.gtf # This shows that gene_ids from the transcriptome file have actually been mapped to LOC ids and XM ids from the annotation. Great! Now I can delete these and not have to worry about re-doing
        $ rm *sam
        # Deleting trim.filter.gz files as well
        $ ls *.filter.gz | wc -l # 84
        # What to do with all the tmp files? Am going to keep for now.

        # Repeating process above for Rubio files
        $ pwd
        /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA
        $ ls *.sam | wc -l # 18
        $ ls *.bam | wc -l # 18
        $ ls *.filter.gz | wc -l # 36
        $ less SRR8551085_1.fastq.gz.clean.trim.filter.gz.bam.gtf # mapped to LOC and XM Id's
        $ rm *.filter.gz
        $ rm *.sam
        # keeping tmp.* files for now

        # Repeating process for Probiotic files
        $ pwd
        /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/2020_Raw_Transcriptome_Data/C_vir_Probiotic_SRA
        $ ls *.sam | wc -l # 6
        $ ls *.bam | wc -l # 6
        $ ls *.filter.gz | wc -l # 12
        $ less SRR5357622_1.fastq.gz.clean.trim.filter.gz.bam.gtf # also was mapped correctly to ref name LOC
        $ rm *.sam
        $ rm *.filter.gz
        # keeping the tmp.* files around
`
* Deleting old gtf files for Zhang_Vibrio because one of them was accidentally overwritten. Will re-do this step
        ` $ ls *.bam.gtf | wc -l # 9
          $ ls *.bam | wc -l # 9
        `
* ROD and Zhang both need the initial assembly steps re-done, while the rest just need to have the merge step fixed and re-abundance re-calculated.

* Inspecting why the Stringtie code failed.

      1. Inspect the mergelist.txt created in each script to check format. Appears to be in correct format with the full path to each file, one file per line.
        `]$ head deLorgeril_mergelist.txt
/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_deLorgeril_OsHV1_SRA/SRR6679052_1.fastq.gz.clean.trim.filter.gz.bam.gtf
/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_deLorgeril_OsHV1_SRA/SRR6679053_1.fastq.gz.clean.trim.filter.gz.bam.gtf
/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_deLorgeril_OsHV1_SRA/SRR6679054_1.fastq.gz.clean.trim.filter.gz.bam.gtf
/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_deLorgeril_OsHV1_SRA/SRR6679055_1.fastq.gz.clean.trim.filter.gz.bam.gtf
/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_deLorgeril_OsHV1_SRA/SRR6679056_1.fastq.gz.clean.trim.filter.gz.bam.gtf
/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_deLorgeril_OsHV1_SRA/SRR6679057_1.fastq.gz.clean.trim.filter.gz.bam.gtf
/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_deLorgeril_OsHV1_SRA/SRR6679058_1.fastq.gz.clean.trim.filter.gz.bam.gtf
/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_deLorgeril_OsHV1_SRA/SRR6679059_1.fastq.gz.clean.trim.filter.gz.bam.gtf
/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_deLorgeril_OsHV1_SRA/SRR6679060_1.fastq.gz.clean.trim.filter.gz.bam.gtf
/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_deLorgeril_OsHV1_SRA/SRR6679061_1.fastq.gz.clean.trim.filter.gz.bam.gtf
[erin_roberts@bluewaves C_gig_deLorgeril_OsHV1_SRA]$`
        2. Inspect output of the merging step. Says no transcripts found in file.
          `Error: no transcripts were found in input file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_deLorgeril_OsHV1_SRA/deLorgeril_mergelist.txt
  53712 reference transcripts loaded.
  22 duplicate reference transcripts discarded.
  0 query transfrags loaded.
  Error: could not any valid reference transcripts in /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_deLorgeril_OsHV1_SRA/deLorgeril_OsHV1_stringtie_merged.gtf (invalid GTF/GFF file?)
`
        3. Deleting the merge step files from deLorgeril OsHV1 and trying to run in an interactive session in the folder where the files are (to answer if this is an issue with the path)
         `$ rm *merged*
          $ module load StringTie/2.1.1-GCCcore-7.3.0
          $ module load gffcompare/0.11.5-foss-2018b
          $ CG=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/Cgig_Genome_and_Indexes
          $ GLO=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_deLorgeril_OsHV1_SRA
          # tested whether not loading gffcompare and module load StringTie/2.1.1-GCCcore-7.3.0  made it difference. It didn't
          # tested whether trying older stringtie version made a difference
          $ module purge
          $ module load StringTie/1.3.3b-foss-2016b
          $ stringtie --merge -A -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o $GLO/deLorgeril_OsHV1_stringtie_merged.gtf deLorgeril_mergelist.txt
          # Error: no transcripts were found in input file deLorgeril_mergelist.txt
          # Trying to test if making new mergelist makes a difference
          $ ls S*.gtf > deLorgeril_mergelist_no_path.txt
          $ cat deLorgeril_mergelist_no_path.txt
          $ module purge
          $ module load StringTie/2.1.1-GCCcore-7.3.0
          $ stringtie --merge -A -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o $GLO/deLorgeril_OsHV1_stringtie_merged.gtf deLorgeril_mergelist_no_path.txt
              Error: no transcripts were found in input file deLorgeril_mergelist_no_path.txt
          # Changing mergelist doesnt make a difference
          # Testing removing -A from the command
          $ stringtie --merge -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o $GLO/deLorgeril_OsHV1_stringtie_merged.gtf deLorgeril_mergelist_no_path.txt
              # doesn't immediately give error that no transcripts were found
          ####### THE ISSUE WAS USING -A DURING THE STRINGTIE MERGE THIS IS NOT ALLOWED #######
          `
* Changing Stringtie script so that the -A option is only used during the re-estimating abundances step. Keeping initial assembly step only for ROD and Zhang.
  * Ran fixed code in the following script on bluewaves "03_Stringtie_Assembly_Quantify_fixed_redo.sh " though it was just fixed and saved in my original file on my computer as "03_Stringtie_Assembly_Quantify.sh"

  * Ran script and got the following error
            `Error: invalid -e usage, GFF reference not given (-G option required).`
  * Testing if switching the order of -e and -A when doing the transcript abundance re-estimation makes a difference. YES THIS FIXED THE ERROR

## 2/11/2020 Checking Stringtie output and converting into DESeq2 format count tables

* Checking script error file `cat Stringtie_ALL_error_2_10_2020` and output file `cat Stringtie_ALL_out_2_10_2020`
  * Error output file just lists the merging statistics for each of the experiments. The statistics look good.
  * Output file
    - deLorgeril completed all steps
    - Rubio Vibrio completed all steps
    - Probiotic completed all steps
    - Dermo completed all steps
    - He OsHV1 completed all steps
    - Zhang Vibrio completed initial stringtie assembly and all subsequent steps
    - ROD completed initial stringtie assembly and all subsequent steps
* Checking the output in each folder
        `# Checking deLorgeril Stringtie output
        $ pwd /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_deLorgeril_OsHV1_SRA
        $ ls *.merge.gtf | wc -l # 42
        $ ls *.bam.gtf | wc -l # 42
        $ mkdir deLorgeril_gffcompare # folder houses all the gffcompare output files
        $ mkdir deLorgeril_HISAT_bam
        $ mkdir deLorgeril_Stringtie_gtf
        # Compress HISAT2 bam files
        $ tar -zcvf deLorgeril_HISAT_bam.archive.tar.gz deLorgeril_HISAT_bam/

        # Checking Rubio Vibrio output
        $ pwd
        /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA
        $ ls *.bam.gtf | wc -l # 18
        $ ls *.merge.gtf | wc -l # 18
        $ mkdir Rubio_Stringtie_gtf
        $ mkdir Rubio_gffcompare
        $ mkdir Rubio_HISAT_bam
        # Compress HISAT2 bam files
        $ tar -zcvf Rubio_HISAT_bam.archive.tar.gz Rubio_HISAT_bam/

        # Checking Probiotic Files
        $ pwd
        /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/2020_Raw_Transcriptome_Data/C_vir_Probiotic_SRA
        $ ls *.merge.gtf | wc -l # 6
        $ ls *.bam.gtf | wc -l # 6
        $ mkdir Probiotic_HISAT_bam
        $ mkdir Probiotic_Stringtie_gtf
        $ mkdir Probiotic_gffcompare
        # Compress HISAT2 bam files
        $ tar -zcvf Probiotic_HISAT_bam.archive.tar.gz Probiotic_HISAT_bam

        # Checking Dermo files
        $ pwd
        /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/2020_Raw_Transcriptome_Data/Dermo_2015_transcriptomes/Fastq
        $ ls *.merge.gtf | wc -l # 97
        $ ls *.bam.gtf | wc -l # 97
        $ mkdir Dermo_gffcompare
        $ mkdir Dermo_HISAT_bam
        $ mkdir Dermo_Stringtie_gtf
        # will compress later

        # Checking He OsHV1
        $ pwd
        /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_He_OsHV1_SRA
        $ ls *merge.gtf | wc -l # 32
        $ ls *bam.gtf | wc -l # 32
        $ mkdir He_gffcompare
        $ mkdir He_HISAT_bam
        $ mkdir He_Stringtie_gtf
        # Will compress later

        # Checking Zhang Vibrio
        $ pwd
        /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Zhang_Vibrio_SRA
        $ ls *.merge.gtf | wc -l # 9
        $ ls *.bam.gtf | wc -l # 9
        $ mkdir Zhang_gffcompare
        $ mkdir Zhang_HISAT_bam
        $ mkdir Zhang_Stringtie_gtf
        # will compress later

        # Checking ROD Vibrio
        $ pwd
        /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/2020_Raw_Transcriptome_Data/C_vir_ROD_SRA
        $ ls *.merge.gtf | wc -l # 12
        $ ls *.bam.gtf | wc -l # 12
        $ mkdir ROD_gffcompare
        $ mkdir ROD_Stringtie_gtf
        $ mkdir ROD_HISAT_bam
        # Compress HISAT BAM
        $ tar -zcvf ROD_HISAT_bam.archive.tar.gz ROD_HISAT_bam
        `
* Creating `04_Prep_Stringtie_DESeq2.sh`
  * had to modify the existing code for creating the sample list files to take into account the underscore in paired end file names. NVM keeping same as original to preserve the underscore
  * Running the script gives the following error in the error file:
        `Traceback (most recent call last):
         File "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2017_MISC_SCRIPTS/prepDE.py", line 254, in <module>
         geneDict[geneIDs[i]][s[0]]+=v[s[0]]
         KeyError: 'SRR6679053'
         Traceback (most recent call last):
         File "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2017_MISC_SCRIPTS/prepDE.py", line 254, in <module>
         geneDict[geneIDs[i]][s[0]]+=v[s[0]]
         KeyError: 'SRR6679053'
         `
  * The output file lists a 1 or 0 next to each sample ID
        `Start Tue Feb 11 12:53:29 EST 2020
        0 SRR6679052
        1 SRR6679053
        0 SRR6679052
        1 SRR6679053
        0 SRR6679052
        1 SRR6679053
        0 SRR6679052
        1 SRR6679053
        0 SRR6679052
        1 SRR6679053
        0 SRR6679052
        1 SRR6679053`
  * Added into the script loading the python module. This didn't fix issue though I think this is necessary
  * Why does SRR6679052 have a 0 and SRR6679053 have a 1?
        `$ less SRR6679052_1.merge.gtf
        # Output shows that transcripts by default are given MSTRG headers rather than the -1 headers with the file name I wanted them to have
        $ less SRR6679053_1.merge.gtf
        # this shows the same thing
        `
        * trying to remove the -s option from prepDE.py running to see if this makes a difference
        * Still received the same error
  * Trying to run the prepDE script just on the command line in deLorgeril folder
        `$ python $P/prepDE.py -i C_gig_deLorgeril_sample_list.txt -g deLorgeril_gene_count_matrix.csv -t deLorgeril_transcript_count_matrix.csv
          0 SRR6679052
          1 SRR6679052
          2 SRR6679053
          Traceback (most recent call last):
          File "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2017_MISC_SCRIPTS/prepDE.py", line 254, in <module>
          geneDict[geneIDs[i]][s[0]]+=v[s[0]]
          KeyError: 'SRR6679053'
          `
  * Is the problem may be that I used stringtie -l when I first aligned, but not in the re-estimation? I'm going to test this on one file from Rubio, SRR8551076_1.fastq.gz.clean.trim.filter.gz.bam, SRR8551077_1.fastq.gz.clean.trim.filter.gz.bam. They do this in the Pertea et al. 2016 tutorial (use -l the first time but not the second time..)
            `$ interactive
            $ pwd
            /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/Rubio_HISAT_bam
            $ CG=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/Cgig_Genome_and_Indexes
            $ GRV=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA
            $ module load StringTie/2.1.1-GCCcore-7.3.0
            $ stringtie -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o SRR8551076_1.test.gtf SRR8551076_1.fastq.gz.clean.trim.filter.gz.bam
            $ stringtie -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o SRR8551077_1.test.gtf SRR8551077_1.fastq.gz.clean.trim.filter.gz.bam
            $ stringtie --merge -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o TEST_stringtie_merged.gtf TEST_mergelist.txt
            $ stringtie -A -e -G TEST_stringtie_merged.gtf -o SRR8551076_1.merge.test.gtf SRR8551076_1.fastq.gz.clean.trim.filter.gz.bam
            $ stringtie -A -e -G TEST_stringtie_merged.gtf -o SRR8551077_1.merge.test.gtf SRR8551077_1.fastq.gz.clean.trim.filter.gz.bam
            $ for i in *.merge.gtf; do
            	# create text file with sample IDs and respective paths
            	echo "$(echo $i |sed "s/\..*//") $GRV/$i" >> test_sample_list.txt
            done
           $ cp /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2017_MISC_SCRIPTS/prepDE.py .

            $ P=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2017_MISC_SCRIPTS
            $ python prepDE.py -i test_sample_list.txt -g test_gene_count_matrix.csv -t test_transcript_count_matrix.csv
            $ python prepDE.py -i test_sample_list.txt -g test_gene_count_matrix.csv -t test_transcript_count_matrix.csv
              0 SRR8551076_1
              1 SRR8551077_1
              Traceback (most recent call last):
              File "prepDE.py", line 252, in <module>
              geneDict.setdefault(geneIDs[i],{}) #gene_id
              KeyError: 'STRG.36898.2'
              # Didn't fix the error running Stringtie without -l didn't make a difference (which is good meaning I don't have to re-run all files through Stringtie)
            # Trying to run without -g
            $ python prepDE.py -i test_sample_list.txt  -t test_transcript_count_matrix.csv # got the same Error

            # trying to redownload prepDE.py to make sure something didn't get messed up in the file prepDE_new.py
            $ python prepDE_new.py -i test_sample_list.txt -t test_transcript_count_matrix.csv
            0 SRR8551076_1
            1 SRR8551077_1
            Traceback (most recent call last):
            File "prepDE_new.py", line 255, in <module>
            geneDict.setdefault(geneIDs[i],{}) #gene_id
            KeyError: 'STRG.36898.2'

            `
    * Noticed that I was running the python prepDE.py as a loop and this is not allowed. Testing removing the loop to see what happens
            - Still got the same error, but now with MSTRG also giving the same error
            `Traceback (most recent call last):
            File "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2017_MISC_SCRIPTS/prepDE.py", line 254, in <module>
            geneDict[geneIDs[i]][s[0]]+=v[s[0]]
            KeyError: 'SRR6679053'
            Traceback (most recent call last):
            File "/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2017_MISC_SCRIPTS/prepDE.py", line 252, in <module>
            geneDict.setdefault(geneIDs[i],{}) #gene_id
            KeyError: 'STRG.36898.2'
              `
  * Found online a Github error post from the developer ppertea, there was an issue with the software but it was fixed by the 2.0.4 V release https://github.com/gpertea/stringtie/issues/238
      *MUST READ THIS POST THIS IS MY EXACT ERROR*
      *also read his second post with the udpated prepde.py script as of Oct. 2019:https://github.com/gpertea/stringtie/issues/234#issuecomment-541494630*
      - Reminder from this second issue 234 from Pertea:
      "To reiterate and clarify: prepDE.py can only be used on a set of stringtie GTF outputs if stringtie was run, for all those outputs:
          - with the -e option
          - with the same file for the -G option.
          Also, make sure that no other GTF files (like the reference annotation file) are present in those sub-directories, only the stringtie output GTF files should be found there, as the default mode of operation for prepDE is to scan all the sub-directories there for .gtf files which are all expected to have been produced by stringtie by following the requirements above (-e option, same -G file)."
      - Link to the latest prepDE.py :https://raw.githubusercontent.com/gpertea/stringtie/master/prepDE.py, need to use with adding the -v option

  * Reran script with the new version of the PrepDE.py script (prepDE_Oct.2019.py)
          `$ python prepDE_Oct.2019.py -v -i test_sample_list.txt -t test_transcript_count_matrix.csv
          >processing sample SRR8551076_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/SRR8551076_1.merge.gtf
          >processing sample SRR8551077_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/SRR8551077_1.merge.gtf
          Error: could not locate transcript STRG.36898.2 entry for sample SRR8551077_1
          Traceback (most recent call last):
          File "prepDE_Oct.2019.py", line 281, in <module>
          geneDict.setdefault(geneIDs[i],{}) #gene_id
            KeyError: 'STRG.36898.2'
        `
  * Even though this error should be fixed on the V 2.1.1 release, trying with the StringTie/1.3.3b-foss-2016b to make sure, without the -l option used
          `$ module load StringTie/1.3.3b-foss-2016b
           $ stringtie -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o SRR8551076_1.test.gtf SRR8551076_1.fastq.gz.clean.trim.filter.gz.bam
           $ stringtie -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o SRR8551077_1.test.gtf SRR8551077_1.fastq.gz.clean.trim.filter.gz.bam
           $ stringtie --merge -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o TEST_stringtie_merged.gtf TEST_mergelist.txt
           $ stringtie -A gene_abund_table -e -G TEST_stringtie_merged.gtf -o SRR8551076_1.merge.test.gtf SRR8551076_1.fastq.gz.clean.trim.filter.gz.bam
           $ stringtie -A -e -G TEST_stringtie_merged.gtf -o SRR8551077_1.merge.test.gtf SRR8551077_1.fastq.gz.clean.trim.filter.gz.bam
           $ for i in *.merge.gtf; do
             # create text file with sample IDs and respective paths
             echo "$(echo $i |sed "s/\..*//") $GRV/$i" >> test_sample_list.txt
           done
           $ python prepDE_Oct.2019.py -i test_sample_list.txt -g test_gene_count_matrix.csv -t test_transcript_count_matrix.csv
           #STOPPING THIS BEFORE FINISHING BECAUSE I BELIEVE THE BULLET POINT BELOW WAS THE ISSUE, GOING TO RERUN WITH NEW VERSION ALL OVER AGAIN TO CHECK
           `
  * POSSIBLY FOUND THE ISSUE. I DIDN'T PUT THE NAME OF THE GENE TABLE AFTER THE -A OPTION IN MY STRINGTIE SCRIPT
          ` pwd
          /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/Rubio_HISAT_bam
          $ module load StringTie/2.1.1-GCCcore-7.3.0
          $ CG=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/Cgig_Genome_and_Indexes
          $ GRV=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/Rubio_HISAT_bam
          $ stringtie -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o SRR8551076_1.test.gtf SRR8551076_1.fastq.gz.clean.trim.filter.gz.bam
          $ stringtie -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o SRR8551077_1.test.gtf SRR8551077_1.fastq.gz.clean.trim.filter.gz.bam
          $ stringtie --merge -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o TEST_stringtie_merged.gtf TEST_mergelist.txt
          $ stringtie -A SRR8551076_1_gene_abd.tab -e -G TEST_stringtie_merged.gtf -o SRR8551076_1.merge.gtf SRR8551076_1.fastq.gz.clean.trim.filter.gz.bam
          $ stringtie -A SRR8551077_1_gene_abd.tab -e -G TEST_stringtie_merged.gtf -o SRR8551077_1.merge.gtf SRR8551077_1.fastq.gz.clean.trim.filter.gz.bam
          $ for i in *.merge.gtf; do
            # create text file with sample IDs and respective paths
            echo "$(echo $i |sed "s/\..*//") $GRV/$i" >> test_sample_list.txt
          done
          $ module load python/2.7.6
          $ python prepDE_Oct.2019.py -i test_sample_list.txt -g test_gene_count_matrix.csv -t test_transcript_count_matrix.csv
          Error: could not locate transcript STRG.36898.2 entry for sample SRR8551077_1
          Traceback (most recent call last):
          File "prepDE_Oct.2019.py", line 281, in <module>
          geneDict.setdefault(geneIDs[i],{}) #gene_id
          KeyError: 'STRG.36898.2'  
          `
  * Nope! Though I will need to fix the -A option and rerun stringtie, removing -A and -l did not fix the issue. Trying again with the older version of Stringtie
          `$ module purge
          $ module load StringTie/1.3.3b-foss-2016b
          $ module load python/2.7.6
          $ stringtie -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o SRR8551076_1.test.gtf SRR8551076_1.fastq.gz.clean.trim.filter.gz.bam
          $ stringtie -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o SRR8551077_1.test.gtf SRR8551077_1.fastq.gz.clean.trim.filter.gz.bam
          $ stringtie --merge -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o TEST_stringtie_merged.gtf TEST_mergelist.txt
          $ stringtie -A SRR8551076_1_gene_abd.tab -e -G TEST_stringtie_merged.gtf -o SRR8551076_1.merge.gtf SRR8551076_1.fastq.gz.clean.trim.filter.gz.bam
          $ stringtie -A SRR8551077_1_gene_abd.tab -e -G TEST_stringtie_merged.gtf -o SRR8551077_1.merge.gtf SRR8551077_1.fastq.gz.clean.trim.filter.gz.bam
          $ for i in *.merge.gtf; do
            # create text file with sample IDs and respective paths
            echo "$(echo $i |sed "s/\..*//") $GRV/$i" >> test_sample_list.txt
          done
          $ python prepDE_Oct.2019.py -i test_sample_list.txt -g test_gene_count_matrix.csv -t test_transcript_count_matrix.csv
          Error: could not locate transcript STRG.36898.2 entry for sample SRR8551077_1
          Traceback (most recent call last):
          File "prepDE_Oct.2019.py", line 281, in <module>
          geneDict.setdefault(geneIDs[i],{}) #gene_id
          KeyError: 'STRG.36898.2'
          ## Got the same issue even with using the old version of Stringtie... not sure what to do now
          `

  * Getting Kevin to download the V 2.0.4 to see if in this supposedly fixed verion the issue stops happening. Is there an issue in how I ran HISAT2?

## 2/12/2020 Bug Checking Stringtie and (hopefully) formatting output for DESeq2

* Kevin installed v2.0.4 so I can compare the output of v1.3.3 and v2.1.0 (which both have the error) with v 2.0.4 to see if I keep getting the same output where STRG values are mis-assigned. Going to stay working in the Rubio folder with those 2 transcriptomes since that has been my test case so far.
  * While this is running I'm still deleting the python loop and the -s in the python line I created in `04_Prep_Stringtie_DESeq2.sh`. Also added in loading the python module at the beginning. Specifying also the new python prepDE.py script made by Pertea in Oct. 2019, called in bluewaves prepDE_Oct.2019.py. Putting a copy of this script in each directory I will run it in to make sure that python correctly "sees" it.
  * In the Stringtie script, `03_Stringtie_Assembly_Quantify.sh` now I'm removing the -l option in the initial assembly and re-estimation and making sure I put in the name of the output gene abundance files after the -A option in the re-estimation. The merge step and gffcompare steps are perfect do not need to re-do.  
        `$ pwd
        /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/Rubio_HISAT_bam
        $ module load python/2.7.6
        $ module load  StringTie/2.0.4-GCCcore-7.3.0
        $ CG=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/Cgig_Genome_and_Indexes
        $ stringtie -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o SRR8551076_1.test.gtf SRR8551076_1.fastq.gz.clean.trim.filter.gz.bam
        $ stringtie -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o SRR8551077_1.test.gtf SRR8551077_1.fastq.gz.clean.trim.filter.gz.bam
        $ ls *test.gtf >> $GRV/TEST_mergelist.txt
        $ stringtie --merge -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o TEST_stringtie_merged.gtf TEST_mergelist.txt
        $ stringtie -A SRR8551076_1_gene_abd.tab -e -G TEST_stringtie_merged.gtf -o SRR8551076_1.merge.gtf SRR8551076_1.fastq.gz.clean.trim.filter.gz.bam
        $ stringtie -A SRR8551077_1_gene_abd.tab -e -G TEST_stringtie_merged.gtf -o SRR8551077_1.merge.gtf SRR8551077_1.fastq.gz.clean.trim.filter.gz.bam
        GRV=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/Rubio_HISAT_bam
        $for i in *.merge.gtf; do
          # create text file with sample IDs and respective paths
          echo "$(echo $i |sed "s/\..*//") $GRV/$i" >> test_sample_list.txt
        done
        $ python prepDE_Oct.2019.py -v -i test_sample_list.txt -g test_gene_count_matrix.csv -t test_transcript_count_matrix.csv
        >processing sample SRR8551076_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/Rubio_HISAT_bam/SRR8551076_1.merge.gtf
        >processing sample SRR8551077_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/Rubio_HISAT_bam/SRR8551077_1.merge.gtf
        ..writing test_transcript_count_matrix.csv
        ..writing test_gene_count_matrix.csv
        All done.
        ###### OMG V 2.0.4 WORKED ######
        `
* I think I may have set my $GRV to the wrong path last time when making my mergelist.txt file when running 2.1.0 and 2.1.3. Going to run these versions one more time here to check error before commenting to pertea et al.
                `$ pwd
                /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/Rubio_HISAT_bam
                $ module purge
                $ module load python/2.7.6
                $ module load StringTie/2.1.1-GCCcore-7.3.0
                $ CG=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/Cgig_Genome_and_Indexes
                $ stringtie -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o SRR8551076_1.test.gtf SRR8551076_1.fastq.gz.clean.trim.filter.gz.bam
                $ stringtie -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o SRR8551077_1.test.gtf SRR8551077_1.fastq.gz.clean.trim.filter.gz.bam
                $ stringtie --merge -G $CG/GCF_000297895.1_oyster_v9_genomic.gff -o TEST_stringtie_merged.gtf TEST_mergelist.txt
                $ stringtie -A SRR8551076_1_gene_abd.tab -e -G TEST_stringtie_merged.gtf -o SRR8551076_1.merge.gtf SRR8551076_1.fastq.gz.clean.trim.filter.gz.bam
                $ stringtie -A SRR8551077_1_gene_abd.tab -e -G TEST_stringtie_merged.gtf -o SRR8551077_1.merge.gtf SRR8551077_1.fastq.gz.clean.trim.filter.gz.bam
                GRV=/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/Rubio_HISAT_bam
                $for i in *.merge.gtf; do
                  # create text file with sample IDs and respective paths
                  echo "$(echo $i |sed "s/\..*//") $GRV/$i" >> test_sample_list.txt
                done
                $ python prepDE_Oct.2019.py -v -i test_sample_list.txt -g test_gene_count_matrix.csv -t test_transcript_count_matrix.csv
                >processing sample SRR8551076_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/Rubio_HISAT_bam/SRR8551076_1.merge.gtf
                >processing sample SRR8551077_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/Rubio_HISAT_bam/SRR8551077_1.merge.gtf
                ..writing test_transcript_count_matrix.csv
                ..writing test_gene_count_matrix.csv
                All done.
                ### WHAT??? THIS WORKED TOO!!!!! YAY ####
              `
* *I have figured out now that using V 2.1.0 works as long as I keep the -l out of it and appropriately put a name after the -A option in the re-estimation.* It may have also helped that have the python script in the same directly. Now I'm going to untar any of the bam folders from HISAT that I tarred yesterday and re-run Stringtie with my updated script. I'm relieved that I don't have to re-run HISAT! Running the Stringtie script separately for each population in the hopes that this will help it run faster.
  * Copying python script into all necessary directories
          `$ cp ./C_gig_Rubio_Vibrio_SRA/prepDE_Oct.2019.py ./C_gig_Zhang_Vibrio_SRA/
           $ cp ./C_gig_Rubio_Vibrio_SRA/prepDE_Oct.2019.py ./C_gig_He_OsHV1_SRA/
           $ cp ./C_gig_Rubio_Vibrio_SRA/prepDE_Oct.2019.py ./C_gig_deLorgeril_OsHV1_SRA/
           $ cp ./C_gig_Rubio_Vibrio_SRA/prepDE_Oct.2019.py ../../C_Vir_subset/2020_Raw_Transcriptome_Data/C_vir_ROD_SRA/
           $ cp ./C_gig_Rubio_Vibrio_SRA/prepDE_Oct.2019.py ../../C_Vir_subset/2020_Raw_Transcriptome_Data/C_vir_Probiotic_SRA/
           $ cp ./C_gig_Rubio_Vibrio_SRA/prepDE_Oct.2019.py ../../C_Vir_subset/2020_Raw_Transcriptome_Data/Dermo_2015_transcriptomes/Fastq/
           `
  * Deleting all `*merge*` files and `*.gtf` files from every folder (including mergelist.txt, all merged gtfs, all initial assembly gtfs), and deleting the original gffcompare files because this is going to be run again. Moving all `*.bam` files back out into main directory for each.
  * untaring any tarred bam file directories
          `$ tar -zxvf deLorgeril_HISAT_bam.archive.tar.gz
           $ tar -zxvf Probiotic_HISAT_bam.archive.tar.gz
           $ tar -zxvf ROD_HISAT_bam.archive.tar.gz
           `
  * Copied separate Stringtie scripts made locally to the cluster
          `$ scp 03_Stringtie_Assembl*Quantify_* erin_roberts@bluewaves:/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/2020_Scripts`
  * Running the Rubio files again first as a test case to make sure my full code works before running all scripts
    * *IT WORKED YAY*
  * RE-running Stringtie and then prepDE_Oct.2019.py for each
    1. Rubio Stringtie complete and prep for DESeq2 complete
          `$  cat Stringtie_Rubio_prepDE_out_2_11_2020
START Wed Feb 12 13:15:59 EST 2020
>processing sample SRR8551076_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/SRR8551076_1.merge.gtf
>processing sample SRR8551077_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/SRR8551077_1.merge.gtf
>processing sample SRR8551078_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/SRR8551078_1.merge.gtf
>processing sample SRR8551079_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/SRR8551079_1.merge.gtf
>processing sample SRR8551080_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/SRR8551080_1.merge.gtf
>processing sample SRR8551081_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/SRR8551081_1.merge.gtf
>processing sample SRR8551082_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/SRR8551082_1.merge.gtf
>processing sample SRR8551083_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/SRR8551083_1.merge.gtf
>processing sample SRR8551084_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/SRR8551084_1.merge.gtf
>processing sample SRR8551085_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/SRR8551085_1.merge.gtf
>processing sample SRR8551086_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/SRR8551086_1.merge.gtf
>processing sample SRR8551087_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/SRR8551087_1.merge.gtf
>processing sample SRR8551088_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/SRR8551088_1.merge.gtf
>processing sample SRR8551089_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/SRR8551089_1.merge.gtf
>processing sample SRR8551090_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/SRR8551090_1.merge.gtf
>processing sample SRR8551091_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/SRR8551091_1.merge.gtf
>processing sample SRR8551092_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/SRR8551092_1.merge.gtf
>processing sample SRR8551093_1 from file /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_Rubio_Vibrio_SRA/SRR8551093_1.merge.gtf
..writing Rubio_transcript_count_matrix.csv
..writing Rubio_gene_count_matrix.csv
All done.
`
    2. deLorgeril started
          `$ sbatch 03_Stringtie_Assembly_Quantify_prepDE_deLorgeril.sh
          $ cat Stringtie_deLorgeril_prepDE_out_2_11_2020
          # shows all samples processed correctly`
    3. Zhang DONE
          `$ sbatch 03_Stringtie_Assemble_Quantify_Zhang_prepDE.sh
          $ cat Stringtie_Zhang_prepDE_out_2_11_2020
          # shows all samples were processed correctly`
    4. He DONE
          `$ sbatch 03_Stringtie_Assemble_Quantify_He_prepDE.sh
          $ cat Stringtie_He_prepDE_out_2_11_2020
          # shows all samples processed correctly`
    5. ROD DONE
          `$ sbatch 03_Stringtie_Assemble_Quantify_ROD_prepDE.sh
          # ROD done
          $ cat Stringtie_ROD_prepDE_out_2_11_2020
          # shows all samples processed correctly`
    6. Probiotic DONE
          `$ sbatch 03_Stringtie_Assemble_Quantify_Probiotic_prepDE.sh
          $ cat Stringtie_Probiotic_prepDE_out_2_11_2020 # all were processed
          # shows all samples were processed correctly
          `
    7. Dermo DONE
          `$ sbatch 03_Stringtie_Assemble_Quantify_Dermo_prepDE.sh
          $ cat Stringtie_Dermo_prepDE_out_2_11_2020`

* While this is running I will set up my new RStudio workspace for analyzing differential expression for each
    - How should I set this up? Using one script and one RStudio workspace may be the best option.
    - Can I partially clear an R studio workspace environment to only save final files that I want to keep?
          ` #  Can remove list of of dataframes that all share a common pattern using the following
            rm(list = ls(pattern = "^tmp"))`
    - Instead of separating directories by species, because I am now analyzing all together, I'm going to put this combined analysis in a new folder
          `/Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/2020_Transcriptome_ANALYSIS`
    - Old scripts and old Dermo 2015 analysis (completed in Fall 2019 from data that had been aligned using CLC) will still be split up by species still.
    - Moving genome files for both species into this new directory
    - Setting up Script based on code format from "Dermo 2015" analysis done in Fall 2019 on Dermo transcriptomes that had been run with CLC.

 * Backing up all current Bluewaves contents on backup hardrive by sshing to data server
        `# Copying data off using following format
        # rsync -rv ./localpath/ remotesystem:/remote/path/
        ssh fs03
        $ cd /data3/marine_diseases_lab/erin/
        # find local path on mac, go to system preferences, sharing, find local IP address under "remote login"
        # To log in to this computer remotely, type “ssh erinroberts@253.59.20.172.s.wireless.uri.edu”.
      $   scp -r ./* erinroberts@253.59.20.172.s.wireless.uri.edu:"/Volumes/EMR\\ Backup/Bluewaves_Backups/2020_Data_Backup/"
        # had to put the path with spaces in it in double quotes!
`

## 2/13/2020 Data Backup finished, download output csvs locally, set up new DESeq2 script

* Created new DESeq2 script which will have all the code combined. Starting to edit.

## 2/14/2020 Editing DESeq2 Script. Creating metadata files for each project

* Starting with Zhang Vibrio to edit DESeq2 code and create template code that will be used and modified for other projects. The analysis for this experiment is relatively straightforward.
  * Creating metadata (coldata) file for experiment using the "Organized_SRA_info.xlsx" spreadsheet to create.
  * ZHANG VIBRIO TRANSCRIPTOME ANALYSIS NOTES
        1. First ran data QC by plotting a PCA of the rlog transformed counts. This PCA indicated that the PBS and control conditions were closely clustered in comparison to the other transcriptomes. The V. tubiashii and the LPS were also somewhat closely clustered. But I decided to keep these separate.
        2. Made DESeq Data set just control vs. treated with the control as PBS and the injected and the treated as LPS, M. lut adn all Vibrio. Then I will use specific contrasts to pull out the differences between controls and each individual challenge type
        3. Re-ran the PCA after calculating the DEseqdataset from matrix, and the PCA didn't change.
        4. Should I remove the MTRG from the beginning?? Testing how this changes my results before keeping in my code.
            * With removing the MSTRG, PBS and the control no longer cluster as closely (though they are still near each other)
            * Now the LPS, PBS and control all cluster
            * `plotPCA(Zhang_dds_rlog, intgroup=c("group", "condition"))` shows that LPS and M. lut cluster most closely. This doesn't at the outset seem biologically relevant because LPS is a control for gram negative bacteria while M. luteus is gram positive 
        5. Now trying to figure out the code to pull out specific contrasts now after the LFC shrinkage with apeglm. Deciding whether or not to subset the list of significant genes for those that had greater than or less than 1 LFC. Doing this drastically decreased the number of genes. Going to compare.  Also loading the apoptosis gene list to find the apoptosis transcript
