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
* C. gigas de Lorgeril et al., 2017 = Download and analyze only the samples from families 11 (susceptible) and 21 (resistant) that were infected with OsHV-1 during a "natural" infection and sequenced. These were all sequenced paired end on hiseq 500. 42 total transcriptomes

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

## 2/14/2020 Analyzing Zhang Vibrio Transcriptome dataset

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
            * Going to proceed with the MSTRG removed
        5. 13 significant apoptosis genes when looking at the overall significant transcription.
        6. Now trying to figure out the code to pull out specific contrasts now after the LFC shrinkage with apeglm. Deciding whether or not to subset the list of significant genes for those that had greater than or less than 1 LFC. Doing this drastically decreased the number of genes. Going to compare. Also loading the apoptosis gene list to find the apoptosis transcript.
            * I could run the lfcsrhink using the "normal" or "ashr" rather than using "apeglm", which would allow me to do the contrasts using the normal procedure? Reading the 2020 vignette extended section on shrinkage estimators says that the "ashr" version is similar to "apeglm" in that it helps preserve the size of large LFcs, but it can also do specific contrasts. Quote from the vignette: "ashr and apeglm can adapt to the scale of the entirety of LFCs, while not over-shrinking the few largest LFCs. The potential for over-shrinking LFC is also why DESeq2’s shrinkage estimator is not recommended for designs with interaction terms.....Finally, normal and ashr can be used with arbitrary specified contrast because normal shrinks multiple coefficients simultaneously (apeglm does not), and because ashr does not estimate a vector of coefficients but models estimated coefficients and their standard errors from upstream methods (here, DESeq2’s MLE). " - USING ASHR DOESN'T MEAN YOU CAN CREATE NEW COEFFICIENTS IF THE VARIABLE WAS NOT INCLUDED IN THE ORIGINAL FORMULA
            * I could also change my initial DESEq2 formula to look at the differences between them all?..Doing more reading about how DESeq1 set up the formulas. Could I possibly include a term called "all", put this in the formula and then pull out specific interactions later? - THIS DOES NOT WORK

* 2/16/2020: comparing the effect of using the original PCA plot to assign the different challenges into similar groups. Based on the oringal rlog transformed PCA plot, I have placed together PBS + control, V. aes and V. alg 2, and LPS M. lut and V. tub, the others are quite far apart on the PCA and I have called them Vibrio. This will allow me specific contrasts for LFC comparison downstream. Assessing how the number of significant DEGs and apoptosis genes differs from having them split into just challenge vs. control. The formula is still `~time + group`, though now group has four different levels.
* 2/17/2020: comparing the results from yesterday to trying a new comparison, where I cluster roughly pathogenic bacteria, non-pathogenic bacteria bacteria (in Zhang_coldata under the path column). Comparing this to splitting up the non-path category into LPS_Mlut and all the non-path vibrios (the group_by_sim column in Zhang_coldata). Also assessing whether adding in a time effect makes sense. There is only one sample at a different time, the original uninjected control. This sample was the only one that was not injected. Controlling for the effect of time would control for the lack of response due to injection.
  * When the formula used is `~time + group_by_sim` vs. `~ group_by_sim` a different list of apoptosis genes comes out as being significantly differentially expressed. Notably, without the effect of time, "diablo homolog" has a very high LFC and is not significantly differentially expressed without this. Conversely, with controlling for the effect of time included, caspase 7 now has a high LFC in all three comparisons with control, and GIMAP4 is highly differentially expressed with non-pathogenic vibrio.  Going to compare the rlog transformed counts of the apoptosis genes between all samples.


* Analysis decisions (for now):
  * The final formula to be used is ~time + group_by_sim. The time effect controls also for the fact that the first control was not injected while the others were and to disregard effects due to this. The exact code is below:
          `Zhang_dds_broken_group <- DESeqDataSetFromMatrix(countData = Zhang_counts,
                                                 colData = Zhang_coldata,
                                                 design = ~time+ group_by_sim)`
  * The groupings used in `group_by_sim` were based off of groupings in the original rlog PCA plot showing clustering of rlog transformed counts. This method of clustering introduces the least bias. The code for this PCA plot is:
          `autoplot(pcZhang,
         data = Zhang_coldata,
         colour="group_by_sim",
         size=5) `
  * The product of DESeq2 differential expression is control+PBS vs Non pathogenic vibrio (V_aes, V_alg1, V_alg2), control+PBS vs pathogenic vibrio (V_ang, V_tub) and LPS and M.lut.
  * Procedural notes: The MSTRG genes are pulled out before the DESeq2 formula and rows with counts less than 10 are removed from analysis. The results function sets alpha to 0.05. For DEG analysis rows with LFC of less than 1 are not filtered out. LFC is calculated using lfcshrink() type "ashr".
  * Plots to be compared between experiments will be the original vst or rlog PCAs, the LFC plots after DESeq2, the plots of rlog transformed count heatmaps where the most variable genes are plotted as the difference from the mean across all samples, and also the PCA plots of the rlog transformed values with just the apoptosis genes included.
        * At how many genes should I cut of this list? Keeping it at the 100 most variable seems to make sense?
        * To allow myself to compare apoptosis expression between all experiments later, I am going to add into the code PCA plotting of the apoptosis rlog transformed data only. Creating an overall `All_coldata.csv` that I can use later to plot a PCA with apoptosis rlog expression data between experiments.
  * The code that includes exploration and multiple testing is now saved as `C_gig_C_vir_transcriptome_DE_GO_Plotting_TESTING.R` while the final code with testing removed is called `C_gig_C_vir_transcriptome_DE_GO_Plotting.R`
* *Talk with Marta about analysis*
  * Marta's thoughts: the variability in the PCA is too high. Grouping the samples other than the control for the sake of increasing the sample size for DESeq2 might not really make sense. Performing an analysis where I compare the rlog transformed counts to the average value for the control may make more sense. Then I could compare this value between experiments. For Marta suggests finishing the other analyses and comparing apoptosis rlog transformed counts between experiments. I could also make heatmaps where I look at the rlog transformed values and subtract the mean across all others from each sample (sort of to normalize)?
  * Marta also thinks I'll need to do some sort of enrichment analysis either enrichment of pathways of interest or gene families of interest
* Ideas from the literature
  * Pillon et al., 2020: chord plot to show number of similarities between groups, network analysis, correlation matrix of fold changes
  * Switonska et al. 2019: network analysis
  * Morrow et al., 2009: correlation matrix of fold changes (looks like heatmap) (Pillon et al., 2020; Morrow et al., 2020)

## 2/18/2020 Starting Analysis of Rubio transcriptomes
* Created `Rubio_coldata.csv`. Had to manually change dashes to underscores because these where in the SRA metadata. R does not like these. Also had to remove `_1` in `Rubio_counts` because this was an artifact of how I named my PE samples. They kep this underscore 1 through to the end of the experiment.
* Added `Rubio_coldata.csv` condition data to the `All_coldata.csv` I'm creating to plot PCAs at the end of my analysis.
* Finished preliminary analysis of Rubio transcriptomes. The formula used to perform differential expression analysis was `~ Condition`. Each treatment type was compared to the `Control_untreated` sample.

## 2/19/2020 Analyzing Probiotic Transcriptomes
* Created `Probiotic_coldata.csv`. As for the Rubio transcriptomes, manually removed the `_1` from all the samples. Added the probiotic experiment Coldata to the `All_coldata.csv` spreadsheet.
* The transcript headers for the `Probiotic_transcript_count_matrix,csv` do not have any XM headers. Was there an issue with the assembly or annotation? Going to check the other C. virginica transcript_count_matrix.csv files before proceeding to make sure I don't have to re-run anything else.
        ` # Checking the Dermo transcriptome count file:
         $ grep 'XM' Dermo_transcript_count_matrix.csv
         # no XMs found
         # Checking ROD XM
         $ grep 'XM'ROD_transcript_count_matrix.csv
         # no XMs found `
  * What happened? Why are there no XM annotations? Going back to my script to check
  * In 2017 when I ran my C_vir pipeline, with Hisat I built my index with `cvir_edited` which had erroneous spaces removed. Then in stringtie I used the normal reference gff file `$F/ref_C_virginica-3.0_top_level.gff3`. This is the same procedure I have followed this year. Error must be caused by another source.
  * Investigating the intermediate merge files, the re-estimated files, and the scripts from Stringtie
        `]$ pwd
/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/2020_Raw_Transcriptome_Data/C_vir_ROD_SRA/ROD_Stringtie_gtf
        $ less SRR1293904.merge.gtf
        # Transcript IDs are all  "rna11" but have LOC Ids
        $ grep 'XM' SRR1293904.merge.gtf
        # no XMs in file

        # No XMs in the HISAT bam output either
        $ less SRR1298711.fastq.gz.clean.trim.filter.gz.bam.gtf
        $ grep 'XM' SRR1298711.fastq.gz.clean.trim.filter.gz.bam.gtf # returns nothing

        # transcript IDs for example in the Cgig files all had XM ids in their merge files
        "rna-XM_011428201.2"
        $ cd /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_deLorgeril_OsHV1_SRA/deLorgeril_Stringtie_gtf
        $ less SRR6679093_1.merge.gtf

        # inspecting HISAT script used to run probiotic Files
        $ nano 02_HISAT2_samtools_sort_Rubio_Pro_deLorgeril.sh
          # index is properly called for..why did it not map
        # inspecting Stringtie Probiotic script
        $ nano 03_Stringtie_Assemble_Quantify_Probiotic_prepDE.sh
          # again no obvious reasons for no mapping of rna to XM

        # Inspecting 2017 R script for analyzing output data and old Stringtie output file from 2017
        $ cd /Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis\ Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA\ ANALYSIS/apoptosis_data_pipeline/DESeq2/C_VIRGINICA_PIPELINE/2017_OLD_Figures_Output/OLD_files
        $ less C_vir_transcript_count_matrix.csv
        $ grep 'XM' C_vir_transcript_count_matrix.csv
        # no XMs in this file either

        # /Users/erinroberts/Documents/PhD_Research/Chapter_1_Apoptosis Paper/Chapter1_Apoptosis_Transcriptome_Analyses_2019/DATA ANALYSIS/apoptosis_data_pipeline/DESeq2/C_VIRGINICA_PIPELINE/SCRIPTS/05_C_Vir_Bac_DESeq2.R

    ##Found at the top of my old DEseq file that I  "Match "rna#" value with the Gene LOC name in the stringtie file". But then how did I get to annotating down to the transcript level?

    ## What I know: there are "rna#" lists in every merge file, the merge.gtf, the merged.annotated.gtf and also in "ref_C_virginica-3.0_top_level.gff3'

    ## Do the `rna#`  match in all of these? Then I can use that match to get the XM values from the `ref_C_virginica-3.0_top_level.gff3` files
    # test example
    $ pwd
    /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/2020_Raw_Transcriptome_Data/C_vir_Probiotic_SRA/Probiotic_Stringtie_gtf
    $ grep 'LOC' Probiotic_stringtie_merged.annotated.gtf
    # test line: NC_035787.1	StringTie	transcript	25840012	25841649	.	+	.	transcript_id "rna52049"; gene_id "MSTRG.23177"; gene_name "LOC111110259"; xloc "XLOC_030060"; ref_gene_id "gene30123"; cmp_ref "rna52049"; class_code "="; tss_id "TSS44959";

    $ pwd
    /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/Cvir_Genome_and_Indexes
    $ grep 'rna52049' ref_C_virginica-3.0_top_level.gff3
    ]$ grep 'rna52049' ref_C_virginica-3.0_top_level.gff3
    NC_035787.1	Gnomon	mRNA	25840012	25841649	.	+	.	ID=rna52049;Parent=gene30123;Dbxref=GeneID:111110259,Genbank:XM_022446679.1;Name=XM_022446679.1;gbkey=mRNA;gene=LOC111110259;model_evidence=Supporting evidence includes similarity to: 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 23 samples with support for all annotated introns;product=uncharacterized LOC111110259;transcript_id=XM_022446679.1
    NC_035787.1	Gnomon	exon	25840012	25840272	.	+	.	ID=id581533;Parent=rna52049;Dbxref=GeneID:111110259,Genbank:XM_022446679.1;gbkey=mRNA;gene=LOC111110259;product=uncharacterized LOC111110259;transcript_id=XM_022446679.1
    NC_035787.1	Gnomon	exon	25840464	25840583	.	+	.	ID=id581534;Parent=rna52049;Dbxref=GeneID:111110259,Genbank:XM_022446679.1;gbkey=mRNA;gene=LOC111110259;product=uncharacterized LOC111110259;transcript_id=XM_022446679.1
    NC_035787.1	Gnomon	exon	25841444	25841649	.	+	.	ID=id581535;Parent=rna52049;Dbxref=GeneID:111110259,Genbank:XM_022446679.1;gbkey=mRNA;gene=LOC111110259;product=uncharacterized LOC111110259;transcript_id=XM_022446679.1
    NC_035787.1	Gnomon	CDS	25840027	25840272	.	+	0	ID=cds47002;Parent=rna52049;Dbxref=GeneID:111110259,Genbank:XP_022302387.1;Name=XP_022302387.1;gbkey=CDS;gene=LOC111110259;product=uncharacterized protein LOC111110259;protein_id=XP_022302387.1
    NC_035787.1	Gnomon	CDS	25840464	25840583	.	+	0	ID=cds47002;Parent=rna52049;Dbxref=GeneID:111110259,Genbank:XP_022302387.1;Name=XP_022302387.1;gbkey=CDS;gene=LOC111110259;product=uncharacterized protein LOC111110259;protein_id=XP_022302387.1
    NC_035787.1	Gnomon	CDS	25841444	25841497	.	+	0	ID=cds47002;Parent=rna52049;Dbxref=GeneID:111110259,Genbank:XP_022302387.1;Name=XP_022302387.1;gbkey=CDS;gene=LOC111110259;product=uncharacterized protein LOC111110259;protein_id=XP_022302387.1

    # YES the LOC ids match here. The merged file transcript ID is in fact just the "Parent=" in the original ref_C_virginica-3.0_top_level.gff3. Why did the software not map down to the XM level? Is there something wrong with how it was interpretting the gff3 file?

    # Comparing C gig and C vir reference gff3 files
    $ less GCF_000297895.1_oyster_v9_genomic.gff
    NW_011934501.1  Gnomon  CDS     3692    3785    .       -       1       ID=cds-XP_011445785.1;Parent=rna-XM_011447483.2;Dbxref=GeneID:105338483,Genbank:XP_011445785.1;Name=XP_011445785.1;gbkey=CDS;gene=LOC105338483;product=DNA-directed RNA polymerase III subunit RPC3;protein_id=XP_011445785.1

    # Conclusion: the mapper is mapping both to the Parent column, it is just that in the C.gig reference genome, the Parent field includes the XM information

    # Does every parent rna only have one XM?
    [erin_roberts@bluewaves Cvir_Genome_and_Indexes]$ grep 'rna1;' ref_C_virginica-3.0_top_level.gff3
    NC_035780.1	Gnomon	mRNA	28961	33324	.	+	.	ID=rna1;Parent=gene1;Dbxref=GeneID:111126949,Genbank:XM_022471938.1;Name=XM_022471938.1;gbkey=mRNA;gene=LOC111126949;model_evidence=Supporting evidence includes similarity to: 3 Proteins%2C and 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 21 samples with support for all annotated introns;product=UNC5C-like protein;transcript_id=XM_022471938.1
    NC_035780.1	Gnomon	exon	28961	29073	.	+	.	ID=id4;Parent=rna1;Dbxref=GeneID:111126949,Genbank:XM_022471938.1;gbkey=mRNA;gene=LOC111126949;product=UNC5C-like protein;transcript_id=XM_022471938.1
    NC_035780.1	Gnomon	exon	30524	31557	.	+	.	ID=id5;Parent=rna1;Dbxref=GeneID:111126949,Genbank:XM_022471938.1;gbkey=mRNA;gene=LOC111126949;product=UNC5C-like protein;transcript_id=XM_022471938.1
    NC_035780.1	Gnomon	exon	31736	31887	.	+	.	ID=id6;Parent=rna1;Dbxref=GeneID:111126949,Genbank:XM_022471938.1;gbkey=mRNA;gene=LOC111126949;product=UNC5C-like protein;transcript_id=XM_022471938.1
    NC_035780.1	Gnomon	exon	31977	32565	.	+	.	ID=id7;Parent=rna1;Dbxref=GeneID:111126949,Genbank:XM_022471938.1;gbkey=mRNA;gene=LOC111126949;product=UNC5C-like protein;transcript_id=XM_022471938.1
    NC_035780.1	Gnomon	exon	32959	33324	.	+	.	ID=id8;Parent=rna1;Dbxref=GeneID:111126949,Genbank:XM_022471938.1;gbkey=mRNA;gene=LOC111126949;product=UNC5C-like protein;transcript_id=XM_022471938.1
    NC_035780.1	Gnomon	CDS	30535	31557	.	+	0	ID=cds0;Parent=rna1;Dbxref=GeneID:111126949,Genbank:XP_022327646.1;Name=XP_022327646.1;gbkey=CDS;gene=LOC111126949;product=UNC5C-like protein;protein_id=XP_022327646.1
    NC_035780.1	Gnomon	CDS	31736	31887	.	+	0	ID=cds0;Parent=rna1;Dbxref=GeneID:111126949,Genbank:XP_022327646.1;Name=XP_022327646.1;gbkey=CDS;gene=LOC111126949;product=UNC5C-like protein;protein_id=XP_022327646.1
    NC_035780.1	Gnomon	CDS	31977	32565	.	+	1	ID=cds0;Parent=rna1;Dbxref=GeneID:111126949,Genbank:XP_022327646.1;Name=XP_022327646.1;gbkey=CDS;gene=LOC111126949;product=UNC5C-like protein;protein_id=XP_022327646.1
    NC_035780.1	Gnomon	CDS	32959	33204	.	+	0	ID=cds0;Parent=rna1;Dbxref=GeneID:111126949,Genbank:XP_022327646.1;Name=XP_022327646.1;gbkey=CDS;gene=LOC111126949;product=UNC5C-like protein;protein_id=XP_022327646.1

    # The RNAID is in the `Parent=` field of exon and CDS lines, but it is in the `ID` section of mRNA lines! I can use the ID line to match to the XMs

    # Moving forward: I will match the rnaID from the  C. virginica to the rna# in the ID column of the reference genome to get the XM ids and product name. The code I created in 2017 for working with the stringtie merged file is not necessary.
`
* While I'm logged in, completing my data backup that was interrupted last week.
        `# Need to only redo the C_gig_deLorgeril_OsHV1 to finish all the C_gig data backup. Then will move on to finish C. virginica
        $ ssh fs03
        $ cd /data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data/C_gig_deLorgeril_OsHV1_SRA/
        $ scp -r ./C_gig_deLorgeril_OsHV1_SRA/ erinroberts@253.59.20.172.s.wireless.uri.edu:"/Volumes/EMR\\ Backup/Bluewaves_Backups/2020_Data_Backup/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/2020_Raw_Transcriptome_Data"
        # finished C_gig data back up

        # C_vir_Probiotic_SRA didn't fully load
        $ scp -r ./C_vir_Probiotic_SRA/ erinroberts@253.59.20.172.s.wireless.uri.edu:"/Volumes/EMR\\ Backup/Bluewaves_Backups/2020_Data_Backup/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/2020_Raw_Transcriptome_Data"

        # Backup complete
        `

* The data set includes samples from two treatments and three timepoints. There are three treatment replicates, but no replicates within treatment for each time point. PCA plots on rlog transformed counts show that a large portion of the variation is determined by the effect of time rather than the effect of treatment. To control for this effect, DESeq2 formula will be `~ Time + Condition` to control for the effect of time.

## 2/19/2020 - 2/20/2020 Analyzing ROD transcriptomes, HE transcriptomes, and DeLorgeril transcriptomes
*ROD Transcriptomes*
  * A note about the data: the susceptible family lacks controls because they died during the experiment unexpectedly. Mcdowell et al., 2014 used PCA plots to group together the "early" 1d and 5d responses and the "late" 15d and 30d responses. My rlog transformed PCA plots support this grouping. Therefore I split the susceptible and resistant family into completely separate analysis. The ROD Susceptible F3L family was compared between Early and Late using the formula `~ Condition`, while the control resistant and ROD challenged Resistant family were compared with the formula `~ Time + Condition` to control for the effect of time on the outcome because PCAs showed a large time component.
* The Resistant family has significantly fewer significant genes overall as compared to the susceptible. The susceptible show a large response while the resistant overall show a lack of response.

*HE transcriptomes*
  * Created `He_coldata.csv` which separates into control and OsHV-1 conditions, with a time column that has the times for sampling. I included the time 0 into the control condition. The formula used for DEseq2 is `~Time + Condition` to control for the effect of time.

*deLorgeril Transcriptomes*
  * Created `deLorgeril_coldata.csv`. Considering whether I should correcting for the effect of time or focusing on the acute response. I could go back and adjust the others to be also looking for the acute response. I have separated the resistant and susceptible families and run DEseq separately.
  * Analyzing the 48hr vs control, 60hr vs control, and 72hr vs control in order to assess the acute response
  * 60hr appears to be the acute response time point for each (most significantly differentially expressed genes), but more so in the susceptible family

## 2/21/2020 Dermo Transcriptome analysis
* Creating `Dermo_coldata.csv`. Flenames and technical replicates with their sample names using the `SRA_metadata_2019_paper_11_15_2019.xlsx` file and then added sample metadata (timepoint, logConc) by matching samplenames to those in `Dermo_trancript_count_matrix.csv` file. Not all the timepoints had units for time in the metadata, added units. Added an h to represent hours for each timepoint. Added Family column. Adding in library prep date also because this was found to be the cause of batch effects in Proestou et al., 2020.
* Using the formula `~Lib_prep_date + Condition + Time` to control for the batch effect and then compare treatment and control and then look at the timepoint with the acute response. The Susceptible family has a greater response

## 2/24/2020 Compare transcriptomes across groups
* 1. Combined apoptosis transcript matrices subset for apoptosis transcripts into a single table for each species. Then removed batch effects using limma to correct for the batch effect of each experiment. Then plotted PCAs and heatmaps for each species. This showed that ROD and probiotic were similar in response for Cvir and for Cgig the OsHV1 challenges were similar.
  * I have been subsetting my genes first and then performing count transformation with rlog or vst after subsetting. Is this appropriate?? Should I perform rlog and then subset for apoptosis genes afterwards?
    * https://support.bioconductor.org/p/109809/: Can plot PCA for a subset of genes using `plotPCA(vsd[idx,])` where IDx is the rows to subset
    * https://www.biostars.org/p/336298/
  * I think I should do the apoptosis subsetting the PCA after the vsd with the full dataset and correction for batch effect. I have added a new set of code that subsets for apoptosis genes after combining and performing vst.  
  * I merge transcript tables by starting with the largest and add in 0 when it is missing. However, I should probably be only working with the consensus set of transcripts.....that way its not like I'm comparing apples to oranges in my heatmaps.
  * I can investigate presence and absence of specific genes and transcripts by comparing the merged.annotated.gtf files from each transcriptome
  *Only work with the consensus set of transcripts when comparing experiments*
    - C virginica experiments have relatively close numbers of transcripts (each row is a transcript)
     `nrow(Dermo_counts) # 67868
      nrow(Probiotic_counts) # 67876
      nrow(ROD_counts) # 67870`
    - C gigas He experiment has a lot more transcripts however. The HE experiment may be affecting results
     `nrow(deLorgeril_counts) # 53701
      nrow(He_counts) # 86859
      nrow(Zhang_counts ) # 53705
      nrow(Rubio_counts) # 53705`
  * I've performed my apoptosis analysis for heatmaps in three different ways, and the results keep coming up the same. Merging vs not merging the gene lists first has little effect.


* 2. In order to combine across species, I think it makes more sense to look at apoptosis gene counts, so I make fewer assumptions about which transcripts are which. Going to try to join by gene name. What do I do for genes with the same name that have multiple copies? I could just choose the one with the highest average counts in each species or the most variance? Could try it both ways?
  * As above, work with the consensus set of genes based on shared gene names (this is the best I can do short of aligning all the genes together and taking their best hits and creating a table). For cases where there are multiple entries per gene
  * Could I do a whole genome alignment to actually match together the most similar genes?

  * Pausing on this while I do some additional reading below

* 3. Additional Reading regarding methods to compare gene expression across species (see notes in "Meta-analysis_of_transcritomes.docx")
  * Martin and Frasier (2018) Nature Communications. Comparative expression profiling reveals widespread coordinated evolution of gene expression across eukaryotes
  * Z hou, Yan, et al. “A Statistical Normalization Method and Differential Expression Analysis for RNA-Seq Data between Different Species.” BMC Bioinformatics, vol. 20, no. 1, BMC Bioinformatics, 2019, pp. 1–10, doi:10.1186/s12859-019-2745-1.
  * Brawand, David, et al. “The Evolution of Gene Expression Levels in Mammalian Organs.” Nature, vol. 478, no. 7369, 2011, pp. 343–48, doi:10.1038/nature10532.
  * Davidson, Rebecca M., et al. “Comparative Transcriptomics of Three Poaceae Species Reveals Patterns of Gene Expression Evolution.” Plant Journal, vol. 71, no. 3, 2012, pp. 492–502, doi:10.1111/j.1365-313X.2012.05005.x.
  * Söllner, Julia F., et al. “An RNA-Seq Atlas of Gene Expression in Mouse and Rat Normal Tissues.” Scientific Data, vol. 4, The Author(s), 2017, pp. 1–11, doi:10.1038/sdata.2017.185.

  * CONCLUSION: I NEED TO IDENTIFY ORTHOLOGS BETWEEN EASTERN AND PACIFIC OYSTER IN ORDER TO MAKE COMPARISONS BETWEEN EXPRESSION. I need to use the normalised counts of 1:1 orthologs. Then I can pipe these into WGCNA for analysis (See Yu et al 2020 for a good example of this). A common tool for this is OrthoFinder.
    * `OrthoFinder/2.3.3-foss-2018b-Python-2.7.15` is available on the cluster!!

* 3. Pathway enrichment analysis? Venn diagrams of the number of transcripts or genes expressed per challenge across gene families

## 2/25/2020 Researching methods for finding orthologs between two species.

1. General methods used in papers
  1. Using the Ensembl database of 1:1 orthologous genes for each pair of species. Also can align cDNA sequences of orthologous gene families using TBA (see Cardoso-Moreira 2019 https://www.nature.com/articles/nature10532#Sec8. They filtered alignments for those that had no gaps) (Sollner et al. 2017, Ferrari 2019)
  2. Using OrthoFinder to identify orthologs (Yu 2020)
  3. Using the UniPropt100 database (Martin and Frasier 2018)

2. OrthoFinder seems to be a widely used, applicable, and respected method for identifying orthologs between species (which is my goal). I am going through the OrthoFinder tutorial now. https://davidemms.github.io

  - Input for Orthofinder is all of protein coding sequences for the species. The input files need to be unzipped in order to be used
  - Downloaded the protein sequences for Cvir `GCF_002022765.2_C_virginica-3.0_protein.faa`
    $ scp GCF_002022765.2_C_virginica-3.0_protein.faa erin_roberts@bluewaves:/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_Vir_subset/Cvir_Genome_and_Indexes/
  - Downloaded all the protein sequnces for Cgig
    $ scp GCF_000297895.1_oyster_v9_protein.faa erin_roberts@bluewaves:/data3/marine_diseases_lab/erin/2017_2020_Transcriptome_Analysis/pipeline_files/C_gig_Bac_Viral_subset/Cgig_Genome_and_Indexes/
