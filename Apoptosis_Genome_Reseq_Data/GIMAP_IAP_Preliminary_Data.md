# Preliminary Gene Family Analysis Data for GIMAP and IAP Gene Families

## Goals:
1. Gather preliminary data on patterns of gene family variation across natural and selected populations
of the eastern oyster for diversified gene families IAP and GIMAP, who both are functionally relevant members
of the apoptosis pathway.

2. Use sequence based approaches to confirm presence of apoptosis genes in oyster

# Data set and Additional Files Necessary to perform the analyses
1. HiSeq Data Sets sequenced on an Illumina HiSeq X Ten PE150, sequenced without PCR concentration. The full sequence
set includes data from 14 different populations from both wild and selected lines of eastern oysters.

| Population Name       | Region         | Location Sampled                           | Ecotype       | Wild/Selected |
|-----------------------|----------------|--------------------------------------------|---------------|---------------|
|1. Laguna Madre (LM)   | Gulf of Mexico | Port Mansfield, TX                         | NA            | Wild          |
|2. Hope Creek (HC)     | Delaware Bay   | Hope Creek, NJ                             | Low salinity  | Wild          |
|3. Cape Shore (CS)     | Delaware Bay   | Cape Shore, NJ                             | High salinity | Wild          |
|4. Sherman Marsh (SM)  | Maine          | Sherman Marsh/Sheepscot River, ME          | Low salinity  | Wild          |    
|5. Hog Island (HI)     | Maine          | Hog Island, Damariscotta River Estuary, ME | High salinity | Wild          |
|6. Sister Lake (SL)    | Louisiana      | Caillou Lake, LA                           | NA            | Wild          |
|7. Calcasieu Lake (CL) | Louisiana      | Grand Isle, LA                             | NA            | Wild          |
|8. Chlora's Point (CLP)| Chesapeake Bay | Choptank River- Chesapeake Bay, VA         | Low salinity  | Wild          |
|9. Hummock Cove (HC-VA)| Chesapeake Bay | Chesapeake Bay, VA                         | High salinity | Wild          |


# Step 1: Download and acquire the Data

All sequence data has been downloaded onto KITT server by JP. To access the data a symbolic link was created using the following commands.

```
# Link to only the natural population files
cd /home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/natural_pop_files

declare -a forward_array=(/home/Shared_Data/LM*.F.fq.gz /home/Shared_Data/HC_?.F.fq.gz /home/Shared_Data/CS*.F.fq.gz /home/Shared_Data/SM*.F.fq.gz /home/Shared_Data/HI*.F.fq.gz /home/Shared_Data/SL*.F.fq.gz /home/Shared_Data/CL_*.F.fq.gz /home/Shared_Data/CLP*.F.fq.gz /home/Shared_Data/HC_VA*.F.fq.gz)

declare -a reverse_array=(/home/Shared_Data/LM*.R.fq.gz /home/Shared_Data/HC_?.R.fq.gz /home/Shared_Data/CS*.R.fq.gz /home/Shared_Data/SM*.R.fq.gz /home/Shared_Data/HI*.R.fq.gz /home/Shared_Data/SL*.R.fq.gz /home/Shared_Data/CL_*.R.fq.gz /home/Shared_Data/CLP*.R.fq.gz /home/Shared_Data/HC_VA*.R.fq.gz)

for i in ${forward_array[@]}; do
  ln -s ${i} /home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/natural_pop_files
done

for i in ${reverse_array[@]}; do
  ln -s ${i} /home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/natural_pop_files
done

#check if I have the correct number of sequence Sets
ls *.fq.gz | wc -l # 106

#Link to only the selected population files
mkdir selected_pop
cd selected_pop
$ ln -s /home/Shared_Data/LOLA* .
$ ln -s /home/Shared_Data/DEBY* .
$ ln -s /home/Shared_Data/OBOYS2* .
$ ln -s /home/Shared_Data/NEH* .
$ ln -s /home/Shared_Data/NG* .
```

# Step 2: Check Sum of Genome Quebec Data
During downloading process JP used wget command supplied by GenomeQuebec to automatically download all data, generate checksum file, and perform checksum.

```
wget -O - "https://genomequebec.mcgill.ca/nanuqMPS/readsetList?projectId=15565&tech=HiSeq" --no-cookies --no-check-certificate --post-data 'j_username=username&j_password=password' | wget --no-cookies --no-check-certificate --post-data 'j_username=username&j_password=password' -ci -
```

#  Step 3: Initial Raw Data Assessment
1. Read counts

Read counts are provided in metadata file and checked using checksum.

2. Read quality of each sample group

Use FASTQC to create FASTQC report with:
  a. Basic statistics
  b. Per base sequence quality scores
  c. Per sequence quality scores
  d. Per base GC content
  e. Per sequence GC content
  f. Per sequence GC content

The following commands were all repeated for the selected_pop files from DEBY, LOLA, HG, NEH, OBOYS2, NG
  ```
# Upload FASTQC environment using conda create
conda create -n natural_pop_files fastqc

# Activate the environment
source activate natural_pop_files

# Generate sequence quality report for each population using command fastqc
declare -a pop_array=(./LM*.F.fq.gz ./HC_?.F.fq.gz ./CS*.F.fq.gz ./SM*.F.fq.gz ./HI*.F.fq.gz ./SL*.F.fq.gz ./CL_*.F.fq.gz ./CLP*.F.fq.gz ./HC_VA*.F.fq.gz)

for i in ${pop_array[@]}; do
  fastqc ${i}
  echo "finished ${i}" $date
done

declare -a pop_R_array=(./LM*.R.fq.gz ./HC_?.R.fq.gz ./CS*.R.fq.gz ./SM*.R.fq.gz ./HI*.R.fq.gz ./SL*.R.fq.gz ./CL_*.R.fq.gz ./CLP*.R.fq.gz ./HC_VA*.R.fq.gz)

for i in ${pop_R_array[@]}; do
  fastqc ${i}
  echo "finished ${i}" $date
done

# Faster approach if all files are in the same folder and can disown the process source

for i in *.fq.qz; do
  fastqc $i
done

# Put process into the background using
^Z
bg
disown -a

# Check that all files have been completed by comparing arrays with bash
declare -a total_array=(./LM*.F.fq.gz ./HC_?.F.fq.gz ./CS*.F.fq.gz ./SM*.F.fq.gz ./HI*.F.fq.gz ./SL*.F.fq.gz ./CL_*.F.fq.gz ./CLP*.F.fq.gz ./HC_VA*.F.fq.gz ./LM*.R.fq.gz ./HC_?.R.fq.gz ./CS*.R.fq.gz ./SM*.R.fq.gz ./HI*.R.fq.gz ./SL*.R.fq.gz ./CL_*.R.fq.gz ./CLP*.R.fq.gz ./HC_VA*.R.fq.gz)
for i in ${total_array[@]}; do
  echo ${i} | sed 's/.fq.gz/_fastqc.html/g' >> total_array.txt
done
declare -a total_array=($(cat total_array.txt))
declare -a html_array=(./*.html)
echo ${total_array[@]} ${html_array[@]} | tr ' ' '\n' | sort | uniq -u #any left over are those that weren't converted to fastqc
# Failed to process file HI_2.R.fq.gz. This file is excluded from FASTQC analysis

#Code used for selected population files
#!/bin/bash
#conda create -n selected_pop fastqc
source activate selected_pop

for i in *.fq.gz; do
	fastqc $i
	echo "done ${i}"
done

# use MultiQC to put together all Files
#make new folder with all Fastqc files
mkdir fastqc_results
mv *.zip ./fastqc_results
mv /home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/natural_pop_files/fastqc_results/*.zip ./fastqc_results
cd fastqc_results
conda install -c bioconda multiqc
multiqc .

# Export file to local folder so you can open it up with .html
cd /Users/erinroberts/Documents/PHD_coursework_TA/Puritz_pop_gen
scp -P 2292 eroberts@kitt.uri.edu:/home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/natural_pop_files/fastqc_results/multiqc_report.html .
```

### MultiQC Results

1. GC content: The GC content was either 35% or 36% across all reads.
2. % failed: None failed the QC report
3. % Duplicates: The highest percent duplicates was 12.5%. Most were lower than this.
4. Mean Quality Score: 105/105 Samples Passed

## MultiQC for all Populations
This graph depicts the mean quality score at each position across a read.
![Mean Quality Score](https://github.com/jpuritz/BIO_594_2018/blob/master/FinalAssignment/EMR_Final_Assignment/Mean_quality_histogram_natural_pops_5_3_18.png "This graph depicts the mean quality score at each position across a read.")

5. Per Sequence Quality Score: 105/105 Samples Passed

This graph depicts the number of reads with average quality scores.
![Per Sequence Quality Score](https://github.com/jpuritz/BIO_594_2018/blob/master/FinalAssignment/EMR_Final_Assignment/Per_sequence_quality_scores_natural_pops_5_3_18.png "This graph depicts the number of reads with average quality scores.")

6. Per Sequence GC Content

This graph depicts the average GC contents of reads and is roughly normally distributed.
![Per Sequence GC Content](https://github.com/jpuritz/BIO_594_2018/blob/master/FinalAssignment/EMR_Final_Assignment/Per_sequence_GC_content_natural_pops_5_3_18.png "This graph depicts the average GC contents of reads and is roughly normally distributed")


Overall the sequence quality is high and we will proceed with further analysis.


#### Step 4: Bioinformatic Processing

1. Adaptor trimming using fastqc-mcf
The user manual for this tool can be found here: https://github.com/ExpressionAnalysis/ea-utils/blob/wiki/FastqMcf.md

The adapters used for this project for multiplexing on a Hiseq 2000 were as follows:

>Hi_seq_Adaptor_Read_1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>Hi_seq_Adaptor_Read_2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

These sequences were input into a file called Hi_seq_adaptors.fa.


Install the software fastqc-mcf (from the ea-utils package)
```
#install ea-utils
conda install -c bioconda ea-utils
```
Run the following script to perform read trimming of HiSeq 2000 adapter, trim low quality ends of reads with a Phred score of less than 20 and remove whole reads with an average score of less than 10.

```
sh -c 'for file in "CL_1" "CL_2" "CL_3" "CL_4" "CL_5" "CL_6" "CLP_1" "CLP_2" "CLP_3" "CLP_4" "CLP_5" "CLP_6" "CS_1" "CS_2" "CS_3" "CS_5" "CS_6" "CS_7" "HC_1" "HC_3" "HC_4" "HC_5" "HC_6" "HC_7" "HC_VA_1" "HC_VA_2" "HC_VA_3" "HC_VA_4" "HC_VA_5" "HC_VA_6" "HI_1" "HI_2" "HI_3" "HI_4" "HI_5" "HI_6" "LM_1_pool" "LM_3" "LM_4" "LM_7" "LM_8" "SL_1" "SL_2" "SL_3" "SL_4" "SL_5" "SL_6" "SM_10" "SM_11" "SM_12" "SM_7" "SM_8" "SM_9"
do
echo "start ${file} $(date)"
/home/eroberts/miniconda3/bin/fastq-mcf \
/home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/natural_pop_files/Hi_seq_adaptors.fa \
/home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/natural_pop_files/${file}.F.fq.gz \
/home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/natural_pop_files/${file}.R.fq.gz \
-l 100 \
-q 20 \
-w 5 \
-x 10 \
-u \
-P 33 \
-o /home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/natural_pop_files/${file}.F.cleaned.fq.gz \
-o /home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/natural_pop_files/${file}.R.cleaned.fq.gz &> /home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/natural_pop_files/${file}.log
echo "${file} done $(date)"
done'

# o =output
# l = minumum remaining sequence Length
# q = quality threshold causing base removal at the end of reads
# w = window size for for quality trimming
# -x = 'N' bad read percentage causing cycle removal

```

# Step 5. Read Mapping to Reference using BWA-MEM
BWA-MEM is recommended for longer sequences that range from 70bp to 1Mbp. It is recommended over BWA-SW and BWA-backtrak because it is faster and more accurate.

The reference file being used is the eastern oyster reference mRNA from NCBI: GCF_002022765.2_C_virginica-3.0_rna.fna. Using this file means that any aligned genes will have the gene information attached.

1. Set up a reference index using samtools faidx

```
# copy genome into folder from bluewaves
scp erin_roberts@bluewaves:/data3/marine_diseases_lab/shared/GCA_002022765.4_C_virginica-3.0_genomic.fna .
# create a reference index in samtools
samtools faidx GCA_002022765.4_C_virginica-3.0_genomic.fna

```
Create a bwa index for use by BWA and then use bwa mem to generate sam records for each read and samtools to sort the bam file.

Put the following code into a bash script and run the bash script to execute.
$ nano bwa.sh
$ bash bwa.sh
```
#!/bin/bash
F=/home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/natural_pop_files
array1=($(ls *.F.cleaned.fq.gz | sed 's/.F.cleaned.fq.gz//g'))

#bwa index GCA_002022765.4_C_virginica-3.0_genomic.fna
echo "done index $(date)"

for i in ${array1[@]}; do
  bwa mem $F/GCA_002022765.4_C_virginica-3.0_genomic.fna ${i}.F.cleaned.fq.gz ${i}.R.cleaned.fq.gz -t 8 -a -M -B 3 -O 5 -R "@RG\tID:${i}\tSM:${i}\tPL:Illumina" 2> bwa.${i}.log | samtools view -@4 -q 1 -SbT $F/GCA_002022765.4_C_virginica-3.0_genomic.fna - > ${i}.bam
  echo "done ${i}"
done

array2=($(ls *.bam | sed 's/.bam//g'))

#now sort the bam files with samtools sort
for i in ${array2[@]}; do
  samtools sort -@8 ${i}.bam -o ${i}.bam && samtools index ${i}.bam
done
```

# STEP 6. Mark and filter out any potential duplicate reads using PICARD
Because libraries were generated without a PCR prep, there should not be duplicate reads. However, this serves as a check of this.

```
# Download the PICARD tool
wget https://github.com/broadinstitute/picard/releases/download/2.17.8/picard.jar

# Mark and output duplicates
#!/bin/bash
F=/home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/natural_pop_files
array1=($(ls *.bam| sed 's/.bam//g'))

for i in ${array1[@]}; do
  java -Xms4g -jar picard.jar MarkDuplicatesWithMateCigar I=${i}.bam O=${i}.md.bam M=${i}_dup_metrics.txt MINIMUM_DISTANCE=302
done

echo "done $(date)"

#initially received error: Exception in thread "main" picard.PicardException: Found a samRecordWithOrdinal with sufficiently large clipping that we may have
 missed including it in an early duplicate marking iteration.  Please increase the minimum distance to at least 302bp
to ensure it is considered (was 300), so the minimum distance was increased.
```
Now we are able to remove duplicates as well as any secondary alignments, mappings with a quality score of less than ten, and reads with more than 80 bp clipped. Finally we can create a BAM index from our fully processed files.

```
#!/bin/bash
F=/home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/natural_pop_files
array1=($(ls *.md.bam | sed 's/.md.bam//g'))

for i in ${array1[@]}; do
  samtools view -@8 -h -F 0x100 -q 10 -F 0x400 ${i}.md.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/'| samtools view -@8 -b > ${i}.F.bam
  echo "done removing dups ${i} $(date)"
done

array2=($(ls *.F.bam | sed 's/.F.bam//g'))
for i in ${array2[@]}; do
        samtools index ${i}.F.bam
        echo "done indexing ${i}"
done

# received error however for CL_1 and HC_5 that says [W::bam_hdr_read] EOF marker is absent. The input is probably truncated. I'm going to see if this works.
```

The command below can then be used to check that reads have been filtered out.

```
F=/home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/natural_pop_files
array1=($(ls *.md.bam | sed 's/.md.bam//g'))

for i in ${array1[@]}; do
  paste <(samtools view -c ${i}.md.bam) <(samtools view -c ${i}.F.bam )
done

# receieved the following error that
[W::bam_hdr_read] EOF marker is absent. The input is probably truncated
[E::bgzf_read] Read block operation failed with error -1 after 219 of 291 bytes
[main_samview] truncated file.

```

Below is a sample output of checking whether duplicates have been removed: This indicates that it was successful.

```
76346792	63466956
45962279	38577443
40905430	34076692
```

# STEP 7. Calculate depth per bp along the reference using samtools

The final filtered bam files have now been generated and we can check the sequencing depth per base pair. The output is a text file with three columns. The first column lists the chromosome, the second column is the base pair and the third column is the depth at that base pair.

```
F=/home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/natural_pop_files
array1=(ls *.F.bam | sed 's/.F.bam//g')
for i in ${array1[@]; do
  samtools depth -aa ${i}.F.bam > ${i}.genome.depth
done
```

# Works Cited

MultiQC: Summarize analysis results for multiple tools and samples in a single report
Philip Ewels, Måns Magnusson, Sverker Lundin and Max Käller
Bioinformatics (2016)
doi: 10.1093/bioinformatics/btw354
PMID: 27312411
