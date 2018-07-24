# Preliminary Gene Family Analysis Data for Gene Families

## Goals:
1. Gather preliminary data on patterns of gene family variation across natural and selected populations
of the eastern oyster for the exapnded genefamily GIMAP, a functionally relevant member
of the apoptosis pathway.

# Data set
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
|10. OBOYS2             | Louisiana      | Calcasieu Lake, LA                         | NA            | Selected      |
|11. DEBY               | Virginia       | York River, VA                             | Selected-High | Selected      |
|12. LOLA               | Virginia       | Unknown                                    | Selected-Low  | Selected      |
|13. NEH93 (NEH)        | Delaware Bay   | Unknown                                    | Selected      | Selected      |
|14. UMFS               | Maine          | Unknown                                    | Selected      | Selected      |

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
3. % Duplicates: The highest percent duplicates was 14.8%. Most were lower than this.
4. Mean Quality Score: 167/167 Samples Passed

## MultiQC for data from all Populations
This graph depicts the mean quality score at each position across a read.
![Mean Quality Score](https://github.com/erinroberts/Apoptosis-Genome-Re-sequencing-Data/blob/master/CV_Gen_Reseq_Sequence_Quality_Histogram.png "This graph depicts the mean quality score at each position across a read.")

5. Per Sequence Quality Score: All Samples Passed

This graph depicts the number of reads with average quality scores.
![Per Sequence Quality Score](https://github.com/erinroberts/Apoptosis-Genome-Re-sequencing-Data/blob/master/CV_Gen_Reseq_Per_Sequence_Quality_Score.png "This graph depicts the number of reads with average quality scores.")

6. Per Sequence GC Content

This graph depicts the average GC contents of reads and is roughly normally distributed.
![Per Sequence GC Content](https://github.com/erinroberts/Apoptosis-Genome-Re-sequencing-Data/blob/master/CV_Gen_Reseq_Per_Sequence_GC_Content.png "This graph depicts the average GC contents of reads and is roughly normally distributed")


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
Get list of unique selected population file prefixes with the following code.
```
cat list | sed 's/\..*$//' | uniq | sed 's/\(.*\)/\"\1\"/'
```

Run the following script to perform read processing as above
```
#!/bin/bash
conda install -c bioconda ea-utils

sh -c 'for file in "DEBY_1" "DEBY_2" "DEBY_3" "DEBY_4" "DEBY_5" "DEBY_6" "HG_HG0F2" "HG_HG2F1" "HG_HG2M5" "LOLA_1" "LOLA_2" "LOLA_3" "LOLA_4" "LOLA_5" "LOLA_6" "NEH_1" "NEH_2" "NEH_3" "NEH_4" "NEH_5" "NEH_6" "NG_NH0H4" "NG_NH2F6" "NG_NH2F8" "NG_NH2M1" "OBOYS2_1" "OBOYS2_2" "OBOYS2_3" "OBOYS2_4" "OBOYS2_5" "OBOYS2_6"
do
echo "start ${file} $(date)"
/home/eroberts/miniconda3/bin/fastq-mcf \
/home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/natural_pop_files/Hi_seq_adaptors.fa \
/home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/selected_pop/${file}.F.fq.gz \
/home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/selected_pop/${file}.R.fq.gz \
-l 100 \
-q 20 \
-w 5 \
-x 10 \
-u \
-P 33 \
-o /home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/selected_pop/${file}.F.cleaned.fq.gz \
-o /home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/selected_pop/${file}.R.cleaned.fq.gz &> /home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/selected_pop/${file}.log
echo "${file} done $(date)"
done'
```

# Step 5. Read Mapping to Reference using BWA-MEM
BWA-MEM is recommended for longer sequences that range from 70bp to 1Mbp.

From here down the codes are only presented for the natural population files. The same analysis was repeated for the selected populations.

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


```

Below is a sample output of checking whether duplicates have been removed: This indicates that it was successful.

```
76346792	63466956
45962279	38577443
40905430	34076692
```

# STEP 7. Mapping statistics with Samtools
1. Calculate depth per bp along the reference using samtools

The final filtered bam files have now been generated and we can check the sequencing depth per base pair. The output is a text file with three columns. The first column lists the chromosome, the second column is the base pair and the third column is the depth at that base pair.

```
F=/home/eroberts/repos/BIO_594_2018/FinalAssignment/EMR_Final_Assignment/natural_pop_files
array1=(ls *.F.bam | sed 's/.F.bam//g')
for i in ${array1[@]}; do
  samtools depth -aa ${i}.F.bam > ${i}.genome.depth
done
```

2.  Simple mapping statistics using SAMtools flagstat
```
for i in *.F.bam; do
  samtools flagstat $i
  echo "done $i flagstat"
done

```

# STEP 8: Identify GIMAP Gene Family Regions of Interest

Preliminary Analysis of annotated GIMAP genes in the Reference Genome (.gff) reveal
53 total GIMAP genes.

Table 1: Summarized GIMAP Gene Data for the Reference Annotation

| Total GIMAP Genes | 53 |
|-------------------|----|
| GIMAP4            | 42 |
| GIMAP7            | 9  |
| GIMAP8            | 2  |

Table 2: GIMAP Genes Per Chromosome in Eastern Oyster Reference Annotation

| Chromosome | Number of Total Genes | Number of GIMAP 4 | Number of GIMAP 7 | Number of GIMAP 8 |
|------------|-----------------------|-------------------|-------------------|-------------------|
| CHR2       | 3                     | 1                 | 0                 | 2                 |
| CHR4       | 5                     | 5                 | 0                 | 0                 |
| CHR5       | 1                     | 1                 | 0                 | 0                 |
| CHR6       | 1                     | 1                 | 0                 | 0                 |
| CHR7       | 10                    | 10                | 0                 | 0                 |
| CHR8       | 25                    | 18                | 7                 | 0                 |
| CHR9       | 8                     | 6                 | 2                 | 0                 |

Next the coordinates for the GIMAP genes in the reference annotation were analyzed. For subsetting the sequences from the chromosomes, a window of 1mb was extracted on either side.

Table 2: Coordinates for GIMAP Gene Extraction on Each Chromosome

| Chromosome ID | Chromosome Number | Start    | End       | Extracted Range    |
|---------------|-------------------|----------|-----------|--------------------|
| NC_035781.1   | CHR2              | 42704802 | 42757768  | 43704802-43757768  |
| NC_035783.1   | CHR4              | 40487217 | 41243082  | 41487217-42243082  |
| NC_035784.1   | CHR5              | 89744459 | 89751167  | 90744459-90751167  |
| NC_035785.1   | CHR6              | 12854692 | 12928162  | 13854692-13928162  |
| NC_035786.1   | CHR7              | 14027963 | 54588883  | 15027963-55588883  |
| NC_035787.1   | CHR8              | 5914739  | 71430550  | 6914739-72430550   |
| NC_035788.1   | CHR9              | 29445429 | 103317686 | 30445429-104317686 |


# STEP 9. Analyze structural variants using LUMPY

LUMPY will be used to detect structural variants and BIC-Seq2 will be used to look at copy number variants.

Analyses were structured based on Greer et al., 2017.

"Using the conventional WGS data as input, tumor SVs were detected using LumPy and somatic copy number variants (CNVs) were detected using BICseq2 [26, 27]. LumPy was run using the lumpyexpress executable with default parameters, and the output VCF file was parsed to bed format for further processing. For copy number calling, BICseq2 first removes potential biases from the se- quencing data (BICseq2-norm v0.2.4) and subsequently calls CNVs from the normalized data (BICseq2-seg v0.7.2). The lambda parameter supplied to BICseq2-seg tunes the smoothness of the resulting CNV profile; a lambda value of 30 was used to call CNVs for the primary tumor and metastatic samples. Amplifications and deletions were called as segments with tumor/normal copy number ratios greater than 1.25 and less than 0.95, respectively.

Tutorial with LUMPY to find structural variants: http://bioinformatics-ca.github.io/bioinformatics_for_cancer_genomics_2016/rearrangement, https://mississippi.snv.jussieu.fr/u/drosofff/w/constructed-lumpy-workflow-imported-from-uploaded-file,
http://ngseasy.readthedocs.io/en/latest/containerized/ngseasy_dockerfiles/ngseasy_lumpy/README/

First we will prepare the data for LUMPY input by making BAM files with only discordant read pairs and split reads. Discordant read pairs are those that do not map as expected (and may be variants). Remembers these reads must first be sorted in a bam file. Next we will create a bam file containing only split reads that were mapped with a large insertion or deletion in the alignment.  

Copy the commands below into a script called LUMPY_prep.sh. Make all the scripts executable first using `chmod u+x script`

```
#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes 1
#SBATCH --exclusive
#SBATCH --mail-user=erin_roberts@my.uri.edu
#SBATCH -o /data3/marine_diseases_lab/erin/CV_Gen_Reseq/discordant_output
#SBATCH -e /data3/marine_diseases_lab/erin/CV_Gen_Reseq/discordant_error
#SBATCH -D /data3/marine_diseases_lab/erin/CV_Gen_Reseq/
cd=/data3/marine_diseases_lab/erin/CV_Gen_Reseq

echo "START $(date)"
module load SAMtools/1.5-foss-2017a

#extract only discordant read pairs
F=/data3/marine_diseases_lab/erin/CV_Gen_Reseq
array1=($(ls *.F.bam | sed 's/.F.bam//g'))
for i in ${array1[@]}; do
  samtools view -b -F 1294 ${i}.F.bam > ${i}.discordants.bam
  samtools view -h ${i}.F.bam | ./extractSplitReads_BwaMem -i stdin | samtools view -Sb - > ${i}.sr.bam
  samtools view ${i}.F.bam |  tail -n+100000 | ./pairend_distro.py -r 150 -X 4 -N 10000 -o ${i}.lib1.histo
  echo "${i} done LUMPY preprocess"
done

#  IF there is a glitch with these, make sure you check that all files are correctly position sorted
#python check_sorting.py \

# Finally, we need to make sure that all of the files are position sorted
array1=($(ls *.F.bam | sed 's/.F.bam//g'))
for i in ${array1[@]}; do
  samtools sort ${i}.discordants.bam -o ${i}.discordants.pe.sort
  samtools sort ${i}.sr.bam -o ${i}.sr.sort
done

echo "DONE $(date)"

```

Now we can finally run the LUMPY command. We must put in the individual mean and stdev for each file so they run correctly. Create a list based on all of this information  Output will be parsed into BED format.

```
# Use these commands to create the arrays
grep "LUMPY" discordant_output | cut -d ' ' -f 1 > mean_stdev_array.txt
grep "mean" discordant_output | sed 's/  /,/g' > mean_lines.txt
```
Put the following commands into a SLURM script called LUMPY.sh
```
#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes 1
#SBATCH --exclusive
#SBATCH --mail-user=erin_roberts@my.uri.edu
#SBATCH -o /data3/marine_diseases_lab/erin/CV_Gen_Reseq/lumpy_output
#SBATCH -e /data3/marine_diseases_lab/erin/CV_Gen_Reseq/lumpy_error
#SBATCH -D /data3/marine_diseases_lab/erin/CV_Gen_Reseq/
cd=/data3/marine_diseases_lab/erin/CV_Gen_Reseq

module load BEDTools/2.26.0-foss-2016b
#need BEDtools to be downloaded in order to view the output
module load LUMPY/0.2.13-foss-2016b
module load SAMtools/1.5-foss-2017a

echo "START $(date)"

#this command will run LUMPY with both paired end and split reads
# pe indicates paired end options
# sr indicates split read options
# default parameters recommended by the writers will be used
# min_non_overlap set to the read length

F=/data3/marine_diseases_lab/erin/CV_Gen_Reseq
array=($(cat mean_stdev_array.txt))
array2=($(cat mean_lines.txt))
for ((i=0;i<${#array[@]};++i)); do
    lumpy \
        -mw 4 \
        -tt 0.0 \
        -pe \
        id:${array[i]},bam_file:${array[i]}.discordants.pe.sort,${array2[i]},histo_file:${array[i]}.lib1.histo,read_length:150,min_non_overlap:150,discordant_z:4,back_distance:20,weight:1,min_mapping_threshold:20\
        -sr \
        bam_file:${array[i]}.sr.sort,back_distance:20,weight:1,id:2,min_mapping_threshold:20 \
        > ${array[i]}.pesr.vcf
        echo "done ${array[i]}"
done

echo "done all $(date)"

# The previous LUMPY command will process all the files separately, we would like to process them
# all together using lumpyexpress
# Lumpyexpress requires comma separated lists of files
# ls *.F.bam | awk -vORS=, '{print $1 }'
# ls *.discordants.pe.sort | awk -vORS=, '{print $1 }'
# ls *.sr.sort | awk -vORS=, '{print $1 }'
# ls -m will output with commas
# created a new lumpyexpress.config file with the path to samtools added

lumpyexpress -B CL_1.F.bam,CL_2.F.bam,CL_3.F.bam,CL_4.F.bam,CL_5.F.bam,CL_6.F.bam,CLP_1.F.bam,CLP_2.F.bam,CLP_3.F.bam,CLP_4.F.bam,CLP_5.F.bam,CLP_6.F.bam,CS_1.F.bam,CS_2.F.bam,CS_3.F.bam,CS_5.F.bam,CS_6.F.bam,CS_7.F.bam,DEBY_1.F.bam,DEBY_2.F.bam,DEBY_3.F.bam,DEBY_4.F.bam,DEBY_5.F.bam,DEBY_6.F.bam,HC_1.F.bam,HC_3.F.bam,HC_4.F.bam,HC_5.F.bam,HC_6.F.bam,HC_7.F.bam,HC_VA_1.F.bam,HC_VA_2.F.bam,HC_VA_3.F.bam,HC_VA_4.F.bam,HC_VA_5.F.bam,HC_VA_6.F.bam,HG_HG0F2.F.bam,HG_HG2F1.F.bam,HG_HG2M5.F.bam,HI_1.F.bam,HI_2.F.bam,HI_3.F.bam,HI_4.F.bam,HI_5.F.bam,HI_6.F.bam,LM_1_pool.F.bam,LM_3.F.bam,LM_4.F.bam,LM_7.F.bam,LM_8.F.bam,LOLA_1.F.bam,LOLA_2.F.bam,LOLA_3.F.bam,LOLA_4.F.bam,LOLA_5.F.bam,LOLA_6.F.bam,NEH_1.F.bam,NEH_2.F.bam,NEH_3.F.bam,NEH_4.F.bam,NEH_5.F.bam,NEH_6.F.bam,NG_NH0H4.F.bam,NG_NH2F6.F.bam,NG_NH2F8.F.bam,NG_NH2M1.F.bam,OBOYS2_1.F.bam,OBOYS2_2.F.bam,OBOYS2_3.F.bam,OBOYS2_4.F.bam,OBOYS2_5.F.bam,OBOYS2_6.F.bam,SL_1.F.bam,SL_2.F.bam,SL_3.F.bam,SL_4.F.bam,SL_5.F.bam,SL_6.F.bam,SM_10.F.bam,SM_11.F.bam,SM_12.F.bam,SM_7.F.bam,SM_8.F.bam,SM_9.F.bam,UMFS_1.F.bam,UMFS_2.F.bam,UMFS_3.F.bam,UMFS_4.F.bam,UMFS_5.F.bam,UMFS_6.F.bam -S CL_1.sr.sort,CL_2.sr.sort,CL_3.sr.sort,CL_4.sr.sort,CL_5.sr.sort,CL_6.sr.sort,CLP_1.sr.sort,CLP_2.sr.sort,CLP_3.sr.sort,CLP_4.sr.sort,CLP_5.sr.sort,CLP_6.sr.sort,CS_1.sr.sort,CS_2.sr.sort,CS_3.sr.sort,CS_5.sr.sort,CS_6.sr.sort,CS_7.sr.sort,DEBY_1.sr.sort,DEBY_2.sr.sort,DEBY_3.sr.sort,DEBY_4.sr.sort,DEBY_5.sr.sort,DEBY_6.sr.sort,HC_1.sr.sort,HC_3.sr.sort,HC_4.sr.sort,HC_5.sr.sort,HC_6.sr.sort,HC_7.sr.sort,HC_VA_1.sr.sort,HC_VA_2.sr.sort,HC_VA_3.sr.sort,HC_VA_4.sr.sort,HC_VA_5.sr.sort,HC_VA_6.sr.sort,HG_HG0F2.sr.sort,HG_HG2F1.sr.sort,HG_HG2M5.sr.sort,HI_1.sr.sort,HI_2.sr.sort,HI_3.sr.sort,HI_4.sr.sort,HI_5.sr.sort,HI_6.sr.sort,LM_1_pool.sr.sort,LM_3.sr.sort,LM_4.sr.sort,LM_7.sr.sort,LM_8.sr.sort,LOLA_1.sr.sort,LOLA_2.sr.sort,LOLA_3.sr.sort,LOLA_4.sr.sort,LOLA_5.sr.sort,LOLA_6.sr.sort,NEH_1.sr.sort,NEH_2.sr.sort,NEH_3.sr.sort,NEH_4.sr.sort,NEH_5.sr.sort,NEH_6.sr.sort,NG_NH0H4.sr.sort,NG_NH2F6.sr.sort,NG_NH2F8.sr.sort,NG_NH2M1.sr.sort,OBOYS2_1.sr.sort,OBOYS2_2.sr.sort,OBOYS2_3.sr.sort,OBOYS2_4.sr.sort,OBOYS2_5.sr.sort,OBOYS2_6.sr.sort,SL_1.sr.sort,SL_2.sr.sort,SL_3.sr.sort,SL_4.sr.sort,SL_5.sr.sort,SL_6.sr.sort,SM_10.sr.sort,SM_11.sr.sort,SM_12.sr.sort,SM_7.sr.sort,SM_8.sr.sort,SM_9.sr.sort,UMFS_1.sr.sort,UMFS_2.sr.sort,UMFS_3.sr.sort,UMFS_4.sr.sort,UMFS_5.sr.sort,UMFS_6.sr.sort -D CL_1.discordants.pe.sort,CL_2.discordants.pe.sort,CL_3.discordants.pe.sort,CL_4.discordants.pe.sort,CL_5.discordants.pe.sort,CL_6.discordants.pe.sort,CLP_1.discordants.pe.sort,CLP_2.discordants.pe.sort,CLP_3.discordants.pe.sort,CLP_4.discordants.pe.sort,CLP_5.discordants.pe.sort,CLP_6.discordants.pe.sort,CS_1.discordants.pe.sort,CS_2.discordants.pe.sort,CS_3.discordants.pe.sort,CS_5.discordants.pe.sort,CS_6.discordants.pe.sort,CS_7.discordants.pe.sort,DEBY_1.discordants.pe.sort,DEBY_2.discordants.pe.sort,DEBY_3.discordants.pe.sort,DEBY_4.discordants.pe.sort,DEBY_5.discordants.pe.sort,DEBY_6.discordants.pe.sort,HC_1.discordants.pe.sort,HC_3.discordants.pe.sort,HC_4.discordants.pe.sort,HC_5.discordants.pe.sort,HC_6.discordants.pe.sort,HC_7.discordants.pe.sort,HC_VA_1.discordants.pe.sort,HC_VA_2.discordants.pe.sort,HC_VA_3.discordants.pe.sort,HC_VA_4.discordants.pe.sort,HC_VA_5.discordants.pe.sort,HC_VA_6.discordants.pe.sort,HG_HG0F2.discordants.pe.sort,HG_HG2F1.discordants.pe.sort,HG_HG2M5.discordants.pe.sort,HI_1.discordants.pe.sort,HI_2.discordants.pe.sort,HI_3.discordants.pe.sort,HI_4.discordants.pe.sort,HI_5.discordants.pe.sort,HI_6.discordants.pe.sort,LM_1_pool.discordants.pe.sort,LM_3.discordants.pe.sort,LM_4.discordants.pe.sort,LM_7.discordants.pe.sort,LM_8.discordants.pe.sort,LOLA_1.discordants.pe.sort,LOLA_2.discordants.pe.sort,LOLA_3.discordants.pe.sort,LOLA_4.discordants.pe.sort,LOLA_5.discordants.pe.sort,LOLA_6.discordants.pe.sort,NEH_1.discordants.pe.sort,NEH_2.discordants.pe.sort,NEH_3.discordants.pe.sort,NEH_4.discordants.pe.sort,NEH_5.discordants.pe.sort,NEH_6.discordants.pe.sort,NG_NH0H4.discordants.pe.sort,NG_NH2F6.discordants.pe.sort,NG_NH2F8.discordants.pe.sort,NG_NH2M1.discordants.pe.sort,OBOYS2_1.discordants.pe.sort,OBOYS2_2.discordants.pe.sort,OBOYS2_3.discordants.pe.sort,OBOYS2_4.discordants.pe.sort,OBOYS2_5.discordants.pe.sort,OBOYS2_6.discordants.pe.sort,SL_1.discordants.pe.sort,SL_2.discordants.pe.sort,SL_3.discordants.pe.sort,SL_4.discordants.pe.sort,SL_5.discordants.pe.sort,SL_6.discordants.pe.sort,SM_10.discordants.pe.sort,SM_11.discordants.pe.sort,SM_12.discordants.pe.sort,SM_7.discordants.pe.sort,SM_8.discordants.pe.sort,SM_9.discordants.pe.sort,UMFS_1.discordants.pe.sort,UMFS_2.discordants.pe.sort,UMFS_3.discordants.pe.sort,UMFS_4.discordants.pe.sort,UMFS_5.discordants.pe.sort,UMFS_6.discordants.pe.sort -t3 -K lumpyexpress.config -o full_lumpy_bam.vcf

```
# STEP 10: Call SVGenotypes using SVTyper

-This is available on Github at https://github.com/hall-lab/svtools/tree/master/svtools/bin/svtyper.
-From their website =
"SVTyper performs breakpoint genotyping of structural variants (SVs) using whole genome sequencing data. Users must supply a VCF file of sites to genotype (which may be generated by LUMPY) as well as a BAM/CRAM file of Illumina paired-end reads aligned with BWA-MEM. SVTyper assesses discordant and concordant reads from paired-end and split-read alignments to infer genotypes at each site. Algorithm details and benchmarking are described in Chiang et al., 2015."

-l flag creates a JSON file with essential metrics on a BAM file. SVTyper will sample the first N reads for the file (1 million by default) to parse the libraries, read groups, and insert size histograms. This can be done in the absence of a VCF file.

Load the following commands into a script called svtyper.sh. The

```
#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes 1
#SBATCH --exclusive
#SBATCH --mail-user=erin_roberts@my.uri.edu
#SBATCH -o /data3/marine_diseases_lab/erin/CV_Gen_Reseq/SVTyper_output
#SBATCH -e /data3/marine_diseases_lab/erin/CV_Gen_Reseq/SVTyper_error
#SBATCH -D /data3/marine_diseases_lab/erin/CV_Gen_Reseq/
cd=/data3/marine_diseases_lab/erin/CV_Gen_Reseq

echo "START $(date)"
module load svtyper/0.6.1-foss-2017a-Python-2.7.13

#Use the following code if you would like to run each file individually with the full bam file just generated
array1=($(ls *.F.bam | sed 's/.F.bam//g'))
for i in ${array1[@]}; do
  svtyper -i full_lumpy_bam.vcf -B ${i}.F.bam -S ${i}.sr.sort > ${i}.sv.gt2.vcf
  echo "done ${i} $(date)"
done

echo "DONE"

#Finally, we can merge the SVTyped files into a single file using the following recommended command: (from https://github.com/arq5x/lumpy-sv/issues/84)
ls *.sv.gt2.vcf > svtyper_vcfs.txt
cat full_lumpy_bam.vcf \
        | vcf_group_multiline.py \
        | scripts/vcf_paste.py \
            -q \
            -m -
             $(<svtyper_vcfs.txt) \
        > final.vcf
echo "done combine $(date)"

```

# STEP 11. Filter VCF file

The following steps were adapted from an excellent protocol developed by Jon Puritz. For more detailed information please see [http://ddocent.com/filtering/](http://ddocent.com/filtering/)
We will perform all of the following commands in single bash script called vcf_filter.sh

1. First we will use VCFTools to filter out any variants that have not been successfully genotyped to more than 50% of individuals ( `--max-missing 0.5`), those with a minor allele count of 3 (`--mac 3)`), and those with a quality score below 30 (`--minQ 30`)  

```
vcftools --vcf final.vcf --max-missing 0.5 --mac 3 --minQ 20 --recode --recode-INFO-all --out final.g5mac3

```

2. Getting rid of these first will help speed up this next command, which applies a minimum mean depth and a minimum depth for a genotype call. Genotypes will be called if they have atleast three reads.

```
vcftools --vcf final.firstfilter.recode.vcf --minDP 3 --recode --recode-INFO-all --out final.g5mac3dp3
```

3. Now we can remove individuals that have a lot of missing data.

```
# The output of this file will be called out.imiss
vcftools --vcf final.g5mac3dp3.recode.vcf --missing-indv

```
4. We can now plot a histogram of individuals that are missing a lot of data using the following command (taken from http://ddocent.com/filtering/ by Jon Puritz).

```
mawk '!/IN/' out.imiss | cut -f5 > totalmissing
gnuplot << \EOF
set terminal dumb size 120, 30
set autoscale
unset label
set title "Histogram of % missing data per individual"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
#set yr [0:100000]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF
```
5. We can create a list with more than 50% missing data using the following mawk command.
```
mawk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv
```
6. This can then be piped into VCFtools so that the low coverage individuals can be removed.
```
vcftools --vcf final.g5mac3dp3.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out final.g5mac3dp3dplm

```
7. Next, we need to restrict the data to those variants that are called in a high percentage of individuals using a genotype call rate of 95% (`--max-missing 0.95`) and then filter based on the mean depth of genotypes of 20% (`--min-meanDP 20`).
```
vcftools --vcf final.g5mac3dp3dplm.recode.vcf --max-missing 0.95 --maf 0.05 --recode --recode-INFO-all --out final.DP3g95maf05 --min-meanDP 20
```

8. Finally, because we are analyzing several individuals, we need to apply a population specific filter. To do this we first need to create what's called a "popmap" file. This file contains two tab separated columns. Lets create one based on our data.

```
array2=($(ls *.F.cleaned.fq. | sed 's/.F.cleaned.fq.gz//g'))

for i in ${array2[@]}; do
  echo -e "${i}:${i}" | awk  '{gsub(":","\t",$0); print;}' >> popmap
done
awk '{split($2,a,/_/);$2=a[1]}1' popmap > popmap_final  
sed 's/ /\t/g' popmap_final > popmap_final_TD

#note: manually added _VA to HC_VA lines using nano
```
8a. Now we need to create 9 lists that have the individual names of each population. We can do this with the following commands that can be used for whatever the unique names of your populations might be.

```
cut -f 2 popmap_final_TD | sort | uniq > unique.txt

#!/bin/bash
while read -r i; do
  echo $i > $i.keep
done < unique.txt

```  
8b. We can use these files to estimate the missing data for loci in each population using vcftools.

```
array3=($(ls *.keep | sed 's'/.keep//')
for i in ${array3[@]}; do
  vcftools --vcf final.DP3g95maf05.recode.vcf --keep ${i}.keep --missing-site --out ${i}
done
```
8c. The output of the above commands outputs files called `*.lmiss` whose last column lists the percentage of missing data for that locus. We can merge all of these files to create a list of loci that have 10% missing data or more to remove.

```
cat CL.lmist CLP.lmiss CS.lmiss HC.lmiss HC_VA.lmiss HI.lmiss LM.lmiss SL.lmiss SM.lmiss | mawk '!/CHR/' | mawk '$6 > 0.1' | cut -f1,2 >> badloci

```
8d. Finally, we can pipe this back into VCFTools to remove any of those bad loci

```
vcftools --vcf final.DP3g95maf05.recode.vcf --exclude-positions badloci --recode --recode-INFO-all --out final.DP3g95p5maf05

```
9. Next we can apply a filter that remove sites that have reads from both strands. The filter command is going to keep loci that have over 100 times more forward alternate reads than reverse alternate reads and 100 times more forward reference reads than reverse reference reads along with the reciprocal.

```
vcffilter -f "SAF / SAR > 100 & SRF / SRR > 100 | SAR / SAF > 100 & SRR / SRF > 100" -s final.DP3g95p5maf05.recode.vcf > final.DP3g95p5maf05.fil1.vcf

# To investigate how many reads were remove we can see how many lines remain
mawk '!/#/' final.DP3g95p5maf05.fil1.vcf | wc -l

```
10. Apply a filter to account for high coverage causing an inflated locus quality score

Heng Li found that in whole genome samples, high coverage can lead to inflated locus quality scores. Based on this, Jon Puritz suggests the following filter to remove any locus that has a quality score below 1/4 of the depth. Because we are not working with RADseq data, we will not implement the second filter he suggest to recalculate mean depth.

```
vcffilter -f "QUAL / DP > 0.25" final.DP3g95p5maf05.fil1.vcf > final.DP3g95p5maf05.fil2.vcf
```

# STEP 12: Use ANNOVAR to Annotate structural variants
1. Download and install ANNOVAR on computer

2. Prepare annotation files for Crassostrea virginica, available at https://www.ncbi.nlm.nih.gov/genome/?term=txid6565[orgn]
-Download reference genomic FASTA file
-Download reference GFF3 format file
-Put both in same directory
```
mkdir ANNOVAR_input
```
3. Use the gtftoGenePhred tool to convert the GFF to a GenePhred file  

```
gtfToGenePred -genePredExt ref_C_virginica-3.0_top_level.gff3 CV_refGene.txt
```

4. Generate a transcript FASTA file with the ANNOVAR provided script
```
perl ../retrieve_seq_from_fasta.pl --format refGene --seqfile
  GCA_002022765.4_C_virginica-3.0_genomic.fna CV_refGene.txt --out CV_refGeneMrna.fa
```
Now the annotation database files are ready for use and can use protocol in Yang 2015 at 2B(ii), using '--buildver' argument set to 'CV'

5. Annotate variants with the C. virginica reference genome annotation

```
perl table_annovar.pl final.DP3g95p5maf05.fil2.vcf ANNOVAR_input/ --vcfinput --outfile CV_final --buildver CV --protocol refGene --operation g
```

The final output is a text file from the previous command.

6. Use table_annovar.pl to directly take VCF as input and generate a VCF_file and a tabular output file  with the INFO column containing various ANNOVAR annotations
(see http://annovar.openbioinformatics.org/en/latest/misc/accessory/#table_annovar-automated-execution-of-multiple-annotation-tasks)

table_annovar.pl final.DP3g95p5maf05.fil2.vcf ANNOVAR_input/ -buildver CV -out CVanno -remove -protocol refGene -operation g -nastring . -vcfinput

7. Subset GIMAP regions of output VCF file with VCFtools

The coordinates used are listed below.

```
$ cat GIMAP_coordinates.txt
NC_035781.1	43704802	43757768
NC_035783.1	41487217	42243082
NC_035784.1	90744459	90751167
NC_035785.1	13854692	13928162
NC_035786.1	15027963	55588883
NC_035787.1	6914739		72430550
NC_035788.1	30445429	104317686
```

```
#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes 1
#SBATCH --exclusive
#SBATCH --mail-user=erin_roberts@my.uri.edu
#SBATCH -o /data3/marine_diseases_lab/erin/CV_Gen_Reseq/vcf_GIMAP_output
#SBATCH -e /data3/marine_diseases_lab/erin/CV_Gen_Reseq/vcf_GIMAP_error
#SBATCH -D /data3/marine_diseases_lab/erin/CV_Gen_Reseq/
cd=/data3/marine_diseases_lab/erin/CV_Gen_Reseq

echo "START $(date)"
module load VCFtools/0.1.14-foss-2017a-Perl-5.24.1

F=/data3/marine_diseases_lab/erin/CV_Gen_Reseq
  vcftools --vcf CV_anno.vcf --chr NC_035781.1 --from-bp 43704802 --to-bp 43757768 --recode --recode-INFO-all --out CV_anno_GIMAP_subset_NC_035781.1
  vcftools --vcf CV_anno.vcf --chr NC_035783.1 --from-bp 41487217 --to-bp 42243082 --recode --recode-INFO-all --out CV_anno_GIMAP_subset_NC_035783.1
  vcftools --vcf CV_anno.vcf--chr NC_035784.1 --from-bp 90744459 --to-bp 90751167 --recode --recode-INFO-all --out CV_anno_GIMAP_subset_NC_035784.1
  vcftools --vcf CV_anno.vcf --chr NC_035785.1 --from-bp 13854692 --to-bp 13928162 --recode --recode-INFO-all --out CV_anno_GIMAP_subset_NC_035785.1
  vcftools --vcf CV_anno.vcf --chr NC_035786.1 --from-bp 15027963 --to-bp 55588883 --recode --recode-INFO-all --out CV_anno_GIMAP_subset_NC_035786.1
  vcftools --vcf CV_anno.vcf --chr NC_035787.1 --from-bp 6914739 --to-bp 72430550 --recode --recode-INFO-all --out CV_anno_GIMAP_subset_NC_035787.1
  vcftools --vcf CV_anno.vcf --chr NC_035788.1 --from-bp 30445429 --to-bp 104317686 --recode --recode-INFO-all --out CV_anno_GIMAP_subset_NC_035788.1


echo "done ${i}"

```

# STEP 13: Identify Consensus Genotypes using SURVIVOR SVTyped files
(https://github.com/fritzsedlazeck/SURVIVOR/wiki)
1. Install SURVIVOR

2. Merge all VCF files
```
ls *vcf > sample_files
```
3. Obtain a consensus call set

The following command will  merge all the vcf files specified in sample_files together using a maximum allowed distance of 1kb. Furthermore we ask SURVIVOR only to report calls supported by 2 callers and they have to agree on the type (1) and on the strand (1) of the SV. Note you can change this behavior by altering the numbers from 1 to e.g. 0. In addition, we told SURVIVOR to only compare SV that are at least 30bp long and print the output in sample_merged.vcf.

```
./SURVIVOR merge sample_files 1000 2 1 1 0 30 sample_merged.vcf
```
4.  Use genomic coordinates mapping to GIMAP genes above to identify GIMAP genes in merged set


# STEP 14: Analyses to characterize structural variants mapped

1. Calculate the total number of variants and the number of variants within each population

We can do this by using vcf-annotate with --fill-type
```
zcat SNP.DP3g95p5maf05.recode.vcf | vcf-annotate --fill-type | grep -oP "TYPE=\w+" | sort | uniq -c > VCF_filltype

grep "snp" VCF_filltype > SNP_VCF_filltype
grep "INDEL" VCF_filltype > indel_VCF_filltype
grep "dup" VCF_filltype > dup_VCF_filltype
```

2. Calculating allele frequencies and the percentage of exclusive variants

2a. Calculate allele frequency
```
# To get the frequency of each allele over all individuals in a VCF file, you can use the -freq argument
vcftools --vcf final.DP3g95p5maf05.recode.vcf --freq --out allele_freq_all

# To calculate for individual populations

array3=($(ls *.keep | sed 's'/.keep//')
for i in ${array3[@]}; do
  vcftools --vcf final.DP3g95p5maf05.recode.vcf--freq --keep ${i}.keep --out ${i}.allele_freq
done

```
2b. Exclusive variants are defined as those that are polymorphic in only one population.
We will do this using the vcftools module vcf-compare

```
#compare the vcf files using the -g option to also output summary stats
vcf-compare -H -g A.vcf.gz B.vcf.gz C.vcf.gz
```

3. Calculate the level of nucleotide diversity using VCFTools

-3a: Calculate for each population separately
In order to reduce the computational intensity the sliding window has been set to 10,000.
```
array3=($(ls *.keep | sed 's'/.keep//')
for i in ${array3[@]}; do
  vcftools --vcf final.DP3g95p5maf05.recode.vcf  --keep ${i}.keep --window-pi 100000 --out = ${i.ND}
done
```

-3b: Use R to calculate a mean pi for all individuals

```
install.packages("fields")
library(fields)
binstats <- stats.bin(dat$POS,dat$PI,N=100)
matplot( binstats$centers, t(binstats$stats[ c("mean", "median","Q1", "Q3"),]), type="l",lty=c(1,2,2,2), col=c('red','blue','green','purple'), ylab="Pi diversity")

```

# STEP 15: Calculate the Fst across different populations

Now we will calculate the Fst statistic between individuals of different populations to get a measure of population differentiation.
```
vcftools --vcf SNP.DP3g95p5maf05.recode.vcf --weir-fst-pop CL.txt --weir-fst-pop CLP.txt --weir-fst-pop CS.txt --weir-fst-pop HC.txt --weir-fst-pop HC_VA.txt --weir-fst-pop HI.txt --weir-fst-pop LM.txt --weir-fst-pop SL.txt --weir-fst-pop SM.txt --out Fst_all_pop


```

# STEP 16: Assess genetic structure by using a PCA through the package adegenet

To carry out the code below, we need to first create a "strata" file that has three columns. The first column is the name of the individual, the second is the population, and the third is the library. We have the first two columns already in our popmap_final_TD file (though with no headings).

We can add on the library type (libraryA) with the following commands

```
#create file with repeated word with same number of lines as popmap file
printf 'LibA\n%.0s' {1..54}

#paste the two together
paste popmap_final_TD libfile > strata

```
Now we can add plot our data using adegenet in R.

```
library(adegenet)
library(vcfR)
library(hierfstat)
library(ape)

#loading in the vcf file
oyster_vcf <- read.vcfR("SNP.DP3g95p5maf05.recode.vcf")

# the genind object store allelic data for individuals as integer allele counts
oyster_genind <- vcfR2genind(oyster_vcf)
strata<- read.table("strata", header=TRUE)
strata_df <- data.frame(strata)
strata(oyster_genind) <- strata_df

setPop(oyster_genind) <- ~Population

#Test Population Structure
fstat(oyster_genind)

oyster.pairwiseFst <- pairwise.fst(oyster_genind, res.type="matrix")

#We can use these genetic distance values to calculate a tree using ape
oyster.tree <- nj(oyster.pairwiseFst)
plot(oyster.tree, type="unr", tip.col=funky(nPop(oyster_genind)), font=2)
annot <- round(oyster.tree$edge.length,2)
edgelabels(annot[annot>0], which(annot>0), frame="n")
add.scale.bar()

# we can test for the presence of population structure using Goudet's *G* statistic.
oyster.gtest <- gstat.randtest(oyster_genind)
oyster.gtest


```

# STEP 17: Use PCAdapt to visualize sample clustering by population location

We will now use PCAdapt to plot clustering of populations in R.

```
R
library(pcadapt)
variants <- read.pcadapt("SNP.DP3g95p5maf05.recode.vcf", type = "vcf")

#In this first step pick a minimum number of possible PCs, we'll go high with 20
PCs <- pcadapt(input = filename, K = 20)

#Now we will plot the percentage of variance explained by each PC using a screeplot
plot(PCs, option = "screeplot")

#plot for the fist ten
plot(PCs, option = "screeplot", K = 10)

# Create list with population designations
poplist.names <- c(rep("CL", 20),rep("CLP", 20),rep("CS", 20), rep("HC,20)),
rep("HC_VA", 20), rep("HI", 20), rep("LM", 20), rep("SL", 20), rep("SM", 20))

#Display individual scores on each of the PCs and produce a score plot
plot(PCs, option = "scores", pop = poplist.names)
#try 2 and 3
plot(PCs, option = "scores", i = 2, j = 3, pop = poplist.names)
#try 3 and 4
plot(PCs, option = "scores", i = 2, j = 3, pop = poplist.names)

# redo with appropriate number of PCs
PC2<- pcadapt(variants, K=3)

#output the summary statistics
summary(PC2)

#putatively investigate outliers with a Manhattan plot
plot(PC2, option="manhattan")

```

# Bibliography and Works Cited
-Website about dDocent filtering
https://github.com/jpuritz/BIO_594_2018/blob/master/Exercises/Week10/EecSeq_code.md
-References about SV detection
https://bpa-csiro-workshops.github.io/btp-manuals-md/modules/cancer-module-snv/snv/
http://bioinformatics-ca.github.io/bioinformatics_for_cancer_genomics_2016/rearrangement
https://bcbio.wordpress.com/tag/lumpy/
http://ngseasy.readthedocs.io/en/latest/containerized/ngseasy_dockerfiles/ngseasy_lumpy/README/
https://mississippi.snv.jussieu.fr/u/drosofff/w/constructed-lumpy-workflow-imported-from-uploaded-file
https://raw.githubusercontent.com/bioinformatics-ca/2015_workshops/master/BiCG_2015/BiCG_2015_Module3_Lab2.txt
Brandler, W. M., D. Antaki, M. Gujral, A. Noor, G. Rosanio, T. R. Chapman, D. J. Barrera, G. N. Lin, D. Malhotra, A. C. Watts, L. C. Wong, J. A. Estabillo, T. E.   Gadomski, O. Hong, K. V. F. Fajardo, A. Bhandari, R. Owen, M. Baughn, J. Yuan, T. Solomon, A. G. Moyzis, M. S. Maile, S. J. Sanders, G. E. Reiner, K. K. Vaux, C.   M. Strom, K. Zhang, A. R. Muotri, N. Akshoomoff, S. M. Leal, K. Pierce, E. Courchesne, L. M. Iakoucheva, C. Corsello, and J. Sebat. 2016. Frequency and Complexity of de Novo Structural Mutation in Autism. Am. J. Hum. Genet. 98:667–679. Available from: http://dx.doi.org/10.1016/j.ajhg.2016.02.018

Greer, S. U., L. D. Nadauld, B. T. Lau, J. Chen, C. Wood-Bouwens, J. M. Ford, C. J. Kuo, and H. P. Ji. 2017. Linked read sequencing resolves complex genomic rearrangements in gastric cancer metastases. Genome Med. 9:1–17.

-vcftools isec code help
https://www.biostars.org/p/140263/
-ANNOVAR help
http://annovar.openbioinformatics.org/en/latest/user-guide/region/
http://www.ngscourse.org/Course_Materials/variant_annotation/tutorial/annovar.html
http://annovar.openbioinformatics.org/en/latest/user-guide/input/#annovar-input-file
http://annovar.openbioinformatics.org/en/latest/user-guide/gene/#what-about-gff3-file-for-new-species
Yang, H., and K. Wang. 2015. Genomic variant annotation and prioritization with ANNOVAR and wANNOVAR. Nat. Protoc. 10:1556–1566. Available from: http://dx.doi.org/10.1038/nprot.2015.105

MultiQC: Summarize analysis results for multiple tools and samples in a single report
Philip Ewels, Måns Magnusson, Sverker Lundin and Max Käller
Bioinformatics (2016)
doi: 10.1093/bioinformatics/btw354
PMID: 27312411


https://www.biostars.org/p/75489/
