#07_C_gigas_RNA_pipeline_Salmon

#this script uses the program Salmon to calculate Transcript abundance

Module load Salmon/0.4.0-foss-2016b-Python-2.7.12  
C=/data3/marine_diseases_lab/erin/Crassostrea_gigas_reference_genome/


#two modes of operation: Indexing, Quantification


#Salmon command to quantify transcripts in a quasi mapping mode)
F=/data3/marine_diseases_lab/erin/Bio_project_SRA/PE_fastq
C=/data3/marine_diseases_lab/erin/Crassostrea_gigas_reference_genome/

array1=($(ls $F/*_1.fq.trim.filter))
array2=($(ls $F/*_2.fq.trim.filter))

#create salmon index for 
#If you want to use Salmon in quasi-mapping-based mode, then you first have to build an Salmon index for your transcriptome. 
#Assume that transcripts.fa contains the set of transcripts you wish to quantify. First, you run the Salmon indexer:

#First create an index out of a normal control file that seems normal and use that to create and index to compare to

salmon index -t transcripts.fa -i transcripts_index --type quasi -k 31


salmon index -t transcript.fa -i C/C_gigas_GTF_index --type quasi -k 31

#Salmon command to quantify paired end reads
salmon quant -i transcripts_index -l <LIBTYPE> -1 reads1.fq -2 reads2.fq -o transcripts_quant

#Salmon command to quantify multiple reads (if I want to put replicates of runs in here. For statistical power purposes I choose not to aggregate reads)

salmon quant -i C/C_gigas_GTF_index -l A -1 lib_1_1.fq lib_2_1.fq -2 lib_1_2.fq lib_2_2.fq -o out

# -l A means salmon will infer the library
-Can quantify Multiple reads together as if there are in one run


#References:
#http://salmon.readthedocs.io/en/latest/salmon.html