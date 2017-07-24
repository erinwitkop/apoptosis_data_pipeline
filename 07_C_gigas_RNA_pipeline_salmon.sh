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
#Choose reference from normal pooled adult tissue samples under no stress, needs to be SE for this analysis, and needs to be long enough
#so that lots can be compared to it

#Sample info:

	#C gigas Adult, mixture tissues, no stress BioProject Accession: PRJNA194084  SRA ID: SRR796589
	#Read type: SE, PCR Length: 3.1G, Sequencer Illumina HiSeq 2000

salmon index -t SRR796589.fastq  -i C/C_gigas_GTF_index --type quasi -k 31

#Salmon command to quantify paired end reads
for i in ${array1[@]}; do
	salmon quant -i C/C_gigas_GTF_index -l A -1 ${i} -2 $(echo ${i}|sed s/_1/_2/) -o ${i}.transcripts_quant
done 

# -l A means salmon will infer the library



#References:
#http://salmon.readthedocs.io/en/latest/salmon.html