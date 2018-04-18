#Commands for Paired End Read Preprocessing using BBTools
#This script runs all of BBTools steps in a single loop for each file and is formatted for zipped files

# Specify current path (just for extra security)
F=/home/eroberts/RNA-seq #or whatever your PATH is

#going to make one array variables and then iterate through them as an index
array1=($(ls $F/*_1.fastq | sed 's/_1.fastq//g'))

for i in ${array1[@]}; do
	gunzip ${i}_1.fastq.gz
  gunzip ${i}_2.fastq.gz
	/usr/local/bin/bbmap/bbduk.sh in1=${i}_1.fastq  out1=${i}_1.fastq.clean in2=${i}_2.fastq out2=${i}_2.fastq.clean ref=$F/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
	echo "adapter trimming ${i}" $(date)
  /usr/local/bin/bbmap/bbduk.sh in1=${i}_1.fastq.clean out1=${i}_1.fastq.clean.trim in2=${i}_2.fastq.clean out2=${i}_2.fastq.clean.trim qtrim=rl trimq=20
  echo "quality trimming ${i}" $(date)
  /usr/local/bin/bbmap/bbduk.sh in1=${i}_1.fastq.clean.trim out1=${i}_1.fastq.clean.trim.filter in2=${i}_2.fastq.clean.trim out2=${i}_2.fastq.clean.trim.filter maq=10
  echo "STOP" $(date)
  echo "quality filtering ${i}" $(date)
  /usr/local/bin/bbmap/bbduk.sh in1=${i}_1.fastq.clean.trim.filter in2=${i}_2.fastq.clean.trim.filter  bhist=${i}.b.hist qhist=${i}.q.hist gchist=${i}.gc.hist lhist=${i}.l.hist gcbins=auto
  echo "histogram DONE" $(date)
	gzip ${i}_1.fastq
	gzip ${i}_2.fastq
done
