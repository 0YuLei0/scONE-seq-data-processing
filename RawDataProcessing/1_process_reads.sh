#!/bin/bash
## The worksapce is a folder you created for scONE-seq
workspace=$1
## In this folder, you need creat a folder named raw to store the fastq files from seqeuncing
mkdir $workspace/catfiles
mkdir $workspace/map_fastq
map_fastq=$workspace/map_fastq
mkdir $map_fastq/DNA_fastq
mkdir $map_fastq/RNA_fastq
mkdir $map_fastq/Un_fastq
## This folder will store the fastq qc summaries for each cell.
mkdir $map_fastq/fastp_report
report=$map_fastq/fastp_report
cd $workspace/raw
for i in `ls *read1.fastq.gz`
do
  sh $PATH/Separation_bash.sh $i ${i%_read1.fastq.gz}_read2.fastq.gz
  ## rename dna rna fastq files and move them to the right folders
  mv fp_r1.fastq.gz $workspace/catfiles/${i%_read1.fastq.gz}_fp_R1.fastq.gz
  mv fp_r2.fastq.gz $workspace/catfiles/${i%_read1.fastq.gz}_fp_R2.fastq.gz
  mv dna_fp_r1.fastq.gz $map_fastq/DNA_fastq/${i%_read1.fastq.gz}_dfp_R1.fastq.gz
  mv dna_fp_r2.fastq.gz $map_fastq/DNA_fastq/${i%_read1.fastq.gz}_dfp_R2.fastq.gz
  mv rna_fp_r1.fastq.gz $map_fastq/RNA_fastq/${i%_read1.fastq.gz}_rfp_R1.fastq.gz
  mv rna_fp_r2.fastq.gz $map_fastq/RNA_fastq/${i%_read1.fastq.gz}_rfp_R2.fastq.gz
  mv trimmed-unknown.R1.fastq.gz $map_fastq/Un_fastq/${i%_read1.fastq.gz}_ufp_R1.fastq.gz
  mv trimmed-unknown.R2.fastq.gz $map_fastq/Un_fastq/${i%_read1.fastq.gz}_ufp_R2.fastq.gz
  rm trimmed*
  mkdir $report/${i%_read1.fastq.gz}
  mv *html *json $report/${i%_read1.fastq.gz}/
done
