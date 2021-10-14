#!/bin/bash
workspace=$1
map_fastq=$workspace/map_fastq
mkdir $map_fastq/Clean_rna
rnadata=$map_fastq/RNA_fastq
cd $rnadata
for i in `ls *R1.fastq.gz`
do
  bowtie2 --threads 12 --sensitive-local -x ~/reference/Ecoli/bowtie_ecoli/dh10b -1 $i -2 ${i%_R1.fastq.gz}_R2.fastq.gz -S temp.sam
  samtools view -@ 12 -b -F 2 temp.sam > unmapped.bam
  samtools fastq -@ 12 -1 paired1.fastq -2 paired2.fastq -n unmapped.bam -s single1.fastq
  rm temp.sam unmapped.bam single1.fastq
  pigz paired1.fastq
  pigz paired2.fastq
  mv paired1.fastq.gz $map_fastq/Clean_rna/${i%_rfp_R1.fastq.gz}_rfp_R1.fastq.gz
  mv paired2.fastq.gz $map_fastq/Clean_rna/${i%_rfp_R1.fastq.gz}_rfp_R2.fastq.gz
done

for i in `ls *R1.fastq.gz`
do
  mv $i ${i%_cl_R1.fastq.gz}_rfp_R1.fastq.gz
  mv ${i%_cl_R1.fastq.gz}_cl_R2.fastq.gz ${i%_cl_R1.fastq.gz}_rfp_R2.fastq.gz
done
