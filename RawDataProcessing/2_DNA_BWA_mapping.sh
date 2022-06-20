#!/bin/bash
## This script will map all DNA fastq with BWA mem
workspace=$1
map_fastq=$workspace/map_fastq
dnadata=$map_fastq/DNA_fastq
mkdir $workspace/MappedResult/bwa_bam
bam_out=$workspace/MappedResult/bwa_bam
reference=$PATH_to_your_DNA_reference
## BE CAURSION: BWA reverse mapping: treat illumia read 1 as read 2; read 2 as read 1. This is because we need deduplicate UMI with read2
cd $dnadata
for f in `ls *_dfp_R1.fastq.gz`
do
  echo "dna bwa_to_human Mapping ${f%_dfp_R1.fastq.gz}"
  bwa mem -t 24 -M $PATH_to_your_DNA_reference ${f%_dfp_R1.fastq.gz}_dfp_R2.fastq.gz $f > $bam_out/temp.sam
  samtools view -@ 24 -S -b -q 59 $bam_out/temp.sam > $bam_out/temp.bam
  samtools sort $bam_out/temp.bam -o $bam_out/${f%_dfp_R1.fastq.gz}_dna.bam -@ 24
  rm $bam_out/temp*
done
