#!/bin/bash
workspace=$1
map_fastq=$workspace/map_fastq
dnadata=$map_fastq/DNA_fastq
bam_out=$workspace/MappedResult/bwa_bam
reference=$PATH_to_your_DNA_reference
umi=$workspace/umidedup
mkdir $umi
mkdir $umi/dedup_fastq
mkdir $umi/dedup_fastq/dedup_dna
dedupdna=$umi/dedup_fastq/dedup_dna
mkdir $umi/umi_mapped
mkdir $umi/umi_mapped/umi_bwabam
mkdir $umi/umi_mapped/umi_bwabed
dna_umibam=$umi/umi_mapped/umi_bwabam
dna_bed=$umi/umi_mapped/umi_bwabed
## 3. umi deduplication and reads extraction
cd $bam_out
for i in `ls *bam`
do
  echo "umitools dedup sample $i"
  ## Extract mapped read2
  samtools view -hbf 64 $i > R2.bam -@ 24
  samtools index R2.bam -@ 24
  ## Use umi_tools to perform deduplication
  umi_tools dedup -I R2.bam -S dedup.bam --umi-separator=":umi_"
  ## convert bam to read2 fastq file
  bedtools bamtofastq -i dedup.bam -fq single.fq
  pigz single.fq
  ## Grep read1 based on read2 fastq head name.
  repair.sh in1=$dnadata/${i%_dna.bam}\_dfp_R1.fastq.gz in2=single.fq.gz out1=$dedupdna/${i%_dna.bam}_dedup_R1.fastq.gz out2=$dedupdna/${i%_dna.bam}_dedup_R2.fastq.gz outsingle=single_dna
  rm dedup.bam single_dna R2.bam* single.fq*
done
## 5. de-dumi mapping
cd $dedupdna
for f in `ls *_dedup_R1.fastq.gz`
do
  echo "dna bwa_to_human Mapping ${f%_dedup_R1.fastq.gz}"
  bwa mem -t 24 -M $PATH_to_your_DNA_reference $f ${f%_dedup_R1.fastq.gz}_dedup_R2.fastq.gz > $dna_umibam/temp.sam
  samtools view -@ 24 -S -b -q 59 $dna_umibam/temp.sam > $dna_umibam/temp.bam
  samtools sort $dna_umibam/temp.bam -o $dna_umibam/${f%_dedup_R1.fastq.gz}_dedup.bam -@ 24
  ## convert bam to bed files
  bamToBed -i $dna_umibam/${f%_dedup_R1.fastq.gz}_dedup.bam > $dna_bed/${f%_dedup_R1.fastq.gz}_dedup.bed
  ## compress bed files
  pigz $dna_bed/${f%_dedup_R1.fastq.gz}_dedup.bed
  rm $dna_umibam/temp*
done
