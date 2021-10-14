#!/bin/bash
workspace=$1
map_fastq=$workspace/map_fastq
dnadata=$map_fastq/DNA_fastq
bam_out=$workspace/MappedResult/bwa_bam
reference=~/reference/human_DNA/liheng_Genome_reference/liheng_human_hg38_bwa
umi=$workspace/umidedup
mkdir $umi/dedup_fastq/dedup_dna
dedupdna=$umi/dedup_fastq/dedup_dna
mkdir $umi/umi_mapped
mkdir $umi/umi_mapped/umi_bwabam
mkdir $umi/umi_mapped/umi_bwabed
dna_umibam=$umi/umi_mapped/umi_bwabam
dna_bed=$umi/umi_mapped/umi_bwabed
## 3. umi deduplication and reads extraction
cd $dnadata
unpigz *gz
cd $bam_out
for i in `ls *bam`
do
  echo "umitools dedup sample $i"
  samtools view -hbf 64 $i > R2.bam -@ 24
  samtools index R2.bam -@ 24
  umi_tools dedup -I R2.bam -S dedup.bam --umi-separator=":umi_"
#  samtools fastq -1 paired1.fq -2 paired2.fq -n dedup.bam -s single.fq
  bedtools bamtofastq -i dedup.bam -fq single.fq
  repair.sh in1=$dnadata/${i%_dna.bam}\_dfp_R1.fastq in2=single.fq out1=$dedupdna/${i%_dna.bam}_dedup_R1.fastq out2=$dedupdna/${i%_dna.bam}_dedup_R2.fastq outsingle=single_rna
  rm dedup.bam single_rna R2.bam* single.fq
done
## 4. gzip to save space
cd $dnadata
pigz *fastq
cd $dedupdna
pigz *fastq
## 5. de-dumi mapping
cd $dedupdna
for f in `ls *_dedup_R1.fastq.gz`
do
  echo "dna bwa_to_human Mapping ${f%_dedup_R1.fastq.gz}"
  bwa mem -t 24 -M $reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna $f ${f%_dedup_R1.fastq.gz}_dedup_R2.fastq.gz > $dna_umibam/temp.sam
  samtools view -@ 24 -S -b -q 59 $dna_umibam/temp.sam > $dna_umibam/temp.bam
  samtools sort $dna_umibam/temp.bam -o $dna_umibam/${f%_dedup_R1.fastq.gz}_dedup.bam -@ 24
  bamToBed -i $dna_umibam/${f%_dedup_R1.fastq.gz}_dedup.bam > $dna_bed/${f%_dedup_R1.fastq.gz}_dedup.bed
  pigz $dna_bed/${f%_dedup_R1.fastq.gz}_dedup.bed
  rm $dna_umibam/temp*
done
