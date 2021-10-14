#!/bin/bash
workspace=$1
map_fastq=$workspace/map_fastq
dnadata=$map_fastq/DNA_fastq
mkdir $workspace/MappedResult
mkdir $workspace/MappedResult/bwa_bam
mkdir $workspace/MappedResult/bwa_bed
bam_out=$workspace/MappedResult/bwa_bam
bed_out=$workspace/MappedResult/bwa_bed
reference=~/reference/human_DNA/liheng_Genome_reference/liheng_human_hg38_bwa
## 2. BWA reverse mapping: treat illumia read 1 as read 2; read 2 as read 1
cd $dnadata
for f in `ls *_dfp_R1.fastq.gz`
do
  echo "dna bwa_to_human Mapping ${f%_dfp_R1.fastq.gz}"
  bwa mem -t 24 -M $reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna ${f%_dfp_R1.fastq.gz}_dfp_R2.fastq.gz $f > $bam_out/temp.sam
  samtools view -@ 24 -S -b -q 40 $bam_out/temp.sam > $bam_out/temp.bam
  samtools sort $bam_out/temp.bam -o $bam_out/${f%_dfp_R1.fastq.gz}_dna.bam -@ 24
  bamToBed -i $bam_out/${f%_dfp_R1.fastq.gz}_dna.bam > $bed_out/${f%_dfp_R1.fastq.gz}.bed
  pigz $bed_out/${f%_dfp_R1.fastq.gz}.bed
  rm $bam_out/temp*
done
