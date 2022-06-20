#!/bin/bash
workspace=$1
map_fastq=$workspace/map_fastq
rnadata=$map_fastq/RNA_fastq
mkdir $workspace/umidedup
umi=$workspace/umidedup
mkdir $umi/umi_mapped
mkdir $umi/umi_mapped/Star_reverse
rnastar=$umi/umi_mapped/Star_reverse
star_from10x=$PATH_to_your_STAR_reference
# dedup part
mkdir $umi/dedup_fastq
mkdir $umi/dedup_fastq/dedup_rna
deduprna=$umi/dedup_fastq/dedup_rna
# umi kquant part
mkdir $umi/umi_mapped/umi_kquant
umikout=$umi/umi_mapped/umi_kquant
kquant_ref=$PATH_to_your_kallisto_reference
# umi salmon intron+exon quant part
mkdir $umi/umi_mapped/umi_salmonIntron
intron_out=$umi/umi_mapped/umi_salmonIntron
intron_ref=$PATH_to_your_salmonIntron_reference
## 1. Star reverse mapping: treat read 1 as read 2; read 2 as read 1
cd $rnadata
for f in `ls *_rfp_R1.fastq.gz`
do
  echo "$f RNA Reads STAR Mapping"
  mkdir $rnastar/${f%_rfp_R1.fastq.gz}
  STAR --runThreadN 16 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --genomeDir $star_from10x --readFilesCommand zcat --readFilesIn ${f%_rfp_R1.fastq.gz}_rfp_R2.fastq.gz $f --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $rnastar/${f%_rfp_R1.fastq.gz}/
done
## 2. umi deduplication and reads extraction
cd $rnastar
for i in `ls`
do
  echo "umitools dedup sample $i"
  samtools view -hbf 64 $i/Aligned.sortedByCoord.out.bam > R2.bam -@ 24
  samtools index R2.bam -@ 24
  umi_tools dedup -I R2.bam -S dedup.bam --umi-separator=":umi_"
  bedtools bamtofastq -i dedup.bam -fq single.fq
  rm paired1.fq paired2.fq dedup.bam
  pigz single.fq
  repair.sh in1=$rnadata/$i\_rfp_R1.fastq.gz in2=single.fq.gz out1=$deduprna/$i\_dedup_R1.fastq.gz out2=$deduprna/$i\_dedup_R2.fastq.gz outsingle=single_rna
  rm single_rna R2.bam* single.fq*
done
## 3 Kallisto quantification exon regions
cd $deduprna
for a in `ls *_dedup_R1.fastq.gz`
do
  echo "Kallisto_Mapping ${a%_dedup_R1.fastq.gz}"
  mkdir $umikout/${a%_dedup_R1.fastq.gz}
  kallisto quant -t 24 --bias -i $kquant_ref -o $umikout/${a%_dedup_R1.fastq.gz}/ $a ${a%_dedup_R1.fastq.gz}_dedup_R2.fastq.gz
done
## 5 Salmon quant with intron+exon reference
cd $deduprna
for a in `ls *_dedup_R1.fastq.gz`
do
  echo "Salmon_intron_Mapping ${a%_dedup_R1.fastq.gz}"
  ~/software/salmon1.8/bin/salmon quant --seqBias --gcBias -p 24 -i $intron_ref -l A -1 $a -2 ${a%_dedup_R1.fastq.gz}_dedup_R2.fastq.gz -o $intron_out/${a%_dedup_R1.fastq.gz}
done
