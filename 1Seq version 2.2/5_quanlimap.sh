#!/bin/bash
workspace=$1
umi=$workspace/umidedup
mkdir $umi/umi_mapped/Star_reverse/bamfiles
mkdir $umi/umi_mapped/Star_reverse/quanlimap
bamfiles=$umi/umi_mapped/Star_reverse/bamfiles
cd $umi/umi_mapped/Star_reverse/
for i in `ls`
do
  mv $i/Aligned.sortedByCoord.out.bam $bamfiles/$i\_starrna.bam
  # samtools index -@ 12 $i\_star.bam
done
cd $bamfiles
for f in `ls *_starrna.bam`
do
  samtools index -@ 12 $f
  mkdir ../quanlimap/${f%_starrna.bam}/
  qualimap rnaseq -bam $f -gtf ~/reference/humanRNA/refdata-gex-GRCh38-2020-A/genes/genes.gtf --java-mem-size=20G -outdir ../quanlimap/${f%_starrna.bam}/
done
bamfiles2=$workspace/MappedResult/bwa_bam
cd $bamfiles2
mkdir report
for i in `ls *dna.bam`
do
  echo "quanlimap sample ${i%_dna.bam}"
  mkdir ./report/${i%_dna.bam}
  qualimap bamqc -bam $i --java-mem-size=20G -outdir ./report/${i%_dna.bam}
done
