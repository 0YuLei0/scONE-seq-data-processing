##########################################################
#!/bin/bash
## 0. working content and file path
# reverse mapping part
workspace=$1
map_fastq=$workspace/map_fastq
rnadata=$map_fastq/Clean_rna
mkdir $workspace/umidedup
umi=$workspace/umidedup
mkdir $umi/umi_mapped
mkdir $umi/umi_mapped/Star_reverse
rnastar=$umi/umi_mapped/Star_reverse
star_from10x=~/reference/humanRNA/star_from10x_ercc
# dedup part
mkdir $umi/dedup_fastq
mkdir $umi/dedup_fastq/dedup_rna
deduprna=$umi/dedup_fastq/dedup_rna
# umi kquant part
mkdir $umi/umi_mapped/umi_kquant
umikout=$umi/umi_mapped/umi_kquant
kquant_ref=~/reference/humanRNA/kallisto_ref
# umi salmon intron+exon quant part
mkdir $umi/umi_mapped/umi_salmonIntron
intron_out=$umi/umi_mapped/umi_salmonIntron
intron_ref=~/reference/humanRNA/salmon_ref/gv34_prerna_s14
## 2. Star reverse mapping: treat illumia read 1 as read 2; read 2 as read 1
cd $rnadata
unpigz *gz
for f in `ls *_rfp_R1.fastq`
do
  echo "$f RNA Reads STAR Mapping"
  mkdir $rnastar/${f%_rfp_R1.fastq}
  STAR --runThreadN 16 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --genomeDir $star_from10x --readFilesIn ${f%_rfp_R1.fastq}_rfp_R2.fastq $f --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $rnastar/${f%_rfp_R1.fastq}/
done
## 3. umi deduplication and reads extraction
cd $rnastar
for i in `ls`
do
  echo "umitools dedup sample $i"
  samtools view -hbf 64 $i/Aligned.sortedByCoord.out.bam > R2.bam -@ 24
  samtools index R2.bam -@ 24
  umi_tools dedup -I R2.bam -S dedup.bam --umi-separator=":umi_"
  # umi_tools dedup -I R2.bam  -S dedup.bam --umi-separator=":ATCGT"
  # samtools fastq -1 paired1.fq -2 paired2.fq -n dedup.bam -s single.fq
  bedtools bamtofastq -i dedup.bam -fq single.fq
  rm paired1.fq paired2.fq dedup.bam
  # repair.sh in1=S961-Tn025-DR01_rfp_R1.fastq in2=aln.end1.fq out1=dedup_R1.fastq out2=dedup_R2.fastq outsingle=single_rna
  repair.sh in1=$rnadata/$i\_rfp_R1.fastq in2=single.fq out1=$deduprna/$i\_dedup_R1.fastq out2=$deduprna/$i\_dedup_R2.fastq outsingle=single_rna
  rm single_rna R2.bam* single.fq
done
## 4. gzip to save space
cd $rnadata
pigz *fastq
cd $deduprna
pigz *fastq
## 5. RNA umi quantification
## 5.1 Kallisto
cd $deduprna
for a in `ls *_dedup_R1.fastq.gz`
do
  echo "Kallisto_Mapping ${a%_dedup_R1.fastq.gz}"
  mkdir $umikout/${a%_dedup_R1.fastq.gz}
  kallisto quant -t 24 --bias -i $kquant_ref/gcv29_rRNA_ebv_ercc.idx -o $umikout/${a%_dedup_R1.fastq.gz}/ $a ${a%_dedup_R1.fastq.gz}_dedup_R2.fastq.gz
done
## 5.2 Salmon quant with intron+exon reference
cd $deduprna
for a in `ls *_dedup_R1.fastq.gz`
do
  echo "Salmon_intron_Mapping ${a%_dedup_R1.fastq.gz}"
  ~/software/salmon1.4/bin/salmon quant --seqBias --gcBias -p 24 -i $intron_ref -l A -1 $a -2 ${a%_dedup_R1.fastq.gz}_dedup_R2.fastq.gz --validateMappings -o $intron_out/${a%_dedup_R1.fastq.gz}
done
