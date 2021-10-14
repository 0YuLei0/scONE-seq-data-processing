#!/bin/bash
workspace=$1
mkdir $workspace/map_fastq
mkdir $workspace/catfiles
map_fastq=$workspace/map_fastq
mkdir $map_fastq/DNA_fastq
mkdir $map_fastq/RNA_fastq
mkdir $map_fastq/Un_fastq
mkdir $map_fastq/fastp_report
report=$map_fastq/fastp_report
cd $workspace/catfiles
for i in `ls *read1.fastq.gz`
do
  sh ~/software/DR_sep/version2.2/allinone_sep.sh $i ${i%_read1.fastq.gz}_read2.fastq.gz
  mv dna_fp_r1.fastq.gz $map_fastq/DNA_fastq/${i%_read1.fastq.gz}_dfp_R1.fastq.gz
  mv dna_fp_r2.fastq.gz $map_fastq/DNA_fastq/${i%_read1.fastq.gz}_dfp_R2.fastq.gz
  mv rna_fp_r1.fastq.gz $map_fastq/RNA_fastq/${i%_read1.fastq.gz}_rfp_R1.fastq.gz
  mv rna_fp_r2.fastq.gz $map_fastq/RNA_fastq/${i%_read1.fastq.gz}_rfp_R2.fastq.gz
  mv unmatch_fp_r1.fastq.gz $map_fastq/Un_fastq/${i%_read1.fastq.gz}_ufp_R1.fastq.gz
  mv unmatch_fp_r2.fastq.gz $map_fastq/Un_fastq/${i%_read1.fastq.gz}_ufp_R2.fastq.gz
  mkdir $report/${i%_read1.fastq.gz}
  mv *html *json $report/${i%_read1.fastq.gz}/
done
## Remove Ecoli contamination in rna reads
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
