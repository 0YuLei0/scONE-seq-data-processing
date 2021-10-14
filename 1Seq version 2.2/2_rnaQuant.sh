#!/bin/bash
workspace=$1
map_fastq=$workspace/map_fastq
rnadata=$map_fastq/Clean_rna
mkdir $workspace/MappedResult
mkdir $workspace/MappedResult/kalli_count
mkdir $workspace/MappedResult/salmon_count
kmap_out=$workspace/MappedResult/kalli_count
sal_out=$workspace/MappedResult/salmon_count
kmap_ref=~/reference/humanRNA/kallisto_ref
sal_intron_ref=~/reference/humanRNA/salmon_ref/gv34_prerna_s14
## Kallisto Transcriptome Mapping with Gencode Reference
cd $rnadata
for a in `ls *_rfp_R1.fastq.gz`
do
  echo "quant ${a%_rfp_R1.fastq.gz}"
  mkdir $kmap_out/${a%_rfp_R1.fastq.gz}
  kallisto quant -t 24 --bias -i $kmap_ref/gcv29_rRNA_ebv_ercc.idx -o $kmap_out/${a%_rfp_R1.fastq.gz}/ $a ${a%_rfp_R1.fastq.gz}_rfp_R2.fastq.gz
  ~/software/salmon1.4/bin/salmon quant --seqBias --gcBias -p 24 -i $sal_intron_ref -l A -1 $a -2 ${a%_rfp_R1.fastq.gz}_rfp_R2.fastq.gz --validateMappings -o $sal_out/${a%_rfp_R1.fastq.gz}
done

for a in `ls *_read1.fastq.gz`
do
  echo "quant ${a%_read1.fastq.gz}"
  kallisto quant -t 24 --bias -i ~/reference/humanRNA/kallisto_ref/gcv29_rRNA_ebv_ercc.idx -o ./kquant/${a%_read1.fastq.gz}/ $a ${a%_read1.fastq.gz}_read2.fastq.gz
done
