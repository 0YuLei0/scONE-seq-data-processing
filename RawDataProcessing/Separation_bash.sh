#!/bin/bash
read1=$1
read2=$2
## 1. use fastp to filter low quality pairs
fastp --thread=24 --disable_adapter_trimming --trim_tail1=1 --trim_poly_g --trim_poly_x --low_complexity_filter -i $read1 -o temp_fp_r1.fastq.gz -I $read2 -O temp_fp_r2.fastq.gz --json=${read1%_read1*}.json --html=${read1%_read1*}.html --report_title="${read1%_read1*}"
## 2. For MGI reads, to avoid later problem in DNA bwa mapping, we should rename the fastq head with Illumina format
zcat temp_fp_r1.fastq.gz |awk '{if(NR%4==1) gsub(/\/1/," 1:N:0"); print $0; }' | pigz > fp_r1.fastq.gz
zcat temp_fp_r2.fastq.gz |awk '{if(NR%4==1) gsub(/\/2/," 2:N:0"); print $0; }' | pigz > fp_r2.fastq.gz
## 3. Use cutadapt to demultiplex DNA/RNA reads, PLEASE REMEMBER TO CHANGE THE PATH TO scONEseq_barcode.fasta
cutadapt --quiet -e 0.4 -g file:$PATH/scONEseq_barcode.fasta -o trimmed-{name}.R2.fastq.gz -p trimmed-{name}.R1.fastq.gz fp_r2.fastq.gz fp_r1.fastq.gz
## 4. Extract RNA umi from fastq files with fastp
fastp --thread=24 --length_required=50 --disable_adapter_trimming --disable_quality_filtering -U --umi_loc=read1 --umi_len=8 --umi_prefix=umi -i trimmed-5RNA.R2.fastq.gz -o rna_fp_r2.fastq.gz -I trimmed-5RNA.R1.fastq.gz -O rna_fp_r1.fastq.gz --json="rna".json --html="rna".html --report_title="rna"
## 5. Extract DNA umi and remvoe Tn5 mosaic sequeces
fastp --thread=24 --length_required=50 --disable_adapter_trimming --disable_quality_filtering -U --umi_loc=read1 --umi_len=8 --umi_prefix=umi -i trimmed-DNA.R2.fastq.gz -o temp -I trimmed-DNA.R1.fastq.gz -O dna_fp_r1.fastq.gz --json="rna".json --html="dna".html --report_title="dna"
seqtk trimfq -b 20 temp | pigz > dna_fp_r2.fastq.gz
## remove temporory files
rm temp
