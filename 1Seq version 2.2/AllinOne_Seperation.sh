#!/bin/bash
read1=$1
read2=$2
barcode=~/software/DR_sep/version2.2
# use fastp to filter low quality pairs
fastp --thread=16 --disable_adapter_trimming --trim_tail1=1 --trim_poly_g --trim_poly_x --low_complexity_filter -i $read1 -o fp_r1.fastq.gz -I $read2 -O fp_r2.fastq.gz --json=${read1%_read1*}.json --html=${read1%_read1*}.html --report_title="${read1%_read1*}"
# Use seqkit to grep dna reads
seqkit grep -s -r -i -f $barcode/dnabarcode_mutation fp_r2.fastq.gz -j 2 -o temp_dna_mut_read2.fastq.gz
seqkit grep -s -r -i -f $barcode/dnabarcode_mutation fp_r2.fastq.gz -j 2 -v -o temp_remain0.fastq.gz
seqkit grep -s -r -i -f $barcode/dnabarcode_insertion temp_remain0.fastq.gz -j 2 -o temp_dna_ins_read2.fastq.gz
seqkit grep -s -r -i -f $barcode/dnabarcode_insertion temp_remain0.fastq.gz -j 2 -v -o temp_remain1.fastq.gz
seqkit grep -s -r -i -f $barcode/dnabarcode_deletion temp_remain1.fastq.gz -j 2 -o temp_dna_del_read2.fastq.gz
seqkit grep -s -r -i -f $barcode/dnabarcode_deletion temp_remain1.fastq.gz -j 2 -v -o not_dna.fastq.gz
# cut the dna barcode with seqtk trimfq
seqtk trimfq -b 5 temp_dna_mut_read2.fastq.gz > dna_read2.fastq
seqtk trimfq -b 6 temp_dna_ins_read2.fastq.gz >> dna_read2.fastq
seqtk trimfq -b 4 temp_dna_del_read2.fastq.gz >> dna_read2.fastq
pigz dna_read2.fastq
rm temp*
# Optional step: incease running time
# Confirm with Mosaic sequencing in DNA reads; AGATGTGTATAAGAGACAG
# mv dna_read2.fastq before_confirm_dna_read2.fastq
# cat before_confirm_dna_read2.fastq | seqkit grep -s -i -p AGATGTGTATAAGAGACAG -m 4 -R 1:27 > dna_read2.fastq
# Use seqkit to grep rna reads
seqkit grep -s -r -i -f $barcode/rnabarcode_mutation not_dna.fastq.gz -j 2 -o temp_rna_mut_read2.fastq.gz
seqkit grep -s -r -i -f $barcode/rnabarcode_mutation not_dna.fastq.gz -j 2 -v -o temp_remain0.fastq.gz
seqkit grep -s -r -i -f $barcode/rnabarcode_insertion temp_remain0.fastq.gz -j 2 -o temp_rna_ins_read2.fastq.gz
seqkit grep -s -r -i -f $barcode/rnabarcode_insertion temp_remain0.fastq.gz -j 2 -v -o temp_remain1.fastq.gz
seqkit grep -s -r -i -f $barcode/rnabarcode_deletion temp_remain1.fastq.gz -j 2 -o temp_rna_del_read2.fastq.gz
seqkit grep -s -r -i -f $barcode/rnabarcode_deletion temp_remain1.fastq.gz -j 2 -v -o unmatch.fastq.gz
# cut the dna barcode with seqtk trimfq
seqtk trimfq -b 5 temp_rna_mut_read2.fastq.gz > rna_read2.fastq
seqtk trimfq -b 6 temp_rna_ins_read2.fastq.gz >> rna_read2.fastq
seqtk trimfq -b 4 temp_rna_del_read2.fastq.gz >> rna_read2.fastq
pigz rna_read2.fastq
rm temp* not_dna.fastq.gz before_confirm_dna_read2.fastq
# Read1 seperation with pair repair algorithm
## Please note that the algorithm can tell the read1 read2; in this case, out1 is real read1!!!
repair.sh in1=dna_read2.fastq.gz in2=fp_r1.fastq.gz out1=dna_pair_read1.fastq.gz out2=dna_pair_read2.fastq.gz outsingle=temp_not_dna_read1.fastq.gz
repair.sh in1=rna_read2.fastq.gz in2=temp_not_dna_read1.fastq.gz out1=rna_pair_read1.fastq.gz out2=rna_pair_read2.fastq.gz outsingle=temp_unmatch_read1.fastq.gz
repair.sh in1=unmatch.fastq.gz in2=temp_unmatch_read1.fastq.gz out1=unm_pair_read1.fastq.gz out2=unm_pair_read2.fastq.gz outsingle=temp
rm temp* dna_read2.fastq.gz rna_read2.fastq.gz unmatch.fastq.gz
# extract umi
fastp --thread=16 --length_required=30 --disable_adapter_trimming --disable_quality_filtering -U --umi_loc=read1 --umi_len=5 --umi_prefix=umi -i rna_pair_read2.fastq.gz -o rna_fp_r2.fastq.gz -I rna_pair_read1.fastq.gz -O rna_fp_r1.fastq.gz --json="rna".json --html="rna".html --report_title="rna"
# dna umi and Mosaic sequence triming
fastp --thread=16 --length_required=30 --disable_adapter_trimming --disable_quality_filtering -U --umi_loc=read1 --umi_len=5 --umi_prefix=umi -i dna_pair_read2.fastq.gz -o temp_dna_fp_r2.fastq.gz -I dna_pair_read1.fastq.gz -O dna_fp_r1.fastq.gz --json="dna".json --html="dna".html --report_title="dna"
seqtk trimfq -b 20 temp_dna_fp_r2.fastq.gz > dna_fp_r2.fastq
pigz dna_fp_r2.fastq
rm temp_dna_fp_r2.fastq.gz
# process unmatch reads with strict criterion:
fastp --thread=16 --length_required=35 --disable_adapter_trimming --trim_poly_g --trim_poly_x --low_complexity_filter --average_qual=20 -i unm_pair_read1.fastq.gz -o unmatch_fp_r1.fastq.gz -I unm_pair_read2.fastq.gz -O unmatch_fp_r2.fastq.gz --json="unmatch".json --html="unmatch".html --report_title="unmatch"
# remove redundant files
rm fp_r[1,2].fastq.gz *pair*fastq.gz
