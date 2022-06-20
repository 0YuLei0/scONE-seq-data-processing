# scONE-seq-data-processing
This repository contains the neccessary scripts to process single-cell RNA and DNA libraries with scONE-seq.

**RawDataProcessing:** this folder includes the neccessary scripts to process single cell DNA and RNA data to matrix-format for downstream analysis. In this section:

1. We QC fastq files and separate DNA/RNA reads. (Use 1_process_reads.sh and Separation_bash.sh)
2. We next performed DNA UMI deduplication and alignment. (Use 2_DNA_BWA_mapping.sh and 3_DNA_UMI_mapping.sh)
3. We next performed RNA UMI deduplication and quantification. (Use 4_RNA_UMI_quantification.sh)

**Barcodes:** this folder includes the cell indexs and DNA/RNA barcode used in scONE-seq.

**DNA_related:** this folder includes the **DNA** data analysis code.

**RNA_related:** this folder includes the **RNA** data analysis code.
