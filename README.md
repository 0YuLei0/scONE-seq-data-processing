# scONE-seq-data-processing
This repository contains the neccessary scripts to process single-cell RNA and DNA libraries with scONE-seq.

**RawDataProcessing:** this folder includes the neccessary scripts to convet single cell DNA and RNA data to matrix-format data for down-stream analysis. In this section:

1. We QC fastq files and separate DNA/RNA reads.
2. We next performed UMI deduplication. 
3. The deduplicated reads could be mapped/quantified. 

**Barcodes:** this folder includes the cell indexs and DNA/RNA barcode used in scONE-seq.

**DNA_related:** this folder includes the **DNA** data analysis code.

**RNA_related:** this folder includes the **RNA** data analysis code.
