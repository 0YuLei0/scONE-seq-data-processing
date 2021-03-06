# scONE-seq-data-processing
This repository contains the necessary scripts to process sequencing data from single-cell RNA and DNA libraries with scONE-seq. To use the scripts in this repository, please create a conda environment and install all necessary packages in that environment. 

```
conda create -n oneseq
conda activate oneseq
conda install fastp=0.20.1 # Necessary packages are listed in Required_packages.txt
```

**RawDataProcessing:** this folder includes the necessary scripts to QC and align/quantify single-cell DNA and RNA data for downstream analysis.
1. We QC fastq files and separate DNA/RNA reads. (Use 1_process_reads.sh and Separation_bash.sh)
2. We next performed DNA UMI deduplication and alignment. (Use 2_DNA_BWA_mapping.sh and 3_DNA_UMI_mapping.sh)
3. We next performed RNA UMI deduplication and quantification. (Use 4_RNA_UMI_quantification.sh)

Briefly, once you have Illumina/MGI sequencing data, you could simply make a directory, for example, let's name it **scONEseq**. Under this directory, you should create a new directory named **Raw**. Then, you can put your fastq files in this **Raw** folder. Before you run the **RawDataProcessing** scripts, please remember to change your Reference, Barcodes, and scripts **PATH** accordingly. For the details about creating references used in the scONE-seq analysis, please check **DNA_related** and **RNA_related** sections.

```
## Directly run the scONE-seq pipline.
nohup sh scONE_pipline.sh $PATH/scONEseq &
## Run different sections seperately
nohup sh 1_process_reads.sh $PATH/scONEseq &
nohup sh 2_DNA_BWA_mapping.sh $PATH/scONEseq &
nohup sh 3_DNA_UMI_mapping.sh $PATH/scONEseq &
nohup sh 4_RNA_UMI_quantification.sh $PATH/scONEseq &
```

**Barcodes:** this folder includes DNA/RNA barcodes used in scONE-seq.

**DNA_related:** this folder includes the instruction to calculate integer CNVs from scONE-seq DNA data.

**RNA_related:** this folder includes the instruction to get preRNA(intron+exon) expression data from single nuclei RNA data.
