#!/bin/bash
## The path of your working content
input=$1
## 
sh $PATH/1_process_reads.sh $input
sh $PATH/2_DNA_BWA_mapping.sh $input
sh $PATH/3_DNA_UMI_mapping.sh $input
sh $PATH/4_RNA_UMI_quantification.sh $input
