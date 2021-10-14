#!/bin/bash
input=$1
sh ~/software/DR_sep/version2.2/0_cat.sh $input
sh ~/software/DR_sep/version2.2/1_process.sh $input
sh ~/software/DR_sep/version2.2/2_rnaQuant.sh $input
sh ~/software/DR_sep/version2.2/3_dnaMap.sh $input
sh ~/software/DR_sep/version2.2/4_umiQuant.sh $input
sh ~/software/DR_sep/version2.2/5_quanlimap.sh $input
sh ~/software/DR_sep/version2.2/3.1_dnaUMI.sh $input

#		~/GBM/1st_batch_gbm 392
1169
#		~/MixSpecies/GBM_Tn5 288
14745
#   ~/GBM/2nd_batch_gbm_1 288
31324
#   ~/GBM/2nd_batch_gbm2
nohup sh ~/software/DR_sep/version2.2/sconeseq.sh ~/GBM/2nd_batch_gbm2 &
39474
# ~/MixSpecies/MixTest20210323
nohup sh ~/software/DR_sep/version2.2/sconeseq.sh ~/MixSpecies/MixTest20210323 &
38882
#   ~/GBM/2nd_batch_gbm_13 288

nohup sh ~/software/DR_sep/version2.2/3.1_dnaUMI.sh ~/GBM/1st_batch_gbm &
nohup sh ~/software/DR_sep/version2.2/3.1_dnaUMI.sh ~/MixSpecies/2nd_batch_gbm0 &

nohup sh ~/software/DR_sep/version2.2/3.1_dnaUMI.sh ~/GBM/2nd_batch_gbm1 &
5388
nohup sh ~/software/DR_sep/version2.2/3.1_dnaUMI.sh ~/GBM/2nd_batch_gbm2 &
4280
nohup sh ~/software/DR_sep/version2.2/3.1_dnaUMI.sh ~/GBM/2nd_batch_gbm3 &
4358
