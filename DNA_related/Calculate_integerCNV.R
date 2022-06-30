## correct integer CNVs by considering the BAF information
baf_cnv_ginkgo <- function(rdr, mbaf, fixed_rdr, binclusters, contral_baf, contral_rdr){

  if(all(colnames(rdr) != colnames(mbaf))){
    message("sample name doesn't match1")
  }
  if(all(colnames(rdr) != colnames(fixed_rdr))){
    message("sample name doesn't match2")
  }
  if(all(colnames(rdr) != colnames(binclusters))){
    message("sample name doesn't match3")
  }
  if(all(rownames(rdr) != rownames(mbaf))){
    message("feature name doesn't match1")
  }
  if(all(rownames(rdr) != rownames(fixed_rdr))){
    message("feature name doesn't match2")
  }
  if(all(rownames(rdr) != rownames(binclusters))){
    message("feature name doesn't match3")
  }
  ## correct all conditions
  return_corrected <- function(condition,merge_region, inputcnv,
                               states,hyper1, hyper2,deep, neu, deci3copy) {
    if(condition == 0){
      return(inputcnv)
    }else if(merge_region == "merging"){
      return(4)
    }else{
      if(states == "del_3" & deep <= 0.1){return(3)}
      else if(states == "del_3" & deep > 0.1){return(inputcnv)}
      else if(states == "del_2" & deep < 0.001){return(4)}
      else if(states == "del_2" & deep > 0.001){return(2)}
      else if(states == "netrual_2"){return(4)}
      else if(states == "netrual_3" & neu >= deci3copy){return(4)}
      else if(states == "netrual_3" & neu < deci3copy){return(3)}
      else if(states == "amp_3"){return(6)}
      else if(states == "amp_2" & hyper2 < 0.05){return(8)}
      #else if(states == "amp_2" & hyper2 < 0.05){return(inputcnv)}
      else if(states == "amp_2" & hyper1 > 0.001){return(inputcnv)}
      else if(states == "amp_2" & hyper1 <= 0.001){return(4)}
    }
  }
  samples <- colnames(rdr)
  features <- rownames(rdr)
  infered_cnv <- matrix(0, nrow = nrow(rdr), ncol = ncol(rdr))
  rownames(infered_cnv) <- features
  colnames(infered_cnv) <- samples

  bin_states <- matrix("a", nrow = nrow(rdr), ncol = ncol(rdr))
  rownames(bin_states) <- features
  colnames(bin_states) <- samples

  merging_states <- matrix("a", nrow = nrow(rdr), ncol = ncol(rdr))
  rownames(merging_states) <- features
  colnames(merging_states) <- samples

  corrected_cnv <- matrix(0, nrow = nrow(rdr), ncol = ncol(rdr))
  rownames(corrected_cnv) <- features
  colnames(corrected_cnv) <- samples

  sample_ploidy_final <- c()
  singlets_states <- c()
  
  ## Transfer Ginkgo segmentation to cluster based annotation
  seg2cluster <- function(segdata){
  samples <- colnames(segdata)
  features <- rownames(segdata)

  outdata <- matrix(0, nrow = length(features), ncol = length(samples))
  rownames(outdata) <- features
  colnames(outdata) <- samples

  for (i in seq(1: length(samples))) {
    mysamples <- samples[i]
    a <- segdata[,mysamples]
    Fs <- 0
    Cs <- 1
    adata <- data.frame(cluster = 1:1000, freq = rep(0,100))
    a <- c(a, 100000000)
    for (ii in seq(1:(length(a)-1))) {
      Fs <- Fs + 1
      if(a[ii] != a[ii+1]){
        if(Fs > 10){
          adata$freq[Cs] <- Fs
          Fs <- 0
          Cs <- Cs+1
        }else if(ii == (length(a)-1)){
          adata$freq[Cs] <- Fs
        }
      }
    }
    adata <- adata[adata$freq > 0,]
    cluster_mat <- rep(adata$cluster, adata$freq)
    outdata[,mysamples] <- cluster_mat
    }
  return(outdata)
  }
  binclusters <- seg2cluster(segdata = fixed_rdr)

  for (i in seq(1: length(samples))) {
    mysample <- samples[i]
    message("Processing Sample ", samples[i], "; processed ", i-1, " samples")
    sample_data <- data.frame(feature = features,
                              type = binclusters[,mysample],
                              rdrvalues = rdr[,mysample],
                              #logrdr = logrdr[,mysample],
                              fixvalues = fixed_rdr[,mysample],
                              bafvalues = mbaf[,mysample]
    )
    #my_rdrmedian <- median(sample_data$rdrvalues)
    perform_test <- sample_data %>% dplyr::group_by(type) %>%
      dplyr::summarise(baf_mean = mean(bafvalues, na.rm = TRUE), baf_median = median(bafvalues, na.rm = TRUE),
                       rdr_mean = mean(rdrvalues, na.rm = TRUE), rdr_median = median(rdrvalues, na.rm = TRUE),
                       fixvalues = mean(fixvalues, na.rm = TRUE),
                       baf_gini = DescTools::Gini(unique(bafvalues)), baf_cv = (sd(unique(bafvalues)) / mean(unique(bafvalues))),
                       rdr_gini = DescTools::Gini(rdrvalues), rdr_cv = (sd(rdrvalues) / mean(rdrvalues)),
                       decision_loh = as.numeric(wilcox.test(unique(bafvalues), mu = 0.35, correct = FALSE, alternative = "less")[3]),
                       decision_3copy = as.numeric(wilcox.test(unique(bafvalues), mu = 1/6, correct = FALSE)[3]),
                       decision_2copy = as.numeric(wilcox.test(unique(bafvalues), contral_baf, alternative = "greater")[3]),
                       decision_2copy_2 = as.numeric(wilcox.test(unique(bafvalues), mu =0.06140351, alternative = "greater")[3]),
                       decision_del = as.numeric(wilcox.test(rdrvalues, contral_rdr, alternative = "less")[3]),
                       decision_3del = as.numeric(wilcox.test(rdrvalues, 3*contral_rdr/4)[3]),
                       decision_deep = as.numeric(wilcox.test(rdrvalues, contral_rdr/2, alternative = "greater")[3]),
                       decision_amp = as.numeric(wilcox.test(rdrvalues, contral_rdr, alternative = "greater")[3]),
                       decision_hyper1 = as.numeric(wilcox.test(rdrvalues, 3*contral_rdr/2, alternative = "less")[3]),
                       decision_hyper2 = as.numeric(wilcox.test(rdrvalues, 3*contral_rdr/2, alternative = "greater")[3]),
                       #decision_neu1 = as.numeric(wilcox.test(rdrvalues, mu = my_rdrmedian, correct = FALSE)[3]),
                       decision_neu = as.numeric(wilcox.test(rdrvalues, contral_rdr)[3])) %>%
      dplyr::mutate(decision_amp_pj = p.adjust(decision_amp, method = "bonferroni")) %>%
      dplyr::mutate(decision_neu_pj = p.adjust(decision_neu, method = "bonferroni")) %>%
      dplyr::mutate(copy1 = ifelse(test = (decision_loh > 0.1 & decision_loh > decision_3copy),yes = 1,no = 10)) %>%
      dplyr::mutate(copy3_1 = ifelse(test = (decision_neu_pj < 0.001 & decision_loh < 0.1 & decision_3del > 0.1 & decision_2copy < 0.2 & baf_gini < 0.45), yes = 3,no = 10)) %>%
      dplyr::mutate(copy3_2 = ifelse(test = (decision_neu_pj < 0.001 & decision_loh < 0.1 & decision_3copy > 0.1 & decision_2copy < 0.1 & baf_gini < 0.45), yes = 3,no = 10)) %>%
      dplyr::mutate(hete_type = ifelse(test = (copy1 != 1 & copy3_1 != 3 & copy3_2 != 3), yes = 2,ifelse(copy1 == 1,yes = 1,no = 3))) %>%
      dplyr::mutate(copystates = ifelse(decision_amp_pj < 10^-8, yes = "amp", ifelse(decision_del < 0.001, yes = "del", "netrual"))) %>%
      dplyr::mutate(alltype = paste0(copystates,"_", hete_type))
    perform_test <- as.data.frame(perform_test)
    perform_test[is.na(perform_test)] <- 0
    perform_test <- perform_test %>%
      dplyr::mutate(merging = ifelse(test = (baf_gini < 0.5 & baf_cv < 1), yes = "good", no = "merging"))

    copy2_rdr <- median(perform_test$rdr_median[perform_test$alltype == "netrual_2"])
    perform_test$cnv <- round(4 * (perform_test$rdr_median / copy2_rdr))
    perform_test$divi <- perform_test$cnv %% perform_test$hete_type
    perform_test[is.na(perform_test)] <- 0

    perform_test_2 <- perform_test %>% dplyr::group_by(type) %>%
      dplyr::summarise(states = alltype, CNV = cnv, Merging = merging,
                       corrected = return_corrected(condition = divi, merge_region = merging, inputcnv = cnv,
                                                    states = alltype, hyper1 = decision_hyper1,hyper2= decision_hyper2,
                                                    deep = decision_deep, neu = decision_neu_pj, deci3copy = decision_3copy))

    sample_data_1 <- sample_data %>% dplyr::inner_join(perform_test,by = "type")
    sample_data_2 <- sample_data %>% dplyr::inner_join(perform_test_2,by = "type")
    infered_cnv[,mysample] <- sample_data_2$CNV
    bin_states[,mysample] <- sample_data_2$states
    merging_states[,mysample] <- sample_data_2$Merging
    corrected_cnv[,mysample] <- sample_data_2$corrected
    ploidy <- sum(sample_data_2$corrected) / length(sample_data_2$corrected)
    sample_ploidy_final <- c(sample_ploidy_final,ploidy)
    dl_states <- all(perform_test$copy1 == 10) & (("del" %in% perform_test$copystates) | ("amp" %in% perform_test$copystates))
    singlets_states <- c(singlets_states, dl_states)
  }
  sample_ploidy <- data.frame(SampleID = samples, ploidy = sample_ploidy_final, doulets = singlets_states)
  estimationCNV <- list(SampleInfo = sample_ploidy,
                        raw_cnv = infered_cnv,
                        bin_anno = bin_states,
                        merge_region = merging_states,
                        F_CNV = corrected_cnv)
  return(estimationCNV)
}
