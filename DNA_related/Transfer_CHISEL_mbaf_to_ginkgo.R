## Transfer CHISEL bins(5Mb) to Ginkgo bins(500kb)
chisel2ginkgo <- function(ginkgo_bins, chisel_bins, mbaf) {
  samples <- colnames(mbaf)
  used_CHR <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
  ginkgo_bins <- ginkgo_bins %>% filter(CHR %in% used_CHR)
  chisel_bins <- chisel_bins %>% filter(CHR %in% used_CHR)
  features <- ginkgo_bins$feature

  ginkgo_mabf <- matrix(0, nrow = length(features), ncol = length(samples))
  rownames(ginkgo_mabf) <- features
  colnames(ginkgo_mabf) <- samples

  for (i in seq(1: length(samples))) {
    mysample <- samples[i]

    ginkgo <- data.frame(CHR = ginkgo_bins$CHR, START = ginkgo_bins$START, END = ginkgo_bins$END,
                         feature = ginkgo_bins$feature) %>%
      makeGRangesFromDataFrame(ignore.strand = T, keep.extra.columns = T)
    chisel <- data.frame(CHR = chisel_bins$CHR, START = chisel_bins$START, END = chisel_bins$END,
                         mbaf = mbaf[,mysample]) %>%
      makeGRangesFromDataFrame(ignore.strand = T, keep.extra.columns = T)

    olaps1 <- findOverlaps(ginkgo, chisel)

    mk_df1 <- tibble(feature = ginkgo$feature[queryHits(olaps1)],
                     mbaf = chisel$mbaf[subjectHits(olaps1)]) %>%
      dplyr::distinct(feature, .keep_all = TRUE)

    ginkgo_mabf <- ginkgo_mabf[mk_df1$feature,]
    ginkgo_mabf[,mysample] <- mk_df1$mbaf
  }
  return(ginkgo_mabf)
}

