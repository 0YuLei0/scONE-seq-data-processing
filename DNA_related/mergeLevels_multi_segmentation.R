## R function to perform segmentation with"copynumber" and merge segmentation with "aCGH"
cytoband_hg38 <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz"
cb <- data.table::fread(cytoband_hg38, header = FALSE, data.table = FALSE)
cb <- cb[cb$V1 %in% c(paste0("chr", 1:22), "chrX", "chrY"), ]
rownames(cb) <- seq(1:nrow(cb))
cb$V1 <- factor(cb$V1)
cb$V4 <- factor(cb$V4)
cb$V5 <- factor(cb$V5)
hg38 <- cb
rm(cytoband_hg38,cb)
mergeLevels_multi_segmentation <- function(smoothed_cnv = smooth_cnv,
                                           postionanno = bin_anno,
                                           gamma.param = 40, hg38 = hg38) {
  ## all chr
  allchr_2 <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X")
  ##
  log.ratio.mats <- smoothed_cnv
  log.ratio.mats[log.ratio.mats <= 0] <- 0.0000001
  log.ratio.mats <- apply(log.ratio.mats, 2, function(x){x / median(x)})
  log.ratio.mats <- log2(log.ratio.mats)
  postionanno <- postionanno %>% dplyr::select(CHR, START)
  postionanno$CHR <- gsub('[chr]', '', postionanno$CHR)
  postionanno$CHR  <- factor(postionanno$CHR, levels = allchr_2)
  sample.matrices <- data.frame(cbind(postionanno, log.ratio.mats))
  sample.matrices <- copynumber::winsorize(sample.matrices,assembly = "hg38")
  segmentations <- multipcf(sample.matrices,gamma = gamma.param, assembly = "hg38", return.est=TRUE, fast = FALSE)
  log.ratio.longs <- segmentations$estimates[,3:ncol(segmentations$estimates)]
  ratio.longs <- 2^log.ratio.longs
  rownames(ratio.longs) <- postionanno$feature
  seg.ml.data <- data.frame(row.names = rownames(smoothed_cnv),ratio.longs)
  for(i in colnames(smoothed_cnv)){
    logratio <- log.ratio.mats[,i]
    logratio <- 2^logratio
    seg.mean <- ratio.longs[,i]
    ml <- aCGH::mergeLevels(logratio,seg.mean)
    seg.ml.data[,i] <- ml$vecMerged
  }
  seg.data <- ratio.longs
  output <- list(seg.data = seg.data, seg.ml.data = seg.ml.data)
  return(output)
}
library(preprocessCore)
densMode <- function(x){
  td <- density(x,bw = 0.001)
  maxDens <- which.max(td$y)
  td$x[maxDens]
}
