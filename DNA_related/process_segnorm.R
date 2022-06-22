## Segment SegNorm data from Ginkgo
## dir: the full path to "SegNorm" file
## samples: Selected samples you want to perfomr segmentation
process_segnorm <- function(dir, samples = c()){
  sc_cnv <- read.delim(dir)
  sc_cnv_anno <- select(sc_cnv,CHR,START,END)
  sc_cnv_anno$feature <- paste0("seg_",seq(1:length(rownames(sc_cnv))))
  sc_cnv_anno <- sc_cnv_anno[sc_cnv_anno$CHR !="chrY",]
  allchr <-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")
  sc_cnv <- sc_cnv %>% filter(CHR !="chrY") %>% select(-CHR,-START,-END)
  rownames(sc_cnv) <- sc_cnv_anno$feature
  sc_cnv <- sc_cnv[,samples]
  col_demode <- 1 / apply(sc_cnv,2,densMode)
  sc_cnv <- sweep(sc_cnv,2,col_demode,"*")
  sc_cnv_seg <- mergeLevels_multi_segmentation(smoothed_cnv = sc_cnv, postionanno = sc_cnv_anno,gamma.param = 40,hg38 = hg38)
  sc_cnv <- data.frame(CHR = sc_cnv_anno$CHR, START = sc_cnv_anno$START, END = sc_cnv_anno$END,
                       feature = sc_cnv_anno$feature, round(sc_cnv_seg$seg.ml.data * 2)) %>%
    makeGRangesFromDataFrame(ignore.strand = T, keep.extra.columns = T)
  out <- list(sc_cnv_seg, sc_cnv)
  return(out)
}
