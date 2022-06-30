library(GenomicRanges)
library(IRanges)
library(tibble)
## Use this function to convert CHISEL output combo file to cell x genome_region matrix
extract_bafANDrdr <- function(combo, cellid, cellid_trim = 10) {
  cellid <- dplyr::select(cellid,X.CELL,BARCODE)
  colnames(cellid) <- c("SampleID", "BARCODE")
  data1 <- inner_join(combo, cellid)
  chisel_binanno <- filter(combo, BARCODE == cellid$BARCODE[1]) %>% select(CHR, START, END, N_reads)
  chisel_binanno$feature <- paste0("chisel_", seq(1:length(chisel_binanno$N_reads)))
  rdr.data <- data.frame(row.names = chisel_binanno$feature, a = rep(1,length(chisel_binanno$CHR)))
  baf.data <- data.frame(row.names = chisel_binanno$feature, a = rep(1,length(chisel_binanno$CHR)))
  mbaf.data <- data.frame(row.names = chisel_binanno$feature, a = rep(1,length(chisel_binanno$CHR)))
  for (i in c(1:length(cellid$SampleID))) {
    running_cell <- cellid$SampleID[i]
    rdr <- data1 %>% filter(SampleID == running_cell) %>% pull(RDR)
    baf <- data1 %>% filter(SampleID == running_cell) %>% pull(BAF)
    mbaf <- abs(0.5 - baf)
    rdr.data[,i] <- rdr
    baf.data[,i] <- baf
    mbaf.data[,i] <- mbaf
  }
  cellid_trim <- -cellid_trim
  SampleID <- cellid$SampleID
  stri_sub(SampleID,cellid_trim, -1) <- ""
  mysample <- SampleID[stri_count_fixed(SampleID, '-') == 3]
  stri_sub(mysample,-3,-3) <- ""
  SampleID[stri_count_fixed(SampleID, '-') == 3] <- mysample
  colnames(rdr.data) <- SampleID
  colnames(baf.data) <- SampleID
  colnames(mbaf.data) <- SampleID
  outputlist <- list(annotation = chisel_binanno, rdr_mat = rdr.data,
                     baf_mat = baf.data, mbaf_mat = mbaf.data)
  return(outputlist)
}

