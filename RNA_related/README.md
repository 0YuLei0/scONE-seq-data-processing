## Genome reference
We used gencode.v34 for RNA aligment and quantification.
To perform snRNA-seq analysis, we followed the pipline from https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/.
We used the followed code to build salmon reference.
```
# In R
# use eisaR and "separate" stratedgy to extract the intron
gtf <- "$PATH/gencode.v34.primary_assembly.annotation.gtf"
grl <- eisaR::getFeatureRanges(
  gtf = gtf,
  featureType = c("spliced", "intron"),
  intronType = "separate",
  flankLength = 40L,
  joinOverlappingIntrons = FALSE,
  verbose = TRUE
)
genome <- Biostrings::readDNAStringSet("$PATH/GRCh38.primary_assembly.genome.fa")
names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)
library(BSgenome) # to enable the "GenomicFeatures::extractTranscriptSeqs", BSgenome has to be loaded
seqs <- GenomicFeatures::extractTranscriptSeqs(genome, grl)
Biostrings::writeXStringSet(seqs, filepath = "gv34.annotation.expanded.fa")
eisaR::exportToGtf(grl, filepath = "gv34.annotation.expanded.gtf")
```
```
# In linux
# Use Salmon to prepare the index
grep ">" GRCh38.primary_assembly.genome.fa | cut -d ">" -f 2 | cut -d " " -f 1 > GRCh38.primary_assembly.genome.chrnames.txt
cat gv34.annotation.expanded.fa GRCh38.primary_assembly.genome.fa > salmon_gv34esiaR_gentrome.fa
nohup salmon index -t salmon_gv34esiaR_gentrome.fa -d GRCh38.primary_assembly.genome.chrnames.txt -p 64 -i gv34_prerna_s14 --gencode &
```
## Summary intron and exon expression
Once you have finished the quantification with salmon, you can import the quantification to R.
1. We first summarise the transcript-level expression to gene-level with tximport R package.
```
# In R
# Import transcript annotation files
library(tximport)
library(data.table)
tx2gene <- read.delim("gcv34.eisaR_seperate.tx2gene.tsv",as.is = T,header = F)
cg <- read.delim("gcv34.eisaR_seperate.intron2exon.tsv",header = T, as.is = TRUE)
havana <- read.csv("esemble2havana.csv", header = F)
# path to quant files
dir <- "d:/YULei/Data/GBM Part/Analyze/CleanRNA_Scount"
samples<- as.data.frame(list.files(dir))
files <- file.path(dir,samples[,1],"quant.sf")
## import data
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
counts<-txi$counts
colnames(counts) <- list.files(dir)
spliced <- counts[rownames(counts) %in% cg$spliced,]
unspliced <- counts[rownames(counts) %in% cg$intron,]
```
2. We next combine intron and exon expression. After running the followed code, you could get 3 matrixs: **prerna_mat**(summarised expression), **spliced_mat**(exon expression), and **unspliced_mat**(intron expression),
```
# Rename intron to thier transcript ID
intron2gene <- unlist(strsplit(rownames(unspliced),split = "-"))
intron2gene <- intron2gene[seq(1, length(intron2gene), by = 2)]
rownames(unspliced) <- intron2gene
# Rename Ensembl ID to havana gene name
transf_havana <- function(x) {
  expr_data <- as.data.frame(x)
  expr_data$V1 <- rownames(expr_data)
  expr_data<- left_join(expr_data,havana)
  expr_data$V1 <- NULL
  expr_data <- as.data.table(expr_data)
  expr_data <- expr_data[,lapply(.SD, sum),by = V2]
  return(expr_data)
}
spliced_counts <- transf_havana(spliced)
unspliced_counts <- transf_havana(unspliced)
## Get PreRNA(intron+exon) expression data from "spliced_counts" + "unspliced_counts"
a <- rbind(spliced_counts, unspliced_counts)
havana <- rownames(a)
a <- as.data.table(a)
a$V2 <- havana
a <- a[,lapply(.SD, sum),by = V2]
prerna_mat <- as.matrix(a[,-1])
rownames(prerna_counts) <- a$V2
## Get "spliced_counts" and "unspliced_counts"
spliced_mat <- as.matrix(spliced_counts[,-1]); rownames(spliced) <- spliced_counts$V2
unspliced_mat <- as.matrix(unspliced_counts[,-1]); rownames(unspliced) <- unspliced_counts$V2
```
