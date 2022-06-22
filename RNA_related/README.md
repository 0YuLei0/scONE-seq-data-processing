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
```
# In R
library(tximport)
tx2gene <- read.delim("d:/YULei/Reference/esiaR_seperate/gcv34.eisaR_seperate.tx2gene.tsv",as.is = T,header = F)
head(tx2gene)
tail(tx2gene)
cg <- read.delim("d:/YULei/Reference/esiaR_seperate/gcv34.eisaR_seperate.intron2exon.tsv",header = T, as.is = TRUE)
head(cg)
## path to quant files
dir <- "d:/YULei/Data/GBM Part/Analyze/CleanRNA_Scount"
list.files(dir)
mysample <- list.files(dir)
samples<-as.data.frame(mysample)
files <- file.path(dir,samples[,1],"quant.sf")
coldata = data.frame(names = mysample, files = files, stringsAsFactors = FALSE)
## import data
txi2 <- tximport(files, type = "salmon",tx2gene = tx2gene)
counts<-txi2$counts;colnames(counts) <- list.files(dir)
spliced <- counts[rownames(counts) %in% cg$spliced,]
unspliced <- counts[rownames(counts) %in% cg$intron,]
spliced <- spliced[rowSums(spliced > 0) > 1,]
unspliced <- unspliced[rowSums(unspliced > 0) > 1,]
dim(spliced)
dim(unspliced)
intron2gene <- unlist(strsplit(rownames(unspliced),split = "-"))
intron2gene <- intron2gene[seq(1, length(intron2gene), by = 2)]
rownames(unspliced) <- intron2gene
dim(unspliced)
## CHECK "spliced, unspliced" counts number
a <- data.frame(type = rep(c("spliced","unspliced"), each =384),
                Counts = c(colSums(spliced),colSums(unspliced)))
ggplot() + geom_boxplot(data = a, aes(y = log10(Counts),fill = type))
## Transfer to havana gene name
havana <- read.csv("../../esemble2havana.csv", header = F)
library(data.table)
transf_havana <- function(x) {
  expr_data <- as.data.frame(x)
  expr_data$V1 <- rownames(expr_data)
  expr_data<- left_join(expr_data,havana)
  #expr_data %>% select(-V1) %>% group_by(V2) %>% summarise_all(funs(sum)) -> expr_data
  expr_data$V1 <- NULL
  expr_data <- as.data.table(expr_data)
  expr_data <- expr_data[,lapply(.SD, sum),by = V2]
  return(expr_data)
}
spliced_counts <- transf_havana(spliced)
unspliced_counts <- transf_havana(unspliced)
## Get PreRNA exprdata from "spliced_counts" + "unspliced_counts"
a <- rbind(spliced_counts, unspliced_counts)
havana <- rownames(a)
a <- as.data.table(a)
a$V2 <- havana
a <- a[,lapply(.SD, sum),by = V2]
prerna_counts <- as.matrix(a[,-1])
rownames(prerna_counts) <- a$V2
## Get "spliced_counts" and "unspliced_counts"
spliced <- as.matrix(spliced_counts[,-1]); rownames(spliced) <- spliced_counts$V2
unspliced <- as.matrix(unspliced_counts[,-1]); rownames(unspliced) <- unspliced_counts$V2
# CHECK if there is any mistakes
a <- data.frame(type = rep(c("spliced","unspliced"), each =384),
                Counts = c(colSums(spliced_counts[,-1]),colSums(unspliced_counts[,-1])))
ggplot() + geom_boxplot(data = a, aes(y = log10(Counts),fill = type))
rm(cg,a,bfc,coldata,prerna_counts,expr_data,havana,newmat,samples,spliced_counts, unspliced_counts, tx2gene, txi, txi2, bfcloc, dir,files, intron2gene, mysample, transf_havana,counts)
## Now, we have spliced, unspliced, and prerna_counts
dim(spliced);dim(unspliced);dim(prerna_counts)
```
