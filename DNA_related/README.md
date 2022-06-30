## Genome reference
We used hg38 suggested by http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use
## Get normalized counts from Ginkgo
We used Ginkgo(https://github.com/robertaboukhalil/ginkgo) to normalize our scDNA data from scONE-seq. To use hg38 human reference, you should setup Ginkgo on your own server and build a new genome reference for Ginkgo following the instruction in (https://github.com/robertaboukhalil/ginkgo/tree/master/genomes/scripts).
Once you have Ginkgo installed, you could 
1. Create a directory under the uploads directory in the ginkgo installation directory
```
cd $Ginkgo
mkdir uploads
mkdir uploads/sconeseq
```
2. Move bed.gz to the uploads directory
```
mv *bed.gz $Ginkgo/uploads/sconeseq
```
3. Create a file with the list of cells. The file must be called "list":
```
ls | grep .bed.gz$ > list
```
4. Create a configuration file with all the options for ginkgo. There is an example file named "config.example" in Ginkgo(https://github.com/robertaboukhalil/ginkgo).
```
nohup bash $Ginkgo/scripts/analyze.sh $Ginkgo/uploads/sconeseq &
```
## Interge CNV calculation (counts-based)
The normlized counts data from Ginkgo stored in **SegNorm** file in $Ginkgo/uploads/sconeseq. We then used R:copynumber and R:aCGH to peform segmentation.
```
# In R
source(mergeLevels_multi_segmentation.R)
source(process_segnorm.R)
processed_counts_data <- process_segnorm($PATH_to_SegNorm)
```
This function output a list object with Segmentation result and CNV result (assuming ploidy to be 2).
Segmentation result is also a list object, which contains Segmentation result from "copynumber" and mergeLevels result from "aCGH".
## Get bin-level allele frequency with CHISEL
We used CHISEL to get scDNA allele frequency information (https://github.com/raphael-group/chisel). 
1. With counts-based CNVs, we could distinguish normal cells and malignant cells.
2. Prepare phased germline SNPs (normal cells were combined for germline SNPs calling). SNPs were identified with your favorite mutation caller (we simply used bcftools). Followed the instruction of CHISEL, these SNPs were phased with Michigan Imputation Server (https://imputationserver.sph.umich.edu/index.html#!pages/home).
```
## From Michigan Imputation Server, you should get phased VCF fow each chromsome
## Filter and combine them
for i in `ls *.vcf`
do
bcftools filter -i 'GT!="0|0"' $i | bcftools filter -i 'GT!="1|1"' -Ov | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' | sed '/^#/d' | awk '{print $1,$2,$10}' | awk -F ":" '{print $1}' > ${i%.dose.vcf}_chisel.tsv
done
cat chr1_chisel.tsv chr2_chisel.tsv chr3_chisel.tsv chr4_chisel.tsv chr5_chisel.tsv chr6_chisel.tsv chr7_chisel.tsv chr8_chisel.tsv chr9_chisel.tsv chr10_chisel.tsv chr11_chisel.tsv chr12_chisel.tsv chr13_chisel.tsv chr14_chisel.tsv chr15_chisel.tsv chr16_chisel.tsv chr17_chisel.tsv chr18_chisel.tsv chr19_chisel.tsv chr20_chisel.tsv chr21_chisel.tsv chr22_chisel.tsv > hg38_phased_chisel.tsv
```
3. Combine tumor cells together.
```
## prepare barcode bam
chisel_prep -r $PATH_to_your_hg38 -j 20 --seed 24 [all your single cell bam files]
```
4. Run chisel main program
```
chisel -t barcodedcells.bam -n normal.bam -r $PATH_to_your_hg38 -l hg38_phased_chisel.tsv
```
## Interge CNV calculation (combining counts and allele frequency)
In this part, we tried to calculate the integer CNVs considering the allele frequency information inferred from CHISEL. This step can help us to remove some miscalling in single cell data.
```
# Load required pakgages and 
source(extract_baf_from_chisel.R)
source(Transfer_CHISEL_mbaf_to_ginkgo.R)
source(Calculate_integerCNV.R)
library(GenomicRanges)
library(IRanges)
library(tibble)
## Import BAF data from CHISEL result, this data is located in combo/combo.tsv of the CHISEL result
combo <- read.delim($PATH_to_combo.tsv,header = F)
colnames(combo) <- c("CHR","START","END","BARCODE","N_reads","Bin_reads","RDR","A_counts","B_Counts","BAF")
## Import Cellid, CHISEL made pseudo-cellid with chisel_prep, which generated a "barcodedcells.info.tsv" files
cellid <- read.delim("$PATH_to_barcodedcells.info.tsv")
## Convert this data to a matrix format
chisel_run <- extract_bafANDrdr(combo, cellid, cellid_trim = 0)
## convert chisel bins (5Mb) to ginkgo bins (500kb)
ginkgo_mbaf <- chisel2ginkgo(ginkgo_bins = bin_anno, chisel_bins = chisel_run$annotation, mbaf = chisel_run$mbaf_mat)
## Only works for our tetraploid sample now.
## For samples from other sources, this code needs to be modified.
## Calling interge CNV with BAF information
ginkgo_bafcnv <- baf_cnv_ginkgo(rdr = ginkgo_rdr, # normalized counts from Ginkgo
                                mbaf = ginkgo_mbaf, # mbaf from CHISEL
                                fixed_rdr = ginkgo_rdr_fixed, # Segmentation result from "processed_counts_data"
                                contral_baf = baf_distr,
                                contral_rdr = rdr_distr
                                )
```
