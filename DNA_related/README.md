## Genome reference
We used hg38 suggested by http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use
## Get normlized counts from Ginkgo
We used Ginkgo(https://github.com/robertaboukhalil/ginkgo) to normalize our scDNA data from scONE-seq. To use hg38 human reference, you should setup Ginkgo on your own server and build new genome reference for Ginkgo forllowing the instruction in (https://github.com/robertaboukhalil/ginkgo/tree/master/genomes/scripts).
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
process_segnorm($PATH_to_SegNorm)
```
This function output a list object with Segmentation result and CNV result (assuming ploidy to be 2).
Segmentation result is also a list object, which contains Segmentation result from "copynumber" and mergeLevels result from "aCGH".
## Get bin-level allele freguency with CHISEL
We used CHISEL to get scDNA allele frequency information (https://github.com/raphael-group/chisel). 
1. With counts-based CNVs, we could distinguish normal cells and malignent cells.
2. Prepare phased germline SNPs (normal cells were combined for germline SNPs calling). SNPs were identified with your favorate mutation caller (we simpliy used bcftools). Followed the instruction of CHISEL, these SNPs were phased with Michigan Imputation Server (https://imputationserver.sph.umich.edu/index.html#!pages/home).
```
## From Michigan Imputation Server, you should get phased VCF fow each chromsome
## Filter and combine them
for i in `ls *dose*`
do
bcftools filter -i 'GT!="0|0"' $i | bcftools filter -i 'GT!="1|1"' -Ov | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' | sed '/^#/d' | awk '{print $1,$2,$10}' | awk -F ":" '{print $1}' > ${i%.dose.vcf}_chisel.tsv
done
cat chr1_chisel.tsv chr2_chisel.tsv chr3_chisel.tsv chr4_chisel.tsv chr5_chisel.tsv chr6_chisel.tsv chr7_chisel.tsv chr8_chisel.tsv chr9_chisel.tsv chr10_chisel.tsv chr11_chisel.tsv chr12_chisel.tsv chr13_chisel.tsv chr14_chisel.tsv chr15_chisel.tsv chr16_chisel.tsv chr17_chisel.tsv chr18_chisel.tsv chr19_chisel.tsv chr20_chisel.tsv chr21_chisel.tsv chr22_chisel.tsv > hg38_phased_chisel.tsv
```
3. Combine tumor cells togather.
```
## prepare barcode bam
chisel_prep -r $PATH_to_your_hg38 -j 20 --seed 24 [all your single cell bam files]
```
4. Run chisel main program
```
chisel -t barcodedcells.bam -n normal.bam -r $PATH_to_your_hg38 -l hg38_phased_snps.tsv
```
## Interge CNV calculation (combining counts and allele freguency)
In this part, we tried to caculated the interge CNVs considering the allele freguency infromation infered from CHISEL. This could be especially useful when working with tumor cells which have 1/3 allele freguency. 
```
combo1 <- read.delim("D:/YULei/Data/GBM Part/Analyze/Chisel/AllCells/Run1/combo/combo.tsv",header = F) ## all cells 1
cellid1 <- read.delim("D:/YULei/Data/GBM Part/Analyze/Chisel/AllCells/Run1/barcodedcells.info.tsv") ## all cells 1
colnames(combo1) <- c("CHR","START","END","BARCODE","N_reads","Bin_reads","RDR","A_counts","B_Counts","BAF")
chisel_run1 <- extract_bafANDrdr(combo1,cellid1,cellid_trim = 10)

```
