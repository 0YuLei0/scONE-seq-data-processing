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

## Get bin-level allele freguency with CHISEL
We used CHISEL to get scDNA allele frequency information (https://github.com/raphael-group/chisel). 
1. With counts-based CNVs, we could distinguish normal cells and malignent cells.
2. Prepare phased germline SNPs (normal cells were combined for germline SNPs calling). SNPs were identified with your favorate mutation caller (we simpliy used bcftools). Followed the instruction of CHISEL, these SNPs were phased with Michigan Imputation Server (https://imputationserver.sph.umich.edu/index.html#!pages/home).
```

```
3. Combine tumor cells togather.
```
## prepare barcode bam
chisel_prep -r $PATH_to_your_hg38 -j 20 --seed 24 [all your single cell bam files]
```
4. Run chisel main program
```
chisel -t barcodedcells.bam -n normal.bam -r $PATH_to_your_hg38 -l hg38_phased_snps.tsv --seed 12 -j 24 -m 100000 -b 5Mb -k 50kb
```
## Interge CNV calculation (combining counts and allele freguency)
