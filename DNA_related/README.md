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
