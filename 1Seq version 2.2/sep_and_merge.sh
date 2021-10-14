## based anotaiton, merge cells based on their cell types
for i in `cat list`
do
  echo "searching $i*.bam"
  ## be careful about the last "\;"
  find `echo ~/GBM/*batch*/umidedup/umi_mapped/umi_bwabam` -name `echo $i*.bam` -exec cp {} /home/lyuah/GBM/clones/normalcells/$i\.bam \;
done

for i in `cat list`
do
  echo "searching $i*.bam"
  ## be careful about the last "\;"
  find `echo ~/GBM/*batch*/umidedup/umi_mapped/umi_bwabam` -name `echo $i*.bam` -exec cp {} /home/lyuah/GBM/clones/Clone2/$i\.bam \;
done

for i in `cat list`
do
  echo "searching $i*.bam"
  ## be careful about the last "\;"
  find `echo ~/GBM/*batch*/umidedup/umi_mapped/umi_bwabam` -name `echo $i*.bam` -exec cp {} /home/lyuah/GBM/clones/Clone1/$i\.bam \;
done

for i in `cat list`
do
  echo "searching $i*.bam"
  ## be careful about the last "\;"
  find `echo ~/GBM/*batch*/umidedup/umi_mapped/umi_bwabam` -name `echo $i*.bam` -exec cp {} /home/lyuah/GBM/clones/Clone3_Cycling/$i\.bam \;
done

for i in `cat list`
do
  echo "searching $i*.bam"
  ## be careful about the last "\;"
  find `echo ~/GBM/*batch*/umidedup/umi_mapped/umi_bwabam` -name `echo $i*.bam` -exec cp {} /home/lyuah/GBM/clones/Clone3_GFAP/$i\.bam \;
done

for i in `cat list`
do
  echo "searching $i*.bam"
  ## be careful about the last "\;"
  find `echo ~/GBM/*batch*/umidedup/umi_mapped/umi_bwabam` -name `echo $i*.bam` -exec cp {} /home/lyuah/GBM/clones/Clone3_VEGFA/$i\.bam \;
done

for i in `cat list`
do
  echo "searching $i*.bam"
  ## be careful about the last "\;"
  find `echo ~/GBM/*batch*/umidedup/umi_mapped/umi_bwabam` -name `echo $i*.bam` -exec cp {} /home/lyuah/GBM/clones/Clone3_MEGF/$i\.bam \;
done

for i in `cat list`
do
  echo "searching $i*.bam"
  ## be careful about the last "\;"
  find `echo ~/GBM/*batch*/umidedup/umi_mapped/umi_bwabam` -name `echo $i*.bam` -exec cp {} /home/lyuah/GBM/clones/highqc_cells/$i\.bam \;
done

for i in `cat list`
do
  echo "searching $i*.bam"
  ## be careful about the last "\;"
  find `echo ~/GBM/*batch*/umidedup/umi_mapped/umi_bwabam` -name `echo $i*.bam` -exec cp {} /home/lyuah/GBM/clones/chisel/CancerCells/$i\.bam \;
done

for i in `ls`
do
  cd ~/GBM/clones/$i
  samtools merge -@ 36 $i.bam *bam
  samtools index $i.bam -@ 36
  cd ~/GBM/clones
done

## CALL SNPs from normal normal cells
## Use freebayes; strelka2 fail to produce usabel chr7 VCF
for i in {1..22}
do
  a=$(echo chr$i)
  echo "Call $a SNPs with freebayes"
  freebayes --no-indels --no-mnps --no-complex -f ~/reference/human_DNA/liheng_Genome_reference/liheng_human_hg38_bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --min-mapping-quality 50 -C 4 -r $a /home/lyuah/GBM/clones/normalcells/normalcells.bam > $a\_free.vcf
done
## filter freebayes
## Assumption:
#  1. Normal Cells all have the same genome;
#  2. Heterozygous SNP should have frequency around 0.5
for i in `ls *vcf`
do
  awk '/#/' $i > temp.vcf
  awk '/TYPE=snp/' $i | awk -F ":" '{if($8 > 0 && $9/$8 > 0.39 && $9/$8 < 0.61) print $0}' >> temp.vcf
  mv temp.vcf ${i%_free.vcf}_snp.vcf
  bgzip ${i%_free.vcf}_snp.vcf
done
## Phasing with Michigan Imputation Server
###### Concat all chr from Phasing VCF files
## rename chr
for i in {1..22}
do
  echo "$i chr$i" >> chr_name_conv.txt
done
## filter hete
for i in `ls *dose*`
do
  bcftools filter -i 'GT!="0|0"' $i | bcftools filter -i 'GT!="1|1"' -Ov > ${i%.dose.vcf.gz}_temp.vcf
done
## concat all VCFs
bcftools concat chr1_temp.vcf chr2_temp.vcf chr3_temp.vcf chr4_temp.vcf chr5_temp.vcf chr6_temp.vcf chr7_temp.vcf chr8_temp.vcf chr9_temp.vcf chr10_temp.vcf chr11_temp.vcf chr12_temp.vcf chr13_temp.vcf chr14_temp.vcf chr15_temp.vcf chr16_temp.vcf chr17_temp.vcf chr18_temp.vcf chr19_temp.vcf chr20_temp.vcf chr21_temp.vcf chr22_temp.vcf --threads 36 -o all_temp.vcf
## replace chr names
bcftools annotate --rename-chrs chr_name_conv.txt all_temp.vcf | bgzip > hg38_phased.vcf.gz
## remove temp files
rm *temp*
## prepare shorter tsv for chisel
zcat hg38_phased.vcf.gz | sed '/^#/d' | awk '{print $1,$2,$10}' | awk -F ":" '{print $1}' > chisel.tsv
###### Liftover Phasing VCF back to hg38
# 1. Liftover hg19 to hg38
awk '{print $1"\t"$2"\t"$2"\t"$3}' chisel.tsv | sort -k1,1 -k2,2n  > temp0
# filtered blacklist and gap(gapped-in-hg19,hg19 overlap gap in hg38, gapped-in-both and blacklist)
bedtools subtract -a temp0 -b ~/software/Lifted/data/gapped-in-hg19.bed -A | sort -k1,1 -k2,2n > temp001
bedtools subtract -a temp001 -b ~/software/Lifted/data/hg19-overlap-gapped-in-hg38.bed -A | sort -k1,1 -k2,2n > temp002
bedtools subtract -a temp002 -b ~/software/Lifted/data/gapped-in-both.bed -A | sort -k1,1 -k2,2n > temp003
bedtools subtract -a temp003 -b ~/software/Lifted/data/blacklist.hg19.bed -A | sort -k1,1 -k2,2n > temp1
rm temp0*
### sort
sort -k1,1 -k2,2n temp1 > temp2
## liftover
sed '492240,493250d' temp2 > temp21
~/software/Lifted/bin/liftOver temp2 ~/software/Lifted/data/hg19ToHg38.over.chain temp3 umapped
#### filter  duplicate, alt chromosome
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$1"."$2"."$3}' temp3 | sort -k1,1 -k2,2n  > temp31
sort -k1,1 -k2,2n temp31 | grep -E "chr(.|..)[[:blank:]]" | sort -k5,5 | uniq -f4 -u | sort -k1,1 -k2,2n > temp002
# our modified bed files
awk '{print $1"\t"$2"\t"$4}' temp002 > lift_chisel.tsv
Phased_SNP = "/home/lyuah/GBM/clones/SNP_strelka/phasing/free_phased/lift_chisel.tsv"

############ Run single cell chisel
rm -rf rdr/ baf/ combo/ calls/ clones/ plots/
##### Psuedo-bulk chisel
## Prepare data
chisel_prep -r ~/reference/human_DNA/liheng_Genome_reference/liheng_human_hg38_bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -j 30 --seed 24 /home/lyuah/GBM/clones/Clone1/Clone1.bam /home/lyuah/GBM/clones/Clone2/Clone2.bam /home/lyuah/GBM/clones/Clone3_Cycling/Clone3_Cycling.bam /home/lyuah/GBM/clones/Clone3_GFAP/Clone3_GFAP.bam /home/lyuah/GBM/clones/Clone3_MEGF/Clone3_MEGF.bam /home/lyuah/GBM/clones/Clone3_VEGFA/Clone3_VEGFA.bam ~/GBM/clones/normalcells/normalcells.bam
## Test Run with Psuedo Bulk
chisel -t barcodedcells.bam -n ~/GBM/clones/normalcells/normalcells.bam -r ~/reference/human_DNA/liheng_Genome_reference/liheng_human_hg38_bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -l ~/GBM/clones/SNP_strelka/phasing/free_phased/lift_chisel.tsv --seed 66 -j 12 -b 5Mb -k 50kb
chisel_plotting -s 7
##### Single cell chisel
## 1. Prepare data
chisel_prep -r ~/reference/human_DNA/liheng_Genome_reference/liheng_human_hg38_bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -j 30 --seed 24 *bam
chisel_prep -r ~/reference/human_DNA/liheng_Genome_reference/liheng_human_hg38_bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -j 30 --seed 24 *bam
## 2. Run chiesl with highqc cells
chisel -t ~/GBM/clones/highqc_cells/barcodedcells.bam -n ~/GBM/clones/normalcells/normalcells.bam -r ~/reference/human_DNA/liheng_Genome_reference/liheng_human_hg38_bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -l /home/lyuah/GBM/clones/SNP_strelka/phasing/free_phased/lift_chisel.tsv --seed 66 -j 36 -b 2.5Mb -k 100kb -m 100000
## Run chiesl with Cancer cells
chisel -t ~/GBM/clones/chisel/CancerCells_1/barcodedcells.bam -n ~/GBM/clones/normalcells/normalcells.bam -r ~/reference/human_DNA/liheng_Genome_reference/liheng_human_hg38_bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -l /home/lyuah/GBM/clones/SNP_strelka/phasing/free_phased/lift_chisel.tsv --seed 66 -j 36 -b 5Mb -k 100kb -m 100000
## No need to rerun the time-conssuming SNP Calling step
python ~/anaconda3/chisel/lib/python2.7/site-packages/chisel/Combiner.py -h
python ~/anaconda3/chisel/lib/python2.7/site-packages/chisel/Combiner.py -r ~/GBM/clones/chisel/CancerCells_1/rdr/rdr.tsv -b ~/GBM/clones/chisel/CancerCells_1/baf/baf.tsv -j 24 -k 200kb > ./combo/combo.tsv
## Adjust
chisel_calling -j 36 --seed 66 -P 6
chisel_cloning --maxdiff 0.4 --minsize 10

## All cells
find `echo ~/GBM/*batch*/umidedup/umi_mapped/umi_bwabam` -name `echo *.bam` -exec cp {} ~/GBM/clones/chisel/allcells/ \;
chisel_prep -r ~/reference/human_DNA/liheng_Genome_reference/liheng_human_hg38_bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -j 30 --seed 24 *bam
chisel -t barcodedcells.bam -n ~/GBM/clones/normalcells/normalcells.bam -r ~/reference/human_DNA/liheng_Genome_reference/liheng_human_hg38_bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -l /home/lyuah/GBM/clones/SNP_strelka/phasing/free_phased/lift_chisel.tsv --seed 66 -j 36 -b 2.5Mb -k 50kb -m 100000
chisel_plotting -s 7

samtools view -H ~/GBM/bulk/DNAbwa/human_dnabam/P1B.bam | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq

samtools addreplacerg --output-fmt BAM -m overwrite_all -r '@RG\tID:P1B\tLB:P1B\tPL:ILLUMINA\tSM:P1B' -@ 36 ~/GBM/bulk/DNAbwa/human_dnabam/P1B.bam -o ~/GBM/bulk/DNAbwa/human_dnabam/P1B_head.bam
samtools addreplacerg --output-fmt BAM -m overwrite_all -r '@RG\tID:T808\tLB:T808\tPL:ILLUMINA\tSM:T808' -@ 36 ~/GBM/bulk/DNAbwa/human_dnabam/T808DNA.bam -o ~/GBM/bulk/DNAbwa/human_dnabam/GBM_T808.bam
samtools view -H ~/GBM/bulk/DNAbwa/human_dnabam/GBM_T808.bam | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq
samtools index ~/GBM/bulk/DNAbwa/human_dnabam/GBM_T808.bam -@ 36

gatk Mutect2 -R ~/reference/human_DNA/liheng_Genome_reference/liheng_human_hg38_bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -I ~/GBM/bulk/DNAbwa/human_dnabam/GBM_P1B.bam -I ~/GBM/bulk/DNAbwa/human_dnabam/GBM_T566.bam -I ~/GBM/bulk/DNAbwa/human_dnabam/GBM_T808.bam -I ~/GBM/bulk/DNAbwa/human_dnabam/GBM_T946.bam -normal P1B -L chrM -O GBM_WES_somatic.vcf.gz

samtools addreplacerg --output-fmt BAM -m overwrite_all -r '@RG\tID:Clone1\tLB:Clone1\tPL:ILLUMINA\tSM:Clone1' -@ 36 Clone1.bam -o Clone1_head.bam
samtools addreplacerg --output-fmt BAM -m overwrite_all -r '@RG\tID:Clone2\tLB:Clone2\tPL:ILLUMINA\tSM:Clone2' -@ 36 Clone2.bam -o Clone2_head.bam
samtools addreplacerg --output-fmt BAM -m overwrite_all -r '@RG\tID:Clone3_GFAP\tLB:Clone3_GFAP\tPL:ILLUMINA\tSM:Clone3_GFAP' -@ 36 Clone3_GFAP.bam -o Clone3_GFAP_head.bam
samtools addreplacerg --output-fmt BAM -m overwrite_all -r '@RG\tID:Clone3_MEGF\tLB:Clone3_MEGF\tPL:ILLUMINA\tSM:Clone3_MEGF' -@ 36 Clone3_MEGF.bam -o Clone3_MEGF_head.bam
samtools addreplacerg --output-fmt BAM -m overwrite_all -r '@RG\tID:Clone3_VEGFA\tLB:Clone3_VEGFA\tPL:ILLUMINA\tSM:Clone3_VEGFA' -@ 36 Clone3_VEGFA.bam -o Clone3_VEGFA_head.bam

samtools addreplacerg --output-fmt BAM -m overwrite_all -r '@RG\tID:scONEseq_Astrocyte\tLB:scONEseq_Astrocyte\tPL:ILLUMINA\tSM:scONEseq_Astrocyte' -@ 36 scONEseq_Astrocyte.bam -o scONEseq_Astrocyte_head.bam
samtools addreplacerg --output-fmt BAM -m overwrite_all -r '@RG\tID:scONEseq_Macrophage\tLB:scONEseq_Macrophage\tPL:ILLUMINA\tSM:scONEseq_Macrophage' -@ 36 scONEseq_Macrophage.bam -o scONEseq_Macrophage_head.bam
samtools addreplacerg --output-fmt BAM -m overwrite_all -r '@RG\tID:scONEseq_Neuron\tLB:scONEseq_Neuron\tPL:ILLUMINA\tSM:scONEseq_Neuron' -@ 36 scONEseq_Neuron.bam -o scONEseq_Neuron_head.bam
samtools addreplacerg --output-fmt BAM -m overwrite_all -r '@RG\tID:scONEseq_Oligodendrocyte\tLB:scONEseq_Oligodendrocyte\tPL:ILLUMINA\tSM:scONEseq_Oligodendrocyte' -@ 36 scONEseq_Oligodendrocyte.bam -o scONEseq_Oligodendrocyte_head.bam

samtools index -@ 40 Clone1_head.bam
samtools index -@ 40 Clone2_head.bam
samtools index -@ 40 Clone3_GFAP_head.bam
samtools index -@ 40 Clone3_MEGF_head.bam
samtools index -@ 40 Clone3_VEGFA_head.bam
samtools index -@ 40 scONEseq_Astrocyte_head.bam
samtools index -@ 40 scONEseq_Macrophage_head.bam
samtools index -@ 40 scONEseq_Neuron_head.bam
samtools index -@ 40 scONEseq_Oligodendrocyte_head.bam

gatk Mutect2 -R ~/reference/human_DNA/liheng_Genome_reference/liheng_human_hg38_bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -I scONEseq_Neuron_head.bam -I scONEseq_Macrophage_head.bam -I Clone1_head.bam -I Clone3_GFAP_head.bam -I Clone3_MEGF_head.bam -I Clone3_VEGFA_head.bam -I scONEseq_Astrocyte_head.bam -I scONEseq_Oligodendrocyte_head.bam -normal scONEseq_Neuron -normal scONEseq_Macrophage -O GBM_oneseq_somatic.vcf.gz

unpigz GBM_oneseq_somatic.vcf.gz
bgzip GBM_oneseq_somatic.vcf
bcftools index -t GBM_oneseq_somatic.vcf.gz

gatk FilterMutectCalls -V GBM_oneseq_somatic.vcf.gz -R ~/reference/human_DNA/liheng_Genome_reference/liheng_human_hg38_bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -O filtered.vcf.gz
bcftools view -v snps filtered.vcf 

table_annovar.pl --vcfinput ~/GBM/clones/Type_Bam/withhead/vcffiles/test.vcf humandb/ -buildver hg38 -out clonal_ --remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation g,r,f,f,f -nastring . -polish

vcf-query ~/software/annovar/clonal_.hg38_multianno.vcf -l > ~/software/annovar/sampleid

vcf-query --list-columns clonal_.hg38_multianno.vcf > sample_ids
for i in `cat ~/software/annovar/sampleid`
do
  echo $i
  cat ~/software/annovar/clonal_.hg38_multianno.vcf | vcf-subset --exclude-ref --columns $i > ~/software/annovar/vcf2maf/$i\.vcf
done
for i in `cat ~/software/annovar/sampleid`
do
  vcf2maf.pl --input-vcf ~/software/annovar/vcf2maf/$i\.vcf --output-maf ~/software/annovar/vcf2maf/$i\.maf --tumor-id $i --ref-fasta ~/reference/human_DNA/liheng_Genome_reference/liheng_human_hg38_bwa/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --inhibit-vep
done

cat *.maf | egrep "^#|^Hugo_Symbol" | head -2 > allsamples.vep.maf
cat *.maf | egrep -v "^#|^Hugo_Symbol" >> allsamples.vep.maf
