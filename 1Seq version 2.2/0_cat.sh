#!/bin/bash
## cat fastq
workspace=$1
mkdir -p $workspace/catfiles
cd $workspace/raw
for i in `ls -d *L001*`
do
  echo "cating sample is $i"
  a=$(echo $i | cut -d '-' -f 2)
  b=$(echo $a | cut -d '_' -f 1)
  c=$(echo $a | cut -d '_' -f 2)
  d=$(echo $a | cut -d '_' -f 3)
  newname=$(echo "$b"-"$c"-"$d")
  find `echo *"$b"_"$c"_"$d"*` -name "*R1_001.fastq.gz" -exec cat {} + > ../catfiles/$newname\_read1.fastq.gz
  find `echo *"$b"_"$c"_"$d"*` -name "*R2_001.fastq.gz" -exec cat {} + > ../catfiles/$newname\_read2.fastq.gz
done

for i in `ls -d *L001*`
do
  echo "cating sample is $i"
  a=$(echo $i | cut -d '-' -f 3)
  b=$(echo $a | cut -d '_' -f 1)
  c=$(echo $a | cut -d '_' -f 2)
  d=$(echo $a | cut -d '_' -f 3)
  e=$(echo $a | cut -d '_' -f 4)
  newname=$(echo "$b"-"$c"-"$d"-"$e")
  find `echo *"$b"_"$c"_"$d"_"$e"*` -name "*R1_001.fastq.gz" -exec cat {} + > ../catfiles/$newname\_read1.fastq.gz
  find `echo *"$b"_"$c"_"$d"_"$e"*` -name "*R2_001.fastq.gz" -exec cat {} + > ../catfiles/$newname\_read2.fastq.gz
done


for i in `ls -d 20200617*L001*`
do
a=$(echo $i | cut -d '_' -f 1)
b=$(echo $a | cut -d '-' -f 2)
c=$(echo $a | cut -d '-' -f 3)
d=$(echo $a | cut -d '-' -f 4)
newname=$(echo "$b"-"$c"-"$d")
find `echo *"$b"-"$c"-"$d"*` -name "*R1_001.fastq.gz" -exec cat {} + > ../catfiles/$newname\_read1.fastq.gz
find `echo *"$b"-"$c"-"$d"*` -name "*R2_001.fastq.gz" -exec cat {} + > ../catfiles/$newname\_read2.fastq.gz
done
