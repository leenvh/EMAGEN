#!/bin/bash/
ref=$(realpath $1)
r1=$(realpath $2)
r2=$(realpath $3)
sample_name=$4
script_dir_rel=$(dirname $0)
script_dir=$(realpath $script_dir_rel)

echo $ref
mkdir $sample_name
cd $sample_name

bwa mem -t 48 -R "@RG\tID:$sample_name\tSM:$sample_name\tPL:Illumina" $ref $r1 $r2 | samtools sort -@ 48 -o alignment.bam
samtools index alignment.bam
python $script_dir/downsample.py --input alignment.bam --bed $script_dir/intervals.bed --output subsampled.bam
gatk HaplotypeCaller -R $ref -I subsampled.bam -O variants.g.vcf.gz -L $script_dir/intervals.bed -ERC GVCF

