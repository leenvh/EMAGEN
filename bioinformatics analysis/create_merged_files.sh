ref=$(realpath $1)

gatk CombineGVCFs -R $ref -O merged.g.vcf.gz `find . | grep variants.g.vcf.gz$ | awk '{print "-V "$1}' | tr '\n' ' '`
gatk GenotypeGVCFs -R $ref -O genotyped.vcf.gz -V merged.g.vcf.gz
source $(conda info --base)/etc/profile.d/conda.sh
conda deactivate
conda activate base
snpEff ann pfalciparum  genotyped.vcf.gz -no-upstream -no-downstream | bcftools view -Oz -o genotyped.ann.vcf.gz
