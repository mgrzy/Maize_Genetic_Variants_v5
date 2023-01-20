#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=2
#SBATCH --nodes=1
#SBATCH --time=5-00:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=Fchr7
#SBATCH --error=../job/filter/chr7.%A_%a.out
#SBATCH --output=../job/filter/chr7.%A_%a.err
#SBATCH --array=1-185

ml java/1.8
ml bcftools
ml gatk4
ml tabix

a="../meta/VCFfiles/chr7.txt"
inPath="../VCF/raw/chr7/"
outPath="../VCF/work/chr7/"

File=`sed -n "$SLURM_ARRAY_TASK_ID"p $a | awk '{print $1}'`

# Filter INDELS
gatk SelectVariants \
-V ${inPath}${File}.vcf.gz \
-select-type INDEL \
-O /scratch/${File}_indels.vcf.gz

gatk VariantFiltration \
-V /scratch/${File}_indels.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "FS > 200.0" --filter-name "FS200" \
-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
-O /scratch/${File}_indels_filtered.vcf.gz

bcftools +setGT /scratch/${File}_indels_filtered.vcf.gz -- -t q -n . -i 'FMT/DP=0' | bcftools view -f PASS -Ov | bcftools annotate --set-id +'%CHROM\_%POS' -O z -o ${outPath}${File}_indels_filtered.vcf.gz
tabix ${outPath}${File}_indels_filtered.vcf.gz

# Filter SNPs

gatk SelectVariants \
-V ${inPath}${File}.vcf.gz \
-select-type SNP \
-O /scratch/${File}_snp.vcf.gz

gatk VariantFiltration \
-V /scratch/${File}_snp.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O /scratch/${File}_snp_filtered.vcf.gz

bcftools +setGT /scratch/${File}_snp_filtered.vcf.gz -- -t q -n . -i 'FMT/DP=0' | bcftools view -f PASS -Ov | bcftools annotate --set-id +'%CHROM\_%POS' -O z -o ${outPath}${File}_snp_filtered.vcf.gz
tabix ${outPath}${File}_snp_filtered.vcf.gz
