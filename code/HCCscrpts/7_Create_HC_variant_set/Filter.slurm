#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=2
#SBATCH --nodes=1
#SBATCH --time=2-10:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --job-name=filter001
#SBATCH --error=../job/final/filter/filter.%A_%a.out
#SBATCH --output=../job/final/filter/filter.%A_%a.err
#SBATCH --array=1-10

ml bcftools

chr=$SLURM_ARRAY_TASK_ID

bcftools annotate ../VCF/unimputed/chr_${chr}_filtered_DP2.vcf.gz -x 'FORMAT' |\
bcftools view -i 'F_MISSING<0.5 & InbreedingCoeff>0 & MAF>0.01 & DP<33550 & DP>1515' -m2 -M2 -Ov |\
bcftools norm -d all -Oz -o ../VCF/filtered/chr_${chr}_filtered_MAF001.vcf.gz

ml tabix
tabix ../VCF/filtered/chr_${chr}_filtered_MAF001.vcf.gz
