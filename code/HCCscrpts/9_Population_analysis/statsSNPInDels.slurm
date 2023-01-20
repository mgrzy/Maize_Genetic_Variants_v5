#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=2-10:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --job-name=stats
#SBATCH --error=../job/stats.%A_%a.out
#SBATCH --output=../job/stats.%A_%a.err
#SBATCH --array=1

ml bcftools
ml tabix

chr=$SLURM_ARRAY_TASK_ID

bcftools +fill-tags ../VCF/filtered/All_filtered_snps.vcf.gz -Ov | bcftools query -f '%POS %DP %ExcessHet %InbreedingCoeff %MAF %NS\n' | gzip -f > ../results/stats/stat_${chr}_stats_snps.txt.gz

bcftools +fill-tags ../VCF/filtered/All_filtered_indels.vcf.gz -Ov | bcftools query -f '%POS %DP %ExcessHet %InbreedingCoeff %MAF %NS\n' | gzip -f > ../results/stats/stat_${chr}_stats_indels.txt.gz
