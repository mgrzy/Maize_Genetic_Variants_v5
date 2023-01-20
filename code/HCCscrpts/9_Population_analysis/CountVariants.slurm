#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=15G
#SBATCH --job-name=count
#SBATCH --error=../job/count.%A_%a.out
#SBATCH --output=../job/count.%A_%a.err
#SBATCH --array=1

ml bcftools
ml tabix
ml bedtools

bcftools query -f'%CHROM\t%POS\t%POS\n' ../VCF/filtered/All_filtered_indels.vcf.gz | bgzip -c > ../meta/indelsPos.bed.gz

tabix ../meta/indelsPos.bed.gz

bedtools coverage -a ../meta/windows.bed.gz -b ../meta/indelsPos.bed.gz -counts > ../results/IndelssPer100kb.txt

bcftools query -f'%CHROM\t%POS\t%POS\n' ../VCF/filtered/All_filtered_snps.vcf.gz | bgzip -c > ../meta/snpsPos.bed.gz

tabix ../meta/snpsPos.bed.gz

bedtools coverage -a ../meta/windows.bed.gz -b ../meta/snpsPos.bed.gz -counts > ../results/SNPsPer100kb.txt
bedtools coverage -a ../meta/windows1MB.bed.gz -b ../meta/snpsPos.bed.gz -counts > ../results/SNPsPer1Mb.txt
