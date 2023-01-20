#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --job-name=pca
#SBATCH --error=../job/analysis/pca.%A_%a.out
#SBATCH --output=../job/analysis/pca.%A_%a.err
#SBATCH --array=1

ml plink

plink --bfile ../plink/All_MAF001 --indep-pairwise 500 100 0.1

plink --bfile ../plink/All_MAF001 --extract plink.prune.in --make-bed --out ../plink/AllPruned

plink --bfile ../plink/All_MAF001 --pca 1515 --out ../results/pca/PCAallTaxa
