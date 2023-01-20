#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=30
#SBATCH --nodes=1
#SBATCH --time=2-10:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --job-name=LDscore
#SBATCH --error=../job/analysis/LDscore.%A_%a.out
#SBATCH --output=../job/analysis/LDscore.%A_%a.err
#SBATCH --array=1

ml anaconda
conda activate gcta

gcta64 --bfile ../plink/All --ld-score --ld-wind 1000 --out ../results/LDscore --thread-num 30
