#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=100G
#SBATCH --job-name=LDwhole
#SBATCH --error=../job/analysis/popLD.%A_%a.out
#SBATCH --output=../job/analysis/popLD.%A_%a.err
#SBATCH --array=1-10

chr=$SLURM_ARRAY_TASK_ID

$HOME/soft/PopLDdecay/bin/PopLDdecay -InVCF ../VCF/filtered/chr_${chr}_filtered_MAF001.vcf.gz ---OutStat ../results/LD/chr_${chr}.stat -MAF 0.05 -MaxDist 200

$HOME/soft/PopLDdecay/bin/PopLDdecay -InVCF ../VCF/filtered/chr_${chr}_filtered_MAF001.vcf.gz ---OutStat ../results/LD/chr_${chr}_Bugeater.stat -MAF 0.05 -MaxDist 200 -SubPop ../meta/Bugeater.txt

$HOME/soft/PopLDdecay/bin/PopLDdecay -InVCF ../VCF/filtered/chr_${chr}_filtered_MAF001.vcf.gz ---OutStat ../results/LD/chr_${chr}_mex.stat -MAF 0.05 -MaxDist 200 -SubPop ../meta/mex.txt

$HOME/soft/PopLDdecay/bin/PopLDdecay -InVCF ../VCF/filtered/chr_${chr}_filtered_MAF001.vcf.gz ---OutStat ../results/LD/chr_${chr}_parv.stat -MAF 0.05 -MaxDist 200 -SubPop ../meta/parviglumis.txt

$HOME/soft/PopLDdecay/bin/PopLDdecay -InVCF ../VCF/filtered/chr_${chr}_filtered_MAF001.vcf.gz ---OutStat ../results/LD/chr_${chr}_trop.stat -MAF 0.05 -MaxDist 200 -SubPop ../meta/trop.txt
