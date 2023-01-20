#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=20
#SBATCH --nodes=1
#SBATCH --time=1-10:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH --job-name=impAll
#SBATCH --error=../job/imp/impAll.%A_%a.out
#SBATCH --output=../job/imp/impAll.%A_%a.err
#SBATCH --array=1-10

chr=$SLURM_ARRAY_TASK_ID

ml java/1.8

java -Xmx1000g -jar $HOME/soft/beagle/beagle.19Apr22.7c0.jar gt=../VCF/filtered/chr_${chr}_filtered_MAF001.vcf.gz out=../VCF/imputed/chr_${chr}_imputed nthreads=20

ml tabix
tabix ../VCF/imputed/chr_${chr}_imputed.vcf.gz
