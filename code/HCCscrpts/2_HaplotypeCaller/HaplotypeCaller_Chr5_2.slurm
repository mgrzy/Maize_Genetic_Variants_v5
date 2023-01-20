#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=4
#SBATCH --nodes=1
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=chr5
#SBATCH --error=../job/gatk/chr5/job.%A_%a.err
#SBATCH --output=../job/gatk/chr5/job.%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mgrzybowski2@unl.edu
#SBATCH --array=1-3

ml gatk4

samplesheet="../meta/BAM2.txt"

BAM=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
samp=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $2}'`
chr="5"

echo "${samp}"

gatk --java-options "-Xmx16g" HaplotypeCaller  \
   -R ../genome/maize.fa \
   -I ${BAM} \
   -O /scratch/${samp}_chr${chr}.g.vcf.gz \
   -L chr${chr} \
   -ERC GVCF

mv /scratch/${samp}_chr${chr}.g.vcf.gz* ../gVCF/chr${chr}
