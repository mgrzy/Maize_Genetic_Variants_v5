#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=2
#SBATCH --nodes=1
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --job-name=jChr5
#SBATCH --error=../job/joint/chr5/chr5.%A_%a.out
#SBATCH --output=../job/joint/chr5/chr5.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mgrzybowski2@unl.edu
#SBATCH --array=1-226

ml anaconda
conda activate gatk4

samplesheet="../meta/ChrReg1Mb/chr5.txt"

chr=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
Start=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $2}'`
End=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $3}'`
Name=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $4}'`

cp -a ../DB/${Name} /scratch

gatk --java-options "-Xmx35g" GenotypeGVCFs \
   -R ../genome/maize.fa \
   -V gendb:///scratch/${Name} \
   -L ${chr}:${Start}-${End} \
   --only-output-calls-starting-in-intervals \
   -O /scratch/Reg_${chr}_${Start}_${End}.vcf.gz \
   --tmp-dir /scratch

mv /scratch/Reg_${chr}_${Start}_${End}.vcf.gz* ../VCF/raw/chr5
