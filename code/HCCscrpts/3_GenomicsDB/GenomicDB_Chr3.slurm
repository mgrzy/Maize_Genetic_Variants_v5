#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=20
#SBATCH --nodes=1
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --job-name=DB3
#SBATCH --error=../job/DB/DB3.%A_%a.out
#SBATCH --output=../job/DB/DB3.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mgrzybowski2@unl.edu
#SBATCH --array=13-17

module load anaconda
conda activate gatk42

T=$SLURM_NTASKS_PER_NODE # n CPU

samplesheet="../meta/ChrReg5Mb.txt"

Chr=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
Start=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $2}'`
End=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $3}'`
Name=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $4}'`

meta="../meta/gVCFPath"

gatk --java-options "-Xmx80g" GenomicsDBImport \
   --genomicsdb-workspace-path /scratch/${Name} \
   --batch-size 50 \
   --sample-name-map ${meta}/AllLines${Chr}.txt \
   --tmp-dir /scratch \
   -L ${Chr}:${Start}-${End} \
   --genomicsdb-shared-posixfs-optimizations true\
   --max-num-intervals-to-import-in-parallel 2 \
   --reader-threads ${T}

mv /scratch/${Name} ../DB/
