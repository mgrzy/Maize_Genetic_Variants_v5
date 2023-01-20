#!/bin/sh
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=16
#SBATCH --nodes=1
#SBATCH --time=5-00:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --job-name=teo2
#SBATCH --error=../job/teo2/job.%A_%a.err
#SBATCH --output=../job/teo2/job.%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mgrzybowski2@unl.edu
#SBATCH --array=1-33

module load anaconda
conda activate geno

samplesheet="../meta/teo2.txt"

fastqPath="../fastq/teo2/"
BAMPath="../BAM/teo2/"

repoPath="../reports/teo2/"

geno="../genome/maizeV5.fa"

samp=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
Taxa=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $2}'`

echo "${samp}"

fastp -i ${fastqPath}/${samp}_1.fastq.gz -I ${fastqPath}/${samp}_2.fastq.gz \
-o /scratch/${samp}_R1.fq.gz -O /scratch/${samp}_R2.fq.gz \
-w 16 -h ${repoPath}/${samp}.html -j ${repoPath}/${samp}.json

../../../speedseq/bin/speedseq align \
-o ${BAMPath}/${samp} \
-R $(echo "@RG\tID:${Taxa}\tSM:${Taxa}\tLB:libX\tPL:ILLUMINA\tPU=unitX") \
-t 16 \
-T ${Taxa} \
${geno} \
/scratch/${samp}_R1.fq.gz \
/scratch/${samp}_R2.fq.gz 
