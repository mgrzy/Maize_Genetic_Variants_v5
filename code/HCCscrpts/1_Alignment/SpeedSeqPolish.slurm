#!/bin/sh
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=16
#SBATCH --nodes=1
#SBATCH --time=5-00:00:00
#SBATCH --mem-per-cpu=6G
#SBATCH --job-name=Pol
#SBATCH --error=../job/Pol/job.%A_%a.err
#SBATCH --output=../job/Pol/job.%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mgrzybowski2@unl.edu
#SBATCH --array=1-11%11

module load anaconda
conda activate geno

samplesheet="../meta/Pol.txt"

fastqPath="../fastq/Pol/"
BAMPath="../BAM/Pol/"

repoPath="../reports/Pol/"

geno="../genome/maizeV5.fa.gz"

samp=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
Taxa=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $2}'`

echo "${samp}"

fastp -i ${fastqPath}/${samp}_RE_R1.fq.gz -I ${fastqPath}/${samp}_RE_R2.fq.gz \
-o /scratch/${samp}_R1.fq.gz -O /scratch/${samp}_R2.fq.gz \
-w 16 -h ${repoPath}/${samp}.html -j ${repoPath}/${samp}.json

../../speedseq/bin/speedseq align \
-o ${BAMPath}/${samp} \
-R $(echo "@RG\tID:${Taxa}\tSM:${Taxa}\tLB:libX\tPL:ILLUMINA\tPU=unitX") \
-t 16 \
-T ${Taxa} \
${geno} \
/scratch/${samp}_R1.fq.gz \
/scratch/${samp}_R2.fq.gz 
