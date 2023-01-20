#!/bin/sh
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=8
#SBATCH --nodes=1
#SBATCH --time=5-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=EU
#SBATCH --error=../job/aligment/eu.%A_%a.err
#SBATCH --output=../job/aligment/eu.%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mgrzybowski2@unl.edu
#SBATCH --array=1-71%71

module load anaconda
conda activate geno

samplesheet="../meta/SampleLists/Eu.txt"

fastqPath="../fastq/eu/"
BAMPath="../BAM/eu/"

repoPath="../reports/eu/"

geno="../genome/maizeV5.fa.gz"

samp=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
Taxa=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $2}'`

echo "${samp}"

fastp -i ${fastqPath}/${samp}_1.fastq.gz -I ${fastqPath}/${samp}_2.fastq.gz \
-o /scratch/${samp}_R1.fq.gz -O /scratch/${samp}_R2.fq.gz \
-w 16 -h ${repoPath}/${samp}.html -j ${repoPath}/${samp}.json

../../speedseq/bin/speedseq align \
-o ${BAMPath}/${Taxa} \
-R $(echo "@RG\tID:${Taxa}\tSM:${Taxa}\tLB:libX\tPL:ILLUMINA\tPU=unitX") \
-t 16 \
-T /scratch/${Taxa} \
${geno} \
/scratch/${samp}_R1.fq.gz \
/scratch/${samp}_R2.fq.gz 

rm ${BAMPath}/${Taxa}.discordants.bam
rm ${BAMPath}/${Taxa}.discordants.bam.bai
rm ${BAMPath}/${Taxa}.splitters.bam
rm ${BAMPath}/${Taxa}.splitters.bam.bai
