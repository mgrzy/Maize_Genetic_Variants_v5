#!/bin/sh
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=8
#SBATCH --nodes=1
#SBATCH --time=5-00:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --job-name=south4
#SBATCH --error=../job/aligment/South.%A_%a.err
#SBATCH --output=../job/aligment/South.%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mgrzybowski2@unl.edu
#SBATCH --array=2,4

module load anaconda
conda activate geno

samplesheet="../meta/SampleLists/South.txt"

fastqPath="../fastq/south/"
BAMPath="../BAM/south/"

repoPath="../reports/south/"

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
