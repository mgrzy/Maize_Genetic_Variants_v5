#!/bin/sh
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=10
#SBATCH --nodes=1
#SBATCH --time=5-00:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --job-name=Nebraska
#SBATCH --error=../job/Nebraska2/job.%A_%a.err
#SBATCH --output=../job/Nebraska2/job.%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mgrzybowski2@unl.edu
#SBATCH --array=1-231

module load anaconda
conda activate geno

samplesheet="../meta/Nebraska2.txt"

fastqPath="../fastq/nebraska/set3/"
BAMPath="../BAM/nebraska/"

repoPath="../reports/nebraska/"

geno="../genome/maizeV5.fa.gz"

samp=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
Taxa=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $2}'`

echo "${samp}"

fastp -i ${fastqPath}/${samp}_1.fastq.gz -I ${fastqPath}/${samp}_2.fastq.gz \
-o /scratch/${samp}_R1.fq.gz -O /scratch/${samp}_R2.fq.gz \
-w 16 -h ${repoPath}/${samp}.html -j ${repoPath}/${samp}.json

../../speedseq/bin/speedseq align \
-o ${BAMPath}/${samp} \
-R $(echo "@RG\tID:${Taxa}\tSM:${Taxa}\tLB:libX\tPL:ILLUMINA\tPU=unitX") \
-t 16 \
-T /scratch/${Taxa} \
${geno} \
/scratch/${samp}_R1.fq.gz \
/scratch/${samp}_R2.fq.gz 

rm ${BAMPath}/${samp}.discordants.bam
rm ${BAMPath}/${samp}.discordants.bam.bai
rm ${BAMPath}/${samp}.splitters.bam
rm ${BAMPath}/${samp}.splitters.bam.bai
