#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=4
#SBATCH --nodes=1
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=30G
#SBATCH --job-name=pi
#SBATCH --error=job/job.%A_%a.out
#SBATCH --output=job/job.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mgrzybowski2@unl.edu
#SBATCH --array=1-1

# 1. Call variant and invariant site
ml anaconda
conda activate gatk4
ml bcftools
ml tabix

samplesheet="reg.txt"

chr=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}'`
Start=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $2}'`
End=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $3}'`
Name=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $4}'`
chrN=`echo ${chr:3}`

cp -a ../DB${chrN}/${Name} /scratch

gatk --java-options "-Xmx115g" GenotypeGVCFs \
   -R ../genome/maize.fa \
   -V gendb:///scratch/${Name} \
   -L ${chr}:${Start}-${End} \
   -all-sites \
   --only-output-calls-starting-in-intervals \
   -O /scratch/Reg_${chr}_${Start}_${End}.vcf.gz \
   --tmp-dir /scratch

# 2. Select and filter SNPs

gatk SelectVariants \
-V /scratch/Reg_${chr}_${Start}_${End}.vcf.gz \
-select-type SNP \
-O /scratch/SNPs.vcf.gz

bcftools +setGT /scratch/SNPs.vcf.gz -- -t q -n . -i 'FMT/DP<2' |\
bcftools view -i 'F_MISSING<0.5 & InbreedingCoeff>0 & MAF>0.01 & INFO/DP<33550 & INFO/DP>1515' -m2 -M2 -Ov |\
bcftools annotate --set-id +'%CHROM\_%POS' -O z -o /scratch/SNPraw.vcf.gz

tabix /scratch/SNPraw.vcf.gz

gatk VariantFiltration \
-V /scratch/SNPraw.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O /scratch/SNPsFiltered.vcf.gz

bcftools view /scratch/SNPsFiltered.vcf.gz -f PASS -Oz -o /scratch/SNPs.vcf.gz
tabix /scratch/SNPs.vcf.gz

# 3. Select and filter invariant site

bcftools +setGT /scratch/Reg_${chr}_${Start}_${End}.vcf.gz -- -t q -n . -i 'FMT/DP<2' |\
bcftools view -i 'N_ALT=0 & F_MISSING<0.5 & INFO/DP<33550 & INFO/DP>1515' -Ov |\
bcftools annotate --set-id +'%CHROM\_%POS' -O z -o /scratch/invariantRaw.vcf.gz

tabix /scratch/invariantRaw.vcf.gz

gatk VariantFiltration \
-V /scratch/invariantRaw.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O /scratch/InvariantFiltered.vcf.gz

bcftools view /scratch/InvariantFiltered.vcf.gz -f PASS -Oz -o /scratch/Invariant.vcf.gz
tabix /scratch/Invariant.vcf.gz

# 4. Merge and sort variant and invariant 

ml picard

picard -Xmx200g SortVcf TMP_DIR=/scratch/ \
I=/scratch/SNPs.vcf.gz \
I=/scratch/Invariant.vcf.gz \
O=/scratch/Variants.vcf.gz

#5. Run Pixy

conda activate pixy

pixy --stats pi fst dxy \
--vcf /scratch/Variants.vcf.gz \
--populations pop.txt \
--window_size 5000 \
--n_cores 4 \
--output_folder results \
--output_prefix Reg_${chr}_${Start}_${End} \
--chromosomes ${chr} \
--interval_start ${Start} \ 
--interval_end ${End}
