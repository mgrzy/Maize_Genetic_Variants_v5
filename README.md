
#  **A Common Resequencing-Based Genetic Marker Dataset for Global Maize Diversity**
![Linux](https://img.shields.io/badge/Linux-FCC624?style=for-the-badge&logo=linux&logoColor=black)
![](https://img.shields.io/badge/shell_script-%23121011.svg?style=for-the-badge&logo=gnu-bash&logoColor=white)
![](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)

### Scripts from paper: https://doi.org/10.1111/tpj.16123

### Variant calling: 
 * *code/HCCscrpts/1_Alignment/* - codes for alignment with SpeedSeeq
 * *code/HCCscrpts/2_HaplotypeCaller/* - codes to run HaplotypeCaller (GATK v. v.4.2.0.0)
 * *code/HCCscrpts/3_GenomicsDB/* - codes for to create  GenomicsDB datastores (GATK v. v.4.2.6.1)
 * *code/HCCscrpts/4_JointCalling/* - codes joint variant calling and producing VCF files (GATK v. v.4.2.6.1)
 * *code/HCCscrpts/5_Filter_VCFs/* - codes for variant filtering with Hard-Filtering criteria
 * *code/HCCscrpts/6_Merge_and_sort_VCFs/* - codes merging and sorting VCFs
 * *code/HCCscrpts/7_Merge_and_sort_VCFs/* - codes subsetting high confidence variant set 
 * *code/HCCscrpts/8_Merge_and_sort_VCFs/* - codes for imputing confidence variant set

### Population analysis:
* *code/HCCscrpts/9_Population_analysis/* - codes for population analysis (MAF, LD, nucleotide diversity, variants number per region, PCA)
* *code/HCCscrpts/10_GWAS/* - codes for GWSA for days to silking

## Data availability:
The additional resequencing data generated as part of this project has been deposited in the European Nucleotide
Archive (ENA) under the study accession numbers: PRJEB56265, PRJEB56295, and PRJEB56320. Raw VCF files
for all 366 million variants identified in this study, imputed VCF files for the 46 million quality and minor allele
frequency filtered variations identified as part of this study and GATK GenomicsDBs files to enable new SNP calling
with additional populations have been deposited are available for download from www.maizegdb.org
