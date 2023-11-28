#!/bin/bash

# POOLED VARIANTS ANALYSIS : detect outliers and candidates on pooled variants set
# Run RDA on pooled set

# srun -p medium -c1 --mem=50G --time=3-00:00:00 -J pooled03_rda -o log/pooled03_rda_%j.log /bin/sh 01_scripts/utils/pooled03_rda.sh &


# VARIABLES
SV_FST_VCF="08_angsd_fst/SVs/merged_SUPP2_MAF0.05_FMISS0.5.SVsFst_RO_PU.vcf.gz"
SNP_FST_VCF="08_angsd_fst/SNPs/SNPs_MAF0.05_FMISS0.5.SNPsFst_RO_PU.vcf.gz"
INDEL_FST_VCF="08_angsd_fst/indels/indels_MAF0.05_FMISS0.5.indelsFst_RO_PU.vcf.gz"

POP1='RO'
POP2='PU'

POOLED_DIR="pooled_analysis"
ID_SEX_POP="02_infos/ID_Sex_Pop_updated.txt"

MAX_QVAL=0.01
SD=3

# LOAD REQUIRED MODULES
module load vcftools/0.1.16
module load bcftools/1.13
module load bedtools/2.30.0
module load R/4.1


# 1. Run RDA
Rscript 01_scripts/utils/pooled_rda.R $POOLED_DIR/"$(basename -s .vcf.gz $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5.vcf.gz)".geno_mat.012 $ID_SEX_POP $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5_CHROM_POS_END_ID $SD $POOLED_DIR