#!/bin/bash

# POOLED VARIANTS ANALYSIS : detect outliers and candidates on pooled variants set
# Run Fisher tests on pooled set

# srun -p small -c1 --mem=50G -J pooled02_fisher -o log/pooled02_fisher_%j.log /bin/sh 01_scripts/utils/pooled02_fisher.sh &


# VARIABLES
SV_FST_VCF="08_angsd_fst/SVs/merged_SUPP2_MAF0.05_FMISS0.5.SVsFst_RO_PU.vcf.gz"
SNP_FST_VCF="08_angsd_fst/SNPs/SNPs_MAF0.05_FMISS0.5.SNPsFst_RO_PU.vcf.gz"
INDEL_FST_VCF="08_angsd_fst/indels/indels_MAF0.05_FMISS0.5.indelsFst_RO_PU.vcf.gz"

POP1='RO'
POP2='PU'

POOLED_DIR="pooled_analysis"
ID_SEX_POP="02_infos/ID_Sex_Pop_updated.txt"

MAX_QVAL=0.01

# LOAD REQUIRED MODULES
module load vcftools/0.1.16
module load bcftools/1.13
module load bedtools/2.30.0
module load R/4.1

# 1. Extract POP and ID from ID_SEX_POP
cut -f1,3 $ID_SEX_POP > 02_infos/ID_POP.txt

# 2. Run Fisher tests
Rscript 01_scripts/utils/fisher_test.R $POOLED_DIR/"$(basename -s .vcf.gz $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5.vcf.gz)".geno_mat.012 02_infos/ID_POP.txt $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5_CHROM_POS_END_ID  $POP1 $POP2 $MAX_QVAL $POOLED_DIR/SVs_SNPs_indels_fisher

# Filter results
#Rscript 01_scripts/utils/filter_fisher.R $POOLED_DIR/"$(basename -s .vcf.gz $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5.vcf.gz)".geno_mat.012 02_infos/ID_POP.txt $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5.table  $POP1 $POP2 $MAX_QVAL $POOLED_DIR/SVs_SNPs_indels_fisher