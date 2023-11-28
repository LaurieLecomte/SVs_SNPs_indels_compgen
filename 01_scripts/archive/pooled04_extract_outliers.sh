#!/bin/bash

# POOLED VARIANTS ANALYSIS : detect outliers and candidates on pooled variants set
# Run RDA on pooled set

# srun -p small -c1 -J pooled02_rda -o log/pooled02_rda_%j.log /bin/sh 01_scripts/utils/pooled02_rda.sh &


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
QUANTILE=0.97

# LOAD REQUIRED MODULES
module load vcftools/0.1.16
module load bcftools/1.13
module load bedtools/2.30.0
module load R/4.1


# 1. Fst outliers 
# Compute quantile cutoff
#Rscript 01_scripts/utils/compute_quantile_Fst_cutoff.R $ANGSD_FST_DIR/"$(basename -s .vcf.gz $ANGSD_FST_VCF)".table $QUANTILE
Rscript 01_scripts/utils/compute_quantile_Fst_cutoff.R $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5.table $QUANTILE 

MIN_FST="$(less $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5_quantile"$QUANTILE"_Fst_cutoff.txt" | head -n1)"

#less $ANGSD_FST_DIR/"$(basename -s .vcf.gz $ANGSD_FST_VCF)".table | awk -v val="$MIN_FST" 'BEGIN{FS="\t"} $5 >= val {print}' > $ANGSD_FST_DIR/SVs_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST".table
#echo "$(less $ANGSD_FST_DIR/SVs_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST".table | wc -l) outlier SVs with Fst >= $MIN_FST"

# Extract Fst outliers
less $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5.table | awk -v val="$MIN_FST" 'BEGIN{FS="\t"} $5 >= val {print}' > $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5_minFst"$MIN_FST".table

#  2. Fisher outliers


# Intersection outliers set
# 1. Find shared outliers between Fst and Fisher and candidates uniques to RDA
Rscript 01_scripts/utils/pooled_compare_outliers_candidates_sites.R $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5_minFst"$MIN_FST".table $POOLED_DIR/RDA_3sd_outliers.table $POOLED_DIR/SVs_SNPs_indels_fisher_outliers_qval"$MAX_QVAL".table $OVERLAP_WIN $QUANTILE $MIN_FST $SD $MAX_QVAL $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared.table $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5_RDA_"$SD"sd_outliers_uniques.table $POOLED_DIR/SVs_SNPs_indels_MAF0.05_FMISS0.5_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_RDA_"$SD"sd_shared.table
