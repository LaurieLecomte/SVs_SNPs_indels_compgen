#!/bin/bash

# Format input SVs VCF for running ANGDS

# manitou
# srun -p small -c 1 -J 01.1_SVs_format_VCF -o log/01.1_SVs_format_VCF_%j.log /bin/sh 01_scripts/01.1_SVs_format_VCF.sh &

# valeria
# srun -p ibis_small -c 1 -J 01.1_SVs_format_VCF -o log/01.1_SVs_format_VCF_%j.log /bin/sh 01_scripts/01.1_SVs_format_VCF.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"

RAW_VCF_DIR="04_vcf/SVs"
RAW_SV_VCF="$RAW_VCF_DIR/merged_SUPP2_MAF0.05_FMISS0.5.vcf.gz"

ANGSD_INPUT_DIR="05_angsd_inputs/SVs"
ANGSD_STATS_DIR="06_angsd_stats/SVs"
ANGSD_BYPOP_DIR="07_angsd_bypop/SVs"
ANGSD_FST_DIR="08_angsd_fst/SVs"

PCA_DIR="09_pca/SVs"
RDA_DIR="10_rda/SVs"
ANNOT_DIR="11_annotation/SVs"
GO_DIR="12_go/SVs"

GENOME_ANNOT="11_annotation/genome_annotation/genome_annotation_table_simplified_1.5.tsv"
ID_SEX_POP="02_infos/ID_Sex_Pop_updated.txt"

# LOAD REQUIRED MODULES
module load bcftools/1.13
module load R/4.1

# 1. Convert GL tags to PL : otherwise, GL values are too small and likely impossible to compare with SNPs/indels
bcftools +tag2tag $RAW_SV_VCF -- -r --gl-to-pl > $RAW_VCF_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)"_PL.vcf

# 2. Recode REF and ALT as dummy SNPs alleles (ANGSD does not work on non-SNPs variants)
Rscript 01_scripts/utils/format_SVs_indels.R $RAW_VCF_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)"_PL.vcf $GENOME $ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".recoded.vcf

