#!/bin/bash

# Format input indels VCF for running ANGDS

# manitou
# srun -p small -c 1 --mem=50G -J 03.1_indels_format_VCF -o log/03.1_indels_format_VCF_%j.log /bin/sh 01_scripts/03.1_indels_format_VCF.sh &

# valeria
# srun -p ibis_small -c 1 --mem=50G -J 03.1_indels_format_VCF -o log/03.1_indels_format_VCF_%j.log /bin/sh 01_scripts/03.1_indels_format_VCF.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"

RAW_VCF_DIR="04_vcf/indels"
RAW_INDEL_VCF="$RAW_VCF_DIR/indels_MAF0.05_FMISS0.5.vcf.gz"

ANGSD_INPUT_DIR="05_angsd_inputs/indels"
ANGSD_STATS_DIR="06_angsd_stats/indels"
ANGSD_BYPOP_DIR="07_angsd_bypop/indels"
ANGSD_FST_DIR="08_angsd_fst/indels"

PCA_DIR="09_pca/indels"
RDA_DIR="10_rda/indels"
ANNOT_DIR="11_annotation/indels"
GO_DIR="12_go/indels"

GENOME_ANNOT="11_annotation/genome_annotation/genome_annotation_table_simplified_1.5.tsv"
ID_SEX_POP="02_infos/ID_Sex_Pop_updated.txt"

# LOAD REQUIRED MODULES
module load bcftools/1.13
module load R/4.1

# 1. Recode REF and ALT as dummy SNPs alleles (ANGSD does not work on non-SNPs variants)
Rscript 01_scripts/utils/format_SVs_indels.R $RAW_INDEL_VCF $GENOME $ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_INDEL_VCF)".recoded.tmp

# 2. Remove tag INFO/INDEL, sort and index
bcftools annotate -x INFO/INDEL $ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_INDEL_VCF)".recoded.tmp | bcftools sort -Oz > $ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_INDEL_VCF)".recoded.vcf.gz
tabix -p vcf $ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_INDEL_VCF)".recoded.vcf.gz -f
