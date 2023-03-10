#!/bin/bash

# Compute SAF and MAF

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

SV_VCF_ANGSD="$ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".recoded.vcf"

N_IND="$(less $ID_SEX_POP | wc -l)"

CPU=4

REGION=$1

# LOAD REQUIRED MODULES
module load angsd/0.937
module load bcftools/1.15

# 1. Calculate the SAF, MAF by regions
angsd -vcf-Gl $SV_VCF_ANGSD -nind $N_IND -P $CPU \
-fai "$GENOME".fai -anc $GENOME -ref $GENOME \
-domaf 1 -dosaf 1 -doMajorMinor 4 -r $REGION \
-out $ANGSD_STATS_DIR/"(basename -s .vcf $SV_VCF_ANGSD)"_chr"$REGION"


# Merge SAF files

