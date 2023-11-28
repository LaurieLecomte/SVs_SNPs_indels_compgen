#!/bin/bash

# Compute SAF and MAF by CHROMOSOME

# manitou
# srun -p small -c 1 -J 01.2.1_SVs_saf_maf_all_merge -o log/01.2.1_SVs_saf_maf_all_merge_%j.log /bin/sh 01_scripts/01.2.1_SVs_saf_maf_all_merge.sh &

# valeria
# srun -p ibis_small -c 1 -J 01.2.1_SVs_saf_maf_all_merge -o log/01.2.1_SVs_saf_maf_all_merge_%j.log /bin/sh 01_scripts/01.2.1_SVs_saf_maf_all_merge.sh &

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

CPU=2

REGION=$1

