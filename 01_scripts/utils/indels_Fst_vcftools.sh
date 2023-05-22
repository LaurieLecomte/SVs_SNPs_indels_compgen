#!/bin/bash

# Estimate genome-wide and per-site Fst by site, using vcftools for validating angsd Fst output

# manitou
# srun -p small -c 1 --mem=20G -J indels_Fst_vcftools -o log/indels_Fst_vcftools_%j.log /bin/sh 01_scripts/utils/indels_Fst_vcftools.sh &

# valeria
# srun -p ibis_small -c 1 --mem=20G -J indels_Fst_vcftools -o log/indels_Fst_vcftools_%j.log /bin/sh 01_scripts/utils/indels_Fst_vcftools.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"

RAW_VCF_DIR="04_vcf/indels"
RAW_INDELS_VCF="$RAW_VCF_DIR/indels_MAF0.05_FMISS0.5.vcf.gz"

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

INDELS_VCF_ANGSD="$ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_INDELS_VCF)".recoded.vcf.gz"

N_IND="$(less $ID_SEX_POP | wc -l)"

CPU=4

#MIN_MAF=0.05
#MAX_MAF=0.95

POPS_FILE="02_infos/pops.txt"

#FILT_ANGSD_VCF="$ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $INDELS_VCF_ANGSD)".maf"$MIN_MAF".vcf.gz"

REALSFS_PATH="/prg/angsd/0.937/misc/realSFS"

# VARIABLES FOR ANGSD
WINDOW=1000000
WIN_STEP=10000 

# LOAD REQUIRED MODULES
module load bcftools/1.13
module load vcftools/0.1.16


# Per site
vcftools --gzvcf $INDELS_VCF_ANGSD --weir-fst-pop 02_infos/RO_IDs.txt --weir-fst-pop 02_infos/PU_IDs.txt --out $ANGSD_FST_DIR/"$(basename -s .vcf.gz $INDELS_VCF_ANGSD)"_vcftools

# Per window
vcftools --gzvcf $INDELS_VCF_ANGSD --weir-fst-pop 02_infos/RO_IDs.txt --weir-fst-pop 02_infos/PU_IDs.txt --fst-window-size $WINDOW --fst-window-step $WIN_STEP --out $ANGSD_FST_DIR/"$(basename -s .vcf.gz $INDELS_VCF_ANGSD)"_vcftools_win
