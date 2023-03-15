#!/bin/bash

# Convert VCF to beagle and normalize genotype likelihoods

# manitou
# srun -p small -c 4 -J 03.4_indels_pca -o log/03.4_indels_pca_%j.log /bin/sh 01_scripts/03.4_indels_pca.sh &

# valeria
# srun -p ibis_small -c 4 -J 03.4_indels_pca -o log/03.4_indels_pca_%j.log /bin/sh 01_scripts/03.4_indels_pca.sh &

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

MIN_MAF=0.05
MAX_MAF=0.95

CHR_LIST="02_infos/chr_list.txt"

INDELS_BEAGLE="$ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $INDELS_VCF_ANGSD)".maf"$MIN_MAF".norm.beagle.gz"
COV_MAT="$PCA_DIR/"$(basename -s .norm.beagle.gz $INDELS_BEAGLE)".cov"


# LOAD REQUIRED MODULES
module load bcftools/1.13
module load R/4.1
module load python/3.7
module load pcangsd/1.10

 
# 1. Analyse covariance matrix
pcangsd --threads $CPU --beagle $INDELS_BEAGLE --out $COV_MAT

# 2. Run PCA in R
Rscript 01_scripts/utils/pca_simple.R "$COV_MAT" $INDELS_BEAGLE $ID_SEX_POP