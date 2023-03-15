#!/bin/bash

# Convert VCF to beagle and normalize genotype likelihoods

# manitou
# srun -p small -c 1 -J 02.2_SNPs_vcf_to_beagle -o log/02.2_SNPs_vcf_to_beagle_%j.log /bin/sh 01_scripts/02.2_SNPs_vcf_to_beagle.sh &

# valeria
# srun -p ibis_small -c 1 -J 02.2_SNPs_vcf_to_beagle -o log/02.2_SNPs_vcf_to_beagle_%j.log /bin/sh 01_scripts/02.2_SNPs_vcf_to_beagle.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"

RAW_VCF_DIR="04_vcf/SNPs"
RAW_SNPS_VCF="$RAW_VCF_DIR/SNPs_MAF0.05_FMISS0.5.vcf.gz"

ANGSD_INPUT_DIR="05_angsd_inputs/SNPs"
ANGSD_STATS_DIR="06_angsd_stats/SNPs"
ANGSD_BYPOP_DIR="07_angsd_bypop/SNPs"
ANGSD_FST_DIR="08_angsd_fst/SNPs"

PCA_DIR="09_pca/SNPs"
RDA_DIR="10_rda/SNPs"
ANNOT_DIR="11_annotation/SNPs"
GO_DIR="12_go/SNPs"

GENOME_ANNOT="11_annotation/genome_annotation/genome_annotation_table_simplified_1.5.tsv"
ID_SEX_POP="02_infos/ID_Sex_Pop_updated.txt"

SNPS_VCF_ANGSD=$RAW_SNPS_VCF

N_IND="$(less $ID_SEX_POP | wc -l)"

CPU=1

MIN_MAF=0.05
MAX_MAF=0.95

CHR_LIST="02_infos/chr_list.txt"

SNPS_BEAGLE="$ANGSD_STATS_DIR/"$(basename -s .vcf.gz $SNPS_VCF_ANGSD)".maf"$MIN_MAF".norm.beagle.gz"
COV_MAT="$PCA_DIR/"$(basename -s .norm.beagle.gz $SNPS_BEAGLE)".cov"

# LOAD REQUIRED MODULES
module load bcftools/1.13
module load R/4.1
module load python/3.7
module load pcangsd/1.10

 
# 1. Analyse covariance matrix
pcangsd --threads $CPU --beagle $SNPS_BEAGLE --out $COV_MAT

# 2. Run PCA in R
Rscript 01_scripts/utils/pca_simple.R $COV_MAT $SNPS_BEAGLE $ID_SEX_POP

