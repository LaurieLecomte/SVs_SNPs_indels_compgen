#!/bin/bash

# Compute LD between SNPs only

# srun -p small -c 3 --mem=50G -J LD_SNPs_only -o log/LD_SNPs_only_%j.log /bin/sh 01_scripts/utils/LD_SNPs_only.sh &


# VARIABLES
VCF_DIR="04_vcf"
SV_VCF="$VCF_DIR/SVs/merged_SUPP2_MAF0.05_FMISS0.5.vcf.gz"
SNP_VCF="$VCF_DIR/SNPs/SNPs_MAF0.05_FMISS0.5.vcf.gz"

LD_DIR="LD/SNPs_vs_SNPs"

MAX_DIST=100000
WIN_KB=100

CPU=3

# LOAD REQUIRED MODULES
module load bcftools/1.13
module load python/3.7
module load plink
module load bedtools

if [[ ! -d $LD_DIR ]]
then
  mkdir $LD_DIR
fi

CHR="OV354430.1"

# Compute LD using VCF of SNPs that DO NOT overlap with SVs, produced in the LD_SVs_SNPs_plink.sh script
plink --r2 --allow-extra-chr --chr $CHR --ld-window 100 --ld-window-kb $WIN_KB --ld-window-r2 0 --vcf LD/SVs_vs_SNPs/non_overlapping_SNPs_10bp.vcf.gz --threads $CPU --out $LD_DIR/SNPs_SNPs_no_overlap_"$WIN_KB"_$CHR