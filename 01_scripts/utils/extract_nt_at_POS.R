# Recode SVs as SNPs for running angsd

# 0. Access to files provided in command line arguments -------------------
argv <- commandArgs(T)
INPUT_VCF <- argv[1]          # input VCF
GENOME <- argv[2]         # reference genome fasta
OUTPUT_VCF <- argv[3]  # output VCF file

# 1. Source recoding function
source('01_scripts/utils/recode_as_SNP.R')

# 2. Process VCF
recode_as_SNPs(INPUT_VCF, GENOME, OUTPUT_VCF)


