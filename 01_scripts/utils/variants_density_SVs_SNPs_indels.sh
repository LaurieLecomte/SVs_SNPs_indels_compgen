#!/bin/bash

# Plot variant density, with SVs, SNPs and indels together

# srun -p small -c 1 -J variants_p -o log/variants_prop_%j.log /bin/sh 01_scripts/utils/variants_prop.sh &

# VARIABLES
GENOME="03_genome/genome.fasta"
VCF_DIR="04_vcf"

DENSITY_DIR="density"
WIN_SIZE=$1
INDELS_DENSITY="$DENSITY_DIR/indels/indels_density_win"$WIN_SIZE"bp.txt"
SV_DENSITY="$DENSITY_DIR/SVs/SVs_density_win"$WIN_SIZE"bp.txt"
SNP_DENSITY="$DENSITY_DIR/SNPs/SNPs_density_win"$WIN_SIZE"bp.txt"

CHR_CONV="02_infos/OV_to_ssa.txt"


#ANGSD_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $FILT_ANGSD_VCF)".SVsFst_"$POP1"_"$POP2".vcf.gz"


# LOAD REQUIRED MODULES
module load bcftools/1.13
module load R/4.1

