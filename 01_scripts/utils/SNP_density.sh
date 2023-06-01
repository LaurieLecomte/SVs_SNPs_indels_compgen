#!/bin/bash

# Compute SNP density for a given window size, specified at positional arg 1
# 02.8_SNPs_rda.sh is a prequisite to this script, which uses the raw 012 matrix and CHROM_POS_END_ID files produced when doing RDA

# manitou
# srun -p small -c 1 -J SNP_density -o log/SNP_density_%j.log /bin/sh 01_scripts/utils/SNP_density.sh 1000000 &

# valeria
# srun -p ibis_small -c 1 -J SNP_density -o log/SNP_density_%j.log /bin/sh 01_scripts/utils/SNP_density.sh 1000000 &

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
#ANNOT_DIR="11_annotation/SNPs"
GO_DIR="12_go/SNPs"

ID_SEX_POP="02_infos/ID_Sex_Pop_updated.txt"

SNPS_VCF_ANGSD=$RAW_SNPS_VCF

N_IND="$(less $ID_SEX_POP | wc -l)"


MIN_MAF=0.05
MAX_MAF=0.95

POPS_FILE="02_infos/pops.txt"

FILT_ANGSD_VCF="$ANGSD_STATS_DIR/"$(basename -s .vcf.gz $SNPS_VCF_ANGSD)".maf"$MIN_MAF".vcf.gz"

POP1='RO'
POP2='PU'

#ANGSD_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $FILT_ANGSD_VCF)".SNPsFst_"$POP1"_"$POP2".vcf.gz" # VCF formatted for angsd, with added Fst values from previous script
RAW_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_SNPS_VCF)".SNPsFst_"$POP1"_"$POP2".vcf.gz" # input VCF, with added Fst values from previous script

GENOME_ANNOT="03_genome/annotation/genome_annotation_table_simplified_1.5.tsv"
ANNOT_TABLE="03_genome/annotation/"$(basename -s .tsv $GENOME_ANNOT)".table"

WIN_SIZE=$1

CHR_BED="02_infos/chrs.bed"

DENSITY_DIR="density/SNPs"


# LOAD REQUIRED MODULES
module load vcftools/0.1.16
module load bcftools/1.13
module load bedtools/2.30.0
module load bedops/2.4.40
module load R/4.1


# 1. Remove header from CHROM_POS_END_ID file (previously used for RDA and Fisher tests)
#less $RDA_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".CHR_POS_END_ID.table | tail -n+2 > $RDA_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".CHR_POS_END_ID.bed

# 2. Split genome by window
bedops --chop $WIN_SIZE -x $CHR_BED > 02_infos/chrs_win"$WIN_SIZE".bed

# 3. Get number of overlap between fixed size windows and sites 
bedtools intersect -a 02_infos/chrs_win"$WIN_SIZE".bed -b $RDA_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".CHR_POS_END_ID.table -c > $DENSITY_DIR/SNPs_density_win"$WIN_SIZE"bp.txt

# 4. Compute mean density by window and plot
Rscript 01_scripts/utils/compute_plot_density.R $DENSITY_DIR/SNPs_density_win"$WIN_SIZE"bp.txt 02_infos/OV_to_ssa.txt $WIN_SIZE