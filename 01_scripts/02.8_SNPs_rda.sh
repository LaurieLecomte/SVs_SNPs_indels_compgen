#!/bin/bash

# Perform RDA on SNPs with population as the only explanatory variable
# Specify window size at positional arg 1

# manitou
# srun -p medium --time=7-00:00:00 -c 2 --mem=50G -J 02.8_SNPs_rda -o log/02.8_SNPs_rda_%j.log /bin/sh 01_scripts/02.8_SNPs_rda.sh 10000 &

# valeria
# srun -p ibis_medium --time=7-00:00:00 -c 2 --mem=50G -J 02.8_SNPs_rda -o log/02.8_SNPs_rda_%j.log /bin/sh 01_scripts/02.8_SNPs_rda.sh 10000 &

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

CPU=2

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

OVERLAP_WIN=$1

# LOAD REQUIRED MODULES
module load vcftools/0.1.16
module load bcftools/1.13
module load R/4.1

# 1. Produce a list of SNPs sites and IDs to make interpretation easier after RDA
bcftools query -f '%CHROM\t%POS\t%END\t%ID\n' $RAW_FST_VCF > $RDA_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".CHR_POS_END_ID.table


# 2. Convert VCF to genotype matrix
vcftools --gzvcf $RAW_FST_VCF --012 --out $RDA_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".geno_mat
echo "imputation done"

# 3. Impute missing data
Rscript 01_scripts/utils/impute_missing.R $RDA_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".geno_mat.012 $ID_SEX_POP

# 4. Run RDA
Rscript 01_scripts/utils/rda.R $RDA_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".geno_mat.012 $ID_SEX_POP $RDA_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".CHR_POS_END_ID.table $RDA_DIR

# 5. Get overlap of outlier sites with known genes
tail -n+2 $RDA_DIR/RDA_outliers.txt > $RDA_DIR/SNPs_RDA_outliers.table
bedtools window -a $ANNOT_TABLE -b $RDA_DIR/SNPs_RDA_outliers.table -w $OVERLAP_WIN > $RDA_DIR/SNPs_"$POP1"_"$POP2"_outliers_RDA_overlap"$OVERLAP_WIN"bp.table