#!/bin/bash

# Perform GO enrichment analysis on a set of genes located near (at 10000 bp) outlier SVs, given their high Fst (0.4)
# Specify window size at positional arg 1, and min Fst threshold at positional arg 2
# Launch in a conda env where goatools is installed

# manitou
# srun -p small -c 1 -J SVs_rda -o log/SVs_rda_%j.log /bin/sh 01_scripts/SVs_rda.sh &

# valeria
# srun -p ibis_small -c 1 -J SVs_rda -o log/SVs_rda_%j.log /bin/sh 01_scripts/SVs_rda.sh &

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

SV_VCF_ANGSD="$ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".recoded.vcf.gz"

N_IND="$(less $ID_SEX_POP | wc -l)"

MIN_MAF=0.05
MAX_MAF=0.95

FILT_ANGSD_VCF="$ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".vcf.gz"

POP1='RO'
POP2='PU'

ANGSD_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $FILT_ANGSD_VCF)".SVsFst_"$POP1"_"$POP2".vcf.gz" # VCF formatted for angsd, with added Fst values from previous script
RAW_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".SVsFst_"$POP1"_"$POP2".vcf.gz" # input VCF (NOT the one formatted for angsd), with added Fst values from previous script

OVERLAP_WIN=$1
MIN_FST=$2

ANNOT_TABLE="11_annotation/genome_annotation/"$(basename -s .tsv $GENOME_ANNOT)".bed"
SV_BED="$ANNOT_DIR/SVs_"$POP1"_"$POP2".bed"


# LOAD REQUIRED MODULES
module load vcftools/0.1.16
module load bcftools/1.13
module load R/4.1


# 1. Produce a list of SV sites and IDs to make interpretation easier after RDA
bcftools query -f '%CHROM\t%POS\t%ID\n' $ANGSD_FST_VCF > $RDA_DIR/"$(basename -s .vcf.gz $ANGSD_FST_VCF)".CHR_POS_ID.txt


# 2. Convert VCF to genotype matrix
vcftools --gzvcf $ANGSD_FST_VCF --012 --out $RDA_DIR/"$(basename -s .vcf.gz $ANGSD_FST_VCF)".geno_mat

# 3. Impute missing data
Rscript 01_scripts/utils/impute_missing.R $RDA_DIR/"$(basename -s .vcf.gz $ANGSD_FST_VCF)".geno_mat.012 $ID_SEX_POP

# 4. Run RDA
Rscript 01_scripts/utils/rda.R $RDA_DIR/"$(basename -s .vcf.gz $ANGSD_FST_VCF)".geno_mat.012 $ID_SEX_POP $RDA_DIR/"$(basename -s .vcf.gz $ANGSD_FST_VCF)".CHR_POS_ID.txt $RDA_DIR