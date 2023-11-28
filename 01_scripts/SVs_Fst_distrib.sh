#!/bin/bash

# 

# manitou
# srun -p small -c 1 -J SVs_Fst_distrib -o log/SVs_Fst_distrib_%j.log /bin/sh 01_scripts/utils/SVs_Fst_distrib.sh &

# valeria
# srun -p ibis_small -c 1 -J SVs_Fst_distrib -o log/SVs_Fst_distrib_%j.log /bin/sh 01_scripts/utils/SVs_Fst_distrib.sh &

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

CPU=4

MIN_MAF=0.05
MAX_MAF=0.95

FILT_ANGSD_VCF="$ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".vcf.gz"


POP1='RO'
POP2='PU'

#ANGSD_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $FILT_ANGSD_VCF)".SVsFst_"$POP1"_"$POP2".vcf.gz" # VCF formatted for angsd, with added Fst values from previous script
RAW_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".SVsFst_"$POP1"_"$POP2".vcf.gz" # input VCF (NOT the one formatted for angsd), with added Fst values from previous script

OVERLAP_WIN=$1
MIN_FST=$2

ANNOT_TABLE="11_annotation/genome_annotation/"$(basename -s .tsv $GENOME_ANNOT)".bed"
SV_BED="$ANNOT_DIR/SVs_"$POP1"_"$POP2".bed"

FST_WIN=100000
FST_STEP=10000

# LOAD REQUIRED MODULES
module load R/4.1
module load bcftools/1.13

# 1. Extract Fst value from VCF
bcftools query -f "%CHROM\t%POS\t%END\t%ID\t%FST_"$POP1"_"$POP2"\n" $RAW_FST_VCF > $ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".table


# 2. Plot Fst distribution and per site Fst
Rscript 01_scripts/utils/plot_Fst_distrib.R $ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".table 02_infos/OV_to_ssa.txt $ANGSD_FST_DIR/"$POP1"_"$POP2"/"$POP1"_"$POP2".SVs.fst $ANGSD_FST_DIR/"$POP1"_"$POP2"/"$POP1"_"$POP2"_win"$FST_WIN"_step"$FST_STEP".txt $FST_WIN


