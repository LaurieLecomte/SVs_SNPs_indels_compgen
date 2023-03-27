#!/bin/bash

# manitou
# srun -p small -c 1 -J indels_Fst_distrib -o log/indels_Fst_distrib_%j.log /bin/sh 01_scripts/utils/indels_Fst_distrib.sh &

# valeria
# srun -p ibis_small -c 1 -J indels_Fst_distrib -o log/indels_Fst_distrib_%j.log /bin/sh 01_scripts/utils/indels_Fst_distrib.sh &

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

POPS_FILE="02_infos/pops.txt"

FILT_ANGSD_VCF="$ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $INDELS_VCF_ANGSD)".maf"$MIN_MAF".vcf.gz"

REALSFS_PATH="/prg/angsd/0.937/misc/realSFS"


POP1='RO'
POP2='PU'

ANGSD_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $FILT_ANGSD_VCF)".indelsFst_"$POP1"_"$POP2".vcf.gz" # VCF formatted for angsd, with added Fst values from previous script
RAW_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_INDELS_VCF)".indelsFst_"$POP1"_"$POP2".vcf.gz" # input VCF (NOT the one formatted for angsd), with added Fst values from previous script


ANNOT_TABLE="11_annotation/genome_annotation/"$(basename -s .tsv $GENOME_ANNOT)".bed"
INDEL_BED="$ANNOT_DIR/indels_"$POP1"_"$POP2".bed"

# LOAD REQUIRED MODULES
module load R/4.1
module load bcftools/1.13

# 1. Extract Fst value from VCF
bcftools query -f "%CHROM\t%POS\t%ID\t%FST_"$POP1"_"$POP2"\n" $RAW_FST_VCF > $ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".table

bcftools query -f "%CHROM\t%POS\t%ID\t%FST_"$POP1"_"$POP2"\n" $ANGSD_FST_VCF > $ANGSD_FST_DIR/"$(basename -s .vcf.gz $ANGSD_FST_VCF)".table

# 2. Plot Fst distribution and per site Fst
Rscript 01_scripts/utils/plot_Fst_distrib.R $ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".table 02_infos/OV_to_ssa.txt $ANGSD_FST_DIR/"$POP1"_"$POP2"/"$POP1"_"$POP2".indels.fst





