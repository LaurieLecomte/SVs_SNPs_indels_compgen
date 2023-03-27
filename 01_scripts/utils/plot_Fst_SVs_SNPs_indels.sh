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

MAX_MISS=0.5

POPS_FILE="02_infos/pops.txt"

FILT_ANGSD_VCF="$ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $INDELS_VCF_ANGSD)".maf"$MIN_MAF".vcf.gz"

REALSFS_PATH="/prg/angsd/0.937/misc/realSFS"


POP1='RO'
POP2='PU'

ANGSD_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $FILT_ANGSD_VCF)".indelsFst_"$POP1"_"$POP2".vcf.gz" # VCF formatted for angsd, with added Fst values from previous script
RAW_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_INDELS_VCF)".indelsFst_"$POP1"_"$POP2".vcf.gz" # input VCF (NOT the one formatted for angsd), with added Fst values from previous script


ANNOT_TABLE="11_annotation/genome_annotation/"$(basename -s .tsv $GENOME_ANNOT)".bed"

SV_BED="11_annotation/SVs/SVs_"$POP1"_"$POP2".bed"
SNP_BED="11_annotation/SNPs/indels_"$POP1"_"$POP2".bed"
INDEL_BED="11_annotation/indels/indels_"$POP1"_"$POP2".bed"

SV_FST_TABLE="08_angsd_fst/SVs/merged_SUPP2_MAF"$MIN_MAF"_FMISS"$MAX_MISS".SVsFst_"$POP1"_"$POP2".table"
SNP_FST_TABLE="08_angsd_fst/SNPs/SNPs_MAF"$MIN_MAF"_FMISS"$MAX_MISS".SNPsFst_"$POP1"_"$POP2".table"
INDEL_FST_TABLE="08_angsd_fst/indels/indels_MAF"$MIN_MAF"_FMISS"$MAX_MISS".indelsFst_"$POP1"_"$POP2".table"

SV_GW_FST="08_angsd_fst/SVs/"$POP1"_"$POP2"/"$POP1"_"$POP2".SVs.fst"
SNP_GW_FST="08_angsd_fst/SNPs/"$POP1"_"$POP2"/"$POP1"_"$POP2".SNPs.fst"
SV_GW_FST="08_angsd_fst/indels/"$POP1"_"$POP2"/"$POP1"_"$POP2".indels.fst"


# 1. Plot Fst distribution and per site Fst for the 3 variant types

Rscript 01_scripts/utils/plot_Fst_SVs_SNPs_indels.R $SV_FST_TABLE $SNP_FST_TABLE $INDEL_FST_TABLE $SV_GW_FST $SNP_GW_FST $INDEL_GW_FST 02_infos/OV_to_ssa.txt
