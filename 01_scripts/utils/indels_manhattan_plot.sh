#!/bin/bash

# Produce Manhattan plot for Fst, with info on candidate variation sites

# Specify window size at positional arg 1, Fst quantile threshold at positional arg 2, number of SDs used for RDA at positional arg 3, and max qvalue used for Fisher tests at arg 4

# Works on ONE population pair at the time, so variables ANGSD_FST_VCF, RAW_FST_VCF and POP_PAIR must be adjusted accordingly - I only have 2 populations (RO and PU), so VCF names will be written as is

# manitou
# srun -p small -c 1 -J indels_manhattan_plot -o log/indels_manhattan_plot_%j.log /bin/sh 01_scripts/utils/indels_manhattan_plot.sh 10000 0.97 3 0.01 &

# valeria
# srun -p ibis_small -c 1 -J indels_manhattan_plot -o log/indels_manhattan_plot_%j.log /bin/sh 01_scripts/utils/indels_manhattan_plot.sh 10000 0.97 3 0.01 &

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
#ANNOT_DIR="11_annotation/indels"
GO_DIR="12_go/indels"
FISHER_DIR="11_fisher_tests/indels"

ID_SEX_POP="02_infos/ID_Sex_Pop_updated.txt"

INDELS_VCF_ANGSD="$ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_INDELS_VCF)".recoded.vcf.gz"

N_IND="$(less $ID_SEX_POP | wc -l)"


#MIN_MAF=0.05
#MAX_MAF=0.95

POPS_FILE="02_infos/pops.txt"

#FILT_ANGSD_VCF="$ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $INDELS_VCF_ANGSD)".maf"$MIN_MAF".vcf.gz"

POP1='RO'
POP2='PU'

ANGSD_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .recoded.vcf.gz $INDELS_VCF_ANGSD)".indelsFst_"$POP1"_"$POP2".vcf.gz" # VCF formatted for angsd, with added Fst values from previous script
RAW_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_INDELS_VCF)".indelsFst_"$POP1"_"$POP2".vcf.gz" # input VCF (NOT the one formatted for angsd), with added Fst values from previous script

GENOME_ANNOT="03_genome/annotation/genome_annotation_table_simplified_1.5.tsv"
ANNOT_TABLE="03_genome/annotation/"$(basename -s .tsv $GENOME_ANNOT)".table"

OVERLAP_WIN=$1
QUANTILE=$2
SD=$3
MAX_QVAL=$4

MIN_FST="$(less "$ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)"_quantile"$QUANTILE"_Fst_cutoff.txt" | head -n1)"
FST_OUTLIERS="$ANGSD_FST_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST".table"

RDA_OUTLIERS="$RDA_DIR/RDA_"$SD"sd_outliers.table"

FISHER_OUTLIERS="$FISHER_DIR/indels_fisher_"$POP1"_"$POP2"_outliers_qval"$MAX_QVAL".table"

GO_DB="12_go/go_db/go-basic.obo"
GO_ANNOT="12_go/go_db/all_go_annotations.csv"

# LOAD REQUIRED MODULES
module load bedtools/2.30.0
module load bcftools/1.13
module load R/4.1

ALL_SITES="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_INDELS_VCF)".indelsFst_"$POP1"_"$POP2".table"
INTERSECT_OUTLIERS="$ANGSD_FST_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared.table"
SHARED_3_OUTLIERS="$ANGSD_FST_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_RDA_"$SD"sd_shared.table"
OV_TO_SSA="02_infos/OV_to_ssa.txt"

Rscript 01_scripts/utils/manhattan_plot_candidates.R $ALL_SITES $FST_OUTLIERS $RDA_OUTLIERS $FISHER_OUTLIERS $INTERSECT_OUTLIERS $SHARED_3_OUTLIERS $QUANTILE $MIN_FST  $SD $MAX_QVAL $OVERLAP_WIN $OV_TO_SSA