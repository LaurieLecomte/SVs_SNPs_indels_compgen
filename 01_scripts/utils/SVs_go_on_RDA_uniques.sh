#!/bin/bash

# Perform GO enrichment analysis on outliers SVs genes identified by RDA but not shared with Fst outliers
# First get overlap between genotyped SVs and known genes and between outlier SVs, while allowing a 10000 bp window around genes
# Specify window size at positional arg 1 and number of sd at arg 2
# Launch in a conda env where goatools is installed

# Works on ONE population pair at the time, so variables ANGSD_FST_VCF, RAW_FST_VCF and POP_PAIR must be adjusted accordingly - I only have 2 populations (RO and PU), so VCF names will be written as is

# manitou
# srun -p small -c 1 -J SVs_go_on_RDA_uniques -o log/SVs_go_on_RDA_uniques_%j.log /bin/sh 01_scripts/utils/SVs_go_on_RDA_uniques.sh 10000 3 &

# valeria
# srun -p ibis_small -c 1 -J SVs_go_on_RDA_uniques -o log/SVs_go_on_RDA_uniques_%j.log /bin/sh 01_scripts/utils/SVs_go_on_RDA_uniques.sh 10000 3 &

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
#ANNOT_DIR="11_annotation/SVs"
GO_DIR="12_go/SVs"

ID_SEX_POP="02_infos/ID_Sex_Pop_updated.txt"

SV_VCF_ANGSD="$ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".recoded.vcf.gz"

N_IND="$(less $ID_SEX_POP | wc -l)"

#MIN_MAF=0.05
#MAX_MAF=0.95

#FILT_ANGSD_VCF="$ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".vcf.gz"
#FILT_ANGSD_VCF="$ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".vcf.gz"

POP1='RO'
POP2='PU'

ANGSD_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".SVsFst_"$POP1"_"$POP2".vcf.gz" # VCF formatted for angsd, with added Fst values from previous script
RAW_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".SVsFst_"$POP1"_"$POP2".vcf.gz" # input VCF (NOT the one formatted for angsd), with added Fst values from previous script

OVERLAP_WIN=$1
SD=$2

GENOME_ANNOT="03_genome/annotation/genome_annotation_table_simplified_1.5.tsv"
ANNOT_TABLE="03_genome/annotation/"$(basename -s .tsv $GENOME_ANNOT)".table"

GO_DB="12_go/go_db/go-basic.obo"
GO_ANNOT="12_go/go_db/all_go_annotations.csv"

FST_OUTLIERS="$GO_DIR/SVs_RO_PU_outliers_minFst0.304_overlap10000bp_outlierIDs.txt"
RDA_OUTLIERS="$GO_DIR/SVs_RO_PU_outliers_RDA_3sd_overlap10000bp_outlierIDs.txt"
FISHER_OUTLIERS="$GO_DIR/SVs_fisher_RO_PU_outliers_qval0.01_overlap10000bp_outlierIDs.txt"

# LOAD REQUIRED MODULES
module load R/4.1


# 1. Extract outliers uniques to RDA
Rscript 01_scripts/utils/compare_outliers.R $FST_OUTLIERS $RDA_OUTLIERS $FISHER_OUTLIERS

# 2. Perform GO enrichment analysis on filtered RDA outliers set
python 12_go/goatools/scripts/find_enrichment.py --pval=0.05 --indent \
  --obo $GO_DB \
  $GO_DIR/SVs_"$POP1"_"$POP2"_outliers_RDA_"$SD"sd_overlap"$OVERLAP_WIN"bp_outlierIDs_uniques.txt \
  $GO_DIR/"$(basename -s .tsv $GENOME_ANNOT)".background.IDs.txt \
  $GO_ANNOT --min_overlap 0.1 \
  --outfile $GO_DIR/SVs_"$POP1"_"$POP2"_RDA_"$SD"sd_outliers_overlap"$OVERLAP_WIN"bp_RDA_uniques_GO.csv
  
# 4. Filter results
MAX_FDR=0.1
MIN_LEVEL=1

Rscript 01_scripts/utils/filter_GO.R $GO_DIR/SVs_"$POP1"_"$POP2"_RDA_"$SD"sd_outliers_overlap"$OVERLAP_WIN"bp_RDA_uniques_GO.csv $MAX_FDR $MIN_LEVEL

#