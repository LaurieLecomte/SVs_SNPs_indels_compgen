#!/bin/bash

# Perform GO enrichment analysis on a set of genes located near (at 10000 bp) outlier SVs, given their high Fst (0.4)
# Specify window size at positional arg 1, and min Fst threshold at positional arg 2
# Launch in a conda env where goatools is installed

# manitou
# srun -p small -c 1 -J 01.8_SVs_fst_go -o log/01.8_SVs_fst_go_%j.log /bin/sh 01_scripts/01.8_SVs_fst_go.sh 10000 0.44 &

# valeria
# srun -p ibis_small -c 1 -J 01.8_SVs_fst_go -o log/01.8_SVs_fst_go_%j.log /bin/sh 01_scripts/01.8_SVs_fst_go.sh 10000 0.44 &

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

CPU=4

MIN_MAF=0.05
MAX_MAF=0.95

FILT_ANGSD_VCF="$ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".vcf.gz"


POP1='RO'
POP2='PU'

ANGSD_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $FILT_ANGSD_VCF)".SVsFst_"$POP1"_"$POP2".vcf.gz" # VCF formatted for angsd, with added Fst values from previous script
RAW_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".SVsFst_"$POP1"_"$POP2".vcf.gz" # input VCF (NOT the one formatted for angsd), with added Fst values from previous script

OVERLAP_WIN=$1
MIN_FST=$2

GENOME_ANNOT="03_genome/annotation/genome_annotation_table_simplified_1.5.tsv"
ANNOT_TABLE="03_genome/annotation/"$(basename -s .tsv $GENOME_ANNOT)".table"

GO_DB="12_go/go_db/go-basic.obo"
GO_ANNOT="12_go/go_db/all_go_annotations.csv"


# 1. Extract all known gene IDs in annotation table = BACKGROUND IDs
less $ANNOT_TABLE | cut -f5 | sort | uniq > $GO_DIR/"$(basename -s .tsv $GENOME_ANNOT)".background.IDs.txt

# 2. Extract gene IDs from outlier SVs table from script 01.7 = OUTLIER IDs
less $ANGSD_FST_DIR/SVs_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_overlap"$OVERLAP_WIN"bp.table | cut -f5 | sort | uniq > $GO_DIR/SVs_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_overlap"$OVERLAP_WIN"bp_outlierIDs.txt

# 3. Run GO enrichment
python 12_go/goatools/scripts/find_enrichment.py --pval=0.05 --indent \
  --obo $GO_DB \
  $GO_DIR/SVs_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_overlap"$OVERLAP_WIN"bp_outlierIDs.txt \
  $GO_DIR/"$(basename -s .tsv $GENOME_ANNOT)".background.IDs.txt \
  $GO_ANNOT --min_overlap 0.1 \
  --outfile $GO_DIR/SVs_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_overlap"$OVERLAP_WIN"bp_GO.csv

#python 15_go/goatools/scripts/find_enrichment.py --pval=0.05 --indent \
#  --obo $GO_DB \
#  $ANNOT_DIR/outlier_SVs_fixed"$FIXED_FST"/outlier_SVs_fixed"$FIXED_FST"_"$DIST".overlap.IDs \
#  03_genome/all_genes.ids \
#  $GO_ANNOT --min_overlap 0.1 \
#  --outfile $GO_DIR/outlier_SVs_fixed"$FIXED_FST"/outlier_SVs_fixed"$FIXED_FST"_win"$DIST"_vs_ALLgenes.csv

# 4. Filter results
MAX_FDR=0.1
MIN_LEVEL=1

Rscript 01_scripts/utils/filter_GO.R $GO_DIR/SVs_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_overlap"$OVERLAP_WIN"bp_GO.csv $MAX_FDR $MIN_LEVEL