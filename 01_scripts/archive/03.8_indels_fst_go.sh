#!/bin/bash

# Perform GO enrichment analysis on a set of genes located near (at 10000 bp) outlier indels, given their high Fst (0.79)
# Specify window size at positional arg 1, and quantile threshold at positional arg 2
# Launch in a conda env where goatools is installed

# manitou
# srun -p small -c 1 -J 03.8_indels_fst_go -o log/03.8_indels_fst_go_%j.log /bin/sh 01_scripts/03.8_indels_fst_go.sh 10000 0.98 &

# valeria
# srun -p ibis_small -c 1 -J 03.8_indels_fst_go -o log/03.8_indels_fst_go_%j.log /bin/sh 01_scripts/03.8_indels_fst_go.sh 10000 0.98 &

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

ID_SEX_POP="02_infos/ID_Sex_Pop_updated.txt"

INDELS_VCF_ANGSD="$ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_INDELS_VCF)".recoded.vcf.gz"

N_IND="$(less $ID_SEX_POP | wc -l)"

CPU=4

#MIN_MAF=0.05
#MAX_MAF=0.95

POPS_FILE="02_infos/pops.txt"

#FILT_ANGSD_VCF="$ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $INDELS_VCF_ANGSD)".maf"$MIN_MAF".vcf.gz"


POP1='RO'
POP2='PU'

ANGSD_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .recoded.vcf.gz $INDELS_VCF_ANGSD)".indelsFst_"$POP1"_"$POP2".vcf.gz" # VCF formatted for angsd, with added Fst values from previous script
RAW_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_INDELS_VCF)".indelsFst_"$POP1"_"$POP2".vcf.gz" # input VCF (NOT the one formatted for angsd), with added Fst values from previous script

OVERLAP_WIN=$1
QUANTILE=$2

GENOME_ANNOT="03_genome/annotation/genome_annotation_table_simplified_1.5.tsv"
ANNOT_TABLE="03_genome/annotation/"$(basename -s .tsv $GENOME_ANNOT)".table"

GO_DB="12_go/go_db/go-basic.obo"
GO_ANNOT="12_go/go_db/all_go_annotations.csv"


# 1. Extract all known gene IDs in annotation table = BACKGROUND IDs
less $ANNOT_TABLE | cut -f5 | sort | uniq > $GO_DIR/"$(basename -s .tsv $GENOME_ANNOT)".background.IDs.txt

# 2. Extract gene IDs from outlier SVs table from script 01.7 = OUTLIER IDs
MIN_FST="$(less "$ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)"_quantile"$QUANTILE"_Fst_cutoff.txt" | head -n1)"

less $ANGSD_FST_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_overlap"$OVERLAP_WIN"bp.table | cut -f5 | sort | uniq > $GO_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_overlap"$OVERLAP_WIN"bp_outlierIDs.txt

# 3. Run GO enrichment
python 12_go/goatools/scripts/find_enrichment.py --pval=0.05 --indent \
  --obo $GO_DB \
  $GO_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_overlap"$OVERLAP_WIN"bp_outlierIDs.txt \
  $GO_DIR/"$(basename -s .tsv $GENOME_ANNOT)".background.IDs.txt \
  $GO_ANNOT --min_overlap 0.1 \
  --outfile $GO_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_overlap"$OVERLAP_WIN"bp_GO.csv

# 4. Filter results
MAX_FDR=0.1
MIN_LEVEL=1

Rscript 01_scripts/utils/filter_GO.R $GO_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_overlap"$OVERLAP_WIN"bp_GO.csv $MAX_FDR $MIN_LEVEL
