#!/bin/bash

# Perform GO enrichment analysis on outliers SVs genes identified by Fisher tests
# First get overlap between genotyped SVs and known genes and between outlier SVs, while allowing a 10000 bp window around genes
# Specify window size at positional arg 1 and max allowed q value at arg 2
# Launch in a conda env where goatools is installed

# manitou
# srun -p small -c 1 -J 01.12_SVs_fisher_go -o log/01.12_SVs_fisher_go_%j.log /bin/sh 01_scripts/01.12_SVs_fisher_go.sh 10000 0.01 &

# valeria
# srun -p ibis_small -c 1 -J 01.12_SVs_fisher_go -o log/01.12_SVs_fisher_go_%j.log /bin/sh 01_scripts/01.12_SVs_fisher_go.sh 10000 0.01 &

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

#MIN_MAF=0.05
#MAX_MAF=0.95

POPS_FILE="02_infos/pops.txt"

#FILT_ANGSD_VCF="$ANGSD_STATS_DIR/"$(basename -s .vcf.gz $SNPS_VCF_ANGSD)".maf"$MIN_MAF".vcf.gz"

POP1='RO'
POP2='PU'

#ANGSD_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $FILT_ANGSD_VCF)".SNPsFst_"$POP1"_"$POP2".vcf.gz" # VCF formatted for angsd, with added Fst values from previous script
RAW_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_SNPS_VCF)".SNPsFst_"$POP1"_"$POP2".vcf.gz" # input VCF, with added Fst values from previous script

GENOME_ANNOT="03_genome/annotation/genome_annotation_table_simplified_1.5.tsv"
ANNOT_TABLE="03_genome/annotation/"$(basename -s .tsv $GENOME_ANNOT)".table"

OVERLAP_WIN=$1
MAX_QVAL=$2

FISHER_DIR="11_fisher_tests/SNPs"

GO_DB="12_go/go_db/go-basic.obo"
GO_ANNOT="12_go/go_db/all_go_annotations.csv"


# 1. Extract all known gene IDs in annotation table = BACKGROUND IDs
less $ANNOT_TABLE | cut -f5 | sort | uniq > $GO_DIR/"$(basename -s .tsv $GENOME_ANNOT)".background.IDs.txt

# 2. Extract gene IDs from outlier indels table from script 01.9 = OUTLIER IDs
less $FISHER_DIR/SNPs_fisher_"$POP1"_"$POP2"_outliers_qval"$MAX_QVAL"_overlap"$OVERLAP_WIN"bp.table | cut -f5 | sort | uniq > $GO_DIR/SNPs_fisher_"$POP1"_"$POP2"_outliers_qval"$MAX_QVAL"_overlap"$OVERLAP_WIN"bp_outlierIDs.txt

# 3. Run GO enrichment
python 12_go/goatools/scripts/find_enrichment.py --pval=0.05 --indent \
  --obo $GO_DB \
  $GO_DIR/SNPs_fisher_"$POP1"_"$POP2"_outliers_qval"$MAX_QVAL"_overlap"$OVERLAP_WIN"bp_outlierIDs.txt \
  $GO_DIR/"$(basename -s .tsv $GENOME_ANNOT)".background.IDs.txt \
  $GO_ANNOT --min_overlap 0.1 \
  --outfile $GO_DIR/SNPs_fisher_"$POP1"_"$POP2"_outliers_qval"$MAX_QVAL"_overlap"$OVERLAP_WIN"bp_GO.csv
  
# 4. Filter results
MAX_FDR=0.1
MIN_LEVEL=1

Rscript 01_scripts/utils/filter_GO.R $GO_DIR/SNPs_fisher_"$POP1"_"$POP2"_outliers_qval"$MAX_QVAL"_overlap"$OVERLAP_WIN"bp_GO.csv $MAX_FDR $MIN_LEVEL
