#!/bin/bash

# Get intersection of Fst and Fisher outliers to produce a single "outliers" set, and extract candidates unique to RDA (for further validation)

# Specify window size at positional arg 1, Fst quantile threshold at positional arg 2, number of SDs used for RDA at positional arg 3, and max qvalue used for Fisher tests at arg 4

# Works on ONE population pair at the time, so variables ANGSD_FST_VCF, RAW_FST_VCF and POP_PAIR must be adjusted accordingly - I only have 2 populations (RO and PU), so VCF names will be written as is

# manitou
# srun -p small -c 1 --time=00:30:00 -J 03.13_indels_compare_outliers -o log/03.13_indels_compare_outliers_%j.log /bin/sh 01_scripts/03.13_indels_compare_outliers.sh 10000 0.97 3 0.01 &

# valeria
# srun -p ibis_small -c 1 --time=00:30:00 -J 03.13_indels_compare_outliers -o log/03.13_indels_compare_outliers_%j.log /bin/sh 01_scripts/03.13_indels_compare_outliers.sh 10000 0.97 3 0.01 &

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


# 1. Find shared outliers between Fst and Fisher and candidates uniques to RDA
Rscript 01_scripts/utils/compare_outliers_candidates_sites.R $FST_OUTLIERS $RDA_OUTLIERS $FISHER_OUTLIERS $OVERLAP_WIN $QUANTILE $MIN_FST $SD $MAX_QVAL $ANGSD_FST_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared.table $RDA_DIR/RDA_"$SD"sd_outliers_uniques.table $ANGSD_FST_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_RDA_"$SD"sd_shared.table

# 2. Get overlap of outliers shared between Fst and Fishers (= confidence highly diffentiated variants set) and known genes
bedtools window -a $ANNOT_TABLE -b $ANGSD_FST_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared.table -w $OVERLAP_WIN > $ANGSD_FST_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared_"$OVERLAP_WIN"bp.table
echo "$(less $ANGSD_FST_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared_"$OVERLAP_WIN"bp.table | cut -f1,5 | sort | uniq | wc -l) unique genes (or duplicated genes on different chromosomes) located at < $OVERLAP_WIN bp of an outlier indel"

# 3. Perform GO enrichment analysis on shared outliers
# Extract all known gene IDs in annotation table = BACKGROUND IDs
less $ANNOT_TABLE | cut -f5 | sort | uniq > $GO_DIR/"$(basename -s .tsv $GENOME_ANNOT)".background.IDs.txt

# Extract gene IDs from shared outliers table
less $ANGSD_FST_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared_"$OVERLAP_WIN"bp.table | cut -f5 | sort | uniq > $GO_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared_"$OVERLAP_WIN"bp_outlierIDs.txt

# Run GO enrichment
python 12_go/goatools/scripts/find_enrichment.py --pval=0.05 --indent \
  --obo $GO_DB \
  $GO_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared_"$OVERLAP_WIN"bp_outlierIDs.txt \
  $GO_DIR/"$(basename -s .tsv $GENOME_ANNOT)".background.IDs.txt \
  $GO_ANNOT --min_overlap 0.1 \
  --outfile $GO_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared_"$OVERLAP_WIN"bp_GO.csv
  
# Filter results
MAX_FDR=0.1
MIN_LEVEL=1
Rscript 01_scripts/utils/filter_GO.R $GO_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared_"$OVERLAP_WIN"bp_GO.csv $MAX_FDR $MIN_LEVEL

# Simplify filtered results
less $GO_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared_"$OVERLAP_WIN"bp_GO.fdr"$MAX_FDR"_depth"$MIN_LEVEL".csv | cut -f1,4-6,8,13 | perl -pe 's/^[.]+(GO\:[0-9\.e\-]+)/\1/' > $GO_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared_"$OVERLAP_WIN"bp_GO.fdr"$MAX_FDR"_depth"$MIN_LEVEL".simpl.txt

# Simplify even further for running REVIGO online
less $GO_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared_"$OVERLAP_WIN"bp_GO.fdr"$MAX_FDR"_depth"$MIN_LEVEL".simpl.txt | tail -n+2 | cut -f1,6 | perl -pe 's/^[.]+(GO\:[0-9\.e\-]+)/\1/' > $GO_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared_"$OVERLAP_WIN"bp_GO.fdr"$MAX_FDR"_depth"$MIN_LEVEL".GO_pval.txt


# 4. Get overlap of outliers shared between all 3 methods and known genes
#bedtools window -a $ANNOT_TABLE -b $ANGSD_FST_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_RDA_"$SD"sd_shared.table -w $OVERLAP_WIN > $ANGSD_FST_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_RDA_"$SD"sd_shared_"$OVERLAP_WIN"bp.table
#echo "$(less $ANGSD_FST_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_RDA_"$SD"sd_shared_"$OVERLAP_WIN"bp.table | cut -f1,5 | sort | uniq | wc -l) unique genes (or duplicated genes on different chromosomes) located at < $OVERLAP_WIN bp of a candidate indel shared between Fst, Fisher and RDA"


# With bedtools intersect
#bedtools intersect -a $ANNOT_TABLE -b $ANGSD_FST_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared.table > $ANGSD_FST_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared_intersect.table
#echo "$(less $ANGSD_FST_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared_intersect.table | cut -f1,5 | sort | uniq | wc -l) unique genes (or duplicated genes on different chromosomes) overlapped by an outlier indel"

#less $ANGSD_FST_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared_intersect.table | cut -f5 | sort | uniq > $GO_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared_intersect_outlierIDs.txt

#python 12_go/goatools/scripts/find_enrichment.py --pval=0.05 --indent \
#  --obo $GO_DB \
#  $GO_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared_intersect_outlierIDs.txt \
#  $GO_DIR/"$(basename -s .tsv $GENOME_ANNOT)".background.IDs.txt \
#  $GO_ANNOT --min_overlap 0.1 \
#  --outfile $GO_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared_intersect_GO.csv

#Rscript 01_scripts/utils/filter_GO.R $GO_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared_intersect_GO.csv $MAX_FDR $MIN_LEVEL

#less $GO_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared_intersect_GO.fdr"$MAX_FDR"_depth"$MIN_LEVEL".csv | cut -f1,4-6,8,13 | perl -pe 's/^[.]+(GO\:[0-9\.e\-]+)/\1/' > $GO_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared_intersect_GO.fdr"$MAX_FDR"_depth"$MIN_LEVEL".simpl.txt

#less $GO_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared_intersect_GO.fdr"$MAX_FDR"_depth"$MIN_LEVEL".simpl.txt | tail -n+2 | cut -f1,6 | perl -pe 's/^[.]+(GO\:[0-9\.e\-]+)/\1/' > $GO_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_shared_intersect_GO.fdr"$MAX_FDR"_depth"$MIN_LEVEL".GO_pval.txt


# 5. FOR CHI SQUARE TESTS : Get overlap of outliers shared between Fst and Fisher but NOT in RDA candidates and known genes (and RDA candidates not in inersection outliers set)
## the files *_outliers_set.table and RDA_cand_set.table" are not user-specified outputs, they are produced from input file names
#bedtools window -a $ANNOT_TABLE -b $ANGSD_FST_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_outliers_set.table -w $OVERLAP_WIN > $ANGSD_FST_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_outliers_set_"$OVERLAP_WIN"bp.table
#echo "$(less $ANGSD_FST_DIR/indels_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_qval"$MAX_QVAL"_outliers_set_"$OVERLAP_WIN"bp.table | cut -f15-18 | sort | uniq | wc -l) intersection outliers that are not RDA candidates are located near a known gene"

#bedtools window -a $ANNOT_TABLE -b "$RDA_DIR/RDA_"$SD"sd_outliers_RDA_cand_set.table" -w $OVERLAP_WIN > $RDA_DIR/RDA_"$SD"sd_outliers_RDA_cand_set_"$OVERLAP_WIN"bp.table
#echo "$(less $RDA_DIR/RDA_"$SD"sd_outliers_RDA_cand_set_"$OVERLAP_WIN"bp.table | cut -f15-18 | sort | uniq | wc -l) RDA candidates that are not in intersection outliers are located near a known gene"

# 4. Get overlap of RDA unique candidates with known genes
#bedtools window -a $ANNOT_TABLE -b $RDA_DIR/RDA_"$SD"sd_outliers_uniques.table -w $OVERLAP_WIN > $RDA_DIR/RDA_"$SD"sd_outliers_uniques_"$OVERLAP_WIN"bp.table

## Extract gene IDs from shared outliers table
#less $RDA_DIR/RDA_"$SD"sd_outliers_uniques_"$OVERLAP_WIN"bp.table | cut -f5 | sort | uniq > $GO_DIR/RDA_"$SD"sd_outliers_uniques_"$OVERLAP_WIN"bp_outlierIDs.txt

## 5. Perform GO enrichment analysis on RDA unique candidates
#python 12_go/goatools/scripts/find_enrichment.py --pval=0.05 --indent \
#  --obo $GO_DB \
#  $GO_DIR/RDA_"$SD"sd_outliers_uniques_"$OVERLAP_WIN"bp_outlierIDs.txt \
#  $GO_DIR/"$(basename -s .tsv $GENOME_ANNOT)".background.IDs.txt \
#  $GO_ANNOT --min_overlap 0.1 \
#  --outfile $GO_DIR/RDA_"$SD"sd_outliers_uniques_"$OVERLAP_WIN"bp_GO.csv
#  
## Filter results
#MAX_FDR=0.1
#MIN_LEVEL=1
#Rscript 01_scripts/utils/filter_GO.R $GO_DIR/RDA_"$SD"sd_outliers_uniques_"$OVERLAP_WIN"bp_GO.csv $MAX_FDR $MIN_LEVEL

## Simplify filtered results
#less $GO_DIR/RDA_"$SD"sd_outliers_uniques_"$OVERLAP_WIN"bp_GO.fdr"$MAX_FDR"_depth"$MIN_LEVEL".csv | cut -f1,4-6,8,13 > $GO_DIR/RDA_"$SD"sd_outliers_uniques_"$OVERLAP_WIN"bp_GO.fdr"$MAX_FDR"_depth"$MIN_LEVEL".simpl.tx#t