#!/bin/bash

# Perform GO enrichment analysis on outliers indels genes identified by RDA
# First get overlap between genotyped indels and known genes and between outlier indels, while allowing a 10000 bp window around genes
# Specify window size at positional arg 1 and number of sd at arg 2
# Launch in a conda env where goatools is installed

# Works on ONE population pair at the time, so variables ANGSD_FST_VCF, RAW_FST_VCF and POP_PAIR must be adjusted accordingly - I only have 2 populations (RO and PU), so VCF names will be written as is

# manitou
# srun -p small -c 1 --time=00:30:00 -J 03.10_indels_rda_go -o log/03.10_indels_rda_go_%j.log /bin/sh 01_scripts/03.10_indels_rda_go.sh 10000 3 &

# valeria
# srun -p ibis_small -c 1 --time=00:30:00 -J 03.10_indels_rda_go -o log/03.10_indels_rda_go_%j.log /bin/sh 01_scripts/03.10_indels_rda_go.sh 10000 3 &


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

CPU=2

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
SD=$2

GO_DB="12_go/go_db/go-basic.obo"
GO_ANNOT="12_go/go_db/all_go_annotations.csv"


# LOAD REQUIRED MODULES
module load vcftools/0.1.16
module load bcftools/1.13
module load R/4.1
module load bedtools/2.30.0


# 0. Get overlap of outlier sites with known genes
tail -n+2 $RDA_DIR/RDA_"$SD"sd_outliers.txt > $RDA_DIR/RDA_"$SD"sd_outliers.table
echo "$(less $RDA_DIR/RDA_"$SD"sd_outliers.table | wc -l) outlier indels"

bedtools window -a $ANNOT_TABLE -b $RDA_DIR/RDA_"$SD"sd_outliers.table -w $OVERLAP_WIN > $RDA_DIR/indels_"$POP1"_"$POP2"_outliers_RDA_"$SD"sd_overlap"$OVERLAP_WIN"bp.table

echo "$(less $RDA_DIR/indels_"$POP1"_"$POP2"_outliers_RDA_"$SD"sd_overlap"$OVERLAP_WIN"bp.table | cut -f1,5 | sort | uniq | wc -l) unique genes (or duplicated genes on different chromosomes) located at < $OVERLAP_WIN bp of an outlier indel"


# 1. Extract all known gene IDs in annotation table = BACKGROUND IDs
less $ANNOT_TABLE | cut -f5 | sort | uniq > $GO_DIR/"$(basename -s .tsv $GENOME_ANNOT)".background.IDs.txt

# 2. Extract gene IDs from outlier indels table from script 03.9 = OUTLIER IDs
less $RDA_DIR/indels_"$POP1"_"$POP2"_outliers_RDA_"$SD"sd_overlap"$OVERLAP_WIN"bp.table | cut -f5 | sort | uniq > $GO_DIR/indels_"$POP1"_"$POP2"_outliers_RDA_"$SD"sd_overlap"$OVERLAP_WIN"bp_outlierIDs.txt

# 3. Run GO enrichment
python 12_go/goatools/scripts/find_enrichment.py --pval=0.05 --indent \
  --obo $GO_DB \
  $GO_DIR/indels_"$POP1"_"$POP2"_outliers_RDA_"$SD"sd_overlap"$OVERLAP_WIN"bp_outlierIDs.txt \
  $GO_DIR/"$(basename -s .tsv $GENOME_ANNOT)".background.IDs.txt \
  $GO_ANNOT --min_overlap 0.1 \
  --outfile $GO_DIR/indels_"$POP1"_"$POP2"_RDA_"$SD"sd_outliers_overlap"$OVERLAP_WIN"bp_GO.csv
  
# 4. Filter results
MAX_FDR=0.1
MIN_LEVEL=1

Rscript 01_scripts/utils/filter_GO.R $GO_DIR/indels_"$POP1"_"$POP2"_RDA_"$SD"sd_outliers_overlap"$OVERLAP_WIN"bp_GO.csv $MAX_FDR $MIN_LEVEL

# Simplify filtered results
less $GO_DIR/indels_"$POP1"_"$POP2"_RDA_"$SD"sd_outliers_overlap"$OVERLAP_WIN"bp_GO.fdr"$MAX_FDR"_depth"$MIN_LEVEL".csv | cut -f1,4-6,8,13 | perl -pe 's/^[.]+(GO\:[0-9\.e\-]+)/\1/' > $GO_DIR/indels_"$POP1"_"$POP2"_RDA_"$SD"sd_outliers_overlap"$OVERLAP_WIN"bp_GO.fdr"$MAX_FDR"_depth"$MIN_LEVEL".simpl.txt
  
# Simplify even further for running REVIGO online
less $GO_DIR/indels_"$POP1"_"$POP2"_RDA_"$SD"sd_outliers_overlap"$OVERLAP_WIN"bp_GO.fdr"$MAX_FDR"_depth"$MIN_LEVEL".simpl.txt | tail -n+2 | cut -f1,6 | perl -pe 's/^[.]+(GO\:[0-9\.e\-]+)/\1/' > $GO_DIR/indels_"$POP1"_"$POP2"_RDA_"$SD"sd_outliers_overlap"$OVERLAP_WIN"bp_GO.fdr"$MAX_FDR"_depth"$MIN_LEVEL".GO_pval.txt


# With bedtools intersect
# with intersect
bedtools intersect -a $ANNOT_TABLE -b $RDA_DIR/RDA_"$SD"sd_outliers.table > $RDA_DIR/indels_"$POP1"_"$POP2"_outliers_RDA_"$SD"sd_intersect.table
echo "$(less $RDA_DIR/indels_"$POP1"_"$POP2"_outliers_RDA_"$SD"sd_intersect.table | cut -f1,5 | sort | uniq | wc -l) unique genes (or duplicated genes on different chromosomes) overlapping a RDA candidate indel"

less $RDA_DIR/indels_"$POP1"_"$POP2"_outliers_RDA_"$SD"sd_intersect.table | cut -f5 | sort | uniq > $GO_DIR/indels_"$POP1"_"$POP2"_outliers_RDA_"$SD"sd_intersect_outlierIDs.txt

python 12_go/goatools/scripts/find_enrichment.py --pval=0.05 --indent \
  --obo $GO_DB \
  $GO_DIR/indels_"$POP1"_"$POP2"_outliers_RDA_"$SD"sd_intersect_outlierIDs.txt \
  $GO_DIR/"$(basename -s .tsv $GENOME_ANNOT)".background.IDs.txt \
  $GO_ANNOT --min_overlap 0.1 \
  --outfile $GO_DIR/indels_"$POP1"_"$POP2"_RDA_"$SD"sd_outliers_intersect_GO.csv
  
Rscript 01_scripts/utils/filter_GO.R $GO_DIR/indels_"$POP1"_"$POP2"_RDA_"$SD"sd_outliers_intersect_GO.csv $MAX_FDR $MIN_LEVEL

less $GO_DIR/indels_"$POP1"_"$POP2"_RDA_"$SD"sd_outliers_intersect_GO.fdr"$MAX_FDR"_depth"$MIN_LEVEL".csv | cut -f1,4-6,8,13 | perl -pe 's/^[.]+(GO\:[0-9\.e\-]+)/\1/' > $GO_DIR/indels_"$POP1"_"$POP2"_RDA_"$SD"sd_outliers_intersect_GO.fdr"$MAX_FDR"_depth"$MIN_LEVEL".simpl.txt

less $GO_DIR/indels_"$POP1"_"$POP2"_RDA_"$SD"sd_outliers_intersect_GO.fdr"$MAX_FDR"_depth"$MIN_LEVEL".simpl.txt | tail -n+2 | cut -f1,6 | perl -pe 's/^[.]+(GO\:[0-9\.e\-]+)/\1/' > $GO_DIR/indels_"$POP1"_"$POP2"_RDA_"$SD"sd_outliers_intersect_GO.fdr"$MAX_FDR"_depth"$MIN_LEVEL".GO_pval.txt