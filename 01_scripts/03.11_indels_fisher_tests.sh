#!/bin/bash

# Run Fisher tests on each site
#  03.9_indels_rda.sh is a prequisite to this script, which uses the raw 012 matrix and CHROM_POS_END_ID files produced when doing RDA
# Specify window size at positional arg 1 and max allowed q value at arg 2

# manitou
# srun -p small -c 1 -J 03.11_indels_fisher_tests -o log/03.11_indels_fisher_tests_%j.log /bin/sh 01_scripts/03.11_indels_fisher_tests.sh 10000 0.01 &

# valeria
# srun -p ibis_small -c 1 -J 03.11_indels_fisher_tests -o log/03.11_indels_fisher_tests_%j.log /bin/sh 01_scripts/03.11_indels_fisher_tests.sh 10000 0.01 &

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
MAX_QVAL=$2

FISHER_DIR="11_fisher_tests/indels"

# LOAD REQUIRED MODULES
module load vcftools/0.1.16
module load bcftools/1.13
module load bedtools/2.30.0
module load R/4.1


# 1. Extract POP and ID from ID_SEX_POP
#cut -f1,3 $ID_SEX_POP > 02_infos/ID_POP.txt

# 2. Run Fisher tests
Rscript 01_scripts/utils/fisher_test.R $RDA_DIR/"$(basename -s .vcf.gz $ANGSD_FST_VCF)".geno_mat.012 02_infos/ID_POP.txt $RDA_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".CHR_POS_END_ID.table $POP1 $POP2 $MAX_QVAL $FISHER_DIR/indels_fisher_"$POP1"_"$POP2"

# 3. Get overlap of outlier sites with known genes
tail -n+2 $FISHER_DIR/indels_fisher_"$POP1"_"$POP2"_outliers_qval"$MAX_QVAL".txt > $FISHER_DIR/indels_fisher_"$POP1"_"$POP2"_outliers_qval"$MAX_QVAL".table

bedtools window -a $ANNOT_TABLE -b $FISHER_DIR/indels_fisher_"$POP1"_"$POP2"_outliers_qval"$MAX_QVAL".table -w $OVERLAP_WIN > $FISHER_DIR/indels_fisher_"$POP1"_"$POP2"_outliers_qval"$MAX_QVAL"_overlap"$OVERLAP_WIN"bp.table
