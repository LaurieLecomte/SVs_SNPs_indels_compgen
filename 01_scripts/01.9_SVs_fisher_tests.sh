#!/bin/bash

# Run Fisher tests on each site
#  01.9_SVs_rda.sh is a prequisite to this script, which uses the raw 012 matrix and CHROM_POS_END_ID files produced when doing RDA 
# Specify window size at positional arg 1 and max allowed q value at arg 2

# manitou
# srun -p small -c 1 -J 01.11_SVs_fisher_tests -o log/01.11_SVs_fisher_tests_%j.log /bin/sh 01_scripts/01.11_SVs_fisher_tests.sh 10000 0.01 &

# valeria
# srun -p ibis_small -c 1 -J 01.11_SVs_fisher_tests -o log/01.11_SVs_fisher_tests_%j.log /bin/sh 01_scripts/01.11_SVs_fisher_tests.sh 10000 0.01 &


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

GENOME_ANNOT="03_genome/annotation/genome_annotation_table_simplified_1.5.tsv"
ANNOT_TABLE="03_genome/annotation/"$(basename -s .tsv $GENOME_ANNOT)".table"

OVERLAP_WIN=$1

MAX_QVAL=$2

FISHER_DIR="11_fisher_tests/SVs"

# LOAD REQUIRED MODULES
module load vcftools/0.1.16
module load bcftools/1.13
module load bedtools/2.30.0
module load R/4.1


# 1. Extract POP and ID from ID_SEX_POP
cut -f1,3 $ID_SEX_POP > 02_infos/ID_POP.txt

# 2. Run Fisher tests
#Rscript 01_scripts/utils/fisher_test.R $RDA_DIR/"$(basename -s .vcf.gz $ANGSD_FST_VCF)".geno_mat.012 02_infos/ID_POP.txt $RDA_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".CHR_POS_END_ID.table $POP1 $POP2 $MAX_QVAL $FISHER_DIR/SVs_fisher_"$POP1"_"$POP2"
Rscript 01_scripts/utils/filter_fisher.R $RDA_DIR/"$(basename -s .vcf.gz $ANGSD_FST_VCF)".geno_mat.012 02_infos/ID_POP.txt $RDA_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".CHR_POS_END_ID.table $POP1 $POP2 $MAX_QVAL $FISHER_DIR/SVs_fisher_"$POP1"_"$POP2" $ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".table

# 3. Get overlap of outlier sites with known genes
tail -n+2 $FISHER_DIR/SVs_fisher_"$POP1"_"$POP2"_outliers_qval"$MAX_QVAL".txt > $FISHER_DIR/SVs_fisher_"$POP1"_"$POP2"_outliers_qval"$MAX_QVAL".table
echo "$(less $FISHER_DIR/SVs_fisher_"$POP1"_"$POP2"_outliers_qval"$MAX_QVAL".table | wc -l) outlier SVs with Fisher qval < $MAX_QVAL"

bedtools window -a $ANNOT_TABLE -b $FISHER_DIR/SVs_fisher_"$POP1"_"$POP2"_outliers_qval"$MAX_QVAL".table -w $OVERLAP_WIN > $FISHER_DIR/SVs_fisher_"$POP1"_"$POP2"_outliers_qval"$MAX_QVAL"_overlap"$OVERLAP_WIN"bp.table
echo "$(less $FISHER_DIR/SVs_fisher_"$POP1"_"$POP2"_outliers_qval"$MAX_QVAL"_overlap"$OVERLAP_WIN"bp.table | cut -f1,5 | sort | uniq | wc -l) unique genes (or duplicated genes on different chromosomes) located at < $OVERLAP_WIN bp of an outlier SV with Fisher qval < $MAX_QVAL"