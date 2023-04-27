#!/bin/bash

# Get overlap between genotyped SVs and known genes and between outlier SVs (with Fst > 0.4), while allowing a 10000 bp window around genes
# Specify window size at positional arg 1, and min Fst threshold at positional arg 2

# Works on ONE population pair at the time, so variables ANGSD_FST_VCF, RAW_FST_VCF and POP_PAIR must be adjusted accordingly - I only have 2 populations (RO and PU), so VCF names will be written as is

# manitou
# srun -p small -c 1 -J 01.7_SVs_fst_outliers_overlap -o log/01.7_SVs_fst_outliers_overlap_%j.log /bin/sh 01_scripts/01.7_SVs_fst_outliers_overlap.sh 10000 0.44 &

# valeria
# srun -p ibis_small -c 1 -J 01.7_SVs_fst_outliers_overlap -o log/01.7_SVs_fst_outliers_overlap_%j.log /bin/sh 01_scripts/01.7_SVs_fst_outliers_overlap.sh 10000 0.44 &

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

#GENOME_ANNOT="11_annotation/genome_annotation/genome_annotation_table_simplified_1.5.tsv"
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
MIN_FST=$2

GENOME_ANNOT="03_genome/annotation/genome_annotation_table_simplified_1.5.tsv"
ANNOT_TABLE="03_genome/annotation/"$(basename -s .tsv $GENOME_ANNOT)".table"


# LOAD REQUIRED MODULES
module load bedtools/2.30.0
module load bcftools/1.13

# 1. Extract CHROM, POS, END and FST from final ANGSD VCF
#bcftools query -f "%CHROM\t%POS\t%END\t%ID\t%FST_"$POP1"_"$POP2"\n" $ANGSD_FST_VCF > $ANGSD_FST_DIR/"$(basename -s .vcf.gz $ANGSD_FST_VCF)".table
bcftools query -f "%CHROM\t%POS\t%END\t%ID\t%FST_"$POP1"_"$POP2"\n" $RAW_FST_VCF > $ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".table

# 2. Find overlap between ALL genotyped SVs and known genes = get number of unique genes near SVs
#bedtools window -a $ANNOT_TABLE -b $ANGSD_FST_DIR/"$(basename -s .vcf.gz $ANGSD_FST_VCF)".table -w $OVERLAP_WIN > $ANGSD_FST_DIR/"$(basename -s .vcf.gz $ANGSD_FST_VCF)"_overlap"$OVERLAP_WIN"bp_allSVs.table
bedtools window -a $ANNOT_TABLE -b $ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".table -w $OVERLAP_WIN > $ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)"_overlap"$OVERLAP_WIN"bp_allSVs.table

#echo "$(less $ANGSD_FST_DIR/"$(basename -s .vcf.gz $ANGSD_FST_VCF)"_overlap"$OVERLAP_WIN"bp_allSVs.table | cut -f1,5 | sort | uniq | wc -l) unique genes (or duplicated genes on different chromosomes) located at < $OVERLAP_WIN bp of a filtered genotyped SV"
echo "$(less $ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)"_overlap"$OVERLAP_WIN"bp_allSVs.table | cut -f1,5 | sort | uniq | wc -l) unique genes (or duplicated genes on different chromosomes) located at < $OVERLAP_WIN bp of a filtered genotyped SV"

# 3. Find overlap between OUTLIER SVs and known genes
## Find outlier SVs with Fst > MIN_FST
#less $ANGSD_FST_DIR/"$(basename -s .vcf.gz $ANGSD_FST_VCF)".table | awk -v val="$MIN_FST" 'BEGIN{FS="\t"} $5 >= val {print}' > $ANGSD_FST_DIR/SVs_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST".table
less $ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".table | awk -v val="$MIN_FST" 'BEGIN{FS="\t"} $5 >= val {print}' > $ANGSD_FST_DIR/SVs_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST".table
echo "$(less $ANGSD_FST_DIR/SVs_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST".table | wc -l) outlier SVs with Fst >= $MIN_FST"

## Get overlap 
bedtools window -a $ANNOT_TABLE -b $ANGSD_FST_DIR/SVs_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST".table -w $OVERLAP_WIN > $ANGSD_FST_DIR/SVs_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_overlap"$OVERLAP_WIN"bp.table

echo "$(less $ANGSD_FST_DIR/SVs_"$POP1"_"$POP2"_outliers_minFst"$MIN_FST"_overlap"$OVERLAP_WIN"bp.table | cut -f1,5 | sort | uniq | wc -l) unique genes (or duplicated genes on different chromosomes) located at < $OVERLAP_WIN bp of an outlier SV with Fst >= $MIN_FST"



