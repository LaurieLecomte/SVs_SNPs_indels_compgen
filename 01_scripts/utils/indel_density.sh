#!/bin/bash

# Compute SV density for a given window size, specified at positional arg 1
# 03.9_indels_rda.sh is a prequisite to this script, which uses the raw 012 matrix and CHROM_POS_END_ID files produced when doing RDA

# manitou
# srun -p small -c 1 -J indel_density -o log/indel_density_%j.log /bin/sh 01_scripts/utils/indel_density.sh 1000000 &

# valeria
# srun -p ibis_small -c 1 -J indel_density -o log/indel_density_%j.log /bin/sh 01_scripts/utils/indel_density.sh 1000000 &

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

MIN_MAF=0.05
MAX_MAF=0.95

POPS_FILE="02_infos/pops.txt"

FILT_ANGSD_VCF="$ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $INDELS_VCF_ANGSD)".maf"$MIN_MAF".vcf.gz"

POP1='RO'
POP2='PU'

ANGSD_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $FILT_ANGSD_VCF)".indelsFst_"$POP1"_"$POP2".vcf.gz" # VCF formatted for angsd, with added Fst values from previous script
RAW_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_INDELS_VCF)".indelsFst_"$POP1"_"$POP2".vcf.gz" # input VCF (NOT the one formatted for angsd), with added Fst values from previous script

GENOME_ANNOT="03_genome/annotation/genome_annotation_table_simplified_1.5.tsv"
ANNOT_TABLE="03_genome/annotation/"$(basename -s .tsv $GENOME_ANNOT)".table"

WIN_SIZE=$1

CHR_BED="02_infos/chrs.bed"

DENSITY_DIR="density/indels"


# LOAD REQUIRED MODULES
module load vcftools/0.1.16
module load bcftools/1.13
module load bedtools/2.30.0
module load bedops/2.4.40
module load R/4.1


# 1. Remove header from CHROM_POS_END_ID file (previously used for RDA and Fisher tests)
#less $RDA_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".CHR_POS_END_ID.table | tail -n+2 > $RDA_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".CHR_POS_END_ID.bed

# 2. Split genome by window
bedops --chop $WIN_SIZE -x $CHR_BED > 02_infos/chrs_win"$WIN_SIZE".bed

# 3. Get number of overlap between fixed size windows and sites 
bedtools intersect -a 02_infos/chrs_win"$WIN_SIZE".bed -b $RDA_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".CHR_POS_END_ID.table -c > $DENSITY_DIR/indels_density_win"$WIN_SIZE"bp.txt

# 4. Compute mean density by window and plot
Rscript 01_scripts/utils/compute_plot_density.R $DENSITY_DIR/indels_density_win"$WIN_SIZE"bp.txt 02_infos/OV_to_ssa.txt $WIN_SIZE