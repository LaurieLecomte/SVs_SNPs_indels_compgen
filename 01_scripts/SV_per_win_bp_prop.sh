#!/bin/bash

# Compute bp covered by SVs for a given window size, specified at positional arg 1
# 01.9_SVs_rda.sh is a prequisite to this script, which uses the raw 012 matrix and CHROM_POS_END_ID files produced when doing RDA

# PREREQUISITE : matching genotyped SVs with known candidates based on POS and ALT seq len in the genotype_SVs_SRLR pipeline (https://github.com/LaurieLecomte/genotype_SVs_SRLR/blob/main/01_scripts/utils/compare_summarize_plot.sh)
# Copy output *_matched_offset5bp.txt to 04_vcf/SVs directory of current pipeline (SVs_SNPs_indels_compgen) 

# manitou
# srun -p small -c 1 -J SV_per_win_bp_prop -o log/SV_per_win_bp_prop_%j.log /bin/sh 01_scripts/utils/SV_per_win_bp_prop.sh 100000 &

# valeria
# srun -p ibis_small -c 1 -J SV_per_win_bp_prop -o log/SV_per_win_bp_prop_%j.log /bin/sh 01_scripts/utils/SV_per_win_bp_prop.sh 100000 &

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

MIN_MAF=0.05
MAX_MAF=0.95

FILT_ANGSD_VCF="$ANGSD_STATS_DIR/"$(basename -s .recoded.vcf.gz $SV_VCF_ANGSD)".maf"$MIN_MAF".vcf.gz"

POP1='RO'
POP2='PU'

ANGSD_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $FILT_ANGSD_VCF)".SVsFst_"$POP1"_"$POP2".vcf.gz" # VCF formatted for angsd, with added Fst values from previous script
RAW_FST_VCF="$ANGSD_FST_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".SVsFst_"$POP1"_"$POP2".vcf.gz" # input VCF (NOT the one formatted for angsd), with added Fst values from previous script

GENOME_ANNOT="03_genome/annotation/genome_annotation_table_simplified_1.5.tsv"
ANNOT_TABLE="03_genome/annotation/"$(basename -s .tsv $GENOME_ANNOT)".table"

WIN_SIZE=$1

CHR_BED="02_infos/chrs.bed"

DENSITY_DIR="density/SVs"

MATCHED_SV="$RAW_VCF_DIR/merged_SUPP2_MAF0.05_FMISS0.5_matched_offset5bp.txt" # this file is produced by the script compare_summarize_plot.sh from the gnotype_SVs_SRLR pipeline (https://github.com/LaurieLecomte/genotype_SVs_SRLR/blob/main/01_scripts/utils/compare_summarize_plot.sh)



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
bedtools intersect -a 02_infos/chrs_win"$WIN_SIZE".bed -b $RDA_DIR/"$(basename -s .vcf.gz $RAW_FST_VCF)".CHR_POS_END_ID.table -wa -wb > $DENSITY_DIR/SVs_intersect_win"$WIN_SIZE"bp.txt

# 4. Compute number of bp and prop
Rscript 01_scripts/utils/SVs_per_win_bp.R $DENSITY_DIR/SVs_intersect_win"$WIN_SIZE"bp.txt 02_infos/OV_to_ssa.txt $WIN_SIZE $MATCHED_SV 