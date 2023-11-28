#!/bin/bash


# VARIABLES
ALL_SV="04_vcf/SVs/merged_SUPP2_MAF0.05_FMISS0.5_CHROM_POS_END.table"
ALL_SNP="04_vcf/SNPs/SNPs_MAF0.05_FMISS0.5_CHROM_POS_END.table"
ALL_INDEL="04_vcf/indels/indels_MAF0.05_FMISS0.5_CHROM_POS_END.table"

SIMPL_TABLE="03_genome/annotation/filtered_elements.table"
OVERLAP_DIR="overlap"
MIN_OVERLAP=0.5 # minimum fraction of genome element that overlaps with SV to be called


# LOAD REQUIRED MODULES
module load bedtools

#
bedtools intersect -a $SIMPL_TABLE -b $ALL_SV -wa -wb -f $MIN_OVERLAP > $OVERLAP_DIR/$"(basename -s CHROM_POS_END.table $ALL_SV)"_overlap.table

bedtools intersect -a $SIMPL_TABLE -b $ALL_SNP -wa -wb -f $MIN_OVERLAP > $OVERLAP_DIR/$"(basename -s CHROM_POS_END.table $ALL_SNP)"_overlap.table

bedtools intersect -a $SIMPL_TABLE -b $ALL_INDEL -wa -wb -f $MIN_OVERLAP > $OVERLAP_DIR/$"(basename -s CHROM_POS_END.table $ALL_INDEL)"_overlap.table
