#!/bin/bash

# Compute SAF and MAF

# manitou
# srun -p small -c 4 -J 02.1_SNPs_saf_maf_all -o log/02.1_SNPs_saf_maf_all_%j.log /bin/sh 01_scripts/02.1_SNPs_saf_maf_all.sh &

# valeria
# srun -p ibis_small -c 4 -J 02.1_SNPs_saf_maf_all -o log/02.1_SNPs_saf_maf_all_%j.log /bin/sh 01_scripts/02.1_SNPs_saf_maf_all.sh &

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
ANNOT_DIR="11_annotation/SNPs"
GO_DIR="12_go/SNPs"

GENOME_ANNOT="11_annotation/genome_annotation/genome_annotation_table_simplified_1.5.tsv"
ID_SEX_POP="02_infos/ID_Sex_Pop_updated.txt"

#SNPS_VCF_ANGSD="$ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SNPS_VCF)".recoded.vcf.gz"
SNPS_VCF_ANGSD=$RAW_SNPS_VCF

N_IND="$(less $ID_SEX_POP | wc -l)"

CPU=4

MIN_MAF=0.05
MAX_MAF=0.95

# LOAD REQUIRED MODULES
module load python/3.7
module load angsd/0.937
module load pcangsd/1.10
module load bcftools/1.13

# 1. Calculate saf and maf
angsd -vcf-Pl $SNPS_VCF_ANGSD -nind $N_IND -P $CPU \
-fai "$GENOME".fai -anc $GENOME -ref $GENOME \
-domaf 1 -dosaf 1 -doMajorMinor 5 \
-out $ANGSD_STATS_DIR/"$(basename -s .vcf.gz $SNPS_VCF_ANGSD)"

#-doMajorMinor	0
#	1: Infer major and minor from GL
#	2: Infer major and minor from allele counts
#	3: use major and minor from a file (requires -sites file.txt)
#	4: Use reference allele as major (requires -ref)
#	5: Use ancestral allele as major (requires -anc)

#-> Inputtype is vcf/bcf
#        -> VCF still beta. Remember that
#           1. SNPs are are discarded
#           2. will use chrom, pos PL columns
#           3. GL tags are interpreted as log10 and are scaled to ln (NOT USED)
#           4. GP tags are interpreted directly as unscaled post probs (spec says phredscaled...) (NOT USED)
#           5. FILTER column is currently NOT used (not sure what concensus is)
#           6. -sites does NOT work with vcf input but -r does
#           7. vcffilereading is still BETA, please report strange behaviour
#           8. Consider adding '-doMajorMinor 1'
#        -> Reading fasta: 03_genome/genome.fasta
#        -> Reading fasta: 03_genome/genome.fasta
#        -> No indel tag in vcf/bcf file, will therefore not be able to filter out SNPs
        
# 2. Filter for maf >0.05 (and <0.95) as this is not a minor all freq but a reference all freq
# If we have -doMajorMinor 4, the maf column in the 7th one, not the 6th one like in original script
#zless $ANGSD_STATS_DIR/"$(basename -s .vcf.gz $SNPS_VCF_ANGSD)".mafs.gz | awk '{ if ($7 >= 0.05 && $7 <= 0.95) { print } }' > $ANGSD_STATS_DIR/"$(basename -s .vcf.gz $SNPS_VCF_ANGSD)".maf0.05.mafs

zless $ANGSD_STATS_DIR/"$(basename -s .vcf.gz $SNPS_VCF_ANGSD)".mafs.gz | awk -v a=$MIN_MAF -v b=$MAX_MAF '{ if ($7 >= a && $7 <= b) { print } }' > $ANGSD_STATS_DIR/"$(basename -s .vcf.gz $SNPS_VCF_ANGSD)".maf"$MIN_MAF".mafs


# 3. Count number of filtered sites
wc -l $ANGSD_STATS_DIR/"$(basename -s .vcf.gz $SNPS_VCF_ANGSD)".maf"$MIN_MAF".mafs

# 4. Convert angsd output to bed
cat $ANGSD_STATS_DIR/"$(basename -s .vcf.gz $SNPS_VCF_ANGSD)".maf"$MIN_MAF".mafs | awk -v OFS='\t' '{print $1,$2-1,$2}' > $ANGSD_STATS_DIR/"$(basename -s .vcf.gz $SNPS_VCF_ANGSD)".maf"$MIN_MAF".bed
head $ANGSD_STATS_DIR/"$(basename -s .vcf.gz $SNPS_VCF_ANGSD)".maf"$MIN_MAF".bed

# 5. Output filtered sites from VCF
#bgzip $SNPS_VCF_ANGSD
#tabix "$SNPS_VCF_ANGSD".gz
bcftools view -R $ANGSD_STATS_DIR/"$(basename -s .vcf.gz $SNPS_VCF_ANGSD)".maf"$MIN_MAF".bed $SNPS_VCF_ANGSD | bcftools sort -Oz > $ANGSD_STATS_DIR/"$(basename -s .vcf.gz $SNPS_VCF_ANGSD)".maf"$MIN_MAF".vcf.gz