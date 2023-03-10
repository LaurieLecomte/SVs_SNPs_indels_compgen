#!/bin/bash

# Compute SAF and MAF

# manitou
# srun -p small -c 1 -J 01.1_SVs_format_VCF -o log/01.1_SVs_format_VCF_%j.log /bin/sh 01_scripts/01.1_SVs_format_VCF.sh &

# valeria
# srun -p ibis_small -c 1 -J 01.1_SVs_format_VCF -o log/01.1_SVs_format_VCF_%j.log /bin/sh 01_scripts/01.1_SVs_format_VCF.sh &

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
ANNOT_DIR="11_annotation/SVs"
GO_DIR="12_go/SVs"

GENOME_ANNOT="11_annotation/genome_annotation/genome_annotation_table_simplified_1.5.tsv"
ID_SEX_POP="02_infos/ID_Sex_Pop_updated.txt"

SV_VCF_ANGSD="$ANGSD_INPUT_DIR/"$(basename -s .vcf.gz $RAW_SV_VCF)".recoded.vcf"

N_IND="$(less $ID_SEX_POP | wc -l)"

CPU=4

# LOAD REQUIRED MODULES
module load angsd/0.937
module load bcftools/1.15

# 1. Calculate the SAF, MAF 
angsd -vcf-Gl $SV_VCF_ANGSD -nind $N_IND -P $CPU \
-fai "$GENOME".fai -anc $GENOME -ref $GENOME \
-domaf 1 -dosaf 1 -doMajorMinor 4 \
-out $ANGSD_STATS_DIR/"(basename -s .vcf $SV_VCF_ANGSD)"


# 2. Filter for maf >0.05 (and <0.95) as this is not a minor all freq but a reference all freq
# If we have -doMajorMinor 4, the maf column in the 7th one, not the 6th one like in original script
gunzip -c $ANGSD_STATS_DIR/"(basename -s .recoded.vcf $SV_VCF_ANGSD)".mafs.gz | awk '{ if ($7 >= 0.05 && $7<=0.95) { print } }' > $ANGSD_STATS_DIR/"(basename -s .recoded.vcf $SV_VCF_ANGSD)".maf0.05.mafs

# 3. Count number of filtered sites
wc -l $ANGSD_STATS_DIR/"(basename -s .recoded.vcf $SV_VCF_ANGSD)".maf0.05.mafs

# 4. Convert angsd output to bed
cat $ANGSD_STATS_DIR/"(basename -s .recoded.vcf $SV_VCF_ANGSD)".maf0.05.mafs | awk -v OFS='\t' '{print $1,$2-1,$2}' > $ANGSD_STATS_DIR/"(basename -s .recoded.vcf $SV_VCF_ANGSD)".maf0.05.bed
head $ANGSD_STATS_DIR/"(basename -s .recoded.vcf $SV_VCF_ANGSD)".maf0.05.bed

# 5. Output filtered sites from VCF
bgzip $SV_VCF_ANGSD
tabix "$SV_VCF_ANGSD".gz
bcftools view -R $ANGSD_STATS_DIR/"(basename -s .recoded.vcf $SV_VCF_ANGSD)".maf0.05.bed "$SV_VCF_ANGSD".gz > $ANGSD_STATS_DIR/"(basename -s .recoded.vcf $SV_VCF_ANGSD)".maf0.05.vcf 



angsd -P $NB_CPU -nQueueSize 50 \
-doMaf 1  -GL 2 -doGlf 2 -doMajorMinor 5 -doHWE 1 \
-anc 02_info/genome.fasta -remove_bads 1 -minMapQ 30 -minQ 20 -skipTriallelic 1 \
-doCounts 1 -minInd $MIN_IND -minMaf $MIN_MAF -setMaxDepth $MAX_DEPTH -setMinDepthInd $MIN_DEPTH \
-b 02_info/bam.filelist \
-r $REGIONS -out 03_saf_maf_gl_all/all_maf"$MIN_MAF"pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"_no_outliers_"$REGIONS"



